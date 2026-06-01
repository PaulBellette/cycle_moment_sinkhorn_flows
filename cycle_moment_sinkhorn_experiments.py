#!/usr/bin/env -S uv run
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "torch>=2.2",
#   "numpy>=1.26",
#   "matplotlib>=3.8",
# ]
# ///
"""
Cycle-moment Sinkhorn flow experiments.

Run:
    uv run cycle_moment_sinkhorn_experiments.py

Outputs plots into:
    ./cycle_moment_outputs/

Experiments:
  A. Friendly circle TSP sanity check with brute force optimum.
  B. Constraint-only primal-dual ringing: rho=0 vs augmented rho>0.
  C. Uniform initialisation degeneracy: trace-gradient norm vs perturbation scale.
  D. Two-cluster TSP-like instance: fixed/scheduled penalty comparison.
  E. Rho-pressure sweep and rho/tau phase diagram.

The code intentionally favours clarity over efficiency. Keep n small.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import torch

Tensor = torch.Tensor

torch.set_default_dtype(torch.float64)


# -----------------------------------------------------------------------------
# Basic utilities
# -----------------------------------------------------------------------------


def set_seed(seed: int) -> None:
    np.random.seed(seed)
    torch.manual_seed(seed)


@dataclass
class SinkhornConfig:
    iters: int = 80
    eps: float = 1e-12


@dataclass
class RunConfig:
    steps: int = 2000
    lr: float = 0.05
    rho: float = 10.0
    tau: float = 1.0
    tau_final: Optional[float] = None
    rho_final: Optional[float] = None
    sinkhorn_iters: int = 80
    seed: int = 0
    log_every: int = 5
    no_self_loops: bool = True
    diag_penalty: float = 1e6
    device: str = "cpu"


@dataclass
class History:
    step: List[int]
    energy: List[float]
    cost: List[float]
    h_norm: List[float]
    h_vals: List[np.ndarray]
    entropy: List[float]
    row_err: List[float]
    col_err: List[float]
    sharpness: List[float]
    tour_cost: List[float]
    best_perm_cost: List[float]

    @staticmethod
    def empty() -> "History":
        return History([], [], [], [], [], [], [], [], [], [], [])


def sinkhorn_from_logits(X: Tensor, tau: float = 1.0, cfg: SinkhornConfig = SinkhornConfig()) -> Tensor:
    """Return Sinkhorn(exp(X/tau)) using log-domain row/column normalisation."""
    Z = X / tau
    # Avoid self-loop by allowing caller to put huge negative values on diagonal.
    for _ in range(cfg.iters):
        Z = Z - torch.logsumexp(Z, dim=1, keepdim=True)
        Z = Z - torch.logsumexp(Z, dim=0, keepdim=True)
    return torch.exp(Z).clamp_min(cfg.eps)


def apply_no_self_loop_mask(X: Tensor, mask_value: float = -50.0) -> Tensor:
    """Softly forbid diagonal entries in logit space."""
    n = X.shape[0]
    Y = X.clone()
    idx = torch.arange(n, device=X.device)
    Y[idx, idx] = mask_value
    return Y


def trace_moments(P: Tensor) -> Tensor:
    """h_k = tr(P^k) - b_k, k=1..n, with b_n=n and b_k=0 otherwise."""
    n = P.shape[0]
    vals = []
    Pk = torch.eye(n, dtype=P.dtype, device=P.device)
    for k in range(1, n + 1):
        Pk = Pk @ P
        target = float(n) if k == n else 0.0
        vals.append(torch.trace(Pk) - target)
    return torch.stack(vals)


def tsp_cost(P: Tensor, C: Tensor) -> Tensor:
    return torch.sum(C * P)


def entropy(P: Tensor) -> Tensor:
    return -torch.sum(P * torch.log(P.clamp_min(1e-12)))


def sharpness(P: Tensor) -> Tensor:
    """Mean row max. 1 means permutation-like; 1/n means uniform-ish."""
    return torch.mean(torch.max(P, dim=1).values)


def greedy_perm_from_P(P: Tensor) -> List[int]:
    """Greedy one-out assignment. Not guaranteed optimal assignment, just diagnostic."""
    n = P.shape[0]
    P_np = P.detach().cpu().numpy()
    remaining_cols = set(range(n))
    perm = [-1] * n
    # Assign rows with most confident max first.
    row_order = list(np.argsort(-P_np.max(axis=1)))
    for i in row_order:
        cols = sorted(remaining_cols, key=lambda j: -P_np[i, j])
        if cols:
            j = cols[0]
            perm[i] = j
            remaining_cols.remove(j)
    # Fill any missing just in case.
    for i in range(n):
        if perm[i] < 0:
            perm[i] = remaining_cols.pop()
    return perm


def assignment_perm_from_P(P: Tensor, max_bruteforce_n: int = 9) -> List[int]:
    """Exact assignment projection: maximise sum_i P[i, perm[i]].

    For the small diagnostic problems in this note, brute force is simpler than
    adding SciPy as a dependency. Falls back to greedy if n is too large.
    """
    n = P.shape[0]
    if n > max_bruteforce_n:
        return greedy_perm_from_P(P)
    P_np = P.detach().cpu().numpy()
    best_score = -math.inf
    best_perm: Tuple[int, ...] = tuple(range(n))
    for perm in itertools.permutations(range(n)):
        score = float(sum(P_np[i, perm[i]] for i in range(n)))
        if score > best_score:
            best_score = score
            best_perm = perm
    return list(best_perm)


def projection_perm_from_P(P: Tensor, method: str = "assignment") -> List[int]:
    if method == "greedy":
        return greedy_perm_from_P(P)
    if method == "assignment":
        return assignment_perm_from_P(P)
    raise ValueError(f"unknown projection method: {method}")


def perm_cost(perm: List[int], C_np: np.ndarray) -> float:
    return float(sum(C_np[i, perm[i]] for i in range(len(perm))))


def count_cycles(perm: List[int]) -> List[List[int]]:
    n = len(perm)
    seen = [False] * n
    cycles = []
    for s in range(n):
        if seen[s]:
            continue
        cur = s
        cyc = []
        while not seen[cur]:
            seen[cur] = True
            cyc.append(cur)
            cur = perm[cur]
        cycles.append(cyc)
    return cycles


def brute_force_atsp(C_np: np.ndarray) -> Tuple[float, Tuple[int, ...]]:
    """Brute force directed Hamiltonian cycles fixing start node 0.

    A tour is 0 -> p1 -> ... -> p_{n-1} -> 0.
    Returns best cost and node order including 0 first.
    """
    n = C_np.shape[0]
    best = math.inf
    best_order: Tuple[int, ...] = tuple(range(n))
    for rest in itertools.permutations(range(1, n)):
        order = (0,) + rest
        cost = 0.0
        for a, b in zip(order, order[1:]):
            cost += C_np[a, b]
        cost += C_np[order[-1], 0]
        if cost < best:
            best = cost
            best_order = order
    return best, best_order


def order_to_perm(order: Iterable[int]) -> List[int]:
    order = list(order)
    n = len(order)
    p = [-1] * n
    for a, b in zip(order, order[1:]):
        p[a] = b
    p[order[-1]] = order[0]
    return p


# -----------------------------------------------------------------------------
# Instance generation
# -----------------------------------------------------------------------------


def circle_instance(n: int, diag_penalty: float = 1e6) -> Tuple[np.ndarray, np.ndarray]:
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    pts = np.stack([np.cos(theta), np.sin(theta)], axis=1)
    C = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=-1)
    np.fill_diagonal(C, diag_penalty)
    return pts, C


def two_cluster_instance(n: int, sep: float = 5.0, sigma: float = 0.35, seed: int = 0, diag_penalty: float = 1e6) -> Tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    n1 = n // 2
    n2 = n - n1
    pts1 = rng.normal(loc=(-sep / 2, 0.0), scale=sigma, size=(n1, 2))
    pts2 = rng.normal(loc=(sep / 2, 0.0), scale=sigma, size=(n2, 2))
    pts = np.concatenate([pts1, pts2], axis=0)
    C = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=-1)
    np.fill_diagonal(C, diag_penalty)
    return pts, C


def zero_cost_instance(n: int, diag_penalty: float = 1e6) -> Tuple[np.ndarray, np.ndarray]:
    pts = np.zeros((n, 2))
    C = np.zeros((n, n), dtype=float)
    np.fill_diagonal(C, diag_penalty)
    return pts, C


# -----------------------------------------------------------------------------
# Optimisers / flows
# -----------------------------------------------------------------------------


def scheduled_value(t: int, steps: int, start: float, final: Optional[float], mode: str = "linear") -> float:
    if final is None:
        return start
    s = t / max(steps - 1, 1)
    if mode == "log":
        return float(math.exp((1 - s) * math.log(start) + s * math.log(final)))
    return float((1 - s) * start + s * final)


def trace_penalty_energy(P: Tensor, C: Tensor, rho: float) -> Tuple[Tensor, Tensor, Tensor]:
    h = trace_moments(P)
    cost = tsp_cost(P, C)
    E = cost + 0.5 * rho * torch.sum(h * h)
    return E, cost, h


def record_history(hist: History, step: int, E: Tensor, cost: Tensor, h: Tensor, P: Tensor, C_np: np.ndarray, best_cost: float) -> None:
    with torch.no_grad():
        perm = greedy_perm_from_P(P)
        hist.step.append(step)
        hist.energy.append(float(E.detach().cpu()))
        hist.cost.append(float(cost.detach().cpu()))
        hist.h_norm.append(float(torch.linalg.vector_norm(h).detach().cpu()))
        hist.h_vals.append(h.detach().cpu().numpy().copy())
        hist.entropy.append(float(entropy(P).detach().cpu()))
        hist.row_err.append(float(torch.max(torch.abs(P.sum(dim=1) - 1)).detach().cpu()))
        hist.col_err.append(float(torch.max(torch.abs(P.sum(dim=0) - 1)).detach().cpu()))
        hist.sharpness.append(float(sharpness(P).detach().cpu()))
        hist.tour_cost.append(perm_cost(perm, C_np))
        hist.best_perm_cost.append(float(best_cost))


def run_penalty_sinkhorn(C_np: np.ndarray, cfg: RunConfig, name: str = "run") -> Tuple[History, Tensor]:
    set_seed(cfg.seed)
    device = torch.device(cfg.device)
    n = C_np.shape[0]
    C = torch.tensor(C_np, device=device)

    # Small random logits are intentional: exact X=0 can be trace-degenerate.
    X = torch.nn.Parameter(0.01 * torch.randn((n, n), device=device))
    opt = torch.optim.Adam([X], lr=cfg.lr)
    hist = History.empty()

    try:
        best_cost, best_order = brute_force_atsp(C_np)
    except Exception:
        best_cost, best_order = math.nan, tuple(range(n))

    sh_cfg = SinkhornConfig(iters=cfg.sinkhorn_iters)

    for t in range(cfg.steps + 1):
        opt.zero_grad(set_to_none=True)
        rho_t = scheduled_value(t, cfg.steps, cfg.rho, cfg.rho_final, mode="log")
        tau_t = scheduled_value(t, cfg.steps, cfg.tau, cfg.tau_final, mode="log")
        X_eff = apply_no_self_loop_mask(X) if cfg.no_self_loops else X
        P = sinkhorn_from_logits(X_eff, tau=tau_t, cfg=sh_cfg)
        E, cost, h = trace_penalty_energy(P, C, rho_t)
        E.backward()
        opt.step()

        if t % cfg.log_every == 0 or t == cfg.steps:
            record_history(hist, t, E, cost, h, P, C_np, best_cost)

    with torch.no_grad():
        X_eff = apply_no_self_loop_mask(X) if cfg.no_self_loops else X
        P = sinkhorn_from_logits(X_eff, tau=scheduled_value(cfg.steps, cfg.steps, cfg.tau, cfg.tau_final, mode="log"), cfg=sh_cfg)
    return hist, P.detach().cpu()


@dataclass
class PrimalDualConfig:
    steps: int = 3000
    lr_x: float = 0.03
    lr_lam: float = 0.02
    rho: float = 0.0
    tau: float = 1.0
    init_scale: float = 0.15
    sinkhorn_iters: int = 80
    seed: int = 0
    log_every: int = 5
    no_self_loops: bool = True
    device: str = "cpu"


def run_primal_dual(C_np: np.ndarray, cfg: PrimalDualConfig) -> Tuple[History, Tensor, Tensor]:
    set_seed(cfg.seed)
    device = torch.device(cfg.device)
    n = C_np.shape[0]
    C = torch.tensor(C_np, device=device)
    X = torch.nn.Parameter(cfg.init_scale * torch.randn((n, n), device=device))
    lam = torch.zeros(n, dtype=torch.float64, device=device)
    opt = torch.optim.SGD([X], lr=cfg.lr_x)
    hist = History.empty()
    sh_cfg = SinkhornConfig(iters=cfg.sinkhorn_iters)

    try:
        best_cost, _ = brute_force_atsp(C_np)
    except Exception:
        best_cost = math.nan

    for t in range(cfg.steps + 1):
        opt.zero_grad(set_to_none=True)
        X_eff = apply_no_self_loop_mask(X) if cfg.no_self_loops else X
        P = sinkhorn_from_logits(X_eff, tau=cfg.tau, cfg=sh_cfg)
        h = trace_moments(P)
        cost = tsp_cost(P, C)
        # Descent in X, ascent in lambda.
        E = cost + torch.sum(lam.detach() * h) + 0.5 * cfg.rho * torch.sum(h * h)
        E.backward()
        opt.step()
        with torch.no_grad():
            lam += cfg.lr_lam * h.detach()

        if t % cfg.log_every == 0 or t == cfg.steps:
            record_history(hist, t, E.detach(), cost.detach(), h.detach(), P.detach(), C_np, best_cost)

    with torch.no_grad():
        X_eff = apply_no_self_loop_mask(X) if cfg.no_self_loops else X
        P = sinkhorn_from_logits(X_eff, tau=cfg.tau, cfg=sh_cfg)
    return hist, P.detach().cpu(), lam.detach().cpu()


# -----------------------------------------------------------------------------
# Diagnostics
# -----------------------------------------------------------------------------


def trace_gradient_norm_at_epsilon(n: int, eps: float, seed: int, sinkhorn_iters: int, no_self_loops: bool = False) -> float:
    set_seed(seed)
    X = torch.nn.Parameter(eps * torch.randn((n, n)))
    if no_self_loops:
        X_eff = apply_no_self_loop_mask(X)
    else:
        X_eff = X
    P = sinkhorn_from_logits(X_eff, tau=1.0, cfg=SinkhornConfig(iters=sinkhorn_iters))
    h = trace_moments(P)
    # Use only k>=2 if desired? Here use full trace penalty exactly as in optimiser.
    E = 0.5 * torch.sum(h * h)
    E.backward()
    return float(torch.linalg.vector_norm(X.grad).detach().cpu())


def finite_difference_jacobian_h(P: Tensor) -> np.ndarray:
    """Numerical Jacobian of h wrt ambient P entries. For small n diagnostics only."""
    P = P.detach().clone().requires_grad_(True)
    h = trace_moments(P)
    rows = []
    for k in range(h.numel()):
        grad, = torch.autograd.grad(h[k], P, retain_graph=True)
        rows.append(grad.detach().cpu().numpy().reshape(-1))
    return np.stack(rows, axis=0)


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------


def ensure_outdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def plot_history(histories: Dict[str, History], outdir: Path, prefix: str, title: str) -> None:
    ensure_outdir(outdir)

    def arr(xs): return np.asarray(xs, dtype=float)

    metrics = [
        ("h_norm", "Trace residual norm ||h(P)||"),
        ("cost", "Relaxed cost <C,P>"),
        ("entropy", "Entropy H(P)"),
        ("sharpness", "Mean row max"),
        ("tour_cost", "Greedy permutation cost"),
    ]

    for attr, ylabel in metrics:
        plt.figure(figsize=(8, 4.5))
        for label, hist in histories.items():
            y = arr(getattr(hist, attr))
            x = arr(hist.step)
            plt.plot(x, y, label=label)
        # best cost line for tour cost plots.
        if attr == "tour_cost":
            first_hist = next(iter(histories.values()))
            if first_hist.best_perm_cost and not math.isnan(first_hist.best_perm_cost[-1]):
                plt.axhline(first_hist.best_perm_cost[-1], linestyle="--", linewidth=1, label="brute force optimum")
        plt.xlabel("step")
        plt.ylabel(ylabel)
        plt.title(title)
        plt.legend()
        plt.tight_layout()
        plt.savefig(outdir / f"{prefix}_{attr}.png", dpi=160)
        plt.close()

    # Individual h_k residuals for first history, and overlay h_norm already done.
    for label, hist in histories.items():
        H = np.stack(hist.h_vals, axis=0) if hist.h_vals else np.zeros((0, 0))
        if H.size == 0:
            continue
        plt.figure(figsize=(8, 4.5))
        for k in range(H.shape[1]):
            plt.plot(hist.step, H[:, k], label=f"h_{k+1}", linewidth=1)
        plt.xlabel("step")
        plt.ylabel("h_k(P)")
        plt.title(f"{title}: trace residual components ({label})")
        plt.legend(ncol=2, fontsize=8)
        plt.tight_layout()
        safe_label = label.replace(" ", "_").replace("/", "_")
        plt.savefig(outdir / f"{prefix}_{safe_label}_h_components.png", dpi=160)
        plt.close()


def plot_points_and_perm(pts: np.ndarray, P: Tensor, C_np: np.ndarray, outdir: Path, prefix: str, title: str, method: str = "assignment") -> None:
    ensure_outdir(outdir)
    perm = projection_perm_from_P(P, method=method)
    cycles = count_cycles(perm)
    plt.figure(figsize=(5, 5))
    plt.scatter(pts[:, 0], pts[:, 1])
    for i, (x, y) in enumerate(pts):
        plt.text(x, y, f" {i}", fontsize=9)
    for i, j in enumerate(perm):
        x0, y0 = pts[i]
        x1, y1 = pts[j]
        plt.arrow(x0, y0, 0.85 * (x1 - x0), 0.85 * (y1 - y0), head_width=0.06, length_includes_head=True, alpha=0.65)
    plt.axis("equal")
    plt.title(f"{title}\ncycles={list(map(len, cycles))}, {method} cost={perm_cost(perm, C_np):.3f}")
    plt.tight_layout()
    plt.savefig(outdir / f"{prefix}_{method}_perm.png", dpi=160)
    plt.close()


def plot_gradient_degeneracy(eps_values: List[float], norms: Dict[str, List[float]], outdir: Path) -> None:
    ensure_outdir(outdir)
    plt.figure(figsize=(7, 4.5))
    for label, vals in norms.items():
        plt.loglog(eps_values, vals, marker="o", label=label)
    plt.xlabel("logit perturbation scale epsilon")
    plt.ylabel("||grad_X 0.5||h(P)||^2||")
    plt.title("Uniform initialisation degeneracy diagnostic")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "uniform_degeneracy_grad_norm.png", dpi=160)
    plt.close()




def trace_jacobian_wrt_logits(X: Tensor, tau: float, sinkhorn_iters: int, no_self_loops: bool) -> Tuple[np.ndarray, np.ndarray]:
    """Return h(X) and Jacobian dh/dvec(X) at a logit point.

    This is a local small-signal object. It intentionally keeps the full logit
    coordinate system, including Sinkhorn gauge null directions.
    """
    X = X.detach().clone().requires_grad_(True)
    X_eff = apply_no_self_loop_mask(X) if no_self_loops else X
    P = sinkhorn_from_logits(X_eff, tau=tau, cfg=SinkhornConfig(iters=sinkhorn_iters))
    h = trace_moments(P)
    rows = []
    for k in range(h.numel()):
        grad, = torch.autograd.grad(h[k], X, retain_graph=True)
        rows.append(grad.detach().cpu().numpy().reshape(-1))
    A = np.stack(rows, axis=0)
    return h.detach().cpu().numpy(), A


def simulate_linear_primal_dual(A: np.ndarray, rho: float, eta: float, dt: float, steps: int, amp: float = 1.0) -> Dict[str, np.ndarray]:
    """Small-signal model around a frozen Sinkhorn point.

    z_dot      = - A^T lambda - rho A^T A z
    lambda_dot = eta A z

    Here y = A z is the linearised trace residual. In a singular direction of
    A with singular value sigma, the scalar mode obeys approximately

        x_ddot + rho*sigma^2*x_dot + eta*sigma^2*x = 0.

    So pure rho=0 is an undamped oscillator and rho>0 damps the trace-visible
    modes. We initialise exactly in the top singular direction so the effect is
    not hidden in gauge/null directions.
    """
    U, S, Vt = np.linalg.svd(A, full_matrices=False)
    if S.size == 0 or S[0] <= 0:
        raise ValueError("Trace Jacobian has no positive singular value; cannot build ringing demo.")

    u0 = U[:, 0].copy()
    z = amp * Vt[0].copy()
    lam = np.zeros(A.shape[0])
    ys = []
    norms = []
    mode_y = []
    lams = []
    for _ in range(steps + 1):
        y = A @ z
        ys.append(y.copy())
        norms.append(float(np.linalg.norm(y)))
        mode_y.append(float(u0 @ y))
        lams.append(lam.copy())

        z_dot = -(A.T @ lam) - rho * (A.T @ y)
        lam_dot = eta * y

        # Semi-implicit Euler is stable for this oscillator at the chosen dt.
        lam = lam + dt * lam_dot
        z = z + dt * z_dot

    return {
        "y": np.asarray(ys),
        "norm": np.asarray(norms),
        "mode_y": np.asarray(mode_y),
        "lambda": np.asarray(lams),
        "singular_values": S,
        "eta": np.asarray([eta]),
        "rho": np.asarray([rho]),
        "dt": np.asarray([dt]),
    }


def plot_linear_ringing(linear_runs: Dict[str, Dict[str, np.ndarray]], outdir: Path, prefix: str, title: str, dt: float) -> None:
    ensure_outdir(outdir)

    # Norm is useful, but it hides sign changes; the top singular-mode coordinate
    # below is the clearest ringing plot.
    plt.figure(figsize=(8, 4.5))
    for label, run in linear_runs.items():
        x = np.arange(run["norm"].shape[0]) * dt
        plt.plot(x, run["norm"], label=label)
    plt.xlabel("linearised time")
    plt.ylabel("||A z||  (linearised trace residual norm)")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / f"{prefix}_linear_h_norm.png", dpi=160)
    plt.close()

    plt.figure(figsize=(8, 4.5))
    for label, run in linear_runs.items():
        x = np.arange(run["mode_y"].shape[0]) * dt
        plt.plot(x, run["mode_y"], label=label)
    plt.axhline(0.0, linewidth=0.8)
    plt.xlabel("linearised time")
    plt.ylabel("top singular trace mode  u₁ᵀ A z")
    plt.title(f"{title}: signed top mode")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / f"{prefix}_linear_top_mode.png", dpi=160)
    plt.close()

    for label, run in linear_runs.items():
        Y = run["y"]
        x = np.arange(Y.shape[0]) * dt
        plt.figure(figsize=(8, 4.5))
        for k in range(Y.shape[1]):
            plt.plot(x, Y[:, k], label=f"h_{k+1}", linewidth=1)
        plt.xlabel("linearised time")
        plt.ylabel("linearised h_k")
        plt.title(f"{title}: components ({label})")
        plt.legend(ncol=2, fontsize=8)
        plt.tight_layout()
        safe_label = label.replace(" ", "_").replace("/", "_").replace("=", "")
        plt.savefig(outdir / f"{prefix}_{safe_label}_linear_h_components.png", dpi=160)
        plt.close()

    # Singular spectrum of A. This tells us how many trace-visible directions exist.
    first = next(iter(linear_runs.values()))
    svals = first["singular_values"]
    plt.figure(figsize=(6, 4))
    plt.semilogy(np.arange(1, len(svals) + 1), svals, marker="o")
    plt.xlabel("index")
    plt.ylabel("singular value of A = Dh(X)")
    plt.title(f"{title}: trace-Jacobian spectrum")
    plt.tight_layout()
    plt.savefig(outdir / f"{prefix}_trace_jacobian_spectrum.png", dpi=160)
    plt.close()

    # Save the chosen oscillator scaling for reproducibility.
    with open(outdir / f"{prefix}_linear_scaling.txt", "w") as f:
        for label, run in linear_runs.items():
            f.write(
                f"{label}: sigma1={run['singular_values'][0]:.8g}, "
                f"eta={float(run['eta'][0]):.8g}, rho={float(run['rho'][0]):.8g}, "
                f"dt={float(run['dt'][0]):.8g}\n"
            )


def plot_rho_pressure(rhos: List[float], summaries: Dict[str, List[float]], cycle_strings: List[str], outdir: Path, prefix: str, title: str) -> None:
    ensure_outdir(outdir)
    x = np.asarray(rhos)
    for key, ylabel in [
        ("h_norm", "final ||h(P)||"),
        ("relaxed_cost", "final relaxed cost <C,P>"),
        ("greedy_cost", "final greedy permutation cost"),
        ("sharpness", "final mean row max"),
        ("entropy", "final entropy H(P)"),
    ]:
        plt.figure(figsize=(7, 4.5))
        plt.semilogx(x, summaries[key], marker="o")
        plt.xlabel("fixed trace penalty rho")
        plt.ylabel(ylabel)
        plt.title(title)
        plt.tight_layout()
        plt.savefig(outdir / f"{prefix}_rho_pressure_{key}.png", dpi=160)
        plt.close()

    # Categorical-ish plot of whether the greedy assignment is one cycle.
    is_single = [1.0 if cs == "[8]" or cs == f"[{len(cycle_strings)}]" else 0.0 for cs in cycle_strings]
    # Better: parse the actual text safely enough for our generated strings.
    is_single = []
    for cs in cycle_strings:
        nums = [int(tok) for tok in cs.strip("[]").replace(",", " ").split() if tok]
        is_single.append(1.0 if len(nums) == 1 else 0.0)
    plt.figure(figsize=(7, 4.5))
    plt.semilogx(x, is_single, marker="o")
    plt.ylim(-0.1, 1.1)
    plt.xlabel("fixed trace penalty rho")
    plt.ylabel("greedy assignment is single cycle")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outdir / f"{prefix}_rho_pressure_single_cycle.png", dpi=160)
    plt.close()

    with open(outdir / f"{prefix}_rho_pressure_summary.txt", "w") as f:
        f.write("rho\th_norm\trelaxed_cost\tgreedy_cost\tsharpness\tentropy\tcycles\n")
        for i, rho in enumerate(rhos):
            f.write(
                f"{rho}\t{summaries['h_norm'][i]}\t{summaries['relaxed_cost'][i]}\t"
                f"{summaries['greedy_cost'][i]}\t{summaries['sharpness'][i]}\t"
                f"{summaries['entropy'][i]}\t{cycle_strings[i]}\n"
            )



def write_summary_csv(rows: List[Dict[str, object]], path: Path) -> None:
    if not rows:
        return
    keys = sorted({k for row in rows for k in row.keys()})
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def median_iqr(values: List[float]) -> Tuple[float, float, float]:
    arr = np.asarray(values, dtype=float)
    return float(np.median(arr)), float(np.quantile(arr, 0.25)), float(np.quantile(arr, 0.75))


def plot_rho_pressure_aggregate(rows: List[Dict[str, object]], outdir: Path, prefix: str, title: str) -> None:
    ensure_outdir(outdir)
    rhos = sorted({float(r["rho"]) for r in rows})
    metrics = [
        ("h_norm", "final ||h(P)||"),
        ("assignment_cost", "assignment-projected cost"),
        ("relaxed_cost", "relaxed cost <C,P>"),
        ("entropy", "entropy H(P)"),
        ("sharpness", "mean row max"),
    ]
    for key, ylabel in metrics:
        med, q1, q3 = [], [], []
        for rho in rhos:
            vals = [float(r[key]) for r in rows if float(r["rho"]) == rho]
            m, lo, hi = median_iqr(vals)
            med.append(m); q1.append(lo); q3.append(hi)
        plt.figure(figsize=(7, 4.5))
        plt.semilogx(rhos, med, marker="o")
        plt.fill_between(rhos, q1, q3, alpha=0.2)
        plt.xlabel("fixed trace penalty rho")
        plt.ylabel(ylabel)
        plt.title(title)
        plt.tight_layout()
        plt.savefig(outdir / f"{prefix}_rho_pressure_{key}_aggregate.png", dpi=160)
        plt.close()

    probs = []
    for rho in rhos:
        vals = [1.0 if bool(r["assignment_single_cycle"]) else 0.0 for r in rows if float(r["rho"]) == rho]
        probs.append(float(np.mean(vals)))
    plt.figure(figsize=(7, 4.5))
    plt.semilogx(rhos, probs, marker="o")
    plt.ylim(-0.05, 1.05)
    plt.xlabel("fixed trace penalty rho")
    plt.ylabel("P(assignment projection is one cycle)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outdir / f"{prefix}_rho_pressure_single_cycle_probability.png", dpi=160)
    plt.close()


def plot_phase_heatmaps(rows: List[Dict[str, object]], outdir: Path, prefix: str, title: str) -> None:
    ensure_outdir(outdir)
    rhos = sorted({float(r["rho"]) for r in rows})
    taus = sorted({float(r["tau"]) for r in rows})
    metrics = [
        ("h_norm", "final ||h(P)||"),
        ("assignment_cost", "assignment-projected cost"),
        ("assignment_single_cycle", "single-cycle indicator"),
        ("entropy", "entropy H(P)"),
        ("sharpness", "mean row max"),
    ]
    for key, label in metrics:
        Z = np.full((len(taus), len(rhos)), np.nan)
        for i, tau in enumerate(taus):
            for j, rho in enumerate(rhos):
                vals = [float(r[key]) for r in rows if float(r["rho"]) == rho and float(r["tau"]) == tau]
                if vals:
                    Z[i, j] = float(np.median(vals))
        plt.figure(figsize=(8, 4.8))
        im = plt.imshow(Z, aspect="auto", origin="lower")
        plt.colorbar(im, label=label)
        plt.xticks(range(len(rhos)), [f"{r:g}" for r in rhos], rotation=45)
        plt.yticks(range(len(taus)), [f"{t:g}" for t in taus])
        plt.xlabel("rho")
        plt.ylabel("tau")
        plt.title(f"{title}: {label}")
        plt.tight_layout()
        plt.savefig(outdir / f"{prefix}_phase_{key}.png", dpi=160)
        plt.close()

# -----------------------------------------------------------------------------
# Experiment driver
# -----------------------------------------------------------------------------


def experiment_circle(outdir: Path, n: int, seed: int) -> None:
    pts, C = circle_instance(n)
    best, order = brute_force_atsp(C)
    print(f"[circle] n={n}, brute force best={best:.6f}, order={order}")

    runs: Dict[str, History] = {}
    cfg1 = RunConfig(steps=1800, lr=0.04, rho=1.0, rho_final=150.0, tau=1.5, tau_final=0.18, seed=seed, log_every=5)
    h1, P1 = run_penalty_sinkhorn(C, cfg1, name="scheduled")
    runs["rho up / tau down"] = h1

    cfg2 = RunConfig(steps=1800, lr=0.04, rho=30.0, rho_final=None, tau=0.7, tau_final=None, seed=seed, log_every=5)
    h2, P2 = run_penalty_sinkhorn(C, cfg2, name="fixed")
    runs["fixed rho/tau"] = h2

    plot_history(runs, outdir, "circle", f"Circle TSP sanity check (n={n})")
    plot_points_and_perm(pts, P1, C, outdir, "circle_scheduled", "Circle scheduled final greedy assignment")
    plot_points_and_perm(pts, P2, C, outdir, "circle_fixed", "Circle fixed final greedy assignment")


def experiment_ringing(outdir: Path, n: int, seed: int) -> None:
    pts, C = zero_cost_instance(n)
    # With diagonal masked in logits, diag cost is mostly irrelevant. Tiny random cost prevents exact symmetry.
    rng = np.random.default_rng(seed)
    C = 0.001 * rng.normal(size=(n, n))
    np.fill_diagonal(C, 1e6)

    # Nonlinear Sinkhorn primal-dual run. This often shows delayed collapse/snap rather
    # than textbook sinusoidal ringing, which is still a useful instability mode.
    runs: Dict[str, History] = {}
    cfg_pd = PrimalDualConfig(steps=2500, lr_x=0.08, lr_lam=0.08, rho=0.0, tau=1.0, init_scale=0.15, seed=seed, log_every=2)
    h_pd, _, _ = run_primal_dual(C, cfg_pd)
    runs["pure primal-dual rho=0"] = h_pd

    cfg_aug = PrimalDualConfig(steps=2500, lr_x=0.08, lr_lam=0.08, rho=25.0, tau=1.0, init_scale=0.15, seed=seed, log_every=2)
    h_aug, _, _ = run_primal_dual(C, cfg_aug)
    runs["augmented rho=25"] = h_aug

    # A softer small-step run tries to keep the dynamics in the interior/linear regime.
    cfg_soft = PrimalDualConfig(steps=6000, lr_x=0.008, lr_lam=0.08, rho=0.0, tau=4.0, init_scale=0.35, seed=seed + 13, log_every=5)
    h_soft, _, _ = run_primal_dual(C, cfg_soft)
    runs["pure primal-dual soft tau=4"] = h_soft

    plot_history(runs, outdir, "ringing", f"Constraint-only ringing/snap demo (n={n})")

    # Linearised small-signal model around a frozen nonuniform Sinkhorn point. This
    # isolates the Jacobian story: pure primal-dual oscillates; rho A^T A damps.
    set_seed(seed + 101)
    X0 = 0.35 * torch.randn((n, n), dtype=torch.float64)
    _, A = trace_jacobian_wrt_logits(X0, tau=3.0, sinkhorn_iters=120, no_self_loops=True)
    # Scale eta and rho from the top singular value so the oscillator is visible
    # on a human-sized time axis. For a singular mode sigma:
    #   x_ddot + rho*sigma^2*x_dot + eta*sigma^2*x = 0.
    # Choose eta for a target period and rho for a target damping ratio.
    svals = np.linalg.svd(A, compute_uv=False)
    sigma1 = float(svals[0])
    target_period = 12.0
    damping_ratio = 0.22
    eta = (2.0 * math.pi / (target_period * sigma1)) ** 2
    rho_damped = 2.0 * damping_ratio * math.sqrt(eta) / sigma1
    dt = target_period / 600.0
    steps = int(5.0 * target_period / dt)
    lin_pure = simulate_linear_primal_dual(A, rho=0.0, eta=eta, dt=dt, steps=steps, amp=1.0 / sigma1)
    lin_damped = simulate_linear_primal_dual(A, rho=rho_damped, eta=eta, dt=dt, steps=steps, amp=1.0 / sigma1)
    plot_linear_ringing(
        {"linear pure rho=0": lin_pure, f"linear damped zeta={damping_ratio:g}": lin_damped},
        outdir,
        "ringing",
        f"Linearised trace-mode ringing around Sinkhorn point (n={n})",
        dt,
    )


def experiment_uniform_degeneracy(outdir: Path, n: int, seed: int) -> None:
    eps_values = [0.0, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1]
    # loglog cannot plot zero x; use tiny placeholder for display but compute exact zero too.
    plot_eps = [1e-10 if e == 0.0 else e for e in eps_values]
    vals_no_diag = []
    vals_with_diag = []
    for e in eps_values:
        vals_no_diag.append(trace_gradient_norm_at_epsilon(n, e, seed=seed, sinkhorn_iters=100, no_self_loops=False))
        vals_with_diag.append(trace_gradient_norm_at_epsilon(n, e, seed=seed, sinkhorn_iters=100, no_self_loops=True))
    print("[uniform] eps vs grad norm, no mask:")
    for e, v in zip(eps_values, vals_no_diag):
        print(f"  eps={e:.1e}: {v:.3e}")
    print("[uniform] eps vs grad norm, diagonal masked:")
    for e, v in zip(eps_values, vals_with_diag):
        print(f"  eps={e:.1e}: {v:.3e}")
    plot_gradient_degeneracy(plot_eps, {"full matrix": vals_no_diag, "diag masked": vals_with_diag}, outdir)

    # Also save singular values of ambient trace Jacobian at uniform for reference.
    P0 = torch.ones((n, n), dtype=torch.float64) / n
    J = finite_difference_jacobian_h(P0)
    s = np.linalg.svd(J, compute_uv=False)
    np.savetxt(outdir / "uniform_ambient_trace_jacobian_singular_values.txt", s)


def experiment_two_cluster(outdir: Path, n: int, seed: int) -> None:
    pts, C = two_cluster_instance(n, sep=6.0, sigma=0.35, seed=seed)
    best, order = brute_force_atsp(C)
    print(f"[two-cluster] n={n}, brute force best={best:.6f}, order={order}")

    runs: Dict[str, History] = {}
    cfg_small = RunConfig(steps=2400, lr=0.035, rho=2.0, tau=0.8, seed=seed, log_every=5)
    h_small, P_small = run_penalty_sinkhorn(C, cfg_small)
    runs["fixed small rho"] = h_small

    cfg_large = RunConfig(steps=2400, lr=0.025, rho=120.0, tau=0.8, seed=seed, log_every=5)
    h_large, P_large = run_penalty_sinkhorn(C, cfg_large)
    runs["fixed large rho"] = h_large

    cfg_sched = RunConfig(steps=2400, lr=0.035, rho=1.0, rho_final=180.0, tau=1.7, tau_final=0.18, seed=seed, log_every=5)
    h_sched, P_sched = run_penalty_sinkhorn(C, cfg_sched)
    runs["rho up / tau down"] = h_sched

    plot_history(runs, outdir, "two_cluster", f"Two-cluster TSP stress test (n={n})")
    plot_points_and_perm(pts, P_small, C, outdir, "two_cluster_small", "Two-cluster fixed small rho")
    plot_points_and_perm(pts, P_large, C, outdir, "two_cluster_large", "Two-cluster fixed large rho")
    plot_points_and_perm(pts, P_sched, C, outdir, "two_cluster_scheduled", "Two-cluster scheduled")



def experiment_rho_pressure(outdir: Path, n: int, seeds: List[int]) -> None:
    """Sweep fixed rho on the two-cluster instance to expose the validity/pressure threshold.

    Uses exact assignment projection for the readout, not the greedy projection.
    Multiple seeds are aggregated into median/IQR plots and a single-cycle probability.
    """
    rhos = [0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 400.0]
    rows: List[Dict[str, object]] = []
    for seed in seeds:
        pts, C = two_cluster_instance(n, sep=6.0, sigma=0.35, seed=seed)
        best, order = brute_force_atsp(C)
        print(f"[rho-pressure] seed={seed} two-cluster brute force best={best:.6f}, order={order}")
        for rho in rhos:
            cfg = RunConfig(steps=300, lr=0.04, rho=rho, tau=0.8, seed=seed, log_every=10, sinkhorn_iters=50)
            hist, P = run_penalty_sinkhorn(C, cfg)
            greedy = greedy_perm_from_P(P)
            assign = assignment_perm_from_P(P)
            assign_cycles = count_cycles(assign)
            row = {
                "experiment": "rho_pressure",
                "seed": seed,
                "n": n,
                "rho": rho,
                "tau": 0.8,
                "h_norm": hist.h_norm[-1],
                "relaxed_cost": hist.cost[-1],
                "greedy_cost": perm_cost(greedy, C),
                "assignment_cost": perm_cost(assign, C),
                "opt_cost": best,
                "assignment_gap": (perm_cost(assign, C) - best) / best,
                "sharpness": hist.sharpness[-1],
                "entropy": hist.entropy[-1],
                "greedy_cycles": str(list(map(len, count_cycles(greedy)))),
                "assignment_cycles": str(list(map(len, assign_cycles))),
                "assignment_single_cycle": len(assign_cycles) == 1,
            }
            rows.append(row)
            print(
                f"  seed={seed:3d} rho={rho:7.2f} h={row['h_norm']:9.3e} "
                f"assign={row['assignment_cost']:8.3f} cycles={row['assignment_cycles']}"
            )
    write_summary_csv(rows, outdir / "two_cluster_rho_pressure_summary.csv")
    plot_rho_pressure_aggregate(rows, outdir, "two_cluster", f"Two-cluster rho-pressure sweep (n={n}, seeds={len(seeds)})")


def experiment_phase_diagram(outdir: Path, n: int, seeds: List[int]) -> None:
    """A small rho/tau phase diagram for the two-cluster instance."""
    rhos = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
    taus = [0.35, 0.5, 0.8, 1.2, 1.8, 3.0]
    rows: List[Dict[str, object]] = []
    # Keep this cheap: use a few seeds, but shorter runs than the main sweep.
    for seed in seeds:
        pts, C = two_cluster_instance(n, sep=6.0, sigma=0.35, seed=seed)
        best, order = brute_force_atsp(C)
        print(f"[phase] seed={seed} brute force best={best:.6f}")
        for tau in taus:
            for rho in rhos:
                cfg = RunConfig(steps=220, lr=0.04, rho=rho, tau=tau, seed=seed, log_every=30, sinkhorn_iters=50)
                hist, P = run_penalty_sinkhorn(C, cfg)
                assign = assignment_perm_from_P(P)
                assign_cycles = count_cycles(assign)
                rows.append({
                    "experiment": "phase_diagram",
                    "seed": seed,
                    "n": n,
                    "rho": rho,
                    "tau": tau,
                    "h_norm": hist.h_norm[-1],
                    "relaxed_cost": hist.cost[-1],
                    "assignment_cost": perm_cost(assign, C),
                    "opt_cost": best,
                    "assignment_gap": (perm_cost(assign, C) - best) / best,
                    "sharpness": hist.sharpness[-1],
                    "entropy": hist.entropy[-1],
                    "assignment_cycles": str(list(map(len, assign_cycles))),
                    "assignment_single_cycle": len(assign_cycles) == 1,
                })
    write_summary_csv(rows, outdir / "two_cluster_phase_diagram_summary.csv")
    plot_phase_heatmaps(rows, outdir, "two_cluster", f"Two-cluster rho/tau phase diagram (n={n}, seeds={len(seeds)})")

def main() -> None:
    parser = argparse.ArgumentParser(description="Cycle-moment Sinkhorn flow experiments")
    parser.add_argument("--outdir", type=Path, default=Path("cycle_moment_outputs"))
    parser.add_argument("--n", type=int, default=8, help="small n; brute force scales factorially")
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--seeds", type=str, default=None, help="comma-separated seeds for aggregate sweeps; default uses --seed")
    parser.add_argument("--skip-circle", action="store_true")
    parser.add_argument("--skip-ringing", action="store_true")
    parser.add_argument("--skip-uniform", action="store_true")
    parser.add_argument("--skip-two-cluster", action="store_true")
    parser.add_argument("--skip-rho-pressure", action="store_true")
    parser.add_argument("--run-phase-diagram", action="store_true", help="run the optional rho/tau grid; slower")
    args = parser.parse_args()

    ensure_outdir(args.outdir)
    seeds = [args.seed] if args.seeds is None else [int(x) for x in args.seeds.split(",") if x.strip()]
    print(f"Writing outputs to {args.outdir.resolve()}")
    print("Tip: n=8 is intentionally small so brute force is feasible.")
    print(f"Aggregate seeds: {seeds}")

    if not args.skip_circle:
        experiment_circle(args.outdir, args.n, args.seed)
    if not args.skip_ringing:
        experiment_ringing(args.outdir, args.n, args.seed)
    if not args.skip_uniform:
        experiment_uniform_degeneracy(args.outdir, args.n, args.seed)
    if not args.skip_two_cluster:
        experiment_two_cluster(args.outdir, args.n, args.seed)
    if not args.skip_rho_pressure:
        experiment_rho_pressure(args.outdir, args.n, seeds)
    if args.run_phase_diagram:
        experiment_phase_diagram(args.outdir, args.n, seeds[: min(len(seeds), 3)])

    print("Done. Generated PNG plots and diagnostics in the output directory.")


if __name__ == "__main__":
    main()
