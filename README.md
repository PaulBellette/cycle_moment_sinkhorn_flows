# Cycle-Moment Sinkhorn Flows

A small research prototype exploring a Brockett-like, trace-constrained Sinkhorn/KL flow on the Birkhoff polytope. The motivating example is the traveling salesperson problem (TSP), but the broader object is a dynamical system on doubly stochastic matrices with structural trace-moment constraints.

The project is **not** intended as a competitive TSP solver. The aim is to understand the geometry and stability of trace-constrained Sinkhorn routing: when trace penalties enforce Hamiltonicity, when they leave cost/optimality unresolved, how symmetric initialisation can kill first-order trace signal, and why augmented primal-dual dynamics damps trace-visible modes.

## Core idea

A directed tour can be represented by a permutation matrix. A Hamiltonian tour is a permutation consisting of a single \(n\)-cycle. Over the Birkhoff polytope,

\[
\mathcal B_n=\{P\ge 0:P\mathbf 1=\mathbf 1,\;P^\top\mathbf 1=\mathbf 1\},
\]

trace powers detect weighted closed walks:

\[
h_k(P)=\operatorname{tr}(P^k)-b_k.
\]

For Hamiltonian-cycle constraints we use

\[
b_k=0\quad (1\le k<n),\qquad b_n=n.
\]

The flow uses KL/Sinkhorn geometry rather than Euclidean projection. Given an energy \(E(P)\) with Euclidean gradient \(G=\nabla_P E(P)\), the Birkhoff-preserving KL-gradient flow has coordinate form

\[
\dot P_{ij}=-P_{ij}(G_{ij}-a_i-b_j),
\]

where \(a,b\) are row/column potentials chosen to preserve the marginal constraints. The multiplicative form preserves positivity and gives a transport-polytope analogue of Brockett-style constrained matrix flows.

## Repository layout

```text
.
├── cycle_moment_sinkhorn_experiments.py   # runnable uv script
├── cycle_moment_outputs/                  # canonical numerical outputs
├── notes/
│   └── cycle_moment_sinkhorn_flows_note.md # working research note
└── paper/
    ├── main.tex                            # draft LaTeX paper
    └── references.bib                      # draft bibliography
```

## Running the experiments

The script is designed to run directly with `uv`:

```bash
uv run cycle_moment_sinkhorn_experiments.py --seeds 0,1,2,3,4
```

The optional phase-diagram sweep is slower:

```bash
uv run cycle_moment_sinkhorn_experiments.py \
  --skip-circle \
  --skip-ringing \
  --skip-uniform \
  --skip-two-cluster \
  --seeds 0,1,2 \
  --run-phase-diagram
```

Outputs are written to `cycle_moment_outputs/`.

## Key numerical results

### 1. Circle sanity check

On points arranged on a circle, a continuation schedule with increasing trace pressure and decreasing Sinkhorn temperature finds the expected geometric tour. A fixed penalty can still produce a valid single cycle, but with bad crossings/chords.

- [circle cost](cycle_moment_outputs/circle_cost.png)
- [circle trace residual norm](cycle_moment_outputs/circle_h_norm.png)
- [scheduled final assignment](cycle_moment_outputs/circle_scheduled_assignment_perm.png)
- [fixed final assignment](cycle_moment_outputs/circle_fixed_assignment_perm.png)

### 2. Two-cluster subtour pressure

A two-cluster instance exposes the subtour failure mode. Low trace pressure favours cheap within-cluster cycle covers such as \([4,4]\); sufficiently large trace pressure forces a single 8-cycle, often at higher assignment cost.

- [scheduled assignment: two subtours](cycle_moment_outputs/two_cluster_scheduled_assignment_perm.png)
- [small-rho assignment: two subtours](cycle_moment_outputs/two_cluster_small_assignment_perm.png)
- [large-rho assignment: single cycle](cycle_moment_outputs/two_cluster_large_assignment_perm.png)
- [two-cluster trace residual norm](cycle_moment_outputs/two_cluster_h_norm.png)

### 3. Rho-pressure transition

Sweeping fixed \(\rho\) reveals a sharp transition from subtour covers to single-cycle assignments.

- [single-cycle probability](cycle_moment_outputs/two_cluster_rho_pressure_single_cycle_probability.png)
- [final trace residual norm](cycle_moment_outputs/two_cluster_rho_pressure_h_norm_aggregate.png)
- [assignment-projected cost](cycle_moment_outputs/two_cluster_rho_pressure_assignment_cost_aggregate.png)
- [entropy](cycle_moment_outputs/two_cluster_rho_pressure_entropy_aggregate.png)
- [sharpness](cycle_moment_outputs/two_cluster_rho_pressure_sharpness_aggregate.png)
- [summary CSV](cycle_moment_outputs/two_cluster_rho_pressure_summary.csv)

### 4. Rho/tau phase diagram

The phase diagram shows that \(\rho\) mostly controls the subtour-to-Hamiltonian transition, while \(\tau\) controls entropy and assignment sharpness.

- [single-cycle phase diagram](cycle_moment_outputs/two_cluster_phase_assignment_single_cycle.png)
- [trace residual phase diagram](cycle_moment_outputs/two_cluster_phase_h_norm.png)
- [assignment cost phase diagram](cycle_moment_outputs/two_cluster_phase_assignment_cost.png)
- [entropy phase diagram](cycle_moment_outputs/two_cluster_phase_entropy.png)
- [sharpness phase diagram](cycle_moment_outputs/two_cluster_phase_sharpness.png)
- [phase diagram CSV](cycle_moment_outputs/two_cluster_phase_diagram_summary.csv)

### 5. Uniform initialisation degeneracy

At the uniform doubly stochastic point, all trace powers with \(k\ge2\) are first-order invisible on the Birkhoff tangent space. With the diagonal/no-self-loop direction removed, the trace penalty is first-order silent at exact uniform initialisation.

- [uniform initialisation degeneracy diagnostic](cycle_moment_outputs/uniform_degeneracy_grad_norm.png)
- [ambient trace-Jacobian singular values](cycle_moment_outputs/uniform_ambient_trace_jacobian_singular_values.txt)

### 6. Linearised ringing and augmented damping

The local primal-dual trace-mode model gives an undamped oscillator when \(\rho=0\) and a damped oscillator when \(\rho>0\):

\[
\ddot q+\rho\sigma^2\dot q+\eta\sigma^2q=0.
\]

- [linearised top-mode ringing and damping](cycle_moment_outputs/ringing_linear_top_mode.png)
- [trace-Jacobian singular spectrum](cycle_moment_outputs/ringing_trace_jacobian_spectrum.png)
- [linear scaling parameters](cycle_moment_outputs/ringing_linear_scaling.txt)

### 7. Nonlinear primal-dual snapping

In the full nonlinear Sinkhorn/logit parameterisation, undamped primal-dual control can appear as delayed high-residual plateaus followed by sharp support changes rather than clean sinusoidal ringing.

- [nonlinear trace residual norm](cycle_moment_outputs/ringing_h_norm.png)
- [pure primal-dual trace components](cycle_moment_outputs/ringing_pure_primal-dual_rho=0_h_components.png)
- [soft-temperature primal-dual trace components](cycle_moment_outputs/ringing_pure_primal-dual_soft_tau=4_h_components.png)
- [augmented trace components](cycle_moment_outputs/ringing_augmented_rho=25_h_components.png)
- [nonlinear entropy](cycle_moment_outputs/ringing_entropy.png)
- [nonlinear sharpness](cycle_moment_outputs/ringing_sharpness.png)

## Draft paper

A LaTeX draft is in [`paper/main.tex`](paper/main.tex). Build it with something like:

```bash
cd paper
pdflatex main
bibtex main
pdflatex main
pdflatex main
```

The draft uses figures directly from `../cycle_moment_outputs/` and citations from `paper/references.bib`.

## Working note

The more informal working derivation is in [`notes/cycle_moment_sinkhorn_flows_note.md`](notes/cycle_moment_sinkhorn_flows_note.md). It contains longer proof sketches and theory-roadmap commentary.

## Current status

The main supported claims are diagnostic rather than solver-performance claims:

1. Trace moments characterise Hamiltonian-cycle permutations over the Birkhoff polytope.
2. KL/Sinkhorn geometry gives a natural positivity-preserving projected flow.
3. Pure penalty flow has a Lyapunov energy.
4. Pure primal-dual trace modes ring in the local linear model.
5. Augmentation damps trace-visible modes.
6. Uniform no-self-loop initialisation is first-order degenerate.
7. Trace pressure creates a sharp transition between subtour covers and single-cycle assignments.

## Reference trail

The draft bibliography currently includes references for Brockett flow, trace acyclicity, trace-guided ATSP relaxations, Sinkhorn/mirror-flow methods, assignment-polytope geometry, and reduced costs/linear programming. The metadata should still be checked before any public preprint-style release.
