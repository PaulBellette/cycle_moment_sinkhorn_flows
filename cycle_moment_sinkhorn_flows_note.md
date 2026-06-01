# Cycle-Moment Sinkhorn Flows: A Brockett-like View of Trace-Constrained Routing on the Birkhoff Polytope

## Abstract

We sketch a geometric flow on the positive interior of the Birkhoff polytope motivated by the traveling salesperson problem, Brockett's double-bracket flow, and modern Sinkhorn/mirror-descent methods. The central object is a KL-gradient flow over doubly stochastic matrices augmented by trace-power constraints of the form

\[
h_k(P)=\operatorname{tr}(P^k)-b_k.
\]

For TSP-like routing, these trace moments distinguish single Hamiltonian cycles from subtour decompositions. The resulting system is not proposed primarily as a competitive TSP solver, but as a mathematical object: a transport-polytope analogue of Brockett-style constrained matrix flows. We discuss its constrained geometry, primal-dual and augmented variants, local linearisation, Lyapunov structure, damping of trace-violating modes, degeneracy near the uniform doubly stochastic matrix, and reduced-cost stability near permutation vertices.

## 1. Motivation

Brockett's double-bracket flow gives a beautiful dynamical-systems view of sorting, diagonalisation, and constrained matrix evolution. The classical setting is typically an isospectral orbit or orthogonal/conjugacy manifold, where commutator structure preserves constraints while a potential function drives the system toward an ordered state.

Here we consider an analogous question on the Birkhoff polytope:

\[
\mathcal B_n=\{P\in \mathbb R_+^{n\times n}:P\mathbf 1=\mathbf 1,\;P^\top\mathbf 1=\mathbf 1\}.
\]

The vertices of \(\mathcal B_n\) are permutation matrices. For directed TSP, a tour can be represented by a permutation matrix corresponding to a single \(n\)-cycle. The usual assignment relaxation allows all permutation matrices, including disjoint subtour decompositions. Trace powers provide a compact way to detect closed walks and subtours.

The motivating question is:

> Can we define a Brockett-like flow on the Birkhoff polytope whose geometry is Sinkhorn/KL rather than orthogonal/commutator, and whose moment constraints favour Hamiltonian cycles?

The practical version of this idea is a differentiable Sinkhorn relaxation for TSP. The more interesting mathematical version is a constrained mirror-gradient flow with trace moments and a stability theory that may apply more generally to differentiable routing and assignment problems.

## 2. The Birkhoff and Trace-Moment Setup

Let

\[
P\in \mathcal B_n^\circ
\]

where \(\mathcal B_n^\circ\) denotes the strictly positive interior of the Birkhoff polytope:

\[
\mathcal B_n^\circ=\{P>0:P\mathbf 1=\mathbf 1,\;P^\top\mathbf 1=\mathbf 1\}.
\]

Let \(C\in \mathbb R^{n\times n}\) be an edge-cost matrix. The linear assignment/TSP-like cost is

\[
\langle C,P\rangle=\operatorname{tr}(C^\top P).
\]

Define trace-moment constraints

\[
h_k(P)=\operatorname{tr}(P^k)-b_k.
\]

For a single directed Hamiltonian cycle, a natural target is

\[
b_k=0,\qquad 1\le k<n,
\]

and

\[
b_n=n.
\]

The constraints \(\operatorname{tr}(P^k)=0\) for \(k<n\) suppress shorter closed cycles, while \(\operatorname{tr}(P^n)=n\) is consistent with a single \(n\)-cycle permutation.

### Remark: why Birkhoff rather than orthogonal?

On an orthogonal or isospectral manifold, trace-power constraints mostly constrain eigenvalues and leave a large continuous orbit. On the nonnegative stochastic/Birkhoff side, trace powers have a direct closed-walk interpretation. This makes the trace constraints substantially more combinatorial.

## 3. KL/Mirror Geometry on the Birkhoff Polytope

A natural flow on \(\mathcal B_n^\circ\) is the KL or entropic mirror-gradient flow. Given a scalar energy \(E(P)\) with Euclidean gradient

\[
G=\nabla_P E(P),
\]

the Birkhoff-preserving mirror flow can be written in coordinates as

\[
\dot P_{ij}=-P_{ij}\left(G_{ij}-a_i-b_j\right),
\]

where the row and column potentials \(a,b\) are chosen so that

\[
\dot P\mathbf 1=0,\qquad \dot P^\top\mathbf 1=0.
\]

Equivalently, the additive potentials remove the forbidden normal directions associated with row and column marginal constraints.

### Multiplier equations for the Birkhoff gauge

Let

\[
r_i=\sum_j P_{ij}G_{ij},\qquad c_j=\sum_i P_{ij}G_{ij}.
\]

Preservation of row and column sums gives

\[
a+Pb=r,
\]

\[
P^\top a+b=c.
\]

Equivalently,

\[
\begin{bmatrix}
I&P\\
P^\top&I
\end{bmatrix}
\begin{bmatrix}
a\\b
\end{bmatrix}
=
\begin{bmatrix}
r\\c
\end{bmatrix}.
\]

There is a gauge freedom \(a\mapsto a+\gamma\mathbf 1\), \(b\mapsto b-\gamma\mathbf 1\), which can be fixed by imposing, for example, \(\mathbf 1^\top a=0\).

### TODO proof sketch

- Show that the coordinate flow preserves strict positivity.
- Show that the chosen \(a,b\) preserve row and column sums.
- Identify the flow as the KL/natural-gradient projection of \(G\) onto the Birkhoff tangent space.
- State precisely the metric/operator \(M_P\) induced by the KL geometry.

## 4. Cycle-Moment Sinkhorn Flow

Define the augmented trace-penalty energy

\[
E_\rho(P)=\langle C,P\rangle+\frac{\rho}{2}\sum_{k=1}^n h_k(P)^2.
\]

Then

\[
G=\nabla_P E_\rho(P)
\]

is

\[
G=C+\rho\sum_{k=1}^n h_k(P)\,k(P^{k-1})^\top.
\]

The penalty flow is

\[
\dot P_{ij}=-P_{ij}\left(G_{ij}-a_i-b_j\right).
\]

This is the continuous-time analogue of exponentiated-gradient descent followed by Sinkhorn projection.

## 5. Primal-Dual and Augmented Primal-Dual Variants

Instead of a pure penalty, introduce live multipliers \(\lambda_k\) and define the Lagrangian

\[
\mathcal L(P,\lambda)=\langle C,P\rangle+\sum_{k=1}^n\lambda_k h_k(P).
\]

The primal-dual mirror flow is

\[
\dot P=-\operatorname{grad}^{\mathrm{KL}}_{\mathcal B}\mathcal L(P,\lambda),
\]

\[
\dot\lambda_k=\eta h_k(P).
\]

In coordinates,

\[
\dot P_{ij}=-P_{ij}\left(G_{ij}-a_i-b_j\right),
\]

with

\[
G=C+\sum_{k=1}^n\lambda_k k(P^{k-1})^\top.
\]

The augmented primal-dual version is

\[
\mathcal L_\rho(P,\lambda)=\langle C,P\rangle+\sum_{k=1}^n\lambda_k h_k(P)+\frac{\rho}{2}\sum_{k=1}^n h_k(P)^2,
\]

with

\[
\dot P=-\operatorname{grad}^{\mathrm{KL}}_{\mathcal B}\mathcal L_\rho(P,\lambda),
\]

\[
\dot\lambda_k=\eta h_k(P).
\]

The augmented term is expected to damp trace-constraint ringing.

### TODO proof sketch

- Derive the primal-dual linearisation.
- Derive the augmented primal-dual linearisation.
- Show how the penalty term contributes a \(\rho A^\top A\)-type damping/stiffness term in trace-violating directions.

## 6. Trace Derivatives and Constraint Jacobian

The first derivative of the trace moments is

\[
D\operatorname{tr}(P^k)[\delta P]
=k\operatorname{tr}(P^{k-1}\delta P).
\]

Equivalently,

\[
\nabla_P \operatorname{tr}(P^k)=k(P^{k-1})^\top.
\]

Let

\[
A=Dh(P)
\]

be the trace-constraint Jacobian, with rows

\[
A_k[\delta P]=k\operatorname{tr}(P^{k-1}\delta P).
\]

The second variation is

\[
D^2\operatorname{tr}(P^k)[\delta P,\delta P]
=
k\sum_{r=0}^{k-2}\operatorname{tr}\left(P^r\delta P P^{k-2-r}\delta P\right).
\]

Thus the local curvature of the trace-constrained Lagrangian is generated entirely by the trace-moment terms; the cost \(\langle C,P\rangle\) is linear and contributes no Hessian.

### TODO proof sketch

- Derive the first derivative using cyclicity of trace.
- Derive the second variation carefully for noncommuting perturbations.
- Express the Hessian as a linear operator on the Birkhoff tangent space.

## 7. Lyapunov Structure

For the pure penalty energy

\[
E_\rho(P)=\langle C,P\rangle+\frac{\rho}{2}\|h(P)\|^2,
\]

and the KL-gradient flow

\[
\dot P=-\operatorname{grad}^{\mathrm{KL}}_{\mathcal B}E_\rho(P),
\]

we expect the dissipation identity

\[
\frac{dE_\rho}{dt}
=-\left\|\operatorname{grad}^{\mathrm{KL}}_{\mathcal B}E_\rho(P)\right\|^2_{\mathrm{KL}}
\le 0.
\]

Thus the penalty flow is Lyapunov stable in the sense of monotone energy descent. It cannot exhibit conservative ringing in the same way as the pure primal-dual flow.

The pure primal-dual system, by contrast, has saddle-like structure and can exhibit oscillatory constraint modes.

### TODO proof sketch

- Define the precise KL metric on \(T_P\mathcal B\).
- Show that the coordinate expression for \(\dot P\) is the negative KL gradient.
- Prove the energy dissipation identity.
- Contrast with primal-dual dynamics by differentiating the linearised constraint residuals.

## 8. Linearisation Around Interior Stationary Points

Let \((P_*,\lambda_*)\) be an interior stationary point. Perturbations satisfy

\[
\delta P\in T_{P_*}\mathcal B,
\]

where

\[
T_{P_*}\mathcal B=\{\delta P:\delta P\mathbf 1=0,\;\delta P^\top\mathbf 1=0\}.
\]

For the primal-dual flow, the linearisation has the schematic block form

\[
\begin{bmatrix}
\delta\dot P\\
\delta\dot\lambda
\end{bmatrix}
=
\begin{bmatrix}
-M_*H_*&-M_*A_*^\top\\
\eta A_*&0
\end{bmatrix}
\begin{bmatrix}
\delta P\\
\delta\lambda
\end{bmatrix}.
\]

For the augmented primal-dual flow, the linearisation becomes

\[
\begin{bmatrix}
\delta\dot P\\
\delta\dot\lambda
\end{bmatrix}
=
\begin{bmatrix}
-M_*(H_*+\rho A_*^\top A_*)&-M_*A_*^\top\\
\eta A_*&0
\end{bmatrix}
\begin{bmatrix}
\delta P\\
\delta\lambda
\end{bmatrix}.
\]

More precisely, the Euclidean expression \(A^\top A\) should be replaced by the corresponding metric-aware operator involving \(M_*\), for example \(A_*M_*A_*^\top\) in the dual constraint space.

The key qualitative consequence is:

\[
\rho A_*^\top A_*
\]

adds damping/stiffness in directions that change the trace moments.

### Validity versus optimality modes

Perturbations split into trace-violating and trace-tangent modes:

\[
A_*\delta P\ne 0,
\]

and

\[
A_*\delta P=0.
\]

The augmented penalty directly controls the first class. The second class is invisible to the trace constraints at first order and must be controlled by constrained curvature, entropy/annealing effects, cost structure, or boundary geometry.

This suggests a conceptual separation:

\[
\text{trace penalties stabilise validity,}
\]

while

\[
\text{cost and reduced-cost structure select optimality among valid cycles.}
\]

### TODO proof sketch

- Define \(M_*\) rigorously on the Birkhoff tangent space.
- Derive the block Jacobian for primal-dual flow.
- Derive the augmented block Jacobian.
- State the mode decomposition in metric-aware form.
- Prove the damping estimate for trace-violating perturbations.

## 9. Degeneracy at the Uniform Doubly Stochastic Point

Let

\[
P_0=\frac{1}{n}\mathbf 1\mathbf 1^\top.
\]

Then

\[
P_0^k=P_0,
\qquad k\ge 1.
\]

For \(k\ge 2\), the trace derivative along a Birkhoff tangent perturbation \(\delta P\) is

\[
D\operatorname{tr}(P^k)[\delta P]
=k\operatorname{tr}(P_0\delta P)
=\frac{k}{n}\sum_{ij}\delta P_{ij}.
\]

But for \(\delta P\in T_{P_0}\mathcal B\), row and column sums vanish, hence

\[
\sum_{ij}\delta P_{ij}=0.
\]

Therefore

\[
D\operatorname{tr}(P^k)[\delta P]=0,
\qquad k\ge 2.
\]

This suggests that most trace-power constraints are first-order invisible at the uniform doubly stochastic point. The optimiser may therefore receive weak or degenerate trace-gradient signal from exactly symmetric initialisation.

A practical implication is that one should not initialise exactly at uniform logits. Small symmetry-breaking perturbations may be necessary.

### TODO proof sketch

- Check the \(k=1\) case separately.
- Account for diagonal/no-self-loop constraints if used.
- Determine the exact rank of \(A=Dh(P_0)\) restricted to \(T_{P_0}\mathcal B\).
- Study whether second-order trace information becomes visible at uniform.

## 10. Boundary Behaviour and Reduced-Cost Stability Near Tours

Permutation tours lie on the boundary of the Birkhoff polytope. The KL flow degenerates at the boundary because

\[
\dot P_{ij}=-P_{ij}\left(G_{ij}-a_i-b_j\right),
\]

so if \(P_{ij}=0\), then \(\dot P_{ij}=0\). Exact permutation vertices therefore require separate tangent-cone or log-coordinate analysis.

Near a permutation support, for an off-support edge with small positive mass,

\[
\dot P_{ij}\approx -P_{ij}\Delta_{ij},
\]

where

\[
\Delta_{ij}=G_{ij}-a_i-b_j.
\]

Thus

\[
\frac{d}{dt}\log P_{ij}\approx -\Delta_{ij}.
\]

If \(\Delta_{ij}>0\), the off-support edge decays exponentially in log coordinates. If \(\Delta_{ij}<0\), that edge tends to grow.

At a permutation vertex, feasible directions of the Birkhoff polytope correspond to alternating-cycle exchange directions. Therefore local stability of a candidate tour should be governed by the sign of the augmented reduced cost along all feasible exchange directions:

\[
\langle G_*,D\rangle>0.
\]

This connects the continuous flow to combinatorial local-search moves such as cycle exchanges and \(k\)-opt-like perturbations.

### TODO proof sketch

- Formalise near-boundary asymptotics using \(P_\epsilon\) with small off-support mass.
- Derive the log-coordinate reduced-cost equation.
- State the tangent cone of the Birkhoff polytope at a permutation vertex.
- Relate tangent-cone generators to alternating cycles.
- Clarify how trace constraints restrict or modify exchange directions.

## 11. Practical Autodiff/Sinkhorn Recipe

The practical differentiable relaxation uses logits \(X\) and a temperature \(\tau\):

\[
P=\operatorname{Sinkhorn}\left(\exp(X/\tau)\right).
\]

Then optimise

\[
E_{\rho,\tau}(X)=\langle C,P\rangle+\frac{\rho}{2}\sum_{k=1}^n h_k(P)^2.
\]

The stability analysis suggests:

1. Do not initialise exactly at \(X=0\), since the trace constraints may be first-order degenerate at the uniform matrix.
2. Use small random logits to break symmetry.
3. Increase \(\rho\) gradually to damp trace-violating modes without making the system excessively stiff.
4. Decrease \(\tau\) gradually to encourage collapse toward permutation vertices.
5. Monitor reduced costs \(\Delta_{ij}=G_{ij}-a_i-b_j\) near convergence as a local stability diagnostic.

A possible continuation schedule is

\[
\rho_t\uparrow,
\qquad
\tau_t\downarrow.
\]

The analysis suggests that \(\rho\) should be tuned in relation to the conditioning of the trace-constraint Jacobian, especially through quantities like

\[
A M A^\top.
\]

## 12. Numerical Investigation

The experiments are diagnostic rather than benchmark-oriented. The aim is not to compete with specialised TSP solvers, but to test the stability predictions of the trace-constrained Sinkhorn-flow picture:

1. trace penalties enforce Hamiltonicity/validity modes;
2. cost and annealing select among valid cycles;
3. symmetric initialisation can suppress first-order trace signal;
4. pure primal-dual trace dynamics has undamped oscillatory modes;
5. augmented penalties damp trace-visible modes.

All experiments use small instances, typically \(n=8\), so that exact Hamiltonian optima and exact assignment projections can be computed by brute force. The script writes machine-readable CSV summaries in addition to plots.

### 12.1 Friendly circle sanity check

For points on a circle, the optimal Euclidean tour is the clockwise or anticlockwise nearest-neighbour cycle. This gives a low-noise sanity check that the relaxation can recover an obvious Hamiltonian cycle.

Two configurations are compared:

- a continuation schedule with \(\rho\uparrow\) and \(\tau\downarrow\);
- a fixed \(\rho,\tau\) run.

The continuation run recovers the geometric optimum. The fixed run can still produce a single Hamiltonian cycle, but with long chords and crossings. This demonstrates the separation between validity and optimality:

\[
\text{trace moments can enforce a single cycle without selecting the best cycle.}
\]

Key files:

- [circle cost](cycle_moment_outputs/circle_cost.png)
- [circle trace residual norm](cycle_moment_outputs/circle_h_norm.png)
- [scheduled final assignment](cycle_moment_outputs/circle_scheduled_assignment_perm.png)
- [fixed final assignment](cycle_moment_outputs/circle_fixed_assignment_perm.png)
- [circle trace residual components, scheduled](cycle_moment_outputs/circle_rho_up___tau_down_h_components.png)
- [circle trace residual components, fixed](cycle_moment_outputs/circle_fixed_rho_tau_h_components.png)
- [circle entropy](cycle_moment_outputs/circle_entropy.png)
- [circle sharpness](cycle_moment_outputs/circle_sharpness.png)

### 12.2 Two-cluster subtour pressure

The two-cluster instance is designed so that the cheap assignment structure is a pair of within-cluster subtours. This stresses the trace constraints because the cost term favours a cycle cover rather than a single Hamiltonian cycle.

Observed behaviour:

- small \(\rho\): cheap subtour covers such as \([4,4]\);
- moderate \(\rho\): transition to a single \(8\)-cycle;
- very large \(\rho\): trace residuals are small, but the selected single cycle can be unnecessarily expensive.

This supports the interpretation of \(\rho\) as a Hamiltonicity pressure. It also shows that excessive pressure can enforce validity before the cost landscape has selected a good tour.

Key files:

- [two-cluster trace residual norm](cycle_moment_outputs/two_cluster_h_norm.png)
- [two-cluster cost](cycle_moment_outputs/two_cluster_cost.png)
- [two-cluster entropy](cycle_moment_outputs/two_cluster_entropy.png)
- [two-cluster sharpness](cycle_moment_outputs/two_cluster_sharpness.png)
- [scheduled assignment, remains two subtours](cycle_moment_outputs/two_cluster_scheduled_assignment_perm.png)
- [small-rho assignment, two subtours](cycle_moment_outputs/two_cluster_small_assignment_perm.png)
- [large-rho assignment, single cycle](cycle_moment_outputs/two_cluster_large_assignment_perm.png)
- [small-rho trace residual components](cycle_moment_outputs/two_cluster_fixed_small_rho_h_components.png)
- [large-rho trace residual components](cycle_moment_outputs/two_cluster_fixed_large_rho_h_components.png)
- [scheduled trace residual components](cycle_moment_outputs/two_cluster_rho_up___tau_down_h_components.png)

### 12.3 Rho-pressure sweep

The fixed-\(\rho\) sweep makes the Hamiltonicity-pressure transition explicit. For each \(\rho\), the optimiser is run on the two-cluster instance and the final matrix is projected to a permutation by exact assignment projection,

\[
\pi^*=\arg\max_\pi \sum_i P_{i,\pi(i)}.
\]

The aggregate plots show median behaviour and interquartile ranges across seeds. The most important curve is the probability that the projected assignment is a single cycle.

The observed qualitative picture is:

\[
\rho\text{ low}\Rightarrow\text{subtour cover},
\]

\[
\rho\text{ intermediate}\Rightarrow\text{single cycle with moderate cost},
\]

\[
\rho\text{ high}\Rightarrow\text{single cycle but possibly poor route}.
\]

Key files:

- [single-cycle probability](cycle_moment_outputs/two_cluster_rho_pressure_single_cycle_probability.png)
- [final trace residual norm](cycle_moment_outputs/two_cluster_rho_pressure_h_norm_aggregate.png)
- [assignment-projected cost](cycle_moment_outputs/two_cluster_rho_pressure_assignment_cost_aggregate.png)
- [relaxed cost](cycle_moment_outputs/two_cluster_rho_pressure_relaxed_cost_aggregate.png)
- [entropy](cycle_moment_outputs/two_cluster_rho_pressure_entropy_aggregate.png)
- [sharpness](cycle_moment_outputs/two_cluster_rho_pressure_sharpness_aggregate.png)
- [rho-pressure summary CSV](cycle_moment_outputs/two_cluster_rho_pressure_summary.csv)

The key empirical observation is a sharp transition: increasing \(\rho\) collapses trace residuals and turns cheap subtour covers into single-cycle assignments, but the assignment-projected cost rises sharply once Hamiltonicity is forced.

### 12.4 Rho/tau phase diagram

The two-parameter phase diagram varies both trace pressure \(\rho\) and Sinkhorn temperature \(\tau\). It shows regions corresponding to:

- soft/fractional matrices;
- cheap subtour covers;
- useful single-cycle solutions;
- brittle high-pressure/high-sharpness solutions.

In the two-cluster example, \(\rho\) controls most of the transition from subtours to single cycles, while \(\tau\) mainly controls entropy and row sharpness. High temperatures require larger trace pressure to force single cycles; low temperatures sharpen assignments but can freeze a subtour structure if trace pressure is insufficient.

Key files:

- [phase diagram: single-cycle indicator](cycle_moment_outputs/two_cluster_phase_assignment_single_cycle.png)
- [phase diagram: trace residual norm](cycle_moment_outputs/two_cluster_phase_h_norm.png)
- [phase diagram: assignment-projected cost](cycle_moment_outputs/two_cluster_phase_assignment_cost.png)
- [phase diagram: entropy](cycle_moment_outputs/two_cluster_phase_entropy.png)
- [phase diagram: sharpness](cycle_moment_outputs/two_cluster_phase_sharpness.png)
- [phase diagram summary CSV](cycle_moment_outputs/two_cluster_phase_diagram_summary.csv)

This experiment motivates a staged stabilisation rule: avoid lowering \(\tau\) so early that subtours freeze before the trace pressure has enough leverage.

### 12.5 Uniform initialisation degeneracy

At the uniform doubly stochastic point

\[
P_0=\frac{1}{n}\mathbf 1\mathbf 1^\top,
\]

we have \(P_0^k=P_0\) for \(k\ge1\). For Birkhoff-tangent perturbations \(\delta P\),

\[
D\operatorname{tr}(P^k)[\delta P]
= k\operatorname{tr}(P_0\delta P)
= \frac{k}{n}\sum_{ij}\delta P_{ij}=0,
\qquad k\ge2.
\]

With diagonal/self-loops masked, the remaining \(k=1\) trace signal is also removed. The numerical diagnostic confirms that the trace-penalty gradient is essentially zero at exact uniform initialisation and grows only after symmetry-breaking logit perturbations.

Key files:

- [uniform initialisation degeneracy diagnostic](cycle_moment_outputs/uniform_degeneracy_grad_norm.png)
- [ambient trace-Jacobian singular values](cycle_moment_outputs/uniform_ambient_trace_jacobian_singular_values.txt)

The practical implication is simple: do not initialise exactly at uniform logits for no-self-loop trace-constrained routing.

### 12.6 Linearised primal-dual ringing and augmented damping

Around a frozen Sinkhorn point, let \(A=Dh(X)\) be the trace-moment Jacobian with respect to logits. The small-signal primal-dual model is

\[
\dot z=-A^\top\lambda-\rho A^\top A z,
\]

\[
\dot\lambda=\eta A z.
\]

Along a singular mode of \(A\) with singular value \(\sigma\), the trace-visible scalar mode obeys

\[
\ddot q + \rho\sigma^2\dot q + \eta\sigma^2 q=0.
\]

Thus \(\rho=0\) gives an undamped oscillator, while \(\rho>0\) gives a damped oscillator. The numerical demo normalises by the top singular value so that the ringing is visible on a fixed time axis.

Key files:

- [linearised top-mode ringing and damping](cycle_moment_outputs/ringing_linear_top_mode.png)
- [linear trace residual norm](cycle_moment_outputs/ringing_linear_h_norm.png)
- [linear pure trace components](cycle_moment_outputs/ringing_linear_pure_rho0_linear_h_components.png)
- [linear damped trace components](cycle_moment_outputs/ringing_linear_damped_zeta0.22_linear_h_components.png)
- [trace-Jacobian singular spectrum](cycle_moment_outputs/ringing_trace_jacobian_spectrum.png)
- [linear ringing scaling parameters](cycle_moment_outputs/ringing_linear_scaling.txt)

This is the cleanest numerical confirmation of the linearised stability analysis.

### 12.7 Nonlinear primal-dual snapping

The full nonlinear Sinkhorn/logit primal-dual system does not always show smooth sinusoidal ringing. In the experiments it can instead sit at high trace residual for many steps, accumulate dual pressure, and then abruptly collapse to a sharp support. This is a nonlinear analogue of the same lack of damping: the trace-control energy is not dissipated smoothly until support hardening occurs.

Key files:

- [nonlinear trace residual norm](cycle_moment_outputs/ringing_h_norm.png)
- [pure primal-dual trace components](cycle_moment_outputs/ringing_pure_primal-dual_rho=0_h_components.png)
- [soft-temperature primal-dual trace components](cycle_moment_outputs/ringing_pure_primal-dual_soft_tau=4_h_components.png)
- [augmented trace components](cycle_moment_outputs/ringing_augmented_rho=25_h_components.png)
- [nonlinear entropy](cycle_moment_outputs/ringing_entropy.png)
- [nonlinear sharpness](cycle_moment_outputs/ringing_sharpness.png)
- [nonlinear relaxed cost](cycle_moment_outputs/ringing_cost.png)

The augmented run removes most of this delayed-collapse behaviour by damping trace-visible modes directly.

### 12.8 Summary of numerical evidence

The experiments support the stability analysis rather than a solver-performance claim:

1. continuation can recover an obvious geometric tour;
2. fixed trace pressure can produce valid but poor tours;
3. clustered instances reveal cheap subtour attractors;
4. increasing \(\rho\) causes a sharp transition from subtours to single cycles;
5. exact uniform initialisation is first-order degenerate under no-self-loop trace constraints;
6. the linearised primal-dual trace modes ring without augmentation;
7. augmented penalties damp the trace-visible modes.

The resulting practical recipe is:

\[
\text{small random logits},
\qquad
\rho\uparrow\text{ carefully},
\qquad
\tau\downarrow\text{ only after trace residuals respond}.
\]

## 13. Relevance Beyond TSP

Although TSP motivates the trace moments, the stability analysis is more general. Many differentiable routing, matching, optimal transport, sparse-attention, and assignment systems use a Sinkhorn or doubly stochastic relaxation:

\[
P=\operatorname{Sinkhorn}(X).
\]

If additional structural constraints are imposed,

\[
h(P)=0,
\]

then the same local analysis applies. The geometry of the optimiser is controlled by:

1. the KL/Sinkhorn metric on the transport polytope;
2. the constraint Jacobian \(A=Dh(P)\);
3. the spectrum of the metric-aware constraint operator \(A M A^\top\);
4. the behaviour near symmetric interior points;
5. reduced-cost stability near boundary/permutation-like supports.

Thus the TSP trace-moment case is a concrete example of a broader phenomenon:

\[
\text{Sinkhorn routing stability is governed by the interaction of mirror geometry and structural-constraint Jacobians.}
\]

## 14. Relation to Brockett Flow

The analogy with Brockett's double-bracket flow is structural rather than literal.

Classical Brockett-style flows:

\[
\text{orthogonal/isospectral geometry}
\]

\[
\text{commutator projection}
\]

\[
\text{sorting/diagonalisation by a potential}
\]

The present flow:

\[
\text{Birkhoff/transport-polytope geometry}
\]

\[
\text{Sinkhorn/KL projection via row-column potentials}
\]

\[
\text{cycle-moment shaping and routing by a cost potential}
\]

In Brockett flow, commutators remove forbidden normal directions on an isospectral orbit. Here, row/column potentials remove forbidden marginal directions on the Birkhoff polytope.

A concise description is:

> A transport-polytope analogue of Brockett flow, where isospectral commutator geometry is replaced by KL/Sinkhorn geometry, and Hamiltonian cycles arise as trace-moment constrained boundary strata.

## 15. Open Questions

1. Give a clean self-contained proof/citation that the trace-power constraints characterise single Hamiltonian cycles over the Birkhoff polytope under nonnegativity and no-self-loop assumptions.
2. What is the exact rank of the trace-constraint Jacobian at symmetric points?
3. How severe is the first-order degeneracy near uniform initialisation in practical finite-temperature Sinkhorn optimisation?
4. Can the spectrum of \(A M A^\top\) be used to adaptively tune augmentation strength \(\rho\)?
5. Can reduced-cost diagnostics predict failure modes of differentiable routing systems more generally?
6. Is there a clean theorem connecting local stability of a tour under this flow to absence of improving alternating-cycle exchanges?
7. Does the augmented primal-dual flow admit a clean global Lyapunov function, or only local damping interpretations?
8. Can this framework be extended to rectangular transport polytopes, partial matchings, capacity-constrained routing, or sparse attention maps?

## 16. Summary

The core object is the trace-constrained KL-gradient flow on the Birkhoff polytope:

\[
\dot P_{ij}=-P_{ij}\left(G_{ij}-a_i-b_j\right),
\]

with

\[
G=C+\rho\sum_k h_k(P)k(P^{k-1})^\top
\]

or, in primal-dual form,

\[
G=C+\sum_k\lambda_k k(P^{k-1})^\top+\rho\sum_k h_k(P)k(P^{k-1})^\top.
\]

The main stability insights are:

1. The pure penalty flow has a Lyapunov energy.
2. Pure primal-dual dynamics can ring.
3. Augmentation adds a \(\rho A^\top A\)-type damping term to trace-violating modes.
4. Trace constraints may be first-order degenerate at the uniform doubly stochastic point.
5. Near permutation vertices, local support stability is governed by reduced augmented edge costs.
6. Feasible boundary perturbations correspond to alternating-cycle exchange directions.
7. These facts suggest practical stabilisation rules for Sinkhorn-based routing optimisers.

The mathematical picture is therefore not simply a TSP relaxation, but a Brockett-inspired modern flow on transport geometry: a cycle-moment Sinkhorn flow.

