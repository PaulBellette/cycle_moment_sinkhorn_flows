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

### Lemma: trace moments characterise Hamiltonian cycles

Let \(P\in\mathcal B_n\). Then \(P\) is the permutation matrix of a single directed \(n\)-cycle if and only if

\[
\operatorname{tr}(P^k)=0,\qquad k=1,\dots,n-1.
\]

In that case, automatically \(\operatorname{tr}(P^n)=n\).

Proof. If \(P\) is a single \(n\)-cycle permutation matrix, then \(P^k\) has no fixed points for \(1\le k<n\). Hence \(\operatorname{tr}(P^k)=0\).

Conversely, suppose \(P\in\mathcal B_n\) and \(\operatorname{tr}(P^k)=0\) for \(1\le k<n\). Since \(P\ge0\),

\[
\operatorname{tr}(P^k)
=
\sum_{i_1,\dots,i_k}
P_{i_1i_2}P_{i_2i_3}\cdots P_{i_ki_1}.
\]

Every term is nonnegative. Therefore \(\operatorname{tr}(P^k)=0\) implies there is no positive-weight closed walk of length \(k\).

By the Birkhoff-von Neumann theorem,

\[
P=\sum_\ell \alpha_\ell \Pi_\ell,
\qquad
\alpha_\ell>0,\quad
\sum_\ell\alpha_\ell=1,
\]

where each \(\Pi_\ell\) is a permutation matrix. If any \(\Pi_\ell\) contained a cycle of length \(m<n\), then that cycle would give a positive term in \(\operatorname{tr}(P^m)\), contradiction. Hence each \(\Pi_\ell\) is a single \(n\)-cycle.

It remains to show that all \(\Pi_\ell\) are the same \(n\)-cycle. Suppose two distinct \(n\)-cycles occur in the decomposition. Relabel the vertices so that one of them is

\[
1\to2\to\cdots\to n\to1.
\]

The other cycle has some edge \(i\to j\) with \(j\ne i+1\) modulo \(n\). Combining this edge with the path in the first cycle from \(j\) back to \(i\) gives a directed cycle of length strictly less than \(n\) in the support of \(P\). This would again make \(\operatorname{tr}(P^m)>0\) for some \(m<n\), contradiction.

Therefore all permutation matrices in the Birkhoff decomposition are the same \(n\)-cycle, and \(P\) equals that permutation matrix. Finally, for an \(n\)-cycle permutation, \(P^n=I\), so \(\operatorname{tr}(P^n)=n\).

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

### Proposition/proof plan: Birkhoff-preserving KL projection

We want to justify three related facts about this coordinate form.

First, the tangent space of the Birkhoff interior is the set of matrices whose row and column sums vanish. This follows by differentiating the constraints \(P(t)\mathbf 1 = \mathbf 1\) and \(P(t)^\top\mathbf 1 = \mathbf 1\). So any feasible velocity \(V\) must satisfy

\[
V\mathbf 1=0,\qquad V^\top\mathbf 1=0.
\]

Second, the KL/Fisher metric on positive matrices is

\[
\langle U,V\rangle_P
=
\sum_{ij}\frac{U_{ij}V_{ij}}{P_{ij}}.
\]

With no row/column constraints, the KL-gradient of an energy \(E\) with Euclidean gradient \(G\) is \(P\odot G\), because

\[
DE(P)[V]=\sum_{ij}G_{ij}V_{ij}
\]

must equal

\[
\langle P\odot G,V\rangle_P
=
\sum_{ij}\frac{P_{ij}G_{ij}V_{ij}}{P_{ij}}
=
\sum_{ij}G_{ij}V_{ij}.
\]

Thus unconstrained KL-gradient descent has the multiplicative form

\[
\dot P_{ij}=-P_{ij}G_{ij}.
\]

Third, the Birkhoff constraints are enforced by subtracting additive row and column potentials before multiplying by \(P\). Set

\[
V_{ij}=-P_{ij}(G_{ij}-a_i-b_j).
\]

Row preservation gives

\[
0=\sum_jV_{ij}
=
-\sum_jP_{ij}G_{ij}
+a_i
+\sum_jP_{ij}b_j.
\]

Writing

\[
r_i=\sum_jP_{ij}G_{ij},
\]

this is

\[
a+Pb=r.
\]

Similarly, column preservation gives

\[
P^\top a+b=c,
\]

where

\[
c_j=\sum_iP_{ij}G_{ij}.
\]

This is the block system above. It is singular only because the velocity depends on \(a_i+b_j\), so the transformation

\[
a\mapsto a+\gamma\mathbf 1,
\qquad
b\mapsto b-\gamma\mathbf 1
\]

changes the potentials but not the flow. A gauge such as \(\sum_i a_i=0\) fixes this.

To verify that this is the KL-projected gradient, define

\[
W_{ij}=P_{ij}(G_{ij}-a_i-b_j).
\]

For any tangent direction \(V\),

\[
\langle W,V\rangle_P
=
\sum_{ij}(G_{ij}-a_i-b_j)V_{ij}.
\]

Expanding,

\[
\langle W,V\rangle_P
=
\sum_{ij}G_{ij}V_{ij}
-
\sum_i a_i\sum_jV_{ij}
-
\sum_j b_j\sum_iV_{ij}.
\]

The last two terms vanish because \(V\) has zero row and column sums. Therefore

\[
\langle W,V\rangle_P
=
DE(P)[V]
\]

for every tangent direction \(V\). This is exactly the defining property of the KL-gradient restricted to the Birkhoff tangent space.

Finally, strict positivity follows directly from the multiplicative form. Each entry satisfies

\[
\dot P_{ij}=P_{ij}R_{ij}(t),
\]

where

\[
R_{ij}(t)=-(G_{ij}-a_i-b_j).
\]

Therefore

\[
P_{ij}(t)
=
P_{ij}(0)
\exp\left(\int_0^t R_{ij}(s)\,ds\right).
\]

Thus positive entries remain positive while the solution exists and the vector field remains finite.

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

### Linearisation reference

The local linearisation and damping calculation are carried out in Section 8. The role of this section is to define the primal-dual and augmented primal-dual systems; the stability mechanism is easier to see after introducing the trace Jacobian and local tangent coordinates.

In brief, near an interior stationary point and in local coordinates \(z\) on the Birkhoff tangent space, the trace constraints linearise as

\[
h(P_*+z)\approx A z.
\]

Ignoring cost curvature for the moment, the pure primal-dual trace subsystem has the normal form

\[
\dot z=-A^\top\lambda,
\qquad
\dot\lambda=\eta A z.
\]

The augmented system adds the local penalty gradient

\[
\rho A^\top A z,
\]

so that

\[
\dot z=-A^\top\lambda-\rho A^\top A z,
\qquad
\dot\lambda=\eta A z.
\]

Along a singular mode \(A v_i=\sigma_i u_i\), this gives

\[
\ddot q+\rho\sigma_i^2\dot q+\eta\sigma_i^2q=0.
\]

Thus the pure primal-dual model has undamped trace-control modes, while the augmented model damps them. Section 8 gives the derivation and the caveat that the full nonlinear Sinkhorn/logit dynamics can instead appear as delayed snapping or support hardening rather than a clean sinusoid.

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

### Proof

Let \(f_k(P)=\operatorname{tr}(P^k)\). Expanding to first order,

\[
(P+\epsilon \delta P)^k
=
P^k
+
\epsilon\sum_{r=0}^{k-1}P^r\delta P P^{k-1-r}
+
O(\epsilon^2).
\]

Taking traces gives

\[
\operatorname{tr}\left((P+\epsilon\delta P)^k\right)
=
\operatorname{tr}(P^k)
+
\epsilon
\sum_{r=0}^{k-1}
\operatorname{tr}(P^r\delta P P^{k-1-r})
+
O(\epsilon^2).
\]

By cyclicity of trace,

\[
\operatorname{tr}(P^r\delta P P^{k-1-r})
=
\operatorname{tr}(P^{k-1}\delta P).
\]

Thus

\[
Df_k(P)[\delta P]
=
k\operatorname{tr}(P^{k-1}\delta P).
\]

Using the Euclidean matrix inner product
\(\langle A,B\rangle=\operatorname{tr}(A^\top B)\), this implies

\[
\nabla_P f_k(P)=k(P^{k-1})^\top.
\]

For the second variation, differentiate the first variation. Since

\[
D(P^{k-1})[\delta P]
=
\sum_{r=0}^{k-2}P^r\delta P P^{k-2-r},
\]

we obtain

\[
D^2 f_k(P)[\delta P,\delta P]
=
k\sum_{r=0}^{k-2}
\operatorname{tr}\left(
P^r\delta P P^{k-2-r}\delta P
\right).
\]

Since \(h_k(P)=f_k(P)-b_k\), the same derivatives hold for \(h_k\).

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

### Proof

Let \(E(P)\) be any smooth energy and let

\[
W=\operatorname{grad}^{\mathrm{KL}}_{\mathcal B}E(P).
\]

The KL-gradient flow is

\[
\dot P=-W.
\]

By the chain rule,

\[
\frac{dE}{dt}=DE(P)[\dot P].
\]

Using the defining property of the KL-gradient,

\[
DE(P)[V]=\langle W,V\rangle_P
\]

for every tangent direction \(V\). Taking \(V=\dot P=-W\), we obtain

\[
\frac{dE}{dt}
=
\langle W,-W\rangle_P
=
-\|W\|_P^2
\le0.
\]

In coordinates,

\[
W_{ij}=P_{ij}(G_{ij}-a_i-b_j),
\]

so

\[
\|W\|_P^2
=
\sum_{ij}\frac{W_{ij}^2}{P_{ij}}
=
\sum_{ij}P_{ij}(G_{ij}-a_i-b_j)^2.
\]

Therefore

\[
\frac{dE}{dt}
=
-\sum_{ij}P_{ij}(G_{ij}-a_i-b_j)^2.
\]

For \(E=E_\rho\), this proves the stated Lyapunov dissipation identity.

The pure primal-dual flow does not satisfy the same identity for the Lagrangian. If

\[
\mathcal L(P,\lambda)=\langle C,P\rangle+\lambda^\top h(P),
\]

then

\[
\frac{d\mathcal L}{dt}
=
D_P\mathcal L[\dot P]+h(P)^\top\dot\lambda.
\]

The first term is nonpositive under KL-gradient descent in \(P\), but with

\[
\dot\lambda=\eta h(P),
\]

the second term is

\[
\eta\|h(P)\|^2\ge0.
\]

Thus the primal-dual system is not simply a dissipative gradient flow of \(\mathcal L\), which is consistent with the oscillatory modes seen in the linearised analysis.

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

Here \(M_*\) denotes the KL mobility/projection operator at \(P_*\). Concretely, for a Euclidean covector \(Q\), \(M_*Q\) is the Birkhoff-tangent vector

\[
(M_*Q)_{ij}=P_{*,ij}(Q_{ij}-\alpha_i-\beta_j),
\]

where \(\alpha,\beta\) are chosen so that \(M_*Q\) has zero row and column sums. Thus \(M_*\) is the metric-dependent analogue of projecting an ordinary Euclidean gradient onto the tangent space. In the simplified local-coordinate model below, this operator is absorbed into the choice of coordinates, which is why the normal form uses \(A^\top A\) directly.

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

### Local trace-mode model: ringing and damping

The full nonlinear Sinkhorn/logit dynamics can snap, saturate, or change support. The following calculation should therefore be read as a local normal-form model for the trace-visible modes near an interior stationary point, not as a complete global convergence theorem.

Let \(z\) denote local coordinates for perturbations in the Birkhoff tangent space, and linearise the trace constraints as

\[
h(P_*+z)\approx A z,
\]

where \(A=Dh(P_*)\) is the trace-constraint Jacobian in these local coordinates. Ignoring cost curvature for the moment, the pure primal-dual trace dynamics has the local form

\[
\dot z=-A^\top\lambda,
\]

\[
\dot\lambda=\eta A z.
\]

Differentiating the second equation gives

\[
\ddot\lambda
=
\eta A\dot z
=
-\eta A A^\top\lambda.
\]

Thus the trace-control subsystem has undamped oscillatory modes. Equivalently, let

\[
A v_i=\sigma_i u_i
\]

be a singular mode of \(A\), and write

\[
z=qv_i,
\qquad
\lambda=pu_i.
\]

Then

\[
\dot q=-\sigma_i p,
\]

\[
\dot p=\eta\sigma_i q.
\]

Hence

\[
\ddot q=-\eta\sigma_i^2 q,
\]

or

\[
\ddot q+\eta\sigma_i^2 q=0.
\]

This is the linearised ringing mechanism.

For the augmented primal-dual system, the local augmented term is

\[
\frac{\rho}{2}\|Az\|^2,
\]

whose gradient is

\[
\rho A^\top A z.
\]

The local model becomes

\[
\dot z=-A^\top\lambda-\rho A^\top A z,
\]

\[
\dot\lambda=\eta A z.
\]

In the same singular mode,

\[
\dot q=-\sigma_i p-\rho\sigma_i^2 q,
\]

\[
\dot p=\eta\sigma_i q.
\]

Differentiating the \(q\)-equation gives

\[
\ddot q
=
-\sigma_i\dot p-\rho\sigma_i^2\dot q
=
-\eta\sigma_i^2q-\rho\sigma_i^2\dot q.
\]

Therefore

\[
\ddot q+\rho\sigma_i^2\dot q+\eta\sigma_i^2q=0.
\]

So the augmentation adds a damping coefficient \(\rho\sigma_i^2\) in each trace-visible singular mode.

Including cost/Lagrangian curvature gives the more general local model

\[
\dot z=-H z-A^\top\lambda-\rho A^\top A z,
\]

\[
\dot\lambda=\eta A z.
\]

The curvature term \(H\) can shift frequencies and add local decay or growth depending on the constrained Hessian. However, the trace-augmentation contribution remains the positive semidefinite damping/stiffness term

\[
\rho A^\top A.
\]

This is the mechanism isolated in the linearised ringing experiment. In the full nonlinear Sinkhorn parameterisation the same undamped trace-control structure may appear not as a clean sinusoid, but as delayed high-residual plateaus followed by sharp support changes.

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

### Proposition: first-order degeneracy at the uniform point

Let

\[
P_0=\frac{1}{n}\mathbf 1\mathbf 1^\top
\]

and let

\[
\delta P\in T_{P_0}\mathcal B_n
=
\{\delta P:\delta P\mathbf 1=0,\;\delta P^\top\mathbf 1=0\}.
\]

Then for every \(k\ge2\),

\[
D\operatorname{tr}(P^k)\big|_{P=P_0}[\delta P]=0.
\]

Proof. Since \(P_0^m=P_0\) for every \(m\ge1\), the trace derivative gives

\[
D\operatorname{tr}(P^k)\big|_{P=P_0}[\delta P]
=
k\operatorname{tr}(P_0^{k-1}\delta P)
=
k\operatorname{tr}(P_0\delta P).
\]

Using \(P_0=(1/n)\mathbf 1\mathbf 1^\top\),

\[
\operatorname{tr}(P_0\delta P)
=
\frac1n\operatorname{tr}(\mathbf 1\mathbf 1^\top\delta P)
=
\frac1n\mathbf 1^\top\delta P\mathbf 1.
\]

But \(\delta P\mathbf 1=0\), so

\[
\mathbf 1^\top\delta P\mathbf 1=0.
\]

Therefore

\[
D\operatorname{tr}(P^k)\big|_{P=P_0}[\delta P]=0,
\qquad k\ge2.
\]

The case \(k=1\) is different:

\[
D\operatorname{tr}(P)[\delta P]=\operatorname{tr}(\delta P).
\]

This is not automatically zero for Birkhoff-tangent perturbations. Thus, if diagonal/self-loop entries are allowed, the \(k=1\) trace term can provide a first-order signal at the uniform point. If the diagonal is hard-masked, or if self-loops are removed from the parameterisation, this \(k=1\) direction is absent, and the trace-moment penalty is first-order silent at \(P_0\).

Consequently, in the no-self-loop setting, exact uniform initialisation is a degenerate point for the trace penalties. Small random logit perturbations are needed to break symmetry and reveal trace-gradient information.

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

### Local log-rate and reduced-cost interpretation

For any strictly positive entry of the KL flow,

\[
\dot P_{ij}
=
-P_{ij}(G_{ij}-a_i-b_j).
\]

Therefore

\[
\frac{d}{dt}\log P_{ij}
=
-(G_{ij}-a_i-b_j).
\]

Define the reduced augmented cost

\[
\Delta_{ij}=G_{ij}-a_i-b_j.
\]

Then

\[
\frac{d}{dt}\log P_{ij}=-\Delta_{ij}.
\]

Near a permutation tour, off-support entries are small but positive if we remain in the Birkhoff interior. For such an off-support edge,

\[
P_{ij}(t)\approx P_{ij}(0)e^{-\Delta_{ij}t}
\]

as long as \(\Delta_{ij}\) varies slowly. Thus positive reduced cost suppresses that edge, while negative reduced cost makes it grow.

At the exact boundary, however, the KL flow degenerates: if \(P_{ij}=0\), then \(\dot P_{ij}=0\). Exact permutation matrices are therefore absorbing supports for the multiplicative flow. Boundary stability is better interpreted through near-boundary limits or log coordinates.

The feasible perturbations of a permutation vertex are not arbitrary single-edge perturbations. They must preserve row and column sums, so they are generated by alternating-cycle exchange directions in the assignment polytope. If \(D\) is such a feasible exchange direction, then

\[
D\mathbf 1=0,\qquad D^\top\mathbf 1=0.
\]

Consequently,

\[
\langle \Delta,D\rangle
=
\sum_{ij}(G_{ij}-a_i-b_j)D_{ij}
=
\sum_{ij}G_{ij}D_{ij}
=
\langle G,D\rangle,
\]

because the row and column potential terms vanish on feasible directions.

Thus a candidate tour is locally stable against a feasible exchange direction \(D\) when

\[
\langle G,D\rangle>0.
\]

Equivalently, every allowed alternating-cycle exchange should have positive augmented reduced cost. If some feasible exchange direction has

\[
\langle G,D\rangle<0,
\]

then the augmented energy decreases along that exchange direction, and the near-boundary flow has a mechanism for leaking mass toward the corresponding non-tour edges.

This should be read as a local reduced-cost diagnostic rather than a global optimality theorem. A full theorem would require a precise tangent-cone statement for the Birkhoff vertex and a careful treatment of the trace constraints on the exchange directions.

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


## 15. Theory Roadmap: Prove, Cite, or Treat as Interpretation

This section is a working plan for turning the current derivation sketch into a more disciplined note. The aim is not to formalise every geometric remark at once, but to separate the claims into:

- results we should prove cleanly in the note;
- standard background we should cite;
- local models we can derive under stated simplifications;
- interpretive claims that should remain narrative or future work.

### 15.1 Claim discipline table

| Claim or ingredient | Treatment | Tactic |
|---|---:|---|
| The vertices of the Birkhoff polytope are permutation matrices. | Cite | Cite Birkhoff-von Neumann. No proof needed. |
| Trace powers of a nonnegative adjacency matrix count weighted closed walks. | Cite + prove briefly | Include the expansion of \(\operatorname{tr}(P^k)\) as a sum over closed walks; cite acyclicity/graphical-model literature if desired. |
| Trace-power constraints over the Birkhoff polytope characterise Hamiltonian-cycle permutations. | Cite + proof sketch | Cite trace-guided ATSP / continuous acyclicity literature, but include a short support-graph proof sketch because this is central to the note. |
| The coordinate flow \(\dot P_{ij}=-P_{ij}(G_{ij}-a_i-b_j)\) preserves positivity. | Prove | Positivity follows from multiplicative form: \(\dot P_{ij}/P_{ij}\) is finite while the solution remains smooth. |
| The same flow preserves row and column sums when \(a,b\) solve the gauge equations. | Prove | Sum over rows and columns to derive \(a+Pb=r\), \(P^\top a+b=c\). |
| The \(a,b\) system has gauge freedom. | Prove | Show \(a\mapsto a+\gamma\mathbf 1,\;b\mapsto b-\gamma\mathbf 1\) leaves \(a_i+b_j\) unchanged. Fix a gauge such as \(\mathbf 1^\top a=0\). |
| The flow is the KL/natural-gradient projection of \(G\) onto the Birkhoff tangent space. | Prove lightly + cite | Define the Fisher/KL metric \(\langle U,V\rangle_P=\sum_{ij}U_{ij}V_{ij}/P_{ij}\), then show the projected gradient has the coordinate form. Cite Sinkhorn/mirror-flow literature for context. |
| First derivative of \(\operatorname{tr}(P^k)\). | Prove | Use cyclicity of trace. |
| Second variation of \(\operatorname{tr}(P^k)\). | Prove | Expand \(P+\epsilon \delta P\) and collect second-order terms, taking care of noncommutation. |
| Pure penalty flow has a Lyapunov energy. | Prove | Once the KL gradient is defined, \(dE/dt=-\|\operatorname{grad}^{\mathrm{KL}}_{\mathcal B}E\|_{\mathrm{KL}}^2\le0\). |
| Full metric-aware Jacobian at a stationary point has the schematic block form in Section 8. | Derive schematically | State this as the principal linearisation in local tangent coordinates. Avoid overclaiming until the metric/projection derivative is written carefully. |
| The local primal-dual trace mode obeys a harmonic oscillator equation when \(\rho=0\). | Prove in local coordinates | Use the simplified small-signal system \(\dot z=-A^\top\lambda,\;\dot\lambda=\eta A z\), then diagonalise by the SVD of \(A\). |
| The augmented term produces damping. | Prove in local coordinates | With \(\dot z=-A^\top\lambda-\rho A^\top A z,\;\dot\lambda=\eta A z\), a singular mode satisfies \(\ddot q+\rho\sigma^2\dot q+\eta\sigma^2q=0\). |
| The split \(A\delta P\neq0\) versus \(A\delta P=0\) separates trace-visible and trace-tangent modes. | Prove as linear algebra; interpret carefully | The penalty contributes first-order stiffness only through \(A\). Calling these “validity” and “optimality” modes is an interpretation. |
| Trace gradients are first-order invisible at the uniform point for \(k\ge2\). | Prove | At \(P_0=J/n\), \(P_0^k=P_0\), and \(D\operatorname{tr}(P^k)[\delta P]=(k/n)\sum_{ij}\delta P_{ij}=0\) on the Birkhoff tangent space. |
| The \(k=1\) trace signal at uniform is special. | Prove / qualify | \(D\operatorname{tr}(P)[\delta P]=\operatorname{tr}(\delta P)\), not automatically zero on the Birkhoff tangent space. With diagonal masking/no-self-loop constraints this signal is removed. |
| Off-support edge masses obey a reduced-cost log-rate near a tour. | Prove local formula | For positive entries, \(d(\log P_{ij})/dt=-(G_{ij}-a_i-b_j)\). Near a permutation support, apply this to small off-support entries. |
| Tangent cone at a permutation vertex is generated by alternating cycles. | Cite / optional proof sketch | Standard assignment-polytope fact. We can cite and use it to interpret boundary exchange modes. |
| Local stability of a tour is equivalent to positive augmented reduced cost on all exchange directions. | Heuristic / future theorem | Plausible and useful, but needs careful tangent-cone and trace-constraint handling before presenting as a theorem. |
| “Modern Brockett” interpretation. | Narrative | Keep as structural analogy, not a theorem. |
| General relevance to Sinkhorn routing. | Discussion + optional generic proposition | The generic \(h(P)\) local model is formal; claims about routing/attention applications should stay as discussion. |
| Adaptive \(\rho\) based on \(A M A^\top\). | Future work | Motivated by analysis and numerics, but not yet established. |

### 15.2 Suggested proof order and current status

The most efficient order is to prove the facts that unlock several sections at once. In this draft, the main proof sketches have now been incorporated: the remaining work is mostly citation polishing, deciding how formal to make the manifold linearisation, and deciding which interpretive claims should remain future work.

#### Step 1: Birkhoff-preserving KL projection

Target sections: 3, 4, 7.

Status: proof plan included in Section 3.

Prove that the coordinate flow

\[
\dot P_{ij}=-P_{ij}(G_{ij}-a_i-b_j)
\]

preserves \(\mathcal B_n^\circ\) when \(a,b\) solve the gauge equations. Then identify it as the KL-projected negative gradient on the Birkhoff polytope.

Deliverables:

- proposition statement;
- proof of row/column preservation;
- proof of positivity;
- gauge-freedom remark;
- metric interpretation.

#### Step 2: Trace calculus

Target sections: 4, 6, 8, 9.

Status: proof included in Section 6.

Prove the first derivative

\[
D\operatorname{tr}(P^k)[\delta P]
=
k\operatorname{tr}(P^{k-1}\delta P),
\]

and the second variation

\[
D^2\operatorname{tr}(P^k)[\delta P,\delta P]
=
k\sum_{r=0}^{k-2}
\operatorname{tr}(P^r\delta P P^{k-2-r}\delta P).
\]

Deliverables:

- lemma statement;
- proof by cyclic trace manipulation;
- note that the linear cost contributes no Hessian.

#### Step 3: Hamiltonian-cycle characterisation

Target sections: 2 and 15.

Status: proof sketch included in Section 2; still needs citation polish. Cite the relevant trace-acyclicity / trace-guided ATSP literature.

Proof route:

1. For \(P\ge0\), \(\operatorname{tr}(P^k)\) is a sum of nonnegative weighted closed walks of length \(k\).
2. Hence \(\operatorname{tr}(P^k)=0\) forbids all positive closed walks of length \(k\).
3. A doubly stochastic support graph must contain directed cycles.
4. If all cycles of length \(<n\) are forbidden, the only possible support cycle length is \(n\).
5. With \(\operatorname{tr}(P^n)=n\), the support and weights collapse to a single \(n\)-cycle permutation.

The last collapse step should be written carefully; it is the only part where handwaving is tempting.

#### Step 4: Lyapunov identity

Target section: 7.

Status: proof included in Section 7.

Once the KL gradient has been defined, show

\[
\frac{dE_\rho}{dt}
=
-\left\|\operatorname{grad}^{\mathrm{KL}}_{\mathcal B}E_\rho(P)\right\|^2_{\mathrm{KL}}
\le0.
\]

This should be a short corollary.

#### Step 5: Local primal-dual damping model

Target sections: 5, 8, 12.6.

Status: local normal-form model included in Section 8; Section 5 now points forward to it.

In local coordinates \(z\) on the Birkhoff tangent space, derive the simplified small-signal model

\[
\dot z=-H z-A^\top\lambda-\rho A^\top A z,
\]

\[
\dot\lambda=\eta A z.
\]

For the clean trace-mode picture, set \(H=0\) or work in a regime where the constraint coupling dominates. If \(A v=\sigma u\), then the trace-visible mode satisfies

\[
\ddot q+\rho\sigma^2\dot q+\eta\sigma^2 q=0.
\]

This exactly matches the linearised ringing/damping figure.

#### Step 6: Uniform degeneracy

Target sections: 9 and 12.5.

Status: proposition included in Section 9.

Prove the first-order degeneracy at \(P_0=J/n\). This is short, central, and strongly supported by the numerical diagnostic.

Deliverables:

- proposition for \(k\ge2\);
- separate remark for \(k=1\);
- no-self-loop/diagonal-mask corollary.

#### Step 7: Boundary reduced-cost analysis

Target section: 10.

Status: local log-rate formula and reduced-cost interpretation included in Section 10; tangent-cone theorem remains future work.

Prove the log-rate formula

\[
\frac{d}{dt}\log P_{ij}=-(G_{ij}-a_i-b_j)
\]

for positive entries. Then interpret near-permutation behaviour using small off-support masses.

For the exchange-direction story:

- cite the alternating-cycle tangent-cone fact;
- state the reduced-cost stability condition as a conjecture or proposition with conditions to be filled in later.

### 15.3 Reference targets

These are the references we probably need, grouped by role. Exact bibliographic details can be filled later.

| Role | Reference target |
|---|---|
| Brockett motivation | Brockett, “Dynamical systems that sort lists, diagonalize matrices and solve linear programming problems.” |
| Birkhoff vertices | Birkhoff-von Neumann theorem / standard convex-polytopes or combinatorial optimisation text. |
| Trace powers and acyclicity | NOTEARS / continuous acyclicity constraints, plus standard graph closed-walk trace facts. |
| Trace-guided TSP relaxation | Recent trace-guided ATSP / doubly stochastic trace-constraint work. |
| Sinkhorn as mirror descent/flow | Sinkhorn mirror descent, entropic OT, continuous-time Sinkhorn flow literature. |
| Assignment-polytope tangent cone | Standard assignment polytope / network-flow / combinatorial optimisation reference. |
| Reduced costs | Linear programming / assignment problem reduced-cost optimality. |

### 15.4 What not to overclaim yet

The following ideas are useful, but should remain carefully worded until proved.

1. **Global convergence to the TSP optimum.**  
   The flow has stationary points and bad local structures; our claim is about geometry and stability diagnostics, not solver optimality.

2. **Full manifold Jacobian without caveats.**  
   The schematic block matrix is correct as a local model, but a fully rigorous manifold linearisation requires careful treatment of the KL metric and projection dependence on \(P\).

3. **Reduced-cost stability as a theorem.**  
   The log-rate formula is easy; the full equivalence between local attraction and no improving exchange cycles requires a precise tangent-cone statement.

4. **General Sinkhorn-routing claims.**  
   The generic local model is real, but application claims should remain discussion unless we test or prove them.

5. **Adaptive \(\rho\) schedules.**  
   The spectrum of \(A M A^\top\) is a plausible control signal, but we have only motivation and diagnostics so far.


## 16. Open Questions

1. Find the cleanest citation trail for the trace-moment Hamiltonian-cycle characterisation, complementing the self-contained proof in Section 2.
2. What is the exact rank of the trace-constraint Jacobian at symmetric points?
3. How severe is the first-order degeneracy near uniform initialisation in practical finite-temperature Sinkhorn optimisation?
4. Can the spectrum of \(A M A^\top\) be used to adaptively tune augmentation strength \(\rho\)?
5. Can reduced-cost diagnostics predict failure modes of differentiable routing systems more generally?
6. Is there a clean theorem connecting local stability of a tour under this flow to absence of improving alternating-cycle exchanges?
7. Does the augmented primal-dual flow admit a clean global Lyapunov function, or only local damping interpretations?
8. Can this framework be extended to rectangular transport polytopes, partial matchings, capacity-constrained routing, or sparse attention maps?

## 17. Summary

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

