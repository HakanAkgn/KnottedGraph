# Exact stationary branches in the retained Schwarz-P to Diamond family

## Source and executed verification

The retained analytic field is

$$h_\lambda(x,y,z)=(1-\lambda)(\cos x+\cos y+\cos z)
+\lambda(\cos x\cos y\cos z-\sin x\sin y\sin z).$$

The solid is its sublevel set inside the original fixed ball. Its radius is
not rounded from a printed decimal; the stored binary64 value is
`0x1.45b8674e54d7cp+2`.

`dev/verify_exact_tpms_branches.py` differentiates this actual source expression,
substitutes the candidate sets below, and checks every residual symbolically.
It also proves strict inclusion of the stated discrete representatives in the
ball using the independent rational Machin-series enclosure of pi.

The tested commit is `34e243236cccb246a4a4f19ab941203cae30068c`.
GitHub Actions run `35414941316`, job `105821614486`, completed successfully.
The full suite reports **647 passed in 38.55 seconds**. Scoped Ruff checks,
the repository consistency check, and the symbolic verification all passed.
The checked-out tracked source was unchanged by these tests.

The result artifact is `exact-tpms-branch-identities`, ID `10575456891`,
12,711 bytes, SHA-256
`6539f9c5d07ee8d08b7747b8e9a13e91887ccf17c890fbbbc7153afae33c76d8`.
It contains `branch_identities.json` and the complete JUnit result. This
artifact was also exported as `KnottedGraph_exact_TPMS_critical_branches.zip`.

## First exact critical-value branch

At the point $(\pi,0,0)$, direct substitution gives

$$\nabla h_\lambda(\pi,0,0)=0,\qquad
h_\lambda(\pi,0,0)=1-2\lambda,$$

and

$$D^2h_\lambda(\pi,0,0)=\operatorname{diag}
(1,2\lambda-1,2\lambda-1).$$

The same critical value occurs at the six coordinate/sign variants with one
coordinate equal to $\pm\pi$ and the other two zero. For $\lambda<1/2$,
the Hessian has two negative eigenvalues and one positive eigenvalue. The
stationary points are nondegenerate away from $\lambda=1/2$.

The branch is present in the nominal plotted interval $0\le c\le0.3$ for
$0.35\le\lambda\le0.5$. These printed endpoints describe the mathematical
window; cell membership in the numerical certificates uses exact binary
parameter endpoints rather than silently replacing their values by decimals.

## Second exact critical-value branch

At $(\pi,\pi,0)$,

$$\nabla h_\lambda(\pi,\pi,0)=0,\qquad
h_\lambda(\pi,\pi,0)=2\lambda-1,$$

and

$$D^2h_\lambda(\pi,\pi,0)=\operatorname{diag}
(1-2\lambda,1-2\lambda,-1).$$

The twelve coordinate/sign variants with two coordinates equal to $\pm\pi$
have this same critical value. For $\lambda>1/2$ all three eigenvalues are
negative. The stationary points are nondegenerate away from $\lambda=1/2$.
The critical-value branch lies in the nominal displayed level interval for
$0.5\le\lambda\le0.65$.

The rational-domain check establishes $2\pi^2<R^2<3\pi^2$. Therefore both
sets of discrete representatives above lie strictly inside the original ball;
they are not artificial critical points outside the sampled design domain.

## Why the central singularity defeats an isolated-root test

At $\lambda=1/2$, the stationary set is not limited to the discrete points.
For every $t$ satisfying $\pi^2+t^2<R^2$,

$$h_{1/2}(\pi,t,0)=0,\qquad
\nabla h_{1/2}(\pi,t,0)=0.$$

Thus the zero level contains an entire critical line, together with its
coordinate/sign images. The tangent direction is flat. On the normal
$(x,z)$ plane the Hessian is

$$\frac12\begin{pmatrix}
1+\cos t & \sin t\\
\sin t & \cos t-1
\end{pmatrix},\qquad
\det=-\frac12\sin^2t.$$

Away from $\sin t=0$, the two normal eigenvalues have opposite signs, while
the full spatial Hessian retains the zero tangent eigenvalue. At the listed
special points the normal determinant also vanishes.

A local verifier that proves an *isolated unique* zero of the spatial-gradient
system cannot certify a box containing a segment of this critical line. A
regular-level continuation test likewise cannot prove that its gradient is
nonzero. Increasing a search budget does not turn these true singularities
into regular points.

## Consequence for interpreting unresolved parameter cells

The regularity solver and the stationary-event solver establish different
positive statements. A regularity certificate excludes all bulk and boundary
critical events throughout a parameter rectangle. A local critical-point
certificate proves that at least one event exists in a specified part of that
rectangle. The exact branches above give additional symbolic event evidence,
including a non-isolated singular set outside the isolated-root verifier's
scope.

A parameter cell intersecting either exact critical branch cannot be filled
with a regularity label for its entire closed area. It must instead be split
along a justified critical-set enclosure or retain an event/unknown label.
No small-region filtering can supply the missing justification.

These identities do **not** enumerate the entire critical set. Other bulk
stationary points and constrained stationary points on the spherical wall
are not excluded. The number of symmetry-related critical points is not the
number of distinct isotopy classes, and local Hessian signatures alone do not
establish the complete global transition sequence.

No all-word formula, cavity-processing implementation, original manuscript or
main figure was revised in deriving or verifying these statements. This is a
separate source-derived result available for a subsequent manuscript revision.
