# Direct continuum equivalence of compact TPMS solids

## Result and scope

This experimental extension compares the analytic source solids directly,
without requiring two extracted graph spines to have identical coordinates.
The branch is research/continuum-isotopy-fab2673, based on audit commit
fab2673f665796c65a32ec42efdc1bf9b6c37678. Existing production files, default
behavior and the original manuscript are unchanged.

All nine corridors chosen by the earlier contraction experiment pass the new
whole-domain regularity test and certificate replay. All nine also pass replay
with independent rational-Taylor trigonometric bounds. The cavity-bearing
P-to-D corridor is certified as a family of solids, not as the invalid legacy
singleton graph. The exact selection matches the archived
`tpms_audit/embedding_contraction_pair_probe.json`.

The second experiment tests all 1,200 closed rectangles between the original
21-by-21 samples in three families. It is not a 1,200-point voxel scan:

| Family | Certified rectangles | Unknown rectangles | Total |
|---|---:|---:|---:|
| Schwarz-P to Diamond | 332 | 68 | 400 |
| Gyroid to Schwarz-P | 317 | 83 | 400 |
| Gyroid to Diamond | 265 | 135 | 400 |
| Total | 914 | 286 | 1200 |

Every successful serialized certificate is replayed. All nine corridors and
all 914 successful grid rectangles also pass replay with independent rational
trigonometric bounds. Unknown does not mean inequivalent. A singular event, an overly
conservative interval enclosure or a finite resource budget can prevent this
sufficient test from accepting a rectangle.

## Same fields, same finite domain

For spatial position r=(x,y,z), define

    G = sin(x)cos(y) + sin(y)cos(z) + sin(z)cos(x),
    P = cos(x) + cos(y) + cos(z),
    D = cos(x)cos(y)cos(z) - sin(x)sin(y)sin(z),
    f(r,lambda,c) = (1-lambda)F0(r) + lambda F1(r) - c,
    Omega(lambda,c) = {r : f(r,lambda,c) <= 0 and |r| <= R}.

These match the audited `_tpms.py` definitions. The radius uses precisely the
same binary64 constant, `0.72*(2.25*math.pi) = 5.0893800988154645`, hex
`0x1.45b8674e54d7cp+2`. Parameter bounds retain the actual numpy.linspace binary
values. The statement concerns these analytic trigonometric fields and this
fixed ball, not a voxel interpolant, a printed rounded radius or a periodic
unit cell. The nonsmooth maximum of the field and domain constraint is not
mistaken for a smooth function: the two constraints are checked separately.

## Sufficient isotopy theorem

Let K be a closed parameter rectangle. Suppose f is smooth near B_R x K.
Assume that for every parameter in K:

1. Every point of f=0 in B_R has nonzero spatial gradient.
2. On f=0 intersect the spherical wall, r cross grad_r f is nonzero.

Then any smooth path in K induces an ambient isotopy of the sublevel solids,
preserving B_R setwise. Piecewise smooth paths are treated by concatenation.
The statement permits disconnected solids and cavities; it needs no graph
spine, handlebody assumption or complete polynomial classifier.

### Proof

For a parameter path p(t), put h(r,t)=f(r,p(t)). A velocity X transporting its
zero set satisfies

    grad_r h . X + partial_t h = 0.

In a neighborhood of an interior zero-set point use

    X = -(partial_t h) grad_r h / |grad_r h|^2.

On the spherical wall let

    a_T = grad_r h - ((r . grad_r h)/|r|^2) r.

The second hypothesis makes a_T nonzero. In a sufficiently small collar use
`X = -(partial_t h) a_T / |a_T|^2`. This satisfies the transport equation and
is tangent to the sphere. The moving zero set in the compact space-time
domain has a finite covering by such neighborhoods. Choose bulk neighborhoods
away from the wall and glue their local velocities with a smooth partition of
unity. The transport equation is linear in X, so it survives this gluing,
as does wall tangency. Extend with a cutoff equal to one near the zero set,
then across the wall in a collar, to obtain a compactly supported ambient
vector field.

Its smooth time-dependent flow exists over the compact time interval,
is a diffeomorphism at every time, and preserves the ball because its wall is
invariant. The space-time vector field (X,1) is tangent to h=0. The flow and
its inverse therefore transport exactly that set; a trajectory cannot cross
it. The negative sublevel solid, with all its components, cavities and spatial
placement, is transported to the final solid. This establishes the claim.

The implementation certifies these hypotheses. It does not integrate or store
the flow. The proof is an explicit regular-level argument, not a claim to have
invented a new general isotopy theorem.

### Fixed-box variant

For a rectangular design domain, the code checks all 27 face strata: the
interior, six faces, twelve edges and eight corners. Only free-coordinate
derivatives count on each stratum; corners require f != 0. Local velocities
in the free coordinates can be glued tangent to the incident faces, giving
the analogous face-preserving isotopy. Closed strata are checked redundantly,
which is conservative. Periodic face identification is not implemented.

## Interval certificate and verifier

A shared expression DAG encloses f, its symbolic spatial derivatives,
|r|^2-R^2 and r cross grad_r f on product boxes in space and parameters.
A box is accepted only when it lies outside the ball, excludes f=0, is
strictly inside the ball with a nonzero gradient component, or has a nonzero
cross-product component while it may meet the wall. Strict interval exclusion
is required. Otherwise the box is bisected, including parameter axes when
needed. No sampled-gradient shortcut or endpoint equality is used.

A complete preorder binary tree records every split and leaf. The verifier
reconstructs every box from the caller-supplied trusted field/domain and
recomputes each test. It rejects incomplete covers, extra tree data, invalid
leaves and problem mismatches. Stored leaf verdicts are not trusted. Budget
exhaustion returns unknown and retains a partial tree; it never implies
inequivalence.

The primary backend uses outward-rounded IEEE-754 binary64 operations and an
80-bit private mpmath interval context for trigonometry. Constants are enclosed
against exact rationals. The second backend uses Fraction Taylor sums with
remainder bounds and a Machin-formula enclosure of pi; it uses no mpmath
transcendental values. The backends still share symbolic gradients, partition
logic and elementary outward binary64 arithmetic. These are experimental
numerical certificates under explicit arithmetic/software assumptions, not
formal verification of Python or two completely independent implementations.
Mpmath's own documentation marks its interval facilities experimental.

## Historical corridors and resulting map

At c=0.105, the retained nine paths are:

| Family | Lambda intervals |
|---|---|
| P to D | [0,.05], [.50,.55], [.95,1] |
| G to P | [.05,.10], [.60,.65], [.95,1] |
| G to D | [.05,.10], [.55,.60], [.95,1] |

The printed decimals abbreviate the exact binary endpoints stored in the plans.
The earlier comparison gave eight unknowns and one identical legacy point
pair; the new result concerns the actual analytic solids throughout each path.
It does not retrospectively validate the old extracted graphs.

In the map, connected components of accepted closed rectangles are colored,
and unknown rectangles remain hatched. Shared corners suffice for a certified
path, so eight-neighbor adjacency is appropriate. The accepted cover has 5, 6
and 6 connected components in the three panels. These are NOT counts of
inequivalent topologies. Disconnected colored components may be equivalent or
connected through unresolved regions. Colors are local to each panel and do
not assert cross-family equivalence. An edge of the certified cover is not
necessarily a physical phase boundary. No label filtering is performed.

The symbolic regressions also verify two P-to-D bulk critical-value loci:
at (pi,0,0), c=1-2*lambda; at (pi,pi,0), c=2*lambda-1. Both points lie within
the design ball. These loci are not an exhaustive event set and do not imply
that every point on them separates distinct global topologies.

## Reproduction

On the separate branch with package dependencies installed:

```sh
python -m pytest -q tests/core/test_field_isotopy.py tests/core/test_field_isotopy_rational.py
python dev/certify_tpms_continuation.py --mode historical-nine --out /tmp/kg-nine-new
python dev/replay_nine_rational.py --data /tmp/kg-nine-new --out /tmp/kg-nine-rational-new.json
python dev/certify_tpms_continuation.py --mode grid --samples 21 --jobs 4 --out /tmp/kg-grid-new
python dev/replay_nine_rational.py --data /tmp/kg-grid-new --out /tmp/kg-grid-rational-new.json
python dev/plot_tpms_continuation.py --data /tmp/kg-grid-new --out /tmp/kg-figures-new
```

Use fresh output paths. The downloadable review bundle contains the full
serialized proof trees, plans, records, figures and exact hashes. The GitHub
extension includes source, tests, reproduction scripts and documentation;
the large generated proof-tree bundle is supplied separately, not silently
claimed to be in the code commit. The actual environment is in each plan.
The new tests were run in a partial source tree with PYTHONPATH=src. The full
native build and existing 465-test repository suite were NOT rerun. No
current-code performance claim is inferred from these container timings.

## Remaining problems

This provides direct continuous equivalence without comparing independently
extracted spines. It does not validate any particular volume-to-spine reduction,
make arbitrary polyline smoothing safe, or make Yamada invariant under arbitrary
spine changes. Transporting an independently validated representative spine
along the certified isotopy is a possible next implementation, not completed
here. The 286 unresolved rectangles are not filled or relabeled.

The Hamiltonian 7140-point application, its periodic versus finite-box domain,
the 96 historical seed conflicts, missing timing hardware/model records and
the all-word graph-family proof remain separate tasks. No completion is
claimed for them. Original figures and author/funding placeholders remain
untouched. This is a research extension, not a newly submission-ready paper.

## References and source identity

The source definitions and historical pair plan are in audit commit
fab2673f665796c65a32ec42efdc1bf9b6c37678, respectively
`src/knotted_graph/applications/phase_map_examples/_tpms.py` and
`User_guide/applications/results/arxiv_revision/tpms_audit/embedding_contraction_pair_probe.json`.
The original manuscript already states the bulk and boundary critical equations.

Relevant prior interval-topology work: S. Plantinga and G. Vegter,
*Isotopic Approximation of Implicit Curves and Surfaces*, Symposium on Geometry
Processing (2004), DOI 10.2312/SGP/SGP04/251-260. This work establishes the
relevance of interval reasoning to isotopic approximation; it does not replace
the specific proof and checks above. Mpmath 1.3.0 interval documentation:
https://mpmath.org/doc/current/contexts.html.
