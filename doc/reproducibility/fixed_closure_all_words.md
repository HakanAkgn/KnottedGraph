# Fixed-closure all-word derivation and geometry correction

Research extension of `2dbac64674b39ae607ad6696f92f91b03ca3cbbd`.

The proof, independent implementation, tests, and reviewed prior-attempt ledger
are in `dev/fixed_closure_proof/`. This is an explicitly isolated research
module; existing production code and the original manuscript are unchanged.

The explicit fifteen-state planar skein realization reproduces the complete
archived Hankel input file with Git blob identity
`ed8d20c2c97d7cdcf62d4e87cae153e658d71224`, without reading training coefficients.
The exact factorizations H=OC and H_a=OM_aC, together with a nonsingular H,
give the archived T_a=C^{-1}M_aC. This proves the all-word raw polynomial law
for the specified fixed cubic pure-braid diagram family, not the separate
commuting motif families.

A distinct geometric defect was found in the notebook's external closures:
six unintended same-height segment intersections. A valid rational polygonal
realization with corrected exterior nesting is supplied and checked. The
all-word theorem concerns that valid canonical family; it does not certify
old self-intersecting coordinates or arbitrary floating-point executions.

50 standalone tests pass, including all 3375 basis associativity triples,
RII and the braid relation, polynomial convention hashes, and nine old/new
geometry pairs. The production/native suite was not rerun. The full result
matrices, geometry records, test log and four-page report are distributed in
the downloadable review bundle. The code regenerates them; they are not
claimed to be all committed in this directory.

From the repository root:

```sh
cd dev/fixed_closure_proof
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pytest -q test_boundary_algebra.py
python audit_realization.py --out /tmp/kg_fixed_closure_new
python canonical_geometry.py
```

Use a new algebra output directory. This work does not complete a TPMS phase
classification, validate a volume-to-spine reduction, rerun the Hamiltonian
application, recover missing historical records, or finalize the manuscript.
