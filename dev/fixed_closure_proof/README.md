# Fixed-closure all-word proof and geometry audit

Start with PROOF.md. This standalone research extension uses only SymPy and
NetworkX; pytest runs the exact regressions. It does not import the production
Yamada evaluator, retrieve training coefficients, or edit the manuscript.

From a fresh extracted copy:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pytest -q test_boundary_algebra.py
python audit_realization.py --out /tmp/kg_fixed_closure_new
python canonical_geometry.py
```

Use a new --out path for the algebra audit. The geometry command writes only
this standalone directory's results/geometry_audit.json. Its two modes compare
the legacy exterior control polygons with an explicitly repaired rational PL
realization; it does not assert equivalence to a self-intersecting old curve.

In the downloadable review bundle, final_results contains complete independent
matrices, the byte-matching Hankel payload and the exact audit certificate.
The bundle also includes test output and the contact audit under results/.
The repository contains the source and compact summary; running the commands
above regenerates the detailed evidence. No old native tests, full Hamiltonian
scans, new TPMS maps or historical data recovery are claimed.
