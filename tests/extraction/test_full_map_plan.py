from pathlib import Path
import subprocess
import sys


def test_exact_full_map_plan_and_shard_coverage():
    root = Path(__file__).resolve().parents[2]
    code = '''
from collections import Counter
import numpy as np
from run_full_resolution_maps import plan, HAMILTONIAN_FAMILIES
rows = plan('hamiltonian')
assert len(rows) == 7140
assert len({r['id'] for r in rows}) == 7140
assert Counter(r['family'] for r in rows) == {key:60*n for key,n in HAMILTONIAN_FAMILIES}
assert sorted(r['lambda'] for r in rows[:60]) == list(np.linspace(0,1,60))
assert rows[60]['level'] == np.linspace(.30,5.25,50)[1]
parts = [rows[i::20] for i in range(20)]
assert len({r['id'] for part in parts for r in part}) == 7140
assert sum(map(len,parts)) == 7140
assert len(plan('tpms')) == 1323
print('exact plans verified')
'''
    subprocess.run([sys.executable, "-c", code], cwd=root / "dev", check=True)
