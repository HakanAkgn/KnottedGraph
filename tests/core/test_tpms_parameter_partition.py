from pathlib import Path
import subprocess
import sys


def test_dyadic_parameter_partition_is_exact_and_complete():
    root = Path(__file__).resolve().parents[2]
    code = '''
from refine_tpms_continuation import subdivide
from certify_tpms_continuation import build_plan
for case in build_plan('grid',21):
    parent = {**case,'depth':0,'tile_x':case['lambda_index'],'tile_y':case['c_index']}
    children = subdivide(parent)
    assert len({r['id'] for r in children}) == 4
    assert children[0]['lambda_bounds'][0] == parent['lambda_bounds'][0]
    assert children[1]['lambda_bounds'][1] == parent['lambda_bounds'][1]
    assert children[0]['lambda_bounds'][1] == children[1]['lambda_bounds'][0]
    assert children[0]['c_bounds'][1] == children[2]['c_bounds'][0]
    assert children[2]['c_bounds'][1] == parent['c_bounds'][1]
    assert len({(r['tile_x'],r['tile_y']) for r in children}) == 4
    assert all(r['depth'] == 1 for r in children)
print('all 1200 root partitions checked')
'''
    subprocess.run([sys.executable, "-c", code], cwd=root / "dev", check=True)
