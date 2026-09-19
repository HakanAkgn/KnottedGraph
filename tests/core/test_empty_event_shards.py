from pathlib import Path
import subprocess
import sys


def test_empty_event_selection_always_serializes_records(tmp_path):
    root = Path(__file__).resolve().parents[2]
    code = '''
import json
from pathlib import Path
import subprocess
import sys
source=Path(sys.argv[1])/'source'
source.mkdir()
(source/'records.json').write_text('[]')
(source/'plan.json').write_text('{}')
out=Path(sys.argv[1])/'result'
subprocess.run([sys.executable,'certify_tpms_event_candidates.py','--input',str(source),'--out',str(out)],check=True)
assert json.loads((out/'event_records.json').read_text()) == []
summary=json.loads((out/'event_summary.json').read_text())
assert summary['selected_cells'] == summary['recorded_cells'] == 0
assert summary['status_counts'] == {}
'''
    subprocess.run([sys.executable, '-c', code, str(tmp_path)], cwd=root / 'dev', check=True)


def test_only_verified_empty_legacy_shards_can_omit_records(tmp_path):
    root = Path(__file__).resolve().parents[2]
    code = '''
import json
from pathlib import Path
import sys
from collect_saved_tpms_events import read_saved_rows
p=Path(sys.argv[1])
(p/'plan.json').write_text(json.dumps({'selected_ids':[]}))
(p/'event_summary.json').write_text(json.dumps({'selected_cells':0,'recorded_cells':0,'status_counts':{}}))
rows, legacy, _=read_saved_rows(p,[])
assert rows == [] and legacy
assert not (p/'event_records.json').exists()
(p/'plan.json').write_text(json.dumps({'selected_ids':['a']}))
(p/'event_summary.json').write_text(json.dumps({'selected_cells':1,'recorded_cells':1,'status_counts':{'unknown':1}}))
try:
    read_saved_rows(p,['a'])
except ValueError:
    pass
else:
    raise AssertionError('missing nonempty cohort was silently accepted')
(p/'event_records.json').write_text(json.dumps([{'id':'b','status':'unknown'}]))
try:
    read_saved_rows(p,['a'])
except ValueError:
    pass
else:
    raise AssertionError('wrong record identity accepted')
'''
    subprocess.run([sys.executable, '-c', code, str(tmp_path)], cwd=root / 'dev', check=True)


def test_exact_stationary_branch_membership_without_float_tolerance():
    root = Path(__file__).resolve().parents[2]
    code = '''
from fractions import Fraction
import numpy as np
from collect_saved_tpms_events import exact_branch_witness
row={'family':'schwarz_p_to_diamond','lambda_bounds':[.5,.5],'c_bounds':[0.,0.]}
w=exact_branch_witness(row)
assert w is not None
assert Fraction(*w['lambda_fraction']) == Fraction(1,2)
assert Fraction(*w['c_fraction']) == 0
row['c_bounds']=[float(np.nextafter(0.,1.)),1e-6]
assert exact_branch_witness(row) is None
row.update(lambda_bounds=[.4,.42],c_bounds=[.15,.22])
w=exact_branch_witness(row)
assert w is not None
lam,level=Fraction(*w['lambda_fraction']),Fraction(*w['c_fraction'])
assert Fraction(.4)<=lam<=Fraction(.42)
assert Fraction(.15)<=level<=Fraction(.22)
assert level == 1-2*lam
row['family']='gyroid_to_diamond'
assert exact_branch_witness(row) is None
'''
    subprocess.run([sys.executable, '-c', code], cwd=root / 'dev', check=True)
