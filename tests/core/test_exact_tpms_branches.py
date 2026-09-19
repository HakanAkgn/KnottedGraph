from pathlib import Path
import subprocess
import sys


def test_exact_branches_and_degenerate_critical_lines_from_original_field():
    root = Path(__file__).resolve().parents[2]
    code = '''
from verify_exact_tpms_branches import verify_identities
result = verify_identities()
assert result['first_branch']['symmetry_related_points'] == 6
assert result['second_branch']['symmetry_related_points'] == 12
assert result['all_symbolic_residuals_zero']
assert result['rational_ball_inclusion_proved']
assert result['degenerate_parameter']['all_gradient_components_identically_zero']
assert result['degenerate_parameter']['isolated_root_uniqueness_inapplicable']
assert not result['all_critical_points_enumerated']
print('exact original-field branch identities and strict source-domain inclusion verified')
'''
    subprocess.run([sys.executable, '-c', code], cwd=root / 'dev', check=True)
