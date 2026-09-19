from pathlib import Path
import subprocess
import sys


def test_certified_event_requires_full_level_enclosure_inside_query():
    root = Path(__file__).resolve().parents[2]
    code = '''
from certify_tpms_event_candidates import event_inside_cell, problem_functions
from knotted_graph.core.field_isotopy import tpms_problem
import sympy as sp
row = {'lambda_bounds':[0.,.1], 'c_bounds':[.1,.2]}
def result(low, high):
    return {'critical_value_hex':[float(low).hex(),float(high).hex()]}
assert event_inside_cell(row,.05,result(.12,.13))
assert not event_inside_cell(row,.2,result(.12,.13))
assert not event_inside_cell(row,.05,result(.05,.15))
assert not event_inside_cell(row,.05,result(.15,.25))
assert not event_inside_cell(row,.05,result(.3,.4))
for family in ('schwarz_p_to_diamond','gyroid_to_schwarz_p','gyroid_to_diamond'):
    problem, *_ = problem_functions(family,.05,'wall')
    source = tpms_problem(family,(.05,.05),(0.,0.))
    expression = source.expression.subs({source.variables[3]:sp.Rational(.05),source.variables[4]:0})
    assert sp.expand(problem.field-expression) == 0
    assert problem.radius == source.radius
print('exact fields and strict query membership checked')
'''
    subprocess.run([sys.executable, '-c', code], cwd=root / 'dev', check=True)
