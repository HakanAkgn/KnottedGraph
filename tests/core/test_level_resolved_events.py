from pathlib import Path
import subprocess
import sys


def test_augmented_system_replay_and_tamper_rejection():
    root=Path(__file__).resolve().parents[2]
    code='''
from copy import deepcopy
import sympy as s
from level_resolved_tpms_events import SquareSystem,certify,replay,domain_membership
x,y,z,l=s.symbols('x y z l',real=True)
system=SquareSystem([2*x,2*y,2*z,x*x+y*y+z*z+l-s.Rational(1,8)],(x,y,z,l))
center=[0.,0.,0.,.125]
box=[(v-1e-6,v+1e-6) for v in center]
c=certify(system,center,box)
assert c is not None
assert replay(system,c,'primary')['valid']
assert replay(system,c,'rational')['valid']
bad=deepcopy(c);bad['inverse_hex']=[[0.0.hex()]*4 for _ in range(4)]
assert not replay(system,bad)['valid']
bad=deepcopy(c);bad['system_sha256']='wrong'
assert not replay(system,bad)['valid']
row={'lambda_bounds':[.1,.2],'c_bounds':[0.,.2]}
assert domain_membership(row,.125,'bulk',c,2.)
assert not domain_membership(row,.3,'bulk',c,2.)
row['lambda_bounds']=[.15,.2]
assert not domain_membership(row,.125,'bulk',c,2.)
mu=s.Symbol('mu',real=True)
sphere=SquareSystem([-mu*x,-mu*y,1-mu*z,x*x+y*y+z*z-4,l-s.Rational(1,4)],(x,y,z,mu,l))
center=[0.,0.,2.,.5,.25]
c=certify(sphere,center,[(v-1e-6,v+1e-6) for v in center])
assert c is not None and replay(sphere,c)['valid']
print('four/five-dimensional square roots and adversarial checks passed')
'''
    subprocess.run([sys.executable,'-c',code],cwd=root/'dev',check=True)
