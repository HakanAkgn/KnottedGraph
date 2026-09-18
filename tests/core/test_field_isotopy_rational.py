import math
import random

import mpmath as mp
import pytest

from knotted_graph.core.field_isotopy import certify, tpms_problem
from knotted_graph.core.field_isotopy_rational import _pi_bounds, _trig_interval, verify_rational


def test_machin_enclosure():
    lo, hi = _pi_bounds()
    with mp.workdps(150):
        assert mp.mpf(lo.numerator)/lo.denominator < mp.pi < mp.mpf(hi.numerator)/hi.denominator
        assert mp.mpf(hi.numerator)/hi.denominator-mp.mpf(lo.numerator)/lo.denominator < mp.mpf('1e-85')


@pytest.mark.parametrize('kind', ['sin', 'cos'])
def test_rational_ranges_include_extrema_and_high_precision_points(kind):
    rng = random.Random(1809)
    intervals = [(-8., 8.), (-math.pi, math.pi), (0., 0.),
                 (1.5, 1.6), (-1.6, -1.5), (3.1, 3.2), (-.1, .1)]
    intervals += [tuple(sorted((rng.uniform(-8, 8), rng.uniform(-8, 8)))) for _ in range(60)]
    with mp.workdps(120):
        for a,b in intervals:
            lo, hi = _trig_interval(kind,a,b)
            for i in range(21):
                t = mp.mpf(a)+(mp.mpf(b)-mp.mpf(a))*i/20
                actual = mp.sin(t) if kind=='sin' else mp.cos(t)
                assert mp.mpf(lo) <= actual <= mp.mpf(hi)


def test_out_of_range_is_conservative():
    assert _trig_interval('sin',-9.,9.) == (-1.,1.)
    assert _trig_interval('cos',-math.inf,math.inf) == (-1.,1.)


def test_independent_replay_of_nontrivial_tpms_corridor():
    problem = tpms_problem('schwarz_p_to_diamond',(0.,.05),(.105,.105))
    cert = certify(problem)
    replay = verify_rational(problem,cert)
    assert replay['valid']
    assert replay['mpmath_transcendentals_used'] is False
    cert['trees'][0]['preorder']='.'
    assert not verify_rational(problem,cert)['valid']
