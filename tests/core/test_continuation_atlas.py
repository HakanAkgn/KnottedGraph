from copy import deepcopy

import pytest
import sympy as sp

from knotted_graph.core.continuation_atlas import ContinuationAtlas
from knotted_graph.core.field_isotopy import FieldProblem, certify


def block(low, high, *, shift=0):
    x, y, z, parameter = sp.symbols("x y z parameter", real=True)
    problem = FieldProblem(x - parameter + shift, (x, y, z, parameter),
                           ((-1., 1.),) * 3 + ((low, high),), "box")
    return problem, certify(problem)


def test_overlap_chain_has_replay_identifiers_and_segments():
    atlas = ContinuationAtlas()
    for interval in ((-.6, -.1), (-.2, .3), (.25, .6)):
        atlas.add(*block(*interval))
    result = atlas.compare((-.5,), (.5,))
    assert result["status"] == "equivalent"
    assert len(result["segments"]) == 3
    assert result["segments"][0]["parameter_entry"] == (-.5,)
    assert result["segments"][-1]["parameter_exit"] == (.5,)
    assert atlas.cover_labels([(-.5,), (0.,), (.5,), (.9,)]) == [0, 0, 0, -1]
    assert result["selected_graph_spines_validated"] is False


def test_disconnected_cover_is_unknown_not_inequivalent():
    atlas = ContinuationAtlas()
    atlas.add(*block(-.5, -.1))
    atlas.add(*block(.1, .5))
    assert atlas.compare((-.3,), (.3,))["status"] == "unknown"
    assert atlas.compare((-.3,), (.3,))["inequivalence_established"] is False
    assert atlas.cover_labels([(-.3,), (.3,)]) == [0, 1]


def test_exact_contact_connects_but_a_real_gap_does_not():
    atlas = ContinuationAtlas()
    atlas.add(*block(-.5, 0.))
    atlas.add(*block(0., .5))
    assert atlas.compare((-.3,), (.3,))["status"] == "equivalent"
    gap = ContinuationAtlas()
    gap.add(*block(-.5, 0.))
    gap.add(*block(1e-12, .5))
    assert gap.compare((-.3,), (.3,))["status"] == "unknown"


def test_unknown_and_tampered_certificates_cannot_enter_cover():
    atlas = ContinuationAtlas()
    problem, proof = block(-.5, .5)
    bad = deepcopy(proof)
    bad["status"] = "unknown"
    with pytest.raises(ValueError):
        atlas.add(problem, bad)
    bad = deepcopy(proof)
    bad["trees"] = []
    with pytest.raises(ValueError):
        atlas.add(problem, bad)
    assert atlas.block_count == 0


def test_changed_source_or_domain_cannot_be_merged():
    atlas = ContinuationAtlas()
    atlas.add(*block(-.5, .5))
    with pytest.raises(ValueError):
        atlas.add(*block(-.5, .5, shift=sp.Rational(1, 10)))


def test_external_certificate_mutation_cannot_change_stored_identity():
    atlas = ContinuationAtlas()
    problem, proof = block(-.5, .5)
    index = atlas.add(problem, proof)
    assert atlas.add(problem, proof) == index
    expected = atlas.compare((-.3,), (.3,))
    proof["trees"].clear()
    assert atlas.compare((-.3,), (.3,)) == expected


@pytest.mark.parametrize("point", [(float("nan"),), (float("inf"),), (), (0., 0.)])
def test_invalid_points_rejected(point):
    atlas = ContinuationAtlas()
    atlas.add(*block(-.5, .5))
    with pytest.raises(ValueError):
        atlas.compare(point, (0.,))


def test_empty_atlas_never_implies_equivalence():
    assert ContinuationAtlas().compare((0.,), (0.,))["status"] == "unknown"
