import numpy as np
import pytest

from knotted_graph.periodic import analyze_periodic_mask
from knotted_graph.periodic.fundamental_group import (
    Presentation,
    count_representations,
    periodic_region_presentations,
    s3_table,
)


def _count(mask):
    presentation, = periodic_region_presentations(mask)
    return count_representations(presentation, s3_table())


def _grid(n):
    coordinates = (np.arange(n) + .137) / n - .5
    return np.meshgrid(coordinates, coordinates, coordinates, indexing="ij", sparse=True)


def test_periodic_ambient_group_and_disconnected_ball_components():
    assert _count(np.ones((6, 7, 8), bool)) == 48
    mask = np.zeros((8, 8, 8), bool)
    mask[1, 1, 1] = mask[5, 5, 5] = True
    assert periodic_region_presentations(mask) == (Presentation(0, ()),) * 2


def test_link_sensitive_count_with_equal_incidence_and_winding():
    x, y, z = _grid(32)
    a = (np.hypot(x, y) - .18)**2 + z*z
    b = (np.hypot(x - .18, z) - .14)**2
    linked = np.minimum(a, b + y*y) < .037**2
    unlink = np.minimum(a, b + (y - .29)**2) < .037**2
    assert analyze_periodic_mask(linked).has_same_invariants(analyze_periodic_mask(unlink))
    # Hom(Z^3, S3) = 48. The Hopf and unlink factors contribute 18 and 36.
    assert _count(~linked) == 864
    assert _count(~unlink) == 1728
    assert _count(np.roll(~linked, (11, 19, 27), (0, 1, 2))) == 864


def test_trefoil_count_with_same_coarse_invariants_as_unknot():
    spatial = pytest.importorskip("scipy.spatial")
    x, y, z = _grid(34)
    t = np.linspace(0, 2 * np.pi, 4096, endpoint=False)
    radius = .167 + .066 * np.cos(3*t)
    core = np.column_stack((radius * np.cos(2*t), radius * np.sin(2*t), .066 * np.sin(3*t)))
    samples = np.stack(np.broadcast_arrays(x, y, z), axis=-1)
    trefoil = spatial.cKDTree(core).query(samples)[0] < .03
    unknot = (np.hypot(x, y) - .167)**2 + z*z < .03**2
    assert analyze_periodic_mask(trefoil).has_same_invariants(analyze_periodic_mask(unknot))
    # A trefoil has 12 homomorphisms into S3. An unknot has six.
    assert _count(~trefoil) == 576
    assert _count(~unknot) == 288


def test_count_respects_free_generators_and_search_budget():
    relation = (1, 2, -1, -2)
    table = s3_table()
    assert count_representations(Presentation(3, (relation,)), table) == 108
    with pytest.raises(ValueError, match="36 assignments"):
        count_representations(Presentation(2, (relation,)), table, max_assignments=35)
