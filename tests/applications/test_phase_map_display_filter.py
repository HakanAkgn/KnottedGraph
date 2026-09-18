import numpy as np

from knotted_graph.applications.phase_map_examples._tpms import stable_labels


def test_protected_failure_cell_is_not_hidden():
    labels = np.ones((3, 3), dtype=int)
    labels[1, 1] = 9
    filtered, count = stable_labels(labels, min_cells=4, protected_labels={9})
    np.testing.assert_array_equal(filtered, labels)
    assert count == 0


def test_small_island_changes_only_display_and_count_is_distinct():
    labels = np.ones((3, 3), dtype=int)
    labels[1, 1] = 9
    filtered, count = stable_labels(labels, min_cells=4)
    assert labels[1, 1] == 9
    assert np.all(filtered == 1)
    assert count == 1


def test_filter_is_invariant_under_label_renaming():
    labels = np.array([[1, 1, 2, 2], [1, 3, 2, 2], [1, 1, 1, 4]])
    renaming = {1: 40, 2: 10, 3: 30, 4: 20}
    renamed = np.vectorize(renaming.__getitem__)(labels)
    filtered, count = stable_labels(labels, min_cells=3)
    other, other_count = stable_labels(renamed, min_cells=3)
    np.testing.assert_array_equal(other, np.vectorize(renaming.__getitem__)(filtered))
    assert count == other_count


def test_tied_neighbors_and_no_large_recipient_remain_unresolved():
    labels = np.array([[1, 1, 3, 2, 2]])
    filtered, count = stable_labels(labels, min_cells=2)
    np.testing.assert_array_equal(filtered, labels)
    assert count == 0
    tiny = np.array([[1, 2], [3, 4]])
    np.testing.assert_array_equal(stable_labels(tiny, min_cells=2)[0], tiny)
