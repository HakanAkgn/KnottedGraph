"""Guided material and compact-TPMS phase-map applications.

Read and plot saved records without the surface-extraction extras. Computing
new records is an optional, resource-intensive application workflow; see
``python -m knotted_graph.applications.phase_map_examples --help``.

These helpers do not replace the generic ``make_yamada_phase_map`` API.
"""

from .records import PhaseMapData, load_phase_map, read_phase_map_records
from .plotting import plot_phase_map

__all__ = ["PhaseMapData", "load_phase_map", "plot_phase_map", "read_phase_map_records"]
