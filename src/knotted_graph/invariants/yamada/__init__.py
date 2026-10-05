"""Exact optimized Yamada polynomial APIs."""

from __future__ import annotations

import networkx as nx
import sympy as sp

from .algebra import laurent_y_to_sigma_polynomial
from .compact import CompactYamadaEvaluator
from .polynomial import (
    Yamada,
    compute_yamada_from_states,
)


def compute_graph_yamada_polynomial(
    graph: nx.Graph,
    variable: sp.Symbol,
) -> sp.Expr:
    """Compute the crossing-free Yamada polynomial for an undirected graph.

    Accepts undirected ``networkx.Graph`` and ``networkx.MultiGraph`` inputs.
    Simple graphs are copied into a ``MultiGraph`` for evaluation.
    """
    if not isinstance(graph, nx.Graph):
        raise TypeError("graph must be a networkx.Graph or networkx.MultiGraph")
    if graph.is_directed():
        raise TypeError("graph must be undirected")
    return CompactYamadaEvaluator().compute(nx.MultiGraph(graph), variable)


__all__ = [
    "Yamada",
    "compute_graph_yamada_polynomial",
    "compute_yamada_from_states",
    "laurent_y_to_sigma_polynomial",
]
