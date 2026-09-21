import networkx as nx

__all__ = [
    "is_trivalent",
    "BouquetGraph",
    "ThetaGraph",
    "weisfeiler_lehman_multigraph_hash",
]


def is_trivalent(G):
    """
    Check if a graph is trivalent -- if all vertices have degree <= 3.
    
    Parameters:
      G : networkx.MultiGraph / networkx.Graph
         The undirected graph to check.
    
    Returns:
      True if the graph is trivalent, False otherwise.
    """
    degs = nx.degree(G)
    return all(degree <= 3 for node, degree in degs)


def BouquetGraph(n):
    """Construct the bouquet with one vertex and ``n`` loop edges."""
    if isinstance(n, bool) or not isinstance(n, int):
        raise TypeError("n must be a non-negative integer")
    if n < 0:
        raise ValueError("n must be a non-negative integer")
    G = nx.MultiGraph()
    G.add_node(0)
    for _ in range(n):
        G.add_edge(0, 0)
    return G


def ThetaGraph(n):
    """Construct the two-vertex theta multigraph with ``n`` parallel edges."""
    if isinstance(n, bool) or not isinstance(n, int):
        raise TypeError("n must be a non-negative integer")
    if n < 0:
        raise ValueError("n must be a non-negative integer")
    G = nx.MultiGraph()
    G.add_nodes_from((0, 1))
    for _ in range(n):
        G.add_edge(0, 1)
    return G


def weisfeiler_lehman_multigraph_hash(
        g_multi: nx.MultiGraph,
        iters: int = 3,
    ):
    """WL hash that tolerates MultiGraphs by collapsing them first."""
    def _simple_copy_with_multiplicity(
            g_multi: nx.MultiGraph
        ):
        """Collapse a (multi)graph into a simple Graph, recording multiplicity."""
        if not g_multi.is_multigraph():
            g_simple = g_multi.copy()
            for u, v in g_simple.edges():
                g_simple[u][v]["m"] = 1
            return g_simple

        g_simple = nx.Graph() if isinstance(g_multi, nx.MultiGraph) else nx.DiGraph()
        g_simple.add_nodes_from(g_multi.nodes())

        for u, v, _k in g_multi.edges(keys=True):
            if g_simple.has_edge(u, v):
                g_simple[u][v]["m"] += 1  # bump multiplicity
            else:
                g_simple.add_edge(u, v, m=1)
        return g_simple
    
    g_for_hash = _simple_copy_with_multiplicity(g_multi)
    return nx.weisfeiler_lehman_graph_hash(
        g_for_hash,
        iterations=iters,
        edge_attr="m"  # include multiplicity in the label
    )
