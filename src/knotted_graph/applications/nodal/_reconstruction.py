"""Shared public nodal reconstruction with explicit geometric evidence."""
from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
import json

import networkx as nx
import numpy as np


def geometry_digest(graph):
    h = sha256()
    for node, data in sorted(graph.nodes(data=True), key=lambda v: repr(v[0])):
        h.update(repr(node).encode())
        h.update(np.asarray(data["pos"], dtype="<f8").tobytes())
    for u, v, key, data in sorted(graph.edges(keys=True, data=True), key=lambda e: repr(e[:3])):
        h.update(repr((u, v, key)).encode())
        h.update(np.asarray(data.get("pts", []), dtype="<f8").tobytes())
    return h.hexdigest()


def counts(graph):
    components = nx.number_connected_components(graph)
    return components, graph.number_of_edges() - graph.number_of_nodes() + components


def reconstruct(model, *, simplify, smooth_epsilon, skeleton_image,
                reconstruction, cubical_options, extract, prune, simplify_edges,
                smooth_edges, is_trivalent):
    if reconstruction not in ("guarded", "cubical"):
        raise ValueError("reconstruction must be guarded or cubical")
    epsilon = float(smooth_epsilon)
    if not np.isfinite(epsilon) or epsilon < 0:
        raise ValueError("smooth_epsilon must be finite and nonnegative")
    options = dict(cubical_options or {})
    if set(options) - {"max_cells", "max_collapses"}:
        raise ValueError("unsupported cubical option")
    if reconstruction == "cubical" and (skeleton_image is not None or epsilon != 0):
        raise ValueError("cubical reconstruction requires the source mask and no geometric smoothing")
    if reconstruction != "cubical" and options:
        raise ValueError("cubical_options require cubical reconstruction")

    target = None
    if skeleton_image is None:
        from knotted_graph.applications.phase_maps import volume_topology
        mask = np.asarray(model._interior_mask, dtype=bool)
        topology = volume_topology(mask)
        target = (topology.connected_components, topology.handle_rank)
        if not mask.any():
            raise ValueError("the source volume is empty")
        input_hash = sha256(mask.astype(np.uint8).tobytes()).hexdigest()
    else:
        mask = None
        image = np.asarray(skeleton_image, dtype=bool)
        if image.ndim != 3:
            raise ValueError("skeleton_image must be a three-dimensional array")
        input_hash = sha256(image.astype(np.uint8).tobytes()).hexdigest()
    key = (reconstruction, bool(simplify), epsilon, input_hash,
           None if mask is None else mask.shape, json.dumps(options, sort_keys=True))
    cached = getattr(model, "skeleton_graph_cache", None)
    if cached is not None and getattr(model, "skeleton_graph_cache_args", None) == key:
        if getattr(model, "_skeleton_graph_digest", None) == geometry_digest(cached):
            return cached

    certificate = None
    if reconstruction == "cubical":
        from knotted_graph.extraction.cubical_spine import (
            cubical_retract, verify_cubical_retract, sample_grid_breaks, certificate_json,
        )
        # Keep the public graph's established index-coordinate convention.
        # The existing physical-coordinate map is affine with positive spacing.
        gridlines = tuple(sample_grid_breaks(np.arange(n, dtype=float)) for n in mask.shape)
        result = cubical_retract(mask, gridlines=gridlines, **options)
        replay = verify_cubical_retract(mask, result.certificate, gridlines=gridlines, **options)
        if not replay["valid"]:
            raise RuntimeError("cubical collapse witness failed replay")
        if result.graph is None:
            raise ValueError("no graph endpoint in the witnessed collapse; result remains unresolved")
        graph, certificate = result.graph, result.certificate
        evidence = {
            "method": "replayed_cubical_collapse", "voxel_retraction_certified": True,
            "certificate_sha256": sha256(certificate_json(certificate).encode()).hexdigest(),
            "collapse_count": replay["collapses"], "coordinate_system": "index",
            "digital_realization": "closed_dual_cells_clipped_at_sampling_endpoints",
        }
    else:
        if mask is not None:
            image = np.asarray(model._skeleton_image, dtype=bool)
            skeleton_topology = volume_topology(image)
            if (skeleton_topology.connected_components, skeleton_topology.handle_rank) != target:
                raise ValueError("source-mask and thinned-skeleton Betti counts disagree")
            graph = extract(image, expected_components=target[0], expected_cycle_rank=target[1])
        else:
            graph = extract(image)
        if simplify:
            original = graph
            candidate = prune(deepcopy(graph))
            # Reinsert an original representative of any entirely pruned tree.
            for component in nx.connected_components(original):
                if not any(n in candidate for n in component):
                    node = min(component, key=repr)
                    candidate.add_node(node, **deepcopy(original.nodes[node]))
            graph = simplify_edges(candidate)
        if epsilon:
            graph = smooth_edges(graph, epsilon=epsilon, copy=False)
        evidence = {
            "method": "homology_guarded_lee" if mask is not None else "external_skeleton",
            "voxel_retraction_certified": False,
            "coordinate_system": "index", "smoothing_epsilon": epsilon,
            "geometric_smoothing_certified": False,
        }
    if target is not None and counts(graph) != target:
        raise ValueError("postprocessed graph disagrees with source component/cycle counts")
    evidence.update({
        "source_mask_sha256": input_hash, "component_cycle_check": target is not None,
        "analytic_source_correspondence_certified": False,
        "periodic_identification": False,
        "source_touches_boundary": bool(topology.touches_boundary) if mask is not None else None,
    })
    graph.graph["reconstruction"] = evidence
    graph.graph["is_trivalent"] = is_trivalent(graph)
    model.is_graph_trivalent = graph.graph["is_trivalent"]
    model.spine_certificate = certificate
    model.skeleton_graph_cache = graph
    model.skeleton_graph_cache_args = key
    model._skeleton_graph_digest = geometry_digest(graph)
    model.__dict__.pop("PDCode", None)
    model._pv_data_args = None
    return graph
