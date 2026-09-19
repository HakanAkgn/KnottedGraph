"""Compact, independently replayed witnesses for the original voxel complex.

The optional C++ program implements exactly the reference free-face rules and
queue order. This runner does not alter the production extraction API. Masks,
full move arrays, terminal geometry, and source hashes are retained, not only
summary statistics. No analytic-to-discrete or regular-neighborhood claim.
"""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

from knotted_graph.extraction import cubical_spine as reference


def decode(ids, shape):
    ids = np.asarray(ids, dtype=np.int64)
    return np.stack((ids // (shape[1] * shape[2]),
                     ids // shape[2] % shape[1], ids % shape[2]), axis=-1)


def reconstruct(mask, kernel, *, gridlines=None, archive=None, timeout=180):
    if sys.byteorder != "little":
        raise RuntimeError("the compact witness encoding currently requires little endian")
    mask, axes, source = reference._input(mask, gridlines)
    kernel = Path(kernel).resolve(strict=True)
    with tempfile.TemporaryDirectory(prefix="kg-collapse-") as temporary:
        directory = Path(temporary)
        mask.astype(np.uint8).tofile(directory / "mask.bin")
        prefix = directory / "witness"
        subprocess.run([str(kernel), *map(str, mask.shape), str(directory / "mask.bin"), str(prefix)],
                       check=True, capture_output=True, timeout=timeout)
        metadata = json.loads(prefix.with_suffix(".json").read_text())
        moves = np.fromfile(prefix.with_suffix(".moves"), dtype="<u4").reshape(-1, 2)
        terminal = np.fromfile(prefix.with_suffix(".terminal"), dtype="<u4")
        shape = metadata["cell_shape"]
        if shape != [2 * n + 1 for n in mask.shape] or len(moves) != metadata["collapses"]:
            raise RuntimeError("kernel output shape mismatch")
        cells = {tuple(map(int, p)) for p in decode(terminal, shape)}
        dimension = max((reference._dim(p) for p in cells), default=-1)
        if reference._counts(cells) != metadata["terminal_cell_counts"]:
            raise RuntimeError("kernel terminal counts mismatch")
        graph = reference._graph(cells, axes) if 0 <= dimension <= 1 else None
        metadata.update({"schema": "knottedgraph.compact_cubical_retraction.v1", "source": source,
                         "terminal_dimension": dimension, "stop_reason": "no_more_free_pairs",
                         "kind": "empty" if dimension < 0 else "graph_retract" if graph is not None else "partial_retract",
                         "kernel_sha256": sha256(kernel.read_bytes()).hexdigest(),
                         "moves_sha256": sha256(moves.tobytes()).hexdigest(),
                         "terminal_complex_sha256": reference._digest(cells),
                         "analytic_source_correspondence_certified": False,
                         "manifold_or_regular_neighborhood_certified": False})
        if archive is not None:
            archive = Path(archive)
            archive.parent.mkdir(parents=True, exist_ok=True)
            # Delta coding reduces archival size but is exactly reversible.
            delta = np.diff(moves.astype(np.int64), axis=0, prepend=np.zeros((1, 2), dtype=np.int64)).astype(np.int32)
            np.savez_compressed(archive, mask=mask, x=axes[0], y=axes[1], z=axes[2],
                                move_delta=delta, terminal=terminal,
                                metadata=np.asarray(json.dumps(metadata, sort_keys=True)))
            metadata["archive_sha256"] = sha256(archive.read_bytes()).hexdigest()
        return graph, metadata, moves, cells


def expand_reference_certificate(mask, axes, metadata, moves, cells):
    """Expand compact output for the existing independent Python verifier."""
    mask, axes, source = reference._input(mask, axes)
    initial = reference._complex(mask, max_cells=20000000)
    pairs = decode(moves, metadata["cell_shape"]).tolist()
    dimension = metadata["terminal_dimension"]
    certificate = {
        "schema": reference._SCHEMA, "source": source, "kind": metadata["kind"],
        "stop_reason": "no_more_free_pairs", "initial_cell_counts": reference._counts(initial),
        "initial_complex_sha256": reference._digest(initial), "collapses": pairs,
        "terminal_cells": [list(c) for c in sorted(cells)],
        "terminal_cell_counts": reference._counts(cells), "terminal_complex_sha256": reference._digest(cells),
        "terminal_dimension": dimension, "analytic_source_correspondence_certified": False,
        "manifold_or_regular_neighborhood_certified": False, "legacy_graph_embedding_validated": False,
    }
    if metadata["kind"] == "graph_retract":
        import networkx as nx
        graph = reference._graph(cells, axes)
        count = nx.number_connected_components(graph)
        certificate["graph_summary"] = {
            "vertices": len(graph), "edges": graph.number_of_edges(), "components": count,
            "cycle_rank": graph.number_of_edges() - len(graph) + count,
            "max_degree": max(dict(graph.degree()).values(), default=0),
        }
    return certificate
