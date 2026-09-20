"""Exact low-arity factorization for the Yamada connectivity state sum.

A fixed graph vertex is an equality constraint on all incident ports with one
weight -1. High arity is therefore representational, not mathematical. We
factor each such equality into a deterministic chain of arity <=3 equality
factors joined by identity wires. Exactly one factor carries the original -1
weight; auxiliary equality factors carry +1.

Identity wires are logical vertex identifications, not physical graph edges, so
they carry no edge sign. If an identity quotient closes an already selected
physical path, however, the graph cycle rank increases and the state receives
the usual +q = A^-1 + 2 + A cycle factor. An included physical edge carries
-1 and therefore contributes -q when it closes a cycle.

This is the sole production diagram evaluator. It contains no runtime-, family-,
crossing-count-, or benchmark-dependent dispatch. A missing compiled factorized
extension is treated as an installation error rather than silently selecting an
older diagram algorithm.
"""

from __future__ import annotations

from collections import defaultdict

from .diagram_frontier import (
    _greedy_factor_order,
    _greedy_factor_order_from_first,
)

_FACTORIZED_IMPORT_ERROR: Exception | None = None
try:
    from . import _yamada_factorized_frontier
except Exception as exc:  # pragma: no cover - compiler/platform installation guard
    _yamada_factorized_frontier = None
    _FACTORIZED_IMPORT_ERROR = exc

FACTOR_EQUALITY_NEG = 0
FACTOR_EQUALITY_POS = 1
FACTOR_CROSSING = 2
WIRE_PHYSICAL = 0
WIRE_IDENTITY = 1

# Hard diagrams may benefit enormously from a different elimination start, but
# multi-start search must not tax the small/structured cases whose historical
# greedy order is already excellent.
MULTISTART_MIN_PEAK_PORTS = 12
MULTISTART_MIN_IMPROVEMENT = 2
MULTISTART_MAX_STARTS = 16


def native_factorized_available() -> bool:
    return _yamada_factorized_frontier is not None


def factorized_import_error() -> Exception | None:
    return _FACTORIZED_IMPORT_ERROR


def require_native_factorized() -> None:
    """Fail clearly instead of silently selecting a superseded diagram backend."""
    if native_factorized_available():
        return
    detail = f" Original import error: {_FACTORIZED_IMPORT_ERROR!r}." if _FACTORIZED_IMPORT_ERROR else ""
    raise RuntimeError(
        "The optimized factorized Yamada extension is not available. Rebuild the "
        "current checkout in the active environment, for example with "
        "`python -m pip install -e '.[benchmark]' --force-reinstall`." + detail
    )


def _factor_order_plan(
    factor_order,
    factor_ports,
    port_factor,
    wire_partner,
):
    """Return exact live-port metrics for one factorized elimination order."""
    processed = bytearray(len(factor_ports))
    active: list[int] = []
    peak_ports = 0
    max_boundary_ports = 0
    boundary_area = 0

    for factor in factor_order:
        active.extend(factor_ports[factor])
        peak_ports = max(peak_ports, len(active))
        processed[factor] = 1
        active = [
            port
            for port in active
            if not processed[port_factor[wire_partner[port]]]
        ]
        boundary = len(active)
        max_boundary_ports = max(max_boundary_ports, boundary)
        boundary_area += boundary

    if active:
        raise RuntimeError("factorized frontier planner did not close")
    return {
        "peak_ports": int(peak_ports),
        "max_boundary_ports": int(max_boundary_ports),
        "boundary_area": int(boundary_area),
    }


def _multistart_factor_order(
    adjacency,
    factor_ports,
    port_factor,
    wire_partner,
    initial_order,
    *,
    crossing_count: int,
):
    """Conservatively rescue genuinely wide factorized production diagrams."""
    initial_plan = _factor_order_plan(
        initial_order,
        factor_ports,
        port_factor,
        wire_partner,
    )
    count = len(factor_ports)
    if (
        count <= 1
        or int(crossing_count) < 12
        or initial_plan["peak_ports"] < MULTISTART_MIN_PEAK_PORTS
    ):
        return list(initial_order), initial_plan, False, 1

    weighted_degree = [sum(neighbors.values()) for neighbors in adjacency]
    ranked = sorted(
        range(count),
        key=lambda node: (
            weighted_degree[node],
            len(factor_ports[node]),
            node,
        ),
    )
    if count <= MULTISTART_MAX_STARTS:
        starts = ranked
    else:
        positions = {
            round(
                index * (count - 1) / (MULTISTART_MAX_STARTS - 1)
            )
            for index in range(MULTISTART_MAX_STARTS)
        }
        starts = [ranked[position] for position in sorted(positions)]

    initial_first = int(initial_order[0])
    if initial_first not in starts:
        starts.insert(0, initial_first)

    best_order = list(initial_order)
    best_plan = initial_plan
    best_key = (
        initial_plan["peak_ports"],
        initial_plan["max_boundary_ports"],
        initial_plan["boundary_area"],
        tuple(initial_order),
    )
    candidates = 0
    for first in starts:
        candidates += 1
        if first == initial_first:
            order = list(initial_order)
        else:
            order = _greedy_factor_order_from_first(
                adjacency,
                factor_ports,
                first,
            )
        plan = _factor_order_plan(
            order,
            factor_ports,
            port_factor,
            wire_partner,
        )
        key = (
            plan["peak_ports"],
            plan["max_boundary_ports"],
            plan["boundary_area"],
            tuple(order),
        )
        if key < best_key:
            best_order = order
            best_plan = plan
            best_key = key

    if (
        initial_plan["peak_ports"] - best_plan["peak_ports"]
        < MULTISTART_MIN_IMPROVEMENT
    ):
        return list(initial_order), initial_plan, False, candidates
    return best_order, best_plan, True, candidates


def build_factorized_frontier(prepared):
    vertex_count = len(prepared.vertex_ids)
    crossing_count = len(prepared.crossing_ids)
    original_port_count = len(prepared.arc_partner)

    vertex_ports = [[] for _ in range(vertex_count)]
    crossing_ports = [[] for _ in range(crossing_count)]
    for port in range(original_port_count):
        fixed = int(prepared.fixed_terminal_index[port])
        crossing = int(prepared.crossing_for_port[port])
        if fixed >= 0:
            vertex_ports[fixed].append(port)
        elif crossing >= 0:
            crossing_ports[crossing].append(port)
        else:
            raise RuntimeError("prepared Yamada port belongs to no factor")

    def neighbor_key(port):
        partner = int(prepared.arc_partner[port])
        fixed = int(prepared.fixed_terminal_index[partner])
        if fixed >= 0:
            return (0, fixed, partner)
        return (1, int(prepared.crossing_for_port[partner]), partner)

    wire_partner = [int(value) for value in prepared.arc_partner]
    wire_type = [WIRE_PHYSICAL] * original_port_count
    port_factor = [-1] * original_port_count
    factor_types = []
    factor_ports = []

    for ports in vertex_ports:
        ordered = sorted(ports, key=neighbor_key)
        if not ordered:
            factor_types.append(FACTOR_EQUALITY_NEG)
            factor_ports.append([])
            continue

        previous_right = None
        for segment_index, original_port in enumerate(ordered):
            factor = len(factor_types)
            local_ports = [original_port]
            port_factor[original_port] = factor

            if previous_right is not None:
                left = len(wire_partner)
                wire_partner.append(previous_right)
                wire_type.append(WIRE_IDENTITY)
                wire_partner[previous_right] = left
                wire_type[previous_right] = WIRE_IDENTITY
                port_factor.append(factor)
                local_ports.append(left)

            if segment_index + 1 < len(ordered):
                right = len(wire_partner)
                wire_partner.append(-1)
                wire_type.append(WIRE_IDENTITY)
                port_factor.append(factor)
                local_ports.append(right)
                previous_right = right
            else:
                previous_right = None

            factor_types.append(
                FACTOR_EQUALITY_NEG if segment_index == 0 else FACTOR_EQUALITY_POS
            )
            factor_ports.append(local_ports)

    crossing_factor_by_index = []
    for ports in crossing_ports:
        if len(ports) != 4:
            raise RuntimeError("prepared crossing must have exactly four ports")
        factor = len(factor_types)
        factor_types.append(FACTOR_CROSSING)
        local_ports = sorted(ports)
        factor_ports.append(local_ports)
        crossing_factor_by_index.append(factor)
        for port in local_ports:
            port_factor[port] = factor

    if any(value < 0 for value in port_factor):
        raise RuntimeError("factorized Yamada port has no owner")
    if any(value < 0 for value in wire_partner):
        raise RuntimeError("unpaired factorization identity wire")

    adjacency = [defaultdict(int) for _ in factor_types]
    for port, partner in enumerate(wire_partner):
        if port >= partner:
            continue
        left = port_factor[port]
        right = port_factor[partner]
        if left != right:
            adjacency[left][right] += 1
            adjacency[right][left] += 1
    initial_factor_order = _greedy_factor_order(adjacency, factor_ports)
    factor_order, order_plan, ordering_multistart, ordering_candidates = (
        _multistart_factor_order(
            adjacency,
            factor_ports,
            port_factor,
            wire_partner,
            initial_factor_order,
            crossing_count=crossing_count,
        )
    )
    initial_order_plan = _factor_order_plan(
        initial_factor_order,
        factor_ports,
        port_factor,
        wire_partner,
    )

    plus_partner = [-1] * len(wire_partner)
    minus_partner = [-1] * len(wire_partner)
    for ports in prepared.ordered_ports:
        for port in ports:
            plus_partner[port] = int(prepared.plus_partner[port])
            minus_partner[port] = int(prepared.minus_partner[port])

    return {
        "factor_types": tuple(factor_types),
        "port_factor": tuple(port_factor),
        "wire_partner": tuple(wire_partner),
        "wire_type": tuple(wire_type),
        "plus_partner": tuple(plus_partner),
        "minus_partner": tuple(minus_partner),
        "factor_order": tuple(factor_order),
        "factor_order_peak_ports": int(order_plan["peak_ports"]),
        "factor_order_max_boundary_ports": int(
            order_plan["max_boundary_ports"]
        ),
        "factor_order_initial_peak_ports": int(
            initial_order_plan["peak_ports"]
        ),
        "factor_order_multistart": bool(ordering_multistart),
        "factor_order_candidates": int(ordering_candidates),
        "original_port_count": original_port_count,
        "crossing_factor_by_index": tuple(crossing_factor_by_index),
    }


def _frontier_args(data):
    return (
        list(data["factor_types"]),
        list(data["port_factor"]),
        list(data["wire_partner"]),
        list(data["wire_type"]),
        list(data["plus_partner"]),
        list(data["minus_partner"]),
        list(data["factor_order"]),
    )


def compute_factorized_frontier_laurent(prepared, *, stats=None):
    """Evaluate one prepared Yamada diagram with the sole production DP."""
    require_native_factorized()
    data = build_factorized_frontier(prepared)
    args = _frontier_args(data)
    try:
        compute_int64 = getattr(
            _yamada_factorized_frontier,
            "compute_factorized_frontier_int64",
            None,
        )
        if compute_int64 is None:
            compute_int64 = _yamada_factorized_frontier.compute_factorized_frontier
        value = compute_int64(*args)
        backend = "native-int64"
        overflowed = False
    except OverflowError:
        compute_bigint = getattr(
            _yamada_factorized_frontier,
            "compute_factorized_frontier_bigint",
            None,
        )
        if compute_bigint is None:
            raise RuntimeError(
                "The optimized factorized Yamada extension overflowed but does "
                "not expose the native bigint fallback. Rebuild this checkout."
            ) from None
        value = compute_bigint(*args)
        backend = "native-bigint"
        overflowed = True
    if stats is not None:
        stats["coefficient_backend"] = backend
        stats["int64_overflow"] = overflowed
        stats["factor_order_peak_ports"] = data["factor_order_peak_ports"]
        stats["factor_order_initial_peak_ports"] = data[
            "factor_order_initial_peak_ports"
        ]
        stats["factor_order_multistart"] = data["factor_order_multistart"]
        stats["factor_order_candidates"] = data["factor_order_candidates"]
    return tuple((int(power), int(coefficient)) for power, coefficient in value)
