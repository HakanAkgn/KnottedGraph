
import numpy as np
import networkx as nx
from shapely.ops import substring
from shapely.geometry import Point, LineString
import re
import hashlib
from collections import defaultdict
from typing import Sequence
from numpy.typing import NDArray


__all__ = [
    "get_rotation_matrix",
    "generate_isotopy_angles",
    "cut_line_string",
    "parse_pd_code_string",
    "build_state_graph",
    "multigraph_key",
]


def _validate_rotation_order(order: str) -> str:
    """Return a validated Euler-axis order used by projection helpers."""

    if not isinstance(order, str):
        raise TypeError(
            "rotation_order must be a string such as 'ZYX' (intrinsic) "
            "or 'xyz' (extrinsic)."
        )
    if len(order) != 3 or any(axis.lower() not in "xyz" for axis in order):
        raise ValueError(
            "rotation_order must be a three-character Euler-axis sequence "
            "containing only x, y, and z, for example 'ZYX' or 'xyz'."
        )
    if not (order.isupper() or order.islower()):
        raise ValueError(
            "rotation_order must use one case consistently: uppercase for "
            "intrinsic rotations (for example 'ZYX') or lowercase for "
            "extrinsic rotations (for example 'xyz')."
        )
    return order


def _validate_positive_integer(value: int, *, name: str) -> int:
    """Return *value* as an ``int`` after strict positive-integer validation."""

    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise TypeError(f"{name} must be a positive integer, not {type(value).__name__}.")
    if value <= 0:
        raise ValueError(f"{name} must be a positive integer.")
    return int(value)


def get_rotation_matrix(
    angles: Sequence[float],
    order: str = "xyz",
    use_radians: bool = False,
) -> NDArray:
    """
    Construct a 3‑D rotation matrix from an ordered triple of elementary rotations.

    Parameters
    ----------
    angles : (3,) sequence
        The three rotation angles (α, β, γ) applied in the specified order.
    order : str, default "xyz"
        Three‑letter string giving the rotation axes.
        • Lower‑case  → extrinsic (rotations about *world* axes)
        • Upper‑case  → intrinsic (rotations about *body‑fixed* axes)
        Accepted letters: x/X, y/Y, z/Z.
    use_radians : bool, default False
        Supply angles in radians if True, else in degrees.

    Returns
    -------
    R : (3, 3) ndarray
        Combined rotation matrix.
    
    Examples
    --------
    1) Classic aerospace yaw‑pitch‑roll (intrinsic Z‑Y‑X)
    >>> yaw, pitch, roll = 30, 15, 5    # degrees
    >>> R = get_rotation_matrix((yaw, pitch, roll), order="ZYX")

    2) Extrinsic roll‑pitch‑yaw for camera pose (x‑y‑z in world frame)
    >>> R_cam = get_rotation_matrix((5, 15, 30), order="xyz")

    3) Euler ZXZ (commonly used in molecular crystallography)
    >>> phi, theta, psi = 45, 60, 10
    >>> R_euler = get_rotation_matrix((phi, theta, psi), order="zxz")
    """
    order = _validate_rotation_order(order)

    # Convert to radians if necessary
    a, b, c = angles if use_radians else np.radians(angles)

    def _single_axis(axis: str, theta: float) -> NDArray:
        """Elementary rotation about world axes (lowercase) or body axes (uppercase)."""
        ct, st = np.cos(theta), np.sin(theta)
        if axis.lower() == "x":
            R = np.array([[1, 0, 0],
                          [0, ct, -st],
                          [0, st,  ct]])
        elif axis.lower() == "y":
            R = np.array([[ ct, 0, st],
                          [  0, 1, 0],
                          [-st, 0, ct]])
        else:  # 'z'
            R = np.array([[ct, -st, 0],
                          [st,  ct, 0],
                          [ 0,   0, 1]])
        return R

    # Build the three elementary matrices
    thetas = [a, b, c]
    R_elems = [_single_axis(ax.lower(), th) for ax, th in zip(order, thetas)]

    # Combine them: intrinsic (upper‑case) postmultiplies;
    # extrinsic (lower‑case) premultiplies
    R = np.eye(3)
    for ax, Rk in zip(order, R_elems):
        R = R @ Rk if ax.isupper() else Rk @ R

    return R


def _matrix_to_euler(matrix: NDArray, order: str) -> NDArray:
    """Decompose a rotation using the same convention as get_rotation_matrix.

    Reduce intrinsic sequences to reversed extrinsic sequences.  The usual
    three-axis and repeated-axis decompositions then avoid a SciPy dependency.
    """
    axes = order.lower()[::-1] if order.isupper() else order
    if axes[0] == axes[1] or axes[1] == axes[2]:
        raise ValueError("Sampling requires an Euler sequence with distinct adjacent axes.")
    i = "xyz".index(axes[0])
    parity = int("xyz".index(axes[1]) != (i + 1) % 3)
    j = (i + 1 + parity) % 3
    k = (i + 2 - parity) % 3
    repeated = axes[0] == axes[2]
    m = matrix
    if repeated:
        s = float(np.hypot(m[i, j], m[i, k]))
        if s > 16 * np.finfo(float).eps:
            angles = np.array([np.arctan2(m[i, j], m[i, k]), np.arctan2(s, m[i, i]), np.arctan2(m[j, i], -m[k, i])])
        else:
            angles = np.array([np.arctan2(-m[j, k], m[j, j]), np.arctan2(s, m[i, i]), 0.0])
    else:
        c = float(np.hypot(m[i, i], m[j, i]))
        if c > 16 * np.finfo(float).eps:
            angles = np.array([np.arctan2(m[k, j], m[k, k]), np.arctan2(-m[k, i], c), np.arctan2(m[j, i], m[i, i])])
        else:
            angles = np.array([np.arctan2(-m[j, k], m[j, j]), np.arctan2(-m[k, i], c), 0.0])
    if parity:
        angles = -angles
    return angles[::-1] if order.isupper() else angles


def generate_isotopy_angles(
    N: int,
    order: str = "ZYX",
    use_radians: bool = False,
) -> NDArray:
    """Return deterministic views distributed over the upper hemisphere.

    The last row of each resulting rotation matrix is the intended viewing
    direction in input coordinates.  In-plane spin is fixed by the camera
    basis; setting an arbitrary Euler angle to zero does not in general
    remove that spin.  Uppercase sequences are intrinsic and lowercase
    sequences extrinsic, as in :func:`get_rotation_matrix`.
    """
    N = _validate_positive_integer(N, name="N")
    order = _validate_rotation_order(order)
    i = np.arange(N)
    phi = i * np.pi * (3 - np.sqrt(5))
    z = (i + 0.5) / N
    r = np.sqrt(1.0 - z**2)
    angles = []
    for azimuth, height, radial in zip(phi, z, r):
        cp, sp = np.cos(azimuth), np.sin(azimuth)
        # Orthonormal right/up/view rows, with determinant +1.
        matrix = np.array([[height * cp, height * sp, -radial],
                           [-sp, cp, 0.0],
                           [radial * cp, radial * sp, height]])
        angles.append(_matrix_to_euler(matrix, order))
    result = np.asarray(angles)
    return result if use_radians else np.degrees(result)


def cut_line_string(line, distances, *, tol=1e-12):
    """
    Split a LineString at one or more distances measured from its start.

    Parameters
    ----------
    line : shapely.geometry.LineString
    distances : float | Iterable[float]
        A single distance or an iterable of distances along the line
        (same units as line.length).  Values outside (0, line.length)
        or closer than *tol* to a previous split are ignored.
    tol : float, optional
        Numerical tolerance when comparing distances (default 1e‑12).

    Returns
    -------
    list[LineString]
        Ordered sub‑segments whose concatenation equals *line*.
    """
    # ── 1. normalise the split positions ──────────────────────────────────
    if np.ndim(distances) == 0:
        distances = [float(distances)]
    else:
        distances = [float(d) for d in distances]

    L = line.length
    # keep only unique, in‑range breakpoints and sort them
    cuts = sorted(
        {d for d in distances if tol < d < L - tol},
    )
    if not cuts:
        return [LineString(line)]

    # ── 2. iterate through [0, d1], [d1, d2], …, [dk, L] ───────────────────
    segments = []
    start = 0.0
    for d in cuts + [L]:
        if d - start > tol:           # skip zero‑length chunks
            # Shapely ≥ 2.0: use fast substring
            try:
                seg = substring(line, start, d)
            except ImportError:
                # manual fallback: interpolate the two endpoints and rebuild
                p0 = line.interpolate(start)
                p1 = line.interpolate(d)
                # collect intermediate vertices between the two points
                coords = [p0.coords[0]]
                for x, y, *z in line.coords:
                    pd = line.project(Point(x, y))
                    if start < pd < d:
                        coords.append((x, y, *z))
                coords.append(p1.coords[0])
                seg = LineString(coords)
            segments.append(seg)
        start = d
    return segments


def parse_pd_code_string(pd_str):
    vertices, crossings = [], []
    # allow zero or more digits/commas inside the brackets
    tokpat = re.compile(r'^(V|X)\[\s*([\d,]*)\s*\]$')
    for raw in pd_str.strip().split(';'):
        token = raw.strip()
        if not token:
            continue

        m = tokpat.match(token)
        if not m:
            raise ValueError(f"Bad token: {token!r}")

        kind, nums = m.groups()
        # if nums is empty, we want an empty list
        labels = [int(n) for n in nums.split(',')] if nums else []

        if kind == 'V':
            vertices.append(labels)
        else:
            crossings.append(labels)

    return vertices, crossings


def build_state_graph(vertices, crossings, state):
    G = nx.MultiGraph()
    nV, nX = len(vertices), len(crossings)
    G.add_nodes_from(range(nV),   kind='V')
    G.add_nodes_from(range(nV, nV+nX), kind='X')
    label_ends = defaultdict(list)
    for vidx, v_lbls in enumerate(vertices):
        for lbl in v_lbls:
            label_ends[lbl].append(vidx)
    for xidx, (i0,i1,i2,i3) in enumerate(crossings):
        Xnode = nV + xidx
        for lbl in (i0,i1,i2,i3):
            label_ends[lbl].append(Xnode)
    for lbl, ends in label_ends.items():
        if len(ends) != 2:
            raise ValueError(f"Label {lbl} appears {len(ends)} times; expected 2")
        G.add_edge(*ends, label=lbl)
    for xidx, res in enumerate(state):
        Xnode = nV + xidx
        i0,i1,i2,i3 = crossings[xidx]
        nbr = {}
        for lbl in (i0,i1,i2,i3):
            a,b = label_ends[lbl]
            nbr[lbl] = b if a==Xnode else a
        if res == 0:
            G.add_edge(nbr[i0], nbr[i3])
            G.add_edge(nbr[i1], nbr[i2])
            G.remove_node(Xnode)
        elif res == 1:
            G.add_edge(nbr[i0], nbr[i2])
            G.add_edge(nbr[i1], nbr[i3])
            G.remove_node(Xnode)
        elif res == 2:
            G.nodes[Xnode]['kind'] = 'V'
        else:
            raise ValueError(f"Invalid state {res}")
    return G


def multigraph_key(G):
    """Canonical unlabeled multigraph key."""
    try:
        import igraph as ig
        def _key(G):
            # Map NetworkX node labels to 0..n-1
            idx = {n: i for i, n in enumerate(G.nodes())}
            # Build edge list (parallel edges appear multiple times)
            edges = [(idx[u], idx[v]) for u, v, _ in G.edges(keys=True)]
            # Build igraph (undirected)
            igG = ig.Graph(len(idx), edges=edges, directed=False)
            # Uniform color partition
            igG.vs['color'] = 0
            # Canonical permutation (Bliss)
            perm = igG.canonical_permutation(color=igG.vs['color'])
            canon = igG.permute_vertices(perm)
            return tuple(sorted(canon.get_edgelist()))
    
    except ImportError:
        # Fallback to NetworkX for environments without igraph
        def _key(G):
            sigs = {}
            deg = dict(G.degree())
            loop_mult = {n: 0 for n in G.nodes()}
            for u, v in G.edges():
                if u == v:
                    loop_mult[u] += 1
            for n in G.nodes():
                # multiset of neighbor degrees, counting multiplicities of edges
                neigh_multideg = []
                for nbr in G.neighbors(n):
                    mult = sum(1 for _ in G.get_edge_data(n, nbr).values())
                    neigh_multideg.extend([deg[nbr]] * mult)
                sigs[n] = (
                    deg[n],
                    loop_mult[n],
                    tuple(sorted(neigh_multideg))
                )
            nodes = sorted(G.nodes(), key=lambda n: (sigs[n], n))
            idx = {n:i for i,n in enumerate(nodes)}
            edges = []
            for u,v in G.edges():
                i,j = idx[u], idx[v]
                if i > j:
                    i, j = j, i
                edges.append((i,j))
            return tuple(sorted(edges))
    
    key = _key(G)
    return hashlib.sha256(repr(key).encode()).hexdigest()

