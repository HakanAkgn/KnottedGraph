# Periodic volume and section invariants

Use `knotted_graph.periodic` for band energies sampled on the Brillouin
three-torus. Array axes correspond to the three ordered primitive reciprocal
vectors. Omit the repeated endpoint on each axis. A threshold must differ
from every vertex energy. The result describes the periodic piecewise-linear
interpolant on the Freudenthal triangulation. To transfer it to a continuous
band, check regularity and the topology of the interpolation.

## An occupied region and its section

This dimensionless cubic model has a genus-one boundary that winds along
the first primitive vector. The section at y = 0 has two open trajectories
in the repeated zone.

~~~python
import numpy as np

from knotted_graph.periodic import analyze_periodic_sublevel, periodic_plane_contours

axis = np.arange(24) / 24
x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
field = 1 - np.cos(2*np.pi*y) + 2*(1 - np.cos(2*np.pi*z))
field += 0.1*(1 - np.cos(2*np.pi*x))
level = 0.4

volume = analyze_periodic_sublevel(field, level)
boundary, = volume.boundaries
assert boundary.genus == 1
assert boundary.winding_basis == ((1, 0, 0),)

sections = periodic_plane_contours(field[:, 0, :], level)
assert len(sections) == 2
assert all(abs(section.winding[0]) == 1 and section.winding[1] == 0
           for section in sections)
assert all(section.area is None for section in sections)
~~~

A boundary records its occupied and complementary region indices, genus
and integral winding basis. The volume also gives the Betti numbers of both
signs. Call `has_same_invariants` on two volumes to compare their decorated
incidence graphs modulo labeling. An integer roll along array axes represents
a shift of the periodic cut.

The section routine uses the two ordered axes of the input section. The
last point of a lifted contour is its starting point plus the winding
vector. For a closed lift, multiply `section.area` by the area of the
reciprocal-plane basis parallelogram to obtain a physical section area.

## Complement representations

Pass the Boolean mask of complementary vertices to
`periodic_region_presentations` in `knotted_graph.periodic.fundamental_group`.
It gives a fundamental-group presentation for each component of the PL
complement. Count homomorphisms into the six-element permutation group with
`count_representations` and `s3_table`.

~~~python
from knotted_graph.periodic.fundamental_group import (
    count_representations,
    periodic_region_presentations,
    s3_table,
)

complements = periodic_region_presentations(field > level)
assert len(complements) == 1
representations = count_representations(complements[0], s3_table())
assert representations > 0
~~~

The count is an exact invariant of the input PL complement. Different
complements can yield the same count. Large presentations may exceed the
assignment budget.
