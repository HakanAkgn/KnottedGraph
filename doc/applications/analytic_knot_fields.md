# Analytic knot fields and topological deformations

<div class="kg-hero">
  <p class="kg-lead">Construct a named knot, torus knot/link, or Artin-braid closure as an analytic complex field; diagnose a sampled tubular level set; convert the accepted volume to a spatial graph; and optionally scan a deformation between two fields.</p>
  <div class="kg-link-row">
    <a href="https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/applications/04_analytic_knot_fields.ipynb">Open the application notebook</a>
    <a href="../feature_status.html">Check extras and scaling</a>
  </div>
</div>

KnottedGraph can construct complex scalar fields whose zero sets represent knots or links, form sampled tubular handlebodies from their sublevel sets, skeletonize those volumes with the package's canonical extractor, and analyze the resulting embedded spatial graph.

The generic route is


$$
\text{knot / braid}\to f(u,v,\bar v)\to F(x,y,z)
\to \{|F|\leq\epsilon\}\to\text{spatial graph}\to\text{topological analysis}.
$$


## Mathematical basis

Bode and Dennis give a constructive algorithm that starts from a braid and constructs a semiholomorphic polynomial, polynomial in \(u,v,\bar v\), whose zero set on the unit three-sphere realizes the closed braid for sufficiently small scaling. See Bode & Dennis, *J. Knot Theory Ramifications* 28, 1850082 (2019), DOI `10.1142/S0218216518500827`, and arXiv:1612.06328. Bode later gave a constructive proof that all link types in \(S^3\) arise from semiholomorphic polynomials; see *European Journal of Mathematics* 9, 85 (2023), DOI `10.1007/s40879-023-00678-1`.

The braid constructor uses a finite Fourier representation. Sample the resulting 3-D field and use the diagnostics below to inspect its tubular geometry.

`from_braid(...)` reports the sampled root-set error relative to the minimum strand separation for the finite Fourier realization. Use tubular diagnostics and a resolution comparison to assess the resulting 3-D level set.

## Construct fields

```python
from knotted_graph.inputs import KnotFunction

# Preferred exact/reference constructors where available.
trefoil = KnotFunction.from_name("3_1")
figure8 = KnotFunction.from_name("4_1")

# Universal route: closure of an Artin braid word.
figure8_from_braid = KnotFunction.from_braid(
    [1, -2, 1, -2],
    strands=3,
)

# Force a named example through the generic compiler.
figure8_from_braid_2 = KnotFunction.from_name(
    "4_1",
    construction="braid",
)

# Exact torus-knot/link fast path f(u,v)=u**p-v**q.
trefoil_exact = KnotFunction.torus(2, 3)
```

A positive integer `i` denotes the Artin generator \(\sigma_i\); `-i` denotes \(\sigma_i^{-1}\). Use the built-in catalogue for named examples and `from_braid(...)` to construct arbitrary braid closures.

The preferred `4_1` constructor uses Rudolph's explicit semiholomorphic figure-eight polynomial


$$
f=u^3-3v^2\bar v^2(1+v^2-\bar v^2)u-2(v^2+\bar v^2),
$$


as quoted in the semiholomorphic-link literature.

## From \(S^3\) to \(\mathbb R^3\)

The package uses


$$
u=\frac{2(x+iy)}{1+r^2},\qquad
v=\frac{r^2-1+2iz}{1+r^2},\qquad r^2=x^2+y^2+z^2,
$$


so \(|u|^2+|v|^2=1\). Generic braid-generated fields may apply a fixed unitary rotation of the \(S^3\) coordinates to move the link away from the omitted stereographic projection point. This ambient diffeomorphism of \(S^3\) preserves the link type.

A requested sublevel radius is rejected if it contains the projection pole, because that would become non-compact in the selected \(\mathbb R^3\) chart. With `require_compact=True`, boundary contact raises an error; enlarge the sampling span to contain the sublevel set.

## Handlebody and convergence checks

For a field \(F\), use


$$
H_\epsilon=\{|F|\leq\epsilon\},\qquad
\partial H_\epsilon=\{|F|=\epsilon\}.
$$


For sufficiently small regular \(\epsilon\), this is a tubular neighborhood of the link. Check its sampled geometry with:

```python
figure8 = KnotFunction.from_name("4_1")

diagnostic = figure8.diagnose_level(
    0.55,
    span=((-4, 4),) * 3,
    dimension=160,
)

report = figure8.tubular_convergence(
    0.55,
    span=((-4, 4),) * 3,
    dimensions=(128, 160),
)
assert report.converged
```

The diagnostics include connected volume components, boundary components, Euler characteristic, closed-boundary status, total genus, and sampling-box contact. For a \(c\)-component link represented by disjoint tubes, the expected result is \(c\) solid tori: \(c\) volume components, \(c\) boundary components, and total boundary genus \(c\).

The regression suite includes the non-torus figure-eight case and checks convergence to one connected solid torus with one closed genus-one boundary at 128 and 160 grid points per axis.

## Convert to a spatial graph

```python
graph = figure8.to_spatial_graph(
    radius=0.55,
    span=((-4, 4),) * 3,
    dimension=160,
)
```

This reuses KnottedGraph's optimized `volume -> Lee skeleton -> sparse MultiGraph` path. Node `pos` and edge `pts` coordinates are transformed from voxel indices back to the physical sampling coordinates before the graph is returned.

## Deform two knot fields

```python
import numpy as np
from knotted_graph.inputs import KnotFunction, KnotFunctionPath
from knotted_graph.applications.knot_deformation import KnotDeformationScan

start = KnotFunction.from_name("3_1")
end = KnotFunction.from_name("4_1")
path = KnotFunctionPath(start, end)

scan = KnotDeformationScan(
    path,
    lambdas=np.linspace(0, 1, 31),
    radii=np.linspace(0.1, 0.5, 21),
    span=((-4, 4),) * 3,
    dimension=128,
    invariant="yamada",
)
result = scan.run()
```

Before linear interpolation, the path RMS-normalizes both endpoint fields and aligns their global complex phase. Multiplying \(F\) by a nonzero scale or global phase does not change its zero set, but would otherwise change a naive interpolation.

The deformation is defined by the chosen field representatives and interpolation path, so its intermediate topology depends on those choices.

## Reusable Bloch phase-space calculation

For periodic Bloch models, use the nodal application layer:

```python
import numpy as np
from knotted_graph.applications.nodal import (
    NodalBlochPath,
    NodalPhaseScan,
    hopf_link_bloch_vector,
    unknot_bloch_vector,
)

path = NodalBlochPath(unknot_bloch_vector, hopf_link_bloch_vector)
scan = NodalPhaseScan(
    path,
    lambdas=np.linspace(0, 1, 41),
    gammas=np.linspace(0.1, 1.0, 30),
    dimension=96,
)
result = scan.run()
```

`NodalBlochPath.at_components(...)` additionally supports independent \((\lambda_x,\lambda_y,\lambda_z)\) mixing, matching the component-wise deformation used in the earlier notebook. Use a periodic Bloch model for Brillouin-zone calculations; the generic knot-field constructors describe fields on $S^3$ and $\mathbb R^3$.
