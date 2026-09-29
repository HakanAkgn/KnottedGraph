# Input Handling

Use an input adapter when your starting point is external data rather than an
already constructed embedded graph. The public adapters normalize supported
sources into either:

- a result object containing an embedded `networkx.MultiGraph`; or
- for surface files, a result object containing `PyVista.PolyData`.

Choose a reader below, then follow the preparation notes for your source data.

<div class="kg-link-row">
  <a href="../feature_status.html">Feature-status matrix</a>
  <a href="workflow_overview.html">Continue through the workflow</a>
  <a href="../api/inputs.html">Inputs API</a>
  <a href="https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/01_getting_started.ipynb">Open Getting Started</a>
</div>

## Supported public routes

Import the public functions from `knotted_graph.inputs`:

| Starting data | Public call | Result | Preparation / next step |
| --- | --- | --- | --- |
| Array or CSV/DAT/JSON/NPY/TSV/TXT/XYZ ordered coordinates | `from_coordinate_chain` | `CoordinateInputResult` | Store one ordered curve as an `(N, 3)` array |
| PDB file or four-character RCSB ID | `from_pdb_backbone` | `PDBBackboneInputResult` | Select chain/model/atom explicitly |
| CIF/mmCIF file or RCSB ID | `from_mmcif_backbone` | `MMCIFBackboneInputResult` | Supports the documented RCSB-style `_atom_site` subset |
| GROMACS GRO snapshot | `from_gromacs_gro` | `PolymerInputResult` | Default coordinate conversion is nm to Å |
| LAMMPS dump | `from_lammps_dump` | `PolymerInputResult` | Reads the first frame and unscaled `x/y/z` columns |
| Paired node/edge CSV files | `from_spatial_graph_csv` | `SpatialGraphInputResult` | Two files; preserves parallel edges and optional polylines |
| OBJ/OFF/PLY/STL/VTK/VTP surface | `from_surface_mesh` | `SurfaceInputResult` | Install `surface`, inspect the mesh, then choose an extraction workflow |

Named knots, torus types, and Artin braid words are handled by
`KnotFunction`; see {doc}`../applications/analytic_knot_fields`. They are
analytic constructors rather than generic file parsers.

## The embedded-graph contract

High-level graph adapters return a dataclass with `.graph` and `.issues`.
The graph is an undirected `networkx.MultiGraph`:

- every node has `pos`, a finite NumPy array of shape `(3,)`;
- every edge has `pts`, a finite array of shape `(N, 3)`;
- the first and last `pts` rows match the positions of the edge endpoints;
- parallel geometric edges remain separate MultiGraph edges; and
- source identifiers, selections, and caller metadata remain available on the
  result or graph metadata.

Always inspect issues before continuing:

```python
result = from_coordinate_chain(
    [[0, 0, 0], [1, 0, 0], [1, 1, 0]],
    input_id="demo-curve",
)

if result.issues:
    for issue in result.issues:
        print("input issue:", issue)

graph = result.graph
print(graph.number_of_nodes(), graph.number_of_edges())
```

Schema or selection errors raise an exception. Review the `.issues` list for
row-level or geometry details to address before continuing to projection.

## 1. Ordered coordinate chains

```python
import numpy as np
from knotted_graph.inputs import from_coordinate_chain

coordinates = np.array(
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.2],
        [1.0, 1.0, 0.0],
    ]
)

result = from_coordinate_chain(coordinates, input_id="open-curve")
graph = result.graph
```

For file input, pass a path. CSV uses named `x`, `y`, and `z` columns by
default; TSV/DAT/TXT read the first three fields; conventional XYZ files may
contain the atom count and a blank or nonblank comment line.

### Closure is explicit

- The default is an open chain.
- `closure="direct"` appends a straight final-to-first segment and produces a
  geometrically closed edge.
- `closed=True` without a closure method is valid only when the supplied first
  and last samples already agree.
- `closure="metadata_only"` records intent but does not close the geometry.

The result's `.coords` retain the original source samples. Inspect edge `pts`
when you need the geometry after direct closure.

## 2. PDB and mmCIF backbones

```python
from knotted_graph.inputs import (
    from_mmcif_backbone,
    from_protein_ca_backbone,
)

pdb_result = from_protein_ca_backbone(
    "protein.pdb",
    chain_id="A",
    model_id=1,
)

cif_result = from_mmcif_backbone(
    "structure.cif",
    chain_id="A",
    atom_name="CA",
    model_id=1,
)
```

A four-character PDB ID may be used instead of a local path when downloading
is enabled. For reproducible offline work, keep a local file and record its
provenance.

When multiple chains match, choose `chain_id` explicitly. The mmCIF reader
currently targets RCSB-style atom-site loops with one complete data row per
physical line. Convert other CIF layouts to that form before loading. PDB and
mmCIF backbone extraction creates an ordered curve from the selected atoms.
To study a protein interaction network, supply the desired graph connections
through the node/edge input route.

## 3. Polymer snapshots

```python
from knotted_graph.inputs import from_gromacs_gro, from_lammps_dump

gro = from_gromacs_gro(
    "chain.gro",
    atom_name="BB",
    output_unit_scale=10.0,
)

lammps = from_lammps_dump(
    "chain.dump",
    molecule_id=7,
)
```

The GRO adapter treats the selected atoms as one ordered coordinate curve;
its default scale converts nanometres to ångströms. The LAMMPS adapter reads
the first frame with unscaled `x/y/z` columns. Before loading, select and order
the atoms of the intended chain, unwrap periodic images when needed, and
convert scaled `xs/ys/zs` coordinates to `x/y/z`.

## 4. Paired spatial-graph CSV

The nodes file requires an identifier and coordinates:

```text
node_id,x,y,z
0,0,0,0
1,1,0,0
```

The edges file requires source and target identifiers:

```text
edge_id,source,target,points_json
e0,0,1,"[[0,0,0],[0.5,0.2,0],[1,0,0]]"
```

Load both files together:

```python
from knotted_graph.inputs import from_spatial_graph_csv

result = from_spatial_graph_csv("nodes.csv", "edges.csv")
graph = result.graph
```

`points_json` is optional; without it the edge is straight. When present, the
polyline endpoints must match the corresponding node positions. Extra CSV
columns are preserved as string attributes. For GraphML, SWC, edge-list or
graph-JSON sources, prepare these two tables as described below.

## 5. Surface meshes

```python
from knotted_graph.inputs import from_surface_mesh

result = from_surface_mesh("surface.ply")
mesh = result.mesh
```

Install the optional dependency first:

```bash
uv sync --extra surface
```

Cleaning and triangulation are enabled by default and may change mesh
connectivity. The result contains a `PyVista.PolyData` mesh and reports open
boundaries in `.issues`. Inspect the mesh, choose any boundary treatment, then
follow {doc}`workflow_overview` to select the extraction and cleanup steps
appropriate to your surface.

## Preparing other file formats

The readers above accept the documented schemas. To use another format, first
read it with a tool for that format and map its data to the appropriate input:

| Source data | Preparation | Continue with |
| --- | --- | --- |
| GraphML, SWC, edge lists or graph JSON | Write node IDs and `x/y/z` positions to the node CSV, and connections to the edge CSV. Include curved edge geometry in `points_json`. | `from_spatial_graph_csv` |
| JSON or NPY storing a single curve | Arrange the samples as an ordered `(N, 3)` coordinate array. | `from_coordinate_chain` |
| NPZ scalar/vector volumes or oriented flows | Load the arrays and define their coordinate grid, then select the relevant sampling, tracing or extraction workflow. | {doc}`workflow_overview` |
| Hamiltonian files | Read the model into the in-memory representation used by the chosen material or nodal example. | {doc}`../applications/index` |

For SWC data, each point's ID, parent ID and coordinates supply the node and
edge rows; omit the parent edge for a root point. Preserve the original
positions and any sampled edge curves during
conversion: these describe the spatial embedding whose topology is analysed.
GraphML attribute names and graph-JSON layouts vary, so map their coordinate
fields explicitly to the schema above.

For the complete reader and application reference, see {doc}`../feature_status`.

## Continue through the pipeline

After loading a graph:

```python
from knotted_graph.core import ensure_embedding
from knotted_graph.projection import select_projection

graph = ensure_embedding(result.graph)
projection = select_projection(graph)
print("crossings:", projection.num_crossings)
```

Continue with {doc}`workflow_overview` for cleanup and extraction decisions or
{doc}`projection_yamada` for projection, PD-code, and invariant interpretation.
For exact signatures and return fields, use {doc}`../api/inputs`.
