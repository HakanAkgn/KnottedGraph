# Protein Derived Spatial Graphs

<div class="kg-hero">
  <p class="kg-lead">Load protein coordinates, construct an embedded graph, optionally relax its geometry, then project it and compute an invariant. Follow the protein applications notebook for input examples and the figure below for the geometric workflow.</p>
  <div class="kg-link-row">
    <a href="https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/applications/03_protein_applications.ipynb">Open 03_protein_applications.ipynb</a>
    <a href="../user_guide/input_adapters.html">PDB/mmCIF input guide</a>
  </div>
</div>

<div class="kg-wide-figure">
  <img src="../site_figures/repulsive_curves.png" alt="Repulsive curves workflow for embedded spatial graphs">
</div>

<a id="what-is-implemented-now"></a>

## Load an ordered backbone

The public input layer can extract an **ordered backbone trace** from a local
PDB/mmCIF file or an RCSB identifier:

```python
from knotted_graph.inputs import from_protein_ca_backbone

result = from_protein_ca_backbone(
    "1CRN",
    chain_id="A",
    model_id=1,
)

print(result.coords.shape)
print(result.graph.number_of_nodes(), result.graph.number_of_edges())
print(result.issues)
```

For nucleic acids use the corresponding backbone helper or select the desired
atom explicitly. Select `chain_id` when the structure has multiple matching
chains. Remote IDs
need network access on the first download, while local files work offline.

The result represents the sampled trace as one geometric curve edge (or a
self-loop when explicitly closed), with the selected atoms stored as samples
in edge `pts`. Read {doc}`../user_guide/input_adapters` for atom selection,
closure, metadata and the supported mmCIF atom-site layout.

<a id="what-is-not-yet-a-generic-workflow"></a>

## Build a protein-derived network

Choose the representation for your question: an ordered backbone, contact
graph, residue-interaction network, cavity skeleton or domain graph. Supply
the corresponding node and edge connections, then use the shared spatial-graph
workflow.

Before creating a derived spatial graph, document:

- which atoms/residues and model/chain were selected;
- the distance, contact, domain, or cavity rule that creates topology;
- coordinate units and any periodic/closure treatment;
- what each node and edge represents scientifically; and
- validation and perturbation checks for the derived graph.

With that mapping defined, continue with cleanup, optional repulsive
relaxation, projection and invariant computation. To use native Repulsor
relaxation, follow the {doc}`../user_guide/repulsive_layout` setup after
loading the backbone.

## Recommended next step

If an ordered backbone is the intended object, continue directly with
{doc}`../user_guide/workflow_overview`. If you need a domain-specific contact or
cavity graph, define and validate its construction rule, then continue with
graph inspection, projection and invariant calculation.
