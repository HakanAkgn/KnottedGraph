# Quick Start volume

`quickstart-handlebody.npz` contains a 96 × 96 × 96 boolean `volume`, its
world-coordinate `origin` and uniform `spacing`, and the tube `radius`.
Voxel index `(i, j, k)` maps to `origin + spacing * [i, j, k]`.
The volume surrounds a trefoil with an additional connecting arc, forming
a handlebody with two independent cycles.

The construction is the trefoil-theta control used in the
[handlebody tutorial](../../User_guide/benchmarks/04_thick_handlebody_validation.ipynb),
sampled here on a 96³ grid with radius `0.2564841814517642`.
Run the volume-to-graph example with:

```bash
uv run --extra knot-fields python examples/volume_to_yamada.py
```

To rebuild the volume from its edge curves, run:

```bash
uv run python dev/build_quickstart_volume.py
```

The displayed surface is the `volume=0.5` isosurface. The example extracts
the graph directly from the boolean volume and computes its PD code and
normalized Yamada polynomial in voxel coordinates.
