# Immutable historical handlebody inputs

This directory packages the earlier 5,000-case seed archive, exact recovered graphs, and result CSV for portable replay. It does not replace the repository's separate saved 100-case batch.

`manifest.json` records the original compressed archive SHA256, every shard's SHA256 and count, complete original non-case metadata (including configuration and generator state), and a content digest of the ordered JSON archive. Ten seed shards contain 500 cases each. The largest is 15,882,452 bytes, below the 50 MB packaging bound and GitHub's single-file limit. The recovered archive and CSV are byte-for-byte copies. Packing read all files back and verified equality of the complete reassembled JSON object with the original archive.

The data has **96 known incompatible seed/CSV associations**. They are listed in `seed_association_conflicts.json`; 94 differ structurally and two (`random_0027`, `random_0036`) have the same abstract metadata but different geometric certification. Historical rows lack exact seed hashes, so agreement of metadata is compatibility rather than an independent original-coordinate identity certificate. Do not relabel seeds by finding a graph with a desired polynomial. The corrected projection replay evaluates the exact provided graphs, records their new hashes, matches all 4,904 compatible pairs, and withholds the other 96 comparisons.

From the repository root:

```bash
uv run --no-sync python dev/shard_handlebody_inputs.py verify \
  --manifest User_guide/benchmarks/results/arxiv_revision/handlebody_inputs/manifest.json

uv run --no-sync python dev/audit_handlebody_recovery.py \
  --seed-archive User_guide/benchmarks/results/arxiv_revision/handlebody_inputs \
  --recovered-archive User_guide/benchmarks/results/arxiv_revision/handlebody_inputs/synthetic_ground_truth_recovered_graphs.json.gz \
  --csv User_guide/benchmarks/results/arxiv_revision/handlebody_inputs/synthetic_ground_truth_yamada_preservation.csv \
  --case all --workers 4 --output _build/handlebody_full_replay.jsonl
```

The runner refuses an existing output file. To inspect only the historical mismatch, use `--case random_4920` with a new output filename. A directory or its `manifest.json` can be supplied as `--seed-archive`.

Optional reassembly for other tools:

```bash
uv run --no-sync python dev/shard_handlebody_inputs.py reassemble \
  --manifest User_guide/benchmarks/results/arxiv_revision/handlebody_inputs/manifest.json \
  --output _build/certified_thick_handlebody_embeddings.json.gz
```

Reassembly preserves the JSON structure, numeric values, case ordering, configuration, and generator state. Its gzip header and whitespace differ from the historical file; use the declared content digest to verify it. No graph is regenerated.
