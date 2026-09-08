#!/usr/bin/env python3
"""Publish preserved research demos with on-demand, content-addressed geometry.

Run on a compute node: the upstream TPMS HTML alone is about 82 MB. This does
not regenerate phase labels, smooth a grid, simplify a mesh or run extraction.
It separates each original region object into an unchanged JSON asset and
changes only the HTML loading controller. Source hashes are deliberately pinned.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "User_guide" / "applications" / "NewPhaseMapPlots"
SOURCES = {
    "materials": (
        "RealMaterials/html/07_hamiltonian_yamada_plotly_region_geometry_with_materials.html",
        "68f7153ea995d0f46df6986ed129657eced512d56f5e4b27a037344f8e6fcc24",
        "RealMaterials/figures/07_hamiltonian_yamada_material_phase_maps_3panel_like_porous.png",
    ),
    "tpms": (
        "TPMS/html/tpms_compact_c0_03_dense_stable_yamada_plotly_region_geometry.html",
        "07c256cf9ab2d47d27b31c41b5d0d256e4089356ff110a5f203c79cbeee2ebe2",
        "TPMS/figures/tpms_compact_scaffold_phase_maps_overview.png",
    ),
}

LOADER = """
let regionSelection = 0;
const loadedRegions = new Map();
function fitSmallScreen(id) {
  const plot = document.getElementById(id);
  const narrow = plot.clientWidth < 520;
  const current = String(plot.layout.title.text);
  if (!current.includes("<br>")) plot.dataset.fullTitle = current;
  const title = plot.dataset.fullTitle || current;
  const text = narrow ? title.replace(": ", "<br>").replace(", region", "<br>region") : title;
  return Plotly.relayout(plot, {"title.text": text,
    "title.font.size": narrow ? 14 : (id === "phaseMap" ? 24 : 21),
    "margin.t": narrow ? 80 : (id === "phaseMap" ? 68 : 50)});
}
async function loadRegion(regionKey) {
  if (loadedRegions.has(regionKey)) return loadedRegions.get(regionKey);
  const entry = regions[regionKey];
  if (!entry) return null;
  const response = await fetch(new URL(entry.asset, document.baseURI));
  if (!response.ok) throw new Error(`HTTP ${response.status} loading region`);
  const region = await response.json();
  loadedRegions.set(regionKey, region);
  return region;
}
"""
SELECTOR = """async function selectRegion(regionKey) {
  const request = ++regionSelection;
  const status = document.getElementById("status");
  status.textContent = "Loading the selected region…";
  let region;
  try {
    region = await loadRegion(regionKey);
  } catch (error) {
    if (request === regionSelection) {
      status.textContent = "Could not load the region. Serve this folder over HTTP " +
        "(not file://), check your connection, then select the region again. " + error.message;
    }
    return;
  }
  if (request !== regionSelection) return;
"""


def extract_payload(html: str) -> tuple[dict, str, str]:
    marker = "const payload = "
    offset = html.index(marker) + len(marker)
    payload, length = json.JSONDecoder().raw_decode(html[offset:])
    return payload, html[:offset], html[offset + length :]


def split_demo(source: Path, destination: Path, *, expected_hash: str) -> dict:
    """Split one known artifact, retaining every region object and non-region field."""
    raw = source.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != expected_hash:
        raise ValueError(
            f"Source changed: {source}; review it before updating the pinned hash."
        )
    payload, prefix, suffix = extract_payload(raw.decode("utf-8"))
    if not isinstance(payload.get("regions"), dict) or not payload["regions"]:
        raise ValueError("Expected a non-empty region dictionary.")
    old = "function selectRegion(regionKey) {\n  const region = regions[regionKey];"
    if suffix.count(old) != 1:
        raise ValueError(
            "Unrecognized region-selection controller; review before publishing."
        )
    asset_dir = destination / "regions"
    asset_dir.mkdir(parents=True, exist_ok=True)
    manifest = {"source_sha256": digest, "source_bytes": len(raw), "regions": {}}
    for key, region in payload["regions"].items():
        data = (
            json.dumps(
                region, ensure_ascii=False, sort_keys=True, separators=(",", ":")
            )
            + "\n"
        ).encode()
        region_hash = hashlib.sha256(data).hexdigest()
        name = f"regions/{region_hash}.json"
        (destination / name).write_bytes(data)
        manifest["regions"][key] = {
            "asset": name,
            "sha256": region_hash,
            "bytes": len(data),
        }
    payload["regions"] = {
        key: {"asset": item["asset"]} for key, item in manifest["regions"].items()
    }
    suffix = suffix.replace(old, SELECTOR, 1)
    suffix = suffix.replace(
        "let activeTransition = transitions[0].key;",
        'let activeTransition = transitions.find(item => item.key === "tib2_d6_F")?.key || transitions[0].key;',
        1,
    )
    suffix = suffix.replace(
        "const regions = payload.regions;",
        "const regions = payload.regions;\n" + LOADER,
        1,
    )
    suffix = suffix.replace(
        'window.addEventListener("resize", () => {',
        'window.addEventListener("resize", () => {\n  if (typeof Plotly === "undefined") return;',
        1,
    )
    for chart, data in (
        ("phaseMap", "[trace, boundaryTrace]"),
        ("geometry", "geometryTraces(region)"),
    ):
        call = f'Plotly.react("{chart}", {data}, layout, {{ responsive: true, displaylogo: false }});'
        if call not in suffix:
            raise ValueError(f"Unrecognized {chart} render call.")
        suffix = suffix.replace(
            call, call[:-1] + f'.then(() => fitSmallScreen("{chart}"));', 1
        )
        resize = f'Plotly.Plots.resize(document.getElementById("{chart}"));'
        suffix = suffix.replace(
            resize, resize[:-1] + f'.then(() => fitSmallScreen("{chart}"));', 1
        )
    prefix = prefix.replace(
        "</head>",
        """<style>
      .layout > * { min-width: 0; }
      #phaseMap, #geometry { max-width: 100%; }
      #status { overflow-wrap: anywhere; }
    </style>\n</head>""",
        1,
    )
    notice = (
        '<aside role="note" style="margin:12px 0;padding:12px;border:1px solid #555;font:16px/1.5 sans-serif">'
        "<strong>Saved research result.</strong> Colors include Yamada values and structural signatures; "
        "display smoothing/merges and contraction grouping are not proofs of polynomial equality. "
        "A representative geometry is loaded on selection. "
        '<a href="../../../applications/material_and_tpms_phase_maps.html">Read scope, provenance and usage</a>.'
        "</aside>"
    )
    prefix = prefix.replace("<body>", "<body>\n" + notice, 1)
    # JSON stays data inside a script element, never executable input.
    text = json.dumps(payload, ensure_ascii=False, separators=(",", ":")).replace(
        "</", "<\\/"
    )
    html = prefix + text + suffix
    bootstrap = "buildButtons();\nplotPhaseMap(activeTransition, activeMode);"
    if bootstrap not in html:
        raise ValueError("Unrecognized demo bootstrap.")
    html = html.replace(
        bootstrap,
        """if (typeof Plotly === "undefined") {
  document.getElementById("status").textContent = "Plotly could not load from its pinned CDN. " +
    "Check your connection or use the static figure on the documentation page.";
} else {
  buildButtons();
  plotPhaseMap(activeTransition, activeMode);
}""",
        1,
    )
    output = destination / "index.html"
    output.write_text(html, encoding="utf-8")
    manifest["index_bytes"] = output.stat().st_size
    manifest["index_sha256"] = hashlib.sha256(output.read_bytes()).hexdigest()
    return manifest


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "doc/assets/demos/new_phase_maps"
    )
    args = parser.parse_args(argv)
    manifest = {
        "upstream_commit": "2b2ae6d",
        "processing": "lossless JSON separation; no scientific recalculation",
        "demos": {},
    }
    for name, (relative, expected, preview) in SOURCES.items():
        target = args.output_dir / name
        result = split_demo(REFERENCE / relative, target, expected_hash=expected)
        result["source"] = "User_guide/applications/NewPhaseMapPlots/" + relative
        shutil.copyfile(REFERENCE / preview, target / "preview.png")
        result["preview_sha256"] = hashlib.sha256(
            (target / "preview.png").read_bytes()
        ).hexdigest()
        manifest["demos"][name] = result
        print(
            f"{name}: {result['source_bytes']:,} input bytes -> {result['index_bytes']:,} initial HTML bytes; "
            f"{len(result['regions'])} region attachments"
        )
    (args.output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
