# Editing the documentation

The [published website](https://hakanakgn.github.io/KnottedGraph/) is the main
reading entry point. Build a local preview when contributing changes to its
pages or examples.

## Build and preview

From the repository root:

```bash
uv run --group docs python -m sphinx -b html -W --keep-going doc doc/_build/html
uv run --group docs python tests/contracts/test_documentation_contracts.py --html-root doc/_build/html
uv run python -m http.server 8000 --directory doc/_build/html
```

Open `http://localhost:8000` in a browser. The build treats warnings as errors,
and the second command checks the generated links. For an offline build, set
`KNOTTED_GRAPH_DOCS_OFFLINE=1` before the Sphinx command. Scientific notebooks
are displayed without executing their cells.

## Write for a new user

Keep the README focused on what the library does, source installation, a
copyable 3D example, citations and direct links to figures and saved data.
Link readers to the hosted guides for details. Put contributor tasks such as
building the website on this page.

Explain input preparation, supported schemas and computational settings in
the relevant functionality guide and nearby Python comments. State required
steps plainly and give the reader a way forward. Preserve information needed
to interpret results accurately, including coordinates, curve closure,
projection and normalization choices.

When changing an example, run it as a new user would and update its expected
output everywhere it appears. Check the rendered page for readable code,
figures, spacing and links. The maintained Quick Start figure can be rebuilt
with `uv run python dev/render_quickstart.py`.
