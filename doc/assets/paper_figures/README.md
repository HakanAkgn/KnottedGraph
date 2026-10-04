# Published software-paper figures

The PDF files are unchanged author figures from
https://arxiv.org/src/2609.31152v1, with figure numbers assigned from that
version of the manuscript. `provenance.json` records the original source
paths, the archive digest and each PDF's SHA-256 digest.

The PNGs were rendered from the corresponding PDF with Poppler:

```bash
pdftoppm -png -singlefile -scale-to-x 1500 -scale-to-y -1 figure-3.pdf figure-3
```

These PNG previews are rendered from the author PDFs with Poppler. The figure
and data guide is `doc/paper_results.md`; its website entry is
`paper_results.html`.
