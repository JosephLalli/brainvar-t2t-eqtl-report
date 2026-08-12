# BrainVar GRCh38-vs-T2T cis-eQTL report

A long-form write-up of the four-arm GRCh38 versus T2T-CHM13 cis-eQTL analysis of the
BrainVar developmental cohort: the missing-genotype convention mismatch that invalidated
the first pass, the rebuild and its audit, and the findings that survived.

## Layout

- `index.html` — the report. Self-contained apart from the figures; theme-aware.
- `figures/` — eight figures, each rendered in a light and a dark variant.
- `build_figures.py` — renders every figure from the collected dataset.
- `build_site.py` — assembles `index.html`, injecting all numbers and tables from the
  same dataset so the prose cannot drift from the run tree.
- `.nojekyll` — serve the files as-is, without Jekyll processing.

## Rebuilding

Both scripts read `report_data.json`, collected from the corrected run roots by
`collect_report_data.py` in the analysis workspace.

```
python3 build_figures.py     # writes figures/*.{light,dark}.svg
python3 build_site.py        # writes index.html
```

No number on the page is transcribed by hand.
