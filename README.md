# BrainVar GRCh38-vs-T2T cis-eQTL report

A long-form write-up of the four-arm GRCh38 versus T2T-CHM13 cis-eQTL analysis of the
BrainVar developmental cohort: two references crossed with a linear and a pangenome-graph
aligner, what each factor changes, and where on the genome the two references give
measurably different answers.

## Start here

- `docs/ANALYSIS_PLAN.md` — the question the project answers, the design, and every
  workstream with its status. Read before starting any work.
- `docs/STATE.md` — what is done, what is next, and the known risks. Update it in the same
  session as the work.
- `docs/CONVENTIONS.md` — run-root structure, statistical rules, and the figure and prose
  conventions, which follow Lalli et al. 2025.

## Layout

- `index.html` — the report. Self-contained apart from the figures; theme-aware.
- `figures/` — nine figures, each rendered in a light and a dark variant.
- `build_figures.py` — renders every figure from the collected dataset.
- `build_site.py` — assembles `index.html`, injecting all numbers and tables from the
  same dataset so the prose cannot drift from the run tree.
- `.nojekyll` — serve the files as-is, without Jekyll processing.

## Rebuilding

Both scripts read `report_data.json`, collected from the run roots by
`collect_report_data.py` in the analysis workspace.

```
python3 build_figures.py     # writes figures/*.{light,dark}.svg
python3 build_site.py        # writes index.html
```

No number on the page is transcribed by hand.
