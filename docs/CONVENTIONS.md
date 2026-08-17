# Conventions

How to do the work so that any two contributors produce compatible output, and so the page
reads as one voice. Read alongside `ANALYSIS_PLAN.md`.

## Where things live

| | Path |
|---|---|
| Analysis workspace (not a git repo) | `/mnt/ssd/lalli/nf_stage/brainvar_t2t_neurodevelopment_gain_loss` |
| Run roots | `<workspace>/runs/<semantic_name>_<YYYYMMDD>` |
| Analysis scripts | `<workspace>/scratch/analyze_<semantic_name>_<YYYYMMDD>.py` |
| Logs | `<workspace>/logs/<semantic_name>_<YYYYMMDD>.log` |
| Collected data contract | `<workspace>/scratch/report_data.json` |
| Page repo (git) | `/mnt/ssd/lalli/nf_stage/brainvar_t2t_eqtl_report` |

Name runs and scripts for **what they do**, never for sequence or position. `directional_window_discordance`, not `analysis_3` or `phase_B`.

## Run roots are the unit of resumable work

Every run root must contain `MANIFEST.json` with at least:

```json
{
  "analysis": "semantic_name",
  "status": "complete",
  "source_run": "...",
  "<parameters>": "...",
  "inferential_boundary": "what this statistic does NOT establish",
  "results": {}
}
```

Include SHA-256 for every external track or table read. If a run consumed another run's
output, name that run root in `source_run`. A run without a manifest did not happen.

**Delete the manifest before re-running a run root.** The collector treats a manifest with
`"status": "complete"` as current, so a manifest left over from an earlier or narrower run
will be read as if it described the new one — silently putting stale numbers on the page.
This has already happened once: a pilot restricted to one chromosome was picked up as though
it were the genome-wide run. Clear the run root, or write to a new one.

## Statistics

- **Correct the background.** Enrichment is always tested against the entities actually
  tested in the same contrast, never against the genome. The tested set is strongly
  non-random and a genome background manufactures significance.
- **Control for variant count.** Window and gene statistics are standardised within
  variant-count strata (40 quantile bins, median and 1.4826×MAD) before thresholding. A
  window with 4,000 variants is not comparable to one with 120.
- **Choose thresholds by error rate, not by sigma.** Report the cut that holds the estimated
  false-discovery proportion at 5%, and show what a fixed 1.96 would have cost.
- **State when an estimate is a lower bound.** Normal-null FDR estimates assume a normal tail;
  LD correlates neighbouring windows and genes, so the true null tail is heavier. Say so.
- **Effect sizes over exponents.** Fisher tests treat genes as independent; genes in one
  duplicated block are not. Quote odds ratios as the trustworthy part.
- **No ground truth exists.** Never phrase a result as accuracy. See `ANALYSIS_PLAN.md`.

## Presentation, following Lalli et al. 2025

The 2025 paper's conventions, adopted deliberately so the two read as a series:

- **Quantify as percentage change and fold change**, not as raw differences. "38% fewer
  assembly-discordant genotypes", "ten-fold reduction" — a reader should get the magnitude
  from the sentence without consulting the figure.
- **Stratify everything.** By MAF bin (log scale), by variant type (SNP vs indel), and by
  region class (syntenic / non-syntenic, segmental duplication, blacklist, newly accessible).
  An unstratified genome-wide number hides the entire result.
- **Show dispersion, never bare means.** Box or violin with explicit interquartile range and
  5th–95th whiskers, or bars with confidence intervals. These distributions are heavily
  zero-inflated; a mean of a zero-inflated fraction implies a shift in the middle that is not
  there. Where the difference lives in the tail, quote the tail.
- **Rank regions by magnitude and pair them with disease relevance** (the paper's Table 2
  pattern). This is the format for the discordance-hotspot table.
- **Direct quantitative sentences.** "We find that alignment to T2T-CHM13 resulted in 38%
  fewer assembly-discordant genotypes and 16% fewer switch errors."

## Figures

- Rendered from one code path into paired `<name>.light.svg` and `<name>.dark.svg`; the page
  selects by `prefers-color-scheme`.
- Palette slots `s1`–`s4` are pre-validated for colour-vision separation and contrast. Do not
  introduce new colours.
- Every figure carries direct value labels and is mirrored by a data table on the page — two
  light-mode slots sit below 3:1 contrast and this is the required relief.
- Add a figure by appending a `fig_<name>(c, mode)` function and registering it in `FIGURES`.

## Prose

- Report a negative, then immediately report what *did* move. Never let "no significant
  enrichment" stand as a headline with a real signal underneath it.
- Never call agreement between arms "replication". The reference and the aligner are the
  experiment; they are not expected to agree. Say "holds under both aligners", or name it as
  a reference-by-aligner interaction.
- Difficulty and biology are not alternatives. Segmental duplications, the MHC, chrX and the
  acrocentric arms are hard to align *because they vary so much between people* that one
  linear reference represents them badly — and between-person variation is the substance of
  disease genetics, so the same property makes them disease-relevant. One cause, two
  consequences. Do not write "positional, not biological".
- **Describe, do not adjudicate.** This is a methods paper. Report where arms differ and what
  those regions are like. Do not write that an arm is better, more accurate, or correct. Lower
  background noise and more or larger results (eGenes, credible sets, betas, certainty) may be
  noted as mildly preferable — as an observation carrying a weak prior, never as a conclusion.
- **Error and biology are not in competition.** Never write "artifact rather than a finding",
  or "genotype error versus real biology". Reducing error is how biology becomes visible. A
  region whose measurement depends on representation *is* the finding, and the actionable
  form is: prior results here may be reference-dependent, and future work here should treat
  reference and aligner as analysis variables.
- **Unsigned means unsigned.** An allele-frequency or genotype discrepancy proves
  representation-dependence, not which arm erred. Truth may lie on either side.
- Prefer the weakest true form. If a comparison is n=7 against n=6, give the p-value and call
  it a description.

## Rebuilding the page

```
python3 collect_report_data.py   # in <workspace>/scratch — reads run roots, writes report_data.json
python3 build_figures.py         # writes figures/*.{light,dark}.svg
python3 build_site.py            # writes index.html
```

Before committing figures, discard any that differ only in matplotlib's embedded timestamp
and randomised element ids — otherwise every rebuild churns the whole `figures/` directory
and buries the real change.
