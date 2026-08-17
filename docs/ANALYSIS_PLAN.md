# BrainVar method-comparison eQTL analysis — central plan

**Read this first, and re-read it before starting any workstream.** Every piece of work in
this project exists to answer one question. If a proposed analysis does not obviously serve
it, that is a reason to challenge the analysis, not to widen the question.

## The question

> **Where** do reference genome, aligner and variant caller change cis-eQTL results, and what
> are the properties of those regions?

This is a **methods paper**. The job is to *describe differences*, not to crown a winner. No
workstream should be designed to show that an arm is better, and no result should be written
that way. Where an arm produces less background noise or a more biologically plausible answer,
that is worth noting; more or larger results (more eGenes, more credible sets, tighter
intervals, larger betas) are interesting and mildly preferable — **but no more than mildly.**
They are weak priors for interpretation, never conclusions, and never the point.

### The thesis

Most people believe alignment and reference choice do not matter much for genetic disease
research. **For most of the genome, they are right, and this project's first job is to confirm
it.** 91–95% of genes and 97% of association calls are untouched.

The argument is about the remainder. Regions that are hard to align are hard *because of
population-level variability* — they differ so much between people that a single linear
reference represents them badly. That same variability is what makes them disease-relevant,
because variation between people is the substance of disease genetics. Difficulty and
biological interest are not correlated by accident; they have one cause.

If that holds, then some of the fields with the longest history of difficulty —
**psychiatry, neurodevelopment, immunology** — are difficult partly because their loci are
exactly the loci where representation choices bite. The claim is not that these fields have
been getting wrong answers. It is that they are disproportionately exposed to a variable most
of the field treats as inert.

So: *for most genetics, this does not matter. For a small, interesting, historically
difficult subset, it does.*

### What the project delivers

Not a verdict on references, and not novel biology. A **map of where measurements are
representation-dependent**, with two uses:

- prior results in those regions should be read knowing the result may depend on reference or
  aligner;
- future work in those regions should treat reference and alignment method as analysis
  variables rather than as background.

## Prior result this builds on

Lalli, Bortvin, McCoy & Werling (2025), *A T2T-CHM13 recombination map and globally diverse
haplotype reference panel improves phasing and imputation*, bioRxiv 10.1101/2025.02.24.639687
— found 38% fewer assembly-discordant genotypes and 16% fewer switch errors on T2T-CHM13,
**with the largest gains on chromosome X and in regions flanking disease-causing CNVs**.

This project asks whether that genotype-level improvement propagates to a downstream
quantitative phenotype. The independent recovery of the same two hotspots (chrX; recurrent
CNV regions) from eQTL data is a central result, and the plan is built to test it rather
than to assume it.

**One method from that paper does not transfer.** It evaluated references against per-sample
assemblies as reference-neutral ground truth. BrainVar is short-read WGS with no native
assemblies, so **there is no ground truth in this dataset** and no analysis here can say which
arm is *correct*. This is the project's defining constraint and it must be stated wherever a
result could be read as an accuracy claim.

What substitutes, and — importantly — what each one can and cannot conclude:

| Substitute | What it establishes | What it does **not** establish | Used by |
|---|---|---|---|
| Allele-frequency concordance | The same donors and the same allele must give the same AF, so a difference proves the region is **representation-dependent** | *Which* arm is wrong. Truth could be the higher or the lower AF; the discrepancy is unsigned | Map representation-dependence |
| Hardy-Weinberg / excess heterozygosity | Het excess is a paralog-collapse signature, reference-neutral | That the arm without it is correct | Background-noise properties |
| Het allele balance from `AD` | A true diploid het is balanced regardless of reference | Anything about non-diploid loci, which is where this matters most | Background-noise properties; ancestry |
| Caller axis | Puts reference and aligner effects on the scale of a change nobody considers controversial | An ordering of the three axes by quality | Variant-caller axis |
| External replication | An outside cohort arbitrates | Cleanly — the outside cohort has its own reference and its own biases | not scoped |

**Error and biology are not in competition.** Identifying where measurement is
representation-dependent *is* how real biology becomes visible; the point of reducing error is
to see signal, not to trade one against the other. Do not write "genotype error versus real
biology", and do not describe a discordant region as "an artifact rather than a finding" — the
correct description is that its measurement depends on representation, and that is the finding.

## Design

Four arms, all 225 BrainVar donors, everything else held fixed. Two references (GRCh38,
T2T-CHM13) crossed with two aligners (linear, pangenome-graph surjected). A fifth and sixth
arm add the variant-caller axis (HaplotypeCaller, linear only, both references).

Contrasts are named `<positive>_minus_<negative>_<held fixed>`:

| Contrast key | Axis | Reads as |
|---|---|---|
| `t2t_minus_grch38_linear` | reference | T2T − GRCh38, aligner fixed at linear |
| `t2t_minus_grch38_graph` | reference | T2T − GRCh38, aligner fixed at graph |
| `graph_minus_linear_grch38` | aligner | graph − linear, reference fixed at GRCh38 |
| `graph_minus_linear_t2t` | aligner | graph − linear, reference fixed at T2T |
| `hc_minus_dv_grch38` | caller | HaplotypeCaller − DeepVariant, GRCh38 linear |
| `hc_minus_dv_t2t` | caller | HaplotypeCaller − DeepVariant, T2T linear |

**Factorial reading, not replication.** The arms are not repeated attempts at one experiment;
the reference, aligner and caller *are* the experiment. An effect present at both levels of
another factor is attributable to its own factor. An effect present at only one level is an
interaction, and the two factors are entangled at that locus. Never describe agreement
between arms as "replication".

## Workstreams

Each is independently claimable. Each ends in a run root with a `MANIFEST.json` and a
one-line status update in `STATE.md` — that is the resumable unit. Do not start a workstream
without reading its dependencies' manifests.

### Obtain concordance baseline `[COMPLETE]`
Establishes claim 1. Per-variant Pearson r, median |Δz|, share of variants moving < 0.5 z,
and agreement of association calls at |z| > 4, for every contrast.
Run root: `runs/gene_discordance_disease_20260816`.

### Calibrate: what counts as a large difference
**Add the variant-caller axis.** `[NOT STARTED]` The yardstick. Re-run the genotype→association
path with HaplotypeCaller in place of DeepVariant, linear arms, both references, then score
`hc_minus_dv_*` through the identical window and gene pipelines. If reference-swap discordance
is smaller than caller-swap discordance, that is the strongest possible reassurance statement;
if larger, it is the strongest caution. Inputs verified present and identically processed:
`inputs/genotypes/callset_factorial/linear_{grch38,t2t}_haplotypecaller.vcf.gz` →
`gatk_*_haplotypecaller.joint_called.biallelic.unique_variant_ids.vcf.gz` (64–67 GB).
*Cost: high (full genotype derivation + association). Restrict to flagged windows first.*

**Measure the sampling-noise floor.** `[DEFERRED — low priority]` Same arm, disjoint donor halves, to give
the sampling-noise floor. Deprioritised: unlikely to yield trustworthy FDR estimates and not
central to the paper's argument.

### Anatomy of the difference
Characterises *how* estimates move, not which movement is correct.

**Decompose Δz into Δβ and Δse.** `[NOT STARTED]` A shifted β means the arms estimate a
different underlying effect; a shifted se means the same effect measured with different
precision. Hypothesis to test: aligner discordance is se-dominated (precision) and reference
discordance is β-dominated (different variant/LD). Nearly free — `beta_a/se_a/beta_b/se_b`
are already in the paired tables. *Cost: low.*

**Describe background-noise properties per arm.** `[NOT STARTED]` Describe, per arm and per region
class: median standard error at matched MAF, residual variance, HWE/excess-het rate, and mean
het allele balance. Report where the arms differ and by how much. A consistent reduction in
background noise is worth stating and is a mild point in an arm's favour — write it as an
observation, not a verdict, and always alongside the region class it occurs in. The question of
interest is *where* the arms' noise properties diverge, and whether that is the same place the
effect estimates diverge. *Cost: medium.*

**Map representation-dependence by allele-frequency concordance.** `[NOT STARTED]`
For an identical normalised allele in the same 225 donors, AF must be identical. A difference
therefore proves the measurement at that site **depends on representation** — with no external
truth set required, because the agreement is logically compelled rather than empirically hoped
for. It is **unsigned**: it does not say which arm is wrong, and truth may be the higher or the
lower frequency. The deliverable is a per-region catalogue of representation-dependent
measurement, which is the caution map described at the top of this document, and a comparison
of its extent across all three axes. *Cost: low. High value — it is the cleanest instrument
the project has.*

**Verify allele identity, then count sign flips.** `[NOT STARTED — gated]` Only interpretable after confirming REF/ALT are
identical. The paired tables already match on exact LiftoverIndel-normalised T2T
position/REF/ALT and orient slopes to the T2T ALT allele, so the gate is expected to pass, but
**verify and report the verification** before reporting any flip rate. A flip after that
verification is a genuinely different claim from a magnitude change. *Cost: low.*

### Asymmetric testability: what was never askable
Every matched-variant comparison is blind to this by construction.

**Characterise genes testable in one reference only.** `[NOT STARTED]` Universes are 18,273 (GRCh38)
vs 18,114 (T2T), frozen natively and never intersected. Characterise the non-overlap: which
genes, what biotype, what sequence context, are the T2T-only ones in newly resolved regions.
*Cost: low. High value — nothing done so far can see this.*

**Measure signal in newly accessible sequence.** `[NOT STARTED]` Define "fills in" concretely as: for each
gene, the mean association signal (mean |z|, and count of variants at |z| > 4) among cis
variants that fall in sequence present in T2T but absent or unmappable in GRCh38, versus the
rest of that gene's cis window. The question is whether newly sequenced sequence carries
*signal* or only *variants*. Requires a defensible "newly accessible" interval set — derive
from GRCh38-unmappable / non-syntenic regions rather than asserting one. *Cost: medium.*

**Test whether arm-exclusive variants carry signal.** `[NOT STARTED]` ~388k T2T-only and ~284k GRCh38-only
identities within the matched frame. Do they ever become lead variants or reach significance?
**Stratify by region class** — a GRCh38-exclusive variant in an ENCODE-blacklist or otherwise
known-difficult region means something very different from one in ordinary sequence. Without
that stratification the comparison is not interpretable. *Cost: medium.*

### Downstream consequence: does the conclusion change
**Measure lead-variant switching by LD.** `[NOT STARTED]` For genes that are eGenes in both arms,
report **LD (r²) between arm A's lead and arm B's lead** as the primary metric, plus physical
distance between the two leads and each lead's signed distance from the TSS. A lead that moves
within an LD block is a cosmetic change; one that moves to an independent block changes the
causal hypothesis. *Cost: medium — needs an LD reference; use the arm's own genotypes.*

**Count credible sets per gene.** `[NOT STARTED]` Number of SuSiE credible sets per gene, and
credible-set membership overlap, per arm. Preferred over conditional/stepwise analysis on cost
grounds. **First task: locate and validate existing SuSiE output.** Candidates seen outside
this workspace: `nf_stage/brainvar_eqtl_e36_susieash_complete_e42388e_20260728T054734Z` and
siblings — but these are e36-era, *not* the current v4/k35 arms, so confirm arm correspondence
before use and re-run if they do not match. *Cost: low if reusable, high if not.*

**Test whether GWAS colocalization changes.** `[NOT STARTED]` The ideal endpoint: does a coloc call appear or
disappear between arms? GWAS rsIDs in T2T coordinates are available at
`genome_refs/T2T-CHM13_v2_ncbi110/chm13v2.0_GWASv1.0rsids_*.vcf.gz` (note: GenBank contig
names, needs renaming; carries dbSNP/ClinVar INFO but **no trait field** — a trait source must
be added). *Cost: high. Schedule last; it depends on lead-variant switching and credible sets.*

**Measure top-k rank concordance.** `[NOT STARTED]` Genome-wide r describes the average result,
but nobody acts on the average — people act on the head of the list: which genes get followed
up, fine-mapped, or put in a figure. Report the overlap of the top 100 / 500 genes between
arms. This converts a statistical statement into an operational one: *how many of your
follow-up candidates would change if you switched method.* That is the number a reader
actually needs. *Cost: low.*

### Who is affected
**Separate the chrX effect by sex.** `[NOT STARTED — high priority]` graph-vs-linear on T2T moves 22.3% of
chrX genes (OR 4.47, p = 4×10⁻³⁰) and no other contrast does. This is unexpected and should be
treated as suspect until explained. XY donors are hemizygous, T2T uses a masked-Y assembly, and
X alignment in males is a different problem than in females. Stratified XX and XY maps already
exist. Determine whether the effect is XY-driven, XX-driven, or both, and whether it is an
artifact of ploidy handling. *Cost: low. Do this early — it may invalidate a headline.*

**Test ancestry-dependence of discordance.** `[NOT STARTED]` Haplotype sampling matches donors to a reference
panel, so donors whose ancestry the panel covers poorly receive worse subgraphs. Regress
per-donor missingness and mean het allele balance on genotype PCs, separately per arm. If the
graph arms show PC dependence the linear arms do not, method choice interacts with ancestry —
a generalisability finding, not a technical one. Ancestry assignments exist at
`runs/wp01_ancestry_hwe_20260808_v3/brainvar_ancestry_assignments.tsv`. *Cost: low.*

### Characterise discordant regions and gene classes `[COMPLETE]`
Method-sensitive genes are enriched for ClinGen recurrent-CNV regions (OR 2.9–9.1), human
accelerated regions (OR 1.7–1.9) and segmental duplications (OR 1.7–2.5), in all four
contrasts. MHC OR 2.5–3.8 in all four. Run roots: `runs/discordant_gene_classes_20260816`,
`runs/gene_discordance_disease_20260816`, `runs/window_discordance_exceedance_20260816`,
`runs/directional_window_discordance_20260816`.
*Open sub-task:* re-test these enrichments controlling for expression level and cis-variant
count, since duplicated-region genes are often lowly expressed and the enrichment could be
partly a power artifact (see *Describe background-noise properties per arm*).

### Build the report page `[IN PROGRESS]`
The GitHub page. See `CONVENTIONS.md` for figure and prose conventions, which follow the 2025
paper. Current page: `index.html`, built by `build_site.py` from `report_data.json`.

## Sequencing

```
Concordance baseline ──> Decompose Δz ─┐
                    └──> One-reference-only genes ─┤
                                                   ├──> Top-k concordance ──> Page
Separate chrX by sex ──(may invalidate)────────────┘

Variant-caller axis ──> recalibrates every effect size reported
Lead-variant switching ──> Credible sets ──> GWAS colocalization
```

Recommended order: **separate the chrX effect by sex** (it may invalidate a headline), then
**decompose Δz** and **map representation-dependence** together (cheap, mechanistic, and they
establish what kind of difference each axis produces), then **characterise one-reference-only
genes**, then **add the variant-caller axis** (expensive, but it recalibrates every effect
size reported), then **lead-variant switching** and **credible sets**, then **top-k concordance**,
with **GWAS colocalization** last.

## Rules

- No number reaches the page except through `report_data.json`. Prose cannot be allowed to
  drift from the run tree.
- Every run root carries a `MANIFEST.json` with inputs, SHA-256 of the tracks and tables it
  read, parameters, and an explicit `inferential_boundary` string naming what the statistic
  does *not* establish.
- State the weakest form of every claim. Where a test is underpowered, say so with the number.
- Report negatives immediately followed by whatever *did* move. A null is a result; a null
  presented as the headline while a real signal sits underneath is a reporting failure.
