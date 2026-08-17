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
**Add the variant-caller axis.** `[COMPLETE — allele-frequency stage]` Run root
`runs/caller_axis_af_20260816`; script `scratch/analyze_caller_axis_af_20260816.py`.

*Scoped in two stages.* The full genotype-to-association rerun is expensive and gated. The
allele-frequency instrument gives the yardstick without it: frequencies are recomputed
identically from each callset over the same donors, and because every contrast here holds the
reference fixed, variant identifiers (`chrom_pos_ref_alt`) match directly with no liftover.
That places the caller and the aligner on one footing immediately. The association-level
version remains outstanding and is what would let the caller axis be scored through the same
window and gene pipelines as the other axes.

*Correctness detail worth keeping:* two analysis donors are absent from the HaplotypeCaller
callsets (`5212_D2` from GRCh38, `6085_D2` from T2T). Without intersecting the donor sets
`bcftools` silently skips them in one callset and includes them in another, so frequencies
would be computed over different samples. The run intersects to the 223 donors present in all
six callsets.

*Result.* Run roots `runs/caller_axis_af_20260816` and
`runs/three_axis_af_yardstick_20260816`. Six callsets, 223 donors present in all of them, six
contigs (chr1, chr6, chr17, chr19, chr22, chrX), frequencies computed identically throughout so
no axis has an advantage of stage, donor set or region.

| Axis | Sites differing by more than 0.01 | Relative to an aligner swap |
|---|---|---|
| Reference (T2T − GRCh38) | 2.34%, 2.61% | **0.91×** |
| Aligner (graph − linear) | 2.56%, 2.91% | 1.00× |
| **Caller (HaplotypeCaller − DeepVariant)** | **5.84%, 6.50%** | **2.26×** |

**Swapping the reference genome perturbs allele frequencies slightly less than swapping the
aligner, and both are under half what swapping the variant caller does.** This is the
reassurance statement the project was missing: the reference choice is a smaller perturbation
than a routine tooling decision most groups make without discussion.

*The caveat is essential and changes what the number means.* Cross-reference contrasts can only
be computed on variants carrying a unique normalised identity in both references — about 1.9
million of the roughly 10 million each callset holds. Same-reference contrasts use all of them.
So the reference axis is measured **only among variants both references can represent**, and
the variants only one reference can represent are excluded entirely.

That is not a flaw; it is the finding, stated precisely. Where two references can both see a
variant, they agree about it better than two aligners do. The reference's distinctive effect is
not in *measuring shared variants differently* — it is in **what it can see at all**, which is
what the gene-universe and arm-exclusive-variant results measure. The two halves fit together:
a reference swap is a small perturbation to the shared part of the genome and a large change to
which part is shared.

Original statement of the task: re-run the genotype→association
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

**Decompose Δz into Δβ and Δse.** `[COMPLETE]` Run root
`runs/effect_precision_and_af_20260816`. The split is exact, not approximate:
`z_a - z_b = (b_a - b_b)/se_a + b_b(1/se_a - 1/se_b)`, verified per contrast to a maximum
residual of 3x10^-6.

*The stated hypothesis was wrong.* The plan predicted aligner discordance would be
precision-dominated and reference discordance effect-dominated. **Both axes are
effect-dominated, and to nearly the same degree** — 93.9% of variants for the reference swap,
94.9% for the aligner swap. Precision barely moves: the median se ratio between arms is 1.0017
(reference) and 1.0000 (aligner).

| Contrast | var(Δz) | var(effect term) | var(precision term) | effect-dominated |
|---|---|---|---|---|
| T2T − GRCh38 · linear | 0.2328 | 0.2499 | 0.0132 | 93.9% |
| T2T − GRCh38 · graph | 0.2568 | 0.2700 | 0.0143 | 93.9% |
| graph − linear · GRCh38 | 0.0102 | 0.0107 | 0.0012 | 94.9% |
| graph − linear · T2T | 0.0111 | 0.0112 | 0.0013 | 94.9% |

**Neither method change buys precision. Both change what is being measured.** That bears
directly on the signal-to-noise framing: where these methods differ they are not reducing noise
around a fixed quantity, they are estimating a different one — consistent with the arms
effectively seeing different genotypes or different local linkage at the sites that move.

*One nuance that does track sequence class:* in duplicated sequence the aligner's precision term
roughly doubles relative to its effect term (0.111 to 0.218), while the reference's barely
shifts (0.117 to 0.138). Pangenomic alignment measurably changes genotype certainty in
duplicated sequence specifically, even though the estimate still dominates.

**Describe background-noise properties per arm.** `[NOT STARTED]` Describe, per arm and per region
class: median standard error at matched MAF, residual variance, HWE/excess-het rate, and mean
het allele balance. Report where the arms differ and by how much. A consistent reduction in
background noise is worth stating and is a mild point in an arm's favour — write it as an
observation, not a verdict, and always alongside the region class it occurs in. The question of
interest is *where* the arms' noise properties diverge, and whether that is the same place the
effect estimates diverge. *Cost: medium.*

**Map representation-dependence by allele-frequency concordance.** `[COMPLETE]` Run root
`runs/effect_precision_and_af_20260816`. Frequencies are placed in the common frame first, with
the reported frequency complemented wherever the identity mapping flipped the allele — which is
also the gate that had to pass before any sign-flip statistic could mean anything.

The instrument is clean: **median |ΔAF| is exactly zero in all four contrasts.**

| Contrast | Exactly equal | Within 0.001 | Within 0.01 | Beyond 0.01 |
|---|---|---|---|---|
| T2T − GRCh38 · linear | 70.2% | 85.7% | 96.2% | **3.8%** |
| T2T − GRCh38 · graph | 70.3% | 85.4% | 96.1% | **3.9%** |
| graph − linear · GRCh38 | 83.1% | 87.8% | 98.4% | **1.6%** |
| graph − linear · T2T | 81.6% | 87.0% | 98.3% | **1.7%** |

The 1.6–3.9% beyond 0.01 is the caution map, and it is unsigned: those sites' measurement
depends on representation, but truth may lie on either side. Maximum |ΔAF| reaches 0.84–0.90, so
a minority of sites disagree almost completely.

*Convergent evidence, from an instrument sharing no statistics with the eQTL analysis:* the
chromosomes with the most AF discordance are chr19, chr22, chr17, chrX, chr6, chr21, chr15 and
chr9 — the same set the window and gene analyses independently identified. chr21 is worst in
graph − linear · T2T (3.9%), matching the acrocentric result from that contrast. Genotype
frequencies and association statistics point at the same regions without sharing a denominator.

**Verify allele identity, then count sign flips.** `[NOT STARTED — gated]` Only interpretable after confirming REF/ALT are
identical. The paired tables already match on exact LiftoverIndel-normalised T2T
position/REF/ALT and orient slopes to the T2T ALT allele, so the gate is expected to pass, but
**verify and report the verification** before reporting any flip rate. A flip after that
verification is a genuinely different claim from a magnitude change. *Cost: low.*

### Asymmetric testability: what was never askable
Every matched-variant comparison is blind to this by construction.

**Characterise genes testable in one reference only.** `[COMPLETE]` Run root
`runs/gene_universe_asymmetry_20260816`.

The headline is that the *net* difference hides the real one. GRCh38 tests 18,273 genes and
T2T 18,114 — a net gap of 159 — but the universes are not nested: **667 genes are testable
only on GRCh38 and 508 only on T2T**, a turnover of 1,175 genes that every matched-gene
comparison in this project is structurally blind to.

**278 of those exclusive genes are eGenes** (145 GRCh38-only, 133 T2T-only): real associations
callable on one reference and unavailable on the other. Exclusive genes are *more* likely to be
eGenes than shared ones — 26.2% for the T2T-only set and 21.7% for GRCh38-only, against 20.8%
among shared genes.

| | GRCh38-only | T2T-only | Shared |
|---|---|---|---|
| Genes | 667 | 508 | 17,606 |
| eGenes | 145 (21.7%) | 133 (26.2%) | 20.8% |
| Largest biotype | lncRNA 549, protein-coding 86 | lncRNA 287, protein-coding 159, rRNA 27 | — |
| Absent from other annotation | 561 of 667 | 266 of 508 | — |
| Mean segdup fraction | not measurable in T2T frame | **0.242** | 0.063 |
| Any segdup overlap | — | **42.1%** | 27.8% |

T2T-exclusive genes sit in sequence roughly **four times as duplicated** as shared genes, which
is the same signature every other workstream has found. The 27 rRNA genes testable only on T2T
are the expected consequence of T2T resolving the rDNA arrays GRCh38 leaves as gaps.

*What the numbers cannot settle:* a majority of exclusive genes are absent from the other
reference's annotation entirely (561 of 667; 266 of 508), so annotation release differences —
RefSeq on GRCh38 against RefSeq Liftoff on T2T — account for much of the turnover rather than
sequence accessibility. The residual (106 and 242 genes present in the other annotation but not
testable) is the part attributable to the reference itself. Both figures should be quoted
together; neither alone is honest.

**Measure signal in newly accessible sequence.** `[COMPLETE]` Run root
`runs/newly_accessible_signal_20260816`. Genome-wide, all 23 contigs.

"Newly accessible" is defined operationally rather than asserted: a variant is accessible to
one arm only when it carries a unique normalised identity there and has no counterpart among
the other arm's tested variants. The sharp question is then per gene — is its single best
association at a variant the other arm cannot represent at all?

| Contrast | Genes whose lead is exclusive | Share | Among genes reaching \|z\| ≥ 4 |
|---|---|---|---|
| T2T − GRCh38 · linear | 2,036 of 18,114 | **11.2%** | 8.9% |
| T2T − GRCh38 · graph | 2,079 of 18,114 | **11.5%** | 9.4% |
| graph − linear · T2T | 864 of 18,114 | 4.8% | 5.0% |

**About one gene in nine has its top association at a variant GRCh38 cannot represent.**

This closes the loop the yardstick opened. Taken together the three results say the same thing
from different directions:

- among variants both references can represent, they agree more closely than two aligners do
  (2.34% of frequencies shifted, against 2.56%);
- but 1,175 genes are testable on one reference only, and 278 of those are eGenes;
- and for 11% of *shared* genes, the strongest association sits on a variant the other
  reference cannot represent.

**The reference's distinctive effect is on access, not on measurement.** Where it can see the
same thing it sees it the same way; what changes is what it can see. That is the mechanical
core of the argument, and no yield comparison can reach it.

*Bounds.* Variants with no unique normalised identity are excluded from both sides, so the
exclusive set is a lower bound on what one arm can represent and the other cannot. A gene's
lead being exclusive does not mean the other arm finds nothing for that gene — only that it
cannot find that variant.

**Test whether arm-exclusive variants carry signal.** `[COMPLETE — chr1 and chr19]`
Run root `runs/arm_exclusive_variants_20260816`.

A variant is exclusive to an arm when its normalised common-frame identity has no counterpart
among the other arm's tested variants. Signal is the best |z| it achieves against any gene,
compared against that same arm's shared variants.

| Contrast | Exclusive variants | Reach \|z\| ≥ 4 | Shared variants do | Ratio | Duplication content |
|---|---|---|---|---|---|
| T2T − GRCh38 · linear | 44,996 | 5.5% | 6.0% | **0.92×** | 0.084 vs 0.065 |
| graph − linear · T2T | 22,785 | **12.4%** | 5.9% | **2.09×** | 0.154 vs 0.065 |

**The two axes add different kinds of variant.** Variants only the pangenome graph can call
are *twice as likely* to carry an association as the variants both aligners share, and they sit
in sequence two and a half times as duplicated. Variants only T2T can call are, if anything,
slightly *less* likely to carry one than shared variants.

So the aligner is not merely adding variants — it is adding variants that produce signal, in
exactly the duplicated sequence where a linear aligner has least to work with. The reference
adds variants of ordinary informativeness. That distinction is invisible to any yield count,
which sees both as "more variants tested".

*Bounds on this.* Computed on chr1 and chr19 only, and restricted to variants carrying a
unique normalised identity — so variants that cannot be placed in the common frame at all are
excluded from both groups, and those are the most reference-specific of all. Summarising a
variant by its best |z| against any gene favours variants tested against more genes.

### Downstream consequence: does the conclusion change
**Measure lead-variant switching by LD.** `[COMPLETE]`
Run root `runs/lead_variant_switching_20260816`.

Restricted to genes that are eGenes in **both** arms, so yield differences cannot contribute.
Leads are placed in the common LiftoverIndel-normalised frame before comparison.

| Contrast | Same lead variant | Lead changes | Median move | Crosses the TSS |
|---|---|---|---|---|
| T2T − GRCh38 · linear | **58.4%** | 41.6% | 14.3 kb | 26% |
| T2T − GRCh38 · graph | **59.3%** | 40.7% | 14.4 kb | 23% |
| graph − linear · GRCh38 | 78.5% | 21.5% | 11.2 kb | 20% |
| graph − linear · T2T | 76.2% | 23.8% | 13.0 kb | 19% |

This is the most consequential result for downstream work so far. Two arms can agree that a
gene has an eQTL, agree closely on its effect size, and still **nominate a different variant
four times out of ten** under a reference swap. For fine-mapping and colocalisation the
identity of the lead is not an incidental detail — it is the result.

Only about an eighth of the moves are under a kilobase, so these are not sub-resolution
jitter; the median move is over ten kilobases and a fifth to a quarter cross to the other side
of the transcription start site.

*The LD stage resolves it, and it softens the result.* Run root
`runs/lead_switch_ld_20260816`. r² between the two leads is computed in the positive arm's own
genotype matrix over pairwise-complete donors, so the leads are compared as that arm measured
them rather than through an external panel carrying its own reference.

| Contrast | Median r² | Same signal (r² > 0.8) | Independent (r² < 0.2) |
|---|---|---|---|
| T2T − GRCh38 · linear | 0.852 | **56.0%** | 11.9% |
| T2T − GRCh38 · graph | 0.854 | 56.1% | 13.2% |
| graph − linear · GRCh38 | 0.929 | 70.8% | 4.0% |
| graph − linear · T2T | 0.934 | 69.0% | 6.2% |

**Most lead switches are relabelling within a haplotype block, not a new causal hypothesis.**
Compounding the two stages: a reference swap changes the lead for about 41% of shared eGenes,
and roughly an eighth of those changes are to an effectively independent variant — so on the
order of **5% of shared eGenes** get a genuinely different causal candidate. For an aligner
swap it is closer to 1%.

That is a far weaker claim than the distance stage on its own implied, and the page has been
corrected to state it this way. The distance result was not wrong, but read alone it invited
the conclusion that four in ten genes get a different answer, which the LD says they do not.

*One conservative bias worth keeping:* pairs are restricted to switched leads whose counterpart
is nameable in the positive arm's own variant space. A lead with no counterpart there is
excluded — and those are precisely the cases where the two arms differ most, so the independent
share is understated rather than inflated.

**Count credible sets per gene.** `[BLOCKED — existing SuSiE output is not reusable]`

Checked and rejected. The candidate outputs
(`nf_stage/brainvar_eqtl_e36_susieash_complete_e42388e_20260728T054734Z` and siblings) do
carry credible sets for GRCh38 and T2T, about 49,000 credible-set variants each, with
cross-reference matching already performed. They cannot be used here for two independent
reasons:

- they are **e36**, a different expression-PC parameterisation from the k = 35 the current
  arms use; and
- they were validated on **2026-07-28**, before the genotype rebuild the four current arms
  rest on, so their credible sets were fit on genotypes from which a large share of variants
  was missing.

Either alone would disqualify them. Producing credible sets for the current arms means a fresh
SuSiE run, which is a substantial compute job and has not been scheduled. Until then, no
statement about credible sets belongs on the page.

**Test whether GWAS colocalization changes.** `[NOT STARTED]` The ideal endpoint: does a coloc call appear or
disappear between arms? GWAS rsIDs in T2T coordinates are available at
`genome_refs/T2T-CHM13_v2_ncbi110/chm13v2.0_GWASv1.0rsids_*.vcf.gz` (note: GenBank contig
names, needs renaming; carries dbSNP/ClinVar INFO but **no trait field** — a trait source must
be added). *Cost: high. Schedule last; it depends on lead-variant switching and credible sets.*

**Measure top-k rank concordance.** `[COMPLETE]` Run root
`runs/topk_rank_concordance_20260816`.

The operational translation of the concordance result, and the two do not say the same thing.
Genome-wide the arms correlate at 0.94 under a reference swap and 91% of genes are unmoved —
but at the head of the ranking, where decisions are actually made:

| Contrast | Top-100 overlap | Candidates changed | Rank correlation |
|---|---|---|---|
| T2T − GRCh38 · linear | 78% | **22 of 100** | 0.837 |
| T2T − GRCh38 · graph | 81% | **19 of 100** | 0.830 |
| graph − linear · GRCh38 | 95% | 5 of 100 | 0.963 |
| graph − linear · T2T | 95% | 5 of 100 | 0.964 |

**Switch reference and about a fifth of your top hundred follow-up candidates change.** Switch
aligner and one in twenty does. Five of the reference-swap changes are genes not testable in
the other arm at all, which links this directly to the gene-universe result.

This is the number a reader planning follow-up actually needs, and it is markedly less
reassuring than the genome-wide correlation. Both belong on the page: the bulk of the map is
stable, and the part of it anyone acts on is less so.

### Who is affected
**Separate the chrX effect by sex.** `[COMPLETE]` Run root
`runs/chrx_discordance_by_sex_20260816`; script
`scratch/analyze_chrx_discordance_by_sex_20260816.py`, dosage audit
`scratch/audit_chrx_dosage_coding_20260816.py`.

*Dosage audit complete and it clears the first hurdle:* chrX hemizygous genotypes are encoded
identically in all four arms — XY het-like dosage share 1.7–2.5%, XX 28.9%, arm-to-arm
variation negligible. The two arms of every contrast therefore encode chrX the same way, so a
differential-ploidy-encoding artifact is excluded.

*Result, full autosomal baseline, all 23 contigs:* the effect is **XX-driven, not XY-driven** —
the opposite of the ploidy hypothesis, and stronger than the whole-cohort figure suggested.

| Contrast | XX chrX vs autosomal | XY chrX vs autosomal |
|---|---|---|
| T2T − GRCh38 · linear | 3.6% vs 5.0% (OR 0.71) | 5.3% vs 5.8% (OR 0.91) |
| T2T − GRCh38 · graph | 3.6% vs 5.0% (OR 0.72) | 5.0% vs 5.8% (OR 0.84) |
| graph − linear · GRCh38 | 1.3% vs 6.8% (OR 0.18) | 1.1% vs 6.7% (OR 0.16) |
| **graph − linear · T2T** | **34.6% vs 6.8% (OR 7.27, p = 3.8×10⁻⁸⁵)** | 1.8% vs 7.0% (OR 0.24) |

Seven of eight cells show nothing on chrX. The eighth moves **more than a third of chrX genes**.
It is not a power artifact — XY has 133 donors against XX's 92, so the null cell has *more*
power than the positive one. In XY, chrX is if anything *depleted* relative to autosomes.

*Mechanism, consistent with all four rows:* pangenomic alignment resolves haplotype diversity.
A hemizygous X carries one haplotype and offers nothing to resolve, which is why XY shows
nothing. GRCh38's X is itself poorly resolved, which is why the graph cannot exploit two
haplotypes there either. Only T2T's complete X plus two X haplotypes gives the graph both the
material and the opportunity — a three-way reference × aligner × sex interaction.

*Consequence for reporting:* the whole-cohort figure of 22.3% is a **dilution** of a stronger
XX-specific effect. Any chrX analysis in a mixed-sex cohort should stratify by sex, and the page
must state the result as XX-specific rather than as a property of chrX.

Original statement of the task: graph-vs-linear on T2T moves 22.3% of
chrX genes (OR 4.47, p = 4×10⁻³⁰) and no other contrast does. This is unexpected and should be
treated as suspect until explained. XY donors are hemizygous, T2T uses a masked-Y assembly, and
X alignment in males is a different problem than in females. Stratified XX and XY maps already
exist. Determine whether the effect is XY-driven, XX-driven, or both, and whether it is an
artifact of ploidy handling. *Cost: low. Do this early — it may invalidate a headline.*

**Test ancestry-dependence of discordance.** `[COMPLETE]` Run root
`runs/ancestry_dependence_20260816`. Cohort composition: EUR 104, AFR 62, AMR 43, plus 16
donors in smaller or uncertain groups.

Two per-donor quantities from the genotype matrices the association actually used —
missingness, and alternate-allele load as a proxy for how far a donor sits from the reference.
The comparison is the *change* within a donor between arms, which removes any constant
per-donor effect.

**The reference swap is not ancestry-neutral, and the direction reverses between groups.**
Moving from GRCh38 to T2T changes median alternate load by:

| Ancestry | Change in alternate load |
|---|---|
| EUR (n = 104) | **−0.047** |
| AMR (n = 43) | −0.029 |
| AFR (n = 62) | **+0.017** |

Kruskal-Wallis p = 1.3×10⁻³⁸, and the pattern is identical under both aligners. European-
ancestry donors come to look *more* like the reference on T2T and African-ancestry donors
*less* like it. That is the expected consequence of the two references' provenance —
T2T-CHM13 derives from a single haploid European-ancestry cell line, while GRCh38 is a mosaic
with substantial African-ancestry contribution — but it is worth measuring rather than
assuming, and it means the reference choice interacts with cohort composition.

The aligner axis shows the same test at roughly a fiftieth of the magnitude
(+0.0011 to +0.0007 across groups within GRCh38), and missingness changes on every axis are
small in absolute terms (+0.0002 to +0.0007) even where significant.

*What this does not establish.* Alternate load describes reference bias; it does not show that
African-ancestry donors get worse eQTL estimates on T2T. Linking the shift to per-donor
contribution to discordance is the natural follow-up and has not been done. Ancestry labels
are projected superpopulation assignments rather than self-reported, and the groups outside
EUR, AFR and AMR are too small to test.

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
