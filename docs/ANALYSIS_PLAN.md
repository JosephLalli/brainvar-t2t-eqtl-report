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

**Background noise, standard-error stage.** `[COMPLETE]`
Run root `runs/background_noise_by_arm_20260816`.

Standard error at matched allele frequency is the natural measure: an arm whose standard
errors are smaller at the same frequency is extracting more from the same donors. Reported per
region class, because the question worth asking is not whether an arm is quieter overall but
whether it is quieter where the sequence is hard.

**No arm is measurably less noisy, anywhere.** Median standard error at MAF 0.2–0.35 is 0.04579,
0.04581, 0.04578 and 0.04581 across the four arms — identical to the fourth decimal. Ratios
between arms:

| Comparison | In duplicated sequence | In ordinary sequence |
|---|---|---|
| graph ÷ linear, GRCh38 | 0.999 | 1.000 |
| graph ÷ linear, T2T | 1.000 | 1.001 |
| T2T ÷ GRCh38, linear | **1.008** | 1.003 |

The largest effect in the table is T2T carrying standard errors 0.8% *larger* than GRCh38 in
duplicated sequence. Standard errors are higher in duplicated sequence than in ordinary
sequence in every arm (about 0.0536 against 0.0522), and no method change repairs that.

This confirms the decomposition and extends it: that result showed the *difference* between
paired estimates is carried by the effect and not the standard error; this shows the *level* of
noise is the same in every arm, including in exactly the hard sequence where a gain would be
most expected.


---

**Describe background-noise properties per arm.** `[COMPLETE — all three stages]`
Run roots `runs/background_noise_by_arm_20260816` (standard error) and
`runs/background_noise_hwe_ab_20260817` (Hardy-Weinberg, allele balance).

The standard-error stage is recorded above and found no arm quieter than another anywhere. The
two remaining instruments ask a different question, and give a different answer.

*Why these two.* BrainVar has no per-donor assemblies, so nothing can say which arm is correct.
These come closest, because they test properties a callset must have whatever reference
produced it. Two paralogous copies collapsed onto one locus make every carrier look
heterozygous, so **excess heterozygosity** beyond Hardy-Weinberg expectation is the classic
collapse signature. And a true diploid heterozygote should draw about half its reads from each
allele, so a **skewed heterozygote allele balance** means the reads at that site are not coming
from one diploid locus. Neither says the arm without the signature is right; both say the arm
with it is more often measuring something other than a single diploid locus.

Scanned across chr1, chr6, chr16, chr22 and chrX in all four DeepVariant arms, requiring at
least 10 heterozygotes per site and 10 reads per heterozygote, then **matched on
1,654,230 variants present in all four arms** in the common frame.

| Arm | Excess het, ordinary | Excess het, duplicated | Skewed-het share, ordinary | Skewed-het share, duplicated |
|---|---|---|---|---|
| linear · GRCh38 | 0.204% | **0.283%** | 0.0396 | **0.0648** |
| graph · GRCh38 | 0.195% | **0.087%** | 0.0364 | **0.0388** |
| linear · T2T | 0.200% | **0.262%** | 0.0393 | **0.0561** |
| graph · T2T | 0.192% | **0.082%** | 0.0362 | **0.0388** |

**Pangenomic alignment removes the duplicated-sequence penalty; the reference change does
not.** In the linear arms duplicated sequence carries clearly elevated collapse signatures —
excess heterozygosity 0.283% against
0.204% in ordinary sequence, and a skewed-het
share of 0.0648 against
0.0396. In the graph arms that elevation
is gone: the duplicated skewed-het share
(0.0388) is essentially its ordinary-sequence
value (0.0364).

Paired on the same variants, the aligner swap in duplicated sequence moves excess heterozygosity
by **−0.195 percentage points** on GRCh38 and −0.180 on T2T, and the mean skewed-het share by
−0.026 and −0.017 (p < 1e-150 throughout). The reference swap in the same sequence moves excess
heterozygosity by −0.021 and −0.005 points, and its effect on allele balance is **not
significant at all** on the graph axis (p = 0.75).

**This is the first result in the project where one arm is cleaner, and it does not contradict
the standard-error stage — it separates two things that were being conflated.** A standard
error is the precision of an estimate given the genotypes. Excess heterozygosity and allele
balance ask whether the genotypes describe one diploid locus. Pangenomic alignment buys no
precision and does buy genotype fidelity in duplicated sequence, which is exactly where a
haplotype-aware method should help.

*Bounds.* The direction is partly by construction: haplotype-aware alignment is designed to
stop paralogous reads collapsing, so finding that it does is confirmation the method works as
intended rather than a surprise. What is new is the magnitude, measured on instruments that
need no ground truth. The matched set also requires a variant to be called and to pass the
heterozygote-count filter in *all four* arms, which in duplicated sequence selects the more
tractable sites — 1,654,230 matched against 46,913–53,850 per arm marginally —
so this understates the difference at the hardest sites rather than overstating it. One further
oddity is worth flagging rather than celebrating: in the graph arms excess heterozygosity in
duplicated sequence is *lower* than in ordinary sequence, which that selection may explain and
which is not otherwise expected.

**Consequence for the project's framing.** Signal-to-noise, in the ordinary sense of a smaller
standard error around the same quantity, is not what any of these method changes buy. What they
change is which variants exist to be measured at all. Read alongside the access result, the two
say the same thing from opposite ends.

*Bound.* A smaller standard error is not by itself better; it is an improvement only if the
estimate it surrounds has not also moved, which this cannot establish. Nothing here ranks the
arms. Hardy-Weinberg and allele-balance stages of this workstream remain unrun and would need
passes over the callsets.

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

**Count credible sets per gene.** `[COMPLETE]`
Run root `runs/susie_finemapping_v4_k35_20260817`; script
`scratch/run_susie_finemapping_v4_k35_20260817.py`.

*Correction to what this plan previously said.* It recorded that the existing e36 fine-mapping
was disqualified partly because its credible sets were fit on genotypes missing a large share
of variants. **That was wrong.** Its adapter did its own NaN-aware mean imputation
(`impute_and_filter_genotypes`) rather than relying on TensorQTL's `-9` sentinel, so the
genotypes it fine-mapped were handled correctly and the missing-dosage defect never reached it.

The run is still not reusable, for two reasons that do hold:

- it was fit at **36** expression principal components, not the 35 these arms use; and
- its **eGene selection** came from permutation results produced before the genotype rebuild,
  so the genes chosen for fine-mapping were chosen from a contaminated map even though the
  fine-mapping itself was sound.

*What the new run enforces.* `tensorqtl.susie.map` calls the same `impute_mean` that
recognises a `-9` sentinel and not IEEE NaN, so SuSiE by that path **is** exposed to the defect
the association runner was built to avoid. The new runner therefore imports the association
runner's loader rather than reimplementing it, so the boundary conversion cannot drift between
the two analyses. Each arm-by-chromosome shard writes its own manifest and is skipped when
complete, so an interrupted run resumes.

Parameters match the previous run: L = 10, coverage 0.95, minimum absolute correlation 0.5,
100 iterations, tolerance 1e-3, in-sample MAF 0.05, one-megabase window. Only genes the same
arm called as eGenes are fine-mapped — roughly three thousand per arm.

*The run.* 92 of 92 shards, **12,189 gene fine-mappings** across the four arms, 40 GPU-minutes
on two devices. Run roots `runs/susie_finemapping_v4_k35_20260817` and
`runs/credible_set_comparison_20260817`.

| Arm | Genes with credible sets | Sets per gene | Genes with more than one | Median set size |
|---|---|---|---|---|
| linear · GRCh38 | 2,602 | 1.113 | 283 | 7 |
| graph · GRCh38 | 2,581 | 1.110 | 277 | 7 |
| linear · T2T | 2,562 | 1.121 | 296 | 7 |
| graph · T2T | 2,540 | 1.122 | 300 | 7 |

*The comparison*, restricted to genes fine-mapped in both arms and with variants compared on
normalised identity:

| Contrast | Genes in both | Same number of sets | Sets overlap | Median Jaccard | **Same top variant** |
|---|---|---|---|---|---|
| T2T − GRCh38 · linear | 2,099 | 92.3% | 96.1% | 0.82 | **60.7%** |
| T2T − GRCh38 · graph | 2,088 | 93.2% | 96.0% | 0.82 | **60.8%** |
| graph − linear · GRCh38 | 2,484 | 97.2% | 98.1% | 0.97 | 78.5% |
| graph − linear · T2T | 2,432 | 96.5% | 97.5% | 0.97 | 77.5% |

**Fine-mapping survives the method change at the level it actually reports.** Credible sets
overlap for 96–98% of genes and the median Jaccard is 0.82 under a reference swap. But the
single highest-posterior variant *within* those overlapping sets differs **four times in ten**
under a reference swap.

That is the same phenomenon the lead-variant and LD stages found, seen through the instrument
that matters most: the arms agree about the *region in contention* and disagree about which
member of it to name. The practical reading is that a credible set is a robust thing to report
across method changes and a single "causal variant" is not.

*Bounds.* Credible-set counts are conditional on each arm's own eGene selection and are only
comparable on the shared gene set used here. Overlap is computed on variants placeable in the
common frame, which is 99.2–100% of them, so a set containing reference-unique variants is
compared on the remainder and disagreement is understated.

**Measure reference-allele orientation across the frame.** `[COMPLETE]`
Read directly from `runs/four_arm_variant_identity_comparison_20260811`; no new run root.

Normalising a GRCh38 variant into the common T2T frame swaps which allele is *reference* for
**30.08%** of them (2,661,900 of 8,848,118 in the linear arm; 30.07% in the graph arm). It is
exactly **zero** in both T2T arms, which are already in the frame. Agreement across two
aligners and disappearance on the native reference establish this as a property of the two
genomes, not of the pipeline. LiftoverIndel's own INFO field defines the flag as "REF/ALT were
flipped during liftover. GTs were altered accordingly", so genotypes are recoded and dosages
stay valid.

The magnitude is expected rather than anomalous: CHM13 is a single haploid genome, so at a
common polymorphism it carries the GRCh38-alternate allele about as often as that allele is
frequent, and the variant set here is MAF-filtered.

**"The reference allele" is a property of a choice, not of a variant**, and the choice changes
for roughly a third of common variation. Any quantity polarised on REF — an allele frequency, a
burden direction, an effect sign — carries that choice. Set against it, only **0.38%** of the
same variants needed haplotype realignment to reconcile indel representation. The correction
the field worries about is small; the one it does not discuss is large.

*Practical consequence for this project:* matching GWAS variants to credible-set variants must
compare **unordered** allele pairs. An ordered `(REF, ALT)` join would silently discard close to
a third of genuine matches, and would do so non-randomly.

---

**Place GWAS variants in the common frame, and check the placement.** `[COMPLETE]`
Run root `runs/gwas_variant_placement_crosscheck_20260817`.

The earlier colocalisation placed catalog rsIDs through a T2T-native GWAS VCF released in
March 2022. **That file is a frozen slice of the catalog, not a general rsID map** — 186,904
records in total — so joining on it silently intersected the analysis with a four-year-old
release. The resulting 39.9% placement rate was
recorded as a coverage limit of T2T. It was a limit of the join. The catalog carries its own
GRCh38 coordinate for 99.7% of its genome-wide significant rsIDs.

Two placements were built and compared:

| Placement | rsIDs placed of 341,981 | Share |
|---|---|---|
| GWAS VCF, r2022-03-08 (the source being replaced) | 136,306 | 39.9% |
| dbSNP155 T2T-native (no liftover at all) | 331,862 | 97.0% |
| **LiftoverIndel from the catalog's GRCh38 coordinates** | **334,179** | **97.7%** |

The LiftoverIndel pass uses the same tool, chain, reference-differences VCF and target FASTA
recorded in `grch38_to_t2t_full_variant_liftover_v4_authorization_20260811.json`, so GWAS and
eQTL variants are placed by one identical procedure instead of by two liftovers assumed to
agree. Alleles come from dbSNP155 GRCh38.p13, the build the catalog reports against. Of 570,685
input records, 564,061 lifted (98.8%); 6,552 were unliftable, 70 had multiple overlaps and 2
failed a reference-sequence check.

**The two placements are independent — one is a native build, the other a liftover — and they
agree on the exact base for 99.91% of the
331,177 rsIDs both can place.** The common frame is therefore settled
empirically rather than assumed, which matters because the assumption is weakest in exactly
the divergent regions this project studies.

*Was the 2022 source's dropout random?* Only mildly non-random, and **less so than predicted**.
rsIDs it missed sit in an FDR-flagged discordance window 16.3%
of the time against 15.1% for those it kept — an odds ratio of
1.09. Mean segmental-duplication content is
0.0411 against 0.0337 and mappability is essentially
identical. The prediction that dropout concentrated in divergent regions is weakly supported at
best; the real problem with the old source was that it lost 60% of the data outright.

---

**Two tooling notes for anyone repeating this.**

1. LiftoverIndel's REF/ALT flip path calls `var.genotype.array()` to recode genotypes and
   therefore **crashes on a sites-only VCF**. One placeholder sample carrying `0/0` satisfies
   it; the flip rewrites that genotype and only coordinates and alleles are read back.
2. `--chrom` restricts the target-FASTA load but **not** the input iteration, so a restricted
   run walks into other contigs and fails with `KeyError` on the reference-differences index.
   Subset the input VCF instead.
3. The tool emits records in source order and warns that the output "still requires indel
   normalizing, sorting"; `bcftools index` fails until it is sorted.

---

**Colocalise every cis-eQTL against a GWAS panel, Bayesian.** `[COMPLETE]`
Run root `runs/gwas_coloc_bayesian_20260817`. **Supersedes `runs/gwas_coloc_20260817` and
`runs/gwas_coloc_v2_20260817`**, which scored credible-set membership rather than
colocalisation.

*Why the earlier runs were replaced.* They asked whether a genome-wide significant GWAS
variant fell inside a gene's SuSiE credible set. That is variant overlap: no model of the GWAS
signal, the answer hinging on one variant's membership, and the variant in question being
whichever index SNP a study happened to report — itself a tagging choice, and therefore
entangled with what was being measured. It also inherits the GWAS Catalog's ascertainment.

*Method.* `coloc.abf` (Giambartolomei 2014). Priors p1 = p2 = 1e-4, p12 = 1e-5; prior effect
variances 0.15² (quantitative expression) and 0.2² (log odds). A gene-trait pair is tested
where the GWAS has p ≤ 1e-5 inside the gene's cis window and ≥ 50 variants are shared; a call
is PP4 ≥ 0.8. Implemented vectorised in polars/numpy with log Bayes factors accumulated by
log-sum-exp in float64, and **validated against the R `coloc` package on 40 random regions
including planted shared signals: largest posterior difference 2.0e-9**.

*The single-causal-variant assumption is retained deliberately.* `coloc.susie` relaxes it but
requires an LD reference matched to every GWAS, and LD-reference mismatch is its best-known
failure mode. The question here is whether four arms differ, so an LD-free method is the safer
instrument: the assumption applies identically to all four and cancels from the comparison
even where it is wrong in absolute terms.

*eQTL input* is `runs/four_arm_nominal_v4_k35_nanfix_20260811` — every cis gene-variant pair,
not only fine-mapped eGenes. Both effect vectors are re-signed to the common frame's alternate
allele before any Bayes factor is computed, which is required because REF/ALT swap for ~30% of
variants across that boundary.

*Panel:* 40 studies —
11 comparison, 10 immune, 2 neurodegenerative, 5 neurodevelopmental, 12 psychiatric. The well-powered
psychiatric and neurodevelopmental GWAS **have no summary statistics in the GWAS Catalog**
(PGC3 schizophrenia and the large bipolar meta-analyses return HTTP 404 there); they come from
the consortium's public figshare deposits and are mapped into the common frame through a
T2T-native dbSNP index, so no liftover is applied to them either.

| Arm | Tests | Genes tested | Colocalisations | Genes with one |
|---|---|---|---|---|
| linear_grch38_dv | 259,239 | 17,936 | 926 | 672 |
| graph_grch38_dv | 259,257 | 17,937 | 930 | 673 |
| linear_t2t_dv | 256,463 | 17,745 | 933 | 664 |
| graph_t2t_dv | 256,516 | 17,745 | 928 | 657 |

| Contrast | Both | Only first | Only second | Calls changed | Confident changed | Median \|ΔPP4\| | r |
|---|---|---|---|---|---|---|---|
| T2T − GRCh38 · linear | 730 | 168 | 172 | **31.8%** | **12.2%** | 0.00215 | 0.915 |
| T2T − GRCh38 · graph | 730 | 164 | 175 | **31.7%** | **12.1%** | 0.00215 | 0.914 |
| graph − linear · GRCh38 | 909 | 21 | 17 | **4.0%** | **1.6%** | 0.00005 | 0.995 |
| graph − linear · T2T | 913 | 15 | 20 | **3.7%** | **1.4%** | 0.00005 | 0.996 |

**A reference swap changes about a third of colocalisation calls and an aligner swap about one
in twenty-five — but the calls that move are the marginal ones.** Restricted to pairs already
confident in one arm, the figures are 12.2% and 1.6%. The evidence barely moves: median |ΔPP4|
is 0.0022 and 0.00005, with r = 0.915 and 0.995.

*Against the superseded proxy:*

| Contrast | Credible-set membership | coloc.abf | Confident coloc.abf |
|---|---|---|---|
| T2T − GRCh38 · linear | 56.6% | 31.8% | **12.2%** |
| T2T − GRCh38 · graph | 53.4% | 31.7% | **12.1%** |
| graph − linear · GRCh38 | 21.0% | 4.0% | **1.6%** |
| graph − linear · T2T | 22.8% | 3.7% | **1.4%** |

**The proxy overstated instability by two to five times.** It is biased upward, not merely
noisier: a conjunction hinging on one variant flips far more readily than a posterior computed
over a region.

**The direction confound dissolves.** The four arms give 926–933 colocalisations, a spread
under 1%. The proxy showed GRCh38 arms finding ~22% more than T2T arms, which the page
declined to interpret on ascertainment grounds. That caution was correct: with full summary
statistics the asymmetry disappears entirely.

*By trait area, and the thesis prediction fails here:*

| Trait area | Studies | Pairs tested | Called | Reference swap | Aligner swap |
|---|---|---|---|---|---|
| immune | 10 | 63,021 | 236 | **36.9%** | 8.7% |
| neurodegenerative | 2 | 4,386 | 17 | **35.3%** | 0.0% |
| comparison | 11 | 117,144 | 472 | **32.8%** | 3.1% |
| psychiatric | 12 | 52,503 | 257 | **28.4%** | 2.1% |
| neurodevelopmental | 5 | 13,356 | 88 | **21.6%** | 2.5% |

**Brain traits change least, not most.** Pooling psychiatric with neurodevelopmental gives
26.7% against 34.2% for the rest — **odds ratio 0.70, p = 0.014** (Fisher p = 0.014). Two
rescues were tested and neither holds:

- *Not a source artefact.* Psychiatric is mostly consortium-sourced and immune entirely
  Catalog-sourced, so source and trait area are nearly collinear. Within psychiatric, the two
  sources give 27.8% and 28.6%, **p = 1**. Restricting to Catalog-only leaves the ordering
  unchanged. Recorded in `SOURCE_CONFOUND.json`.
- *Not because brain calls are better determined.* Neurodevelopmental has the **lowest**
  median PP4 among its calls (0.873 against 0.910 immune) and the fewest above 0.95, so it
  sits closer to the threshold and should flip more.

*What this bounds.* It does not overturn the regional results, which are measured directly on
association statistics rather than on a downstream conjunction. It establishes something
narrower: **among loci where a well-powered GWAS has already resolved a clean signal,
brain-trait colocalisations are not more fragile than others.** A gene-trait pair is only
testable where the GWAS already has an association in the window, so if representation choice
bites hardest where GWAS has *not* yet produced clean signals, those loci are excluded by
construction. That is testable and is not tested here.

*Resources.* 40 traits, 1,031,475 tests, peak 10.0 GB RSS on 24 threads,
18 minutes. The first attempt joined all traits against all gene-variant pairs at once and
reached 159 GB on 65 cores; it was killed and restructured to one trait at a time with an
overlap prefilter.

---

**Fine-mapping settings, recorded because they decide what a credible set means.**
SuSiE via `tensorqtl.susie`, per arm, over that arm's own eGenes at
5% FDR with 35 expression PCs:
L = 10, coverage 95%, purity `min_abs_corr` = 0.5,
MAF ≥ 0.05, ±1 Mb window, tolerance 0.001 within
100 iterations.

- **L = 10** caps independent signals per gene; a gene at the cap is censored rather
  than resolved, so sets-per-gene is a between-arm comparison and not a biological quantity.
- **Coverage 95%** decides set size. Raising it would raise the overlap
  statistics and not the top-variant agreement, so the two are not interchangeable.
- **`min_abs_corr` = 0.5** is the purity filter, removing sets that are the
  algorithm splitting one signal rather than a real second one.
- Missing dosages reach tensorqtl as the `-9` sentinel, not IEEE NaN, via the same boundary
  adapter the association runs use; without it variants are dropped silently.

**Is the colocalisation negative a selection effect?** `[COMPLETE — three tests]`
Run roots `runs/gwas_coloc_bayesian_20260817` (gate sweep and window stratification,
`GATE_AND_WINDOW_STRATIFICATION.json`, `GWAS_COVERAGE_IN_FLAGGED_WINDOWS.json`) and
`runs/neuro_gene_sets_20260817`.

The colocalisation comparison found brain traits *less* sensitive to a reference swap than
others, against prediction. Because a pair is only testable where a GWAS already has a signal
in the window, that null is ambiguous between "these regions are ordinary" and "the instrument
cannot see into them". Three tests separate those.

---

### 1. GWAS coverage inside the flagged windows

| Trait group | Studies | Coverage inside flagged windows ÷ outside |
|---|---|---|
| immune | 10 | 0.870 |
| comparison | 11 | 0.862 |
| neurodegenerative | 2 | 0.798 |
| psychiatric | 12 | 0.777 |
| neurodevelopmental | 5 | 0.765 |

**Every trait is thinner inside the discordance windows, and brain traits are the thinnest.**
But the standard-error ratio at matched allele frequency is ≈1.000: conditional on a variant
being measured there is no power loss. The studies do not lose precision where they look —
**they fail to look**. That is a coverage artefact, not quiet biology.

*Unresolved confounds.* Flagged windows are segmental-duplication enriched, so some depletion
is expected for any technology; and the psychiatric studies are mostly older consortium
releases, so the trait difference could be an imputation-panel vintage effect rather than a
trait effect.

### 2. Stratifying by window — **the regional claim survives, and the genome-wide test was diluting it**

The genome-wide comparison averages over the ~80% of tested genes that sit outside a flagged
window. Splitting on the gene's anchor:

| Contrast | Inside flagged windows | Outside | Odds ratio | p |
|---|---|---|---|---|
| reference (T2T − GRCh38, linear) | 109/294 = **37.1%** | 233/778 = 29.9% | 1.38 | 0.028 |
| aligner (graph − linear, GRCh38) | 21/261 = **8.0%** | 17/686 = 2.5% | **3.44** | 0.00028 |

**Colocalisation calls are significantly more likely to change inside the windows this project
flagged from eQTL statistics alone** — 3.4× more likely on the aligner axis. This is the
project's regional claim reproducing at the endpoint that matters, on an instrument that
shares no statistics with the one that defined the windows.

**But the trait-specific claim does not follow.** Within flagged windows brain traits are
still no more affected than others (OR 0.73–0.74, p = 0.18–0.26, not significant), and outside
them they are significantly *less* affected (OR 0.65–0.67, p ≈ 0.017).

### 3. Loosening the testing gate

The gate sweep scored every pair reaching p ≤ 1e-3 and recorded each one's strongest GWAS
p-value, so stricter gates are recovered by filtering. Brain-versus-other odds ratio on the
reference axis:

| Gate | 5e-8 | 1e-5 | 1e-4 | 1e-3 |
|---|---|---|---|---|
| odds ratio | 0.56 | 0.69 | 0.71 | 0.71 |
| p | 0.004 | 0.012 | 0.013 | 0.011 |

**The deficit narrows as weaker GWAS signals are admitted — in the direction the
underpowered-GWAS explanation predicts — but it never reverses.** Partial support at best.

### 4. Sequencing-derived gene sets — **retracted; the instrument carries the same bias**

This test asked whether the flagged windows are enriched for curated disease genes, on the
reasoning that gene sets built from sequencing rather than from arrays would be immune to the
coverage depletion that limits the colocalisation test. **That reasoning was wrong.** SFARI,
DDG2P and gnomAD are all built from short reads aligned to GRCh38 with linear aligners — the
same procedure whose blind spots this project exists to map. Exome capture probes are designed
against GRCh38 and perform poorly in duplicated or divergent sequence, multi-mapping reads are
discarded, and a gene whose variants cannot be called reliably never accumulates the evidence
to be curated as a disease gene at all. The control inherits the confound it was meant to break.

The raw result was a mild depletion of disease genes in flagged windows — SFARI OR 0.84
(p = 0.026), DDG2P OR 0.80 (p = 1.7e-4), and a higher median gnomAD LOEUF (0.906 against
0.868, p = 6.2e-5) implying *less* constrained genes.

**All of it is a power artefact, and the arithmetic shows it exactly.** LOEUF is the upper
bound of a confidence interval on the observed-to-expected ratio of loss-of-function variants,
and a confidence bound widens when there is less information. Inside flagged windows:

| Quantity | Flagged | Unflagged | Ratio | p |
|---|---|---|---|---|
| possible LoF sites (sequence only) | 180.0 | 200.0 | 0.900 | 3.7e-7 |
| expected LoF variants | 36.2 | 40.3 | 0.900 | 1.2e-6 |
| observed LoF variants | 20.0 | 21.0 | 0.952 | 4.3e-4 |
| confidence-interval width | 0.419 | 0.393 | **1.066** | 5.9e-6 |
| LOEUF (the upper bound) | 0.906 | 0.868 | 1.044 | 6.2e-5 |
| **o/e point estimate** | 0.620 | 0.603 | 1.029 | **0.052 (n.s.)** |

The point estimate of constraint does not differ significantly; only the *bound* does, and the
interval is 6.6% wider. Genes in flagged windows offer **10% fewer possible loss-of-function
sites** — that is the reference and its annotation, not biology.

Stratifying genes by expected LoF count, which matches them on how much information gnomAD
has, **the LOEUF gap disappears in every quintile**:

| Expected LoF quintile | ~11 | ~25 | ~40 | ~61 | ~112 |
|---|---|---|---|---|---|
| p | 0.454 | 0.469 | 0.23 | 0.131 | 0.789 |

The DDG2P depletion largely dissolves under the same stratification (odds ratios 0.88, 1.29,
0.71, 0.81, 1.13; only one quintile nominally significant, and two above 1).

**Nothing in this test bears on the hypothesis.** It measures where short-read GRCh38 pipelines
can see, which is what the flagged windows already encode. Recorded in
`runs/neuro_gene_sets_20260817/CONSTRAINT_POWER_ARTEFACT.json`.

*What would be a valid instrument.* One whose ascertainment does not pass through short-read
GRCh38 alignment: assembly-based benchmarks (the GIAB Q100 truth set and the Challenging
Medically Relevant Genes benchmark, both available on CHM13v2.0 locally), cytogenetically
ascertained disease regions, or long-read studies. It is worth noting that the project's own
ClinGen result — flagged genes enriched 2.9–9.1× for recurrent genomic-disorder regions — comes
from exactly such a route, since those regions were defined by karyotype and chromosomal
microarray rather than by sequencing.

### What the three tests establish together

- **The regional claim is confirmed at the colocalisation endpoint** and was being hidden by
  genome-wide dilution. This is a new positive result.
- **The trait-specific claim is untestable with these data, not refuted.** The apparent
  brain-trait deficit survives source, confidence and coverage checks but not adjustment for
  GWAS signal strength: brain studies are about three orders of magnitude weaker inside a
  gene's window (median 10⁻⁶·² against 10⁻⁹·⁴), and matched on signal the difference is not
  significant (OR 0.80, 95% CI 0.59–1.08, p = 0.15; joint with coverage 0.78, p = 0.12). The
  first three checks all left the deficit standing, so stopping there would have produced a
  confident published negative that the fourth removes.

The honest summary is that *where* a method change matters is now established on two
independent instruments, and *which diseases it matters for* is not established at all.

---

**Decompose the reference contrast into genotype and phenotype terms.** `[COMPLETE]`
Run root `runs/crossed_reference_association_20260821`.

*Why.* Expression is quantified once per **reference**, not once per arm. An aligner swap
therefore holds the phenotype exactly fixed while a reference swap changes the genotypes, the
RNA alignment and annotation, and the expression covariates together. The published reference
contrast cannot say which of those moves the result. This runs the missing cell.

| Cell | genotypes | expression |
|---|---|---|
| published `linear_grch38_dv` | GRCh38 | GRCh38 |
| **crossed** (new) | GRCh38 | T2T |
| **control** (new) | T2T | T2T |
| published `linear_t2t_dv` | T2T | T2T |

Both new cells run through identical code, so control-versus-crossed differs *only* in the
genotype matrix. GRCh38 variants are relabelled into the common frame and dosages recoded to
2 − d on the 30.60% that swap reference and alternate;
8,700,370 of 8,848,118 GRCh38 variants carry a unique common identity.

*The gate.* Control against the published `linear_t2t_dv` — the same cell by two code paths,
differing only in a seed derived from the cell name:

| | genes | eGenes | turnover | r | top variant |
|---|---|---|---|---|---|
| control vs published | 18,114 | 3,043 vs 3,041 | **0.9%** | **0.9999** | **100.0%** |

The relabelling machinery reproduces production exactly, so the crossed number is interpretable.

*The result.*

| Contrast | differs in | eGene turnover | r | top variant |
|---|---|---|---|---|
| control vs crossed | **genotypes only** | **6.6%** | 0.9927 | 72.2% |
| crossed vs published GRCh38 | **expression only** | **20.1%** | 0.9409 | — |
| published T2T vs GRCh38 | both, as published | **21.7%** | 0.9380 | — |

**The published reference contrast is dominated by RNA re-quantification, not by genotype
representation.** Swapping expression alone reproduces 20.1%
of the 21.7% confounded turnover; swapping
genotypes alone gives 6.6%. The correlations agree: a genotype
swap barely moves gene-level significance
(r = 0.9927) while an expression swap accounts
for essentially the whole drop
(0.9409 against the confounded
0.9380).

*How this sits with the magnitude result.* The expression-stratification analysis
(`runs/expression_confound_decomposition_20260819`) found the roughly nine-fold
reference-to-aligner ratio in median |Δz| survives holding expression stable. Both are true and
they measure different things: that is variant-level effect *magnitude*, this is gene-level
call *identity*. The gene-class and disease enrichments are computed on call identity, which is
the quantity dominated by expression.

*Bounds.* eGene turnover counts calls crossing a threshold rather than measuring how far
effects move, and is conditional on genes in the T2T annotation. The crossed cell is not a
biologically meaningful pipeline — nobody calls variants on one reference and quantifies RNA on
another — and exists only to separate the two terms. Top-variant agreement is interpretable
only for the genotype term; it reads 0.0% elsewhere because the two sides use different
identifier coordinate systems, which is a formatting artefact.

*What still needs doing.* The gene-class and disease enrichments should be recomputed on the
genotype term alone. Whichever survive are attributable to variant representation; whichever do
not belong to RNA quantification and must be described that way.

---

**Which mechanism carries the disease-region enrichment?** `[COMPLETE]`
Run roots `runs/crossed_reference_nominal_20260821` and
`runs/genotype_term_gene_classes_20260821`.

The crossed association showed the reference contrast is mostly RNA re-quantification by
volume. That raised an obvious worry about the project's strongest disease claim: the
gene-class enrichments are computed on genes flagged by the *reference contrast*, so they might
be a property of where RNA quantification disagrees rather than of where variant representation
does. This tests it by recomputing them on the genotype term alone — control (T2T genotypes)
against crossed (GRCh38 genotypes), with expression, covariates and annotation held identical.

Nominal cis surfaces were run for both cells (3.3 min each on GPU) so the published gene-level
discordance statistic could be reproduced exactly, then genes were flagged by the same recipe
and passed through the same class-assignment functions rather than a reimplementation.

| Class | Published reference axis | **Genotype term** | genes in class | flagged in vs out |
|---|---|---|---|---|
| genomic disorder | 2.91 (p = 1e-23) | **4.02** (p = 4e-24) | 717 | 12.8% vs 3.5% |
| human accelerated | 1.84 (p = 1e-08) | **1.54** (p = 0.007) | 875 | 5.8% vs 3.9% |
| segmental duplication | 1.65 (p = 2e-16) | **2.40** (p = 3e-24) | 4,413 | 6.6% vs 2.9% |

**The enrichments are not an artefact of RNA re-quantification. Two of the three are stronger
on the genotype term.** Isolating variant representation raises genomic disorder from
2.91x to 4.02x and segmental
duplication from 1.65x to 2.40x;
human accelerated regions weaken from 1.84x to
1.54x but hold.

*Why the two results fit together.* The genotype term flags a **smaller** set —
579 genes (4.0% of
14,507) against the published 1,256
(8.8%) — but a **more concentrated** one:
12.8% of genomic-disorder genes are flagged against 3.5%
elsewhere, a ratio the reference contrast does not reach. Variant representation differs
specifically in duplicated and divergent sequence, which is where recurrent-CNV regions sit.
RNA quantification differs more broadly and less selectively. **The reference effect is mostly
RNA by volume and genotype by mechanism** — separable claims with different remedies:
annotation and RNA alignment for the first, genotyping reference for the second.

*Bounds.* The two flagged sets differ in size and use different z-cuts (genotype term
z > 2.88), so "discordant" does not mean the same thing on both sides and the
odds ratios are not strictly commensurable. The genotype term is measured only on variants both
references can represent, so it excludes the access effect and is a **lower bound**. Neither
flagged set is a controlled discovery set — the normal-tail false-discovery estimate is
uncalibrated in both, as this project has already shown. Enrichment is unmatched on
mappability, gene density and expression level, so it shows disease genes are present in these
regions rather than why.

---

**Validating the dosage re-signing.** `[COMPLETE]`
Recorded in `runs/crossed_reference_nominal_20260821/FLIP_VALIDATION.json`.

The gate for the crossed association compared the control cell against the published T2T arm,
which exercises the relabelling code only where it does nothing — T2T variants have zero flips.
It said nothing about the 2,661,896 GRCh38 dosages recoded as 2 − d, and a sign error there
would manufacture a genotype effect indistinguishable from a real one.

Allele frequency settles it. Both cells report the frequency of the common-frame *alternate*
allele, so re-signed and untouched variants must agree equally well; an inverted sign would
show re-signed variants matching 1 − af instead. On chr22:

| Variants | r(af control, af crossed) | median abs diff | had the sign been inverted |
|---|---|---|---|
| re-signed (34,284) | 0.9351 | **0.0000** | r = -0.9351, median 0.3911 |
| untouched (79,272) | 0.9990 | **0.0000** | r = -0.9990, median 0.6933 |

The re-signing is correct. One observation left open: re-signed variants correlate at
0.9351 against 0.9990 for untouched ones,
despite identical median absolute agreement. The benign reading is a restricted-range effect,
since flipped variants have a different frequency distribution — visible in the inverted-sign
column, where their median distance to 1 − af is 0.39
against 0.69. The substantive reading is that
sites where the two references disagree about which allele is *reference* are also sites where
genotyping disagrees more — the project's thesis in miniature. Not resolved here.

---

**Retraction (internal).** An earlier note in this session claimed the gene-class analysis
intersected hg38 ClinGen intervals against T2T gene positions without lifting, and called it a
critical flaw in the project's strongest disease claim. **That was wrong.**
`analyze_discordant_gene_classes_20260816.py` loads both annotations and matches each track to
its own frame: GRCh38 gene positions against the hg38 ClinGen track, CHM13 positions against
the CHM13 HAR and segdup tracks. No liftover is needed because nothing crosses frames. The
error came from reading the manifest's track paths and grepping for "lift" rather than reading
how the tracks were used. The mismatch existed only in a first draft of the genotype-term
script, which parsed the T2T annotation alone; that draft was replaced with one that delegates
class assignment to the original functions.

---

**Is the enrichment about the class, or about everything that travels with it?**
`[COMPLETE]`
Run root `runs/matched_gene_class_enrichment_20260821`.

The enrichments compare flagged genes against every gene tested, matched on nothing. Genes in
these classes carry more cis variants, are longer, sit in denser neighbourhoods, are expressed
differently and are harder to align -- any of which could produce the enrichment with no
contribution from the class itself. This holds them fixed, in two tiers, with three estimators
that fail differently: a logistic regression of flagged on class plus covariates and their
squares; a propensity for class membership cut into 20 strata and
pooled by Mantel-Haenszel; and a 1:1 nearest-neighbour match with a
0.2 SD caliper, tested
by McNemar. Balance is reported as standardised mean differences against a 0.10 target.

*Tier one -- statistical power.* Cis-variant count, gene length, gene density and expression
level. These change how many chances a gene has to look discordant and how precisely each is
measured. Nothing about them is biological.

| Axis | Class | Crude | Regression | Stratified | Matched (95% CI) |
|---|---|---|---|---|---|
| genotype term | recurrent genomic-disorder region | 4.02 | 3.57 | 3.50 | 5.87 (3.39-10.14) |
| genotype term | segmental duplication | 2.40 | 2.94 | 2.78 | 2.68 (2.14-3.37) |
| genotype term | human accelerated region | 1.54 | 1.61 | 1.61 | 1.39 (0.90-2.13) |
| reference axis | recurrent genomic-disorder region | 2.91 | 2.83 | 2.81 | 2.77 (1.98-3.86) |
| reference axis | segmental duplication | 1.65 | 1.80 | 1.78 | 1.74 (1.49-2.03) |
| reference axis | human accelerated region | 1.84 | 1.96 | 1.90 | 1.64 (1.22-2.21) |

**None of the enrichment is a counting artefact.** Segmental duplication and human
accelerated regions come out *stronger* than crude on both axes, meaning the power covariates
were masking them rather than manufacturing them. Only genomic disorder attenuates, and barely
-- 4.02x to 3.57x on the genotype term and 2.91x to 2.83x on the reference axis, both still far
from one. Balance is good throughout: the worst standardised mean
difference on this tier is
0.028.

*Tier two -- adding mappability.* This is **not** a stricter version of tier one, and the
manifest records why per class. A segmental duplication is by definition sequence that recurs
elsewhere, which is what makes it unmappable; conditioning on mappability there removes the
mechanism rather than a confounder. Recurrent-CNV regions are bounded by segdups, so the same
applies in weaker form. Only for human accelerated regions -- defined by substitution rate
rather than copy structure -- is mappability cleanly a confounder.

| Axis | Class | Crude | Regression | Stratified | Matched (95% CI) |
|---|---|---|---|---|---|
| genotype term | recurrent genomic-disorder region | 4.02 | 1.37 | 1.81 (unbalanced) | 1.17 (0.81-1.68) |
| genotype term | segmental duplication | 2.40 | 2.06 | 1.91 | 1.87 (1.51-2.32) |
| genotype term | human accelerated region | 1.54 | 1.15 | 1.10 | 0.92 (0.62-1.37) |
| reference axis | recurrent genomic-disorder region | 2.91 | 1.31 | 1.70 (unbalanced) | 1.32 (0.96-1.80) |
| reference axis | segmental duplication | 1.65 | 1.39 | 1.37 | 1.39 (1.19-1.62) |
| reference axis | human accelerated region | 1.84 | 1.61 | 1.47 | 1.54 (1.15-2.07) |

Three things follow.

1. **Segmental duplication survives conditioning on its own mechanism**, on both axes and all
   three estimators, with balance at or under 0.03. Genotype term
   2.06x,
   reference axis
   1.39x.
   The segdup annotation carries information that 100-mer uniqueness over the gene span does
   not -- unsurprising, since the discordance statistic is computed over a 1 Mb window while
   mappability here is measured over the gene body.

2. **The human-accelerated enrichment on the genotype term is entirely mappability.** It falls
   from 1.54x to
   1.15x
   (p = 0.41), and
   the matched estimator puts it at
   0.92. This
   is the tier where mappability is a genuine confounder for this class, so **this one is a
   negative result and is reported as such.** The reference axis behaves differently and holds
   (1.61x,
   p = 2.2e-05).

3. **The genomic-disorder enrichment is mostly carried by mappability.** It falls to
   1.37x on
   the genotype term and
   1.31x on
   the reference axis, with matched intervals that include one. Since mappability is partly
   downstream of the class here, that is not evidence the enrichment was spurious -- but the
   page can no longer claim it is independent of alignment difficulty. On this cell the
   stratified estimator should be ignored: its mappability balance is
   0.207, a
   common-support failure, because genomic-disorder genes are systematically unmappable and
   there are too few comparable controls. The matched estimator solves that by discarding the
   8 cases
   it cannot match rather than extrapolating, and balances to
   0.040.

*What this does and does not settle.* It removes the "these genes just had more chances"
explanation entirely, which was the live one. It replaces the claim that flagged genes are
enriched for disease sequence *independently* of alignment difficulty with the narrower and
better-supported claim that the enrichment **is** alignment difficulty, concentrated in the
duplicated sequence where disease genes sit. That is the mechanism the design proposed, so it
is not a retreat -- but it is a different sentence from the one the page was making, and the
distinction is not rhetorical: the first would survive better alignment, the second is what
better alignment is for. Mappability is still measured over the gene span rather than the cis
window, so the adjustment is partial on the one covariate that matters most.

---

**Does the matched result survive measuring mappability where the statistic is?**
`[COMPLETE]`
Run root `runs/window_mappability_20260821`.

The matched enrichment adjusted for mappability across each gene's *span*, while the statistic
it explains is a mean over that gene's cis variants, which sit anywhere within a megabase of
the transcription start site. Mappability carries most of the attenuation, so measuring it in
the wrong place was the weakest joint in that result. Three definitions were computed and the
whole analysis re-run under each.

| Definition | What it measures | Exactly 1.0 | r against the variant-level definition |
|---|---|---|---|
| gene span | the gene's own footprint; what the matched analysis used | 65% | 0.761 |
| TSS +/- 1 Mb | the nominal cis window, positionally | 0% | 0.439 |
| anchored footprint | lowest to highest anchored variant | 45% | 0.930 |
| at the variants | share of the gene's anchored cis variants in mappable sequence | 65% | 1 |

Transcription start sites were recovered from the nominal output as variant position minus
`start_distance`, which is the coordinate tensorqtl actually used and needs no strand
assumption.

**The nominal cis window is a poor stand-in for where the evidence is**
(r = 0.44), and this was the definition the plan had intended to use
for the reference axis. A 2 Mb window always contains some unmappable sequence, so it is smooth
and well-behaved and measures broad regional mappability rather than the mappability of the
loci that produced the statistic. The gene span does much better
(r = 0.76) and the anchored footprint better still
(r = 0.93).

*Adjusted odds ratios, power and mappability both held fixed, genotype term:*

| Class | Unmatched | Mappability at the gene span | at the TSS +/- 1 Mb | at the variants |
|---|---|---|---|---|
| recurrent genomic-disorder region | 4.02x | 1.37x (p = 0.04) | 1.70x (p = 0.0001) | 1.27x (p = 0.1) |
| segmental duplication | 2.40x | 2.06x (p = 7e-14) | 2.50x (p = 1e-23) | 2.08x (p = 2e-14) |
| human accelerated region | 1.54x | 1.15x (p = 0.4) | 1.38x (p = 0.04) | 1.29x (p = 0.1) |

**The published matched result is confirmed by the faithful definition.** Gene-span and
variant-level mappability give the same answers to within noise: segmental duplication
2.06x against
2.08x, genomic disorder
1.37x against
1.27x, human accelerated
1.15x against
1.29x. Segmental duplication
survives under all three definitions; the other two sit at or near one under both faithful
ones. **The caveat is therefore downgraded rather than merely closed**: the mismatch existed but
did not matter.

*Why the window column is the outlier, and why that is reassuring.* Adjusting for a covariate
measured with error under-adjusts, so a noisy proxy for the confounder leaves an estimate
closer to the crude one. That is the ordering observed: in every row the window column is the
adjustment nearest the crude value and the two faithful definitions sit further from it. It is
not strictly bracketed -- segmental duplication comes out at 2.50x against a crude 2.40x, a
slight overshoot -- but the ranking holds in all three, and the window definition is the one
with r = 0.44. The pattern the correlations predict is the pattern the
odds ratios show.

*Bounds.* The reference axis has no variant-level definition: its run retained only the gene
table, not the anchored variant set, so a variant-level check there would mean reproducing the
common-frame identity join over the full four-arm nominal scan. Its gene-span result is
reported as-is, supported by the genotype-term agreement rather than by a direct check. Both
faithful definitions are also conditioned on variant existence -- variants are called where
reads align, which is why 65% of genes sit at exactly 1.0
on the variant-level scale. And every caveat from the matched analysis still applies:
mappability is a mediator for segmental duplication and partly so for recurrent-CNV regions.

---

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
