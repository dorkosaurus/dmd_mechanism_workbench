# Hypothesis assembly — a scientific walkthrough

> How the workbench gathers evidence and forms mechanism hypotheses
> that connect a patient's pathogenic variant to their clinical phenotype.

## The scientific question

Given a rare-disease patient with a known pathogenic variant, we want to
answer: **how does this specific variant produce the clinical signs and
symptoms we observe in this specific patient?**

The answer is a *mechanism* — a causal chain of biological events that
begins with the variant and terminates in a phenotype. The mechanism is
the *why* connecting genotype to clinic. Every therapeutic decision
implicitly rests on one.

The mechanism itself is not directly observable. What we have is
evidence at various points along the chain — sequence data, protein
predictions, pathway memberships, expression atlases, clinical labs. The
job of the workbench is to **assemble that scattered evidence into a
coherent chain**, then score how well the assembled chain holds
together. That assembled, evidence-linked chain is what we call a
*hypothesis*. Multiple hypotheses coexist per patient; ranking them is
a separate step.

## The scaffold: biological organization

Rare-disease biology is naturally hierarchical. A variant does not
cause a phenotype in one step — it causes protein changes, which
disrupt pathways, which affect subcellular structures, which impair
specific cell types, which weaken tissues, which manifest as
observable signs. Each level is a scale of organization at which the
consequences of the variant can be *observed and reasoned about
independently*.

We adopt the following seven-layer hierarchy as the assembly scaffold:

```
   variant → protein → pathway → subcellular → cell type → tissue → phenotype
```

Each layer is a **node** in a directed acyclic graph. Adjacent layers
are connected by **edges** that represent causal transitions:
*"variant produces this protein change"*, *"this protein change
disrupts this pathway"*, and so on.

Every hypothesis is a **walk through this graph** — a specific chain
of node-claims and edge-claims that together tell a mechanism story
from variant to phenotype.

## Formation vs. ranking

We treat these as two distinct operations:

**Formation** = gathering the evidence for each node and each edge in
the chain. This is empirical work: consulting data sources and running
predictive models to accumulate the premises that inform each step.
Formation happens per patient because the patient's context (variant,
labs, phenotype) determines which evidence applies.

**Ranking** = comparing whole assembled chains to decide which one best
explains the patient's observations. Ranking is a downstream operation
on the finished chains; it does not itself gather evidence.

Keeping these separate matters because:

- A weakly-ranked hypothesis can still have strong evidence at
  specific steps — worth surfacing in the trace even if not the winner.
- Adding a new data source affects formation (more evidence available);
  ranking then re-scores over the enriched set without changing rules.
- We can compare hypotheses cell-by-cell along the chain, not just
  by a single scalar score.

## Two kinds of evidence: node premises and edge premises

At every position along the chain, evidence takes one of two forms:

**Node evidence** — what we know about the biological state at a
single layer of organization. Examples: the variant record itself
(node: variant); the ESM3-predicted 3D structure of the affected
isoform (node: protein); the tissue-scale imaging of muscle
composition (node: tissue).

**Edge evidence** — what we know about the causal transition from one
layer to the next. Examples: NMD prediction linking a variant to
transcript loss (edge: variant → protein); pathway-component tissue
expression linking a pathway disruption to cell-type dysfunction
(edge: pathway → subcellular); functional labs linking tissue
dysfunction to a clinical measurement (edge: tissue → phenotype).

Both kinds matter. A chain with strong node evidence but weak edges
is a chain with unexplained gaps — the layers are populated but
their causal relationship is not established. A chain with strong
edges but empty nodes has coupling without content. The formation
process aims to accumulate both.

## Data sources and which layer each informs

The workbench currently integrates seven data or model sources. Each
source produces typed evidence rows (**premises**) that are then
attributed to node and edge positions in the biological chain:

| Source | Type | Informs | Nature of the evidence |
|---|---|---|---|
| **Zhang 2024 patient cohort** | data | variant (node), phenotype (node) | Per-patient variant record + clinical phenotype label from the published supplementary tables (N=418) |
| **ClinVar NMD cross-tab** | model | variant→protein (edge) | Predicted NMD classification (triggering / escaping / transcript-dependent) applied to 11,790 ClinVar DMD variants; establishes cohort-level base rates for variant → transcript-loss transitions |
| **DMD isoform architecture** | data | variant→protein (edge), protein (node) | Which of the seven DMD isoforms is impacted given the variant's exon position (first-shared-exon rule); a structural argument for the variant → protein transition |
| **ESM3 protein fold** | model | protein (node), variant→protein (edge) | Per-patient ESM3-predicted 3D structure of the wild-type and truncated protein product; visualizes what protein content is preserved and what is lost |
| **Reactome pathway memberships** | data | pathway (node), pathway→subcellular (edge) | Which biological pathways include dystrophin; a claim about downstream molecular consequences of protein loss |
| **Human Protein Atlas cell-type expression** | data | cell type (node), cellType→tissue (edge) | Which cell types express which DMD isoforms; a claim about which cells will be functionally affected by isoform-specific loss |
| **Per-patient clinical labs** | model / data | cell type (node), tissue (node), phenotype (node) | Fifteen assays per patient (CK, LVEF, FVC, MRI fat fraction, 6MWT, NSAA, IQ, ERG, UACR, etc.); each lab occupies its natural layer of the chain and provides the concrete empirical anchor for that node |

Additional sources — AbSplice for aberrant-splicing predictions,
ESM2 log-likelihood ratios for variant-level functional-disruption
scoring, and downstream network sources (Open Targets, STRING, ARCHS4,
Geneformer) — are called out in the rare-disease-platform deck as
priority additions. Each would enter as a new premise producer,
attributing to the appropriate node or edge without disturbing sources
already integrated.

## How a hypothesis emerges

For each patient, hypothesis formation proceeds in four moves:

**1. Assemble the patient's premise bundle.** All patient-scoped
premises (labs, variant record, isoform impact, protein fold) plus all
gene-scoped or cohort-scoped premises (HPA expression, Reactome
memberships, ClinVar NMD cross-tab) are gathered into a single set.
For Patient 5 in our current roster this is approximately ten patient
premises + three cohort premises = thirteen total.

**2. Instantiate the candidate chain templates.** For DMD, four
canonical mechanism templates (H01–H04) are well-established:
out-of-frame truncation with sarcolemmal fragility (H01);
in-frame rescue with partial function (H02); nonsense/splice-driven
NMD with tissue-graded loss (H03); distal-promoter variants with
isoform-specific loss (H04). Each template defines a specific walk
through the layer graph — which nodes it emphasizes and which
transitions it invokes. For less well-characterized diseases, an
LLM-driven generator would propose additional candidate walks.

**3. Attribute premises to chain positions.** Each premise in the
patient's bundle is routed to the node or edge it informs, for each
candidate hypothesis. The routing is source-driven: a Zhang record
always informs the variant and phenotype nodes; an ESM3 fold always
informs the protein node and the variant → protein edge; a lab
informs its own layer of the hierarchy (CK → cell type; LVEF →
tissue; 6MWT → phenotype). Each attribution is stored as a signed
weight — positive if the premise supports the hypothesis, negative
if it argues against.

**4. Compute the score vector for the assembled chain.** The
attribution counts + weights per layer and per edge roll up into
four dimensions:

- **Aggregate** — total evidence weight summed across nodes and edges
- **Coverage** — fraction of the seven biological layers that carry
  any evidence for this hypothesis
- **Consistency** — degree to which the premises agree with one
  another (v0 always 1.0; a later pass will detect contradictions)
- **Parsimony** — inverse of the number of empty layers, penalizing
  chains with unexplained gaps

The score is a vector, not a scalar. A hypothesis with a high
aggregate but low coverage is a chain that is well-supported at a few
layers but breaks in the middle. A hypothesis with modest aggregate
but 100% coverage is a chain that touches every layer of biology
even if none of the individual pieces are dramatic. The clinician
or the AAV design agent decides which dimension to prioritize.

## Illustration: Patient 5, hypothesis H01

Patient 5 carries a DMD frameshift variant at exon 75
(`p.(Ser3552Lysfs*6)`) and presents as a non-ambulant sixteen-year-old
with the classical severe motor phenotype. The workbench forms four
candidate hypotheses; H01 (out-of-frame → truncated dystrophin →
sarcolemmal fragility) is the top-ranked.

The assembly trace for H01 for this patient looks like:

| Layer | Evidence gathered | Source |
|---|---|---|
| **variant** | frameshift at exon 75, ACMG likely-pathogenic | Zhang 2024 record |
| variant → protein | patient's variant hits Dp427m/c/p through Dp71 in the isoform architecture; NMD prediction at cohort scale | isoform architecture, ClinVar NMD cross-tab |
| **protein** | truncated Dp71 folded via ESM3 (493 aa); WW domain preserved, C-terminal syntrophin-binding lost | ESM3 fold, isoform architecture |
| pathway → subcellular | dystrophin-glycoprotein complex assembly pathway | Reactome |
| **pathway** | five DMD-relevant Reactome pathways with specificity scores | Reactome |
| **subcellular** | *(no direct premise — a gap)* | — |
| cellType → tissue | HPA cell-type expression identifies which cells rely on the affected isoforms | Human Protein Atlas |
| **cell type** | myofiber damage biomarkers: CK 79× ULN; aldolase, LDH elevated | synthetic labs, HPA |
| **tissue** | LVEF 42% (cardiac impairment); FVC 17% predicted (respiratory failure); MRI fat fraction 91% (skeletal muscle replaced) | synthetic labs |
| **phenotype** | 6MWT 0m (non-ambulant); NSAA 6/34; time-to-stand 30s; patient phenotype label DMD | synthetic labs, Zhang record |

Coverage: six of seven layers evidenced (subcellular alone is empty
because we do not yet have a source that directly probes dystrophin
localization at the sarcolemma). Aggregate weight: 8.4. Parsimony: 0.5
(one gap penalizes it modestly).

The chain reads coherently top to bottom: the variant produces a
truncated protein that lacks the C-terminal DGC anchor, the DGC fails
to assemble, myofibers become mechanically fragile, muscle is
progressively replaced by fibro-fatty tissue, and the patient loses
ambulation. Each step has empirical support from at least one source.
The single gap (subcellular) is honestly flagged rather than papered
over with narrative.

## Compare to a hypothesis that scores lower

For the same Patient 5, H02 (in-frame rescue → BMD-like partial
function) is the lowest-ranked hypothesis (rank 4, aggregate 4.0,
coverage 3/7). Its chain is broken in the middle:

- Variant node: has evidence (Zhang record), but the variant is a
  frameshift — physically inconsistent with an in-frame rescue claim
- Protein node: some evidence (isoform impact)
- Pathway, subcellular, cell type, tissue nodes: no evidence
- Phenotype node: strong evidence, but the patient's clinical
  presentation (non-ambulant, severe motor decline) actively argues
  *against* a BMD-like partial-function mechanism

The chain has holes precisely where the mechanism story would need to
be strongest — the middle layers where the "partial function" claim
would need to be established. The score vector surfaces this: high
aggregate at a few nodes, but coverage is only 43%. A single scalar
score would have obscured the reason H02 doesn't fit.

## The pipeline

Assembly is a linear pipeline, run once per patient during hydration:

```
   raw sources        →   typed premises   →   chain-attributed hypotheses
   ─────────────         ──────────────       ────────────────────────────
   Zhang cohort           patient premises     H01 with chain links + score
   ClinVar variants       cohort premises      H02 with chain links + score
   HPA expression                              H03 with chain links + score
   Reactome pathways                           H04 with chain links + score
   Synthetic labs
   ESM3 folds
```

Each stage is idempotent and re-runnable. Adding a new source (an
AbSplice bake, an ESM2 LLR bake, a STRING PPI ingest) means:

1. Register the source in the `premise_source` registry with an ID
   and a version.
2. Emit its rows into the `premise` table.
3. Declare which layer or edge it informs via the source-to-position
   mapping.
4. Re-run the hypothesis bake — the world model consumes the enriched
   premise set automatically. Rankings shift where the new evidence
   fires; hypotheses re-form with richer chains.

No scorer code is rewritten. The scientific pipeline scales by adding
sources, not by adding rules.

## What ranking gives us that a scalar score does not

A vector-valued ranking supports different clinical and research
questions with the same underlying data:

- **Which hypothesis is most complete?** Sort by coverage. The chain
  with the fewest gaps is the one the model has evidence about at every
  layer of biology.
- **Which hypothesis has the strongest total evidence?** Sort by
  aggregate. Useful when a well-established but narrower mechanism is
  the right answer.
- **Which chain is most parsimonious?** Sort by parsimony. Prefer
  chains with fewer weak links over chains with many.
- **Where are the evidence gaps I should investigate?** For any
  hypothesis, layers with zero premises identify targeted experiments
  or data acquisitions that would strengthen the chain.
- **Which hypotheses share the same evidence at specific steps?**
  Cross-hypothesis analysis of node premises. If H01 and H03 both
  cite the same PTC premise at the variant → protein edge, that
  evidence is *overdetermined* — supporting either hypothesis with
  equal force.

## Known gaps and roadmap

**Subcellular-layer evidence.** No source currently informs the
subcellular localization node directly. Immunohistochemistry of
dystrophin at the sarcolemma would fill this — requires muscle
biopsy, uncommon in current clinical practice.

**Aberrant splicing predictions.** No AbSplice premise producer yet.
Splice-site patients (P10 in the current roster) rely on the coarse
ClinVar NMD cross-tab for the variant → protein edge; a per-patient
splice-outcome prediction would sharpen the evidence.

**ESM2 log-likelihood ratios.** Structural fold (ESM3) is baked for
two patients; per-variant functional disruption scoring (ESM2 LLR) is
not yet integrated. Would add a variant-level severity premise
distinct from the ACMG classification.

**Downstream network sources.** Open Targets molecular interactors,
STRING PPI, ARCHS4 co-expression, Geneformer perturbation predictions,
and CMap drug signatures are all called out in the rare-disease
platform deck as priorities. Each would extend the chain into the
molecular-neighborhood space beyond the direct pathway node.

**Contradiction detection.** The consistency dimension of the score
vector is a placeholder (always 1.0 in v0). A real implementation
would flag pairs of premises within a hypothesis that argue in
opposite directions (e.g., a phenotype label supporting H02 but
motor labs supporting H01), lowering consistency accordingly.

**LLM-refined claims.** The mechanism-claim text for each hypothesis
is currently template-filled. An LLM refinement pass conditioned on
the assembled premise bundle would produce per-patient customized
claims, "considered but discarded" transparency, and testable
predictions — layered on top of the deterministic template-scored
backbone.

**Novel hypothesis discovery.** For DMD the four canonical templates
cover essentially all documented mechanisms. For less-well-studied
rare diseases, novel walks through the layer graph — mechanisms not
in any template catalog — would be proposed by an LLM operating on
the assembled premise bundle, held in a review queue for human
curation before promotion.

## Summary

The workbench treats hypothesis formation as **evidence assembly along
a biological hierarchy**. Sources contribute typed premises; each
premise is attributed to the specific node or edge of the chain it
informs; each candidate hypothesis is a walk through the graph that
collects those attributions into a chain-decomposed evidence bundle;
ranking is a separate, vector-valued operation that compares whole
chains along multiple dimensions.

The result is a hypothesis object that carries not just a score but a
full audit trail: which sources argued for it, at which layer, with
what weight, and where the chain has gaps. That trace is what
distinguishes a real world-model output from an opaque ranking — and
it is what makes the pipeline extensible to new evidence sources
without rewriting the scoring logic.
