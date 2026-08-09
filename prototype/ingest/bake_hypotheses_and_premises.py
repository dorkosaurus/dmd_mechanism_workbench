"""Populate the premise registry + patient-scoped hypotheses/therapeutics.

Phase 1a of the world-model refactor:
- Register each data source / model as a `premise_source`
- Emit `premise` rows from existing baked substrate (labs, isoforms, HPA,
  Reactome, ClinVar NMD, ESM3 folds, Zhang cohort)
- Generate ≥3 (up to 4) `patient_hypothesis` rows per patient, each with
  a stored rationale + a `hypothesis_premise` audit trail
- Generate `patient_therapeutic` rows per top-ranked hypothesis
  (currently curated from the existing therapy nodes; will be replaced
  by a real AAV design model in later phases)

Migration-safe: creates tables if missing, idempotent (wipes prior rows
with matching generator_id before insert).

Run:
    python3 -m prototype.ingest.bake_hypotheses_and_premises
"""
from __future__ import annotations

import hashlib
import html
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DB = REPO / "data" / "mechanism.sqlite"

HYP_GENERATOR = ("HYP-MODEL", "v0-scored-2026-08")
AAV_GENERATOR = ("AAV-MODEL", "v0-curated-2026-08")
NOW = datetime.now(timezone.utc).isoformat(timespec="seconds")

# Migration DDL — safe to re-run on an already-built DB.
SCHEMA_MIGRATION = """
CREATE TABLE IF NOT EXISTS premise_source (
  source_id     TEXT PRIMARY KEY,
  source_type   TEXT NOT NULL,
  version       TEXT,
  description   TEXT,
  reference_url TEXT
);

CREATE TABLE IF NOT EXISTS premise (
  premise_id    TEXT PRIMARY KEY,
  source_id     TEXT NOT NULL,
  scope         TEXT NOT NULL,
  scope_key     TEXT NOT NULL,
  evidence      TEXT NOT NULL,
  confidence    REAL,
  provenance    TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS ix_premise_scope  ON premise(scope, scope_key);
CREATE INDEX IF NOT EXISTS ix_premise_source ON premise(source_id);

CREATE TABLE IF NOT EXISTS patient_hypothesis (
  hypothesis_id      TEXT PRIMARY KEY,
  patient_id         TEXT NOT NULL,
  variant_key        TEXT NOT NULL,
  mechanism_template TEXT,
  rank               INTEGER NOT NULL,
  score              REAL NOT NULL,
  confidence         REAL NOT NULL,
  claim              TEXT NOT NULL,
  rationale          TEXT NOT NULL,
  score_vector       TEXT,
  generator_id       TEXT NOT NULL,
  generator_version  TEXT NOT NULL,
  generated_at       TEXT NOT NULL,
  input_context_hash TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS ix_ph_patient ON patient_hypothesis(patient_id, rank);
CREATE INDEX IF NOT EXISTS ix_ph_variant ON patient_hypothesis(variant_key);

CREATE TABLE IF NOT EXISTS hypothesis_premise (
  hypothesis_id TEXT NOT NULL,
  premise_id    TEXT NOT NULL,
  weight        REAL NOT NULL,
  rationale     TEXT,
  PRIMARY KEY (hypothesis_id, premise_id)
);

CREATE TABLE IF NOT EXISTS hypothesis_chain_link (
  hypothesis_id TEXT NOT NULL,
  link_type     TEXT NOT NULL,       -- 'node' | 'edge'
  layer_from    TEXT NOT NULL,
  layer_to      TEXT NOT NULL,
  premise_id    TEXT NOT NULL,
  weight        REAL NOT NULL,
  rationale     TEXT,
  PRIMARY KEY (hypothesis_id, link_type, layer_from, layer_to, premise_id)
);
CREATE INDEX IF NOT EXISTS ix_hcl_hyp   ON hypothesis_chain_link(hypothesis_id);
CREATE INDEX IF NOT EXISTS ix_hcl_layer ON hypothesis_chain_link(layer_from, layer_to);

CREATE TABLE IF NOT EXISTS patient_therapeutic (
  therapeutic_id     TEXT PRIMARY KEY,
  patient_id         TEXT NOT NULL,
  hypothesis_id      TEXT NOT NULL,
  rank               INTEGER NOT NULL,
  score              REAL NOT NULL,
  confidence         REAL NOT NULL,
  modality           TEXT NOT NULL,
  design             TEXT NOT NULL,
  rationale          TEXT NOT NULL,
  eligibility_status TEXT,
  generator_id       TEXT NOT NULL,
  generator_version  TEXT NOT NULL,
  generated_at       TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS ix_pt_patient ON patient_therapeutic(patient_id, rank);
CREATE INDEX IF NOT EXISTS ix_pt_hyp     ON patient_therapeutic(hypothesis_id);
"""

# The 10-patient roster (mirrors hydrate_patient_view.ROSTER).
ROSTER_KEYS = [
    ("S1_novel",    "2"),
    ("S1_novel",    "30"),
    ("S2_reported", "258"),
    ("S1_novel",    "57"),
    ("S1_novel",    "5"),
    ("S1_novel",    "11"),
    ("S2_reported", "202"),
    ("S2_reported", "225"),
    ("S1_novel",    "49"),
    ("S2_reported", "266"),
]

# The four canonical mechanism templates (kept in mechanism.sqlite.hypotheses).
HYP_TEMPLATES = [
    ("01", "Out-of-frame deletions → truncated dystrophin → sarcolemmal fragility"),
    ("02", "In-frame deletions → partial-function dystrophin → BMD phenotype"),
    ("03", "Nonsense / splice variants → NMD → tissue-graded transcript loss"),
    ("04", "Distal-promoter variants → tissue-specific isoform loss (Dp140/Dp71)"),
]


# ----------------------------------------------------------------------
# Premise source registry
# ----------------------------------------------------------------------
PREMISE_SOURCES = [
    ("zhang_2024",     "data",  "supp S1+S2 (2024)", "Per-patient genotype/phenotype records (Zhang et al. 2024, PMC11344408)", "https://pmc.ncbi.nlm.nih.gov/articles/PMC11344408/"),
    ("clinvar_nmd",    "model", "aenmd-rule v1",     "NMD prediction (PTC ≥50nt upstream of last EEJ) applied to ClinVar DMD variants", None),
    ("hpa_expression", "data",  "HPA 2024",          "Human Protein Atlas single cell type specific nCPM", "https://www.proteinatlas.org/"),
    ("reactome",       "data",  "Reactome v96",      "Reactome DMD pathway memberships + specificity scoring", "https://reactome.org/"),
    ("isoform_arch",   "data",  "UniProt+RefSeq",    "DMD isoform architecture (first-shared-exon per isoform)", "https://www.uniprot.org/uniprotkb/P11532"),
    ("synthetic_labs", "model", "synthetic_v1",      "Per-patient synthetic clinical labs (Birnkrant 2018 ranges)", None),
    ("esm3",           "model", "esm3-open-2024-03", "ESM3 protein fold + InterPro function annotation", "https://forge.evolutionaryscale.ai/"),
]


def register_sources(conn):
    conn.execute("DELETE FROM premise_source")
    for row in PREMISE_SOURCES:
        conn.execute("INSERT INTO premise_source VALUES (?,?,?,?,?)", row)


# ----------------------------------------------------------------------
# Premise emission helpers
# ----------------------------------------------------------------------
def _pid(*parts: str) -> str:
    """Stable ID from parts."""
    return "|".join(str(p) for p in parts)


def _hash(*parts: str) -> str:
    return hashlib.sha1("|".join(str(p) for p in parts).encode()).hexdigest()[:12]


def _prov(cache_path: str | None = None) -> str:
    return json.dumps({"fetched_at": NOW, "cache": cache_path})


def emit_zhang_premise(conn, cohort: str, pid: str, phen: str, exon: str,
                        consequence: str, aa_change: str | None, acmg: str | None) -> str:
    premise_id = _pid("zhang", cohort, pid)
    ev = {
        "cohort": cohort, "patient_id": pid, "phenotype": phen,
        "exon": exon, "consequence": consequence, "aa_change": aa_change, "acmg": acmg,
    }
    conn.execute(
        "INSERT OR REPLACE INTO premise VALUES (?,?,?,?,?,?,?)",
        (premise_id, "zhang_2024", "patient", pid, json.dumps(ev), 1.0, _prov()),
    )
    return premise_id


def emit_lab_premises(conn, cohort: str, pid: str) -> list[str]:
    """Emit one premise per abnormal lab for this patient."""
    ids = []
    rows = conn.execute(
        "SELECT assay_key, label, layer, tissue, value, ref_low, ref_high, flag, unit "
        "FROM patient_labs WHERE cohort=? AND patient_id=? AND flag != 'normal'",
        (cohort, pid),
    ).fetchall()
    for (key, label, layer, tissue, value, ref_lo, ref_hi, flag, unit) in rows:
        premise_id = _pid("labs", cohort, pid, key)
        fold = (value / ref_hi) if (flag == "high" and ref_hi > 0) else \
               (ref_lo / value) if (flag == "low" and value > 0) else 1.0
        ev = {"assay": key, "label": label, "layer": layer, "tissue": tissue,
              "value": value, "unit": unit, "flag": flag,
              "ref_low": ref_lo, "ref_high": ref_hi, "fold": round(fold, 2)}
        conn.execute(
            "INSERT OR REPLACE INTO premise VALUES (?,?,?,?,?,?,?)",
            (premise_id, "synthetic_labs", "patient", pid, json.dumps(ev), 0.85, _prov()),
        )
        ids.append(premise_id)
    return ids


def emit_isoform_premise(conn, cohort: str, pid: str, exon_n: int | None) -> str:
    """Per-patient isoform-impact projection (which isoforms are hit)."""
    if exon_n is None: return ""
    hit = []; spared = []
    for iso in conn.execute(
        "SELECT isoform_id, first_shared_exon FROM isoforms ORDER BY rank"
    ):
        (hit if iso[1] <= exon_n else spared).append(iso[0])
    ev = {"variant_exon": exon_n, "hit_isoforms": hit, "spared_isoforms": spared}
    premise_id = _pid("iso_arch", cohort, pid)
    conn.execute(
        "INSERT OR REPLACE INTO premise VALUES (?,?,?,?,?,?,?)",
        (premise_id, "isoform_arch", "patient", pid, json.dumps(ev), 0.9, _prov()),
    )
    return premise_id


def emit_hpa_premise(conn) -> str:
    """One cohort-scope premise: DMD expression by HPA cell type."""
    rows = conn.execute(
        "SELECT cell_type, score, tissue FROM celltype_expression "
        "WHERE gene_symbol='DMD' ORDER BY score DESC LIMIT 7"
    ).fetchall()
    ev = {"cell_types": [{"name": r[0], "score": r[1], "tissue": r[2]} for r in rows]}
    premise_id = _pid("hpa", "DMD")
    conn.execute(
        "INSERT OR REPLACE INTO premise VALUES (?,?,?,?,?,?,?)",
        (premise_id, "hpa_expression", "cohort", "DMD", json.dumps(ev), 0.95, _prov()),
    )
    return premise_id


def emit_reactome_premise(conn) -> str:
    rows = conn.execute(
        "SELECT pathway_name, score FROM pathway_enrichment "
        "WHERE gene_symbol='DMD' ORDER BY score DESC LIMIT 5"
    ).fetchall()
    ev = {"pathways": [{"name": r[0], "specificity": r[1]} for r in rows]}
    premise_id = _pid("reactome", "DMD")
    conn.execute(
        "INSERT OR REPLACE INTO premise VALUES (?,?,?,?,?,?,?)",
        (premise_id, "reactome", "cohort", "DMD", json.dumps(ev), 0.9, _prov()),
    )
    return premise_id


def emit_nmd_cohort_premise(conn) -> str:
    """Cohort-scope NMD × ACMG cross-tab (built at hydration time; here
    we materialize it as a premise so hypothesis rationales can cite it)."""
    # Grid from the hydrator's build_nmd_cohort — replicate the summary here
    # for storage. We just re-count from clinvar_phenotype.
    from prototype.ingest.hydrate_mechanism_view import build_nmd_cohort
    ev = build_nmd_cohort(conn)
    premise_id = _pid("nmd_cohort", "DMD")
    conn.execute(
        "INSERT OR REPLACE INTO premise VALUES (?,?,?,?,?,?,?)",
        (premise_id, "clinvar_nmd", "cohort", "DMD", json.dumps(ev), 0.9, _prov()),
    )
    return premise_id


def emit_esm3_premise(conn, cohort: str, pid: str) -> str | None:
    """ESM3 fold premise for patients whose variant we've folded."""
    known = {
        ("S2_reported", "258"): {"isoform": "Dp71", "cache": "P11532_Dp71",
                                  "variant_resi_full": 3083, "variant_label": "Trp3083*"},
        ("S1_novel",    "5"):   {"isoform": "Dp71", "cache": "P11532_Dp71",
                                  "variant_resi_full": 3552, "variant_label": "Ser3552Lysfs*6",
                                  "truncated_cache": "P11532_Dp71_truncP5"},
    }
    k = known.get((cohort, pid))
    if not k: return None
    premise_id = _pid("esm3", cohort, pid)
    ev = dict(k)
    ev["mean_plddt"] = 0.69
    conn.execute(
        "INSERT OR REPLACE INTO premise VALUES (?,?,?,?,?,?,?)",
        (premise_id, "esm3", "patient", pid, json.dumps(ev), 0.7, _prov(k.get("cache"))),
    )
    return premise_id


# ----------------------------------------------------------------------
# Hypothesis generation (world-model stub v0)
# ----------------------------------------------------------------------
# For now this replicates the scoring rules from hydrate_patient_view.py,
# but as a "generator" that emits stored hypothesis rows with rationale
# + a hypothesis_premise audit trail.

def parse_exon(s: str | None) -> int | None:
    if not s: return None
    t = s.strip()
    if t.lower().startswith("int"): t = t[3:]
    try: return int(t)
    except ValueError: return None


# (weight, rationale) per (template, feature) — mirrors the current
# hydrate_patient_view.score_hypotheses rules.
def score_templates(consequence: str, phenotype: str, exon_n: int | None
                    ) -> dict[str, tuple[float, str]]:
    cons = (consequence or "").lower()
    phen = phenotype or ""
    out: dict[str, tuple[float, str]] = {}

    # H01 — out-of-frame → truncation
    if cons == "frameshift":
        out["01"] = (9.5, "frameshift → out-of-frame → truncated dystrophin (Monaco rule)")
    elif cons == "nonsense":
        out["01"] = (6.0, "premature stop → truncated protein (before NMD)")
    elif cons == "splice-site":
        out["01"] = (5.0, "splice defect can produce out-of-frame transcript")
    elif cons == "missense":
        out["01"] = (1.5, "missense does not truncate — poor fit")
    else:
        out["01"] = (0.5, "not a truncating variant")

    # H02 — in-frame / partial-function / BMD
    if phen == "BMD":
        out["02"] = (9.0, "BMD phenotype label — partial-function dystrophin consistent")
    elif phen == "IMD":
        out["02"] = (6.5, "intermediate phenotype — some partial rescue plausible")
    elif cons == "missense":
        out["02"] = (5.5, "missense can preserve some dystrophin function")
    else:
        out["02"] = (1.5, "no evidence of in-frame rescue")

    # H03 — NMD / PTC
    if cons == "nonsense":
        out["03"] = (9.5, "PTC → NMD-mediated transcript loss (direct fit)")
    elif cons == "splice-site":
        out["03"] = (9.0, "splice defect → downstream PTC → NMD-eligible")
    elif cons == "frameshift":
        out["03"] = (8.5, "frameshift creates downstream PTC → NMD-eligible")
    else:
        out["03"] = (1.5, "no PTC generated")

    # H04 — distal-isoform loss
    if exon_n is None:
        out["04"] = (1.0, "variant position not parseable")
    elif exon_n >= 63:
        out["04"] = (5.0, f"exon {exon_n} — hits Dp71 (last ubiquitous isoform)")
    elif exon_n >= 56:
        out["04"] = (4.0, f"exon {exon_n} — hits Dp116 (Schwann)")
    elif exon_n >= 45:
        out["04"] = (3.0, f"exon {exon_n} — hits Dp140 (brain/kidney)")
    else:
        out["04"] = (1.5, f"exon {exon_n} — proximal, doesn't selectively affect distal isoforms")
    return out


# Rough claim templates per hypothesis, filled in per patient.
CLAIM_TEMPLATES = {
    "01": "Patient {p}'s {cons} at exon {exon} produces an out-of-frame transcript → truncated dystrophin lacking the C-terminal DGC anchor → sarcolemmal fragility in skeletal + cardiac muscle.",
    "02": "Patient {p}'s {cons} at exon {exon} may retain partial function via in-frame rescue → BMD-like presentation with slower progression.",
    "03": "Patient {p}'s {cons} at exon {exon} generates a PTC upstream of the last exon-exon junction → NMD degrades the transcript → near-total Dp427m loss.",
    "04": "Patient {p}'s {cons} at exon {exon} selectively ablates distal isoforms (Dp140 / Dp116 / Dp71) → tissue-specific dysfunction beyond muscle.",
}


# Which premise sources support which templates (for the audit trail).
# Weight interpretation: score contribution / 10 (relative unit).
# Chain layers (must match the order used in the layered biological view).
LAYERS = ["variant", "protein", "pathway", "subcellular",
          "cellType", "tissue", "phenotype"]


# Attribute a premise (by its source_id) to a position in the biological
# chain. A premise informs either a NODE (layer_from == layer_to) or an
# EDGE (adjacent layers). Some sources inform multiple positions.
def premise_chain_positions(source_id: str, evidence: dict) -> list[tuple[str, str, str]]:
    """Return [(link_type, layer_from, layer_to)] positions for this premise."""
    if source_id == "zhang_2024":
        # Variant record itself + patient phenotype label
        return [("node", "variant", "variant"), ("node", "phenotype", "phenotype")]
    if source_id == "clinvar_nmd":
        # Cohort NMD cross-tab argues for the variant→protein transition
        return [("edge", "variant", "protein")]
    if source_id == "isoform_arch":
        # Per-patient isoform-hit projection: variant→protein edge + protein node
        return [("edge", "variant", "protein"), ("node", "protein", "protein")]
    if source_id == "esm3":
        # ESM3 fold: protein node + evidence for the variant→protein transition
        return [("node", "protein", "protein"), ("edge", "variant", "protein")]
    if source_id == "reactome":
        # Pathway memberships: pathway node + edge to cellType
        return [("node", "pathway", "pathway"), ("edge", "pathway", "subcellular")]
    if source_id == "hpa_expression":
        # HPA cell-type expression: cellType node + edge to tissue
        return [("node", "cellType", "cellType"), ("edge", "cellType", "tissue")]
    if source_id == "synthetic_labs":
        # Each lab has its own layer (evidence.layer field: cellType|tissueType|phenotype)
        layer = evidence.get("layer", "phenotype")
        # Our labs use "tissueType" in the labs table but LAYERS uses "tissue"
        norm = "tissue" if layer == "tissueType" else layer
        return [("node", norm, norm)]
    return []


def premise_weights_for_template(hyp_id: str, patient_premises: dict, cohort_premises: dict
                                  ) -> list[tuple[str, float, str]]:
    """Return [(premise_id, signed_weight, rationale)] linking premises to this hyp."""
    out = []
    zhang_pid = patient_premises.get("zhang")
    iso_pid   = patient_premises.get("iso")
    esm3_pid  = patient_premises.get("esm3")
    lab_pids  = patient_premises.get("labs", [])
    nmd_pid   = cohort_premises.get("nmd_cohort")
    hpa_pid   = cohort_premises.get("hpa")
    reactome_pid = cohort_premises.get("reactome")

    # Every hypothesis draws on the Zhang variant record.
    if zhang_pid:
        out.append((zhang_pid, 1.0, "variant + phenotype + exon from Zhang 2024 record"))

    # Isoform architecture projection is the substrate for isoform-impact reasoning.
    if iso_pid:
        out.append((iso_pid, 0.7, "per-isoform hit/spared projection from variant exon position"))

    # NMD cohort premise informs H03 primarily.
    if nmd_pid and hyp_id == "03":
        out.append((nmd_pid, 0.8, "ClinVar 11,790-variant NMD × ACMG cross-tab supports NMD mechanism at cohort scale"))

    # HPA expression informs H01/H03/H04 (which cell types express which isoforms).
    if hpa_pid and hyp_id in ("01", "03", "04"):
        out.append((hpa_pid, 0.5, "HPA cell-type expression identifies which tissues lose dystrophin"))

    # Reactome informs H01 (DGC pathway).
    if reactome_pid and hyp_id == "01":
        out.append((reactome_pid, 0.4, "Reactome DGC-formation pathway is a first-order consequence of dystrophin loss"))

    # ESM3 fold informs H01/H03 by making the truncated protein claim structural.
    if esm3_pid and hyp_id in ("01", "03"):
        out.append((esm3_pid, 0.6, "ESM3 fold of WT + truncated dystrophin visualizes the DGC-anchor loss"))

    # Abnormal labs contribute per-lab weight for the mechanism they support.
    if hyp_id == "01":
        # elevated CK, low LVEF, low FVC, high MRI ff → all support muscle fragility
        for lp in lab_pids:
            if any(k in lp for k in ("CK", "LVEF", "FVC", "MRI_ff", "m6MWT")):
                out.append((lp, 0.4, "abnormal lab consistent with H01 sarcolemmal-fragility mechanism"))
    elif hyp_id == "02":
        for lp in lab_pids:
            if any(k in lp for k in ("NSAA", "m6MWT")):
                out.append((lp, 0.3, "motor function retention weakly supports partial-rescue hypothesis"))
    elif hyp_id == "03":
        for lp in lab_pids:
            if any(k in lp for k in ("CK", "MRI_ff")):
                out.append((lp, 0.5, "muscle-damage lab consistent with total protein loss from NMD"))
            if "IQ" in lp:
                out.append((lp, 0.2, "tissue-graded NMD escape — IQ status is diagnostic"))
    elif hyp_id == "04":
        for lp in lab_pids:
            if any(k in lp for k in ("IQ", "ERG", "UACR")):
                out.append((lp, 0.6, "distal-tissue lab abnormality points at distal isoform loss"))

    return out


def bake_hypotheses_for_patient(conn, cohort: str, pid: str, phen: str, age: float | None,
                                 amb: str | None, exon_str: str | None, nuc: str | None,
                                 aa: str | None, cons: str | None, acmg: str | None
                                 ) -> tuple[int, int]:
    """Emit premises + patient_hypothesis + patient_therapeutic rows for one patient.
    Returns (n_hyps, n_therapies)."""
    exon_n = parse_exon(exon_str)
    variant_key = f"{cohort}#{pid}:{nuc or '?'}"

    # 1. Emit patient-scope premises
    zhang_pid = emit_zhang_premise(conn, cohort, pid, phen, exon_str, cons, aa, acmg)
    iso_pid   = emit_isoform_premise(conn, cohort, pid, exon_n)
    esm3_pid  = emit_esm3_premise(conn, cohort, pid)  # None for most
    lab_pids  = emit_lab_premises(conn, cohort, pid)

    patient_premises = {"zhang": zhang_pid, "iso": iso_pid, "esm3": esm3_pid, "labs": lab_pids}
    cohort_premises = {
        "nmd_cohort": _pid("nmd_cohort", "DMD"),
        "hpa":        _pid("hpa", "DMD"),
        "reactome":   _pid("reactome", "DMD"),
    }

    # 2. Score templates and emit patient_hypothesis rows
    scores = score_templates(cons, phen, exon_n)
    ranked = sorted(scores.items(), key=lambda kv: -kv[1][0])
    input_hash = _hash(variant_key, phen, cons or "", str(exon_n), acmg or "")

    # Look up each premise's source_id + evidence so we can attribute
    # it to a chain position. Cache the lookup for this patient.
    def _get_premise_source_and_ev(premise_id: str) -> tuple[str, dict]:
        row = conn.execute(
            "SELECT source_id, evidence FROM premise WHERE premise_id=?",
            (premise_id,),
        ).fetchone()
        if not row: return ("unknown", {})
        try: ev = json.loads(row[1])
        except Exception: ev = {}
        return (row[0], ev)

    n_hyps = 0
    for rank, (tmpl_id, (score, fit)) in enumerate(ranked, start=1):
        hyp_id = f"P_{cohort}#{pid}:h{tmpl_id}:v1"
        claim = CLAIM_TEMPLATES[tmpl_id].format(
            p=f"{cohort}#{pid}",
            cons=(cons or "variant").lower(),
            exon=exon_str or "?",
        )
        rationale = f"Rule fired: {fit}. Score {score:.2f}/10 (rank #{rank})."

        # Emit hypothesis_premise + hypothesis_chain_link audit rows
        conn.execute("DELETE FROM hypothesis_premise    WHERE hypothesis_id=?", (hyp_id,))
        conn.execute("DELETE FROM hypothesis_chain_link WHERE hypothesis_id=?", (hyp_id,))
        prem_weights = premise_weights_for_template(tmpl_id, patient_premises, cohort_premises)

        # Track per-layer contributions to compute the score vector.
        layer_evidence: dict[str, float] = {layer: 0.0 for layer in LAYERS}
        edge_evidence:  dict[tuple[str, str], float] = {}
        n_chain_links = 0
        for (pr_id, w, rat) in prem_weights:
            if not pr_id: continue
            conn.execute(
                "INSERT OR REPLACE INTO hypothesis_premise VALUES (?,?,?,?)",
                (hyp_id, pr_id, w, rat),
            )
            source_id, ev = _get_premise_source_and_ev(pr_id)
            for (link_type, lf, lt) in premise_chain_positions(source_id, ev):
                conn.execute(
                    "INSERT OR REPLACE INTO hypothesis_chain_link VALUES (?,?,?,?,?,?,?)",
                    (hyp_id, link_type, lf, lt, pr_id, w, rat),
                )
                n_chain_links += 1
                if link_type == "node":
                    layer_evidence[lf] = layer_evidence.get(lf, 0) + abs(w)
                else:
                    edge_evidence[(lf, lt)] = edge_evidence.get((lf, lt), 0) + abs(w)

        # Compute score vector:
        #   aggregate   = sum of all evidence weights (nodes + edges)
        #   coverage    = fraction of layers with any evidence
        #   consistency = 1.0 (no contradiction detection yet — placeholder)
        #   parsimony   = 1 / (1 + number of empty layers) — chains with
        #                 fewer gaps score higher
        n_layers_covered = sum(1 for v in layer_evidence.values() if v > 0)
        n_layers_empty   = len(LAYERS) - n_layers_covered
        aggregate = sum(layer_evidence.values()) + sum(edge_evidence.values())
        coverage  = n_layers_covered / len(LAYERS)
        consistency = 1.0
        parsimony = 1.0 / (1.0 + n_layers_empty)
        score_vector = {
            "aggregate":   round(aggregate, 3),
            "coverage":    round(coverage, 3),
            "consistency": round(consistency, 3),
            "parsimony":   round(parsimony, 3),
            "layerScores": {layer: round(v, 3) for layer, v in layer_evidence.items()},
            "edgeScores":  {f"{lf}->{lt}": round(v, 3) for (lf, lt), v in edge_evidence.items()},
            "chainLinks":  n_chain_links,
        }

        conn.execute(
            "INSERT OR REPLACE INTO patient_hypothesis "
            "(hypothesis_id, patient_id, variant_key, mechanism_template, rank, score, "
            " confidence, claim, rationale, score_vector, generator_id, generator_version, "
            " generated_at, input_context_hash) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (hyp_id, pid, variant_key, tmpl_id, rank, score,
             min(1.0, score / 10.0), claim, rationale,
             json.dumps(score_vector),
             HYP_GENERATOR[0], HYP_GENERATOR[1], NOW, input_hash),
        )
        n_hyps += 1

    # 3. Emit patient_therapeutic rows for the top hypothesis
    top_hyp_id = f"P_{cohort}#{pid}:h{ranked[0][0]}:v1"
    n_therapies = emit_therapeutics(conn, cohort, pid, ranked[0][0], top_hyp_id)
    return n_hyps, n_therapies


# Curated therapy templates per mechanism template. Real AAV design model
# would generate these on the fly conditioned on patient immune profile;
# for v0 we use the therapy_node text from the hypotheses table plus
# structured design blobs.
THERAPY_DESIGNS = {
    "01": [
        {"modality": "AAV",         "design": {"capsid": "AAVrh74", "promoter": "MHCK7", "transgene": "µ-dystrophin"},
         "rationale": "Delivers a truncated but functional dystrophin under a muscle-specific promoter; rescues DGC assembly.",
         "eligibility": "screening_required"},
        {"modality": "ASO",         "design": {"target_exon": "51", "chemistry": "PMO"},
         "rationale": "Exon-skipping restores in-frame Dp427m for eligible deletions.",
         "eligibility": "screening_required"},
    ],
    "02": [
        {"modality": "small_molecule", "design": {"agent": "utrophin_upregulator"},
         "rationale": "Utrophin upregulation compensates for reduced dystrophin function in BMD-like patients.",
         "eligibility": "unknown"},
    ],
    "03": [
        {"modality": "readthrough",    "design": {"agent": "ataluren", "target": "PTC"},
         "rationale": "Promotes ribosomal readthrough of PTCs — restores partial protein from NMD-degraded transcripts.",
         "eligibility": "screening_required"},
    ],
    "04": [
        {"modality": "AAV",            "design": {"capsid": "AAV9", "promoter": "tissue_restricted", "transgene": "Dp140_or_Dp71_mini"},
         "rationale": "Isoform-selective AAV delivers a Dp140/Dp71 mini-gene under a CNS/renal promoter to rescue distal-isoform loss.",
         "eligibility": "unknown"},
    ],
}


def emit_therapeutics(conn, cohort: str, pid: str, tmpl_id: str, hyp_id: str) -> int:
    designs = THERAPY_DESIGNS.get(tmpl_id, [])
    for i, d in enumerate(designs, start=1):
        tid = f"P_{cohort}#{pid}:h{tmpl_id}:t{i}"
        conn.execute(
            "INSERT OR REPLACE INTO patient_therapeutic VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (tid, pid, hyp_id, i, 0.7, 0.7, d["modality"],
             json.dumps(d["design"]), d["rationale"], d["eligibility"],
             AAV_GENERATOR[0], AAV_GENERATOR[1], NOW),
        )
    return len(designs)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def _ensure_column(conn, table: str, column: str, decl: str) -> None:
    cols = {r[1] for r in conn.execute(f"PRAGMA table_info({table})")}
    if column not in cols:
        conn.execute(f"ALTER TABLE {table} ADD COLUMN {column} {decl}")


def main() -> None:
    conn = sqlite3.connect(DB)
    conn.executescript(SCHEMA_MIGRATION)
    _ensure_column(conn, "patient_hypothesis", "score_vector", "TEXT")

    # 1. Register premise sources
    register_sources(conn)
    print(f"[sources]  {len(PREMISE_SOURCES)} premise sources registered")

    # 2. Emit cohort-scope premises (shared across patients)
    emit_hpa_premise(conn)
    emit_reactome_premise(conn)
    emit_nmd_cohort_premise(conn)
    print("[cohort]   HPA + Reactome + NMD cohort premises emitted")

    # 3. Wipe & regenerate per-patient hypotheses + therapeutics for the roster
    #    (we scope by generator_id so re-runs are idempotent across versions).
    conn.execute("DELETE FROM patient_hypothesis WHERE generator_id=? AND generator_version=?",
                 HYP_GENERATOR)
    conn.execute("DELETE FROM patient_therapeutic WHERE generator_id=? AND generator_version=?",
                 AAV_GENERATOR)

    total_h, total_t = 0, 0
    for (cohort, pid) in ROSTER_KEYS:
        row = conn.execute(
            "SELECT phenotype_label, age_years, ambulatory, exon, "
            "       nucleotide, aa_change, consequence, acmg "
            "FROM patient_phenotype WHERE cohort=? AND patient_id=?",
            (cohort, pid),
        ).fetchone()
        if not row:
            print(f"[warn] {cohort}#{pid} not in patient_phenotype")
            continue
        (phen, age, amb, exon, nuc, aa, cons, acmg) = row
        aa = html.unescape(aa) if aa else aa
        nuc = html.unescape(nuc) if nuc else nuc
        n_h, n_t = bake_hypotheses_for_patient(
            conn, cohort, pid, phen, age, amb, exon, nuc, aa, cons, acmg)
        total_h += n_h
        total_t += n_t

    conn.commit()

    # 4. Sanity summary
    n_prem = conn.execute("SELECT COUNT(*) FROM premise").fetchone()[0]
    n_hp   = conn.execute("SELECT COUNT(*) FROM patient_hypothesis").fetchone()[0]
    n_hpp  = conn.execute("SELECT COUNT(*) FROM hypothesis_premise").fetchone()[0]
    n_pt   = conn.execute("SELECT COUNT(*) FROM patient_therapeutic").fetchone()[0]
    print(f"\n[premises] {n_prem:>5} rows (cohort + patient scope)")
    print(f"[patient_hypothesis]    {n_hp:>5} rows ({total_h} freshly baked)")
    print(f"[hypothesis_premise]    {n_hpp:>5} audit-trail rows")
    print(f"[patient_therapeutic]   {n_pt:>5} rows ({total_t} freshly baked)")
    print()

    # Spot-check
    print("Sample: top hypothesis + supporting premises for Patient 5 (S1_novel#5):")
    for row in conn.execute("""
        SELECT h.mechanism_template, h.rank, h.score,
               COUNT(hp.premise_id) as n_premises
        FROM patient_hypothesis h
        LEFT JOIN hypothesis_premise hp ON h.hypothesis_id = hp.hypothesis_id
        WHERE h.patient_id = '5'
        GROUP BY h.hypothesis_id ORDER BY h.rank
    """):
        print(f"  H{row[0]} rank #{row[1]} score {row[2]:.1f} · {row[3]} premises consulted")

    conn.close()


if __name__ == "__main__":
    main()
