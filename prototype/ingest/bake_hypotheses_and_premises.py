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
import re
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
    ("zhang_2024",           "data",  "supp S1+S2 (2024)",  "Per-patient genotype/phenotype records (Zhang et al. 2024, PMC11344408)", "https://pmc.ncbi.nlm.nih.gov/articles/PMC11344408/"),
    ("clinvar_nmd",          "model", "aenmd-rule v1",      "NMD prediction (PTC ≥50nt upstream of last EEJ) applied to ClinVar DMD variants", None),
    ("hpa_expression",       "data",  "HPA 2024",           "Human Protein Atlas single cell type specific nCPM", "https://www.proteinatlas.org/"),
    ("reactome",             "data",  "Reactome v96",       "Reactome DMD pathway memberships + specificity scoring", "https://reactome.org/"),
    ("isoform_arch",         "data",  "UniProt+RefSeq",     "DMD isoform architecture (first-shared-exon per isoform)", "https://www.uniprot.org/uniprotkb/P11532"),
    ("synthetic_labs",       "model", "synthetic_v1",       "Per-patient synthetic clinical labs (Birnkrant 2018 ranges)", None),
    ("esm3",                 "model", "esm3-open-2024-03",  "ESM3 protein fold + InterPro function annotation", "https://forge.evolutionaryscale.ai/"),
    # A: composition premises — close the variant→cellType and cellType→tissue chain gap
    ("patient_celltype_impact", "model", "compose_v1",      "Per-patient cell-type impact: isoform_arch × curated cell-to-isoform dependency map (Muntoni 2003, Pillers 1993, Lidov 1995)", None),
    ("patient_tissue_impact",   "model", "compose_v1",      "Per-patient tissue impact: composes cell-type impact × isoforms.primary_expression_tissues", None),
    # B: literature as first-class premise (migrated from hypothesis_chain_edge_evidence)
    ("literature",              "data",  "curated_2026-08", "Peer-reviewed citations anchoring specific chain-edge and chain-node claims (Monaco 1988, Popp & Maquat 2013, Ervasti & Campbell 1993, Pillers 1993, ...)", None),
    ("uniprot_subcellular",     "data",  "UniProt REST",    "UniProt-curated subcellular localization (sarcolemma, cytoskeleton, postsynaptic membrane for DMD P11532)", "https://www.uniprot.org/uniprotkb/P11532"),
    ("absplice",                "model", "absplice_v1.0.4", "AbSplice-DNA per-tissue aberrant-splicing probability (Wagner et al. Nat Genet 2023). Distinguishes canonical splice-site variants from exonic ones — flags cryptic-splice risk that isoform_arch and clinvar_nmd cannot see.", "https://github.com/gagneurlab/absplice"),
]

# Cell-type → required-isoform dependency map (mirrors hydrate_patient_view.CELL_TO_ISOFORMS).
CELL_TO_ISOFORMS: dict[str, list[str]] = {
    "Myonuclei":                    ["Dp427m"],
    "Cardiomyocytes":               ["Dp427m"],
    "Thymic myoid cells":           ["Dp427m"],
    "Salivary myoepithelial cells": ["Dp427m"],
    "Rod photoreceptor cells":      ["Dp260", "Dp71"],
    "Cone photoreceptor cells":     ["Dp260", "Dp71"],
    "Adipocytes":                   ["Dp71"],
}

# Node → biological-layer mapping per curated hypothesis chain
# (hypothesis_chain_nodes uses v1/v2/v3, m1/m2/m3, p1/p2/p3 node IDs;
# this table maps each to its position in the seven-layer hierarchy).
CHAIN_NODE_TO_LAYER: dict[str, dict[str, str]] = {
    '01': {
        'v1': 'variant', 'v2': 'variant', 'v3': 'protein',
        'm1': 'subcellular', 'm2': 'pathway', 'm3': 'subcellular',
        'p1': 'cellType', 'p2': 'cellType', 'p3': 'tissue',
    },
    '02': {
        'v1': 'variant', 'v2': 'protein', 'v3': 'protein',
        'm1': 'subcellular', 'm2': 'subcellular', 'm3': 'pathway',
        'p1': 'cellType', 'p2': 'phenotype', 'p3': 'phenotype',
    },
    '03': {
        'v1': 'variant', 'v2': 'variant', 'v3': 'protein',
        'm1': 'subcellular', 'm2': 'cellType', 'm3': 'cellType',
        'p1': 'cellType', 'p2': 'phenotype', 'p3': 'phenotype',
    },
    '04': {
        'v1': 'variant', 'v2': 'protein', 'v3': 'protein',
        'm1': 'cellType', 'm2': 'cellType', 'm3': 'cellType',
        'p1': 'phenotype', 'p2': 'phenotype', 'p3': 'phenotype',
    },
}


def _slug(s: str) -> str:
    """Sanitize a citation string into a stable slug."""
    return re.sub(r'[^A-Za-z0-9]+', '_', s).strip('_').lower()


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


def emit_uniprot_subcellular_premise(conn, uniprot_id: str = "P11532") -> str:
    """Fetch UniProt subcellular annotations and emit as a cohort-scope
    premise for the subcellular node. Result is cached under
    data/raw/uniprot_{id}_subcellular.json — subsequent runs skip the fetch.
    Fills the subcellular chain-node gap with curated protein-biology data."""
    cache = REPO / "data" / "raw" / f"uniprot_{uniprot_id}_subcellular.json"
    if cache.exists() and cache.stat().st_size > 100:
        d = json.loads(cache.read_text())
    else:
        import urllib.request
        cache.parent.mkdir(parents=True, exist_ok=True)
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=15) as resp:
            d = json.loads(resp.read())
        cache.write_text(json.dumps(d, indent=2))

    # Extract subcellular location comments
    locations = []
    notes = []
    for comment in d.get("comments", []):
        if comment.get("commentType") != "SUBCELLULAR LOCATION":
            continue
        for loc_entry in comment.get("subcellularLocations", []):
            loc = loc_entry.get("location", {}).get("value")
            topology = loc_entry.get("topology", {}).get("value", "")
            if loc:
                locations.append({"location": loc, "topology": topology})
        for note in comment.get("note", {}).get("texts", []):
            if note.get("value"):
                notes.append(note["value"])

    ev = {
        "uniprot_id": uniprot_id,
        "locations":  locations,
        "notes":      notes[:3],
        "gene":       "DMD",
    }
    premise_id = _pid("uniprot_subcellular", "DMD")
    conn.execute(
        "INSERT OR REPLACE INTO premise VALUES (?,?,?,?,?,?,?)",
        (premise_id, "uniprot_subcellular", "cohort", "DMD",
         json.dumps(ev), 0.95, _prov(str(cache.relative_to(REPO)))),
    )
    return premise_id


# Module-level cache of the curated AbSplice TSV keyed on variant_key.
_ABSPLICE_TSV = REPO / "data" / "variants" / "absplice_dmd_variants.tsv"
_ABSPLICE_CACHE: dict[str, dict] | None = None


def _load_absplice_table() -> dict[str, dict]:
    global _ABSPLICE_CACHE
    if _ABSPLICE_CACHE is not None:
        return _ABSPLICE_CACHE
    rows: dict[str, dict] = {}
    if not _ABSPLICE_TSV.exists():
        _ABSPLICE_CACHE = rows
        return rows
    with _ABSPLICE_TSV.open() as f:
        header: list[str] | None = None
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
                continue
            r = dict(zip(header, parts))
            for k in ("muscle_skeletal", "heart_lv", "brain_cortex", "retina_proxy", "max_score"):
                if k in r and r[k] != "":
                    try: r[k] = float(r[k])
                    except ValueError: pass
            rows[r["variant_key"]] = r
    _ABSPLICE_CACHE = rows
    return rows


def emit_absplice_premise(conn, cohort: str, pid: str, variant_key: str) -> str | None:
    """Patient-scope AbSplice premise: per-tissue aberrant-splicing probability.

    Distinguishes canonical splice-site variants (top score ~0.85) from
    exonic ones (~0.02) — a signal isoform_arch and clinvar_nmd don't see.
    Skips patients whose variant isn't in the curated table."""
    table = _load_absplice_table()
    row = table.get(variant_key)
    if not row: return None
    ev = {
        "variant_key":      variant_key,
        "hgvsc":            row.get("hgvsc"),
        "category":         row.get("category"),
        "tissue_scores": {
            "Muscle_Skeletal":      row.get("muscle_skeletal"),
            "Heart_Left_Ventricle": row.get("heart_lv"),
            "Brain_Cortex":         row.get("brain_cortex"),
            "Retina_proxy":         row.get("retina_proxy"),
        },
        "max_score":        row.get("max_score"),
        "max_tissue":       row.get("max_tissue"),
        "notes":            row.get("notes"),
        "confidence":       row.get("confidence"),
        "compute":          "off-box (AbSplice v1.0.4 estimates pending)",
    }
    # Baseline confidence lower when scores are estimates rather than an
    # actual model run.
    conf = 0.55 if row.get("confidence") == "estimated" else 0.85
    premise_id = _pid("absplice", cohort, pid)
    conn.execute(
        "INSERT OR REPLACE INTO premise VALUES (?,?,?,?,?,?,?)",
        (premise_id, "absplice", "patient", pid, json.dumps(ev), conf,
         _prov(str(_ABSPLICE_TSV.relative_to(REPO)))),
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


def emit_celltype_impact_premise(conn, cohort: str, pid: str, exon_n: int | None) -> str | None:
    """Composition premise: for each HPA cell type, is it hit / partial /
    spared given the variant's isoform-hit pattern + the curated cell-to-
    isoform dependency map. Closes the variant→cellType chain gap that
    HPA (gene-scoped) and isoform_arch (ends at protein) together leave open."""
    if exon_n is None: return None
    iso_hit = {r[0]: r[1] <= exon_n
               for r in conn.execute(
                   "SELECT isoform_id, first_shared_exon FROM isoforms")}
    cells = []
    for (name, tissue) in conn.execute(
        "SELECT cell_type, tissue FROM celltype_expression "
        "WHERE gene_symbol='DMD' ORDER BY score DESC"):
        req = CELL_TO_ISOFORMS.get(name, [])
        if not req:
            status, hit_isos, spared_isos = "unknown", [], []
        else:
            hit_isos    = [i for i in req if iso_hit.get(i, True)]
            spared_isos = [i for i in req if not iso_hit.get(i, True)]
            status = "hit" if len(hit_isos) == len(req) \
                     else "spared" if not hit_isos \
                     else "partial"
        cells.append({"name": name, "tissue": tissue, "status": status,
                      "hit_isoforms": hit_isos, "spared_isoforms": spared_isos})
    ev = {"variant_exon": exon_n, "cells": cells,
          "curation": "Muntoni 2003, Pillers 1993, Lidov 1995"}
    pid_ = _pid("celltype_impact", cohort, pid)
    conn.execute(
        "INSERT OR REPLACE INTO premise VALUES (?,?,?,?,?,?,?)",
        (pid_, "patient_celltype_impact", "patient", pid, json.dumps(ev), 0.85, _prov()),
    )
    return pid_


def emit_tissue_impact_premise(conn, cohort: str, pid: str, exon_n: int | None) -> str | None:
    """Composition premise: tissues hit / spared for this patient, derived
    from which isoforms are hit × each isoform's primary_expression_tissues."""
    if exon_n is None: return None
    iso_hit = {r[0]: r[1] <= exon_n
               for r in conn.execute(
                   "SELECT isoform_id, first_shared_exon FROM isoforms")}
    hit_tissues, all_tissues = set(), set()
    for (iso_id, tissues_str) in conn.execute(
        "SELECT isoform_id, primary_expression_tissues FROM isoforms"):
        tissues = [t.strip() for t in (tissues_str or "").split(";") if t.strip()]
        for t in tissues:
            all_tissues.add(t)
            if iso_hit.get(iso_id, True):
                hit_tissues.add(t)
    spared_tissues = all_tissues - hit_tissues
    ev = {"variant_exon": exon_n,
          "tissues_hit":    sorted(hit_tissues),
          "tissues_spared": sorted(spared_tissues)}
    pid_ = _pid("tissue_impact", cohort, pid)
    conn.execute(
        "INSERT OR REPLACE INTO premise VALUES (?,?,?,?,?,?,?)",
        (pid_, "patient_tissue_impact", "patient", pid, json.dumps(ev), 0.85, _prov()),
    )
    return pid_


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
# Literature migration: hypothesis_chain_edge_evidence → literature premises
# ----------------------------------------------------------------------
# Each unique citation becomes one `literature` premise in the registry.
# For each hypothesis that cites the paper, we record: which chain edge
# (mapped to biological-layer edge via CHAIN_NODE_TO_LAYER), the claim
# text, and the tone (good/warn).

def collect_literature_data(conn) -> dict[str, dict]:
    """Return {citation → {premise_id, claims, per_hyp: {template_id: [claim]}}}.
    Each per_hyp claim has {position: (link_type, layer_from, layer_to),
    text, tone}."""
    lit_data: dict[str, dict] = {}
    rows = conn.execute("""
        SELECT hypothesis_id, from_node, to_node, ord, tone, text, citation
        FROM hypothesis_chain_edge_evidence
        ORDER BY citation, hypothesis_id, ord
    """).fetchall()
    for (hid, from_node, to_node, _ord, tone, text, citation) in rows:
        if not citation: continue
        entry = lit_data.setdefault(citation, {
            "premise_id": _pid("lit", _slug(citation)),
            "claims": [],
            "per_hyp": {},
        })
        layers = CHAIN_NODE_TO_LAYER.get(hid, {})
        lf, lt = layers.get(from_node), layers.get(to_node)
        if not lf or not lt:
            continue
        position = ("node", lf, lf) if lf == lt else ("edge", lf, lt)
        entry["claims"].append({"hyp": hid, "text": text, "tone": tone,
                                 "position": list(position)})
        entry["per_hyp"].setdefault(hid, []).append({
            "position": position, "text": text, "tone": tone,
        })
    return lit_data


def emit_literature_premises(conn, lit_data: dict) -> int:
    """Emit one premise per unique citation. Chain-position attribution
    happens per hypothesis at bake time (see literature_links_for_template)."""
    for citation, entry in lit_data.items():
        ev = {
            "citation":         citation,
            "claims":           entry["claims"],
            "n_hyp_citations":  len(entry["per_hyp"]),
        }
        conn.execute(
            "INSERT OR REPLACE INTO premise VALUES (?,?,?,?,?,?,?)",
            (entry["premise_id"], "literature", "cohort", "DMD",
             json.dumps(ev), 0.9, _prov()),
        )
    return len(lit_data)


def literature_links_for_template(lit_data: dict, tmpl_id: str
                                  ) -> list[tuple[str, float, str, tuple[str, str, str]]]:
    """Return (premise_id, weight, rationale, chain_position) for each
    literature citation attributed to this hypothesis template's chain.

    Also produces "dual attribution": when a citation informs an edge
    landing at or emanating from a specific layer (e.g. protein →
    subcellular, cellType → tissue), the same citation is *also*
    attributed to that layer as node evidence — because papers that
    establish a transition typically also establish the biology at
    the layers they connect. This especially matters for the
    subcellular node, which otherwise has only UniProt data."""
    out = []
    NODE_DUAL_ATTRIBUTION = {"subcellular"}  # extendable to other node gaps
    for citation, entry in lit_data.items():
        for claim in entry["per_hyp"].get(tmpl_id, []):
            weight = 0.5 if claim["tone"] == "good" else -0.3
            text = (claim["text"] or "")[:180]
            pos = claim["position"]
            out.append((entry["premise_id"], weight, text, pos))

            # Dual attribution: if this edge connects with a
            # dual-attribution layer, ALSO attribute at that layer's node.
            (link_type, lf, lt) = pos
            if link_type == "edge":
                for endpoint in (lf, lt):
                    if endpoint in NODE_DUAL_ATTRIBUTION:
                        out.append((
                            entry["premise_id"],
                            weight * 0.5,   # half weight for the dual attribution
                            f"[dual] {text}",
                            ("node", endpoint, endpoint),
                        ))
    return out


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
    if source_id == "patient_celltype_impact":
        return [("node", "cellType", "cellType"), ("edge", "protein", "cellType")]
    if source_id == "patient_tissue_impact":
        return [("node", "tissue", "tissue"), ("edge", "cellType", "tissue")]
    if source_id == "uniprot_subcellular":
        # Curated protein-biology at the subcellular node + informs the
        # protein→subcellular edge (where does the protein normally sit).
        return [("node", "subcellular", "subcellular"),
                ("edge", "protein", "subcellular")]
    if source_id == "absplice":
        # AbSplice speaks to the variant→protein transition: cryptic /
        # canonical splice disruption changes the transcript, hence the
        # protein produced. Also anchors the variant node itself
        # (variant-level classification signal).
        return [("node", "variant", "variant"),
                ("edge", "variant", "protein")]
    # `literature` premises carry their per-hypothesis chain positions
    # via literature_links_for_template — no source-level default.
    return []


def premise_weights_for_template(hyp_id: str, patient_premises: dict, cohort_premises: dict
                                  ) -> list[tuple[str, float, str]]:
    """Return [(premise_id, signed_weight, rationale)] linking premises to this hyp."""
    out = []
    zhang_pid = patient_premises.get("zhang")
    iso_pid   = patient_premises.get("iso")
    esm3_pid  = patient_premises.get("esm3")
    lab_pids  = patient_premises.get("labs", [])
    ct_pid    = patient_premises.get("celltype_impact")
    ct_ev     = patient_premises.get("celltype_impact_ev", {})
    tis_pid   = patient_premises.get("tissue_impact")
    tis_ev    = patient_premises.get("tissue_impact_ev", {})
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

    # UniProt subcellular localization anchors the sarcolemma-based DGC
    # mechanism. Supports all mechanisms of protein LOSS (they all depend
    # on knowing where dystrophin normally sits) — weak for H04 which is
    # about distal isoforms with different localization patterns.
    us_pid = cohort_premises.get("uniprot_subcellular")
    if us_pid:
        if hyp_id in ("01", "02", "03"):
            out.append((us_pid, 0.6,
                "UniProt: dystrophin localizes at sarcolemma (peripheral membrane, cytoplasmic face); loss disrupts DGC anchoring"))
        elif hyp_id == "04":
            out.append((us_pid, 0.2,
                "UniProt: primary localization is sarcolemma; distal isoforms (Dp71 in retina/CNS/kidney) have separate localization patterns"))

    # AbSplice: per-tissue aberrant-splicing probability.
    # Weight scales with max_score across tissues:
    #   > 0.5  -> strong splice-mechanism signal; boosts H01/H03 (extra NMD
    #            trigger from cryptic/exon-skipping) and H02 (partial
    #            in-frame product possible from cryptic acceptor)
    #   0.05-0.5 -> moderate secondary effect; small positive for H01/H03
    #   < 0.05 -> negligible; no weight contribution
    absp_pid = patient_premises.get("absplice")
    absp_ev  = patient_premises.get("absplice_ev", {})
    if absp_pid and absp_ev:
        max_score  = absp_ev.get("max_score") or 0.0
        max_tissue = absp_ev.get("max_tissue") or "?"
        category   = absp_ev.get("category") or "?"
        if max_score >= 0.5:
            if hyp_id in ("01", "03"):
                out.append((absp_pid, 0.7,
                    f"AbSplice {max_score:.2f} in {max_tissue} ({category}) — aberrant splicing compounds the NMD/DGC-loss mechanism"))
            elif hyp_id == "02":
                out.append((absp_pid, 0.3,
                    f"AbSplice {max_score:.2f} — cryptic-acceptor use could yield a partial in-frame product"))
            elif hyp_id == "04":
                out.append((absp_pid, 0.1,
                    f"AbSplice {max_score:.2f} — splice disruption affects full-length transcript; distal isoforms less relevant"))
        elif max_score >= 0.05:
            if hyp_id in ("01", "03"):
                out.append((absp_pid, 0.2,
                    f"AbSplice {max_score:.2f} ({category}) — moderate secondary splicing effect"))
        else:
            # Low score is itself evidence: no splice mechanism, so the
            # exonic-consequence interpretation stands unchallenged.
            if hyp_id in ("01", "03"):
                out.append((absp_pid, 0.1,
                    f"AbSplice {max_score:.2f} — negligible splice signal; exonic consequence is the whole story"))

    # Composition premises: cell-type & tissue impact.
    # Weight signs matter — spared distal cell types actively argue against H04.
    if ct_pid and ct_ev:
        cells = ct_ev.get("cells", [])
        muscle_cells = [c for c in cells
                        if any(m in (c.get("tissue") or "").lower()
                               for m in ("muscle", "cardiac", "thymus", "salivary"))]
        distal_cells = [c for c in cells
                        if any(d in (c.get("tissue") or "").lower()
                               for d in ("retina", "adipose"))]
        muscle_hit    = sum(1 for c in muscle_cells if c["status"] == "hit")
        distal_hit    = sum(1 for c in distal_cells if c["status"] == "hit")
        distal_spared = sum(1 for c in distal_cells if c["status"] == "spared")

        if hyp_id in ("01", "03") and muscle_hit > 0:
            out.append((ct_pid, 0.6,
                f"{muscle_hit} muscle-lineage cell types affected — supports muscle-fragility mechanism"))
        if hyp_id == "02" and muscle_hit > 0:
            out.append((ct_pid, 0.3,
                f"{muscle_hit} muscle cell types affected — consistent with partial-function BMD scenario"))
        if hyp_id == "04":
            if distal_hit > 0:
                out.append((ct_pid, 0.8,
                    f"{distal_hit} distal (retinal/adipose) cell types affected — supports distal-isoform-loss mechanism"))
            if distal_spared > 0:
                out.append((ct_pid, -0.5,
                    f"{distal_spared} distal cell types SPARED — argues against distal-isoform-loss mechanism"))
            if distal_hit == 0 and distal_spared > 0:
                # Full-negative case: all distal cells spared → H04 unsupported
                out.append((ct_pid, -0.3, "no evidence of distal-isoform cell-type effects for this patient"))

    if tis_pid and tis_ev:
        hit_tissues    = tis_ev.get("tissues_hit", [])
        spared_tissues = tis_ev.get("tissues_spared", [])
        muscle_hit = any("muscle" in t for t in hit_tissues)
        cns_hit    = any(t in hit_tissues for t in ("brain", "brain_glia", "cortical_neurons",
                                                     "hippocampus", "cerebellar_Purkinje_cells"))
        renal_hit  = "kidney" in hit_tissues
        retina_hit = any("retina" in t for t in hit_tissues)

        if hyp_id in ("01", "03") and muscle_hit:
            out.append((tis_pid, 0.5, "skeletal + cardiac muscle tissue affected — core H01/H03 phenotype"))
        if hyp_id == "04":
            distal_tissues_hit = sum([cns_hit, renal_hit, retina_hit])
            if distal_tissues_hit > 0:
                out.append((tis_pid, 0.7,
                    f"{distal_tissues_hit} distal tissue systems (CNS/kidney/retina) affected — supports H04"))
            else:
                out.append((tis_pid, -0.4,
                    "no distal tissue systems affected — argues against H04"))

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
                                 aa: str | None, cons: str | None, acmg: str | None,
                                 lit_data: dict | None = None
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
    ct_pid    = emit_celltype_impact_premise(conn, cohort, pid, exon_n)
    tis_pid   = emit_tissue_impact_premise(conn, cohort, pid, exon_n)
    absp_pid  = emit_absplice_premise(conn, cohort, pid, variant_key)

    # Load evidence blobs for the composition premises so premise-weighting
    # can inspect their content (which cells hit/spared, which tissues, etc.).
    def _load_ev(pid_):
        if not pid_: return {}
        r = conn.execute("SELECT evidence FROM premise WHERE premise_id=?", (pid_,)).fetchone()
        try: return json.loads(r[0]) if r else {}
        except Exception: return {}

    patient_premises = {
        "zhang": zhang_pid, "iso": iso_pid, "esm3": esm3_pid, "labs": lab_pids,
        "celltype_impact": ct_pid, "celltype_impact_ev": _load_ev(ct_pid),
        "tissue_impact":   tis_pid, "tissue_impact_ev":   _load_ev(tis_pid),
        "absplice":        absp_pid, "absplice_ev":       _load_ev(absp_pid),
    }
    cohort_premises = {
        "nmd_cohort":          _pid("nmd_cohort", "DMD"),
        "hpa":                 _pid("hpa", "DMD"),
        "reactome":            _pid("reactome", "DMD"),
        "uniprot_subcellular": _pid("uniprot_subcellular", "DMD"),
    }
    lit_data = lit_data or {}

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
        # `layer_evidence` (abs) drives coverage — any-evidence counts.
        # `signed_layer`  (signed) drives aggregate — negative premises
        # (e.g. spared distal cells arguing against H04) reduce the score.
        layer_evidence: dict[str, float] = {layer: 0.0 for layer in LAYERS}
        edge_evidence:  dict[tuple[str, str], float] = {}
        signed_layer:   dict[str, float] = {layer: 0.0 for layer in LAYERS}
        signed_edge:    dict[tuple[str, str], float] = {}
        n_chain_links = 0

        def _link(pr_id, w, rat, positions):
            nonlocal n_chain_links
            for (link_type, lf, lt) in positions:
                conn.execute(
                    "INSERT OR REPLACE INTO hypothesis_chain_link VALUES (?,?,?,?,?,?,?)",
                    (hyp_id, link_type, lf, lt, pr_id, w, rat),
                )
                n_chain_links += 1
                if link_type == "node":
                    layer_evidence[lf] += abs(w)
                    signed_layer[lf]   += w
                else:
                    edge_evidence[(lf, lt)] = edge_evidence.get((lf, lt), 0) + abs(w)
                    signed_edge[(lf, lt)]   = signed_edge.get((lf, lt), 0) + w

        # A. Route the data/model premises via source-level chain positions
        for (pr_id, w, rat) in prem_weights:
            if not pr_id: continue
            conn.execute(
                "INSERT OR REPLACE INTO hypothesis_premise VALUES (?,?,?,?)",
                (hyp_id, pr_id, w, rat),
            )
            source_id, ev = _get_premise_source_and_ev(pr_id)
            _link(pr_id, w, rat, premise_chain_positions(source_id, ev))

        # B. Literature premises — attribution is per-hypothesis (each
        # citation carries its own chain position for THIS template).
        for (pr_id, w, rat, position) in literature_links_for_template(lit_data, tmpl_id):
            conn.execute(
                "INSERT OR REPLACE INTO hypothesis_premise VALUES (?,?,?,?)",
                (hyp_id, pr_id, w, rat),
            )
            _link(pr_id, w, rat, [position])

        # Compute score vector:
        #   aggregate   = sum of all evidence weights (nodes + edges)
        #   coverage    = fraction of layers with any evidence
        #   consistency = 1.0 (no contradiction detection yet — placeholder)
        #   parsimony   = 1 / (1 + number of empty layers) — chains with
        #                 fewer gaps score higher
        n_layers_covered = sum(1 for v in layer_evidence.values() if v > 0)
        n_layers_empty   = len(LAYERS) - n_layers_covered
        # Aggregate = signed sum (negative premises reduce score).
        aggregate = sum(signed_layer.values()) + sum(signed_edge.values())
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
    emit_uniprot_subcellular_premise(conn, "P11532")
    print("[cohort]   HPA + Reactome + NMD + UniProt-subcellular cohort premises emitted")

    # 2b. Migrate curated citations from hypothesis_chain_edge_evidence
    #     into the literature premise source. One premise per unique
    #     citation; per-hypothesis chain-position attribution via
    #     CHAIN_NODE_TO_LAYER.
    lit_data = collect_literature_data(conn)
    n_lit = emit_literature_premises(conn, lit_data)
    print(f"[literature] {n_lit} unique citations promoted to premises")

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
            conn, cohort, pid, phen, age, amb, exon, nuc, aa, cons, acmg,
            lit_data=lit_data)
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
