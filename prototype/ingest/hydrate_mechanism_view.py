"""Read data/mechanism.sqlite → write workbench/mechanism_data.json.

Emits a JSON payload matching the shape of the DATA object in
workbench/mechanism.html. Every field is derived from a SELECT; the
query lives next to the assignment so it's obvious what changes when
the substrate tables change.

Run:
    python3 -m prototype.ingest.hydrate_mechanism_view
"""
from __future__ import annotations

import json
import re
import sqlite3
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DB = REPO / "data" / "mechanism.sqlite"
OUT = REPO / "workbench" / "mechanism_data.json"

# Palette hint → CSS var. Kept out of SQL so DB stays presentation-neutral.
CVAR = {
    "accent": "var(--accent)", "violet": "var(--violet)", "teal": "var(--teal)",
    "bad": "var(--bad)", "warn": "var(--warn)", "good": "var(--good)",
    "sky": "var(--sky)", "pink": "var(--pink)", "slate": "var(--slate)",
}


def q_one(conn, sql, *args):
    r = conn.execute(sql, args).fetchone()
    return r[0] if r else None


def build_gene(conn):
    r = conn.execute(
        "SELECT symbol, full_name, uniprot, locus, n_exons, locus_size_mb, isoform_names "
        "FROM gene_meta WHERE symbol='DMD'"
    ).fetchone()
    return {
        "symbol": r[0], "fullName": r[1], "uniprot": r[2], "locus": r[3],
        "nExons": r[4], "locusSizeMb": r[5], "isoformNames": json.loads(r[6]),
    }


def build_header(conn):
    return {
        "variantsAnalyzed":    q_one(conn, "SELECT COUNT(*) FROM lovd_variants"),
        "uniqueVariants":      q_one(conn, "SELECT COUNT(DISTINCT dbid) FROM lovd_variants"),
        "phenotyped":          q_one(conn, "SELECT SUM(n_patients) FROM phenotype_summary "
                                            "WHERE cohort='total_2097'"),
        "mechanismConfidence": int(q_one(conn, "SELECT value FROM settings WHERE key='mechanism_confidence'")),
    }


def build_genetic_evidence(conn):
    total = q_one(conn, "SELECT COUNT(*) FROM lovd_variants")
    by_class = dict(conn.execute(
        "SELECT mut_type_class, COUNT(*) FROM lovd_variants GROUP BY mut_type_class"
    ).fetchall())
    struct_pct = round(100 * by_class.get("structural", 0) / total)
    snv_pct    = round(100 * by_class.get("snv", 0)        / total)
    other_pct  = 100 - struct_pct - snv_pct   # keep sum = 100
    top = conn.execute(
        "SELECT mut_type, ROUND(100.0*COUNT(*)/(SELECT COUNT(*) FROM lovd_variants)) "
        "FROM lovd_variants GROUP BY mut_type ORDER BY 2 DESC LIMIT 4"
    ).fetchall()
    footer = " · ".join(f"{k} {int(v)}%" for k, v in top)
    return {
        "total": total,
        "big": f"{total:,}",
        "bigSub": "LOVD-DMD variant reports",
        "breakdown": [
            {"label": "Structural", "pct": struct_pct, "color": CVAR["accent"]},
            {"label": "SNV",        "pct": snv_pct,    "color": CVAR["violet"]},
            {"label": "Other",      "pct": other_pct,  "color": CVAR["teal"]},
        ],
        "footer": footer,
    }


def build_phenotype_dist(conn):
    # Source: Zhang et al. 2024 (Orphanet J Rare Dis, PMC11344408) reports
    # observed clinical distribution across a 2,097-patient cohort. This
    # is the right substrate for the phenotype tile: ClinVar's per-variant
    # labels are submission-biased (~98% DMD), which misrepresents the
    # actual distribution of dystrophinopathy severity in patients.
    rows = conn.execute(
        "SELECT phenotype_label, n_patients FROM phenotype_summary "
        "WHERE cohort='total_2097' ORDER BY n_patients DESC"
    ).fetchall()
    total = sum(n for _, n in rows)
    color = {"DMD": CVAR["bad"], "BMD": CVAR["teal"], "IMD": CVAR["warn"],
             "pending": CVAR["slate"]}
    label_full = {"DMD": "DMD (Duchenne)", "BMD": "BMD (Becker)",
                  "IMD": "IMD (intermediate)", "pending": "Undetermined"}
    segs = [{"label": label_full.get(l, l), "pct": round(100 * n / total),
             "color": color.get(l, CVAR["slate"])} for l, n in rows]
    top_label, top_n = rows[0]
    return {
        "segments": segs,
        "center": {
            "value": f"{round(100 * top_n / total)}%",
            "label": label_full.get(top_label, top_label),
        },
        "footer": f"Zhang et al. 2024 · N={total:,} patients",
    }


def build_isoform_impact(conn):
    total_exons = q_one(conn, "SELECT n_exons FROM gene_meta WHERE symbol='DMD'")
    rows = conn.execute(
        f"SELECT eu.isoform_id, ROUND(100.0 * SUM(eu.used) / {total_exons}) "
        "FROM exon_usage eu JOIN isoforms i USING (isoform_id) "
        "GROUP BY eu.isoform_id ORDER BY 2 DESC, i.rank"
    ).fetchall()

    def color_for(pct: float) -> str:
        if pct >= 90:  return CVAR["bad"]
        if pct >= 40:  return CVAR["warn"]
        return CVAR["teal"]

    return {
        "rows": [{"label": iso, "value": int(pct), "unit": "%", "color": color_for(pct)}
                 for iso, pct in rows],
        "max": 100,
        "footer": f"Exon coverage of NM_004006.2 ({total_exons} exons)",
    }


def build_hbar(conn, sql: str, footer: str, hard_max: float | None = None):
    rows = conn.execute(sql).fetchall()
    max_v = hard_max if hard_max is not None else (max((r[1] for r in rows), default=1))
    return {
        "rows": [{"label": r[0], "value": r[1], "color": CVAR.get(r[2], CVAR["slate"])}
                 for r in rows],
        "max": max_v,
        "footer": footer,
    }


def build_hypotheses(conn):
    out = []
    for h in conn.execute(
        "SELECT id, rank, name, subtitle, supporting_variants, odds_ratio, "
        "evidence_score, druggability, therapeutic, selected, lede "
        "FROM hypotheses ORDER BY rank"
    ).fetchall():
        hid = h[0]
        row = {
            "id": h[0], "name": h[2], "subtitle": h[3],
            "supporting": h[4], "oddsRatio": h[5], "evidence": h[6],
            "druggability": h[7], "therapeutic": h[8], "selected": bool(h[9]),
        }
        evi = conn.execute(
            "SELECT tone, text, citation FROM hypothesis_evidence "
            "WHERE hypothesis_id=? ORDER BY ord", (hid,)
        ).fetchall()
        nodes = conn.execute(
            "SELECT node_id, col, row, tier, label1, label2, meta "
            "FROM hypothesis_chain_nodes WHERE hypothesis_id=?", (hid,)
        ).fetchall()
        edges = conn.execute(
            "SELECT from_node, to_node FROM hypothesis_chain_edges WHERE hypothesis_id=?",
            (hid,)
        ).fetchall()
        # Group edge evidence by (from, to). Empty list if none curated.
        edge_evi: dict[tuple[str, str], list[dict]] = {}
        for (a, b, tone, text, cite) in conn.execute(
            "SELECT from_node, to_node, tone, text, citation "
            "FROM hypothesis_chain_edge_evidence WHERE hypothesis_id=? "
            "ORDER BY from_node, to_node, ord", (hid,)
        ):
            edge_evi.setdefault((a, b), []).append(
                {"tone": tone, "text": text, "cite": cite})
        therapy = conn.execute(
            "SELECT label1, label2 FROM hypothesis_therapeutic_node WHERE hypothesis_id=?",
            (hid,)
        ).fetchone()

        if h[10] or evi or nodes:  # detail present iff we have lede/evidence/chain
            row["detail"] = {
                "lede": h[10],
                "evidence": [{"tone": t, "text": x, "cite": c} for (t, x, c) in evi],
                "chain": ({
                    "nodes": [{"id": n[0], "col": n[1], "row": n[2], "tier": n[3],
                               "label1": n[4], "label2": n[5], "meta": n[6]} for n in nodes],
                    "edges": [{"from": e[0], "to": e[1],
                               "id": f"{e[0]}-{e[1]}",
                               "evidence": edge_evi.get((e[0], e[1]), [])} for e in edges],
                    "therapeutic": ({"label1": therapy[0], "label2": therapy[1]}
                                    if therapy else None),
                } if nodes else None),
            }
        out.append(row)
    return out


# ======================================================================
# NMD × ACMG cross-tab from ClinVar (for the Marimekko widget).
# ----------------------------------------------------------------------
# For each ClinVar DMD variant, classify:
#   - ACMG: benign | uncertain | pathogenic
#   - NMD:  triggering | escaping | transcript_dependent
#
# NMD classification rules (approximation of Popp & Maquat 2013):
#   - Nonsense (p.XxxNNNTer or p.XxxNNN*)
#       → PTC at codon NNN. NMD-eligible if NNN < 3600 (roughly ex ≤ 78),
#         NMD-escape if NNN ≥ 3600 (last exon of Dp427m at ~aa 3580-3685).
#   - Frameshift (fs) → creates downstream PTC → NMD-triggering
#   - Splice-site (c.*+1/2, c.*-1/2 near splice site) → exon skipping →
#         PTC → NMD-triggering
#   - Missense, synonymous, in-frame indel, promoter deletion, gross del
#         → no PTC → transcript_dependent (protein-level or dose-level effect)
# ======================================================================

DMD_PROTEIN_LEN = 3685              # Dp427m
NMD_ESCAPE_AA_THRESHOLD = 3600      # PTCs beyond this are in the last exon


def classify_variant_nmd(clin_sig: str | None, variant_name: str | None) -> tuple[str, str]:
    cs = (clin_sig or '').lower()
    if 'conflicting' in cs or 'not provided' in cs or 'no classification' in cs or cs.strip() in ('', '-'):
        acmg = 'uncertain'
    elif 'uncertain' in cs:
        acmg = 'uncertain'
    elif 'benign' in cs:
        acmg = 'benign'
    elif 'pathogenic' in cs:
        acmg = 'pathogenic'
    else:
        acmg = 'uncertain'

    name = variant_name or ''
    # nonsense: p.(XxxNNNTer) or p.XxxNNN* or p.(XxxNNN*)
    m = re.search(r'p\.\(?[A-Za-z]{3}(\d+)(?:Ter|\*)\)?', name)
    if m:
        pos = int(m.group(1))
        return acmg, ('escape' if pos >= NMD_ESCAPE_AA_THRESHOLD else 'triggering')
    # frameshift
    if re.search(r'fs\*?\d*', name):
        return acmg, 'triggering'
    # canonical splice-site (c.NNN+1..+2 or -1..-2)
    if re.search(r'c\.\d+[+-][12](?![0-9])', name):
        return acmg, 'triggering'
    # everything else — missense, synonymous, in-frame indel, gross del,
    # promoter changes: no PTC per se; effect is transcript-context-dependent
    return acmg, 'transcript_dependent'


def build_nmd_cohort(conn) -> dict:
    NMD_CATS = ('triggering', 'transcript_dependent', 'escape')
    ACMG_CATS = ('benign', 'uncertain', 'pathogenic')
    grid = {a: {n: 0 for n in NMD_CATS} for a in ACMG_CATS}
    total = 0
    for (cs, nm) in conn.execute(
        "SELECT clin_sig, variant_name FROM clinvar_phenotype"
    ):
        acmg, nmd = classify_variant_nmd(cs, nm)
        grid[acmg][nmd] += 1
        total += 1
    acmg_totals = {a: sum(grid[a].values()) for a in ACMG_CATS}
    nmd_totals  = {n: sum(grid[a][n] for a in ACMG_CATS) for n in NMD_CATS}
    return {
        "source": "ClinVar DMD variant subset",
        "total": total,
        "acmgOrder": list(ACMG_CATS),
        "nmdOrder":  list(NMD_CATS),
        "acmgTotals": acmg_totals,
        "nmdTotals":  nmd_totals,
        "grid": grid,
        "rules": {
            "nmd_escape_aa_threshold": NMD_ESCAPE_AA_THRESHOLD,
            "protein_length": DMD_PROTEIN_LEN,
            "notes": ("Approximation of Popp & Maquat 2013: PTCs upstream of "
                      "the last exon-exon junction trigger NMD; late PTCs escape. "
                      "Splice-site variants are assumed to trigger NMD via "
                      "exon skipping. In-frame and missense variants are "
                      "'transcript-dependent' — no PTC generated."),
        },
    }


def main() -> None:
    conn = sqlite3.connect(DB)
    data = {
        "gene":   build_gene(conn),
        "header": build_header(conn),
        "tiles": {
            "geneticEvidence": build_genetic_evidence(conn),
            "phenotypeDist":   build_phenotype_dist(conn),
            "isoformImpact":   build_isoform_impact(conn),
            "cellTypes": build_hbar(conn,
                "SELECT cell_type, score, color_hint FROM celltype_expression "
                "WHERE gene_symbol='DMD' ORDER BY score DESC LIMIT 7",
                "Human Protein Atlas · nCPM (log-scaled)",
                hard_max=10),
            "pathways": build_hbar(conn,
                "SELECT pathway_name, score, color_hint FROM pathway_enrichment "
                "WHERE gene_symbol='DMD' ORDER BY score DESC LIMIT 7",
                "Reactome membership · specificity = 10 − log(n_genes)",
                hard_max=10),
        },
        "hypotheses": build_hypotheses(conn),
        "hypothesesMeta": {
            "top":   q_one(conn, "SELECT COUNT(*) FROM hypotheses"),
            "total": q_one(conn, "SELECT COUNT(*) FROM lovd_variants"),
            "unit":  "hypotheses shown",
        },
        "nmdCohort": build_nmd_cohort(conn),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(data, indent=2, ensure_ascii=False))
    print(f"[wrote] {OUT} ({OUT.stat().st_size / 1024:.1f} KB)")
    conn.close()


if __name__ == "__main__":
    main()
