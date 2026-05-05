#!/usr/bin/env python3
"""
merge_annotations_v2.py
=======================
Unisce annotazioni OrthoFinder e Swiss-Prot in un unico file completo.

UTILIZZO: python merge_annotations_v2.py
"""

import re
import sys
from pathlib import Path
import pandas as pd

ORTHOFINDER_CSV = "symbol_to_human_ortholog.csv"
SWISSPROT_CSV   = "unannotated_genes_swissprot.csv"
OUTPUT_CSV      = "annotation_complete.csv"

METAZOA = [
    "sapiens", "musculus", "rattus", "familiaris", "taurus", "scrofa",
    "gallus", "danio", "xenopus", "drosophila", "caenorhabditis",
    "strongylocentrotus", "ciona", "branchiostoma", "petromyzon",
    "latimeria", "oryzias", "takifugu", "tetraodon", "anolis",
    "monodelphis", "macaca", "pan troglodytes", "gorilla", "pongo"
]

def is_metazoan(organism):
    org = str(organism).lower()
    return any(k in org for k in METAZOA)

def normalize_gene(gene, organism):
    if pd.isna(gene) or str(gene) == "nan" or str(gene) == "":
        return None
    if not is_metazoan(organism):
        return None
    return str(gene).upper()

def main():
    print("Carico OrthoFinder...")
    ortho = pd.read_csv(ORTHOFINDER_CSV, dtype=str)
    print(f"  Totale: {len(ortho):,} | Con ortologo: {ortho['human_gene_name'].notna().sum():,}")

    print("Carico Swiss-Prot...")
    sp = pd.read_csv(SWISSPROT_CSV, dtype=str)
    print(f"  Totale: {len(sp):,} | Con hit: {sp['swissprot_hit'].notna().sum():,}")

    # Aggiungi query_proteins a sp
    qp_map = ortho.set_index("symbol")["query_proteins"].to_dict()
    sp["query_proteins"] = sp["symbol"].map(qp_map)

    # Normalizza gene name Swiss-Prot
    sp["human_gene_name"] = sp.apply(
        lambda r: normalize_gene(r.get("gene_name_sp"), r.get("organism", "")), axis=1
    )

    n_transfer = sp["human_gene_name"].notna().sum()
    print(f"  Hit trasferibili da metazoi: {n_transfer:,}")

    # ── OrthoFinder annotati ──────────────────────────────────────────────────
    oa = ortho[ortho["human_gene_name"].notna()].copy()
    oa["annotation_source"] = "OrthoFinder"
    oa["swissprot_hit"] = None
    oa["protein_name"]  = None
    oa["organism_sp"]   = None
    oa["pident"]        = None
    oa["evalue"]        = None

    # ── Swiss-Prot annotati ───────────────────────────────────────────────────
    sa = sp[sp["human_gene_name"].notna()].copy()
    sa["annotation_source"] = "SwissProt_homolog"
    sa["human_protein_id"]  = None
    sa["orthogroup"]        = None
    sa = sa.rename(columns={"organism": "organism_sp"})

    # ── Unknown ───────────────────────────────────────────────────────────────
    ortho_missing = set(ortho[ortho["human_gene_name"].isna()]["symbol"])
    sp_missing    = set(sp[sp["human_gene_name"].isna()]["symbol"])
    unknown_syms  = ortho_missing & sp_missing

    ua = pd.DataFrame({"symbol": sorted(unknown_syms)})
    ua["annotation_source"] = "unknown"
    ua["human_gene_name"]   = None
    ua["human_protein_id"]  = None
    ua["orthogroup"]        = None
    ua["swissprot_hit"]     = None
    ua["protein_name"]      = None
    ua["organism_sp"]       = None
    ua["pident"]            = None
    ua["evalue"]            = None
    ua["query_proteins"]    = ua["symbol"].map(qp_map)

    # ── Colonne finali ────────────────────────────────────────────────────────
    cols = [
        "symbol", "human_gene_name", "annotation_source",
        "human_protein_id", "orthogroup",
        "swissprot_hit", "protein_name", "organism_sp",
        "pident", "evalue", "query_proteins"
    ]

    # Assicura che tutte le colonne esistano
    for df in [oa, sa, ua]:
        for c in cols:
            if c not in df.columns:
                df[c] = None

    result = pd.concat([oa[cols], sa[cols], ua[cols]], ignore_index=True)
    result = result.sort_values("symbol").reset_index(drop=True)

    n_o = (result["annotation_source"] == "OrthoFinder").sum()
    n_s = (result["annotation_source"] == "SwissProt_homolog").sum()
    n_u = (result["annotation_source"] == "unknown").sum()
    total = len(result)

    print(f"\n{'='*50}")
    print(f" RISULTATI FINALI")
    print(f"{'='*50}")
    print(f"  OrthoFinder       : {n_o:,} ({n_o/total*100:.1f}%)")
    print(f"  SwissProt_homolog : {n_s:,} ({n_s/total*100:.1f}%)")
    print(f"  Unknown           : {n_u:,} ({n_u/total*100:.1f}%)")
    print(f"  Totale            : {total:,}")

    result.to_csv(OUTPUT_CSV, index=False)
    print(f"\n[OK] Salvato: {OUTPUT_CSV}")

    print("\n-- Anteprima SwissProt_homolog (prime 5 righe) --")
    preview = result[result["annotation_source"] == "SwissProt_homolog"][
        ["symbol", "human_gene_name", "protein_name", "organism_sp", "pident", "evalue"]
    ].head(5)
    print(preview.to_string(index=False))

if __name__ == "__main__":
    main()
