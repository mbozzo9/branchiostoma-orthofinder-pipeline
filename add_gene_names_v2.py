#!/usr/bin/env python3
"""
add_gene_names_v2.py
====================
Aggiunge le colonne gene_name e annotation_source a un file Excel di counts
usando l'annotazione completa prodotta da merge_annotations_v2.py.

La colonna annotation_source indica l'origine dell'annotazione:
  - OrthoFinder       : ortologo umano formale
  - SwissProt_homolog : omologo da metazoi via DIAMOND vs Swiss-Prot
  - unknown           : nessuna annotazione trovata

La colonna annotation_source viene posizionata DOPO le colonne dei counts
in modo da non interferire con lo script R (DESeq2 legge solo le colonne
numeriche dei campioni).

INPUT
-----
  - <input>.xlsx               file Excel dei counts
  - annotation_complete.csv    output di merge_annotations_v2.py

OUTPUT
------
  - <input>_annotated.xlsx
  - <input>_annotated.csv

UTILIZZO
--------
  python add_gene_names_v2.py --input gene_count_larvae_DBS.xlsx
  python add_gene_names_v2.py --input gene_count_adults_DBS.xlsx

  # Specificare un file di annotazione diverso:
  python add_gene_names_v2.py --input gene_count_larvae_DBS.xlsx \
                               --annotation annotation_complete.csv
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

DEFAULT_ANNOTATION = "annotation_complete.csv"


def main():
    parser = argparse.ArgumentParser(
        description="Aggiunge gene_name e annotation_source a un file Excel di counts"
    )
    parser.add_argument("--input", required=True,
                        help="File Excel di input (es. gene_count_larvae_DBS.xlsx)")
    parser.add_argument("--annotation", default=DEFAULT_ANNOTATION,
                        help=f"File di annotazione completa (default: {DEFAULT_ANNOTATION})")
    args = parser.parse_args()

    input_path      = Path(args.input)
    annotation_path = Path(args.annotation)

    # ── Verifica file ─────────────────────────────────────────────────────────
    for f in [input_path, annotation_path]:
        if not f.exists():
            print(f"[ERRORE] File non trovato: {f}")
            sys.exit(1)

    # ── Leggi Excel ───────────────────────────────────────────────────────────
    print(f"  Lettura {input_path}...")
    df = pd.read_excel(input_path, dtype=str)

    if "gene_id" not in df.columns:
        print(f"[ERRORE] Colonna 'gene_id' non trovata.")
        sys.exit(1)

    print(f"  Geni nel file: {len(df):,}")

    # Estrai symbol da gene_id (parte prima del "|")
    df["_symbol"] = df["gene_id"].str.split("|").str[0].str.strip()

    # ── Leggi annotazione completa ────────────────────────────────────────────
    print(f"  Lettura {annotation_path}...")
    annot = pd.read_csv(annotation_path, dtype=str)

    symbol_to_name   = annot.set_index("symbol")["human_gene_name"].to_dict()
    symbol_to_source = annot.set_index("symbol")["annotation_source"].to_dict()

    print(f"  Geni annotati nel file: {annot['human_gene_name'].notna().sum():,}")

    # ── Popola colonne ────────────────────────────────────────────────────────
    df["gene_name"]         = df["_symbol"].map(symbol_to_name)
    df["annotation_source"] = df["_symbol"].map(symbol_to_source).fillna("unknown")
    df = df.drop(columns=["_symbol"])

    # ── Ordine colonne ────────────────────────────────────────────────────────
    # gene_id | gene_name | <counts> | annotation_source
    # annotation_source viene messa in fondo per non disturbare R
    count_cols = [c for c in df.columns
                  if c not in ("gene_id", "gene_name", "annotation_source")]

    df = df[["gene_id", "gene_name"] + count_cols + ["annotation_source"]]

    # ── Converti counts in numerici ───────────────────────────────────────────
    for col in count_cols:
        df[col] = pd.to_numeric(df[col], errors="ignore")

    # ── Statistiche ───────────────────────────────────────────────────────────
    n_ortho = (df["annotation_source"] == "OrthoFinder").sum()
    n_sp    = (df["annotation_source"] == "SwissProt_homolog").sum()
    n_unk   = (df["annotation_source"] == "unknown").sum()
    total   = len(df)

    print(f"\n  OrthoFinder       : {n_ortho:,} ({n_ortho/total*100:.1f}%)")
    print(f"  SwissProt_homolog : {n_sp:,} ({n_sp/total*100:.1f}%)")
    print(f"  Unknown           : {n_unk:,} ({n_unk/total*100:.1f}%)")
    print(f"  Totale            : {total:,}")

    # ── Salva output ──────────────────────────────────────────────────────────
    stem     = input_path.stem
    out_xlsx = input_path.parent / f"{stem}_annotated.xlsx"
    out_csv  = input_path.parent / f"{stem}_annotated.csv"

    df.to_excel(out_xlsx, index=False)
    df.to_csv(out_csv, index=False)

    print(f"\n[OK] Output salvato:")
    print(f"  Excel : {out_xlsx}")
    print(f"  CSV   : {out_csv}")

    # ── Anteprima ─────────────────────────────────────────────────────────────
    print("\n-- Anteprima prime 5 righe con gene_name --")
    preview = df[df["gene_name"].notna()][
        ["gene_id", "gene_name", "annotation_source"]
    ].head(5)
    print(preview.to_string(index=False))


if __name__ == "__main__":
    main()
