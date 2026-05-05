#!/usr/bin/env python3
"""
annotate_unannotated.py
=======================
Annota i geni di B. lanceolatum senza ortologo umano usando DIAMOND
blastp contro Swiss-Prot (tutti gli organismi).

WORKFLOW
--------
  1. Scarica Swiss-Prot da UniProt (~250MB)
  2. Estrae le sequenze proteiche dei geni senza ortologo da klBraLanc5.faa
  3. Crea database DIAMOND di Swiss-Prot
  4. Lancia DIAMOND blastp
  5. Produce tabella con best hit per ogni gene

INPUT
-----
  - symbol_to_human_ortholog.csv  (output di OrthoFinder)
  - klBraLanc5.faa                (proteoma di B. lanceolatum)
  - GCF_035083965_1_klBraLanc5_hap2_feature_table.txt

OUTPUT
------
  - unannotated_genes_swissprot.csv
    Colonne: symbol, query_protein, swissprot_hit, identity, coverage,
             evalue, bitscore, protein_name, organism

UTILIZZO
--------
  python annotate_unannotated.py

  # Per specificare percorsi diversi:
  python annotate_unannotated.py \
      --ortholog symbol_to_human_ortholog.csv \
      --faa ~/OrthoFinderbis/orthofinder_run/proteomes/klBraLanc5.faa \
      --feature GCF_035083965_1_klBraLanc5_hap2_feature_table.txt

DIPENDENZE
----------
  conda activate orthofinder_env
  (DIAMOND è già installato nell'ambiente)
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

# ──────────────────────────────────────────────────────────────────────────────
# CONFIGURAZIONE DEFAULT
# ──────────────────────────────────────────────────────────────────────────────

DEFAULT_ORTHOLOG  = "symbol_to_human_ortholog.csv"
DEFAULT_FAA       = os.path.expanduser("~/OrthoFinderbis/orthofinder_run/proteomes/klBraLanc5.faa")
DEFAULT_FEATURE   = "GCF_035083965_1_klBraLanc5_hap2_feature_table.txt"
DEFAULT_WORKDIR   = Path("swissprot_annotation")
SWISSPROT_URL     = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
N_THREADS         = 10
TOP_HITS          = 1    # best hit per gene (aumenta a 3 per avere più opzioni)
EVALUE_THRESHOLD  = 1e-5

# ──────────────────────────────────────────────────────────────────────────────
# UTILITÀ
# ──────────────────────────────────────────────────────────────────────────────

def run(cmd, check=True):
    print(f"\n$ {cmd}")
    result = subprocess.run(cmd, shell=True, text=True)
    if check and result.returncode != 0:
        print(f"[ERRORE] Comando uscito con codice {result.returncode}")
        sys.exit(result.returncode)
    return result


def require(tool):
    if shutil.which(tool) is None:
        print(f"[ERRORE] '{tool}' non trovato nel PATH.")
        print(f"  Attiva l'ambiente conda: conda activate orthofinder_env")
        sys.exit(1)


# ──────────────────────────────────────────────────────────────────────────────
# STEP 1 — Trova geni senza ortologo umano
# ──────────────────────────────────────────────────────────────────────────────

def get_unannotated_genes(ortholog_csv, feature_table):
    print("\n[STEP 1] Identifico geni senza ortologo umano...")

    ortho = pd.read_csv(ortholog_csv, dtype=str)
    unannotated = ortho[ortho["human_gene_name"].isna()].copy()
    print(f"  Geni senza ortologo umano: {len(unannotated):,}")

    # Ogni gene può avere più proteine (isoforme) — prendiamo la prima XP_
    # che è quella con cui OrthoFinder ha lavorato
    unannotated["best_xp"] = unannotated["query_proteins"].str.split(",").str[0].str.strip()

    # Mappa symbol -> best XP_
    symbol_to_xp = dict(zip(unannotated["symbol"], unannotated["best_xp"]))
    print(f"  Proteine da cercare: {len(symbol_to_xp):,}")

    return unannotated, symbol_to_xp


# ──────────────────────────────────────────────────────────────────────────────
# STEP 2 — Estrai sequenze da klBraLanc5.faa
# ──────────────────────────────────────────────────────────────────────────────

def extract_sequences(faa_path, xp_set, output_faa):
    print(f"\n[STEP 2] Estraggo sequenze da {faa_path}...")

    faa_path = Path(faa_path)
    if not faa_path.exists():
        print(f"[ERRORE] File non trovato: {faa_path}")
        sys.exit(1)

    found = 0
    writing = False
    with open(faa_path) as fin, open(output_faa, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                # Header formato: >XP_066272809.1 protein name [Branchiostoma lanceolatum]
                acc = line.split()[0][1:]  # rimuovi ">"
                acc_nov = re.sub(r"\.\d+$", "", acc)
                if acc in xp_set or acc_nov in xp_set:
                    writing = True
                    fout.write(line)
                    found += 1
                else:
                    writing = False
            elif writing:
                fout.write(line)

    print(f"  Sequenze estratte: {found:,} / {len(xp_set):,}")
    if found == 0:
        print("[ERRORE] Nessuna sequenza trovata. Controlla il percorso del file .faa")
        sys.exit(1)
    return found


# ──────────────────────────────────────────────────────────────────────────────
# STEP 3 — Scarica Swiss-Prot e crea database DIAMOND
# ──────────────────────────────────────────────────────────────────────────────

def prepare_swissprot(workdir):
    print(f"\n[STEP 3] Preparo Swiss-Prot...")

    sp_gz   = workdir / "uniprot_sprot.fasta.gz"
    sp_faa  = workdir / "uniprot_sprot.fasta"
    sp_db   = workdir / "swissprot_db"

    # Download
    if not sp_faa.exists() and not sp_gz.exists():
        print(f"  Download Swiss-Prot (~250MB)...")
        run(f"curl -L '{SWISSPROT_URL}' -o {sp_gz}")
    else:
        print(f"  [SKIP] Swiss-Prot già presente.")

    # Decomprimi
    if not sp_faa.exists():
        print(f"  Decompressione Swiss-Prot...")
        run(f"gzip -d {sp_gz}")
    else:
        print(f"  [SKIP] Swiss-Prot già decompresso.")

    # Crea database DIAMOND
    if not Path(str(sp_db) + ".dmnd").exists():
        print(f"  Creazione database DIAMOND...")
        run(f"diamond makedb --in {sp_faa} -d {sp_db} --threads {N_THREADS}")
    else:
        print(f"  [SKIP] Database DIAMOND già presente.")

    return sp_db


# ──────────────────────────────────────────────────────────────────────────────
# STEP 4 — Lancia DIAMOND blastp
# ──────────────────────────────────────────────────────────────────────────────

def run_diamond(query_faa, sp_db, output_tsv):
    print(f"\n[STEP 4] Lancio DIAMOND blastp...")
    print(f"  Query: {query_faa}")
    print(f"  Database: {sp_db}")
    print(f"  Soglia e-value: {EVALUE_THRESHOLD}")
    print(f"  Thread: {N_THREADS}")

    run(
        f"diamond blastp "
        f"--query {query_faa} "
        f"--db {sp_db} "
        f"--out {output_tsv} "
        f"--outfmt 6 qseqid sseqid pident length qcovhsp evalue bitscore stitle "
        f"--max-target-seqs {TOP_HITS} "
        f"--evalue {EVALUE_THRESHOLD} "
        f"--more-sensitive "
        f"--threads {N_THREADS} "
        f"--quiet"
    )
    print(f"  Output: {output_tsv}")


# ──────────────────────────────────────────────────────────────────────────────
# STEP 5 — Parsing risultati e produzione tabella finale
# ──────────────────────────────────────────────────────────────────────────────

def parse_results(diamond_tsv, unannotated_df, symbol_to_xp, output_csv):
    print(f"\n[STEP 5] Parsing risultati DIAMOND...")

    cols = ["query", "subject", "pident", "length", "qcovhsp",
            "evalue", "bitscore", "stitle"]

    hits = pd.read_csv(diamond_tsv, sep="\t", header=None, names=cols)
    print(f"  Hit trovati: {len(hits):,}")

    # Estrai nome proteina e organismo dal campo stitle
    # Formato Swiss-Prot: sp|P12345|GENE_HUMAN Protein name OS=Homo sapiens OX=9606 GN=GENE PE=1 SV=1
    def parse_stitle(stitle):
        protein_name = stitle
        organism = ""
        gn = ""

        os_match = re.search(r"OS=(.+?)(?:\s+OX=|\s+GN=|\s+PE=|$)", stitle)
        gn_match = re.search(r"GN=(\S+)", stitle)

        if os_match:
            organism = os_match.group(1).strip()
        if gn_match:
            gn = gn_match.group(1).strip()

        # Rimuovi il codice Swiss-Prot dall'inizio
        name_clean = re.sub(r"^sp\|\S+\|\S+\s+", "", stitle)
        name_clean = re.sub(r"\s+OS=.+$", "", name_clean).strip()

        return pd.Series([name_clean, organism, gn])

    hits[["protein_name", "organism", "gene_name_sp"]] = hits["stitle"].apply(parse_stitle)

    # Crea mapping XP_ -> hit
    xp_to_hit = {}
    for _, row in hits.iterrows():
        qp     = row["query"].strip()
        qp_nov = re.sub(r"\.\d+$", "", qp)
        if qp not in xp_to_hit and qp_nov not in xp_to_hit:
            xp_to_hit[qp] = row
            xp_to_hit[qp_nov] = row

    # Costruisci tabella finale
    records = []
    for _, row in unannotated_df.iterrows():
        symbol  = row["symbol"]
        xp      = row["best_xp"]
        xp_nov  = re.sub(r"\.\d+$", "", xp)
        all_xp  = row["query_proteins"]

        hit = xp_to_hit.get(xp) or xp_to_hit.get(xp_nov)

        if hit is not None:
            records.append({
                "symbol":         symbol,
                "query_protein":  xp,
                "swissprot_hit":  hit["subject"],
                "pident":         hit["pident"],
                "qcovhsp":        hit["qcovhsp"],
                "evalue":         hit["evalue"],
                "bitscore":       hit["bitscore"],
                "protein_name":   hit["protein_name"],
                "gene_name_sp":   hit["gene_name_sp"],
                "organism":       hit["organism"],
                "all_xp":         all_xp,
            })
        else:
            records.append({
                "symbol":         symbol,
                "query_protein":  xp,
                "swissprot_hit":  None,
                "pident":         None,
                "qcovhsp":        None,
                "evalue":         None,
                "bitscore":       None,
                "protein_name":   None,
                "gene_name_sp":   None,
                "organism":       None,
                "all_xp":         all_xp,
            })

    result = pd.DataFrame(records)
    n_hit    = result["swissprot_hit"].notna().sum()
    n_nohit  = result["swissprot_hit"].isna().sum()

    print(f"\n  Geni con hit Swiss-Prot   : {n_hit:,} ({n_hit/len(result)*100:.1f}%)")
    print(f"  Geni senza hit            : {n_nohit:,} ({n_nohit/len(result)*100:.1f}%)")
    print(f"  Totale                    : {len(result):,}")

    result.to_csv(output_csv, index=False)
    print(f"\n[OK] Output salvato in: {output_csv}")

    # Anteprima
    print("\n-- Anteprima prime 10 righe con hit --")
    preview_cols = ["symbol", "protein_name", "gene_name_sp", "organism", "pident", "evalue"]
    print(result[result["swissprot_hit"].notna()].head(10)[preview_cols].to_string(index=False))

    return result


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Annota geni senza ortologo umano con DIAMOND vs Swiss-Prot"
    )
    parser.add_argument("--ortholog",  default=DEFAULT_ORTHOLOG,
                        help="File CSV di OrthoFinder (default: symbol_to_human_ortholog.csv)")
    parser.add_argument("--faa",       default=DEFAULT_FAA,
                        help="Proteoma klBraLanc5.faa")
    parser.add_argument("--feature",   default=DEFAULT_FEATURE,
                        help="Feature table NCBI")
    args = parser.parse_args()

    require("diamond")

    # Crea cartella di lavoro
    workdir = DEFAULT_WORKDIR
    workdir.mkdir(exist_ok=True)

    query_faa   = workdir / "unannotated_queries.faa"
    diamond_tsv = workdir / "diamond_swissprot_results.tsv"
    output_csv  = "unannotated_genes_swissprot.csv"

    # Step 1: trova geni senza ortologo
    unannotated_df, symbol_to_xp = get_unannotated_genes(args.ortholog, args.feature)

    # Step 2: estrai sequenze
    xp_set = set(symbol_to_xp.values())
    extract_sequences(args.faa, xp_set, query_faa)

    # Step 3: prepara Swiss-Prot
    sp_db = prepare_swissprot(workdir)

    # Step 4: lancia DIAMOND
    run_diamond(query_faa, sp_db, diamond_tsv)

    # Step 5: parsing e output
    parse_results(diamond_tsv, unannotated_df, symbol_to_xp, output_csv)

    print("\n[PIPELINE COMPLETATA]")


if __name__ == "__main__":
    main()
