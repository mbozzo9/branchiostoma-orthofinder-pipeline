#!/usr/bin/env python3
"""
build_NP_mapping.py
===================
Legge gene2refseq (decompresso) e Homo_sapiens.gene_info.gz
e produce un file NP_ -> gene_symbol per umani.

UTILIZZO:
    python build_NP_mapping.py
"""

import gzip
import csv
import os
from pathlib import Path

GENE2REFSEQ  = Path("orthofinder_run/gene2refseq")
GENE_INFO    = Path("orthofinder_run/Homo_sapiens.gene_info.gz")
OUTPUT       = Path("orthofinder_run/human_NP_to_symbol.txt")

# ── Step 1: carica GeneID → Symbol da gene_info ──────────────────────────────
print("Carico gene_info...")
geneid_to_symbol = {}
with gzip.open(GENE_INFO, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        geneid  = parts[1]
        symbol  = parts[2]
        geneid_to_symbol[geneid] = symbol

print(f"  GeneID umani caricati: {len(geneid_to_symbol):,}")

# ── Step 2: leggi gene2refseq e filtra NP_ umane ─────────────────────────────
print("Leggo gene2refseq e filtro NP_ umane...")
print("  (questo può richiedere qualche minuto)")

count = 0
with open(GENE2REFSEQ, "r") as fin, open(OUTPUT, "w") as fout:
    for i, line in enumerate(fin):
        if i == 0:
            continue  # header
        if i % 1_000_000 == 0:
            print(f"  Righe processate: {i:,}...")
        
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        
        tax_id  = parts[0]
        geneid  = parts[1]
        np_acc  = parts[5]  # protein accession
        
        if tax_id != "9606":
            continue
        if not np_acc.startswith("NP_"):
            continue
        
        symbol = geneid_to_symbol.get(geneid)
        if symbol:
            fout.write(f"{np_acc}\t{symbol}\n")
            count += 1

print(f"\n[OK] Scritte {count:,} coppie NP_ -> symbol in {OUTPUT}")
print("Prime 5 righe:")
with open(OUTPUT) as f:
    for _ in range(5):
        print(" ", f.readline().strip())
