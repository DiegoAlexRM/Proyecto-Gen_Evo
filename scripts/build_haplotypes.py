#!/usr/bin/env python3

from Bio import SeqIO
from collections import defaultdict
import sys

if len(sys.argv) != 3:
    print("Uso: python build_haplotypes.py <input.fasta> <prefix_salida>")
    sys.exit(1)

fasta_in = sys.argv[1]
prefix = sys.argv[2]

records = list(SeqIO.parse(fasta_in, "fasta"))

haplos = defaultdict(list)

for rec in records:
    seq = str(rec.seq).upper()
    haplos[seq].append(rec.id)

fasta_out = f"{prefix}_haplotypes.fasta"
ids_out = f"{prefix}_haplotypes_ids.tsv"

# FASTA con representativos
with open(fasta_out, "w") as f:
    for i, (seq, ids) in enumerate(haplos.items(), 1):
        hap_id = f"{prefix}_hap{i}"
        f.write(f">{hap_id}\n{seq}\n")

# Tabla IDs
with open(ids_out, "w") as f:
    f.write("hap_id\tsequence_ids\n")
    for i, (seq, ids) in enumerate(haplos.items(), 1):
        hap_id = f"{prefix}_hap{i}"
        f.write(f"{hap_id}\t{','.join(ids)}\n")

print(f" Haplotipos generados: {len(haplos)}")
print(f" FASTA: {fasta_out}")
print(f" IDs:   {ids_out}")
