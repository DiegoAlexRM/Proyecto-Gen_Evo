#!/usr/bin/env python3
"""
build_haplotypes_version2.py
--------------------------------------
Construye haplotipos usando SOLO las mutaciones dominantes (≥10 %).

Entrada:
  1) FASTA limpio (sin gaps), solo secuencias relevantes
  2) TSV de mutaciones (mutaciones_variante.tsv)
  3) FASTA alineado con referencia (para leer aminoácidos reales)

Salida:
  - <prefix>_haplotypes_dom.fasta
  - <prefix>_haplotypes_dom_summary.tsv
  - <prefix>_haplotypes_dom_ids.tsv
  - <prefix>_haplotypes_dom_profiles.tsv
"""

import sys
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

if len(sys.argv) != 4:
    print("Uso: python3 build_haplotypes_version2.py <clean.fasta> <mutaciones.tsv> <prefix>")
    sys.exit(1)

clean_fasta = sys.argv[1]
mut_tsv     = sys.argv[2]
prefix      = sys.argv[3]

# ==============================
# 1. Leer mutaciones y obtener posiciones ≥10%
# ==============================
mut = pd.read_csv(mut_tsv, sep="\t")
dom = mut[mut["frecuencia(%)"] >= 10]

if dom.empty:
    print("No hay mutaciones ≥10%")
    sys.exit(1)

dom_positions = sorted(dom["pos"].unique())
print(f"Posiciones dominantes (≥10%): {dom_positions}")

# Diccionario pos → (refAA →?) (porque en mut TSV ya viene ref_aa en mut ex: T19R)
pos_to_ref = {}
for _, row in dom.iterrows():
    mut_str = row["mut"]  # ejemplo: "T19R"
    ref_aa = mut_str[0]
    pos    = int(row["pos"])
    pos_to_ref[pos] = ref_aa

# ==============================
# 2. Leer FASTA limpio
# ==============================
records = list(SeqIO.parse(clean_fasta, "fasta"))
print(f"📄 Secuencias cargadas: {len(records)}")

# ==============================
# 3. Construir perfiles dominantes
# ==============================
haplos = defaultdict(list)
haplo_profiles = {}

for rec in records:
    seq = str(rec.seq)
    perfil = []

    for pos in dom_positions:
        ref_aa = pos_to_ref[pos]
        q_aa = seq[pos - 1]

        if q_aa == ref_aa:
            mut_str = f"{ref_aa}{pos}{ref_aa}"   # sin mutación dominante
        else:
            mut_str = f"{ref_aa}{pos}{q_aa}"

        perfil.append(mut_str)

    key = ";".join(perfil)
    haplos[key].append(rec.id)

# ==============================
# 4. Escribir FASTA de representantes
# ==============================
fasta_out = f"{prefix}_haplotypes_dom.fasta"
with open(fasta_out, "w") as out:
    for i, (perfil, ids) in enumerate(haplos.items(), 1):
        hap_id = f"{prefix}_hapD{i}"
        rep_seq = str(SeqIO.to_dict(records)[ids[0]].seq)

        out.write(f">{hap_id}\n{rep_seq}\n")
        haplo_profiles[hap_id] = perfil

# ==============================
# 5. Resumen TSV
# ==============================
summary_out = f"{prefix}_haplotypes_dom_summary.tsv"
with open(summary_out, "w") as out:
    out.write("hap_id\tn_seq\tfrecuencia(%)\tmutaciones_dominantes\n")
    total = len(records)
    for i, (perfil, ids) in enumerate(haplos.items(), 1):
        hap_id = f"{prefix}_hapD{i}"
        freq = (len(ids) / total) * 100
        out.write(f"{hap_id}\t{len(ids)}\t{freq:.2f}\t{haplo_profiles[hap_id]}\n")

# ==============================
# 6. IDs por haplotipo
# ==============================
ids_out = f"{prefix}_haplotypes_dom_ids.tsv"
with open(ids_out, "w") as out:
    out.write("hap_id\tsequence_ids\n")
    for i, (perfil, ids) in enumerate(haplos.items(), 1):
        hap_id = f"{prefix}_hapD{i}"
        out.write(f"{hap_id}\t{','.join(ids)}\n")

# ==============================
# 7. Perfiles completos (opcional)
# ==============================
profiles_out = f"{prefix}_haplotypes_dom_profiles.tsv"
with open(profiles_out, "w") as out:
    out.write("hap_id\tperfil_dominante\n")
    for i, (perfil, ids) in enumerate(haplos.items(), 1):
        hap_id = f"{prefix}_hapD{i}"
        out.write(f"{hap_id}\t{perfil}\n")

print("Haplotipos dominantes creados")
print(f"FASTA: {fasta_out}")
print(f"Summary: {summary_out}")
print(f"IDs: {ids_out}")
print(f"Profiles: {profiles_out}")
