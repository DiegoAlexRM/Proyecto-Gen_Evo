import sys
import os
from Bio import SeqIO

if len(sys.argv) != 4:
    print("Uso: python split_haplotypes_for_modeling.py <haplotypes_dom.fasta> <carpeta_salida> <variante>")
    sys.exit(1)

fasta_in = sys.argv[1]
out_dir = sys.argv[2]
variant = sys.argv[3]

# Crear carpeta de salida
os.makedirs(out_dir, exist_ok=True)

records = list(SeqIO.parse(fasta_in, "fasta"))

print(f"Total de secuencias detectadas: {len(records)}")

for i, rec in enumerate(records, start=1):
    out_name = f"{variant}_hap_{i}.fasta"        # ⬅ nombre final
    out_path = os.path.join(out_dir, out_name)
    SeqIO.write(rec, out_path, "fasta")
    print(f" Exportado: {out_path}")

print("\nListo. Secuencias individuales generadas para modelado.")
