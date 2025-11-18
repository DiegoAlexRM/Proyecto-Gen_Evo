from Bio import SeqIO
import re
import sys

if len(sys.argv) != 3:
    print("Uso: python crear_sitios_glyc_ref.py S_reference.fasta sitios_glyc_ref.txt")
    sys.exit(1)

fasta_file = sys.argv[1]
output_file = sys.argv[2]

# Cargar referencia
ref = next(SeqIO.parse(fasta_file, "fasta"))
seq = str(ref.seq)

# Buscar motivos N[^P][ST]
with open(output_file, "w") as out:
    for m in re.finditer(r'N[^P][ST]', seq):
        out.write(f"{m.start()+1} {m.group()}\n")