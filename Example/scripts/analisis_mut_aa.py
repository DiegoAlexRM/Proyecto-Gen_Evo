import sys
from Bio import SeqIO
import csv
from collections import defaultdict

if len(sys.argv) != 3:
    print("Uso: python analizar_mut_aa.py archivo.fasta salida.tsv")
    sys.exit(1)

fasta_file = sys.argv[1]
archivo_salida = sys.argv[2]

# Diccionarios de polaridad y carga
aa_polaridad = {
    'A': 'apolar', 'V': 'apolar', 'L': 'apolar', 'I': 'apolar', 'M': 'apolar',
    'F': 'apolar', 'W': 'apolar', 'P': 'apolar', 'G': 'apolar',
    'S': 'polar', 'T': 'polar', 'C': 'polar', 'N': 'polar', 'Q': 'polar', 'Y': 'polar',
    'D': 'polar', 'E': 'polar', 'K': 'polar', 'R': 'polar', 'H': 'polar',
    'U': 'polar', 'O': 'polar', 'B': 'polar', 'Z': 'polar',
    'X': 'unknown', '*': 'unknown', '-': 'unknown'
}

aa_carga = {
    'D': 'negativa', 'E': 'negativa',
    'K': 'positiva', 'R': 'positiva', 'H': 'positiva',
    'A': 'neutra', 'V': 'neutra', 'L': 'neutra', 'I': 'neutra', 'M': 'neutra',
    'F': 'neutra', 'Y': 'neutra', 'W': 'neutra',
    'P': 'neutra', 'G': 'neutra',
    'S': 'neutra', 'T': 'neutra', 'C': 'neutra', 'N': 'neutra', 'Q': 'neutra',
    'U': 'neutra', 'O': 'neutra', 'B': 'neutra', 'Z': 'neutra',
    'X': 'desconocida', '*': 'desconocida', '-': 'desconocida'
}

# Leer alineamiento
alineadas = list(SeqIO.parse(fasta_file, "fasta"))

# Detectar referencia
try:
    ref_record = next(rec for rec in alineadas if 'QHD43416.1' in rec.id)
except StopIteration:
    print("⚠️ No se encontró la referencia 'QHD43416.1' en el archivo.")
    sys.exit(1)

ref_seq = str(ref_record.seq)
ref_id = ref_record.id

# Diccionario para contar mutaciones globales
mutacion_global = defaultdict(int)

# Lista para almacenar filas
filas = []

total_secuencias = 0

for rec in alineadas:
    if rec.id == ref_id:
        continue
    total_secuencias += 1
    query_seq = str(rec.seq)
    for i, (ref_aa, q_aa) in enumerate(zip(ref_seq, query_seq), start=1):

        # 👉 saltar posiciones donde la referencia tiene un gap
        if ref_aa == '-':
            continue

        if ref_aa != q_aa and ref_aa not in ['-', 'X'] and q_aa not in ['-', 'X']:
            ref_polar = aa_polaridad.get(ref_aa, 'unknown')
            mut_polar = aa_polaridad.get(q_aa, 'unknown')
            ref_carga = aa_carga.get(ref_aa, 'desconocida')
            mut_carga = aa_carga.get(q_aa, 'desconocida')

            polaridad_cambio = f"{ref_polar} , {mut_polar}"
            carga_cambio = f"{ref_carga} , {mut_carga}"
            cambio_fq = 'si' if (ref_polar != mut_polar or ref_carga != mut_carga) else 'no'

            mut = f"{ref_aa}{i}{q_aa}"
            mutacion_global[mut] += 1
            filas.append([rec.id, i, mut, polaridad_cambio, carga_cambio, cambio_fq])

# Abrir archivo de salida
with open(archivo_salida, "w", newline="") as out_file:
    writer = csv.writer(out_file, delimiter="\t")
    writer.writerow(["id", "pos", "mut", "polaridad_cambio", "carga_cambio", "cambio_fisicoquimico", "conteo", "frecuencia(%)"])

    for fila in filas:
        mut = fila[2]
        conteo = mutacion_global[mut]
        frecuencia = (conteo / total_secuencias) * 100
        writer.writerow(fila + [conteo, f"{frecuencia:.2f}"])

print(f"✔️ Mutaciones exportadas a: {archivo_salida}")
