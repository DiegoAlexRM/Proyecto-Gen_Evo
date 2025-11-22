import sys
from Bio import SeqIO

if len(sys.argv) != 4:
    print("Uso: python analizar_sitios_glyc_tsv.py alineado.fasta sitios_glyc_ref.txt salida.tsv")
    sys.exit(1)

fasta_file = sys.argv[1]
sitios_file = sys.argv[2]
output_file = sys.argv[3]

# Leer sitios de glicosilación de referencia
sitios = []
with open(sitios_file) as f:
    for line in f:
        pos, motivo = line.strip().split()
        sitios.append((int(pos), motivo))

# Leer alineamiento
secuencias = list(SeqIO.parse(fasta_file, "fasta"))

# Buscar referencia
ref_record = next((s for s in secuencias if "QHD43416.1" in s.id), None)
if not ref_record:
    print("No se encontró la referencia QHD43416.1 en el archivo.")
    sys.exit(1)

ref_seq = str(ref_record.seq)

# Crear archivo de salida
with open(output_file, "w") as out:
    out.write("seq_id\tsitio_pos\tmotivo_ref\tmotivo_obs\tconservado\tmutacion_directa\tmutaciones_ventana\n")
    
    for s in secuencias:
        if s.id == ref_record.id:
            continue  # saltar referencia

        for pos, motivo_ref in sitios:
            try:
                motivo_query = s.seq[pos-1:pos+2]
                mut_directa = motivo_query != ref_seq[pos-1:pos+2]
                destruido = (
                    motivo_query[0] != 'N' or
                    motivo_query[1] == 'P' or
                    motivo_query[2] not in ['S', 'T']
                )
                conservado = not destruido and 'X' not in motivo_query

                # Mutaciones en ventana de ±10aa
                window_ref = ref_seq[pos-11:pos+10]
                window_seq = s.seq[pos-11:pos+10]
                mutaciones = sum(
                    1 for a, b in zip(window_ref, window_seq)
                    if a != b and b not in ['-', 'X']
                )

                out.write(f"{s.id}\t{pos}\t{motivo_ref}\t{motivo_query}\t{'Si' if conservado else 'No'}\t{'Si' if mut_directa else 'No'}\t{mutaciones}\n")

            except:
                out.write(f"{s.id}\t{pos}\t{motivo_ref}\tERROR\tNo\tNo\tNA\n")
