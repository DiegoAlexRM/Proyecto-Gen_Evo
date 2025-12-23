import pandas as pd
import sys

if len(sys.argv) != 4:
    print("Uso: python correlate_mut_glyc.py resume_mutaciones.tsv resumen_glyc.tsv salida.tsv")
    sys.exit(1)

archivo_mut = sys.argv[1]
archivo_glyc = sys.argv[2]
archivo_salida = sys.argv[3]

# Leer archivos
mut = pd.read_csv(archivo_mut, sep="\t")
glyc = pd.read_csv(archivo_glyc, sep="\t")

# Extraer solo las columnas necesarias de sitios de glicosilación
sitios = glyc[["seq_id", "sitio_pos"]].drop_duplicates().rename(columns={
    "seq_id": "secuencia",
    "sitio_pos": "posicion"
})

# Añadir columna de marcador
mut["en_ventana_glicosilada"] = "no"

# Comparar posiciones dentro de ±10
for i, row in mut.iterrows():
    id_ = row["id"]
    pos = row["pos"]
    sitios_locales = sitios[sitios["secuencia"] == id_]

    for _, sitio in sitios_locales.iterrows():
        if abs(pos - sitio["posicion"]) <= 10:
            mut.at[i, "en_ventana_glicosilada"] = "si"
            break

# Añadir columnas de conteo y frecuencia para cada mutación (mut)
conteo_mut = mut["mut"].value_counts().to_dict()
total_secuencias = mut["id"].nunique()

mut["conteo"] = mut["mut"].map(conteo_mut)
mut["frecuencia(%)"] = mut["conteo"].apply(lambda x: (x / total_secuencias) * 100)

# Guardar
mut.to_csv(archivo_salida, sep="\t", index=False)
print(f" Resultado guardado en: {archivo_salida}")
