#!/bin/bash

# ===============================
# Script Name: glycotestV2.sh
# Description: Prepara alineamientos y sitios N-glicosilados
#              para modelado estructural de Spike
# Author: Dante
# Version: 2.0
# ===============================

usage() {
    echo "==============================================================================="
    echo "        _                 _            _                                    "
    echo "       | |               | |          | |                                   "
    echo "   __ _| | _  _  ___ ___ | |_ ___  ___| |_                                  "
    echo "  / _ \ | | || |/ __/ _ \| __/ _ \/ __| __|                                 "
    echo " | (_| | | |_| | (_| (_) | ||  __/\__ \ |_     v2.0 - Glycoanalysis Prep Script"
    echo "  \__, |_|\__, |\___\___/ \__\___||___/\__|                                 "
    echo "   __/ |   __/ |                                                            "
    echo "  |___/   |___/                                                             "
    echo "-------------------------------------------------------------------------------"
    echo "Uso: $0 -i <archivo.fasta>"
    echo 
    echo "Requisitos (de preferencia en entorno Pixi):"
    echo "  - nextclade"
    echo "  - seqkit"
    echo "  - mafft"
    echo "  - Biopython"
    echo 
    echo "Ejemplo:"
    echo "  ./glycotestV2.sh -i genomas.fasta"
    echo "=========================================="
    exit 0
}

# ===== Lectura de argumentos =====
FASTA_FILE=""
while [[ $# -gt 0 ]]; do
    case $1 in
        -i)
            FASTA_FILE="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo " Opción desconocida: $1"
            usage
            ;;
    esac
done

if [[ -z "$FASTA_FILE" ]]; then
    usage
fi

# ===== Directorios =====
mkdir -p rawdata results scripts results/aa_analysis

echo "=========================================="
echo "1️ Descargando dataset Nextclade..."
pixi run nextclade dataset get --name 'sars-cov-2' --output-dir nextclade_dataset

echo "=========================================="
echo "2️ Ejecutando Nextclade..."
pixi run nextclade run -D nextclade_dataset -O results/nextclade_output "$FASTA_FILE"

echo "=========================================="
echo "3️ Filtrando secuencias sin frameshift en Spike..."
awk -F'\t' 'NR==1 || $30 !~ /S:/' results/nextclade_output/nextclade.tsv > results/nextclade_sinframeshift_S.tsv

echo "=========================================="
echo "4️ Extrayendo secuencias Spike válidas..."
cut -f2 results/nextclade_sinframeshift_S.tsv > results/ids_spike_validos.txt
pixi run seqkit grep -f results/ids_spike_validos.txt results/nextclade_output/nextclade.cds_translation.S.fasta > results/spike_sin_frameshift.fasta

echo "=========================================="
echo "4️.5️ Filtrando secuencias con aminoácidos desconocidos (X)..."

# Conteo inicial de secuencias Spike sin frameshift
TOTAL=$(grep -c "^>" results/spike_sin_frameshift.fasta)

# Eliminar secuencias que tengan al menos una X
pixi run seqkit grep -v -s -r -p "X" \
    results/spike_sin_frameshift.fasta \
    > results/spike_sin_frameshift_noX.fasta

# Conteo después del filtrado
TOTAL_NOX=$(grep -c "^>" results/spike_sin_frameshift_noX.fasta)

ELIMINADAS=$((TOTAL - TOTAL_NOX))

echo "   Secuencias totales (sin frameshift):  $TOTAL"
echo "   Secuencias sin X:                     $TOTAL_NOX"
echo "   Secuencias eliminadas por X:          $ELIMINADAS"

echo "=========================================="
echo "5️ Descargando referencia Spike (Wuhan-Hu-1)..."
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=QHD43416.1&rettype=fasta" > nextclade_dataset/S_reference.fasta

echo "=========================================="
echo "6️ Separando por variantes y alineando..."
variants=("Lambda" "Gamma" "Delta" "Omicron")
for variant in "${variants[@]}"; do
    echo "→ Procesando $variant..."
    awk -F'\t' -v var="$variant" '$5 == var {print $2}' results/nextclade_sinframeshift_S.tsv > results/${variant}_ids.txt
    mkdir -p results/${variant}
    pixi run seqkit grep -f results/${variant}_ids.txt results/spike_sin_frameshift_noX.fasta > results/${variant}/spike_${variant}.fasta
    cat nextclade_dataset/S_reference.fasta results/${variant}/spike_${variant}.fasta > results/${variant}/${variant}_con_ref.fasta
    pixi run mafft --auto results/${variant}/${variant}_con_ref.fasta > results/${variant}/${variant}_alineado.fasta
done

echo "=========================================="
echo "7️ Creando sitios N-glicosilados de referencia..."
python3 scripts/crear_sitios_glyc_ref.py \
    nextclade_dataset/S_reference.fasta \
    scripts/sitios_glyc_ref.txt

echo "=========================================="
echo "8️ Extrayendo tabla de mutaciones..."
for variant in "${variants[@]}"; do
    echo "=== Variante: $variant ==="
    python3 scripts/analisis_mut_aa.py \
        results/$variant/${variant}_alineado.fasta \
        results/$variant/mutaciones_${variant}.tsv
done

echo "=========================================="
echo "8.5️ NxS/T potenciales + ventana ±10 (ALL secuencias) por variante..."
for variant in "${variants[@]}"; do
    echo "=== Variante: $variant ==="
    python3 scripts/analisis_nglyc.py \
      "results/$variant/${variant}_alineado.fasta" \
      "scripts/sitios_glyc_ref.txt" \
      "results/$variant/nglyc_${variant}_ALL.tsv"

    # Copia la tabla ALL a summary (útil para plots globales si quieres)
    cp "results/$variant/nglyc_${variant}_ALL.tsv" "results/summary_tables/"
done

echo "=========================================="
echo "9️ Filtrando IDs con mutaciones ≥10%..."

for variant in "${variants[@]}"; do
    echo "→ Variante $variant"

    # Extraer IDs con frecuencia ≥10%
    awk -F'\t' '$8 >= 10 {print $1}' \
        results/$variant/mutaciones_${variant}.tsv | sort -u \
        > results/$variant/ids_muts10.txt

    echo "   IDs con mutaciones ≥10%: $(wc -l < results/$variant/ids_muts10.txt)"
done

echo "=========================================="
echo "10️ Extrayendo secuencias relevantes (≥1 mutación ≥10%)..."

for variant in "${variants[@]}"; do
    seqkit grep -f results/$variant/ids_muts10.txt \
        results/$variant/spike_${variant}.fasta \
        > results/$variant/spike_${variant}_muts10.fasta
done

echo "=========================================="
echo "1️1️ Limpiando secuencias (removiendo gaps)..."

for variant in "${variants[@]}"; do
    seqkit seq -u results/$variant/spike_${variant}_muts10.fasta \
        > results/$variant/spike_${variant}_muts10_clean.fasta
done

echo "=========================================="
echo "1️2️ Construyendo haplotipos reales (identidad 100%)..."

for variant in "${variants[@]}"; do
    python3 scripts/build_haplotypes.py \
        results/$variant/spike_${variant}_muts10_clean.fasta \
        results/$variant/${variant}
done

echo "=========================================="
echo "1️3️ Construyendo haplotipos dominantes (versión 2)..."

for variant in "${variants[@]}"; do
    python3 scripts/build_haplotypes_version2.py \
        results/$variant/spike_${variant}_muts10_clean.fasta \
        results/$variant/mutaciones_${variant}.tsv \
        results/$variant/${variant}
done

echo "=========================================="
echo "1️3️ Dividiendo haplotipos en FASTAs individuales..."

for variant in "${variants[@]}"; do
    python3 scripts/split_haplotypes_for_modeling.py \
        results/$variant/${variant}_haplotypes_dom.fasta \
        results/$variant/Seq_model \
        $variant
done

echo "=========================================="
echo "13.5️ Correlate FINAL (HAPDOM): mutaciones en ventana NxS/T (tabla lista para R)..."

for variant in "${variants[@]}"; do
    echo "=== Variante: $variant ==="

    IDS_FILE="results/$variant/${variant}_haplotypes_dom_ids.tsv"

    if [[ ! -f "$IDS_FILE" ]]; then
        echo "⚠️  No existe $IDS_FILE. Revisa outputs de build_haplotypes_version2.py"
        continue
    fi

    # 1) Expandir la columna sequence_ids (separada por comas) -> 1 ID por línea
    #    Formato TSV: hap_id <tab> sequence_ids
    awk -F'\t' 'NR==1{next} {print $2}' "$IDS_FILE" \
      | tr ',' '\n' \
      | sed 's/^[[:space:]]*//; s/[[:space:]]*$//' \
      | grep -v '^$' \
      | sort -u \
      > "results/$variant/ids_haplotipos_dom.txt"

    echo "   IDs HAPDOM: $(wc -l < results/$variant/ids_haplotipos_dom.txt)"

    # 2) Filtrar mutaciones a solo esos IDs (ID en columna 1 de mutaciones_*.tsv)
    awk -F'\t' '
      NR==FNR { ids[$1]=1; next }
      FNR==1 { print; next }
      ($1 in ids) { print }
    ' "results/$variant/ids_haplotipos_dom.txt" \
      "results/$variant/mutaciones_${variant}.tsv" \
      > "results/$variant/mutaciones_${variant}_HAPDOM.tsv"

    echo "   Filas mutaciones HAPDOM: $(($(wc -l < results/$variant/mutaciones_${variant}_HAPDOM.tsv) - 1))"

    # 3) Correlate mutaciones (HAPDOM) vs ventanas NxS/T (usa nglyc ALL)
    python3 scripts/correlate_mut_glyc.py \
      "results/$variant/mutaciones_${variant}_HAPDOM.tsv" \
      "results/$variant/nglyc_${variant}_ALL.tsv" \
      "results/$variant/mutaciones_con_ventana_nglyc_${variant}_HAPDOM.tsv"

    # 4) Copias para R
    mkdir -p results/summary_tables
    cp "results/$variant/mutaciones_${variant}_HAPDOM.tsv" "results/summary_tables/"
    cp "results/$variant/mutaciones_con_ventana_nglyc_${variant}_HAPDOM.tsv" "results/summary_tables/"
done

echo "=========================================="
echo "14️ Dividiendo haplotipos dominantes en FASTAs individuales..."
for variant in "${variants[@]}"; do
    python3 scripts/split_haplotypes_for_modeling.py \
      results/$variant/${variant}_haplotypes_dom.fasta \
      results/$variant/Seq_model \
      $variant
done

echo "=========================================="
echo "✅ Listo, Dante."
echo "Tablas listas para R en: results/summary_tables/"
echo "  - nglyc_<var>_ALL.tsv"
echo "  - mutaciones_<var>_HAPDOM.tsv"
echo "  - mutaciones_con_ventana_nglyc_<var>_HAPDOM.tsv"
































