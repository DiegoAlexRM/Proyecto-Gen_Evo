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
            echo "❌ Opción desconocida: $1"
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
echo "1️⃣ Descargando dataset Nextclade..."
pixi run nextclade dataset get --name 'sars-cov-2' --output-dir nextclade_dataset

echo "=========================================="
echo "2️⃣ Ejecutando Nextclade..."
pixi run nextclade run -D nextclade_dataset -O results/nextclade_output "$FASTA_FILE"

echo "=========================================="
echo "3️⃣ Filtrando secuencias sin frameshift en Spike..."
awk -F'\t' 'NR==1 || $30 !~ /S:/' results/nextclade_output/nextclade.tsv > results/nextclade_sinframeshift_S.tsv

echo "=========================================="
echo "4️⃣ Extrayendo secuencias Spike válidas..."
cut -f2 results/nextclade_sinframeshift_S.tsv > results/ids_spike_validos.txt
pixi run seqkit grep -f results/ids_spike_validos.txt results/nextclade_output/nextclade.cds_translation.S.fasta > results/spike_sin_frameshift.fasta

echo "=========================================="
echo "4️⃣.5️⃣ Filtrando secuencias con aminoácidos desconocidos (X)..."

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
echo "5️⃣ Descargando referencia Spike (Wuhan-Hu-1)..."
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=QHD43416.1&rettype=fasta" > nextclade_dataset/S_reference.fasta

echo "=========================================="
echo "6️⃣ Separando por variantes y alineando..."
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
echo "7️⃣ Creando sitios N-glicosilados de referencia..."
python3 scripts/crear_sitios_glyc_ref.py \
    nextclade_dataset/S_reference.fasta \
    scripts/sitios_glyc_ref.txt

echo "=========================================="
echo "8️⃣ Extrayendo tabla de mutaciones..."
for variant in "${variants[@]}"; do
    echo "=== Variante: $variant ==="
    python3 scripts/analisis_mut_aa.py \
        results/$variant/${variant}_alineado.fasta \
        results/$variant/mutaciones_${variant}.tsv
done

echo "=========================================="
echo "9️⃣ Filtrando IDs con mutaciones ≥10%..."

for variant in "${variants[@]}"; do
    echo "→ Variante $variant"

    # Extraer IDs con frecuencia ≥10%
    awk -F'\t' '$8 >= 10 {print $1}' \
        results/$variant/mutaciones_${variant}.tsv | sort -u \
        > results/$variant/ids_muts10.txt

    echo "   IDs con mutaciones ≥10%: $(wc -l < results/$variant/ids_muts10.txt)"
done

echo "=========================================="
echo "🔟 Extrayendo secuencias relevantes (≥1 mutación ≥10%)..."

for variant in "${variants[@]}"; do
    seqkit grep -f results/$variant/ids_muts10.txt \
        results/$variant/spike_${variant}.fasta \
        > results/$variant/spike_${variant}_muts10.fasta
done

echo "=========================================="
echo "1️⃣1️⃣ Limpiando secuencias (removiendo gaps)..."

for variant in "${variants[@]}"; do
    seqkit seq -u results/$variant/spike_${variant}_muts10.fasta \
        > results/$variant/spike_${variant}_muts10_clean.fasta
done

echo "=========================================="
echo "1️⃣2️⃣ Construyendo haplotipos reales (identidad 100%)..."

for variant in "${variants[@]}"; do
    python3 scripts/build_haplotypes.py \
        results/$variant/spike_${variant}_muts10_clean.fasta \
        results/$variant/${variant}
done

echo "=========================================="
echo "1️⃣3️⃣ Construyendo haplotipos dominantes (versión 2)..."

for variant in "${variants[@]}"; do
    python3 scripts/build_haplotypes_version2.py \
        results/$variant/spike_${variant}_muts10_clean.fasta \
        results/$variant/mutaciones_${variant}.tsv \
        results/$variant/${variant}
done

echo "=========================================="
echo "1️⃣3️⃣ Dividiendo haplotipos en FASTAs individuales..."

for variant in "${variants[@]}"; do
    python3 scripts/split_haplotypes_for_modeling.py \
        results/$variant/${variant}_haplotypes_dom.fasta \
        results/$variant/Seq_model \
        $variant
done

echo "=========================================="
echo "✅ ¡Listo para modelado 3D!"
echo "Archivos generados:"
echo "  • results/<variante>/<variante>_alineado.fasta"
echo "  • scripts/sitios_glyc_ref.txt"
echo "  • results/<variante>/mutaciones_<variante>.tsv"
echo "Usar estos archivos como entrada para modelado por homología y simulación de glicosilación."
echo "=========================================="

































