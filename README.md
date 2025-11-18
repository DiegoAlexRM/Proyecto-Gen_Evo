📘 Proyecto Gen_Evo – Pipeline Spike SARS-CoV-2

Pipeline para procesar secuencias Spike, identificar mutaciones relevantes, generar haplotipos reales y preparar secuencias listas para modelado estructural 3D.

🔧 Requisitos

nextclade

seqkit

mafft

python3 + biopython

Idealmente usar Pixi o Conda.

🚀 Uso
bash scripts/glycotestV2.sh -i secuencias.fasta

📂 Qué hace el pipeline

Corre Nextclade

Filtra Spike sin frameshift

Alinea por variante (MAFFT + referencia Wuhan)

Detecta mutaciones AA

Filtra secuencias con mutaciones ≥10%

Limpia gaps

Construye haplotipos (100% identidad)

Exporta FASTA listo para modelado

📁 Salidas principales

Para cada variante:

mutaciones_<variante>.tsv

ids_muts10.txt

spike_<variante>_muts10_clean.fasta

haplotypes_<variante>.tsv

hap_<n>.fasta (secuencias representativas)

🎯 Objetivo

Obtener haplotipos representativos basados en mutaciones frecuentes, adecuados para:

modelado por homología

análisis 3D de proximidad a glicanos

comparación estructural entre variantes

