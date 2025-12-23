# Análisis evolutivo de sitios de N glicosilacion y mutaciones en Spike de SARS-CoV-2

---
⚠️ Nota importante

En Example/ se incluye un mini-dataset (2 genomas) que muestra cómo luce la estructura de resultados final del pipeline.
Para asegurar reproducibilidad, el pipeline debe ejecutarse con Pixi, usando el pixi.toml del repositorio.
Este archivo fija las versiones exactas de algunas dependencias, evitando errores o comportamientos distintos entre entornos.

---

## Hipótesis

Las mutaciones de alta frecuencia en la proteína Spike coevolucionan con los sitios de N-glicosilación, mostrando patrones posicionales no aleatorios en términos de carga eléctrica y propiedades fisicoquímicas dentro de ventanas ±10 aminoácidos alrededor de los secuones N-X-S/T.

---

## Objetivo

Caracterizar la relación evolutiva entre mutaciones de alta frecuencia y sitios N-glicosilados en Spike, comparando la secuencia de referencia Wuhan con variantes posteriores, mediante análisis posicionales y fisicoquímicos.

---

---
## Objetivos específicos

* Detectar sitios N-X-S/T en la proteína Spike de referencia (Wuhan).
* Identificar mutaciones de alta frecuencia en variantes (Gamma, Delta, Lambda, Omicron).
* Comparar antes (Wuhan) vs después (variantes) en términos de:
  -carga eléctrica
  -polaridad
  -cambios fisicoquímicos
* Evaluar patrones posicionales de mutaciones dentro de ventanas ±10 aa alrededor de sitios glicosilados.
* Identificar posiciones relativas con comportamiento conservado entre variantes.
---

## Flujo de trabajo

![Brutal ](https://github.com/user-attachments/assets/1657fd95-6936-440c-9b83-369de9868391)

---

## Lenguajes y herramientas

### Lenguajes
- Bash (script base)
- Python 3 (análisis de mutaciones, haplotipos y secuencias)  
- R (visualización)

### Herramientas externas
- **Nextclade** — clasificación y QC  
- **SeqKit** — filtrado y extracción  (Modificar variantes en https://github.com/DiegoAlexRM/Proyecto-Gen_Evo/blob/main/glycotestV1.2.sh)
- **MAFFT** — alineamiento por variante  
- **Modeller** *(in progress)*  
- **GlycoSHIELD** *(futuro análisis glicanos)*

---

## Cómo usar el pipeline

### Ejecutar

```bash
./glycotestV2.sh -i genomas.fasta
```
* Para visualización de datos ejecutar Rcode_graphs_code.Rmd en Rstudio

## Scripts incluidos

| Script | Función |
|--------|---------|
| **analisis_mut_aa.py** | Mutaciones AA + frecuencia + cambios fisicoquímicos |
| **crear_sitios_glyc_ref.py** | Detecta N-X-S/T en Spike Wuhan |
| **analisis_nglyc.py** | conservación de sitios potenciales N-X-S/T |
| **correlate_mut_glyc.py** | Correlaciona mutaciones aminoacídicas con su proximidad a sitios N-X-S/T, cuantificando mutaciones dentro de ventanas definidas de 10aa |
| **build_haplotypes.py** | Haplotipos reales |
| **build_haplotypes_version2.py** | Haplotipos dominantes (≥10%) |
| **split_haplotypes_for_modeling.py** | FASTA individuales para modelaje |
| **Rcode_graphs_code.Rmd** | Visualización de datos |
---

## Estado del proyecto

### Análisis de secuencias  
### Construcción de haplotipos  
### Limpieza y filtrado  
### Generación de secuencias finales 
### Visualización de datos
### Proximidad mutación-glicano (análisis posicional)
### Modelado estructural *(in progress)*  
### Dinámica de glicosilación *(in progress)*  
### Proximidad mutación-glicano *(in progress)*  

---

## Contacto

**Diego Rivas Montani**  
Laboratorio de Epidemiología Molecular y Genética (CITBM)  



