# 🧬 Análisis estructural y evolutivo de glicanos y mutaciones en Spike de SARS-CoV-2

---
⚠️ Nota importante

En Example/ se incluye un mini-dataset (2 genomas) que muestra cómo luce la estructura de resultados final del pipeline.
Para asegurar reproducibilidad, el pipeline debe ejecutarse con Pixi, usando el pixi.toml del repositorio.
Este archivo fija las versiones exactas de todas las dependencias, evitando errores o comportamientos distintos entre entornos.

---

## 📌 Hipótesis

Las mutaciones de alta frecuencia en la proteína Spike del SARS-CoV-2 coevolucionan con los sitios de N-glicosilación, modificando el blindaje glicano, la accesibilidad antigénica y la estructura tridimensional de la proteína a lo largo de la evolución de las variantes.

---

## 🎯 Objetivo

Evaluar cómo las mutaciones frecuentes de Spike afectan la estructura y el entorno tridimensional de los sitios N-glicosilados mediante haplotipos representativos y modelado estructural.

## Flujo de trabajo

![Brutal ](https://github.com/user-attachments/assets/1657fd95-6936-440c-9b83-369de9868391)

---

## 🧪 Datos utilizados

- Secuencias de GISAID
- Variantes analizadas: **Lambda, Gamma, Delta y Ómicron**.

### Criterios de inclusión
✔ Sin frameshift en Spike  
✔ Sin aminoácidos ambiguos (“X”)  
✔ Con ≥1 mutación de frecuencia ≥10%

---

## 🛠️ Lenguajes y herramientas

### Lenguajes
- Bash  
- Python 3 (BioPython, pandas)  
- R (para visualización)

### Herramientas externas
- **Nextclade** — clasificación y QC  
- **SeqKit** — filtrado y extracción  
- **MAFFT** — alineamiento por variante  
- **Modeller** *(in progress)*  
- **GlycoSHIELD** *(futuro análisis glicanos)*

---

## ▶️ Cómo usar el pipeline

### 1️⃣ Ejecutar

```bash
./glycotestV2.sh -i genomas.fasta
```

## 🧱 Scripts incluidos

| Script | Función |
|--------|---------|
| **analisis_mut_aa.py** | Mutaciones AA + frecuencia + cambios fisicoquímicos |
| **crear_sitios_glyc_ref.py** | Detecta N-X-S/T en Spike Wuhan |
| **build_haplotypes.py** | Haplotipos reales |
| **build_haplotypes_version2.py** | Haplotipos dominantes (≥10%) |
| **split_haplotypes_for_modeling.py** | FASTA individuales para modelaje |

---

## 🧬 Estado del proyecto

### ✔ Análisis de secuencias  
### ✔ Construcción de haplotipos  
### ✔ Limpieza y filtrado  
### ✔ Generación de secuencias finales  
### 🔧 Modelado estructural *(in progress)*  
### 🔧 Dinámica de glicosilación *(in progress)*  
### 🔧 Proximidad mutación-glicano *(in progress)*  

---

## 📚 Contacto

**Diego Rivas Montani**  
Laboratorio de Epidemiología Molecular y Genética (CITBM)  



