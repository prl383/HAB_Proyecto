# Proyecto de Análisis Funcional y Propagación en Redes

Este proyecto contiene un pipeline bioinformático en Python para el análisis funcional y la priorización de genes basada en redes a partir de una lista de entrada (p. ej., genes diferencialmente expresados). El pipeline automatiza la consulta a bases de datos biológicas, ejecuta algoritmos de propagación en red y genera visualizaciones para una fácil interpretación de los resultados.

## Objetivo
Desarrollar un script en Python que realice análisis funcional con propagación en redes a partir de una lista de genes diferencialmente expresados.

## Características

*   **Análisis de Enriquecimiento Funcional:**
    *   Integración con **g:Profiler** (GO, KEGG, Reactome).
    *   Integración con **Enrichr** (GO Biological Process, KEGG Human).
    *   Integración con la base de datos **STRING** para enriquecimiento funcional.
*   **Propagación en Red (Priorización de Genes):**
    *   Implementación de **Random Walk with Restart (RWR)** para puntuar nodos según su proximidad a las semillas.
    *   Implementación del algoritmo **DIAMOnD** para la expansión iterativa del módulo de genes.
*   **Modular y Configurable:** Control total a través de argumentos de línea de comandos (CLI).
*   **Visualización de Resultados:** Generación automática de gráficos de barras en formato PNG para los resultados más significativos.
*   **Salidas Organizadas:** Exportación de resultados en formatos tabulares (`.csv`) y unificado (`.xlsx`).

## Flujo de Trabajo

El pipeline se divide en dos tareas principales que se pueden ejecutar de forma independiente:

1.  **Análisis Funcional:**
    *   **Entrada:** `genes_input.txt`
    *   **Proceso:** `analisis_funcional.py` consulta g:Profiler, Enrichr y STRING.
    *   **Salida (Datos):** `results/*.csv` y `results/resultados_enriquecimiento.xlsx`.
    *   **Visualización:** `plot_results.py` lee los CSV y genera los gráficos `results/plots/*.png`.

2.  **Propagación en Red:**
    *   **Entrada:** `genes_seed.txt` y un fichero de red PPI (p. ej., `string_network.tsv`).
    *   **Proceso:** `network_propagation.py` ejecuta RWR o DIAMOnD.
    *   **Salida (Datos):** `results/rwr_scores.csv` o `results/diamond_ranking.csv`.
    *   **Visualización:** `plot_propagation.py` lee los CSV de ranking y genera los gráficos `results/plots/*.png`.

## Instalación y Requisitos

### Prerrequisitos

*   Python 3.8 o superior.

### Instalación de Librerías

Se recomienda crear un entorno virtual para instalar las dependencias.

```bash
# Crear un entorno virtual (opcional pero recomendado)
python -m venv venv
source venv\Scripts\activate  # En Linux: venv/bin/activate

# Instalar las librerías necesarias
pip install pandas matplotlib requests gprofiler-official networkx openpyxl mygene
```

También puedes crear un fichero `requirements.txt` con el siguiente contenido y ejecutar `pip install -r requirements.txt`:

```
pandas
matplotlib
requests
gprofiler-official
networkx
openpyxl
mygene
```

## Uso

A continuación se muestran ejemplos de cómo ejecutar cada parte del pipeline desde la línea de comandos.

### Tarea 1: Análisis Funcional

**1. Preparar el fichero de genes**

Crea un fichero `data/genes_input.txt` con tus genes de interés, uno por línea o separados por comas.

```
# data/genes_input.txt
COX4I2
ND1
ATP6
# También se puede poner en una línea: COX4I2, ND1, ATP6
```

**2. Ejecutar el script de análisis**

Este comando ejecutará el análisis con g:Profiler, Enrichr y STRING, y guardará los resultados en la carpeta `results/`, incluyendo un Excel unificado.

```bash
python scripts\analisis_funcional.py -i data\genes_input.txt -o results\ --excel
```

**3. Visualizar los resultados del enriquecimiento**

Este comando leerá los CSV generados en el paso anterior y creará los gráficos de barras en `results/plots/`.

```bash
python scripts\plot_results.py -d results\ -o results\plots\ --top 15
```

### Tarea 2: Propagación en Red

**1. Preparar los ficheros de entrada**

*   **Semillas:** Un fichero de texto con los genes semilla, similar al anterior (`data/genes_seed.txt`).
*   **Red:** Un fichero con la red de interacciones. Por ejemplo, un fichero TSV descargado de STRING con las columnas `protein1_hugo`, `protein2_hugo`, `combined_score`.

**2. Ejecutar el algoritmo de propagación**

*   **Ejemplo con RWR:**
    ```bash
    python scripts\network_propagation.py --network data\string_network_filtered_hugo-400.tsv --format string --seeds data\genes_input.txt --algo rwr --outdir results --min-score 700
    ```

*   **Ejemplo con DIAMOnD:**
    ```bash
    python scripts\network_propagation.py --network data\string_network_filtered_hugo-400.tsv --format string --seeds data\genes_input.txt --algo diamond --outdir results --steps 20 --min-score 700
    ```

**3. Visualizar los resultados de la propagación**

Este comando buscará `rwr_scores.csv` y/o `diamond_ranking.csv` en la carpeta `results/` y generará los gráficos correspondientes.

```bash
python scripts\plot_propagation.py --dir results\ --outdir results\plots\ --top 20
```

## Descripción de los Scripts

*   `analisis_funcional.py`: Script principal para el análisis de enriquecimiento. Contacta las APIs de g:Profiler, Enrichr y STRING.
*   `plot_results.py`: Visualiza los resultados del enriquecimiento. Genera gráficos de barras con los términos más significativos.
*   `network_propagation.py`: Script principal para la propagación en red. Implementa los algoritmos RWR y DIAMOnD.
*   `plot_propagation.py`: Visualiza los resultados de la propagación. Genera gráficos de barras con los genes mejor clasificados.
*   `example_gene_conversion.py`: Script de utilidad para convertir identificadores de genes (p. ej., de Símbolo HUGO a UniProt ID) usando `MyGene.info`.

## Salidas Esperadas

Después de ejecutar el pipeline completo, la carpeta `results/` contendrá:

#### Ficheros de Datos (`.csv`, `.xlsx`)

*   `gprofiler_enrichment.csv`: Términos enriquecidos de GO, KEGG y Reactome (vía g:Profiler).
*   `enrichr_GO_Biological_Process_2023.csv`: Resultados de enriquecimiento para GO BP (vía Enrichr).
*   `enrichr_KEGG_2021_Human.csv`: Resultados de enriquecimiento para KEGG (vía Enrichr).
*   `string_enrichment.csv`: Resultados de enriquecimiento funcionales de STRING.
*   `resultados_enriquecimiento.xlsx`: Fichero Excel con todos los resultados de enriquecimiento en pestañas separadas.
*   `rwr_scores.csv`: Ranking de todos los nodos de la red según su score RWR.
*   `diamond_ranking.csv`: Ranking de los nodos añadidos por el algoritmo DIAMOnD.

#### Figuras (`.png`)

La subcarpeta `results/plots/` contendrá:

*   `plot_gprofiler_top.png`: Top términos de g:Profiler.
*   `plot_enrichr_..._top.png`: Top términos de Enrichr (GO y KEGG).
*   `plot_string_top.png`: Top términos de STRING.
*   `plot_rwr_top.png`: Top genes según el score de RWR.
*   `plot_diamond_top.png`: Top genes según la significancia en DIAMOnD (-log10(p-valor)).

## Estructura del repositorio

```
├── data/                                                  # Archivos de entrada (lista de genes, redes, etc.)
│   └── genes\_input.txt                                   # Genes semilla
│   └── string\_network\_filtered\_hugo-400.tsv            # Red de interacciones proteicas filtrada de STRING
│   └── uniprot\_output.tsv                                # Lista de genes con IDs UniProt
│   └── Allcontrasts\_GLM-Treat\_P-0.1\_FC-1.25\_2025-10-14\_16.57.27.tsv
├── results/                                               # Resultados generados por el script
│   └── plots/
│   └── diamond\_ranking.csv
│   └── enrichr\_GO\_Biological\_Process\_2023.csv
│   └── enrichr\_KEGG\_2021\_Human.csv
│   └── gprofiler\_enrichment.csv
│   └── rwr\_scores.csv
│   └── string\_enrichment.csv
│   └── uniprot\_output.tsv
├── scripts/                                               # Código fuente del proyecto
│   └── analisis\_funcional.py
│   └── example\_gene_conversion.py
│   └── network\_propagation.py
│   └── plot\_propagation.py
│   └── plot\_results.py
├── LICENSE
├── README.md
└── requirements.txt

```
---

## Referencias

### Herramientas Científicas y Bases de Datos

*   **g:Profiler:** Para el análisis de enriquecimiento funcional.
    > Kolberg, L., Raudvere, U., Kuzmin, I., Vilo, J., & Peterson, H. (2023). g:Profiler—a web server for functional enrichment analysis and conversions of gene lists (2023 update). *Nucleic Acids Research*, 51(W1), W277-W282.

*   **Enrichr:** Para el análisis de enriquecimiento utilizando las librerías del Ma'ayan Lab.
    > Xie, Z., Bailey, A., Kuleshov, M. V., Clarke, D. J., Evangelista, J. E., Jenkins, S. L., ... & Ma'ayan, A. (2021). The Gene Set Library knowledge base for gene set enrichment analysis. *Bioinformatics*, 37(19), 3381-3383.

*   **STRING Database:** Para la red de interacción proteína-proteína y su servicio de enriquecimiento. La versión 11.5 fue utilizada en este pipeline.
    > Szklarczyk, D., Gable, A. L., Nastou, K. C., Lyon, D., Kirsch, R., Pyl, P. T., ... & Jensen, L. J. (2021). The STRING database in 2021: customizable protein–protein networks, and functional characterization of user-uploaded gene/measurement sets. *Nucleic Acids Research*, 49(D1), D605-D612.

*   **MyGene.info:** Para la conversión de identificadores de genes.
    > Xin, J., Mark, A., Afrasiabi, C., Tsueng, G., Juchli, M., Gopalakrishnan, K., ... & Wu, C. (2016). High-performance python web services for querying gene and variant annotation. *bioRxiv*, 095143.

### Algoritmos de Propagación en Red

*   **Random Walk with Restart (RWR):** El enfoque de "difusión" implementado se inspira en métodos como el utilizado en GUILD.
    > Guney, E., Menche, J., Vidal, M., & Barábasi, A. L. (2016). Network-based in silico drug efficacy screening. *Nature communications*, 7(1), 1-13.

*   **DIAMOnD (Disease Module Detection):** Para la identificación iterativa de proteínas relevantes basada en significancia conectiva.
    > Ghiassian, S. D., Menche, J., & Barabási, A. L. (2015). A DIseAse MOdule Detection (DIAMOnD) algorithm derived from a systematic analysis of connectivity patterns of disease proteins in the human interactome. *PLoS computational biology*, 11(4), e1004120.

### Software y Librerías

El desarrollo de este pipeline no habría sido posible sin el ecosistema de computación científica de Python:

*   **Pandas:** Para la manipulación y análisis de datos.
*   **NetworkX:** Para la creación, manipulación y estudio de redes complejas.
*   **Matplotlib:** Para la generación de visualizaciones y gráficos.
*   **gprofiler-official, requests, openpyxl:** Para la interacción con APIs y la manipulación de ficheros.

### Desarrollo y Asistencia con IA

Parte del código, la estructura de los scripts y la documentación de este proyecto se desarrollaron con la asistencia de Inteligencia Artificial generativa. Se consultaron Modelos de Lenguaje Grandes (LLMs) como **ChatGPT** de OpenAI y **Gemini** de Google para:

*   Generar código repetitivo (`boilerplate`) para funciones y clases.
*   Refactorizar y optimizar bloques de código para mejorar su legibilidad y eficiencia.
*   Sugerir estrategias para la depuración de errores.
*   Redactar y formatear la documentación, incluyendo este fichero `README.md`.

Los autores del proyecto ha revisado, validado y adaptado todo el contenido generado por la IA, y asume la total responsabilidad por la corrección y funcionalidad del código final.

### Autores

Este trabajo lo ha realizado:

1. Anabel Yu Flores Moral
2. Achraf Ousti El Moussati
3. Patricia Rodríguez Lidueña
