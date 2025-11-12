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
source venv/bin/activate  # En Windows: venv\Scripts\activate

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
python scripts/analisis_funcional.py -i data/genes_input.txt -o results/ --excel
```

**3. Visualizar los resultados del enriquecimiento**

Este comando leerá los CSV generados en el paso anterior y creará los gráficos de barras en `results/plots/`.

```bash
python scripts/plot_results.py -d results/ -o results/plots/ --top 15
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
python scripts/plot_propagation.py --dir results/ --outdir results/plots/ --top 20
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
├── data/                # Archivos de entrada (lista de genes, redes, etc.)
├── results/             # Resultados generados por el script
├── scripts/             # Código fuente del proyecto
├── docs/                # Documentación adicional (opcional)
├── README.md            # Este archivo
└── requirements.txt     # Dependencias del proyecto
```
