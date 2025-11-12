#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
analisis_funcional.py — Análisis funcional (Tarea 1)
---------------------------------------------------
Implementa un análisis funcional reproducible para los genes de entrada 
(por defecto: COX4I2, ND1, ATP6 leídos desde data/genes_input.txt). 
El script acepta rutas de entrada/salida por CLI (compatible con Windows) 
y exporta resultados de:

1) g:Profiler (GO:BP, GO:CC, GO:MF, Reactome, KEGG)
   - Servicio: https://biit.cs.ut.ee/gprofiler/gost
   - Bases de datos: Gene Ontology, Reactome, KEGG (vía g:Profiler)
   - Ajuste de múltiples contrastes: g:SCS (propio de g:Profiler)
   - Filtrado por p ≤ 0.05.

2) Enrichr (vía API de Ma’ayan Lab)
   - Librerías usadas: GO_Biological_Process_2023, KEGG_2021_Human
   - Devuelve p-value y combined score. Filtrado por p ≤ 0.05.

3) STRING (v11.5) — enriquecimiento funcional / PPIs
   - Endpoint: https://version-11-5.string-db.org/api
   - Reporta categorías (incl. “Process”) y FDR.
   - Filtrado por FDR < 0.05.

Salida:
- CSVs en la carpeta results/ (creada automáticamente si no existe).
- Opcionalmente, XLSX unificado con el flag --excel.

Uso (Windows):
  python scripts\analisis_funcional.py -i data\genes_input.txt -o results\

Requisitos: ver requirements.txt

"""

from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
from typing import List, Dict, Any, Optional

import pandas as pd
import requests
from gprofiler import GProfiler


# ------------------------------- utilidades IO -------------------------------

def read_gene_list(input_path: Path) -> List[str]:
    """Lee símbolos génicos desde una línea (separados por comas) o varias líneas.
    - Ignora vacíos y comentarios (#).
    - Normaliza a mayúsculas.
    """
    text = input_path.read_text(encoding="utf-8")
    # Soporte para ambos formatos: por líneas o separados por comas
    # 1) quita comentarios por línea
    cleaned_lines = []
    for raw in text.splitlines():
        raw = raw.strip()
        if not raw or raw.startswith("#"):
            continue
        cleaned_lines.append(raw)
    # 2) une y separa por coma o espacios
    joined = ",".join(cleaned_lines)
    parts = [p.strip().upper() for p in joined.replace("\t", " ").split(",")]
    genes = [p for p in parts if p]  # quita vacíos
    if not genes:
        raise ValueError(f"No se encontraron genes en {input_path}.")
    return genes



def ensure_results_dir(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)


def safe_write_csv(df: pd.DataFrame, path: Path) -> None:
    if df.empty:
        # Crear CSV con encabezado mínimo para evitar errores al abrir
        path.write_text("", encoding="utf-8")
        return
    df.to_csv(path, index=False, quoting=csv.QUOTE_MINIMAL)


# ------------------------------- g:Profiler ----------------------------------

def run_gprofiler(
    genes: List[str],
    organism: str = "hsapiens",
    sources: Optional[List[str]] = None,
    p_threshold: float = 0.05
) -> pd.DataFrame:
    """
    Ejecuta enriquecimiento con g:Profiler y devuelve un DataFrame filtrado por p≤p_threshold.

    Notas:
      - g:Profiler aplica por defecto corrección de múltiples contrastes (g:SCS).
      - Campos típicos: source, native (ID del término), name, p_value, intersection_size, term_size, etc.
    """
    if sources is None:
        # GO:BP/MF/CC + Reactome + KEGG
        sources = ["GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG"]

    gp = GProfiler(return_dataframe=True)
    res = gp.profile(
        organism=organism,
        query=genes,
        sources=sources,
        user_threshold=1.0  # recuperamos todo y filtramos localmente
    )

    if res is None or len(res) == 0:
        return pd.DataFrame()

    # Filtrado de significancia
    res = res.loc[res["p_value"] <= p_threshold].copy()

    # Renombrados amigables
    rename_map = {
        "native": "term_id",
        "name": "term_name",
        "p_value": "p_adj",
        "source": "source_db",
        "intersection_size": "overlap",
        "term_size": "term_size",
        "query_size": "query_size",
    }
    for k, v in rename_map.items():
        if k in res.columns:
            res.rename(columns={k: v}, inplace=True)

    # Orden por significancia
    if "p_adj" in res.columns:
        res.sort_values(by="p_adj", inplace=True)

    return res


# -------------------------------- Enrichr ------------------------------------

ENRICHR_ADD_URL = "https://maayanlab.cloud/Enrichr/addList"
ENRICHR_ENRICH_URL = "https://maayanlab.cloud/Enrichr/enrich"
ENRICHR_LIBS = ["GO_Biological_Process_2023", "KEGG_2021_Human"]


def run_enrichr(genes: List[str], libraries: Optional[List[str]] = None) -> Dict[str, pd.DataFrame]:
    """
    Ejecuta enriquecimiento con Enrichr para las librerías indicadas.
    Devuelve diccionario {library: DataFrame (p<=0.05)} con columnas principales.
    """
    libs = libraries or ENRICHR_LIBS

    # Primero registramos la lista
    payload = {
        "list": (None, "\n".join(genes)),
        "description": (None, "Functional analysis gene list"),
    }
    r = requests.post(ENRICHR_ADD_URL, files=payload, timeout=60)
    r.raise_for_status()
    user_list_id = r.json()["userListId"]

    results = {}
    for lib in libs:
        params = {"userListId": user_list_id, "backgroundType": lib}
        resp = requests.get(ENRICHR_ENRICH_URL, params=params, timeout=120)
        resp.raise_for_status()
        data = resp.json()

        if lib not in data:
            results[lib] = pd.DataFrame()
            continue

        # Enrichr retorna listas; mapeamos a un DataFrame
        # Índices (según documentación): 0:Rank, 1:Term, 2:P-value, 3:Z-score, 4:Combined score, 5:Overlapping genes, 6:Adjusted p-value, etc.
        rows = []
        for row in data[lib]:
            rows.append(
                {
                    "term_name": row[1],
                    "p_value": float(row[2]),
                    "z_score": row[3],
                    "combined_score": row[4],
                    "overlapping_genes": row[5],
                    "p_adj": float(row[6]) if len(row) > 6 else None,
                    "lib": lib,
                }
            )
        df = pd.DataFrame(rows)
        # Filtramos por p<=0.05 (preferimos p ajustado si está disponible)
        if "p_adj" in df.columns and df["p_adj"].notna().any():
            df = df.loc[df["p_adj"] <= 0.05].copy()
            df.sort_values(by="p_adj", inplace=True)
        else:
            df = df.loc[df["p_value"] <= 0.05].copy()
            df.sort_values(by="p_value", inplace=True)
        results[lib] = df.reset_index(drop=True)

    return results


# -------------------------------- STRING -------------------------------------

STRING_API = "https://version-11-5.string-db.org/api"
STRING_FORMAT = "json"
STRING_METHOD = "enrichment"


def run_string_enrichment(genes: List[str], species: int = 9606) -> pd.DataFrame:
    """
    Llama al endpoint de enriquecimiento de STRING (v11.5) y devuelve un DataFrame.
    Filtra por FDR<0.05 y ordena por FDR.
    """
    url = "/".join([STRING_API, STRING_FORMAT, STRING_METHOD])
    params = {
        "identifiers": "%0d".join(genes),  # %0d = CRLF
        "species": species,
        "caller_identity": "analisis_funcional_cli",
    }
    resp = requests.post(url, data=params, timeout=120)
    resp.raise_for_status()
    data = json.loads(resp.text)

    if not isinstance(data, list):
        return pd.DataFrame()

    rows = []
    for row in data:
        # Campos típicos: term, preferredNames, fdr, category, description
        rows.append(
            {
                "term_id": row.get("term"),
                "term_name": row.get("description"),
                "preferred_names": ",".join(row.get("preferredNames", [])),
                "category": row.get("category"),
                "fdr": float(row.get("fdr")) if row.get("fdr") is not None else None,
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    # Filtramos por significancia
    if "fdr" in df.columns:
        df = df.loc[df["fdr"].notna() & (df["fdr"] < 0.05)].copy()
        df.sort_values(by="fdr", inplace=True)
    return df.reset_index(drop=True)


# --------------------------------- CLI main ----------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Análisis funcional de genes (GO/KEGG/Reactome, Enrichr y STRING)."
    )
    p.add_argument(
        "-i", "--input",
        type=str,
        required=True,
        help=r"Ruta al archivo de entrada (p. ej. data\genes_input.txt)"
    )
    p.add_argument(
        "-o", "--outdir",
        type=str,
        required=True,
        help=r"Carpeta de salida (p. ej. results\ )"
    )
    p.add_argument(
        "--organism",
        type=str,
        default="hsapiens",
        help="Organismo para g:Profiler (por defecto: hsapiens)."
    )
    p.add_argument(
        "--no-enrichr",
        action="store_true",
        help="Desactiva el análisis adicional con Enrichr."
    )
    p.add_argument(
        "--no-string",
        action="store_true",
        help="Desactiva el enriquecimiento STRING."
    )
    p.add_argument(
        "--excel",
        action="store_true",
        help="Además de CSV, escribe un Excel unificado (resultados.xlsx). Requiere openpyxl."
    )
    return p


def main():
    parser = build_parser()
    args = parser.parse_args()

    input_path = Path(args.input)
    out_dir = Path(args.outdir)
    ensure_results_dir(out_dir)

    # 1) Leer genes
    genes = read_gene_list(input_path)
    print(f"[INFO] Genes ({len(genes)}): {', '.join(genes)}")

    # 2) g:Profiler
    print("[INFO] Ejecutando g:Profiler (GO/Reactome/KEGG)…")
    gprof_df = run_gprofiler(genes, organism=args.organism)
    gprof_csv = out_dir / "gprofiler_enrichment.csv"
    safe_write_csv(gprof_df, gprof_csv)
    print(f"[OK] g:Profiler -> {gprof_csv} ({len(gprof_df)} términos significativos)")

    # 3) Enrichr (opcional)
    enrichr_results: Dict[str, pd.DataFrame] = {}
    if not args.no_enrichr:
        print("[INFO] Ejecutando Enrichr (GO BP y KEGG)…")
        try:
            enrichr_results = run_enrichr(genes)
        except Exception as e:
            print(f"[WARN] Enrichr falló: {e}")
        else:
            for lib, df in enrichr_results.items():
                path = out_dir / f"enrichr_{lib}.csv"
                safe_write_csv(df, path)
                print(f"[OK] Enrichr {lib} -> {path} ({len(df)} términos significativos)")

    # 4) STRING (opcional)
    string_df = pd.DataFrame()
    if not args.no_string:
        print("[INFO] Ejecutando STRING (enrichment)…")
        try:
            string_df = run_string_enrichment(genes)
        except Exception as e:
            print(f"[WARN] STRING falló: {e}")
        else:
            string_csv = out_dir / "string_enrichment.csv"
            safe_write_csv(string_df, string_csv)
            print(f"[OK] STRING -> {string_csv} ({len(string_df)} términos significativos)")

    # 5) Excel unificado (opcional)
    if args.excel:
        try:
            with pd.ExcelWriter(out_dir / "resultados_enriquecimiento.xlsx") as xw:
                if not gprof_df.empty:
                    gprof_df.to_excel(xw, sheet_name="gProfiler", index=False)
                for lib, df in enrichr_results.items():
                    if not df.empty:
                        df.to_excel(xw, sheet_name=f"Enrichr_{lib[:28]}", index=False)
                if not string_df.empty:
                    string_df.to_excel(xw, sheet_name="STRING", index=False)
            print(f"[OK] Excel -> {out_dir / 'resultados_enriquecimiento.xlsx'}")
        except Exception as e:
            print(f"[WARN] No se pudo escribir Excel (instala openpyxl): {e}")

    print("[DONE] Análisis funcional completado.")


if __name__ == "__main__":
    main()
