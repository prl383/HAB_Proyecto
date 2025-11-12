#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
plot_results.py — Visualización de resultados de enriquecimiento
----------------------------------------------------------------
Lee los CSV generados por `tu_script.py` (g:Profiler, Enrichr y STRING) y produce
gráficos de barras (horizontal) con los términos/top pathways más significativos.

Métricas usadas por fuente:
- g:Profiler: usa p_adj -> score = -log10(p_adj)
- Enrichr: usa p_adj si existe, si no p_value -> score = -log10(p)
- STRING: usa fdr -> score = -log10(fdr)

Salidas:
- PNG por cada fuente seleccionada (en la carpeta indicada con -o/--outdir).

Uso (Windows):
  python scripts\\plot_results.py -d results\\ -o results\\plots\\ --top 15

  Análisis funcional con salida combinada en Excel (.xlsx)
  python scripts\analisis_funcional.py -i data\genes_input.txt -o results\ --excel

Requisitos:
  - pandas
  - matplotlib
"""

from __future__ import annotations
import argparse
from pathlib import Path
from typing import Optional, Tuple
import math

import pandas as pd
import matplotlib.pyplot as plt


# ------------------------------- utilidades ----------------------------------

def _safe_read_csv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path)


def _compute_score_from_p(p: float) -> float:
    if p is None or pd.isna(p) or p <= 0:
        return float("nan")
    return -math.log10(p)


def _prep_terms_for_plot(df: pd.DataFrame, term_col: str, score_col: str,
                         top: int) -> pd.DataFrame:
    """Ordena por score descendente y acorta nombres muy largos."""
    if df.empty:
        return df
    df = df.copy()
    df = df.dropna(subset=[score_col])
    df = df.sort_values(by=score_col, ascending=False).head(top)

    # Acorta nombres muy largos para el eje Y (opcional)
    def _shorten(name: str, max_len: int = 80) -> str:
        s = str(name)
        return s if len(s) <= max_len else s[:max_len - 1] + "…"

    df["term_for_plot"] = df[term_col].apply(_shorten)
    return df


def _plot_barh(df: pd.DataFrame, title: str, term_col: str, score_col: str,
               outfile: Path, dpi: int = 150, width: float = 10, height: float = 6):
    if df.empty:
        print(f"[WARN] Nada que plotear para {title}.")
        return

    # Configurar figura (sin estilos ni colores específicos)
    plt.figure(figsize=(width, height))
    plt.barh(df["term_for_plot"], df[score_col])
    plt.xlabel("-log10(p)")  # etiqueta genérica (aplica también a FDR)
    plt.title(title)
    plt.tight_layout()
    # Hacer que el término más significativo aparezca arriba
    plt.gca().invert_yaxis()
    outfile.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(outfile, dpi=dpi)
    plt.close()
    print(f"[OK] Guardado: {outfile}")


# ------------------------------- cargadores ----------------------------------

def load_gprofiler(results_dir: Path, top: int) -> Tuple[pd.DataFrame, str, str]:
    path = results_dir / "gprofiler_enrichment.csv"
    df = _safe_read_csv(path)
    if df.empty:
        return df, "", ""
    # Espera columnas: term_name, p_adj, source_db, etc. (según tu_script.py)
    if "p_adj" not in df.columns:
        return pd.DataFrame(), "", ""
    df["score"] = df["p_adj"].apply(_compute_score_from_p)
    df = _prep_terms_for_plot(df, "term_name", "score", top)
    title = "g:Profiler — Top términos (por -log10(p_adj))"
    fname = "plot_gprofiler_top.png"
    return df, title, fname


def load_enrichr(results_dir: Path, lib: str, top: int) -> Tuple[pd.DataFrame, str, str]:
    # lib esperado: "GO_Biological_Process_2023" o "KEGG_2021_Human"
    path = results_dir / f"enrichr_{lib}.csv"
    df = _safe_read_csv(path)
    if df.empty:
        return df, "", ""
    # p_adj si existe; si no, p_value
    p_col = "p_adj" if "p_adj" in df.columns and df["p_adj"].notna().any() else "p_value"
    if p_col not in df.columns:
        return pd.DataFrame(), "", ""
    df["score"] = df[p_col].apply(_compute_score_from_p)
    df = _prep_terms_for_plot(df, "term_name", "score", top)
    pretty = lib.replace("_", " ")
    title = f"Enrichr ({pretty}) — Top términos (por -log10({p_col}))"
    fname = f"plot_enrichr_{lib}_top.png"
    return df, title, fname


def load_string(results_dir: Path, top: int) -> Tuple[pd.DataFrame, str, str]:
    path = results_dir / "string_enrichment.csv"
    df = _safe_read_csv(path)
    if df.empty:
        return df, "", ""
    if "fdr" not in df.columns:
        return pd.DataFrame(), "", ""
    df["score"] = df["fdr"].apply(_compute_score_from_p)
    df = _prep_terms_for_plot(df, "term_name", "score", top)
    title = "STRING — Top términos (por -log10(FDR))"
    fname = "plot_string_top.png"
    return df, title, fname


# --------------------------------- CLI ---------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Visualiza resultados de enriquecimiento en gráficos de barras."
    )
    p.add_argument(
        "-d", "--results-dir",
        type=str,
        required=True,
        help=r"Carpeta con CSVs de salida (p. ej. results\ )"
    )
    p.add_argument(
        "-o", "--outdir",
        type=str,
        required=True,
        help=r"Carpeta donde guardar los PNG (p. ej. results\plots\ )"
    )
    p.add_argument(
        "--top",
        type=int,
        default=15,
        help="Cuántos términos mostrar (default: 15)."
    )
    p.add_argument(
        "--which",
        type=str,
        nargs="+",
        choices=["gprofiler", "enrichr_go_bp", "enrichr_kegg", "string", "all"],
        default=["all"],
        help="Qué fuentes graficar. Por defecto: all."
    )
    p.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Resolución de las figuras (default: 150)."
    )
    p.add_argument(
        "--size",
        type=float,
        nargs=2,
        metavar=("WIDTH", "HEIGHT"),
        default=(10.0, 6.0),
        help="Tamaño de figura en pulgadas, ej. --size 11 7 (default: 10 6)."
    )
    return p


def main():
    args = build_parser().parse_args()
    results_dir = Path(args.results_dir)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    width, height = args.size

    targets = set(args.which)
    if "all" in targets:
        targets = {"gprofiler", "enrichr_go_bp", "enrichr_kegg", "string"}

    # g:Profiler
    if "gprofiler" in targets:
        df, title, fname = load_gprofiler(results_dir, args.top)
        if not df.empty:
            _plot_barh(df, title, "term_name", "score", out_dir / fname, dpi=args.dpi,
                       width=width, height=height)
        else:
            print("[INFO] g:Profiler: sin datos para graficar (¿CSV vacío o no existe?).")

    # Enrichr GO BP
    if "enrichr_go_bp" in targets:
        df, title, fname = load_enrichr(results_dir, "GO_Biological_Process_2023", args.top)
        if not df.empty:
            _plot_barh(df, title, "term_name", "score", out_dir / fname, dpi=args.dpi,
                       width=width, height=height)
        else:
            print("[INFO] Enrichr GO BP: sin datos para graficar.")

    # Enrichr KEGG
    if "enrichr_kegg" in targets:
        df, title, fname = load_enrichr(results_dir, "KEGG_2021_Human", args.top)
        if not df.empty:
            _plot_barh(df, title, "term_name", "score", out_dir / fname, dpi=args.dpi,
                       width=width, height=height)
        else:
            print("[INFO] Enrichr KEGG: sin datos para graficar.")

    # STRING
    if "string" in targets:
        df, title, fname = load_string(results_dir, args.top)
        if not df.empty:
            _plot_barh(df, title, "term_name", "score", out_dir / fname, dpi=args.dpi,
                       width=width, height=height)
        else:
            print("[INFO] STRING: sin datos para graficar.")

    print("[DONE] Gráficos generados.")


if __name__ == "__main__":
    main()
