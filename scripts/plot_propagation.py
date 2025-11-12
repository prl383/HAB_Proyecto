#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
plot_propagation.py — Visualización de RWR y DIAMOnD
----------------------------------------------------
Genera gráficos de barras a partir de:
- rwr_scores.csv (columnas esperadas: node, score; score alto = más relevante)
- diamond_ranking.csv (columnas esperadas: node, score=pvalue, step_added; p bajo = más relevante)

Por defecto:
- RWR: muestra Top-N por 'score' (desc).
- DIAMOnD: convierte p a -log10(p) para que barras grandes = más significativo.

Uso (Windows, una línea):
  python scripts\\plot_propagation.py --dir results --outdir results\\plots --top 20
"""

from __future__ import annotations
import argparse
from pathlib import Path
import math
import pandas as pd
import matplotlib.pyplot as plt


def _safe_read_csv(p: Path) -> pd.DataFrame:
    if not p or not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(p)
    except Exception:
        # segundo intento por si tiene separador raro
        return pd.read_csv(p, sep=None, engine="python")


def _ensure_outdir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _shorten(text: str, max_len: int = 70) -> str:
    s = str(text)
    return s if len(s) <= max_len else s[: max_len - 1] + "…"


def _plot_barh(df: pd.DataFrame, term_col: str, score_col: str,
               title: str, outfile: Path, width: float, height: float, dpi: int):
    if df.empty:
        print(f"[INFO] No hay datos para '{title}', no se genera figura.")
        return
    plt.figure(figsize=(width, height))
    plt.barh(df[term_col], df[score_col])
    plt.xlabel(score_col)
    plt.title(title)
    plt.tight_layout()
    plt.gca().invert_yaxis()
    _ensure_outdir(outfile.parent)
    plt.savefig(outfile, dpi=dpi)
    plt.close()
    print(f"[OK] Guardado: {outfile}")


def _compute_neglog10(x):
    if pd.isna(x) or x is None:
        return None
    try:
        x = float(x)
    except Exception:
        return None
    if x <= 0:
        return None
    return -math.log10(x)


def plot_rwr(rwr_path: Path, outdir: Path, top: int, figsize, dpi: int):
    df = _safe_read_csv(rwr_path)
    if df.empty:
        print(f"[INFO] rwr_scores vacío o no encontrado: {rwr_path}")
        return
    # columnas esperadas: node, score
    if "node" not in df.columns or "score" not in df.columns:
        print(f"[WARN] rwr_scores no tiene columnas esperadas (node, score): {rwr_path}")
        return
    df = df.sort_values("score", ascending=False).head(top).copy()
    df["label"] = df["node"].map(_shorten)
    _plot_barh(df, "label", "score",
               title=f"RWR — Top {len(df)} nodos por score",
               outfile=outdir / "plot_rwr_top.png",
               width=figsize[0], height=figsize[1], dpi=dpi)


def plot_diamond(diamond_path: Path, outdir: Path, top: int, figsize, dpi: int, raw_p: bool):
    df = _safe_read_csv(diamond_path)
    if df.empty:
        print(f"[INFO] diamond_ranking vacío o no encontrado: {diamond_path}")
        return
    # columnas esperadas: node, score(=pvalue), step_added
    if "node" not in df.columns or "score" not in df.columns:
        print(f"[WARN] diamond_ranking no tiene columnas esperadas (node, score[, step_added]): {diamond_path}")
        return

    df = df.copy()
    if raw_p:
        # barras por p-value (pequeño = mejor) -> orden ascendente
        df = df.sort_values("score", ascending=True).head(top)
        score_col = "score"
        x_label = "p-value"
        title = f"DIAMOnD — Top {len(df)} nodos por p-value (menor es mejor)"
    else:
        # -log10(p) (grande = mejor)
        df["neglog10_p"] = df["score"].apply(_compute_neglog10)
        df = df.dropna(subset=["neglog10_p"])
        df = df.sort_values("neglog10_p", ascending=False).head(top)
        score_col = "neglog10_p"
        x_label = "-log10(p)"
        title = f"DIAMOnD — Top {len(df)} nodos por -log10(p)"

    # Etiquetas: si tenemos step_added, lo mostramos
    if "step_added" in df.columns:
        df["label"] = df.apply(
            lambda r: _shorten(f"{r['node']} (step {int(r['step_added'])})"),
            axis=1
        )
    else:
        df["label"] = df["node"].map(_shorten)

    # Renombrar columna de puntuación solo para el eje X
    df = df.rename(columns={score_col: x_label})

    _plot_barh(df, "label", x_label,
               title=title,
               outfile=outdir / ("plot_diamond_top_rawp.png" if raw_p else "plot_diamond_top.png"),
               width=figsize[0], height=figsize[1], dpi=dpi)


def build_parser():
    p = argparse.ArgumentParser(
        description="Visualiza resultados de propagación en red (RWR y DIAMOnD)."
    )
    # Fuentes
    p.add_argument("--dir", type=str, default=None,
                   help=r"Carpeta con resultados (autodetecta rwr_scores.csv y diamond_ranking.csv).")
    p.add_argument("--rwr", type=str, default=None, help="Ruta a rwr_scores.csv.")
    p.add_argument("--diamond", type=str, default=None, help="Ruta a diamond_ranking.csv.")
    # Salida
    p.add_argument("--outdir", type=str, required=True, help=r"Carpeta para PNGs (p. ej. results\plots).")
    # Opciones de plot
    p.add_argument("--top", type=int, default=20, help="Cantidad de términos/genes a mostrar (default: 20).")
    p.add_argument("--dpi", type=int, default=150, help="Resolución de las figuras (default: 150).")
    p.add_argument("--size", type=float, nargs=2, metavar=("WIDTH", "HEIGHT"),
                   default=(10.0, 6.0), help="Tamaño de figura en pulgadas (default: 10 6).")
    p.add_argument("--diamond-raw-p", action="store_true",
                   help="Para DIAMOnD, graficar p-value crudo (en vez de -log10(p)).")
    return p


def main():
    args = build_parser().parse_args()
    outdir = Path(args.outdir)
    _ensure_outdir(outdir)

    # Detectar archivos si se pasó --dir
    rwr_path = Path(args.rwr) if args.rwr else None
    diamond_path = Path(args.diamond) if args.diamond else None
    if args.dir:
        base = Path(args.dir)
        if rwr_path is None:
            cand = base / "rwr_scores.csv"
            if cand.exists():
                rwr_path = cand
        if diamond_path is None:
            cand = base / "diamond_ranking.csv"
            if cand.exists():
                diamond_path = cand

    if rwr_path:
        print(f"[INFO] Graficando RWR desde: {rwr_path}")
        plot_rwr(rwr_path, outdir, args.top, args.size, args.dpi)
    else:
        print("[INFO] No se proporcionó rwr_scores.csv (omite gráfica RWR).")

    if diamond_path:
        print(f"[INFO] Graficando DIAMOnD desde: {diamond_path}")
        plot_diamond(diamond_path, outdir, args.top, args.size, args.dpi, args.diamond_raw_p)
    else:
        print("[INFO] No se proporcionó diamond_ranking.csv (omite gráfica DIAMOnD).")

    print("[DONE] Figuras generadas.")


if __name__ == "__main__":
    main()
