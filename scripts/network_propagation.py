#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
network_propagation.py — Tarea 2: Network Propagation
-----------------------------------------------------
Ejecuta propagación en red usando:
  1) RWR (Random Walk with Restart) — enfoque GUILD-like (NetRank/difusión).
  2) DIAMOnD — selección iterativa basada en p-valor hipergeométrico.

Redes soportadas (formato de entrada):
  - STRING (TSV con cabeceras: protein1_hugo, protein2_hugo, combined_score)
  - GUILD  (texto: 'u w v' separados por espacio; w = peso numérico)
  - DIAMOnD (texto: 'u,v' separados por coma; sin peso)

Semillas:
  - Archivo de texto con símbolos/IDs separados por comas en una sola línea
    o una por línea. Ej.: "ENO1, PGK1, HK2".

Salidas:
  - Para RWR: CSV con score de cada nodo (probabilidades estacionarias) → results\rwr_scores.csv
  - Para DIAMOnD: CSV con ranking iterativo (nodo, score, step_added) → results\diamond_ranking.csv

CLI (Windows, una línea):
  # RWR sobre STRING (HUGO), filtrando interacciones fuertes (min-score 700)
  python scripts\network_propagation.py --network data\string_network_filtered_hugo-400.tsv --format string --seeds data\genes_seed.txt --algo rwr --outdir results --min-score 700

  # DIAMOnD sobre STRING (HUGO), 20 pasos y min-score 700
  python scripts\network_propagation.py --network data\string_network_filtered_hugo-400.tsv --format string --seeds data\genes_seed.txt --algo diamond --outdir results --steps 20 --min-score 700

Requisitos: ver requirements.txt (pandas, networkx).
Notas:
- El parámetro --min-score solo aplica a redes STRING (usa la columna combined_score).
- Asegúrate de que el espacio de IDs de las semillas coincide con el de la red cargada.

"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import pandas as pd
import networkx as nx


# =============================== Utilidades IO ===============================

def read_seeds(path: Path) -> List[str]:
    """
    Lee genes semilla desde un archivo con:
      - una sola línea separada por comas, o
      - varias líneas (una por gen).
    Ignora líneas vacías y comentarios (#). Normaliza espacios.
    """
    text = path.read_text(encoding="utf-8")
    lines: List[str] = []
    for raw in text.splitlines():
        raw = raw.strip()
        if not raw or raw.startswith("#"):
            continue
        lines.append(raw)
    if not lines:
        raise ValueError(f"No se encontraron semillas en {path}")
    # Unificar y separar por coma
    joined = ",".join(lines)
    seeds = [p.strip() for p in joined.split(",") if p.strip()]
    # Mantener mayúsculas (para HUGO suele ser útil). No forzamos uppercase
    # porque hay redes con IDs numéricos/mixtos y no queremos transformarlos.
    return seeds


# ============================ Carga de redes (G) =============================

def load_string_network(path: Path, min_score: int = 0, undirected: bool = True) -> nx.Graph:
    """
    Carga red STRING filtrada (HUGO), TSV con cabeceras:
      protein1_hugo, protein2_hugo, combined_score
    Aplica umbral min_score si se indica (por ejemplo, 400, 700, 800).
    """
    df = pd.read_csv(path, sep="\t", header=0, dtype={"protein1_hugo": str,
                                                      "protein2_hugo": str,
                                                      "combined_score": int})
    if "combined_score" in df.columns and min_score > 0:
        df = df.loc[df["combined_score"] >= min_score].copy()

    edges = list(
        zip(df["protein1_hugo"].astype(str),
            df["protein2_hugo"].astype(str),
            df["combined_score"].astype(float))
    )

    G = nx.Graph() if undirected else nx.DiGraph()
    for u, v, w in edges:
        # STRING suele considerarse no-direccional
        G.add_edge(u, v, weight=w)
    return G


def load_guild_network(path: Path, undirected: bool = True) -> nx.Graph:
    """
    Carga red GUILD: líneas 'u w v' separadas por espacio.
    """
    G = nx.Graph() if undirected else nx.DiGraph()
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            u, w, v = parts[0], parts[1], parts[2]
            try:
                weight = float(w)
            except ValueError:
                # si el segundo campo no es numérico, asumimos peso=1 y tratamos como 'u v x'
                weight = 1.0
                v = parts[-1]
            G.add_edge(u, v, weight=weight)
    return G


def load_diamond_network(path: Path, undirected: bool = True) -> nx.Graph:
    """
    Carga red DIAMOnD: líneas 'u,v' (sin peso).
    """
    G = nx.Graph() if undirected else nx.DiGraph()
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "," not in line:
                continue
            u, v = [tok.strip() for tok in line.split(",", 1)]
            G.add_edge(u, v, weight=1.0)
    return G


def build_graph(network_path: Path, fmt: str,
                min_score: int = 0,
                undirected: bool = True) -> nx.Graph:
    fmt = fmt.lower()
    if fmt == "string":
        return load_string_network(network_path, min_score=min_score, undirected=undirected)
    if fmt == "guild":
        return load_guild_network(network_path, undirected=undirected)
    if fmt == "diamond":
        return load_diamond_network(network_path, undirected=undirected)
    raise ValueError(f"Formato de red no soportado: {fmt}")


# ================================ Algoritmos =================================
# --------------------------- 1) Random Walk w/ Restart -----------------------

def rwr_scores(G: nx.Graph,
               seeds: Iterable[str],
               alpha: float = 0.85,
               tol: float = 1e-9,
               max_iter: int = 100) -> Dict[str, float]:
    """
    Random Walk with Restart (RWR):
      r_{t+1} = (1-alpha) * P^T * r_t + alpha * r0
    donde P es la matriz de transición (normalizada por grado saliente).
    r0 es la distribución inicial (uniforme sobre semillas).

    Devuelve dict {nodo: score} (distribución estacionaria aproximada).
    """
    nodes = list(G.nodes())
    if not nodes:
        return {}

    # Asegurar que todas las semillas estén en la red
    seeds = [s for s in seeds if s in G]
    if not seeds:
        raise ValueError("Ninguna semilla está presente en la red. "
                         "Asegúrate de que el ID de semillas coincide con el ID de nodos.")

    # Índices
    idx = {n: i for i, n in enumerate(nodes)}

    # r0: uniforme sobre semillas
    r = [0.0] * len(nodes)
    for s in seeds:
        r[idx[s]] = 1.0 / len(seeds)

    # Precomputar vecinos y grados (para normalizar)
    neighbors = {n: list(G.neighbors(n)) for n in nodes}
    degrees = {n: len(neighbors[n]) for n in nodes}

    def multiply_Pt(vec: List[float]) -> List[float]:
        # y = P^T * vec
        y = [0.0] * len(nodes)
        for u in nodes:
            du = degrees[u]
            if du == 0:
                continue
            vu = vec[idx[u]] / du
            for v in neighbors[u]:
                y[idx[v]] += vu
        return y

    # Iteraciones hasta convergencia
    r0 = r[:]
    for _ in range(max_iter):
        Pt_r = multiply_Pt(r)
        new_r = [(1.0 - alpha) * val + alpha * r0_i for val, r0_i in zip(Pt_r, r0)]
        # Criterio de parada
        diff = sum(abs(a - b) for a, b in zip(new_r, r))
        r = new_r
        if diff < tol:
            break

    return {n: r[idx[n]] for n in nodes}


# --------------------------- 2) DIAMOnD (hipergeométrico) --------------------

def _comb(n: int, k: int) -> int:
    # Combinatoria exacta usando Python 3.8+ (maneja enteros grandes)
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)


def _hypergeom_pmf(i: int, N: int, K: int, n: int) -> float:
    # P(X = i) para X ~ Hypergeom(N, K, n)
    # N = total nodos, K = # semilla, n = grado del candidato, i = enlaces del candidato con semillas
    num = _comb(K, i) * _comb(N - K, n - i)
    den = _comb(N, n)
    return 0.0 if den == 0 else num / den


def _hypergeom_sf(k_minus_1: int, N: int, K: int, n: int) -> float:
    # Cola superior P(X >= k) = sum_{i=k}^{min(n,K)} PMF(i)
    k = k_minus_1 + 1
    max_i = min(n, K)
    if k > max_i:
        return 0.0
    s = 0.0
    for i in range(k, max_i + 1):
        s += _hypergeom_pmf(i, N, K, n)
    return s


def diamond_iterative(G: nx.Graph,
                      seed_set: Set[str],
                      steps: int = 50) -> pd.DataFrame:
    """
    Implementación de DIAMOnD (básica):
      - Iterativamente, para cada nodo candidato c (no en seed_set):
          * k = # de vecinos de c que están en seed_set
          * n = grado de c
          * N = # total de nodos de la red
          * K = # actual de semillas
        Compute p = P[X >= k] con X~Hypergeom(N,K,n) (cola superior).
      - Añadir el nodo con p más bajo (más significativo).
      - Repetir 'steps' veces (o hasta agotar candidatos).
    Devuelve un DataFrame con columnas: node, score(=p), step_added.
    """
    N = G.number_of_nodes()
    if N == 0:
        return pd.DataFrame(columns=["node", "score", "step_added"])

    seeds = set(n for n in seed_set if n in G)
    if not seeds:
        raise ValueError("Ninguna semilla está presente en la red para DIAMOnD.")

    ranking: List[Tuple[str, float, int]] = []
    current_seeds = set(seeds)

    for step in range(1, steps + 1):
        K = len(current_seeds)
        candidate_scores: List[Tuple[str, float]] = []
        for c in G.nodes():
            if c in current_seeds:
                continue
            neigh = set(G.neighbors(c))
            k = len(neigh & current_seeds)
            n = len(neigh)
            if n == 0:
                continue
            pval = _hypergeom_sf(k - 1, N=N, K=K, n=n)
            candidate_scores.append((c, pval))

        if not candidate_scores:
            break

        # Seleccionar el mejor (p mínimo); si empata, usa el de mayor k implícitamente por p
        best_node, best_p = min(candidate_scores, key=lambda t: t[1])
        ranking.append((best_node, best_p, step))
        current_seeds.add(best_node)

    return pd.DataFrame(ranking, columns=["node", "score", "step_added"])


# ================================ CLI / Main =================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Network Propagation (RWR y DIAMOnD) sobre redes STRING/GUILD/DIAMOnD."
    )
    p.add_argument("--network", "-n", type=str, required=True,
                   help=r"Ruta a la red (p. ej. data\string_network_filtered_hugo-400.tsv)")
    p.add_argument("--format", "-f", type=str, required=True,
                   choices=["string", "guild", "diamond"],
                   help="Formato de la red de entrada.")
    p.add_argument("--seeds", "-s", type=str, required=True,
                   help=r"Archivo con genes semilla (p. ej. data\genes_seed.txt)")
    p.add_argument("--algo", "-a", type=str, required=True,
                   choices=["rwr", "diamond"],
                   help="Algoritmo a ejecutar: rwr (GUILD-like) o diamond.")
    p.add_argument("--outdir", "-o", type=str, required=True,
                   help=r"Carpeta de salida (p. ej. results\ )")
    # Opciones específicas
    p.add_argument("--min-score", type=int, default=0,
                   help="(STRING) Umbral mínimo de combined_score para cargar aristas.")
    p.add_argument("--steps", type=int, default=50,
                   help="(DIAMOnD) Nº de pasos (nodos a añadir).")
    p.add_argument("--alpha", type=float, default=0.85,
                   help="(RWR) Probabilidad de restart (default: 0.85).")
    p.add_argument("--tol", type=float, default=1e-9,
                   help="(RWR) Tolerancia de convergencia (default: 1e-9).")
    p.add_argument("--max-iter", type=int, default=100,
                   help="(RWR) Máximo de iteraciones (default: 100).")
    p.add_argument("--directed", action="store_true",
                   help="Tratar la red como dirigida (por defecto no).")
    return p


def main():
    args = build_parser().parse_args()

    network_path = Path(args.network)
    seeds_path = Path(args.seeds)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) Cargar red
    print(f"[INFO] Cargando red '{network_path}' (formato: {args.format})…")
    G = build_graph(network_path, fmt=args.format,
                    min_score=args.min_score,
                    undirected=not args.directed)
    print(f"[OK] Nodos: {G.number_of_nodes()}  Aristas: {G.number_of_edges()}")

    # 2) Leer semillas
    seeds = read_seeds(seeds_path)
    print(f"[INFO] Semillas ({len(seeds)}): {', '.join(seeds)}")
    seeds_in = [s for s in seeds if s in G]
    if not seeds_in:
        raise SystemExit(
            "[ERROR] Ninguna semilla coincide con los nodos de la red.\n"
            "Verifica que el tipo de ID de semillas corresponde al de la red.\n"
            "Ejemplos:\n"
            " - STRING (HUGO): semillas como ENO1, PGK1, HK2\n"
            " - GUILD/DIAMOnD (IDs numéricos): semillas deben ser esos IDs."
        )
    if len(seeds_in) < len(seeds):
        missing = set(seeds) - set(seeds_in)
        print(f"[WARN] Semillas no encontradas en la red: {', '.join(map(str, missing))}")
    print(f"[OK] Semillas utilizadas: {', '.join(seeds_in)}")

    # 3) Ejecutar algoritmo
    if args.algo == "rwr":
        print("[INFO] Ejecutando RWR (GUILD-like)…")
        scores = rwr_scores(G, seeds_in, alpha=args.alpha, tol=args.tol, max_iter=args.max_iter)
        df = pd.DataFrame({"node": list(scores.keys()), "score": list(scores.values())})
        df.sort_values("score", ascending=False, inplace=True)
        outfile = outdir / "rwr_scores.csv"
        df.to_csv(outfile, index=False)
        print(f"[OK] Guardado: {outfile} (top 10)\n{df.head(10).to_string(index=False)}")

    elif args.algo == "diamond":
        print("[INFO] Ejecutando DIAMOnD (iterativo con p-valor hipergeométrico)…")
        df = diamond_iterative(G, set(seeds_in), steps=args.steps)
        outfile = outdir / "diamond_ranking.csv"
        df.to_csv(outfile, index=False)
        print(f"[OK] Guardado: {outfile} (top 10)\n{df.head(10).to_string(index=False)}")

    print("[DONE] Propagación completada.")


if __name__ == "__main__":
    main()
