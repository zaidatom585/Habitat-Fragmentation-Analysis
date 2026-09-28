"""
Run the habitat fragmentation analysis for all species in data/species.csv.

Outputs (in results/):
  fragmentation_metrics.csv   one row of metrics per species
  temporal_summary.csv        averages by designation period
  taxon_summary.csv           averages by taxonomic group
  comparison_with_original.csv  script results next to the original
                                Google Earth Pro measurements
  *.png                       charts and habitat maps
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fragmentation import analyze_species

ROOT = Path(__file__).parent
PERIODS = [(1960, 2000, "1960-2000"), (2001, 2012, "2001-2012"), (2013, 2025, "2013-2025")]


def period_of(year: int) -> str:
    for start, end, label in PERIODS:
        if start <= year <= end:
            return label
    return "other"


def plot_maps(species: pd.DataFrame, details: dict, out: Path, cluster_km: float) -> None:
    """Small-multiple maps of every species' habitat, colored by cluster."""
    fig, axes = plt.subplots(2, 5, figsize=(22, 9))
    cmap = plt.get_cmap("tab20")
    for ax, (_, row) in zip(axes.ravel(), species.iterrows()):
        d = details[row.common_name]
        patch_of_polygon = {}
        for patch_idx, members in enumerate(d["patch_members"]):
            for poly_idx in members:
                patch_of_polygon[poly_idx] = patch_idx
        lats = []
        for poly_idx, ring in enumerate(d["polygons_lonlat"]):
            if poly_idx not in patch_of_polygon:
                continue
            step = max(1, len(ring) // 400)          # thin very detailed rings for plotting
            r = ring[::step]
            color = cmap(d["cluster_labels"][patch_of_polygon[poly_idx]] % 20)
            ax.fill(r[:, 0], r[:, 1], color=color, linewidth=0.3, edgecolor=color)
            lats.append(r[:, 1].mean())
        ax.set_aspect(1 / np.cos(np.radians(np.mean(lats))))
        ax.set_title(f"{row.common_name}\n{d['metrics'][f'clusters_{cluster_km:g}km']} clusters, "
                     f"{d['metrics']['patches']} patches", fontsize=10)
        ax.tick_params(labelsize=7)
    fig.suptitle(f"Critical habitat by species (colors = clusters, patches within {cluster_km:g} km)", fontsize=14)
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    plt.close(fig)


def plot_metrics(metrics: pd.DataFrame, out: Path, cluster_km: float) -> None:
    m = metrics.sort_values("largest_patch_km2")
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
    axes[0].barh(m.common_name, m.largest_patch_km2, color="#2a9d8f")
    axes[0].set_xscale("log")
    axes[0].set_xlabel("Largest patch (km², log scale)")
    axes[1].barh(m.common_name, m[f"clusters_{cluster_km:g}km"], color="#e9c46a")
    axes[1].set_xlabel(f"Number of clusters ({cluster_km:g} km rule)")
    axes[2].barh(m.common_name, m.max_nn_km, color="#e76f51")
    axes[2].set_xscale("log")
    axes[2].set_xlabel("Most isolated patch: distance to nearest neighbor (km, log)")
    for ax in axes:
        ax.grid(axis="x", alpha=0.3)
    fig.suptitle("Fragmentation metrics by species", fontsize=14)
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    plt.close(fig)


def plot_temporal(metrics: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
    for _, r in metrics.iterrows():
        ax.scatter(r.designation_year, r.median_patch_km2, s=80, color="#264653")
        ax.annotate(r.common_name, (r.designation_year, r.median_patch_km2),
                    textcoords="offset points", xytext=(6, 4), fontsize=8)
    ax.set_yscale("log")
    ax.set_xlabel("Year critical habitat was designated")
    ax.set_ylabel("Median patch size (km², log scale)")
    ax.set_title("Median patch size by designation year (n = 10, descriptive only)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--kml-dir", default=ROOT / "data" / "kml", type=Path)
    parser.add_argument("--out-dir", default=ROOT / "results", type=Path)
    parser.add_argument("--merge-km", default=0.5, type=float, help="merge polygons closer than this (km)")
    parser.add_argument("--cluster-km", default=10.0, type=float, help="cluster patches closer than this (km)")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    species = pd.read_csv(ROOT / "data" / "species.csv")
    species["designation_year"] = pd.to_datetime(species.designation_date).dt.year

    rows, details = [], {}
    for _, sp in species.iterrows():
        path = args.kml_dir / sp.kml_file
        if not path.exists():
            print(f"Skipping {sp.common_name}: {path.name} not found (see README for download link)")
            continue
        start = time.time()
        res = analyze_species(str(path), merge_km=args.merge_km, cluster_km=args.cluster_km)
        scalar = {k: v for k, v in res.items() if not isinstance(v, (list, np.ndarray))}
        details[sp.common_name] = {**res, "metrics": scalar}
        rows.append({**sp.to_dict(), **scalar})
        print(f"{sp.common_name:<32} {res['patches']:>4} patches  {time.time() - start:5.1f}s")

    metrics = pd.DataFrame(rows)
    metrics["period"] = metrics.designation_year.apply(period_of)
    num_cols = ["total_area_km2", "largest_patch_km2", "smallest_patch_km2", "median_patch_km2", "median_nn_km", "max_nn_km", "extent_km"]
    metrics[num_cols] = metrics[num_cols].round(2)
    metrics["largest_patch_share"] = metrics.largest_patch_share.round(3)
    metrics.drop(columns=["designation_date"]).to_csv(args.out_dir / "fragmentation_metrics.csv", index=False)

    ck = f"clusters_{args.cluster_km:g}km"
    summary_cols = [ck, "largest_patch_km2", "smallest_patch_km2", "median_patch_km2", "median_nn_km", "max_nn_km", "extent_km"]
    temporal = metrics.groupby("period")[summary_cols].mean().round(2)
    temporal.insert(0, "species", metrics.groupby("period").size())
    temporal.to_csv(args.out_dir / "temporal_summary.csv")

    taxon = metrics.groupby("taxon")[summary_cols].mean().round(2)
    taxon.insert(0, "species", metrics.groupby("taxon").size())
    taxon.to_csv(args.out_dir / "taxon_summary.csv")

    original = pd.read_csv(ROOT / "data" / "original_measurements.csv")
    comp = original.merge(metrics[["common_name", ck, "largest_patch_km2", "smallest_patch_km2", "extent_km"]], on="common_name")
    comp = comp.rename(columns={
        "clusters": "clusters_original", ck: f"clusters_script_{args.cluster_km:g}km",
        "largest_fragment": "largest_original", "largest_patch_km2": "largest_script_km2",
        "smallest_fragment": "smallest_original", "smallest_patch_km2": "smallest_script_km2",
        "max_separation_km": "max_separation_original_km", "extent_km": "extent_script_km",
        "unit": "original_area_unit",
    })
    comp.to_csv(args.out_dir / "comparison_with_original.csv", index=False)

    plot_metrics(metrics, args.out_dir / "metrics_by_species.png", args.cluster_km)
    plot_temporal(metrics, args.out_dir / "patch_size_by_year.png")
    if len(details) == len(species):
        plot_maps(species, details, args.out_dir / "habitat_maps.png", args.cluster_km)
    print(f"\nResults written to {args.out_dir}/")


if __name__ == "__main__":
    main()
