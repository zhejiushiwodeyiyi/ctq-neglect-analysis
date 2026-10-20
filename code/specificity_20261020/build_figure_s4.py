#!/usr/bin/env python3
"""Build Supplementary Figure S4 from aggregate specificity reanalysis outputs.

This script does not fit models or read participant-level records. It copies the
reported aggregate contrasts into figure-specific source-data tables, then
renders the same numbers as a two-panel forest plot.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import sys

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "code/figure_qa_tools"))
from audit_panel_alignment import require_matplotlib_panel_alignment


MODEL_INPUT: Path
ANCILLARY_INPUT: Path
SOURCE: Path
FIGURES: Path
QA: Path
STEM = "FigureS4_content_specificity_sensitivity"

MODEL_ORDER = [
    ("OLS", "linear", "Linear OLS", "#315f86", "o"),
    ("fractional_logit", "linear", "Linear fractional logit", "#b65a32", "s"),
    ("OLS", "natural_cubic_spline_3df", "Spline OLS", "#315f86", "o"),
    ("fractional_logit", "natural_cubic_spline_3df", "Spline fractional logit", "#b65a32", "s"),
]
ANCILLARY_ORDER = [
    ("wave2_to_3_unweighted_identity_check", "W2–3 unweighted reference", "#737b82", "o"),
    ("wave2_to_3_IPW_p99_common_prior_linear_WLS", "W2–3, IPW (p99 cap)", "#315f86", "D"),
    ("grade_ge7_existing_common_prior_linear_OLS_rescaled", "W2–3, grade ≥7", "#168277", "^"),
    ("wave1_to_2_common_prior_linear_OLS", "W1–2, early interval", "#b65a32", "s"),
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def assert_contrast(row: dict[str, str]) -> None:
    estimate = float(row["estimate"])
    low = float(row["ci_low"])
    high = float(row["ci_high"])
    if not all(map(math.isfinite, (estimate, low, high))):
        raise ValueError(f"Non-finite estimate or interval in {row['label']}")
    if not low <= estimate <= high:
        raise ValueError(f"Estimate falls outside interval in {row['label']}")
    if int(row["n"]) <= 0 or int(row["schools"]) <= 0:
        raise ValueError(f"Invalid sample count in {row['label']}")


def figure_source_data() -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    models_raw = read_csv(MODEL_INPUT)
    ancillary_raw = read_csv(ANCILLARY_INPUT)
    models_by_key = {(row["model"], row["form"]): row for row in models_raw}
    ancillary_by_key = {row["check"]: row for row in ancillary_raw}
    if len(models_by_key) != 4 or len(ancillary_by_key) != 4:
        raise ValueError("Expected exactly four pre-specified rows in each aggregate file")

    models = []
    for model, form, label, color, marker in MODEL_ORDER:
        row = models_by_key[(model, form)]
        models.append({
            "panel": "a",
            "source_key": f"{model}:{form}",
            "label": label,
            "estimate": row["delta"],
            "ci_low": row["school_low"],
            "ci_high": row["school_high"],
            "n": row["n_participants"],
            "schools": row["schools"],
            "color": color,
            "marker": marker,
            "contrast": "Difference between adjusted follow-up item-mean contrasts at neglect scores 17 vs 10",
            "unit": "DASS-21 per-item mean points (0 to 3)",
            "interval_method": "95% school-cluster CR1, t(df=schools-1), approximate",
        })

    ancillary = []
    for key, label, color, marker in ANCILLARY_ORDER:
        row = ancillary_by_key[key]
        ancillary.append({
            "panel": "b",
            "source_key": key,
            "label": label,
            "estimate": row["delta_anhedonic_minus_autonomic"],
            "ci_low": row["school_cr1_t_low"],
            "ci_high": row["school_cr1_t_high"],
            "n": row["n"],
            "schools": row["schools"],
            "color": color,
            "marker": marker,
            "contrast": "Difference between adjusted follow-up item-mean contrasts at neglect scores 17 vs 10",
            "unit": "DASS-21 per-item mean points (0 to 3)",
            "interval_method": "95% school-cluster CR1, t(df=schools-1), approximate",
        })

    for row in models + ancillary:
        assert_contrast(row)
    if abs(float(models[0]["estimate"]) - float(ancillary[0]["estimate"])) > 1e-12:
        raise ValueError("Unweighted identity check differs from primary OLS linear result")
    if not {r["schools"] for r in models} == {"10"}:
        raise ValueError("Four-model grid must share the same 10-school sample")

    fields = list(models[0])
    for path, rows in [
        (SOURCE / "figS4_content_specificity_model_grid.csv", models),
        (SOURCE / "figS4_content_specificity_ancillary.csv", ancillary),
    ]:
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
    return models, ancillary


def forest_panel(ax, rows: list[dict[str, str]], title: str, show_counts: bool) -> None:
    y = list(range(len(rows), 0, -1))
    ax.axvline(0, color="#68747e", lw=0.8, zorder=0)
    ax.grid(axis="x", color="#e7ebee", lw=0.6, zorder=0)
    for ypos, row in zip(y, rows):
        estimate = float(row["estimate"])
        low = float(row["ci_low"])
        high = float(row["ci_high"])
        ax.plot([low, high], [ypos, ypos], color=row["color"], lw=1.7,
                solid_capstyle="round", zorder=2)
        ax.plot([low, low], [ypos - 0.085, ypos + 0.085], color=row["color"], lw=1.2, zorder=2)
        ax.plot([high, high], [ypos - 0.085, ypos + 0.085], color=row["color"], lw=1.2, zorder=2)
        ax.scatter([estimate], [ypos], s=34, color=row["color"], marker=row["marker"],
                   edgecolor="white", linewidth=0.5, zorder=3)
    labels = []
    for row in rows:
        if show_counts:
            labels.append(f"{row['label']}  (n={int(row['n']):,}; {row['schools']} schools)")
        else:
            labels.append(row["label"])
    ax.set_yticks(y, labels)
    ax.set_ylim(0.5, len(rows) + 0.5)
    ax.set_xlim(-0.032, 0.080)
    ax.set_xticks([-0.02, 0, 0.02, 0.04, 0.06, 0.08])
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.tick_params(axis="y", length=0, pad=8, labelsize=7.6)
    ax.tick_params(axis="x", length=3, labelsize=7.4)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#69757e")
    ax.spines["bottom"].set_linewidth(0.6)
    ax.set_title(title, loc="left", fontsize=9.4, fontweight="bold", pad=11)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build supplementary Figure S4 from local aggregate outputs")
    parser.add_argument("--model-grid", type=Path, required=True, help="four_model_grid.csv")
    parser.add_argument("--ancillary", type=Path, required=True, help="ancillary_contrasts.csv")
    parser.add_argument("--output", type=Path, required=True, help="New output directory outside this repository")
    opts = parser.parse_args()
    global MODEL_INPUT, ANCILLARY_INPUT, SOURCE, FIGURES, QA
    MODEL_INPUT, ANCILLARY_INPUT = opts.model_grid, opts.ancillary
    for path in (MODEL_INPUT, ANCILLARY_INPUT):
        if not path.is_file():
            raise FileNotFoundError(path)
    root = opts.output.resolve()
    if root.is_relative_to(REPO):
        raise ValueError("Keep generated figure and aggregate source data outside the public code repository")
    root.mkdir(parents=True, exist_ok=False)
    SOURCE, FIGURES, QA = root / "source_data", root / "figures", root / "qa"
    SOURCE.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    models, ancillary = figure_source_data()
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7.6,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.linewidth": 0.6,
        "savefig.facecolor": "white",
    })
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.1), sharex=True)
    fig.subplots_adjust(left=0.435, right=0.965, top=0.86, bottom=0.20, hspace=0.42)
    forest_panel(axes[0], models, "a   Model-form grid · W2–3 (n=4,219; 10 schools)", False)
    forest_panel(axes[1], ancillary, "b   Interval, grade and attrition checks", True)
    axes[0].tick_params(axis="x", labelbottom=True)
    axes[1].set_xlabel("Adjusted content contrast for neglect 17 vs 10 (0–3 item-mean points)",
                       fontsize=8.5, labelpad=9)
    fig.text(0.435, 0.94, "Anhedonic-item contrast minus autonomic-item contrast",
             ha="left", va="center", fontsize=10.2, fontweight="bold")
    fig.text(0.435, 0.085,
             "Points and bars: adjusted contrasts and approximate 95% school-cluster CR1 intervals.",
             ha="left", va="center", fontsize=6.6, color="#505c64")
    fig.text(0.435, 0.058,
             "Exploratory, post-result analyses; no multiplicity adjustment.",
             ha="left", va="center", fontsize=6.6, color="#505c64")
    fig.canvas.draw()
    require_matplotlib_panel_alignment(
        fig,
        json_out=QA / f"{STEM}.alignment.json",
        overlay_svg=QA / f"{STEM}.alignment.svg",
        tolerance_pt=1.5,
        gutter_tolerance_pt=1.5,
        strict=True,
    )
    fig.savefig(FIGURES / f"{STEM}.pdf", facecolor="white")
    fig.savefig(FIGURES / f"{STEM}.svg", facecolor="white")
    fig.savefig(FIGURES / f"{STEM}.png", dpi=600, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
