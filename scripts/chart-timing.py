#!/usr/bin/env python3
"""Generate a stacked bar chart from pipeline timing.json files.

Usage: chart-timing.py OUTPUT_PNG LABEL1=timing1.json [LABEL2=timing2.json ...]

Each argument after the output path is a label=path pair pointing to a
timing.json produced by `maxbin-rs pipeline`.
"""

import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


STAGES = [
    ("genecall_s", "Gene calling"),
    ("hmmer_s", "HMMER"),
    ("recursive_em_s", "EM + recursion"),
]

COLORS = ["#b2df8a", "#33a02c", "#fb9a99"]


def load_timing(path):
    with open(path) as f:
        return json.load(f)


def main():
    if len(sys.argv) < 3:
        print(__doc__, file=sys.stderr)
        sys.exit(1)

    out_path = sys.argv[1]
    runs = []
    for arg in sys.argv[2:]:
        label, path = arg.rsplit("=", 1)
        runs.append((label, load_timing(path)))

    labels = [r[0] for r in runs]
    x = np.arange(len(labels))
    width = 0.5

    fig, ax = plt.subplots(figsize=(max(4, len(runs) * 1.5 + 1), 5))

    bottoms = np.zeros(len(runs))
    for (key, stage_label), color in zip(STAGES, COLORS):
        values = [r[1].get(key, 0.0) for r in runs]
        bars = ax.bar(x, values, width, bottom=bottoms, label=stage_label, color=color)
        # Label segments that are wide enough to read.
        for bar, val in zip(bars, values):
            if val > 0.5:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_y() + bar.get_height() / 2,
                    f"{val:.1f}s",
                    ha="center", va="center", fontsize=8,
                )
        bottoms += values

    ax.set_ylabel("Wall-clock time (seconds)")
    ax.set_title("maxbin-rs pipeline timing")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
