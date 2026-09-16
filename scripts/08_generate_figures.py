#!/usr/bin/env python3
"""
08_generate_figures.py

STATUS: NOT EXECUTED against real results (no real candidate/score data
exists yet in this environment). This script is provided so that, once
real results exist, figures can be produced without inventing data.

Produces the four figures required by the project spec:
  Fig 1: workflow diagram (static; can be generated with no result data —
         this is the only one of the four that does not depend on
         real numbers, so it is the only one safe to render now)
  Fig 2: host vs. negative-control bitscore distributions
         (requires real *_blastp.out files from BOTH host and control runs)
  Fig 3: candidate counts surviving each filter stage
         (requires real counts logged at each pipeline stage)
  Fig 4: empirical FDR vs. bitscore-difference threshold
         (requires re-running 06_calculate_fdr.py at several threshold
         values, i.e. several real mimicDetector runs)

Figures 2-4 will refuse to run without real input files, to avoid
plotting placeholder/fabricated numbers.
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def fig1_workflow(outdir):
    """The only figure that can legitimately be drawn without real
    result data: it depicts the fixed pipeline structure, not results."""
    steps = [
        "Parasite proteome",
        "12-aa overlapping k-mers",
        "BLASTP / PAM30 / wordsize 2\n(vs host & control)",
        "E-value \u2264 0.001",
        "Host-control bitscore diff \u2265 1",
        "Merge overlapping hits",
        "Solvent accessibility > 75%\n(POPSCOMP)",
        "Low-complexity filter\n(Segmasker, \u226450%)",
        "Target-decoy FDR analysis",
        "Final candidate mimics",
    ]
    fig, ax = plt.subplots(figsize=(5, 12))
    ax.axis("off")
    y = 1.0
    dy = 1.0 / (len(steps))
    for i, step in enumerate(steps):
        ax.text(0.5, 1 - i * dy, step, ha="center", va="center",
                 fontsize=9, wrap=True,
                 bbox=dict(boxstyle="round", fc="#e8f0fe", ec="#4472c4"))
        if i < len(steps) - 1:
            ax.annotate("", xy=(0.5, 1 - (i + 0.55) * dy),
                        xytext=(0.5, 1 - (i + 0.15) * dy),
                        arrowprops=dict(arrowstyle="->"))
    ax.set_title("mimicDetector workflow (paper Fig. 1 structure)")
    out_path = Path(outdir) / "figure1_workflow.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Wrote {out_path}")


def fig2_score_distributions(host_scores_path, control_scores_path, outdir):
    if not host_scores_path or not control_scores_path:
        print("SKIPPED Figure 2: real host/control bitscore files not supplied.")
        return
    if not Path(host_scores_path).exists() or not Path(control_scores_path).exists():
        print("SKIPPED Figure 2: real host/control bitscore files do not exist.")
        return
    host_scores = [float(x) for x in Path(host_scores_path).read_text().split()]
    control_scores = [float(x) for x in Path(control_scores_path).read_text().split()]
    fig, ax = plt.subplots()
    ax.hist(host_scores, bins=30, alpha=0.6, label="pathogen vs host")
    ax.hist(control_scores, bins=30, alpha=0.6, label="pathogen vs control")
    ax.set_xlabel("BLASTP bitscore")
    ax.set_ylabel("Count")
    ax.set_title("Host vs. negative-control alignment bitscores")
    ax.legend()
    out_path = Path(outdir) / "figure2_host_vs_control_scores.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Wrote {out_path}")


def fig3_filtering_funnel(counts_path, outdir):
    if not counts_path or not Path(counts_path).exists():
        print("SKIPPED Figure 3: real per-stage candidate counts file not supplied.")
        return
    stages, counts = [], []
    for line in Path(counts_path).read_text().splitlines():
        if not line.strip():
            continue
        stage, count = line.split("\t")
        stages.append(stage)
        counts.append(int(count))
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(stages[::-1], counts[::-1])
    ax.set_xlabel("Candidates remaining")
    ax.set_title("Candidate filtering funnel")
    out_path = Path(outdir) / "figure3_filtering_funnel.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Wrote {out_path}")


def fig4_fdr_curve(fdr_sweep_path, outdir):
    if not fdr_sweep_path or not Path(fdr_sweep_path).exists():
        print("SKIPPED Figure 4: real FDR-sweep file not supplied "
              "(requires multiple real mimicDetector + FDR runs at "
              "different bitscore-difference thresholds).")
        return
    xs, ys = [], []
    for line in Path(fdr_sweep_path).read_text().splitlines():
        if not line.strip():
            continue
        b, fdr = line.split("\t")
        xs.append(float(b))
        ys.append(float(fdr))
    fig, ax = plt.subplots()
    ax.plot(xs, ys, marker="o")
    ax.axhline(0.05, color="red", linestyle="--", label="FDR = 0.05")
    ax.set_xlabel("Bitscore-difference threshold")
    ax.set_ylabel("Empirical FDR")
    ax.set_title("Empirical FDR vs. bitscore-difference threshold")
    ax.legend()
    out_path = Path(outdir) / "figure4_fdr_curve.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Wrote {out_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="../results/figures")
    ap.add_argument("--host-scores", help="whitespace-separated real bitscores, pathogen vs host")
    ap.add_argument("--control-scores", help="whitespace-separated real bitscores, pathogen vs control")
    ap.add_argument("--stage-counts", help="TSV: stage_name<TAB>count, one per pipeline stage")
    ap.add_argument("--fdr-sweep", help="TSV: bitscore_diff<TAB>empirical_FDR")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fig1_workflow(outdir)
    fig2_score_distributions(args.host_scores, args.control_scores, outdir)
    fig3_filtering_funnel(args.stage_counts, outdir)
    fig4_fdr_curve(args.fdr_sweep, outdir)


if __name__ == "__main__":
    main()
