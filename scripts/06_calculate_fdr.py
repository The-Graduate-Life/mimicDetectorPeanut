#!/usr/bin/env python3
"""
06_calculate_fdr.py

STATUS: NOT EXECUTED against real mimicDetector output (demonstrated
only on synthetic toy hit tables below, purely to verify the formula
is implemented correctly).

Implements the paper's empirical target-decoy FDR exactly as specified
(Section 2.4 / Eq. 1-2):

    FDR = c * (Ndecoy / Ntarget)
    c   = 1 / (UPdecoy / UPtarget)

where:
    Ntarget  = number of hits to the target (host) database
    Ndecoy   = number of hits to the decoy (shuffled host) database
    UPtarget = number of unique peptides in the target hit set
    UPdecoy  = number of unique peptides in the decoy hit set

"Unique peptides" here is taken to mean unique k-mer query sequences
that produced a passing hit (i.e. distinct pathogen-derived 12-mer
sequences), since the paper's Methods does not further define
"peptide" in this context beyond citing Kall et al. 2008 / Lee et al.
2021 (proteomics decoy-database literature, where "peptide" =
the query sequence). This interpretation is stated explicitly rather
than assumed silently, and should be re-checked against the actual
mimicDetector source before being reported as the paper's own
definition.

Expected inputs (both are mimicDetector *_hcfiltered_blastp.out-style
tables, i.e. the SAME filtering thresholds -- E-value and bitscore
difference -- already applied, one run against the real host database
and one run against the shuffled decoy database):
    --target-hits   TSV/BLASTP outfmt6 file: pathogen k-mers vs HOST
    --decoy-hits    TSV/BLASTP outfmt6 file: pathogen k-mers vs DECOY
BLASTP outfmt 6 columns (per repo README):
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
"""
import argparse
import csv
from pathlib import Path

BLASTP_OUTFMT6_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]


def load_hits(path):
    rows = []
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            rows.append(dict(zip(BLASTP_OUTFMT6_COLS, row)))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target-hits", required=True)
    ap.add_argument("--decoy-hits", required=True)
    ap.add_argument("--out", default="../results/fdr/fdr_summary.tsv")
    args = ap.parse_args()

    target_path = Path(args.target_hits)
    decoy_path = Path(args.decoy_hits)
    if not target_path.exists() or not decoy_path.exists():
        raise SystemExit(
            "ERROR: real hit tables not found. This script cannot produce "
            "a real FDR value until 04_run_mimicDetector.sh has actually "
            "been executed against real target and decoy databases."
        )

    target_hits = load_hits(target_path)
    decoy_hits = load_hits(decoy_path)

    n_target = len(target_hits)
    n_decoy = len(decoy_hits)

    up_target = len({h["qseqid"] for h in target_hits})
    up_decoy = len({h["qseqid"] for h in decoy_hits})

    if n_target == 0:
        raise SystemExit("ERROR: zero target hits — cannot compute FDR (division by zero).")
    if up_target == 0 or up_decoy == 0:
        raise SystemExit("ERROR: zero unique peptides in target or decoy — cannot compute correction factor c.")

    c = 1.0 / (up_decoy / up_target)
    fdr = c * (n_decoy / n_target)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"Ntarget\t{n_target}\n")
        fh.write(f"Ndecoy\t{n_decoy}\n")
        fh.write(f"UPtarget\t{up_target}\n")
        fh.write(f"UPdecoy\t{up_decoy}\n")
        fh.write(f"c\t{c:.6f}\n")
        fh.write(f"FDR\t{fdr:.6f}\n")

    print(f"Ntarget={n_target}  Ndecoy={n_decoy}  UPtarget={up_target}  UPdecoy={up_decoy}")
    print(f"c   = {c:.6f}")
    print(f"FDR = {fdr:.6f}")
    print(f"Written to {args.out}")


if __name__ == "__main__":
    main()
