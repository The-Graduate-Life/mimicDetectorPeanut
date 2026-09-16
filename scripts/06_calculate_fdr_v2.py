#!/usr/bin/env python3
"""
06_calculate_fdr_v2.py

Corrected version of 06_calculate_fdr.py.

THE PROBLEM WITH THE ORIGINAL: it computes UPtarget/UPdecoy as the
number of unique qseqid (query) values. But this project's
*_hcfiltered_blastp.out-style inputs (both the real target file and
the decoy-equivalent file built by build_decoy_hcfiltered.py) already
keep exactly ONE row per unique query k-mer -- that's what
filter_blast_bitscore's groupby('query').idxmax() does upstream. So
every qseqid in these files is ALREADY unique by construction, meaning
UPtarget == Ntarget and UPdecoy == Ndecoy will ALWAYS hold for any
input built this way, which makes:
    c = UPtarget/UPdecoy = Ntarget/Ndecoy
    FDR = c * (Ndecoy/Ntarget) = 1.0
...algebraically, unconditionally, regardless of the real data.
Confirmed by running the original script against this project's real
files: FDR=1.000000 exactly, which is a tautology, not a finding.

THE FIX: count unique peptides by SUBJECT (sseqid -- the host/decoy
protein actually hit) instead of query. sseqid is NOT deduplicated
anywhere upstream, so genuine redundancy differences between target
and decoy hit patterns (e.g. many query k-mers piling onto the same
few subjects vs. spread evenly across many) are preserved and can
meaningfully differ between UPtarget and UPdecoy. This still follows
the Kall et al. 2008 / Lee et al. 2021 target-decoy logic the paper
cites (N = total hits, UP = unique underlying identifications among
possibly-redundant hits) -- it just applies "unique" to the dimension
that actually carries redundancy in this pipeline's output, rather
than the dimension that was already deduplicated before this script
ever saw the data.

This is a documented reinterpretation, not a silent guess: the
original script's own docstring already flagged its query-based
"unique peptide" definition as something that "should be re-checked
against the actual mimicDetector source before being reported as the
paper's own definition" -- this is that re-check, done empirically
against this project's real, degenerate FDR=1.0 result.
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
    ap.add_argument("--out", default="../results/fdr/fdr_summary_v2.tsv")
    args = ap.parse_args()

    target_path = Path(args.target_hits)
    decoy_path = Path(args.decoy_hits)
    if not target_path.exists() or not decoy_path.exists():
        raise SystemExit("ERROR: real hit tables not found.")

    target_hits = load_hits(target_path)
    decoy_hits = load_hits(decoy_path)

    n_target = len(target_hits)
    n_decoy = len(decoy_hits)

    # FIX: unique by subject (sseqid), not query (qseqid) -- see module docstring.
    up_target = len({h["sseqid"] for h in target_hits})
    up_decoy = len({h["sseqid"] for h in decoy_hits})

    print(f"[for comparison] query-based UP (original, degenerate): "
          f"UPtarget={len({h['qseqid'] for h in target_hits})} "
          f"UPdecoy={len({h['qseqid'] for h in decoy_hits})} (always == N)")

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
        fh.write(f"UPtarget_by_subject\t{up_target}\n")
        fh.write(f"UPdecoy_by_subject\t{up_decoy}\n")
        fh.write(f"c\t{c:.6f}\n")
        fh.write(f"FDR\t{fdr:.6f}\n")

    print(f"\nNtarget={n_target}  Ndecoy={n_decoy}  "
          f"UPtarget(by subject)={up_target}  UPdecoy(by subject)={up_decoy}")
    print(f"c   = {c:.6f}")
    print(f"FDR = {fdr:.6f}")
    print(f"Written to {args.out}")


if __name__ == "__main__":
    main()
