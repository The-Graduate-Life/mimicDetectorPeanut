#!/usr/bin/env python3
"""
02_qc.py

STATUS: NOT EXECUTED (no real FASTA files exist yet in this environment).

Performs the input QC required before any mimicDetector run:
  - number of proteins
  - total amino-acid count
  - min / max / mean / median protein length
  - number of duplicate sequences
  - number of duplicate identifiers
  - count of ambiguous residues (B, J, O, U, X, Z)
  - malformed records (non-standard characters, no sequence, no header)
  - empty sequences
  - unusually short sequences (< 5 aa, since mimicDetector motifs are
    5-20 aa; shorter than the min mimicry-relevant length)

Does NOT delete anything. Prints a report and writes a QC TSV per file.
Run once per proteome (parasite / host / control) after 01_download_proteomes.sh
has produced real FASTA files.
"""
import argparse
import hashlib
import statistics
import sys
from collections import Counter
from pathlib import Path

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS_AA = set("BJOUXZ")
VALID_AA = STANDARD_AA | AMBIGUOUS_AA


def parse_fasta(path):
    """Minimal FASTA parser; yields (header, sequence)."""
    header, seq_chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("fasta", help="Path to input proteome FASTA")
    ap.add_argument("--label", required=True, help="parasite / host / control")
    ap.add_argument("--outdir", default="../data/processed")
    args = ap.parse_args()

    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        sys.exit(
            f"ERROR: {fasta_path} does not exist. This script cannot run "
            f"until 01_download_proteomes.sh has actually been executed "
            f"on a machine with network access."
        )

    lengths = []
    seq_counter = Counter()
    id_counter = Counter()
    n_empty = 0
    n_malformed = 0
    n_ambiguous_residue_seqs = 0
    n_ambiguous_residues_total = 0
    n_short = 0  # < 5 aa
    total_aa = 0

    for header, seq in parse_fasta(fasta_path):
        seq_id = header.split()[0] if header else ""
        id_counter[seq_id] += 1

        if not seq:
            n_empty += 1
            continue

        seq_upper = seq.upper()
        invalid_chars = set(seq_upper) - VALID_AA
        if invalid_chars:
            n_malformed += 1

        n_ambig = sum(seq_upper.count(a) for a in AMBIGUOUS_AA)
        if n_ambig:
            n_ambiguous_residue_seqs += 1
            n_ambiguous_residues_total += n_ambig

        lengths.append(len(seq_upper))
        total_aa += len(seq_upper)
        seq_counter[seq_upper] += 1

        if len(seq_upper) < 5:
            n_short += 1

    n_proteins = len(lengths)
    n_dup_seqs = sum(c - 1 for c in seq_counter.values() if c > 1)
    n_dup_ids = sum(c - 1 for c in id_counter.values() if c > 1)

    report = {
        "label": args.label,
        "file": str(fasta_path),
        "n_proteins": n_proteins,
        "total_aa": total_aa,
        "min_length": min(lengths) if lengths else None,
        "max_length": max(lengths) if lengths else None,
        "mean_length": round(statistics.mean(lengths), 2) if lengths else None,
        "median_length": statistics.median(lengths) if lengths else None,
        "n_duplicate_sequences": n_dup_seqs,
        "n_duplicate_identifiers": n_dup_ids,
        "n_sequences_with_ambiguous_residues": n_ambiguous_residue_seqs,
        "n_ambiguous_residues_total": n_ambiguous_residues_total,
        "n_malformed_records": n_malformed,
        "n_empty_sequences": n_empty,
        "n_unusually_short_sequences_lt5aa": n_short,
        "file_sha256": hashlib.sha256(fasta_path.read_bytes()).hexdigest(),
    }

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / f"{args.label}_qc.tsv"
    with open(out_path, "w") as fh:
        fh.write("field\tvalue\n")
        for k, v in report.items():
            fh.write(f"{k}\t{v}\n")

    print(f"QC report written to {out_path}")
    for k, v in report.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
