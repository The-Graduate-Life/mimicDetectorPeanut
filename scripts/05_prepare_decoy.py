#!/usr/bin/env python3
"""
05_prepare_decoy.py

STATUS: NOT EXECUTED against real data (demonstrated only on a tiny
synthetic toy FASTA — see results/raw/toy_example/).

Builds a shuffled-sequence decoy database from the HOST proteome
(peanut, *Arachis hypogaea*, for this project), per the paper's
target-decoy FDR strategy (Section 2.4), generalized from the paper's
own human-proteome case:
  "target = host proteome; decoy = shuffled host sequences"

Design choices made explicit here (the paper's Methods text does not
give shuffling implementation detail beyond "shuffled human sequences",
so these choices are documented rather than silently assumed):
  - Per-sequence shuffling (each host sequence is shuffled independently),
    which PRESERVES per-sequence length and overall amino-acid
    composition (it's a permutation of the same residues) but does NOT
    preserve any local motif/domain structure — that is the point of a
    decoy.
  - Identifiers are changed to "<original_id>_decoy" so decoy hits can
    be distinguished unambiguously in downstream BLASTP output.
  - A fixed random seed is used and recorded, so the decoy is
    reproducible from the same input FASTA.
  - Unique-peptide counting (UPtarget, UPdecoy) is handled downstream
    in 06_calculate_fdr.py, not here; this script only builds the decoy
    sequence file.

This does NOT claim to reproduce the paper's own decoy exactly, since
the paper does not specify a per-residue vs per-sequence vs global
shuffling algorithm, nor whether a fixed seed was used. That is an
explicit, documented limitation (see report/final_report.md, Limitations).
"""
import argparse
import random
from pathlib import Path


def parse_fasta(path):
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


def shuffle_sequence(seq, rng):
    chars = list(seq)
    rng.shuffle(chars)
    return "".join(chars)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("host_fasta")
    ap.add_argument("--out", required=True, help="Output decoy FASTA path")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    host_path = Path(args.host_fasta)
    if not host_path.exists():
        raise SystemExit(f"ERROR: {host_path} does not exist.")

    rng = random.Random(args.seed)
    n_written = 0
    with open(args.out, "w") as out_fh:
        for header, seq in parse_fasta(host_path):
            seq_id = header.split()[0] if header else f"seq{n_written}"
            shuffled = shuffle_sequence(seq.upper(), rng)
            out_fh.write(f">{seq_id}_decoy\n{shuffled}\n")
            n_written += 1

    print(f"Wrote {n_written} shuffled decoy sequences to {args.out}")
    print(f"Random seed used: {args.seed} (record this in software_versions.txt)")
    print("Shuffling method: per-sequence residue permutation")
    print("Length preserved: YES (permutation of the same residues)")
    print("Amino-acid composition preserved: YES, per-sequence")
    print("Identifiers: '<original_id>_decoy'")


if __name__ == "__main__":
    main()
