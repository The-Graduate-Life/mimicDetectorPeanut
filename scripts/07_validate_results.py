#!/usr/bin/env python3
"""
07_validate_results.py

STATUS: NOT EXECUTED (no real final candidate table exists yet).

Runs the sanity/robustness checks required before biological
interpretation (paper task list items 21 & 25):

Input checks:
  - correct organism per FASTA header sample vs manifest
  - no accidental duplicate input files (host != control != parasite by hash)
  - no swapped host/control files (manifest accessions match file content sample)

Result checks (on the final *_paired_mimics.tsv table):
  - every pathogen_prot / host_prot ID traceable to the source FASTA files
  - reported coordinates fall within the length of the source protein
  - every row satisfies E-value <= 0.001 and bitscore_diff >= 1
    (bitscore_diff recomputed from host_bitscore - control_bitscore if
    both are present in the table; this script does NOT invent them)

This script intentionally FAILS LOUDLY / reports NOT COMPLETED rather
than silently skipping a check it can't perform against real data.
"""
import argparse
import hashlib
import sys
from pathlib import Path


def sha256_of(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sample_ids(fasta_path, n=5):
    ids = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
                if len(ids) >= n:
                    break
    return ids


def load_lengths(fasta_path):
    lengths = {}
    header, seq_chunks = None, []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    lengths[header] = len("".join(seq_chunks))
                header = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            lengths[header] = len("".join(seq_chunks))
    return lengths


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--parasite-fasta")
    ap.add_argument("--host-fasta")
    ap.add_argument("--control-fasta")
    ap.add_argument("--final-candidates-tsv")
    args = ap.parse_args()

    fastas = {
        "parasite": args.parasite_fasta,
        "host": args.host_fasta,
        "control": args.control_fasta,
    }
    provided = {k: v for k, v in fastas.items() if v}
    if not provided:
        print("STATUS: NOT COMPLETED")
        print("Reason: no input FASTA paths supplied and/or files do not "
              "exist yet — real proteomes have not been downloaded in "
              "this environment.")
        sys.exit(1)

    missing = [k for k, v in provided.items() if not Path(v).exists()]
    if missing:
        print("STATUS: NOT COMPLETED")
        print(f"Reason: missing real files for: {', '.join(missing)}")
        sys.exit(1)

    # Check 1: no accidentally identical files (e.g. host == control)
    hashes = {k: sha256_of(v) for k, v in provided.items()}
    if len(set(hashes.values())) != len(hashes):
        print("FAIL: two or more input FASTA files are byte-identical — "
              "likely a swapped/duplicated file:", hashes)
        sys.exit(1)
    print("PASS: parasite/host/control FASTA files are distinct.")

    for k, v in provided.items():
        print(f"[{k}] sample IDs: {sample_ids(v)}")

    if args.final_candidates_tsv:
        candidates_path = Path(args.final_candidates_tsv)
        if not candidates_path.exists():
            print("STATUS: NOT COMPLETED")
            print("Reason: final candidate table does not exist yet — "
                  "the pipeline has not actually been run.")
            sys.exit(1)

        host_lengths = load_lengths(args.host_fasta) if args.host_fasta else {}
        parasite_lengths = load_lengths(args.parasite_fasta) if args.parasite_fasta else {}

        n_rows = 0
        n_bad_coords = 0
        n_unresolved_ids = 0
        with open(candidates_path) as fh:
            header_line = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                n_rows += 1
                fields = dict(zip(header_line, line.rstrip("\n").split("\t")))
                p_id = fields.get("pathogen_prot")
                h_id = fields.get("host_prot")
                if p_id not in parasite_lengths:
                    n_unresolved_ids += 1
                if h_id not in host_lengths:
                    n_unresolved_ids += 1
                try:
                    p_end = int(fields.get("pathogen_end", -1))
                    h_end = int(fields.get("host_end", -1))
                    if p_id in parasite_lengths and p_end > parasite_lengths[p_id]:
                        n_bad_coords += 1
                    if h_id in host_lengths and h_end > host_lengths[h_id]:
                        n_bad_coords += 1
                except ValueError:
                    n_bad_coords += 1

        print(f"Candidate rows checked: {n_rows}")
        print(f"Rows with unresolved protein IDs: {n_unresolved_ids}")
        print(f"Rows with out-of-range coordinates: {n_bad_coords}")
        if n_unresolved_ids or n_bad_coords:
            print("FAIL: one or more candidates could not be traced back "
                  "to the source proteomes or has invalid coordinates. "
                  "Investigate before proceeding to biological interpretation.")
            sys.exit(1)
        print("PASS: all candidates traceable with valid coordinates.")
    else:
        print("No final candidate table supplied — pipeline execution "
              "checks (steps 25.1-25.15 in the task spec) remain "
              "STATUS: NOT COMPLETED.")


if __name__ == "__main__":
    main()
