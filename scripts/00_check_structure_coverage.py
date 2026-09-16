#!/usr/bin/env python3
"""
00_check_structure_coverage.py

STATUS: NOT EXECUTED (no network in the environment this was written in).
This is a real, runnable script for YOU to run before investing time in
steps 4-8 of the execution guide.

Purpose: check what fraction of your parasite and host FASTA accessions
already have an AlphaFold DB structure, WITHOUT downloading any actual
structure files yet. Two methods are implemented; use Method 1 first.

METHOD 1 (recommended, fast, bulk): AlphaFold DB publishes a flat list
of every UniProt accession it covers:
    https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/accession_ids.txt
Download it once (~gigabytes of text, one accession per line -- do this
with `wget`/`curl` and `gzip -d` if compressed, not this script), then
point this script at it with --coverage-list. This gives you an exact
answer for every single protein in your proteome in seconds, no rate
limits, no guessing from a sample.

METHOD 2 (fallback, slow, sampled): if you don't want to download the
full accession list, this script can instead query the AlphaFold API
per-accession for a random SAMPLE of your proteome:
    GET https://alphafold.ebi.ac.uk/api/prediction/{accession}
Returns a JSON array (non-empty = structure exists, empty = it doesn't).
This only gives you an estimated coverage percentage with a confidence
interval, not an exact per-protein answer, and is slow / rate-limit-
sensitive for tens of thousands of proteins -- use it only as a quick
gut-check before committing to Method 1's bulk download.

Either way: run this on BOTH the parasite (Aspergillus flavus) and host
(Arachis hypogaea / peanut) FASTA files. The peanut proteome is the one
flagged as uncertain in STATUS.md / README.md -- if its coverage comes
back low, that tells you to plan for AlphaFold3/ColabFold structure
prediction for the uncovered proteins (or to scope the analysis down to
only the subset with existing structures) BEFORE running the rest of
the pipeline.
"""
import argparse
import random
import sys
import time
from pathlib import Path

try:
    import requests
except ImportError:
    requests = None


def extract_accessions(fasta_path):
    """Pull the bare UniProt accession from each FASTA header.
    Handles both '>sp|ACCESSION|NAME...' (UniProt-style) and a bare
    '>ACCESSION ...' header (e.g. after bioawk renaming in 03_prepare_inputs.sh).
    """
    accessions = []
    with open(fasta_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            first_field = header.split()[0] if header else ""
            if "|" in first_field:
                parts = first_field.split("|")
                acc = parts[1] if len(parts) >= 2 else parts[0]
            else:
                acc = first_field
            if acc:
                accessions.append(acc)
    return accessions


def method1_bulk_list(accessions, coverage_list_path):
    """Exact check against the full AlphaFold-covered accession list."""
    covered = set()
    with open(coverage_list_path) as fh:
        for line in fh:
            covered.add(line.strip())
    n_total = len(accessions)
    n_covered = sum(1 for a in accessions if a in covered)
    missing = [a for a in accessions if a not in covered]
    return n_total, n_covered, missing


def method2_api_sample(accessions, sample_size, sleep_seconds):
    if requests is None:
        sys.exit(
            "ERROR: the 'requests' package is required for Method 2 "
            "(pip install requests --break-system-packages)."
        )
    sample = random.sample(accessions, min(sample_size, len(accessions)))
    n_covered = 0
    missing = []
    for acc in sample:
        url = f"https://alphafold.ebi.ac.uk/api/prediction/{acc}"
        try:
            resp = requests.get(url, timeout=15)
        except Exception as exc:
            print(f"  WARN: request failed for {acc}: {exc}", file=sys.stderr)
            continue
        if resp.status_code == 200 and resp.json():
            n_covered += 1
        else:
            missing.append(acc)
        time.sleep(sleep_seconds)  # be polite to the API
    return len(sample), n_covered, missing


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", help="Path to a real proteome FASTA (parasite or host)")
    ap.add_argument("--label", required=True, help="e.g. parasite / host")
    ap.add_argument("--coverage-list",
                     help="Path to a locally downloaded accession_ids.txt "
                          "from the AlphaFold FTP site (Method 1, recommended)")
    ap.add_argument("--sample-api", type=int, default=0,
                     help="If set (and --coverage-list is not), sample this "
                          "many accessions and check each via the live "
                          "AlphaFold API (Method 2, fallback)")
    ap.add_argument("--sleep", type=float, default=0.34,
                     help="Seconds to sleep between API calls in Method 2 "
                          "(default keeps you under ~3 requests/sec)")
    args = ap.parse_args()

    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        sys.exit(f"ERROR: {fasta_path} does not exist.")

    accessions = extract_accessions(fasta_path)
    if not accessions:
        sys.exit(f"ERROR: no FASTA headers/accessions found in {fasta_path}.")

    print(f"[{args.label}] {len(accessions)} accessions found in {fasta_path}")

    if args.coverage_list:
        cov_path = Path(args.coverage_list)
        if not cov_path.exists():
            sys.exit(f"ERROR: {cov_path} does not exist.")
        n_total, n_covered, missing = method1_bulk_list(accessions, cov_path)
        pct = 100 * n_covered / n_total
        print(f"[{args.label}] METHOD 1 (exact, bulk list)")
        print(f"[{args.label}] Covered: {n_covered} / {n_total} ({pct:.1f}%)")
        if missing:
            out_path = fasta_path.parent / f"{args.label}_missing_structures.txt"
            with open(out_path, "w") as fh:
                fh.write("\n".join(missing) + "\n")
            print(f"[{args.label}] Wrote {len(missing)} missing accessions to {out_path}")
    elif args.sample_api:
        n_sampled, n_covered, missing = method2_api_sample(
            accessions, args.sample_api, args.sleep)
        pct = 100 * n_covered / n_sampled if n_sampled else 0.0
        print(f"[{args.label}] METHOD 2 (estimated, API sample of {n_sampled})")
        print(f"[{args.label}] Estimated coverage: {n_covered} / {n_sampled} ({pct:.1f}%)")
        print(f"[{args.label}] NOTE: this is a point estimate from a random "
              f"sample, not an exact count. Re-check with Method 1 before "
              f"relying on it for a go/no-go decision.")
    else:
        sys.exit(
            "ERROR: provide either --coverage-list (Method 1, recommended) "
            "or --sample-api N (Method 2, fallback)."
        )


if __name__ == "__main__":
    main()
