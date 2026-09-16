#!/usr/bin/env python3
"""
build_decoy_hcfiltered.py

Replicates the exact filtering logic in filter_blast_bitscore()
(psychoscope.py, lines 139-184) that produces the real
*_hcfiltered_blastp.out target-hits file -- but with the decoy BLASTP
output standing in for the real host hits, and the SAME pathogen-vs-
control (A. oryzae) hits reused unchanged, since the control side of
the comparison doesn't change for a decoy/null model of the host.

For each side (decoy, control), keeps only the single top-bitscore hit
per query k-mer. Joins decoy top-hits with control top-hits on query
(missing side filled with 0). Computes bits_diff = bitscore_decoy -
bitscore_control. Filters to bits_diff >= bit_diff AND evalue_decoy <=
min_e (same thresholds as config.yaml, matching the real target file).

Output: same 12-column structure as raw BLASTP outfmt6 (no header),
identical in form to the real *_hcfiltered_blastp.out -- ready to pass
as --decoy-hits to 06_calculate_fdr.py.

Usage:
    python3 build_decoy_hcfiltered.py \\
        --decoy-blast /path/to/..._decoy.12mers_blastp.out \\
        --control-blast /path/to/..._blastp.out \\
        --out /path/to/decoy_hcfiltered_blastp.out \\
        --bit-diff 1 --min-e 0.001
"""
import argparse
import pandas as pd
from pathlib import Path

COLNAMES = ['query', 'subject', 'pident', 'length', 'mismatch',
            'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']


def _top_hit_per_query(blastfile):
    df = pd.read_csv(blastfile, sep=r'\s+', header=None)
    df.columns = COLNAMES
    return df.loc[df.groupby('query')['bitscore'].idxmax()]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--decoy-blast', required=True,
                    help='Raw BLASTP outfmt6 output: pathogen k-mers vs DECOY db')
    p.add_argument('--control-blast', required=True,
                    help='Existing raw BLASTP outfmt6 output: pathogen k-mers vs control (real file, unchanged)')
    p.add_argument('--out', required=True)
    p.add_argument('--bit-diff', type=float, default=1, help='config.yaml: bit_diff')
    p.add_argument('--min-e', type=float, default=0.001, help='config.yaml: min_e')
    args = p.parse_args()

    print(f"Loading control (top hit per query)...")
    cont_df = _top_hit_per_query(args.control_blast)
    print(f"  {len(cont_df)} unique control queries")

    print(f"Loading decoy (top hit per query)...")
    decoy_df = _top_hit_per_query(args.decoy_blast)
    print(f"  {len(decoy_df)} unique decoy queries")

    topbits_df = decoy_df.set_index('query').join(
        cont_df.set_index('query'), lsuffix='_decoy', rsuffix='_control')
    topbits_df = topbits_df.fillna(0)
    topbits_df['bits_diff'] = topbits_df['bitscore_decoy'] - topbits_df['bitscore_control']

    filtered_df = topbits_df[
        (topbits_df['bits_diff'] >= float(args.bit_diff)) &
        (topbits_df['evalue_decoy'] <= float(args.min_e))
    ]
    print(f"Passed bits_diff >= {args.bit_diff} and evalue <= {args.min_e}: {len(filtered_df)} rows")

    dcols = filtered_df.columns.str.contains(r'_decoy')
    out_df = filtered_df.iloc[:, dcols]
    out_df.to_csv(Path(args.out), sep='\t', header=False, index=True)
    print(f"Wrote decoy-equivalent hcfiltered file to {args.out}")
    print(f"(Same 12-column structure as raw BLASTP outfmt6 -- pass this as "
          f"--decoy-hits to 06_calculate_fdr.py)")


if __name__ == "__main__":
    main()
