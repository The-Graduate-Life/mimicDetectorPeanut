#!/usr/bin/env python3
"""
diagnose_lcr_evalue.py

Read-only diagnostic. Does NOT modify any pipeline files or code.
Replicates the exact LCR-overlap filter and final E-value re-scoring
logic from lcr_filter() in psychoscope.py, against your real
q75_filtered.tsv and LCR bed file, so you can see:
  1. How many candidates survive the LCR-overlap check (self.lcr = 0.5)
  2. For those survivors, the actual bitscore/E-value from the same
     PAM30 local re-alignment the pipeline itself uses -- so you can
     see how close (or far) they came to the min_e = 0.001 cutoff,
     rather than just a final zero.

Usage:
    python3 diagnose_lcr_evalue.py \\
        --qsasa-filtered /path/to/..._q75_filtered.tsv \\
        --lcr-coords /path/to/host-segmask-intervals.bed \\
        --host-fasta /path/to/host.fasta \\
        --pathogen-fasta /path/to/pathogen.fasta

All parameter defaults below are taken directly from this project's
real config.yaml and the hardcoded constants in psychoscope.py's
MimicDetectionII.__init__ (lines 467-471): matrix=PAM30, lamb=0.339,
kconst=0.28, gopen=-15, gext=-3. dbsize/min_e/lcr match config.yaml.
"""
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--qsasa-filtered', required=True)
    p.add_argument('--lcr-coords', required=True)
    p.add_argument('--host-fasta', required=True)
    p.add_argument('--pathogen-fasta', required=True)
    p.add_argument('--lcr', type=float, default=0.5, help='config.yaml: lcr')
    p.add_argument('--min-e', type=float, default=0.001, help='config.yaml: min_e')
    p.add_argument('--dbsize', type=int, default=42904241, help='config.yaml: dbsize')
    p.add_argument('--matrix', default='PAM30')
    p.add_argument('--lamb', type=float, default=0.339)
    p.add_argument('--kconst', type=float, default=0.28)
    p.add_argument('--gopen', type=float, default=-15)
    p.add_argument('--gext', type=float, default=-3)
    args = p.parse_args()

    qfilt_df = pd.read_csv(args.qsasa_filtered, sep='\t', skip_blank_lines=True,
                            header=None, names=['pathogen_prot', 'pathogen_start',
                                                 'pathogen_end', 'pathogen_qsasa',
                                                 'pathogen_length', 'host_prot',
                                                 'host_start', 'host_end',
                                                 'host_qsasa', 'host_length'])
    print(f"[1] QSASA-filtered candidates loaded: {len(qfilt_df)}")

    seg_df = pd.read_csv(args.lcr_coords, sep='\t', header=None,
                          names=['host_prot', 'seg_start', 'seg_end'])
    print(f"    LCR interval rows loaded: {len(seg_df)}")

    # Exact same drop-based filter as lcr_filter()
    to_drop_idx = (qfilt_df.reset_index()
                   .merge(seg_df, how='left', on=['host_prot'])
                   .query('(host_start < seg_end) & (host_end > seg_start) & '
                          '(((host_end - seg_start)/(host_end - host_start + 1)) > @args.lcr)')
                   )['index']
    lcrfilt_df = qfilt_df.drop(to_drop_idx.unique())
    print(f"[2] Survived LCR-overlap filter (> {args.lcr*100:.0f}% overlap removed): {len(lcrfilt_df)}")

    if lcrfilt_df.empty:
        print("    -> Zero survivors. The LCR filter itself removed everything;")
        print("       the final E-value re-scoring step never even runs on real candidates.")
        return

    def _seq_records(fasta_path, prot_list):
        record_dict = {}
        prot_set = set(prot_list)
        for record in SeqIO.parse(fasta_path, 'fasta'):
            if record.id in prot_set:
                record_dict[record.id] = record.seq
        return record_dict

    precords = _seq_records(args.pathogen_fasta, lcrfilt_df['pathogen_prot'].to_list())
    hrecords = _seq_records(args.host_fasta, lcrfilt_df['host_prot'].to_list())

    def _subseq(prot, start, end, record_dict):
        return str(record_dict[prot][start - 1:end])

    lcrfilt_df = lcrfilt_df.copy()
    lcrfilt_df['pathogen_sequence'] = lcrfilt_df.apply(
        lambda x: _subseq(x.pathogen_prot, x.pathogen_start, x.pathogen_end, precords), axis=1)
    lcrfilt_df['host_sequence'] = lcrfilt_df.apply(
        lambda x: _subseq(x.host_prot, x.host_start, x.host_end, hrecords), axis=1)

    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load(args.matrix)
    aligner.open_gap_score = args.gopen
    aligner.extend_gap_score = args.gext

    def _aln(seq1, seq2):
        bits = aligner.score(seq1, seq2)
        bitsc = ((args.lamb * bits) - np.log(args.kconst)) / np.log(2)
        evalue = args.kconst * len(seq1) * args.dbsize * (2.71828 ** (-args.lamb * bits))
        return bitsc, evalue

    lcrfilt_df[['bitscore', 'e_value']] = lcrfilt_df.apply(
        lambda x: _aln(x.pathogen_sequence, x.host_sequence), axis=1, result_type='expand')

    print(f"\n[3] Re-alignment results for all {len(lcrfilt_df)} LCR-survivors "
          f"(threshold: e_value <= {args.min_e}):\n")
    show_cols = ['pathogen_prot', 'pathogen_start', 'pathogen_end',
                 'host_prot', 'host_start', 'host_end', 'bitscore', 'e_value']
    sorted_df = lcrfilt_df[show_cols].sort_values('e_value')
    pd.set_option('display.width', 200)
    pd.set_option('display.max_rows', None)
    print(sorted_df.round({'bitscore': 1, 'e_value': 5}).to_string(index=False))

    n_pass = (lcrfilt_df['e_value'] <= args.min_e).sum()
    print(f"\n[4] Of these, {n_pass} clear the e_value <= {args.min_e} threshold.")
    if n_pass == 0:
        closest = sorted_df.iloc[0]
        print(f"    Closest miss: e_value={closest['e_value']:.5f} "
              f"({closest['pathogen_prot']} vs {closest['host_prot']}, "
              f"bitscore={closest['bitscore']:.1f})")


if __name__ == "__main__":
    main()
