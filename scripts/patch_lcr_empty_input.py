#!/usr/bin/env python3
"""
patch_lcr_empty_input.py

STATUS: NOT EXECUTED (run this once, directly on the machine with the
real psychoscope.py, e.g.:
    python3 patch_lcr_empty_input.py ~/mimicDetector/mimicDetector/scripts/psychoscope.py

Fixes a real crash confirmed on this project's actual data:
    TypeError: unsupported operand type(s) for -: 'object' and 'int64'

Root cause: the QSASA-filtered input file
(...12mers_b1_e001_q75_filtered.tsv) had 0 data rows -- confirmed via
`wc -l` on the real file (36 raw candidate pairs, 0 survived the QSASA
>0.75 threshold on both pathogen and host sides simultaneously). This
is a genuine, valid result (this project's documentation already
anticipated that the A. oryzae control, being extremely close to
A. flavus, could filter out nearly everything), not a data-corruption
bug.

But reading a fully empty file with pd.read_csv(header=None,
names=[...]) makes pandas default every column to object dtype, since
there's no data to infer a numeric type from. That crashes the
arithmetic in the LCR-overlap .query() call a few lines later.

lcr_filter() already has a graceful "empty result" branch (writes a
header-only output file and returns) -- but it's positioned AFTER the
crashing merge/query, only catching emptiness that arises from LCR
filtering itself, not emptiness that was already present in the input.
This patch adds the identical check immediately after reading the
input, before the merge/query ever runs.
"""
import sys
from pathlib import Path

OLD = """        qfilt_df = pd.read_csv(qsasa_filtered, sep='\\t', skip_blank_lines=True,
                               header=None, names=['pathogen_prot', 'pathogen_start',
                                                   'pathogen_end', 'pathogen_qsasa', 
                                                   'pathogen_length', 'host_prot', 
                                                   'host_start', 'host_end', 
                                                   'host_qsasa', 'host_length'])
        # load all host LCR coordinates"""

NEW = """        qfilt_df = pd.read_csv(qsasa_filtered, sep='\\t', skip_blank_lines=True,
                               header=None, names=['pathogen_prot', 'pathogen_start',
                                                   'pathogen_end', 'pathogen_qsasa', 
                                                   'pathogen_length', 'host_prot', 
                                                   'host_start', 'host_end', 
                                                   'host_qsasa', 'host_length'])
        if qfilt_df.empty:
            # No candidates survived QSASA filtering -- a genuine, valid result
            # (confirmed on this project's real data: 0/36 pairs cleared the
            # QSASA threshold on both sides). Reading a fully empty file makes
            # every column default to object dtype (nothing to infer a numeric
            # type from), which would crash the arithmetic in the LCR overlap
            # query below. Handle it the same way the code already handles an
            # empty RESULT after LCR filtering (see the lcrfilt_df.empty branch
            # further down) -- write the same header-only output and return.
            with open(Path(outfile), 'w') as f:
                wrt_str = '\\t'.join(['pathogen_prot', 'pathogen_start', 'pathogen_end', 
                                     'pathogen_sequence', 'pathogen_qsasa', 'host_prot', 
                                     'host_start', 'host_end', 'host_sequence', 'host_qsasa',
                                     'bitscore', 'e_value'])
                f.writelines(wrt_str)
            return
        # load all host LCR coordinates"""

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 patch_lcr_empty_input.py /path/to/psychoscope.py", file=sys.stderr)
        sys.exit(1)

    path = Path(sys.argv[1])
    text = path.read_text()

    n = text.count(OLD)
    if n == 0:
        print("ERROR: expected code block not found. The file may already be patched, "
              "or may not match what this script expects. No changes made.", file=sys.stderr)
        sys.exit(1)
    if n > 1:
        print(f"ERROR: expected exactly 1 occurrence, found {n}. Ambiguous match -- "
              f"no changes made, needs manual review.", file=sys.stderr)
        sys.exit(1)

    backup = path.with_suffix(path.suffix + ".bak4")
    backup.write_text(text)
    print(f"Backed up original to {backup}")

    patched = text.replace(OLD, NEW)
    path.write_text(patched)
    print(f"Patched {path} (1 occurrence, empty-input handling in lcr_filter).")

if __name__ == "__main__":
    main()
