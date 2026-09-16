#!/usr/bin/env python3
"""
patch_qsasa_h_counter_init.py

STATUS: NOT EXECUTED (run this once, directly on the machine with the
real psychoscope.py, e.g.:
    python3 patch_qsasa_h_counter_init.py ~/mimicDetector/mimicDetector/scripts/psychoscope.py

Fixes a real crash confirmed on this project's actual data, only
reachable now that the get_hpops_files path bug is fixed (i.e. real
host pops data is finally being read):
    UnboundLocalError: cannot access local variable 'h' where it is
    not associated with a value

Root cause: `h` is used as a simple counter ("h += 1") to track how
many times a host protein's aligned region falls outside its pops
file's range, logged via write_log a line above. But `h` is never
initialized anywhere in region_avg_qsasa() -- a genuine, simple
pre-existing bug. It was never triggered before because the entire
host-data branch containing this line was unreachable: get_hpops_files
pointed at a directory with no real files, so every host lookup
silently failed to a FileNotFoundError-caught path before ever reaching
this counter.

Fix: initialize h = 0 at the top of the function, alongside the
existing region_vals_lists initialization.
"""
import sys
from pathlib import Path

OLD = '''    def region_avg_qsasa(self, kmer_coords_dict, ppops_files, hpops_files):
        colnames = ['ResidNe', 'Chain', 'ResidNr', 'iCode', 'Phob/A^2', 'Phil/A^2', 'SASA/A^2', 'Q(SASA)', 'N(overl)', 'Surf/A^2']
        region_vals_lists = []'''

NEW = '''    def region_avg_qsasa(self, kmer_coords_dict, ppops_files, hpops_files):
        colnames = ['ResidNe', 'Chain', 'ResidNr', 'iCode', 'Phob/A^2', 'Phil/A^2', 'SASA/A^2', 'Q(SASA)', 'N(overl)', 'Surf/A^2']
        region_vals_lists = []
        # FIX: h is incremented later ("h += 1") as an "out of pops file range"
        # event counter but was never initialized -- a genuine pre-existing bug,
        # only reachable (and only now confirmed via a live crash) once real host
        # pops data is actually being read (see the get_hpops_files path fix).
        h = 0'''

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 patch_qsasa_h_counter_init.py /path/to/psychoscope.py", file=sys.stderr)
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

    backup = path.with_suffix(path.suffix + ".bak5")
    backup.write_text(text)
    print(f"Backed up original to {backup}")

    patched = text.replace(OLD, NEW)
    path.write_text(patched)
    print(f"Patched {path} (1 occurrence, h counter initialization in region_avg_qsasa).")

if __name__ == "__main__":
    main()
