#!/usr/bin/env python3
"""
patch_qsasa_hfile_match.py

STATUS: NOT EXECUTED (run this once, directly on the machine with the
real psychoscope.py, e.g.:
    python3 patch_qsasa_hfile_match.py ~/mimicDetector/mimicDetector/scripts/psychoscope.py

Fixes a real, pre-existing bug in region_avg_qsasa(), confirmed via a
live crash on this project's real data:
    FileNotFoundError: [Errno 2] No such file or directory:
    '.../arachis_hypogaea/pops/pops_A0A445DSU1.out'

Root cause: line "hfiles = [file for file in hpops_files if h_prot_name
in file]" does a SUBSTRING match, not an exact match. If a failing
accession (e.g. A0A445DSU1, which has no popscomp output because its
structure failed AlphaFold coverage or POPScomp's own checks) happens
to be a substring of a different, successfully-processed accession's
filename (e.g. A0A445DSU12), hfiles is wrongly non-empty. The code then
throws that (wrong) match away and manually reconstructs the path for
the ORIGINAL, genuinely-missing accession, and crashes trying to open
it.

This bug was always latent but harmless before: the pipeline's old
all-or-nothing popscomp gate guaranteed every host protein had a pops
file, so this path was never exercised with a real gap. Relaxing that
gate (a separate, correct fix applied earlier in this project) is what
first exposed it.

Fix: use the actual matched file (found via an EXACT basename match,
not substring containment) instead of reconstructing a path by hand,
and treat "no exact match found" the same way an already-handled empty
pops file is treated -- skip that host protein's QSASA values (NaN),
log nothing further needed since this is expected/common now, and
continue -- instead of crashing.
"""
import sys
from pathlib import Path

OLD = '''                        if hfiles:
                            hfile = f"{self.outdir}/{self.host_name}/pops/pops_{h_prot_name}.out"
                            try: 
                                h_prot_sasa = pd.read_csv(f'{hfile}', sep="\\s+", header=None, engine='python', skiprows=3, skipfooter=3)
                            except pd.errors.EmptyDataError:
                                h_region_vals = 0
                                h_region_vals_lists.append(h_region_vals)
                                h_region_avg = "NaN"
                                continue'''

NEW = '''                        # Use an EXACT basename match, not substring containment (the
                        # original "h_prot_name in file" check could false-positive on
                        # accessions that are substrings of other, longer accessions,
                        # e.g. A0A445DSU1 matching A0A445DSU12's filename). Confirmed via
                        # a live crash on this project's real data before this fix.
                        exact_hfiles = [f for f in hfiles if Path(f).name == f'pops_{h_prot_name}.out']
                        if exact_hfiles:
                            hfile = exact_hfiles[0]
                            try: 
                                h_prot_sasa = pd.read_csv(f'{hfile}', sep="\\s+", header=None, engine='python', skiprows=3, skipfooter=3)
                            except (pd.errors.EmptyDataError, FileNotFoundError):
                                h_region_vals = 0
                                h_region_vals_lists.append(h_region_vals)
                                h_region_avg = "NaN"
                                continue'''

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 patch_qsasa_hfile_match.py /path/to/psychoscope.py", file=sys.stderr)
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

    backup = path.with_suffix(path.suffix + ".bak")
    backup.write_text(text)
    print(f"Backed up original to {backup}")

    patched = text.replace(OLD, NEW)
    path.write_text(patched)
    print(f"Patched {path} (1 occurrence, host pops file matching in region_avg_qsasa).")

if __name__ == "__main__":
    main()
