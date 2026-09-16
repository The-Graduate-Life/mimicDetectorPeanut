#!/usr/bin/env python3
"""
patch_qsasa_pfile_missing.py

STATUS: NOT EXECUTED (run this once, directly on the machine with the
real psychoscope.py, e.g.:
    python3 patch_qsasa_pfile_missing.py ~/mimicDetector/mimicDetector/scripts/psychoscope.py

Fixes a real, pre-existing bug, symmetric to the one already fixed on
the host side (patch_qsasa_hfile_match.py). Confirmed via a live crash
on this project's real data:
    FileNotFoundError: [Errno 2] No such file or directory:
    '.../aspergillus_flavus/pops/pops_A0A7U2MMS4.out'

Root cause (traced to mimicDetector.smk's get_ppops_files/
get_hpops_files): the list of "pops files" passed into
region_avg_qsasa() is NOT an actual glob of existing popscomp output.
It's mechanically derived from the INPUT STRUCTURES directory listing,
renaming each structure filename into an assumed corresponding output
path -- with no check that popscomp actually succeeded for that
accession. A0A7U2MMS4 is one of two pathogen accessions confirmed
earlier in this project to fail POPScomp's own internal geometry
checks (exits 0 but writes no output); its structure file exists, so
it appears in this list, but its pops output file does not.

The pathogen-side read loop here only catches
pd.errors.EmptyDataError, not FileNotFoundError -- so it crashes hard
on any such accession instead of skipping it gracefully, the same way
an empty pops file is already skipped.

Fix: add FileNotFoundError to the caught exceptions, mirroring the fix
already applied to the host-side read.
"""
import sys
from pathlib import Path

OLD = '''                # open pops file, make dataframe
                try: 
                    prot_sasa = pd.read_csv(pfile, sep="\\s+", header=None, engine='python', skiprows=3, skipfooter=3)
                except pd.errors.EmptyDataError:
                    continue'''

NEW = '''                # open pops file, make dataframe
                try: 
                    prot_sasa = pd.read_csv(pfile, sep="\\s+", header=None, engine='python', skiprows=3, skipfooter=3)
                except (pd.errors.EmptyDataError, FileNotFoundError):
                    # FileNotFoundError added: ppops_files (see get_ppops_files in
                    # mimicDetector.smk) is mechanically derived from the structures
                    # directory, not validated against actual popscomp output -- an
                    # accession whose structure exists but whose popscomp run
                    # internally failed (confirmed real case: A0A7U2MMS4) appears
                    # here as if its output exists. Skip it, same as an empty file.
                    continue'''

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 patch_qsasa_pfile_missing.py /path/to/psychoscope.py", file=sys.stderr)
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

    backup = path.with_suffix(path.suffix + ".bak3")
    backup.write_text(text)
    print(f"Backed up original to {backup}")

    patched = text.replace(OLD, NEW)
    path.write_text(patched)
    print(f"Patched {path} (1 occurrence, pathogen pops file reading in region_avg_qsasa).")

if __name__ == "__main__":
    main()
