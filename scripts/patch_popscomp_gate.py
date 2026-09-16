#!/usr/bin/env python3
"""
patch_popscomp_gate.py

STATUS: NOT EXECUTED (run this once, directly on the machine with the
real mimicDetector.smk, e.g.:
    python3 patch_popscomp_gate.py ~/mimicDetector/mimicDetector/mimicDetector.smk

Fixes a real bug found during troubleshooting: run_pathogen_popscomp and
run_host_popscomp each only write their flag output file if EVERY input
structure produced a matching popscomp output file
(`sorted(pops_prots) == sorted(struct_prots)`). But POPScomp's own
internal geometry sanity checks (e.g. "Too short atom distance ... A ;
system exit = 0") can cause it to exit 0 -- i.e. "successfully" as far
as Snakemake's shell() is concerned -- WITHOUT writing that structure's
output file. Confirmed on this project's real data: 2 of 13318
A. flavus structures hit this (both due to implausible atom distances,
almost certainly low-confidence/clashing regions of the AlphaFold
model), which silently blocked the pathogen flag file forever, with no
error Snakemake could see.

This patch: still runs popscomp on every structure exactly as before,
but writes the flag as soon as popscomp has been attempted on every
structure, regardless of how many individual structures failed
POPScomp's internal check. Failed accessions are written to
<outdir>/popscomp_failed_structures.txt so they can be excluded from
downstream QSASA filtering and documented as a real limitation -- the
same treatment already given to accessions with no AlphaFold structure
at all, not silently dropped.

Applies to BOTH run_pathogen_popscomp and run_host_popscomp, since both
rules contain byte-identical gate logic (confirmed by inspection).
str.replace() below replaces every occurrence, i.e. both rules, in one
pass.
"""
import sys
from pathlib import Path

OLD = '''        if sorted(pops_prots) == sorted(struct_prots):
            with open(output.flag, "w") as f:
                f.write(f"Number of protein structures processed: {len(pops_prots)}")'''

NEW = '''        missing = sorted(set(struct_prots) - set(pops_prots))
        if missing:
            with open(f'{params.outdir}/../popscomp_failed_structures.txt', 'w') as f:
                f.write('\\n'.join(missing) + '\\n')
        # Flag is written once popscomp has been attempted on every
        # structure, even if some individual structures failed POPScomp's
        # own internal geometry checks (it exits 0 even then, so this
        # cannot be detected via shell() exit code). Failed accessions
        # are logged to popscomp_failed_structures.txt for downstream
        # exclusion + documentation, not silently dropped.
        with open(output.flag, "w") as f:
            f.write(f"Number of protein structures processed: {len(pops_prots)} / {len(struct_prots)} ({len(missing)} failed POPScomp's internal checks -- see popscomp_failed_structures.txt)")'''

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 patch_popscomp_gate.py /path/to/mimicDetector.smk", file=sys.stderr)
        sys.exit(1)

    path = Path(sys.argv[1])
    text = path.read_text()

    n = text.count(OLD)
    if n == 0:
        print("ERROR: expected gate pattern not found. The file may already be patched, "
              "or may not match what this script expects. No changes made.", file=sys.stderr)
        sys.exit(1)

    backup = path.with_suffix(path.suffix + ".bak")
    backup.write_text(text)
    print(f"Backed up original to {backup}")

    patched = text.replace(OLD, NEW)
    path.write_text(patched)
    print(f"Patched {n} occurrence(s) of the gate in {path} (expected 2: pathogen + host rules).")
    if n != 2:
        print(f"WARNING: expected exactly 2 occurrences (pathogen + host rules), found {n}. "
              f"Please double check the patched file before running Snakemake.", file=sys.stderr)

if __name__ == "__main__":
    main()
