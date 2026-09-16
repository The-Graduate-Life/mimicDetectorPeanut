#!/usr/bin/env python3
"""
patch_hpops_output_dir.py

STATUS: NOT EXECUTED (run this once, directly on the machine with the
real mimicDetector.smk, e.g.:
    python3 patch_hpops_output_dir.py ~/mimicDetector/mimicDetector/mimicDetector.smk

THE REAL ROOT CAUSE, found after tracing why every single row of
qsasa_averages.tsv showed h_avg = "NaN" with zero exceptions (36/36),
even for host proteins whose popscomp run genuinely succeeded (of the
~95,501 that did, out of ~95,753 with AlphaFold coverage).

get_hpops_files() builds the list of expected host popscomp OUTPUT
file paths using the literal string 'host':
    outstr=f'{wildcards.outdir}/host/pops/pops_'
But the REAL output directory, written by the actual run_host_popscomp
rule (confirmed at line ~145 of this same file:
    outdir=f'{config["outdir"]}/{config["host_name"]}/pops'
), is named after config['host_name'] (in this project: arachis_hypogaea),
NOT the literal string 'host'. 'host' is only correct for the INPUT
structures directory (data/mimic_input/host/structures/ -- a real,
separate requirement of this pipeline, set up earlier in this project),
not the output.

This meant every path in the list get_hpops_files() returns pointed to
a directory that has never contained any real files for this run --
so EVERY host pops-file lookup during QSASA averaging silently failed
(caught by the earlier FileNotFoundError patch to psychoscope.py,
which was correct in isolation but ended up masking this much larger,
systemic bug) and defaulted to "NaN", regardless of whether that
structure's popscomp run actually succeeded. The pathogen-side
equivalent (get_ppops_files) does NOT have this bug -- it consistently
uses wildcards.pathogen for both input and output, since pathogen
input/output directories share the same real species-name convention.

Fix: use config["host_name"] instead of the literal string 'host' when
building the OUTPUT path guess, matching the real run_host_popscomp
rule's own convention. The INPUT indir line is correctly left
unchanged (it must stay 'host', matching this pipeline's real
structures-directory requirement).

IMPORTANT: after applying this fix, filter_qsasa must be RE-RUN from
scratch (its existing qsasa_averages.tsv / q75_filtered.tsv outputs
were computed entirely from the broken path and are not valid --
delete or let Snakemake detect the rule change and rebuild them).
"""
import sys
from pathlib import Path

OLD = '''def get_hpops_files(wildcards):
    indir=f'{config["indir"]}/host/structures/'
    outstr=f'{wildcards.outdir}/host/pops/pops_'
    hstruct=[''.join([indir, f]) for f in os.listdir(indir) if Path(''.join([indir, f])).exists()]
    hpops=[p.replace(indir, outstr).replace('.pdb', '.out').replace('.gz', '') for p in hstruct]
    return hpops'''

NEW = '''def get_hpops_files(wildcards):
    indir=f'{config["indir"]}/host/structures/'
    # FIX: was f'{wildcards.outdir}/host/pops/pops_' -- the literal 'host' is only
    # correct for the INPUT structures dir. The real run_host_popscomp rule writes
    # output to {outdir}/{host_name}/pops (see its own outdir param above), so the
    # OUTPUT path guess must use config["host_name"], not the literal string 'host'.
    # Confirmed real bug: every host QSASA lookup silently failed and defaulted to
    # NaN for the entire run because of this mismatch, regardless of whether that
    # structure's popscomp run actually succeeded.
    outstr=f'{wildcards.outdir}/{config["host_name"]}/pops/pops_'
    hstruct=[''.join([indir, f]) for f in os.listdir(indir) if Path(''.join([indir, f])).exists()]
    hpops=[p.replace(indir, outstr).replace('.pdb', '.out').replace('.gz', '') for p in hstruct]
    return hpops'''

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 patch_hpops_output_dir.py /path/to/mimicDetector.smk", file=sys.stderr)
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
    print(f"Patched {path} (1 occurrence, get_hpops_files output directory).")
    print("")
    print("IMPORTANT: filter_qsasa's existing outputs were computed from the broken")
    print("path and must be re-run. Delete them before relaunching:")
    print("  rm -f <outdir>/<pathogen>/<pathogen>.<pathogen>.*qsasa_averages.tsv")
    print("  rm -f <outdir>/<pathogen>/<pathogen>.<pathogen>.*q75_filtered.tsv")
    print("  rm -f <outdir>/<pathogen>/<pathogen>.<pathogen>.*hcfiltered_blastp.out")
    print("(or just let --rerun-incomplete / a rule-code-change trigger handle it,")
    print("since this file's content changed, Snakemake should detect and rebuild.)")

if __name__ == "__main__":
    main()
