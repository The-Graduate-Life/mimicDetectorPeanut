#!/usr/bin/env python3
"""
patch_popscomp_shell.py

STATUS: NOT EXECUTED (run this once, directly on the machine with the
real mimicDetector.smk, e.g.:
    python3 patch_popscomp_shell.py ~/mimicDetector/mimicDetector/mimicDetector.smk

Fixes a real crash found during troubleshooting: the per-structure call
`shell(f'~/POPScomp/POPSC/src/pops --pdb {file} --residueOut --popsOut
{outfile} {z[0]}')` uses Snakemake's shell() helper, which raises
CalledProcessError on a nonzero exit code -- and that exception
propagates up and kills the ENTIRE run_pathogen_popscomp /
run_host_popscomp rule, not just the one bad structure.

Confirmed real trigger: host accession A0A444YTB4's AlphaFold structure
file is a genuine, non-corrupt, well-formed 50-line PDB stub with ZERO
ATOM records (just header/title/compound/source metadata) -- our
download script's "200 + non-empty file" check was too loose to catch
this. pops exits with a real error code on it ("Error: Could not find
atoms in input file"), which crashed the whole host popscomp rule.

This patch: replaces the shell() call with subprocess.run(), which does
NOT raise on a nonzero exit code by default. The existing post-loop
diff logic (from patch_popscomp_gate.py) already correctly classifies
any structure that didn't produce output as "missing" -- this just
stops a single bad structure's nonzero exit from crashing the whole
rule before that diff logic ever runs. Also adds `import subprocess`
at the top of the file if not already present.

Applies to BOTH run_pathogen_popscomp and run_host_popscomp, since both
rules contain the identical shell() call pattern (confirmed by
inspection).
"""
import sys
from pathlib import Path

OLD = "shell(f'~/POPScomp/POPSC/src/pops --pdb {file} --residueOut --popsOut {outfile} {z[0]}')"

NEW = ("subprocess.run(f'~/POPScomp/POPSC/src/pops --pdb {file} --residueOut "
       "--popsOut {outfile} {z[0]}', shell=True)  # subprocess.run does NOT "
       "raise on nonzero exit (unlike shell()), so one bad structure can't "
       "crash the whole rule -- the post-loop diff already classifies it as missing")

IMPORT_OLD = "import pandas as pd"
IMPORT_NEW = "import pandas as pd\nimport subprocess"

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 patch_popscomp_shell.py /path/to/mimicDetector.smk", file=sys.stderr)
        sys.exit(1)

    path = Path(sys.argv[1])
    text = path.read_text()

    n = text.count(OLD)
    if n == 0:
        print("ERROR: expected shell() call not found. The file may already be patched, "
              "or may not match what this script expects. No changes made.", file=sys.stderr)
        sys.exit(1)

    if "import subprocess" in text:
        print("Note: 'import subprocess' already present, not adding again.")
        import_n = 0
    else:
        import_n = text.count(IMPORT_OLD)
        if import_n != 1:
            print(f"ERROR: expected exactly 1 occurrence of '{IMPORT_OLD}' to anchor the "
                  f"import addition, found {import_n}. No changes made.", file=sys.stderr)
            sys.exit(1)

    backup = path.with_suffix(path.suffix + ".bak2")
    backup.write_text(text)
    print(f"Backed up original to {backup}")

    patched = text.replace(OLD, NEW)
    if import_n == 1:
        patched = patched.replace(IMPORT_OLD, IMPORT_NEW)

    path.write_text(patched)
    print(f"Patched {n} occurrence(s) of the shell() call in {path} (expected 2: pathogen + host rules).")
    if n != 2:
        print(f"WARNING: expected exactly 2 occurrences (pathogen + host rules), found {n}. "
              f"Please double check the patched file before running Snakemake.", file=sys.stderr)

if __name__ == "__main__":
    main()
