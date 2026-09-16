#!/usr/bin/env bash
# 03_prepare_inputs.sh
#
# STATUS: NOT EXECUTED.
#
# Arranges the downloaded proteomes into the directory layout
# mimicDetector's mimic_configuration.py expects, renames headers to
# bare UniProt accessions (as the repo README recommends, using
# bioawk), and generates the config.yaml via mimic_configuration.py.
#
# PREREQUISITE NOT ADDRESSED HERE: mimicDetector's QSASA/POPSCOMP step
# additionally requires a `structures/` subdirectory under each of
# host/ and each pathogen/ directory, containing one .pdb(.gz) file per
# protein, named to match the FASTA accession. Bulk-downloading/curating
# AlphaFold structures for the full A. flavus and (large, allotetraploid)
# peanut proteomes is a separate, substantial acquisition step and is
# NOT scripted here — see STATUS.md and README.md.

set -euo pipefail

MIMIC_REPO="../mimicDetector"          # git clone https://github.com/Kayleerich/mimicDetector.git
INDIR="../data/mimic_input"
OUTDIR="../mimicDetector_results/aspergillus_flavus_run"

# Directory names below are NOT arbitrary -- verified directly against the
# real mimicDetector source (scripts/utils.py make_config_file, and
# mimicDetector.smk) rather than assumed. The host directory in particular
# MUST be literally named "host" (not the host species name): utils.py
# hardcodes db_file=f"{indir}/host/{host_name}.fasta" and the Snakemake
# rules hardcode f'{config["indir"]}/host/structures/'. Pathogen and
# controls directory names (species-named / "controls") were also checked
# and are correct as originally written.
mkdir -p "$INDIR/aspergillus_flavus/structures"
mkdir -p "$INDIR/host/structures"
mkdir -p "$INDIR/controls"

# Rename FASTA headers to bare UniProt accessions, per repo README:
#   >sp|P04004|VTNC_HUMAN ...  ->  >P04004
# Using plain awk instead of bioawk: bioawk isn't a standard package (it's
# not in apt; normally installed via bioconda/from source), and this
# rewrite doesn't need bioawk's FASTA-aware line-joining -- it only ever
# touches header lines and passes sequence lines through untouched, so a
# plain per-line awk filter is sufficient and has no extra dependency.
rename_headers() {
  awk '/^>/ { split($0, a, "|"); print ">"a[2]; next } { print }' "$1" > "$2"
}

rename_headers ../data/raw/parasite/aspergillus_flavus.UP000001875.fasta \
  "$INDIR/aspergillus_flavus/aspergillus_flavus.fasta"

rename_headers ../data/raw/host/arachis_hypogaea.UP000289738.fasta \
  "$INDIR/host/arachis_hypogaea.fasta"

rename_headers ../data/raw/control/aspergillus_oryzae.UP000006564.fasta \
  "$INDIR/controls/aspergillus_oryzae.fasta"

# Generate the mimicDetector config file using its own configuration script.
python "$MIMIC_REPO/mimic_configuration.py" \
  -p aspergillus_flavus \
  -s arachis_hypogaea \
  -c aspergillus_oryzae \
  -i "$INDIR" \
  -o "$OUTDIR" \
  -f aspergillus_flavus \
  -k 12 \
  -b 1 \
  -e 0.001 \
  -q 0.75 \
  -l 0.50 \
  -t 8

echo "Config written under $OUTDIR/config.yaml"
echo "NOTE: -b/-q here are set to the paper's RECOMMENDED optimal values"
echo "(bitscore diff = 1, min QSASA = 0.75), which differ from the CLI's"
echo "own defaults (-b default 2, -q default 0.50). This substitution is"
echo "intentional and documented per the paper's Section 2.4 findings,"
echo "not a silent deviation. -l is set to 0.50 to match the CLI's"
echo "'max_lcr' semantics for the paper's '>50% low-complexity is removed'"
echo "criterion; this mapping should be re-verified against the actual"
echo "pipeline source at execution time before trusting it."
echo ""
echo "NOTE: these thresholds were benchmarked by the paper only on"
echo "animal-host systems targeting the human proteome. Their validity"
echo "for this fungus-vs-plant system is NOT assumed and must be checked"
echo "empirically via 06_calculate_fdr.py on real output (see README.md)."
