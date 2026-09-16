#!/usr/bin/env bash
# 04_run_mimicDetector.sh
#
# STATUS: NOT EXECUTED.
#
# Runs the actual mimicDetector Snakemake pipeline against the
# config.yaml produced in step 03, using the primary (optimized)
# parameters specified by the 2026 paper:
#   k = 12, BLASTP, PAM30, wordsize = 2, E-value <= 0.001,
#   host-control bitscore diff >= 1, mean QSASA > 0.75,
#   low-complexity fraction <= 0.50, homologues retained (no
#   full-length homologue-removal step — this is the pipeline's
#   updated default behaviour per the paper; no override is applied).

set -euo pipefail

MIMIC_REPO="../mimicDetector"
CONFIG="../mimicDetector_results/aspergillus_flavus_run/config.yaml"
LOG="../logs/04_run_mimicDetector.log"

if [ ! -f "$CONFIG" ]; then
  echo "ERROR: $CONFIG not found. Run 03_prepare_inputs.sh first (which" \
       "itself requires 01_download_proteomes.sh + real structure files)." >&2
  exit 1
fi

echo "Snakemake dry run (sanity check of the DAG before real execution):"
snakemake -s "$MIMIC_REPO/mimicDetector.smk" --configfile "$CONFIG" -n --rerun-incomplete \
  | tee -a "$LOG"

echo "Executing pipeline..."
snakemake -s "$MIMIC_REPO/mimicDetector.smk" --configfile "$CONFIG" \
  --cores 8 --keep-going --printshellcmds --rerun-incomplete \
  2>&1 | tee -a "$LOG"

echo "Pipeline run complete. Outputs expected under:"
echo "  ../mimicDetector_results/aspergillus_flavus_run/aspergillus_flavus/"
echo "Key intermediate/final files (per repo README naming convention):"
echo "  aspergillus_flavus.aspergillus_oryzae.12mers_blastp.out"
echo "  aspergillus_flavus.arachis_hypogaea.12mers_blastp.out"
echo "  *_hcfiltered_blastp.out   (E-value + bitscore-diff filtered)"
echo "  *_qsasa_averages.tsv      (all candidates + unfiltered mean QSASA)"
echo "  *_q75_filtered.tsv        (QSASA filtered)"
echo "  *_q75_l50_paired_mimics.tsv  (final LCR-filtered candidates)"
