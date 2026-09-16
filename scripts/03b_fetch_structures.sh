#!/usr/bin/env bash
# 03b_fetch_structures.sh
#
# REVISION 2: the original version of this script used the Google Cloud
# public-datasets-deepmind-alphafold-v4 bucket via gsutil. That bucket
# allows anonymous LISTING but not anonymous GET (confirmed empirically:
# "Anonymous caller does not have storage.objects.get access" on every
# single object, for both taxids, even though gsutil successfully listed
# and matched all the shard names first). Rather than chase Google Cloud
# auth setup, this version uses the EBI-hosted per-accession endpoint
# instead: https://alphafold.ebi.ac.uk/files/AF-<ACCESSION>-F1-model_v4.pdb
# This is the SAME domain/mechanism already verified working with plain
# curl, no auth, earlier in this project's troubleshooting (the
# /api/prediction/<accession> coverage-check calls that returned clean
# 200/404s). Trade-off: ~111k individual HTTPS requests instead of a
# handful of bulk tars -- slower, but it's the path actually confirmed
# to work, not a third guess.
#
# STATUS: NOT EXECUTED (no network access to alphafold.ebi.ac.uk from the
# sandbox this was written in -- run this on your own machine).

set -euo pipefail

INDIR="../data/mimic_input"
PARALLEL_JOBS=4   # lowered from 8: the host role's first attempt returned
                  # MISSING for all 97687 accessions, most plausibly a
                  # transient EBI rate-limit/throttle triggered right as
                  # the script burst from finishing ~13.7k parasite
                  # requests straight into launching 97.7k host requests
                  # (a manual single-accession retest immediately after
                  # succeeded cleanly, so the endpoint itself is fine).
                  # Lower, gentler concurrency reduces the chance of
                  # re-triggering that.

declare -A ROLES=(
  [aspergillus_flavus]="$INDIR/aspergillus_flavus/aspergillus_flavus.fasta $INDIR/aspergillus_flavus/structures"
  [host]="$INDIR/host/arachis_hypogaea.fasta $INDIR/host/structures"
)

fetch_one() {
  acc="$1"
  outdir="$2"
  out="$outdir/${acc}.pdb.gz"
  [ -s "$out" ] && return 0   # already downloaded (resumable across re-runs)

  # Fast path: AlphaFold DB's current release is v6 database-wide (a
  # coordinated version bump, not staggered per-protein -- confirmed via
  # a live /api/prediction/<acc> response during troubleshooting), so
  # guessing "-model_v6.pdb" directly succeeds for the vast majority of
  # accessions without needing a JSON lookup first. An earlier version of
  # this script hardcoded "_v4", which is now stale (confirmed via a
  # 404 NoSuchKey on a live request) -- that's what this v6 guess fixes.
  url="https://alphafold.ebi.ac.uk/files/AF-${acc}-F1-model_v6.pdb"
  code=$(curl -s -o "$outdir/${acc}.pdb.tmp" -w "%{http_code}" "$url" || echo "000")

  # Fallback path: only for accessions where the v6 guess actually misses
  # (a handful might be on an older/different version, or the guessed
  # naming scheme could shift again in the future). Ask the JSON API for
  # its own authoritative pdbUrl rather than guess a second version
  # number, and use exactly that.
  if [ "$code" != "200" ] || [ ! -s "$outdir/${acc}.pdb.tmp" ]; then
    rm -f "$outdir/${acc}.pdb.tmp"
    real_url=$(curl -s "https://alphafold.ebi.ac.uk/api/prediction/${acc}" \
      | grep -o '"pdbUrl":"[^"]*"' | head -1 | sed -E 's/"pdbUrl":"([^"]*)"/\1/')
    if [ -n "$real_url" ]; then
      code=$(curl -s -o "$outdir/${acc}.pdb.tmp" -w "%{http_code}" "$real_url" || echo "000")
    fi
  fi

  if [ "$code" = "200" ] && [ -s "$outdir/${acc}.pdb.tmp" ]; then
    gzip -c "$outdir/${acc}.pdb.tmp" > "$out"
    rm -f "$outdir/${acc}.pdb.tmp"
    echo "OK $acc"
  else
    rm -f "$outdir/${acc}.pdb.tmp"
    echo "MISSING $acc"
  fi
}
export -f fetch_one

for role in "${!ROLES[@]}"; do
  read -r fasta structdir <<< "${ROLES[$role]}"
  mkdir -p "$structdir"
  echo "[$role] Fetching structures for accessions in $fasta ..."

  grep '^>' "$fasta" | tr -d '>' | sort -u > "../data/metadata/${role}_all_accessions.txt"
  n_total=$(wc -l < "../data/metadata/${role}_all_accessions.txt")
  echo "[$role] $n_total accessions to check."

  # xargs -P for parallelism; -I{} passes each accession to fetch_one.
  # Progress is written to a log so you can watch it with:
  #   tail -f ../logs/03b_${role}_fetch.log
  mkdir -p ../logs
  cat "../data/metadata/${role}_all_accessions.txt" \
    | xargs -P "$PARALLEL_JOBS" -I{} bash -c 'fetch_one "$@"' _ {} "$structdir" \
    | stdbuf -oL tee "../logs/03b_${role}_fetch.log" \
    | awk '{c[$1]++} END{for (k in c) print k, c[k]}'

  ls "$structdir" | sed -E 's/\.pdb(\.gz)?$//' | sort -u > "../data/metadata/${role}_have_structures.txt"
  comm -23 "../data/metadata/${role}_all_accessions.txt" "../data/metadata/${role}_have_structures.txt" \
    > "../data/metadata/${role}_missing_structures.txt"

  n_missing=$(wc -l < "../data/metadata/${role}_missing_structures.txt")
  n_have=$((n_total - n_missing))
  pct=$(awk -v h="$n_have" -v t="$n_total" 'BEGIN{printf "%.1f", (t>0)? 100*h/t : 0}')
  echo "[$role] Coverage: $n_have / $n_total accessions have a structure ($pct%)"
  echo "[$role] Missing-accession list: ../data/metadata/${role}_missing_structures.txt"
  echo "[$role] These proteins cannot be mimicry candidates on either side of the comparison -- document as a limitation, not a bug."
done

echo "Done. Review coverage numbers and *_missing_structures.txt before running 04_run_mimicDetector.sh."
