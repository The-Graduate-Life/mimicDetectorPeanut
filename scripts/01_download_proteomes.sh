#!/usr/bin/env bash
# 01_download_proteomes.sh
#
# STATUS: NOT EXECUTED — the sandbox this project was built in has no
# outbound network access. This script is real and should work as-is
# on a machine with internet access, but every claim about its output
# (protein counts, checksums, release version) must be filled in from
# ACTUAL command output after running it — never assumed.
#
# Downloads the three UniProt reference proteomes selected for this
# project (see data/metadata/proteome_manifest.tsv for accessions and
# justification) and records checksums + provenance.

set -euo pipefail

RAW_DIR="../data/raw"
META_DIR="../data/metadata"
LOG_DIR="../logs"
mkdir -p "$RAW_DIR/parasite" "$RAW_DIR/host" "$RAW_DIR/control" "$LOG_DIR"

DATE_STAMP=$(date -u +%Y-%m-%dT%H:%M:%SZ)
echo "Download run started: $DATE_STAMP" | tee -a "$LOG_DIR/01_download.log"

declare -A PROTEOMES=(
  [parasite]="UP000001875 aspergillus_flavus"
  [host]="UP000289738 arachis_hypogaea"
  [control]="UP000006564 aspergillus_oryzae"
)

# Per-role override for the uniprotkb/stream *query string* (not the
# accession itself). Normally "proteome:<accession>" is correct and
# sufficient. UP000001875 (A. flavus, this strain) is a documented
# exception: it is a legitimate, current, non-reassigned accession
# (confirmed live against /proteomes/UP000001875.json — id, strain,
# and taxid all match, proteinCount=13485) but the "proteome:" search
# filter itself returns zero hits for it, apparently because it is a
# UniProt "Non Reference proteome" built only from an unassembled WGS
# assembly, and non-reference proteomes are not reliably exposed to
# that particular search filter. The individual protein records DO
# exist in searchable UniProtKB, in standard UniProtKB format (real
# accessions like B8N8Q9, full headers) — confirmed live via
# organism_id:332952, which returned 13485 matches, exactly matching
# the proteome's declared protein count. So we query by strain taxid
# instead of by proteome accession for this one role. Re-verify this
# workaround is still needed if UniProt reprocesses this proteome.
declare -A QUERY_OVERRIDE=(
  [parasite]="organism_id:332952"
)

for role in "${!PROTEOMES[@]}"; do
  read -r accession name <<< "${PROTEOMES[$role]}"
  out_fasta="$RAW_DIR/$role/${name}.${accession}.fasta"
  query="${QUERY_OVERRIDE[$role]:-proteome:${accession}}"
  url="https://rest.uniprot.org/uniprotkb/stream?query=${query}&format=fasta&compressed=false"
  entry_url="https://rest.uniprot.org/proteomes/${accession}.json"

  echo "[$role] Downloading $accession ($name) ..." | tee -a "$LOG_DIR/01_download.log"
  if [ -n "${QUERY_OVERRIDE[$role]:-}" ]; then
    echo "[$role] NOTE: using query override '${QUERY_OVERRIDE[$role]}' instead of 'proteome:${accession}' (see comment above PROTEOMES for why)." | tee -a "$LOG_DIR/01_download.log"
  fi

  # Idempotency check: skip a role ONLY if a live sha256 of the current
  # on-disk file matches the checksum recorded for it in provenance.
  # Checking for a matching provenance ROW is not enough on its own — a
  # later failed/partial download can silently overwrite a previously-
  # good file (see the atomic .partial write below for why that's now
  # prevented going forward), leaving a stale provenance row that no
  # longer describes what's actually on disk. Comparing checksums is
  # what actually catches that instead of trusting a stale row.
  if [ -s "$out_fasta" ] && [ -f "$META_DIR/download_provenance.tsv" ]; then
    recorded_checksum=$(awk -F'\t' -v r="$role" -v a="$accession" -v n="$name" -v f="$out_fasta" \
      '$1==r && $2==a && $3==n && $4==f {print $6}' "$META_DIR/download_provenance.tsv" | tail -1)
    live_checksum=$(sha256sum "$out_fasta" | awk '{print $1}')
    if [ -n "$recorded_checksum" ] && [ "$recorded_checksum" = "$live_checksum" ]; then
      echo "[$role] Already downloaded and verified ($out_fasta, checksum matches provenance) — skipping. Delete the file and its provenance row to force a re-download." | tee -a "$LOG_DIR/01_download.log"
      continue
    elif [ -n "$recorded_checksum" ]; then
      echo "[$role] WARNING: $out_fasta exists but its checksum ($live_checksum) does NOT match the last recorded provenance checksum ($recorded_checksum)." | tee -a "$LOG_DIR/01_download.log"
      echo "[$role] The on-disk file was modified/corrupted since it was last successfully recorded (e.g. a later failed download overwrote it). Re-downloading to fix this." | tee -a "$LOG_DIR/01_download.log"
    fi
  fi

  # Pre-flight check against the proteome *entry* endpoint (not the
  # uniprotkb/stream search endpoint). This matters because a stale or
  # reassigned accession does NOT 404 on /uniprotkb/stream — it just
  # silently matches zero sequences (HTTP 200, Content-Length: 0), which
  # is indistinguishable from a transient network hiccup until you check
  # the accession itself. Querying /proteomes/{accession}.json directly
  # tells us whether the ID exists, is redundant/obsolete, or has been
  # superseded, instead of us guessing from an empty search result.
  entry_json=$(curl --http1.1 --fail --retry 3 --retry-delay 5 -sSL "$entry_url" || true)
  if [ -z "$entry_json" ]; then
    echo "[$role] ERROR: proteome accession $accession not found at $entry_url" | tee -a "$LOG_DIR/01_download.log"
    echo "[$role] This usually means the accession was merged, redacted, or reassigned by UniProt." | tee -a "$LOG_DIR/01_download.log"
    echo "[$role] Do NOT guess a replacement ID. Find the current one with, e.g.:" | tee -a "$LOG_DIR/01_download.log"
    echo "[$role]   curl -s 'https://rest.uniprot.org/proteomes/stream?query=organism_id:<TAXID>&format=tsv&fields=upid,organism,organism_id,proteome_type,strain'" | tee -a "$LOG_DIR/01_download.log"
    echo "[$role] Confirm the returned upid against the strain in data/metadata/proteome_manifest.tsv, then update PROTEOMES here, 03_prepare_inputs.sh, the manifest, and README.md together." | tee -a "$LOG_DIR/01_download.log"
    echo "[$role] Not recording provenance for a missing accession." | tee -a "$LOG_DIR/01_download.log"
    continue
  fi
  returned_upid=$(echo "$entry_json" | grep -o '"id"[[:space:]]*:[[:space:]]*"[^"]*"' | head -1 | sed -E 's/.*"([^"]+)"$/\1/')
  if [ -n "$returned_upid" ] && [ "$returned_upid" != "$accession" ]; then
    echo "[$role] ERROR: $accession redirects to a different current accession: $returned_upid" | tee -a "$LOG_DIR/01_download.log"
    echo "[$role] Update PROTEOMES here (and 03_prepare_inputs.sh, the manifest, README.md) to $returned_upid, then re-run." | tee -a "$LOG_DIR/01_download.log"
    echo "[$role] Not recording provenance for a superseded accession." | tee -a "$LOG_DIR/01_download.log"
    continue
  fi

  # --http1.1: work around observed 'HTTP/2 stream reset' errors from
  #            UniProt's streaming endpoint on some networks (a known,
  #            generic curl/nghttp2 interaction issue, not UniProt-specific).
  # --retry 5 --retry-delay 5 --retry-all-errors: retry transient network
  #            failures. Plain --retry only covers a narrow default set
  #            (timeouts, 5xx, some FTP codes) — it does NOT cover a
  #            mid-stream connection drop (curl exit 56 / CURLE_RECV_ERROR,
  #            "SSL_read: unexpected eof while reading"), which is exactly
  #            what was observed on the large host (peanut, ~97k sequence)
  #            download. --retry-all-errors widens retry coverage to that
  #            class of failure too. Bumped retry count from 3 to 5 since
  #            larger responses are more exposed to a mid-stream drop.
  # --fail: make curl return a nonzero exit code on HTTP error responses
  #         (4xx/5xx), instead of silently writing an error page as if it
  #         were the FASTA body.
  #
  # Written to a .partial file first, moved into place only after a full
  # success + non-zero sequence count. This is deliberate: curl's -o opens
  # (and truncates) its output file immediately, before the transfer
  # completes. Writing straight to $out_fasta means a download that fails
  # or resets mid-stream overwrites/corrupts whatever good file was
  # already there — which is exactly what happened to the host proteome
  # here (a good 97687-sequence file got clobbered down to a 1992-sequence
  # partial one by a later failed attempt, silently, with no error at the
  # time it happened). Writing to a temp path and only mv'ing on success
  # means a failed attempt can never touch a previously-good file.
  tmp_fasta="${out_fasta}.partial"
  if ! curl --http1.1 --fail --retry 5 --retry-delay 5 --retry-all-errors -sSL "$url" -o "$tmp_fasta"; then
    echo "[$role] ERROR: download failed (see curl exit code above)." | tee -a "$LOG_DIR/01_download.log"
    echo "[$role] Not recording provenance for a failed/partial download. $out_fasta (if it exists) is untouched." | tee -a "$LOG_DIR/01_download.log"
    rm -f "$tmp_fasta"
    continue
  fi

  n_seqs=$(grep -c '^>' "$tmp_fasta" || true)

  if [ "$n_seqs" -eq 0 ]; then
    echo "[$role] ERROR: downloaded file has 0 sequences (empty/truncated response)." | tee -a "$LOG_DIR/01_download.log"
    echo "[$role] Not recording provenance for an empty file. $out_fasta (if it exists) is untouched." | tee -a "$LOG_DIR/01_download.log"
    rm -f "$tmp_fasta"
    continue
  fi

  mv "$tmp_fasta" "$out_fasta"
  checksum=$(sha256sum "$out_fasta" | awk '{print $1}')

  echo "[$role] $accession -> $out_fasta" | tee -a "$LOG_DIR/01_download.log"
  echo "[$role] sequences: $n_seqs" | tee -a "$LOG_DIR/01_download.log"
  echo "[$role] sha256: $checksum" | tee -a "$LOG_DIR/01_download.log"

  # Append actual provenance to metadata (do not hand-edit protein counts
  # or checksums elsewhere; regenerate this file from real output).
  echo -e "${role}\t${accession}\t${name}\t${out_fasta}\t${n_seqs}\t${checksum}\t${DATE_STAMP}" \
    >> "$META_DIR/download_provenance.tsv"

  # Be polite to UniProt's API and reduce the chance of triggering
  # further stream resets/rate limiting on back-to-back large requests.
  sleep 3
done

echo "Download run complete." | tee -a "$LOG_DIR/01_download.log"
echo "Next: run 02_qc.py against the files listed in $META_DIR/download_provenance.tsv"
