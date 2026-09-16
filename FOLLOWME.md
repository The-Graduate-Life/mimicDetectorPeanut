# FOLLOWME: Step-by-Step Reproduction Guide

This guide reproduces the *A. flavus* vs. peanut molecular mimicry analysis from a clean environment. It reflects the actual commands, scripts, and fixes used to produce the real results in `report/final_report.md` — including patches to the vendored mimicDetector code that are **required** for correct results, not optional cleanup.

**Read Section 16 (Required Patches) before running the primary analysis (Section 8).** The unpatched, upstream mimicDetector code will run to completion without crashing on some inputs, but will silently produce incorrect host-side QSASA values (all `NaN`) — see Section 16 for why.

---

## 1. Prerequisites

- Linux or WSL2 (this analysis was run on WSL2, kernel 6.6.87.2-microsoft-standard-WSL2).
- conda or mamba.
- ~15 GB free disk space (proteome downloads, AlphaFold structures, BLASTP output, and popscomp output for ~111,000 proteins together are substantial).
- Network access to UniProt REST API and the AlphaFold DB EBI file server.
- **At least 8 GB of RAM available**, more than the ~3.7 GB this analysis had access to. The QSASA-filtering step (`filt_qsasa.py`) peaked at ~3.4 GB RSS on a ~13.7k x ~97.7k protein comparison; this analysis required a 6 GB swapfile to complete without being killed. See Section 15.
- If running under WSL2: disable the VM idle-timeout before starting any multi-hour step (Section 15) — the default behavior can silently kill detached background jobs.

## 2. Repository setup

```bash
git clone https://github.com/The-Graduate-Life/mimicDetectorPeanut mimicDetector_project
cd mimicDetector_project
git clone https://github.com/Kayleerich/mimicDetector.git mimicDetector
cd mimicDetector
git checkout c540df77788d8296f62ba662f5e3b673d9c0ec34
cd ..
```

## 3. Environment / software installation

```bash
conda env create -f environment/environment.yml
conda activate mimics
blastp -version   # confirm 2.16.0
snakemake --version   # confirm 7.32.4
```

Install POPScomp separately (not distributed via conda):

```bash
git clone https://github.com/Fraternalilab/POPScomp.git ~/POPScomp
cd ~/POPScomp
# follow the POPScomp README's build instructions (bootstrap/configure/make);
# confirm the resulting binary at ~/POPScomp/POPSC/src/pops
```

Apply the required code patches now (Section 16) before proceeding — the primary run in Section 8 depends on them.

## 4. Directory preparation

mimicDetector requires the host input directory to be named literally `host/` (not the species name) — this is a real constraint of the tool's input-discovery logic, not a naming choice made in this project.

```bash
mkdir -p data/mimic_input/{host,controls}/structures
mkdir -p data/mimic_input/aspergillus_flavus/structures
mkdir -p data/raw/{parasite,host,control}
mkdir -p data/metadata data/processed
mkdir -p mimicDetector_results/aspergillus_flavus_run
```

## 5. Input data acquisition

```bash
cd scripts
bash 01_download_proteomes.sh
```

Downloads the three UniProt proteomes (pathogen UP000001875, host UP000289738, control UP000006564) via the UniProt REST API, with checksum-based idempotency. Real counts from this run: 13,729 pathogen sequences, 97,687 host sequences, 12,061 control sequences.

Optional, before investing time in structure fetching: check AlphaFold coverage in advance without downloading structures yet.

```bash
python3 00_check_structure_coverage.py --label pathogen --coverage-list <path-to-alphafold-accession-list> data/raw/parasite/aspergillus_flavus.UP000001875.fasta
```

## 6. Data quality control

```bash
python3 02_qc.py
```

Checks for duplicate sequences and malformed records in each downloaded proteome. Real result: 0 duplicates, 0 malformed records in all three; 63 host sequences (238 residue positions) contain ambiguous (`X`) residues.

## 7. Preprocessing

```bash
bash 03_prepare_inputs.sh
```

Builds `data/mimic_input/` from the raw downloads, including the literal `host/` directory rename (Section 4) and `config.yaml`. Verify `dbsize` in the generated config equals host total AA + control total AA.

```bash
bash 03b_fetch_structures.sh
```

Fetches AlphaFold structures per-accession from `https://alphafold.ebi.ac.uk/files/AF-{ACCESSION}-F1-model_v6.pdb`. Real coverage: 97.0% pathogen (13,318/13,729), 98.0% host (95,753/97,687) — the gap in both cases is fully accounted for by AlphaFold's exclusion of sequences >2,700 aa.

## 8. Primary analysis

**Confirm the patches in Section 16 are applied to `mimicDetector/mimicDetector.smk` and `mimicDetector/scripts/psychoscope.py` before running this step.**

```bash
bash 04_run_mimicDetector.sh
```

This runs the full Snakemake pipeline (BLASTP, POPScomp, QSASA filtering, LCR filtering) via `mimicDetector/mimicDetector.smk`. Real runtime notes:
- The two POPScomp stages (pathogen + host structures) are strictly sequential (one structure at a time), not parallelized despite `threads: 8` being declared, and took several hours each for ~13k and ~97k structures respectively.
- The QSASA-filtering step is memory-intensive (Section 1, Section 15).
- Run detached for anything beyond a few minutes: `setsid nohup bash 04_run_mimicDetector.sh > ../logs/04_overall.log 2>&1 < /dev/null &` — see Section 15 for why plain `nohup`/`disown` were insufficient on WSL2.

Real final counts from this run: 4,109 raw hits → 289 QSASA survivors → 7 LCR survivors → 0 final candidates.

## 9. Statistical analysis: FDR estimation

Build the decoy database:

```bash
python3 05_prepare_decoy.py data/mimic_input/host/arachis_hypogaea.fasta \
  --out mimicDetector_results/aspergillus_flavus_run/arachis_hypogaea_decoy.fasta --seed 42
```

BLASTP the pathogen k-mers against the decoy (same parameters as the real host/control BLASTP runs in `mimicDetector.smk`):

```bash
makeblastdb -in mimicDetector_results/aspergillus_flavus_run/arachis_hypogaea_decoy.fasta -dbtype prot -parse_seqids

blastp -db mimicDetector_results/aspergillus_flavus_run/arachis_hypogaea_decoy.fasta \
  -query mimicDetector_results/aspergillus_flavus_run/aspergillus_flavus/aspergillus_flavus-12mers.fasta \
  -num_threads 8 -outfmt 6 \
  -out mimicDetector_results/aspergillus_flavus_run/aspergillus_flavus/aspergillus_flavus.arachis_hypogaea_decoy.12mers_blastp.out \
  -evalue 1 -comp_based_stats 0 -matrix PAM30 -word_size 2 -ungapped -dbsize 42904241
```

Build the decoy-equivalent filtered hit table (replicates the same bits_diff/E-value filter used for the real target hits):

```bash
python3 build_decoy_hcfiltered.py \
  --decoy-blast mimicDetector_results/aspergillus_flavus_run/aspergillus_flavus/aspergillus_flavus.arachis_hypogaea_decoy.12mers_blastp.out \
  --control-blast mimicDetector_results/aspergillus_flavus_run/aspergillus_flavus/aspergillus_flavus.aspergillus_oryzae.12mers_blastp.out \
  --out mimicDetector_results/aspergillus_flavus_run/aspergillus_flavus/aspergillus_flavus.arachis_hypogaea_decoy.12mers_b1_e001_hcfiltered_blastp.out \
  --bit-diff 1 --min-e 0.001
```

Calculate FDR (use `06_calculate_fdr_v2.py`, not the original `06_calculate_fdr.py` — see Section 16, item 7):

```bash
python3 06_calculate_fdr_v2.py \
  --target-hits mimicDetector_results/aspergillus_flavus_run/aspergillus_flavus/aspergillus_flavus.aspergillus_flavus.12mers_b1_e001_hcfiltered_blastp.out \
  --decoy-hits mimicDetector_results/aspergillus_flavus_run/aspergillus_flavus/aspergillus_flavus.arachis_hypogaea_decoy.12mers_b1_e001_hcfiltered_blastp.out \
  --out mimicDetector_results/aspergillus_flavus_run/aspergillus_flavus/fdr_result_v2.tsv
```

Real result: Ntarget=4,109, Ndecoy=359, UPtarget(by subject)=1,104, UPdecoy(by subject)=172, c=6.42, **FDR ≈ 0.56**.

## 10-11. Intermediate and final files

See README.md, "Key output files" table.

## 12. Visualization

`matplotlib` is required but not in `environment/environment.yml` (a real gap — add it: `conda install -n mimics -c conda-forge matplotlib`, then pin the installed version in `environment.yml`/`environment/software_versions.txt`; not yet recorded).

Figure 1 (workflow diagram) does not depend on real numbers and exists at `results/figures/figure1_workflow.png`.

**Figure 2 (bitscore distributions): attempted, not completed.** The command below is correct and uses the real file (note: the "host" BLASTP output is misleadingly named `aspergillus_flavus.aspergillus_flavus.12mers_blastp.out`, not `...arachis_hypogaea...` — the pipeline's `id` wildcard happens to equal the pathogen name in this single-pathogen run):

```bash
python3 08_generate_figures.py --outdir results/figures \
  --host-scores mimicDetector_results/aspergillus_flavus_run/aspergillus_flavus/aspergillus_flavus.aspergillus_flavus.12mers_blastp.out \
  --control-scores mimicDetector_results/aspergillus_flavus_run/aspergillus_flavus/aspergillus_flavus.aspergillus_oryzae.12mers_blastp.out
```

Required patch first (Section 16, item 8) — the original code loads the entire ~1GB host file into memory and attempts to parse every column as a float, not just bitscore; both a correctness bug and a memory risk.

Even patched, this step could not be completed in this environment: the host machine (7.7GB total RAM, ~3.9GB allocated to the WSL2 VM) repeatedly became unstable while processing the 1GB host file — confirmed via Windows Event Viewer to involve full machine restarts (`Kernel-Power`/`USER32` initiated, not a WSL-internal event), not resolved by `.wslconfig` idle-timeout settings or added swap. Suspected but not confirmed cause: system-wide memory pressure from the WSL VM's allocation competing with Windows itself on a small-RAM machine. **Deferred as a follow-up** — recommend running this one step on a machine with more available RAM (any machine with more than ~8GB free should complete it in seconds) rather than continuing to push this environment. No change to the scientific result — Figure 2 is a supplementary visualization, not part of the underlying analysis.

Figure 3 (stage-survival counts) has not yet been attempted; unlike Figure 2 it should be lightweight (a handful of integers, not a 1GB file) and is a reasonable next step:

```bash
python3 08_generate_figures.py --outdir results/figures --stage-counts <path-per-script's-actual-argument-format-check-with---help>
```

(Check `python3 08_generate_figures.py --help` for the exact expected input format for `--stage-counts` before running — not yet verified against the real script.)

Figure 4 (FDR vs. threshold sweep) requires re-running the FDR calculation at multiple bitscore-difference thresholds, which has not been done in this analysis.

## 13. Validation / quality checks

**Not yet run** (`07_validate_results.py`). Should be run against the real final candidate table before any downstream use, even though it has 0 rows (confirms no malformed/untraceable data slipped through undetected):

```bash
python3 07_validate_results.py \
  --parasite-fasta data/raw/parasite/aspergillus_flavus.UP000001875.fasta \
  --host-fasta data/raw/host/arachis_hypogaea.UP000289738.fasta \
  --control-fasta data/raw/control/aspergillus_oryzae.UP000006564.fasta \
  --final-candidates-tsv mimicDetector_results/aspergillus_flavus_run/aspergillus_flavus/aspergillus_flavus.aspergillus_flavus.12mers_b1_e001_q75_l50_paired_mimics.tsv
```

## 14. Expected outputs

| Stage | File | Real row count |
|---|---|---|
| Bitscore-diff/E-value filter | `*_hcfiltered_blastp.out` | 4,109 |
| QSASA filter | `*_q75_filtered.tsv` | 289 |
| Final (LCR + E-value re-scoring) | `*_q75_l50_paired_mimics.tsv` | 0 (header only) |
| FDR estimate | `fdr_result_v2.tsv` | FDR ≈ 0.56 |

## 15. Troubleshooting

Issues actually encountered during this analysis, in case they recur:

- **WSL2 VM disappearing mid-run, killing detached jobs.** `nohup`/`disown`/`setsid` all protect against your *terminal* closing, but WSL2 can tear down the entire underlying VM after an idle period, killing everything regardless. Fix: create `C:\Users\<you>\.wslconfig` with:
  ```
  [wsl2]
  vmIdleTimeout=-1
  ```
  **Watch for leading whitespace** if editing in Notepad via copy-paste — a leading space before `[wsl2]` silently breaks parsing, with no error, and the setting is simply ignored. Verify with `Get-Content C:\Users\<you>\.wslconfig` from PowerShell before trusting it.
- **`ls *.out | wc -l` reporting 0 files in a directory that actually has tens of thousands.** Large directories can exceed the shell's argument-list expansion limit, causing a silent empty result. Use `find <dir> -type f | wc -l` instead for any directory expected to hold more than a few thousand files.
- **Stale Snakemake lock after an ungraceful kill** (VM restart, segfault, Ctrl+C mid-run): `LockException: Directory cannot be locked`. Fix: `snakemake -s mimicDetector.smk --configfile <config> --unlock`, then re-run normally.
- **Memory pressure during QSASA filtering** (Section 1): add swap if the machine has less than ~8GB RAM: `sudo fallocate -l 6G /swapfile2 && sudo chmod 600 /swapfile2 && sudo mkswap /swapfile2 && sudo swapon /swapfile2`. Note swap does not survive a VM restart and must be re-enabled (`sudo swapon /swapfile2`) after one.
- **`Segmentation fault (core dumped)` from Snakemake itself** (not a rule's script) immediately after a fresh WSL VM restart, while building the DAG: encountered once, did not recur after applying the `vmIdleTimeout` fix and re-running; may have been coincidental to the VM's cold-start state, not a Snakemake defect.

## 16. Reproducibility notes: required patches to the vendored mimicDetector code

The following bugs were found in `mimicDetector` commit `c540df7` during this analysis and patched (patch scripts in `scripts/patch_*.py`). **All are required for a correct run**, not optional cleanup — several will silently produce wrong results (not crashes) if left unpatched.

1. **`patch_popscomp_gate.py`** (`mimicDetector.smk`) — the original popscomp completion-flag logic required every structure to succeed before writing the flag, meaning a single structure with a POPScomp geometry failure would permanently block the pipeline. Relaxed to write the flag after attempting all structures, logging failures to `popscomp_failed_structures.txt`.
2. **`patch_popscomp_shell.py`** (`mimicDetector.smk`) — a structure that makes `pops` exit with a nonzero code (not just silently produce no output) crashed the entire rule via Snakemake's `shell()` helper. Switched to `subprocess.run()`, which does not raise on nonzero exit.
3. **`patch_qsasa_hfile_match.py`** (`psychoscope.py`) — a substring match (`accession in filename`) could false-positive-match a different, longer accession's filename (e.g. `A0A445DSU1` matching `A0A445DSU12`). Fixed to an exact basename match.
4. **`patch_qsasa_pfile_missing.py`** (`psychoscope.py`) — the pathogen-side pops-file read did not catch `FileNotFoundError` for accessions whose popscomp run failed internally; added, symmetric with the host-side fix.
5. **`patch_hpops_output_dir.py`** (`mimicDetector.smk`) — **the most significant fix.** `get_hpops_files()` built the host popscomp *output* path using the literal string `"host"`, but the real output directory (per the actual `run_host_popscomp` rule) is named after `config["host_name"]`. Left unpatched, **every host QSASA lookup silently fails and defaults to "NaN" for the entire run**, regardless of whether that structure's popscomp run actually succeeded, invalidating QSASA filtering completely without any error being raised.
6. **`patch_qsasa_h_counter_init.py`** (`psychoscope.py`) — a logging counter variable (`h`) is incremented but never initialized, causing `UnboundLocalError`. Only reachable once fix #5 is applied and real host data starts flowing through this code path.
7. **`patch_lcr_empty_input.py`** (`psychoscope.py`) — reading a genuinely empty QSASA-filtered input (a valid outcome, not an error) crashes on a dtype mismatch in a later arithmetic step. Added an early empty-input check with graceful header-only output.

Apply in the order listed (patches 1-2 to `mimicDetector.smk`; patches 3, 4, 6, 7 to `psychoscope.py`; patch 5 to `mimicDetector.smk`):

```bash
cd scripts
python3 patch_popscomp_gate.py ../mimicDetector/mimicDetector.smk
python3 patch_popscomp_shell.py ../mimicDetector/mimicDetector.smk
python3 patch_hpops_output_dir.py ../mimicDetector/mimicDetector.smk
python3 patch_qsasa_hfile_match.py ../mimicDetector/scripts/psychoscope.py
python3 patch_qsasa_pfile_missing.py ../mimicDetector/scripts/psychoscope.py
python3 patch_qsasa_h_counter_init.py ../mimicDetector/scripts/psychoscope.py
python3 patch_lcr_empty_input.py ../mimicDetector/scripts/psychoscope.py
```

Additionally, **use `06_calculate_fdr_v2.py`, not `06_calculate_fdr.py`**, for FDR calculation — the original's "unique peptide" definition (unique query k-mer) is mathematically degenerate against this pipeline's already-deduplicated hit tables and always returns FDR=1.0 regardless of the real data (see `report/final_report.md`, Section 7, for the full explanation).

8. **`patch_fig2_bitscore_parsing.py`** (`08_generate_figures.py`, this project's own script, not vendored) — `fig2_score_distributions()` loaded an entire BLASTP outfmt6 file into memory via `.read_text().split()` and attempted to parse every whitespace-separated token as a float, including the non-numeric `qseqid`/`sseqid` ID columns. Confirmed to cause an out-of-memory kill on the real ~1GB host-vs-pathogen file before it could even reach the point of raising a `ValueError` on a non-numeric column. Fixed to stream the file line-by-line and parse only the last (bitscore) column.

```bash
python3 patch_fig2_bitscore_parsing.py ../scripts/08_generate_figures.py
```

The decoy database uses a fixed random seed (42) for reproducibility (`05_prepare_decoy.py`).
