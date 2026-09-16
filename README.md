# *Aspergillus flavus* → Peanut Molecular Mimicry Detection

Applies the [mimicDetector](https://github.com/Kayleerich/mimicDetector) pipeline (Rich & Wasmuth, *Bioinformatics* 42(2):btag012, 2026) to a species pair not examined in the original paper, to test for candidate molecular mimicry between a fungal pathogen and a plant host.

## Research objective

mimicDetector's original validation surveyed 17 pathogens against the human proteome only. This project applies the same pipeline, parameters, and thresholds to a genuinely new system: a fungal crop pathogen and its plant host, to test whether the method's animal/human-calibrated thresholds transfer to a cross-kingdom (fungus-vs-plant) comparison, and whether any candidate mimicry signal exists between this specific pathogen-host pair.

## Biological system

| Role | Organism | Accession | Why |
|---|---|---|---|
| Pathogen | *Aspergillus flavus* (strain NRRL 3357) | [UniProt UP000001875](https://www.uniprot.org/proteomes/UP000001875) | Major pre- and post-harvest pathogen of peanut, causes yellow mold disease and aflatoxin (carcinogenic) contamination |
| Host | *Arachis hypogaea* / peanut (cv. Tifrunner) | [UniProt UP000289738](https://www.uniprot.org/proteomes/UP000289738) | Allotetraploid crop host |
| Negative control | *Aspergillus oryzae* (strain RIB40) | [UniProt UP000006564](https://www.uniprot.org/proteomes/UP000006564) | Domesticated, non-pathogenic close relative (~99.5% homology to *A. flavus*); used for the pipeline's bitscore-difference filter |

Full provenance, including URLs and selection rationale, in [data/metadata/proteome_manifest.tsv](data/metadata/proteome_manifest.tsv).

## Analytical workflow

1. Download and QC the three proteomes (pathogen, host, control).
2. Fetch AlphaFold structures for all three proteomes.
3. Fragment the pathogen proteome into overlapping 12-aa k-mers; BLASTP each against host and control.
4. Retain hits where E-value ≤0.001 and bitscore(host) − bitscore(control) ≥1.
5. Score retained candidates by mean solvent accessibility (QSASA, via POPScomp) on both pathogen and host sides; retain >0.75 on both.
6. Remove candidates with >50% overlap with a low-complexity/repetitive region (Segmasker).
7. Re-score final candidates with an independent local realignment; retain E-value ≤0.001.
8. Estimate false discovery rate via a shuffled-peanut decoy database, using the target-decoy method.

Full step-by-step commands: **[FOLLOWME](FOLLOWME.md)**. Full results and their interpretation: **[Report](report/final_report.md)**.

## Software and versions

Pinned in `environment/environment.yml` (conda env `mimics`): Python 3.11.12, Snakemake 7.32.4, BLAST+ 2.16.0 (includes Segmasker), BEDTools 2.31.1, Biopython 1.85, pandas 2.2.3, pybedtools 0.12.0.

**Not included in the conda environment** — must be installed separately: [POPScomp](https://github.com/Fraternalilab/POPScomp) (structural solvent-accessibility calculation).

Full environment record, including real execution dates and known compute constraints: `environment/software_versions.txt`.

## Repository organization

```
.
├── README.md                    This file
├── FOLLOWME.md                  Step-by-step reproduction guide
├── COPYING.md                   License
├── mimicDetector/                Vendored mimicDetector tool (git clone,
│                                  commit c540df7). This is the ONLY copy
│                                  of the pipeline code that is actually
│                                  executed — scripts/04_run_mimicDetector.sh
│                                  invokes mimicDetector/mimicDetector.smk
│                                  directly. Includes this project's
│                                  fixes to real bugs found during
│                                  execution (see FOLLOWME.md, Section 16).
├── scripts/                      This project's own analysis scripts
│   ├── 00-04                     Data acquisition, QC, and pipeline execution
│   ├── 05-08                     Decoy generation, FDR calculation,
│   │                              validation, and figure generation
│   ├── build_decoy_hcfiltered.py Builds the decoy-equivalent filtered
│   │                              hit table needed for FDR calculation
│   ├── diagnose_lcr_evalue.py     Read-only diagnostic: re-runs the
│   │                              LCR-overlap and final E-value scoring
│   │                              steps standalone, for inspecting
│   │                              borderline candidates
│   └── patch_*.py                Applied, one-off patches to the vendored
│                                  mimicDetector code (see FOLLOWME.md,
│                                  Section 16, for what each one fixes)
├── data/
│   ├── metadata/                  Proteome provenance, accession lists,
│   │                               structure-coverage records
│   ├── mimic_input/                BLAST databases built from the
│   │                               downloaded proteomes
│   └── processed/                  QC output tables
├── environment/                   Environment/version/provenance records
├── mimicDetector_results/
│   ├── aspergillus_flavus_run/     Real analysis output (see below)
│   └── test/                       Output of a smoke-test run against
│                                   mimicDetector's own bundled test
│                                   fixtures (mimicDetector/test_data/),
│                                   used to confirm the pipeline and
│                                   environment were correctly configured
│                                   before committing to the real run.
├── report/final_report.md         Full results, methodology, and limitations
├── results/                        Figures and a toy formula-verification
│                                   example (not real biological results)
└── logs/                           Execution logs
```

## Key output files

Under `mimicDetector_results/aspergillus_flavus_run/aspergillus_flavus/`:

| File | Content |
|---|---|
| `*_hcfiltered_blastp.out` | Hits surviving the bitscore-difference/E-value filter (4,109 rows) |
| `*_qsasa_averages.tsv` | Mean QSASA per candidate, pre-threshold |
| `*_q75_filtered.tsv` | Candidates surviving QSASA >0.75 on both sides (289 rows) |
| `*_q75_l50_paired_mimics.tsv` | Final candidates after LCR filtering + E-value re-scoring (**0 rows** — see below) |
| `fdr_result_v2.tsv` | Corrected target-decoy FDR estimate (~56% at the initial filtering stage) |

## Main finding

**No candidates survive the complete pipeline.** Of 4,109 initial hits, 289 survive QSASA filtering, 7 survive low-complexity filtering, and 0 survive final statistical re-scoring (closest miss: E=0.00153 against a 0.001 cutoff). A corrected FDR estimate of ~56% at the initial filtering stage indicates the candidate pool is substantially noise-dominated before structural filtering is even applied. This is a genuine negative result for this species pair and control choice — full interpretation and caveats in `report/final_report.md`.

## Important limitations (see `report/final_report.md`, Section 11, for full list)

- The paper's filtering thresholds were calibrated on animal/human host-pathogen systems only; their transfer to a fungus-vs-plant comparison is not independently validated beyond this project's own FDR check.
- The target-decoy FDR estimate covers only the initial sequence-similarity stage; it cannot be extended to validate the QSASA (structural) filtering step, since shuffled decoy sequences have no real 3D structure to score.
- *A. oryzae*, while satisfying the pipeline's control-selection criteria, is unusually close to *A. flavus* (~99.5% homology), likely making this control more stringent than those used in the original paper.

## Reproducibility

See **`FOLLOWME.md`** for the complete, tested, step-by-step procedure, including the real bugs found in the vendored mimicDetector code and the exact fixes applied (necessary for the pipeline to run to completion — the unpatched upstream code will not produce correct QSASA values for the host side, among other issues).

## Citation

This project uses the mimicDetector pipeline. If reporting results from this analysis, cite:

> Rich, K. & Wasmuth, J. (2026). mimicDetector: a pipeline for protein motif mimicry detection in host-pathogen interactions. *Bioinformatics*, 42(2), btag012. https://academic.oup.com/bioinformatics/article/42/2/btag012/8423033
