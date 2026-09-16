# Candidate Molecular Mimicry Between *Aspergillus flavus* and Peanut (*Arachis hypogaea*) Host Proteins: A mimicDetector Analysis

## Abstract
This analysis applied mimicDetector (Rich & Wasmuth 2026) to a species pair not examined in the original paper: *Aspergillus flavus* (fungal pathogen) vs. *Arachis hypogaea* / peanut (host), with *Aspergillus oryzae* (domesticated, non-pathogenic close relative) as negative control. Method: k=12 BLASTP/PAM30/wordsize=2 k-mer search, E-value <=0.001, host-control bitscore difference >=1, QSASA >75%, low-complexity overlap <=50%, final PAM30 re-alignment E-value <=0.001. FDR: target-decoy against a shuffled peanut proteome, per Kall et al. 2008-style target-decoy correction.

**Result: no high-confidence candidates.** Of 4,109 raw hits surviving the initial bitscore-difference/E-value filter, 289 survived QSASA filtering (>0.75 on both pathogen and host sides), 7 survived low-complexity-region filtering, and 0 survived the final statistical re-scoring step (closest miss: E=0.00153 against a 0.001 cutoff). A corrected target-decoy FDR estimate for the initial filtering stage was ~56%, indicating substantial expected noise even before structural filtering. This is a genuine negative result for this species pair and control, not an artifact of overly strict thresholds (see Section 8-9).

## 1. Introduction
Molecular mimicry — where pathogen molecules structurally or functionally resemble host molecules — is best documented in animal/human host-pathogen systems, but the same principle applies in principle to plant pathogens, which must also evade or subvert host immune signaling (e.g. PAMP-triggered immunity, effector-triggered immunity) to establish infection. *Aspergillus flavus* is of particular interest as a case study because it is simultaneously (i) a major agricultural pathogen of peanut, causing direct yield loss and carcinogenic aflatoxin contamination of a globally important food crop, and (ii) a species with an unusually well-matched, real, non-pathogenic close relative (*A. oryzae*) available as a negative control. Applying the same pipeline, parameters, and thresholds used in the original human-targeted paper to this new, taxonomically very different host (a plant) is itself informative: any mimicry signal found (or not found) should be interpreted in light of the fact that the FDR benchmarking underlying the paper's recommended thresholds was performed only on animal-host systems (Section 9).

## 2. Biological system
- **Pathogen**: *Aspergillus flavus* (UniProt UP000001875, strain ATCC 200026/FGSC A1120/NRRL 3357) — causes yellow mold disease and aflatoxin contamination of peanut and maize.
- **Host**: *Arachis hypogaea* / peanut (UniProt UP000289738, cv. Tifrunner) — allotetraploid (AABB genome).
- **Negative control**: *Aspergillus oryzae* (UniProt UP000006564, strain RIB40) — domesticated industrial "koji mold" (sake, soy sauce, miso fermentation), not a plant pathogen, reported ~99.5% gene homology to *A. flavus* NRRL3357. Selected because it is phylogenetically very close to the pathogen while lacking its plant-pathogenic and aflatoxigenic phenotype, matching the pipeline authors' own data-selection guidance as closely as any real, available proteome allows.

## 3. Data
See `data/metadata/proteome_manifest.tsv` for full accession, release, and URL records.
- Pathogen proteome: 13,729 sequences downloaded (13,318 with AlphaFold structure coverage, 97.0%).
- Host proteome: 97,687 sequences downloaded (95,753 with AlphaFold structure coverage, 98.0%).
- The ~2-3% coverage gap in both proteomes is entirely accounted for by AlphaFold's documented exclusion of sequences >2,700 aa from proteome-scale predictions — not a pipeline error.

## 4. Computational environment
See `environment/environment.yml` and `environment/software_versions.txt`.

## 5. Methods
Pathogen proteins fragmented into overlapping 12-aa k-mers; BLASTP (PAM30, wordsize 2) against host and control proteomes; retain hits with E-value <=0.001 AND bitscore(host) - bitscore(control) >=1; merge overlapping k-mer alignments into preliminary candidates; retain candidates with mean QSASA (POPSCOMP) >75% on both pathogen and host sides; remove candidates with >50% low-complexity overlap (Segmasker intervals); final candidates re-scored via independent PAM30 local realignment (E-value <=0.001). No full-length homologue-removal step is applied (paper's updated/retained-homologue methodology), applied here without modification to a plant host.

## 6. Negative control
The *A. oryzae* proteome is searched under identical BLASTP parameters to the host search. A candidate is only retained if its pathogen-host bitscore exceeds its best pathogen-control bitscore by >=1 bit — this is the mechanism by which the negative control functions, and is distinct from the FDR decoy in Section 7. Because *A. oryzae* is unusually close to *A. flavus* (Section 2), this control is expected to behave more stringently than the controls used in the original paper's benchmark species.

## 7. FDR analysis
Target = real peanut (host) proteome hits at the bitscore-difference/E-value stage; decoy = per-sequence shuffled peanut proteome (seed=42, `scripts/05_prepare_decoy.py`), run through the identical filter (`scripts/build_decoy_hcfiltered.py`).

`FDR = c * (Ndecoy/Ntarget)`, `c = UPtarget/UPdecoy`

**Note on "unique peptide" (UP):** the paper's Methods do not define this term precisely enough to resolve directly. An initial implementation counting unique *query* k-mers (`scripts/06_calculate_fdr.py`) is mathematically degenerate for this pipeline's hit tables, since both target and decoy hit tables are already deduplicated to one row per query upstream (top-bitscore-per-query selection) — this makes UP identically equal to N by construction, forcing FDR=1.0 regardless of the real data. `scripts/06_calculate_fdr_v2.py` corrects this by counting unique **subject** (host/decoy protein) hits instead, which is not deduplicated upstream and can genuinely differ between target and decoy hit patterns. This is a documented methodological choice, not an assumption silently carried over from the paper.

**Real result:**

| Metric | Target | Decoy |
|---|---|---|
| N (hits, bits_diff/E-value filtered) | 4,109 | 359 |
| UP (unique subject proteins hit) | 1,104 | 172 |

c = 6.42, **FDR ≈ 0.56 (56%)** at the initial sequence-similarity stage.

**Scope of this estimate:** it applies only to the first-stage sequence-similarity filter. It cannot be extended to validate the downstream QSASA filtering step, since QSASA requires a real folded 3D structure (AlphaFold + POPScomp), and shuffled decoy sequences have no biologically real structure to fold. This is an unresolved gap in what target-decoy FDR can validate for this pipeline design, not an oversight to be silently assumed away.

## 8. Results
- Raw BLASTP hits surviving bits_diff/E-value filter: **4,109**
- Surviving QSASA filtering (>0.75 both sides): **289**
- Surviving low-complexity-region filtering (<50% overlap): **7**
- Surviving final E-value re-scoring (<=0.001): **0**

The 7 LCR-survivors were borderline misses, not far outliers: E-values ranged 0.00153-0.00654 against the 0.001 cutoff (closest: *A0A7U2MNP9* vs *A0A445BFV6*, bitscore 38.0, E=0.00153).

## 9. Biological interpretation
No candidates survive the full pipeline at the stated thresholds. Combined with the ~56% FDR at the initial filtering stage, this indicates the raw candidate pool is already substantially noise-dominated before structural filtering is even applied — the zero final count is not simply an artifact of an overly conservative downstream threshold. Given this, no biological interpretation of specific candidate proteins is offered here; there are none that meet the pipeline's own significance criteria.

## 10. Controls and reproducibility
Negative control: Section 6. FDR: Section 7. Full environment and command provenance: `environment/software_versions.txt`. See `FOLLOWME.md` for the complete, step-by-step reproduction procedure, including the code fixes required to run this pipeline correctly (Section 16).

## 11. Limitations
- Sequence similarity detected by this pipeline is NOT evidence of functional mimicry, structural mimicry, or biochemical interaction — only of sequence-level candidate mimicry meeting the stated statistical criteria.
- **Cross-kingdom threshold transfer**: the paper's recommended filtering thresholds (b=1, e=0.001, q=0.75) were benchmarked only against a human host using three animal-infecting pathogens. Their suitability for a plant host has not been independently established by the paper.
- **Control stringency**: *A. oryzae*'s very close relationship to *A. flavus* may make the negative-control filter unusually conservative, potentially suppressing true candidates rather than only false ones. A less-stringent control was considered as a sensitivity check but was not run in this analysis due to the computational cost of repeating the structural (POPScomp) stage against a new control species.
- **FDR scope**: as noted in Section 7, the target-decoy FDR estimate covers only the initial sequence-similarity stage and cannot validate the downstream QSASA filtering step.
- Proteome completeness: the peanut reference proteome, being a recently sequenced allotetraploid crop genome, may have more incomplete or ambiguous gene models than a long-curated model organism proteome. A ~244-accession discrepancy between the downloaded pathogen proteome count and the official UniProt count was not resolved (ruled out: duplicates, fragments, isoform splitting).
- Annotation limitations: many *A. flavus* and peanut proteins are poorly annotated, which would limit biological interpretation of any candidates, had any survived filtering.
- Decoy assumption: the per-sequence shuffling algorithm used here (`scripts/05_prepare_decoy.py`) is a documented, reasonable interpretation of the paper's "shuffled [host] sequences" but is not verified to be identical to the original authors' (human-proteome-specific) implementation.
- Potential evolutionary conservation: had candidates survived filtering, some could reflect ancient shared ancestry/domain conservation (e.g., core metabolic enzymes conserved across the fungal-plant divide) rather than convergent mimicry — a general risk of this method's no-homologue-removal design, not specific to this null result.
- Structural coverage gaps: ~2-3% of both proteomes lack AlphaFold structure coverage (Section 3); a small number of additional structures (2 confirmed pathogen-side, an unquantified small number host-side) failed POPScomp's internal geometry checks and were excluded from QSASA scoring (`popscomp_failed_structures.txt` in each species' output directory).

## 12. Conclusion
This analysis does not support a high-confidence molecular mimicry hit between *A. flavus* and peanut, given *A. oryzae* as the bitscore-difference control. Zero of 4,109 initial candidates survive the complete filtering pipeline, and a corrected target-decoy FDR estimate of ~56% at the initial stage indicates this is consistent with the candidate pool being substantially noise-dominated from the outset, not merely a product of downstream threshold strictness. This is a genuine negative result for this specific species pair, control choice, and threshold set — not a claim that no molecular mimicry exists between *A. flavus* and peanut under any possible analysis design.
