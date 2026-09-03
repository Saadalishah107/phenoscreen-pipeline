# PhenoScreen: Genotype-Guided Natural Product Virtual Screening Pipeline

**A containerized, Nextflow-orchestrated workflow that bridges statistical genomics, variant annotation, structural biology, and deep learning-based drug-target interaction screening.**

[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A526.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/Docker-Containerized-blue.svg)](https://www.docker.com/)
[![PyTorch](https://img.shields.io/badge/PyTorch-DeepPurpose-EE4C2C.svg)](https://pytorch.org/)

---

## Project Overview

PhenoScreen accepts raw, multi-sample genomic variant calls and a phenotype of interest, and processes them through four decoupled, containerized stages to produce a ranked list of natural-product compounds predicted to bind the proteins implicated by those variants. The pipeline is implemented as a Directed Acyclic Graph (DAG) using Nextflow DSL2, with each stage isolated in its own Docker image so that the PLINK/R, Biopython/REST API, and RDKit/PyTorch dependency trees do not conflict with one another.

The pipeline is designed to prioritize scientific honesty over cosmetic completeness. Every candidate locus is explicitly labeled as `significant_fdr` or `exploratory_top_n` rather than being silently combined into a single category. Every predicted structure carries a confidence and source field so that mock data and predicted data are never ambiguous. Every entry in the final results table is documented as a computational prediction rather than a validated experimental finding.

---

## Scientific Reasoning and Design Philosophy

This pipeline is best described as genotype-guided virtual screening rather than classical phenotypic drug discovery, since it begins from genetic variation rather than from a cellular phenotype assay and proceeds forward toward candidate compounds. Several design decisions follow directly from this framing.

- **Statistical honesty over convenience.** The association module never promotes a non-significant SNP to candidate status without clear labeling. Every locus recorded in `candidate_loci.tsv` carries a `selection_status` value of either `significant_fdr` (a locus that passed the FDR threshold) or `exploratory_top_n` (a locus that did not pass the threshold and is included only so that downstream stages have input during pilot runs).
- **Real biological grounding rather than placeholder translation.** Candidate loci are annotated through the Ensembl VEP REST API and filtered to protein-altering consequences only, specifically missense and stop-gained variants. Mutated protein sequences are reconstructed from the true reference sequence combined with the annotated amino acid change, rather than substituted with a generic placeholder.
- **Sequence-based drug-target interaction (DTI) prediction, labeled accurately.** DeepPurpose scores protein sequence and ligand SMILES pairs using the `Morgan_CNN_BindingDB` model. Predicted three-dimensional structures are retained as structural evidence, recorded in `structure_manifest.tsv` along with a mean Cα B-factor confidence proxy, but these structures are not passed into the DTI model itself. This pipeline does not perform structure-based docking, and the documentation does not claim otherwise.
- **Mock execution modes as a core feature.** The `--struct_backend mock` and `--dti_mode mock` options allow the entire seven-stage DAG to be exercised end to end without contacting external APIs or downloading multi-gigabyte PyTorch and DeepPurpose model weights. This capability is essential for fast continuous integration and for smoke testing within GitHub Codespaces.

---

## System Architecture and Orchestration

The pipeline entry point, `main.nf`, validates input parameters, constructs four input channels (VCF, phenotype, COCONUT compound library, and an optional annotation fixture), and delegates execution to a nested subworkflow defined in `workflows/phenoscreen.nf`, which wires together seven processes. The pipeline uses the Nextflow DSL2 `output {}` block, which requires Nextflow version 26.04.0 or later, to declare the final publish destination for every named output. This approach replaces per-process `publishDir` directives with a single, auditable mapping from each Nextflow channel name to its corresponding results subdirectory.

```
VCF -> LD_PRUNING -> PCA_ASSOCIATION -> ANNOTATE_VARIANTS -> PREDICT_STRUCTURE
                             |                                        |
                       candidate_loci.tsv                 structure_manifest.tsv
                                                                       |
COCONUT -> PROCESS_CHEM_LIBRARY ----------------------------> DTI_SCREENING -> RANK_HITS
                                                                                    |
                                                                 top_predicted_candidates.tsv
```

### Infrastructure: Three Phase-Isolated Images

| Image | Base | Purpose |
|---|---|---|
| `phenoscreen/genomics:2.0` | `condaforge/mambaforge` | PLINK 1.9/2, R 4.3 (SNPRelate, ggplot2, data.table) |
| `phenoscreen/structure:2.0` | `python:3.11-slim` | Biopython and `requests`, used for both variant annotation and structure retrieval |
| `phenoscreen/cheminformatics:2.0` (or `-test:2.0` for mock DTI) | `condaforge/mambaforge` | RDKit, PyTorch (CPU only, installed via pip to avoid MKL conflicts with conda-installed RDKit), DeepPurpose, and `descriptastorus` |

These images are not published to any container registry. They must be built locally before running the pipeline with `-profile docker`:

```bash
bash docker/build_images.sh          # builds genomics, structure, and cheminformatics-test
bash docker/build_images.sh --full   # also builds the full cheminformatics image, required for real --dti_mode deepurpose runs
```

---

## Detailed Module Specifications

### Module 1: Genomics (`LD_PRUNING` and `PCA_ASSOCIATION`)

**Objective:** Remove redundant, correlated SNPs, correct for population stratification, and identify variant-phenotype associations.

`LD_PRUNING` runs PLINK's `--indep-pairwise` procedure (default window, step, and r-squared values of `50 5 0.2`, configurable through `--indep_pairwise`) on the input VCF, and emits the pruned `.bed`, `.bim`, and `.fam` files together as a single tuple channel so that downstream steps always receive all three files as a unit.

`PCA_ASSOCIATION` computes principal components using SNPRelate, then performs a principal-component-adjusted regression for each SNP. A linear model is used when `--phenotype_type` is set to `continuous`, and a logistic regression model, implemented through `glm()` with a binomial family, is used when `--phenotype_type` is set to `binary`. Benjamini-Hochberg false discovery rate correction is applied using the `--fdr_threshold` parameter, which defaults to `0.05`.

Outputs include `candidate_loci.tsv`, `association_results.tsv`, PCA eigenvectors, a scree plot, a scatter plot, and Manhattan and quantile-quantile plots.

### Module 2: Variant Annotation (`ANNOTATE_VARIANTS`)

**Objective:** Resolve each candidate SNP to a gene, transcript, and reconstructed mutant protein sequence.

In live mode, the module queries the Ensembl VEP REST API using the endpoint `/vep/{species}/id/{rsid}`, filters transcript consequences to `missense_variant` and `stop_gained` categories, prefers the canonical transcript when multiple candidates exist, and retrieves the reference protein sequence through `/sequence/id/{protein_id}`. The annotated amino acid substitution or truncation is then applied to construct the variant-specific protein sequence.

In offline or fixture mode, enabled by setting `--offline_annotation true` or by supplying `--annotation_tsv <fixture.tsv>`, the module substitutes a pre-curated TSV file indexed by rsID. This mode is used by both the synthetic `test` profile and the curated NCI-60 integration fixture.

The number of live API calls is bounded by the `--max_annotation_variants` parameter during development and testing.

Outputs include `variant_annotation.tsv`, which contains gene, transcript, and protein identifiers, the annotated consequence, and both the reference and mutant protein sequences, as well as `protein_sequences.fasta`.

### Module 3: Structure Prediction (`PREDICT_STRUCTURE`)

**Objective:** Obtain a three-dimensional structural model for each annotated protein target without requiring local GPU inference.

Three backends are available through the `--struct_backend` parameter: `alphafold`, which retrieves a precomputed model from the AlphaFold Database by UniProt accession; `esmfold`, which queries the single-sequence, alignment-free ESMFold REST API; and `mock`, which generates deterministic placeholder structures for rapid continuous integration testing. Transient failures are retried with exponential backoff, configured through the `--max_retries` and `--retry_delay` parameters.

Outputs include the predicted structures in `structures/*.pdb` format, along with `structure_manifest.tsv`, which records the status, source, and mean Cα B-factor for each target. For structures generated by ESMFold or AlphaFold, this B-factor value also serves as a per-residue confidence proxy corresponding to predicted local distance difference test (pLDDT) scores, so structural confidence can be audited rather than assumed.

### Module 4: Cheminformatics and DTI Screening (`PROCESS_CHEM_LIBRARY`, `DTI_SCREENING`, and `RANK_HITS`)

`PROCESS_CHEM_LIBRARY` standardizes the input SMILES or SDF library using RDKit, computes Morgan fingerprints according to the `--fingerprint_radius` and `--fingerprint_bits` parameters, and generates an expanded set of two-dimensional molecular descriptors, including molecular weight, LogP, topological polar surface area, hydrogen bond donor and acceptor counts, rotatable bond count, aromatic ring count, heavy atom count, ring count, formal charge, quantitative estimate of drug-likeness (QED), and fraction of sp3-hybridized carbons.

`DTI_SCREENING` runs DeepPurpose, using the model specified by `--dti_model` (default `Morgan_CNN_BindingDB`), across every combination of target protein sequence and compound SMILES string, processed in batches according to `--dti_batch_size`. Setting `--dti_mode mock` substitutes a fast, deterministic scoring function for pipeline validation purposes. The `--max_compounds` parameter limits the size of the screening library during iterative testing.

`RANK_HITS` sorts the complete affinity matrix, joins the results with the corresponding molecular descriptors, and writes the final output table, `top_predicted_candidates.tsv`, containing the top candidates as specified by `--top_n` (default `50`), along with a summary visualization.

---

## Configuration and Environment Setup

### Host Requirements

- Nextflow version 26.04.0 or later, required for the DSL2 `output {}` publish block syntax used in `main.nf`
- Docker Engine, or alternatively Singularity or Conda, as described in the profiles section below

### Key Parameters (`nextflow.config`)

| Parameter | Default | Purpose |
|---|---|---|
| `--vcf`, `--phenotype`, `--coconut` | Required, no default | Core pipeline inputs |
| `--annotation_tsv` | `null` | Optional offline annotation fixture |
| `--phenotype_type` | `continuous` | Specifies `continuous` (linear model) or `binary` (logistic model) |
| `--fdr_threshold` | `0.05` | Benjamini-Hochberg adjusted significance cutoff |
| `--candidate_top_n` | `20` | Size of the exploratory fallback panel |
| `--max_annotation_variants` | `20` | Limit on live Ensembl VEP lookups |
| `--struct_backend` | `alphafold` | One of `alphafold`, `esmfold`, or `mock` |
| `--dti_mode` | `deepurpose` | One of `deepurpose` or `mock` |
| `--dti_model` | `Morgan_CNN_BindingDB` | DeepPurpose pretrained encoder pairing |
| `--max_compounds` | `0` (unlimited) | Limit on the size of the screening library |
| `--top_n` | `50` | Number of candidates included in the final ranked table |

### Execution Profiles

- **`test`**: Loads the synthetic fixtures located in `test_data/`, and sets `--offline_annotation true`, `--struct_backend mock`, and `--dti_mode mock`. This profile makes no external API calls and provides the fastest possible full-pipeline smoke test.
- **`real_small`**: A reduced-scale configuration, using fewer principal components, fewer annotation variants, and a capped compound count, intended for validating the live-data execution path without incurring the runtime cost of a full production run.
- **`docker`, `singularity`, and `conda`**: Select the container or environment backend used for execution.

---

## Execution Guide

### Option A: Synthetic Smoke Test

This is the recommended first run for verifying that the pipeline is correctly installed and configured.

```bash
bash docker/build_images.sh
nextflow run main.nf -profile test,docker
```

This command runs entirely offline against the fixtures in `test_data/`, and verifies the DAG structure, channel wiring, and container handoffs without contacting any external API.

### Option B: NCI-60 Real-Data Integration Test

This test uses a curated fixture located in `real_data/nci60/`, generated by the scripts `scripts/dataset_curation/curate_nci60.py` and `generate_fixture.py`. The fixture contains real oncogenic variants, including BRAF V600E, TP53 R273H, and KRAS G12S, along with actual Paclitaxel drug-sensitivity phenotype data.

```bash
nextflow run main.nf \
  -profile docker \
  --vcf real_data/nci60/nci60_mutations.vcf \
  --phenotype real_data/nci60/phenotype.tsv \
  --coconut test_data/sample_coconut.smi \
  --annotation_tsv real_data/nci60/nci60_annotation.tsv \
  --phenotype_type continuous \
  --struct_backend alphafold \
  --dti_mode deepurpose \
  --candidate_top_n 30 \
  --max_annotation_variants 30 \
  --max_compounds 8 \
  --top_n 10 \
  --outdir results_nci60_production
```

**Note regarding the compound library used in this run.** The command above uses the small, synthetic `test_data/sample_coconut.smi` file rather than a full COCONUT database export. This configuration was chosen specifically to validate that the annotation, structure prediction, and DeepPurpose inference stages function correctly on real protein targets, and it does not represent a production-scale screen. A complete COCONUT SDF or SMILES export should be substituted before the contents of `top_predicted_candidates.tsv` are treated as an actual screening result.

---

## Current Limitations

- **Statistical power.** The NCI-60 integration fixture, with a sample size of ten, functions as a computational integration test rather than a statistically powered association study.
- **Compound library scale.** Production runs to date have relied on the small synthetic SMILES set to validate the architecture on real protein targets, rather than on a complete COCONUT export. Results from these runs should be interpreted as validation of the computational architecture, not as the output of a genuine screening effort.
- **Somatic versus germline variation.** The genomics module, built on PLINK and SNPRelate, is designed for germline association testing. Adapting these methods for somatic tumor mutation burden would be necessary before applying this pipeline to true oncology cohorts.
- **Sequence-based rather than structure-based DTI prediction.** Predicted protein structures are retained as supporting evidence and quality control, including a confidence proxy, but they are not provided as input to the DeepPurpose model. This pipeline does not perform molecular docking.
- **Compute environment.** The cheminformatics container runs PyTorch in CPU-only mode. Large-scale screening or GPU-accelerated inference would require substituting a CUDA-enabled image.

---

## Repository Layout

```
.
├── main.nf                          # Entry point: parameter validation, channels, and output {} publish map
├── workflows/phenoscreen.nf         # Subworkflow wiring all seven processes
├── nextflow.config                  # Parameters, profiles, and per-label container and resource mapping
├── modules/
│   ├── genomics/                    # LD_PRUNING, PCA_ASSOCIATION
│   ├── annotation/                  # ANNOTATE_VARIANTS
│   ├── structure/                   # PREDICT_STRUCTURE
│   └── cheminformatics/             # PROCESS_CHEM_LIBRARY, DTI_SCREENING, RANK_HITS
├── bin/                             # Python and R scripts invoked by each process
├── docker/                          # Per-phase Dockerfiles and build_images.sh
├── real_data/nci60/                 # Curated real-variant integration fixture
├── scripts/dataset_curation/        # Fixture generation scripts (curate_nci60.py, generate_fixture.py)
├── test_data/                       # Synthetic and offline annotation fixtures for the test profile
├── docs/                            # Architecture notes and example result artifacts
└── CITATION.cff, LICENSE, Makefile, CODESPACES.md
```

---

## References

- Cock, P. J., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics.
- Huang, K., et al. (2020). DeepPurpose: a deep learning library for drug-target interaction prediction.
- Landrum, G., et al. RDKit: Open-source cheminformatics software.
- Zheng, X., et al. (2012). A high-performance computing toolset for relatedness and principal component analysis of SNP data (SNPRelate).
- Purcell, S., et al. (2007). PLINK: a tool set for whole-genome association and population-based linkage analyses.
