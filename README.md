# PhenoScreen

**An End-to-End Containerized Pharmacogenomics and Drug-Target Interaction Pipeline**

[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A522.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/Docker-Containerized-blue.svg)](https://www.docker.com/)
[![PyTorch](https://img.shields.io/badge/PyTorch-DeepPurpose-EE4C2C.svg)](https://pytorch.org/)

## 1. Overview & Research Scope

PhenoScreen is a reproducible, HPC-ready bioinformatics workflow designed to bridge population genomics, structural biology, and deep learning-based cheminformatics. Orchestrated via **Nextflow** and fully isolated within **Docker** containers, the pipeline automates the discovery of drug-target interactions directly from raw genomic variation.

It ingests multi-sample genomic variants and clinical phenotypes, identifies statistically significant loci, maps them to 3D protein structures, and executes neural network-based screening against natural compound libraries. This repository serves as a robust proof-of-concept for scalable, multi-omics machine learning architectures in computational oncology and pharmacogenomics, engineered specifically to demonstrate production-grade software development for high-performance computing (HPC) environments.

## 2. Pipeline Architecture

The workflow is constructed as a Directed Acyclic Graph (DAG) divided into four highly modular, containerized execution environments:

- **Module 1: Genomics (PLINK & R)**
  - Performs linkage disequilibrium (LD) pruning.
  - Executes Principal Component Analysis (PCA) to correct for population stratification.
  - Conducts variant-phenotype association testing to isolate statistically significant single nucleotide polymorphisms (SNPs).

- **Module 2: Variant Annotation (Ensembl REST API)**
  - Maps candidate SNPs to genomic coordinates.
  - Filters variants for protein-altering consequences (missense/nonsense).
  - Retrieves both wild-type and mutated amino acid sequences.

- **Module 3: Structural Bioinformatics (AlphaFold / ESMFold)**
  - Translates candidate mutated sequences into 3D structural manifolds.
  - Generates high-fidelity structural data essential for downstream binding analysis.

- **Module 4: Deep Learning DTI Screening (DeepPurpose / PyTorch)**
  - Computes Morgan fingerprints for the provided chemical library.
  - Utilizes convolutional neural networks (CNNs) trained on the BindingDB dataset to predict binding affinities (pIC50/Kd) between the mutated protein targets and candidate ligands.

## 3. Tech Stack & Concepts

This pipeline was designed to adhere to modern bioinformatics and MLOps best practices:

- **Workflow Orchestration:** Advanced Nextflow DSL2 scripting for parallel execution, process caching (`-resume`), and channel-based dataflow management.
- **Environment Virtualization:** Strict Docker containerization for every module, eliminating "works on my machine" dependencies and ensuring cross-platform reproducibility.
- **MLOps & Dependency Management:** Custom Dockerfile engineering to resolve complex C++ compiler conflicts between PyTorch Intel MKL bindings and cheminformatics libraries (`descriptastorus`, RDKit).
- **Multi-Omics Integration:** Seamless programmatic handoffs between statistical genomics (VCF/PLINK), web APIs (Ensembl), and deep learning tensors (PyTorch).

## 4. Inputs and Outputs

### Required Inputs

The pipeline accepts standardized biological file formats:

- `--vcf`: A multi-sample Variant Call Format (VCFv4.2) file containing genotype data.
- `--phenotype`: A tab-separated values (TSV) file mapping samples to continuous or binary traits (format: `sample_id\tphenotype`).
- `--coconut`: A SMILES (`.smi`) file containing the chemical library for screening (e.g., COCONUT database extracts).
- `--annotation_tsv` *(optional/testing)*: A pre-curated TSV fixture for offline variant-to-protein mapping, used to bypass live API rate limits.

### Generated Outputs

Results are systematically organized into a structured `--outdir`:

- `01_genomics/` — Filtered candidate loci (TSV), PCA scatter plots, Manhattan plots, and QQ plots (PNG).
- `02_annotation/` — Consequence mappings (TSV) and extracted protein sequences (FASTA).
- `03_structure/` — AlphaFold PDB predictions and structural linkage manifests.
- `04_cheminformatics/` — Computed ligand descriptor matrices.
- `05_dti/` — Raw neural network affinity prediction matrices.
- `06_results/` — The final `top_predicted_candidates.tsv` detailing ranked ligand-target pairs, accompanied by summary visualizations.

## 5. Dependencies

PhenoScreen abstracts underlying software dependencies away from the host system. The only host requirements are:

- **Nextflow** (≥ 22.04.0)
- **Docker** Engine

Internal containerized environments include:

- `phenoscreen/genomics` — PLINK 1.9, R (ggplot2, data.table).
- `phenoscreen/structure` — Python 3.10, BioPython, Requests.
- `phenoscreen/cheminformatics` — PyTorch (CPU-optimized), DeepPurpose, RDKit, pandas, descriptastorus.

## 6. Testing & Usage Instructions

Ensure Nextflow and Docker are installed and running on your host machine before executing the commands below.

### Option A: Run the Synthetic Smoke Test (Fast CI/CD Verification)

This mode uses a lightweight synthetic dataset to verify the integrity of the Nextflow DAG and container handoffs without requiring external API lookups.

```bash
nextflow run main.nf \
  -profile docker \
  --vcf test_data/sample.vcf \
  --phenotype test_data/sample_phenotype.tsv \
  --coconut test_data/sample_coconut.smi \
  --struct_backend mock \
  --dti_mode mock \
  --outdir results_test
```

### Option B: Run the NCI-60 "Gold-Standard" Production Test

This mode executes the pipeline using a curated integration fixture containing real oncogenic variants (BRAF V600E, TP53 R273H, KRAS G12S) and actual Paclitaxel drug-sensitivity metrics, testing live structural modeling and PyTorch DTI inference.

```bash
nextflow run main.nf \
  -profile docker \
  --vcf real_data/nci60/nci60_mutations.vcf \
  --phenotype real_data/nci60/phenotype.tsv \
  --coconut test_data/sample_coconut.smi \
  --annotation_tsv real_data/nci60/nci60_annotation.tsv \
  --struct_backend alphafold \
  --dti_mode deepurpose \
  --candidate_top_n 30 \
  --max_annotation_variants 30 \
  --max_compounds 8 \
  --top_n 10 \
  --outdir results_nci60_production
```

## 7. Current Limitations

Maintaining scientific rigor requires distinguishing between computational validation and biological discovery:

- **Statistical Power:** The NCI-60 validation fixture (N=10 samples) is engineered strictly as a computational integration test. It lacks the statistical power required for a true Genome-Wide Association Study (GWAS).
- **Somatic vs. Germline Assumptions:** The current genomics module uses PLINK, which is optimized for germline inherited variance. Adapting the mathematical models specifically for somatic tumor mutation burdens is required for true oncology cohort deployments.
- **Compute Constraints:** The cheminformatics Docker image is currently optimized for CPU execution to allow continuous integration and local Codespace testing. Large-scale structural predictions and DTI screening require GPU-accelerated (CUDA) PyTorch images.

## 8. Future Applications & Scalability

PhenoScreen is built with the explicit intention of scaling from local development environments to institutional research clusters.

- **HPC Deployment:** By swapping the Nextflow `-profile docker` flag for a Singularity/SLURM profile, the pipeline can be immediately deployed on university supercomputers.
- **Clinical Cohort Analysis:** The architecture is designed to ingest large-scale, real-world pharmacogenomic datasets (e.g., DepMap, CCLE, or TCGA) to systematically identify novel drug repurposing candidates for treatment-resistant cancer sub-clones.

## Disclaimer

This pipeline produces **computational predictions only**, statistical associations, predicted protein structures, and predicted binding affinities. Results require experimental validation before any biological or pharmacological conclusions are drawn.
