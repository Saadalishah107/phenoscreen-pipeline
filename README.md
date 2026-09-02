# PhenoScreen v2 — Genotype-Guided Natural-Product Virtual Screening

PhenoScreen v2 is a modular Nextflow DSL2 workflow for a research-software portfolio:

```text
WGS/VCF
  │
  ├── PLINK conversion + LD pruning
  │
  └── R/SNPRelate PCA + phenotype association
             │
             ▼
       candidate variants
             │
             ▼
       variant annotation
       gene/transcript/consequence
             │
             ▼
       protein sequence
          ┌──┴───────────────┐
          │                  │
          ▼                  ▼
   AlphaFold DB/ESMFold   DeepPurpose
   structural artifact    sequence + SMILES DTI
          │                  │
          └──────┬───────────┘
                 │
          ranked predicted
             candidates
```

The workflow deliberately separates **structural prediction** from **sequence-based DTI**. A PDB is not silently treated as an input to DeepPurpose.

## Scientific status

This repository is a computational workflow and virtual-screening prototype. Its outputs are predictions, not experimental validation, binding measurements, or clinical recommendations.

The included synthetic data are for integration testing only.

## Major v2 corrections

- Correct PLINK `.bed/.bim/.fam` tuple wiring into R.
- Real Ensembl VEP REST annotation path.
- Explicit protein-altering consequence filtering.
- Protein sequence resolution before structure/DTI.
- AlphaFold DB lookup with ESMFold fallback.
- Structural prediction is a parallel branch, not a false DTI input.
- Continuous and binary phenotype association modes.
- FDR correction and explicit `significant_fdr` versus `exploratory_top_n` labels.
- Expanded RDKit descriptors including QED, fraction Csp3, heavy atoms, rings and charge.
- DeepPurpose model names use the full pretrained model identifier.
- Batch-aware, sequence/SMILES DTI implementation.
- Final output is named `top_predicted_candidates.tsv`.
- Lightweight deterministic test profile.
- Centralized workflow-level output publishing using current Nextflow output syntax.
- GitHub Actions smoke test that actually executes the pipeline.

## Quick start

Install Nextflow and Docker, then:

```bash
./docker/build_images.sh
nextflow run . -profile test,docker -resume
```

For GitHub Codespaces, see [CODESPACES.md](CODESPACES.md).

## Production example

```bash
nextflow run .   --vcf cohort.vcf.gz   --phenotype phenotype.tsv   --coconut coconut.sdf   --phenotype_type continuous   --species human   --struct_backend alphafold   --dti_mode deepurpose   --dti_model Morgan_CNN_BindingDB   --outdir results
```

For very large COCONUT libraries, use a filtered subset or add an explicit screening budget before launching a full cross-product.

## Input phenotype

Continuous:

```text
sample_id    phenotype
SAMPLE000    39.16
SAMPLE001    50.77
```

Binary:

```text
sample_id    phenotype
SAMPLE000    0
SAMPLE001    1
```

## Annotation behavior

Production annotation calls Ensembl VEP REST for candidate variants and retrieves the Ensembl protein sequence for a selected protein-altering transcript consequence.

The test profile uses `test_data/sample_annotation.tsv` so the smoke test is deterministic and does not depend on external annotation availability.

## Structure behavior

`alphafold` first attempts AlphaFold DB when a UniProt accession is available and falls back to ESMFold for unresolved targets.

`esmfold` directly uses the ESMFold API.

`mock` is only for CI/integration testing.

## DTI behavior

`deepurpose` uses a pretrained DeepPurpose model with protein sequence + compound SMILES. The default is `Morgan_CNN_BindingDB`.

`mock` is a deterministic non-biological surrogate used only to test Nextflow integration.

DeepPurpose's own documentation lists pretrained model names such as `Morgan_CNN_BindingDB` and explains that its DTI inputs are drug SMILES and protein amino-acid sequences. The model should therefore not be described as consuming PDB coordinates directly.

## Repository layout

```text
.
├── main.nf
├── workflows/
│   └── phenoscreen.nf
├── modules/
│   ├── genomics/
│   ├── annotation/
│   ├── structure/
│   └── cheminformatics/
├── bin/
├── conf/
├── docker/
├── test_data/
├── .devcontainer/
├── .github/workflows/ci.yml
├── CODESPACES.md
└── nextflow.config
```

## Reproducibility

Nextflow manages task isolation and caching. Use `-resume` to restart a failed run without repeating successful tasks.

The pipeline emits timeline, report, trace and DAG metadata under `pipeline_info/`.

## License

MIT.
