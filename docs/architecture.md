# Architecture

## Dataflow

```text
VCF + phenotype
      |
      v
LD_PRUNING
      |
      | tuple(bed,bim,fam)
      v
PCA_ASSOCIATION
      |
      v
ANNOTATE_VARIANTS
      |
      +----------------------+
      |                      |
      v                      v
PREDICT_STRUCTURE      PROCESS_CHEM_LIBRARY
                             |
                             v
                       DTI_SCREENING
                             |
                             v
                         RANK_HITS
```

## Important scientific boundaries

### Variant annotation

A genomic identifier is not a protein identifier. Candidate variants are annotated to transcript consequences, then protein sequences are retrieved. Only protein-altering consequences are passed downstream.

### Structure

The structure branch is independent of the DTI model. Its purpose is structural characterization and a reusable artifact for later docking/structural analysis.

### DTI

DeepPurpose receives amino-acid sequences and compound SMILES. The model output is a predicted score from the selected pretrained model. It is not an experimental binding constant.

### Candidate selection

Association results contain beta, standard error, p-value and FDR. A variant is marked `significant_fdr` only when it passes the configured FDR threshold. If none pass, the pipeline can select the configured top-N variants as `exploratory_top_n`; those are not significant findings.

## Nextflow design

The entry workflow uses the current DSL2 workflow `publish:`/`output {}` model. Process modules do not contain `publishDir`. This keeps dataflow modules reusable and centralizes output organization.

A named subworkflow (`PHENOSCREEN`) exposes its outputs through `emit:` so the entry workflow controls publication.

## Test design

The test profile is deliberately deterministic:

- synthetic VCF;
- synthetic phenotype;
- local annotation fixture;
- mock structure generator;
- mock DTI score;
- small SMILES library.

This tests workflow correctness without claiming biological evidence.
