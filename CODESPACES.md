# GitHub Codespaces test guide

## 1. Open the repository

Push this repository to GitHub, open it, then choose **Code → Codespaces → Create codespace on main**.

The repository contains `.devcontainer/devcontainer.json`, so the Codespace installs Docker and Java 17 and installs Nextflow automatically.

Check:

```bash
nextflow -version
docker version
git status
```

Nextflow 26.04 uses the current v2 language parser by default; the pipeline is written for current DSL2 syntax rather than the older per-process publishing style.

## 2. Build only the lightweight test containers

Do this first. The smoke test intentionally uses:

- a deterministic annotation fixture instead of live Ensembl;
- `mock` structure generation;
- a deterministic DTI surrogate instead of DeepPurpose.

This proves the **Nextflow wiring, containers, channels, files and reports** without downloading a large ML environment.

```bash
./docker/build_images.sh
```

Check:

```bash
docker images | grep phenoscreen
```

You should see:

- `phenoscreen/genomics:2.0`
- `phenoscreen/structure:2.0`
- `phenoscreen/cheminformatics-test:2.0`

## 3. Run the complete smoke test

```bash
nextflow run . -profile test,docker -resume
```

A successful run should finish every process with `✔`.

Then inspect:

```bash
find results_test -maxdepth 3 -type f | sort
```

Important files:

```text
results_test/
├── 01_genomics/
│   ├── candidates/candidate_loci.tsv
│   ├── association/association_results.tsv
│   └── pca/
├── 02_annotation/
│   ├── variant_annotation.tsv
│   └── protein_sequences.fasta
├── 03_structure/
│   └── structure_manifest.tsv
├── 04_cheminformatics/
│   └── molecular_descriptors.tsv
├── 05_dti/
│   └── affinity_matrix.tsv
├── 06_results/
│   └── top_predicted_candidates.tsv
└── pipeline_info/
```

## 4. Validate the outputs

Run:

```bash
test -s results_test/01_genomics/candidates/candidate_loci.tsv
test -s results_test/01_genomics/association/association_results.tsv
test -s results_test/02_annotation/variant_annotation.tsv
test -s results_test/02_annotation/protein_sequences.fasta
test -s results_test/03_structure/structure_manifest.tsv
test -s results_test/04_cheminformatics/molecular_descriptors.tsv
test -s results_test/05_dti/affinity_matrix.tsv
test -s results_test/06_results/top_predicted_candidates.tsv
```

Inspect the biological hand-off:

```bash
column -t -s $'\t' results_test/02_annotation/variant_annotation.tsv | head -10
```

Inspect predicted candidates:

```bash
column -t -s $'\t' results_test/06_results/top_predicted_candidates.tsv | head -15
```

## 5. Inspect the DAG and provenance

Open the generated HTML files under:

```text
results_test/pipeline_info/
```

Also inspect:

```bash
nextflow log
```

For a failed run, use:

```bash
nextflow log
```

Then rerun with:

```bash
nextflow run . -profile test,docker -resume
```

`-resume` is important: Nextflow can reuse successful cached tasks and rerun only the failed portion.

## 6. Run the live annotation path

The production annotation module uses Ensembl VEP REST. Do not use the offline fixture.

For example:

```bash
nextflow run .   --vcf /absolute/path/cohort.vcf.gz   --phenotype /absolute/path/phenotype.tsv   --coconut /absolute/path/coconut.sdf   --species human   --offline_annotation false   --struct_backend alphafold   --dti_mode mock   --outdir results_live
```

For a real cohort, use a real annotation reference and expect some variants to have no protein-altering consequence. Those variants are correctly excluded from the structure/DTI branch.

## 7. Run real DeepPurpose

The real model environment is deliberately separate because it is much heavier.

Build it:

```bash
./docker/build_images.sh --full
```

Then run:

```bash
nextflow run .   -profile docker   --vcf /absolute/path/cohort.vcf.gz   --phenotype /absolute/path/phenotype.tsv   --coconut /absolute/path/coconut.sdf   --struct_backend alphafold   --dti_mode deepurpose   --dti_model Morgan_CNN_BindingDB   --outdir results_real
```

DeepPurpose is a sequence/SMILES DTI model here. The PDB files are produced in the parallel structure branch for structural characterization and downstream structural work; the DTI model does **not** consume PDB coordinates.

## 8. Test only one module when debugging

Use Nextflow's process selectors. Examples:

```bash
nextflow run . -profile test,docker -stub-run
```

and, after a full run, inspect the task directory from the Nextflow log.

The safest debugging loop is:

```bash
nextflow run . -profile test,docker -resume
nextflow log
```

Then inspect the failing task directory and its `.command.sh`, `.command.out`, and `.command.err`.

## 9. Clean and repeat from scratch

```bash
make clean
make test
```

Do not delete `work/` while you are investigating a failure; it contains the exact staged inputs and command scripts needed for debugging.

## 10. What counts as a successful test?

For the test profile, success means:

1. PLINK creates a valid pruned dataset.
2. R/SNPRelate consumes all three PLINK files (`.bed`, `.bim`, `.fam`).
3. PCA and association outputs are generated.
4. Candidate variants are explicitly labelled significant or exploratory.
5. The candidate-to-protein transition is exercised using a local fixture.
6. Structure artifacts are generated by the mock backend.
7. RDKit standardizes the small test library and calculates descriptors.
8. DTI and ranking run without using structural coordinates.
9. The final table is called **predicted candidates**, not validated hits.
10. Nextflow produces timeline/report/trace/DAG metadata.

The mock DTI score is only an integration-test surrogate and must never be interpreted as a biological prediction.
