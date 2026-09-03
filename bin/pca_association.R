#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(SNPRelate)
  library(data.table)
  library(ggplot2)
})

# ============================================================
# COMMAND-LINE OPTIONS
# ============================================================

opt <- parse_args(
  OptionParser(
    option_list = list(
      make_option("--bfile_prefix", type = "character"),
      make_option("--phenotype", type = "character"),
      make_option("--phenotype_type", type = "character", default = "continuous"),
      make_option("--n_pcs", type = "integer", default = 10),
      make_option("--top_n", type = "integer", default = 20),
      make_option("--fdr", type = "double", default = 0.05)
    )
  )
)

# ============================================================
# CREATE GDS FILE FROM PLINK BED/BIM/FAM
# ============================================================

dir.create("pca_tmp", showWarnings = FALSE)

gds_file <- "pca_tmp/genotypes.gds"

snpgdsBED2GDS(
  paste0(opt$bfile_prefix, ".bed"),
  paste0(opt$bfile_prefix, ".fam"),
  paste0(opt$bfile_prefix, ".bim"),
  gds_file,
  verbose = FALSE
)

g <- snpgdsOpen(gds_file)

on.exit(
  snpgdsClose(g),
  add = TRUE
)

# ============================================================
# READ GENOTYPE INFORMATION
# ============================================================

sample_ids <- read.gdsn(
  index.gdsn(g, "sample.id")
)

snp_id <- read.gdsn(
  index.gdsn(g, "snp.id")
)

# Try to obtain rsIDs
rsid_node <- index.gdsn(
  g,
  "snp.rs.id",
  silent = TRUE
)

if (!is.null(rsid_node)) {
  snp_rsid <- read.gdsn(rsid_node)
} else {
  snp_rsid <- as.character(snp_id)
}

snp_chr <- read.gdsn(
  index.gdsn(g, "snp.chromosome")
)

snp_pos <- read.gdsn(
  index.gdsn(g, "snp.position")
)

cat(sprintf(
  "Loaded %d samples and %d SNPs.\n",
  length(sample_ids),
  length(snp_id)
))

# ============================================================
# PCA
# ============================================================

n_pcs <- min(
  opt$n_pcs,
  max(1, length(sample_ids) - 1),
  max(1, length(snp_id) - 1)
)

pca <- snpgdsPCA(
  g,
  num.thread = 2,
  eigen.cnt = n_pcs,
  verbose = FALSE
)

pc <- as.data.table(
  pca$eigenvect[, seq_len(n_pcs), drop = FALSE]
)

colnames(pc) <- paste0(
  "PC",
  seq_len(n_pcs)
)

pc$sample_id <- as.character(
  pca$sample.id
)

# Save PCA eigenvectors

fwrite(
  pc[, c(
    "sample_id",
    paste0("PC", seq_len(n_pcs))
  ), with = FALSE],
  "pca_eigenvectors.tsv",
  sep = "\t"
)

# ============================================================
# PCA SCREE PLOT
# ============================================================

vp <- data.frame(
  PC = seq_len(n_pcs),
  VarianceExplained =
    100 * pca$varprop[seq_len(n_pcs)]
)

ggsave(
  "pca_scree_plot.png",
  ggplot(
    vp,
    aes(PC, VarianceExplained)
  ) +
    geom_col() +
    labs(
      y = "% variance explained",
      title = "PCA scree plot"
    ),
  width = 6,
  height = 4,
  dpi = 150
)

# ============================================================
# READ PHENOTYPE
# ============================================================

ph <- fread(
  opt$phenotype
)

if (!all(
  c("sample_id", "phenotype") %in% names(ph)
)) {
  stop(
    "Phenotype TSV needs columns: sample_id and phenotype."
  )
}

ph[, sample_id := as.character(sample_id)]

if (opt$phenotype_type == "binary") {

  ph$phenotype <- as.factor(
    ph$phenotype
  )

} else {

  ph$phenotype <- as.numeric(
    ph$phenotype
  )
}

# Check binary phenotype

if (
  opt$phenotype_type == "binary" &&
  length(unique(na.omit(ph$phenotype))) != 2
) {
  stop(
    "Binary phenotype must contain exactly two observed classes."
  )
}

# ============================================================
# MATCH PHENOTYPE WITH GENOTYPE SAMPLES
# ============================================================

m <- merge(
  ph,
  pc,
  by = "sample_id"
)

cat(sprintf(
  "Matched %d phenotype samples to %d genotype samples.\n",
  nrow(m),
  nrow(pc)
))

if (nrow(m) < 10) {
  stop(
    "Too few phenotype samples matched genotype sample IDs. ",
    "Check sample_id naming between genotype and phenotype files."
  )
}

# ============================================================
# SELECT PCA COVARIATES
# ============================================================

# Use up to 3 PCs for small datasets such as NCI-60

n_cov_pcs <- min(
  3,
  n_pcs
)

pcs <- paste0(
  "PC",
  seq_len(n_cov_pcs)
)

terms <- paste(
  pcs,
  collapse = " + "
)

# ============================================================
# LOAD GENOTYPES
# ============================================================

geno <- snpgdsGetGeno(
  g,
  verbose = FALSE
)

cat("Starting association analysis...\n")

res <- vector(
  "list",
  length(snp_id)
)

# ============================================================
# SNP ASSOCIATION ANALYSIS
# ============================================================

for (i in seq_along(snp_id)) {

  # Match genotype values to phenotype samples

  m$genotype <- geno[, i][
    match(
      m$sample_id,
      pc$sample_id
    )
  ]

  # Remove missing genotype/phenotype values

  analysis_data <- m[
    is.finite(genotype) &
      !is.na(phenotype)
  ]

  # Require enough samples

  if (nrow(analysis_data) < 10) {
    next
  }

  # Skip monomorphic SNPs

  if (
    length(
      unique(analysis_data$genotype)
    ) < 2
  ) {
    next
  }

  # Build association model

  f <- as.formula(
    paste(
      "phenotype ~ genotype +",
      terms
    )
  )

  # Run regression

  fit <- tryCatch(

    if (opt$phenotype_type == "binary") {

      glm(
        f,
        data = analysis_data,
        family = binomial()
      )

    } else {

      lm(
        f,
        data = analysis_data
      )
    },

    error = function(e) NULL
  )

  # Skip failed models

  if (is.null(fit)) {
    next
  }

  # Extract coefficients

  co <- summary(fit)$coefficients

  # Check genotype term exists

  if (
    !"genotype" %in% rownames(co)
  ) {
    next
  }

  # Save result

  res[[i]] <- data.table(

    snp_id = i,

    rsid = snp_rsid[i],

    chromosome = snp_chr[i],

    position = snp_pos[i],

    beta = co["genotype", "Estimate"],

    se = co["genotype", "Std. Error"],

    p_value = co[
      "genotype",
      ncol(co)
    ]
  )
}

# ============================================================
# COMBINE ASSOCIATION RESULTS
# ============================================================

valid_results <- Filter(
  Negate(is.null),
  res
)

# ============================================================
# HANDLE NO VALID ASSOCIATIONS
# ============================================================

if (length(valid_results) == 0) {

  cat(
    "WARNING: No valid SNP association models were produced.\n"
  )

  empty_results <- data.table(

    snp_id = integer(),

    rsid = character(),

    chromosome = character(),

    position = numeric(),

    beta = numeric(),

    se = numeric(),

    p_value = numeric(),

    fdr = numeric(),

    effect_size = numeric(),

    selection_status = character()
  )

  fwrite(
    empty_results,
    "association_results.tsv",
    sep = "\t"
  )

  fwrite(
    empty_results,
    "candidate_loci.tsv",
    sep = "\t"
  )

  a <- empty_results

  sig <- empty_results

  cand <- empty_results

} else {

  # Combine results

  a <- rbindlist(
    valid_results,
    fill = TRUE
  )

  # Remove invalid p-values

  a <- a[
    is.finite(p_value) &
      !is.na(p_value)
  ]

  # ==========================================================
  # HANDLE RESULTS AFTER FILTERING
  # ==========================================================

  if (nrow(a) == 0) {

    cat(
      "WARNING: No valid p-values after filtering.\n"
    )

    empty_results <- data.table(

      snp_id = integer(),

      rsid = character(),

      chromosome = character(),

      position = numeric(),

      beta = numeric(),

      se = numeric(),

      p_value = numeric(),

      fdr = numeric(),

      effect_size = numeric(),

      selection_status = character()
    )

    fwrite(
      empty_results,
      "association_results.tsv",
      sep = "\t"
    )

    fwrite(
      empty_results,
      "candidate_loci.tsv",
      sep = "\t"
    )

    a <- empty_results

    sig <- empty_results

    cand <- empty_results

  } else {

    # ========================================================
    # FDR CORRECTION
    # ========================================================

    a[, fdr := p.adjust(
      p_value,
      method = "BH"
    )]

    a[, effect_size := beta]

    setorder(
      a,
      p_value
    )

    a[, selection_status := "not_selected"]

    a[
      fdr < opt$fdr,
      selection_status := "significant_fdr"
    ]

    # Save all association results

    fwrite(
      a,
      "association_results.tsv",
      sep = "\t"
    )

    # ========================================================
    # SELECT CANDIDATES
    # ========================================================

    sig <- a[
      selection_status == "significant_fdr"
    ]

    expl <- a[
      selection_status == "not_selected"
    ][
      seq_len(
        min(
          opt$top_n,
          .N
        )
      )
    ]

    # Use significant SNPs if available.
    # Otherwise use exploratory top SNPs.

    if (nrow(sig) > 0) {

      cand <- sig

    } else {

      cand <- copy(expl)

      cand[
        ,
        selection_status := "exploratory_top_n"
      ]
    }

    fwrite(
      cand,
      "candidate_loci.tsv",
      sep = "\t"
    )
  }
}

# ============================================================
# PCA SCATTER PLOT
# ============================================================

if (
  nrow(pc) >= 2 &&
  "PC1" %in% names(pc) &&
  "PC2" %in% names(pc)
) {

  p <- ggplot(
    pc,
    aes(
      PC1,
      PC2
    )
  ) +
    geom_point() +
    labs(
      title = "Population-structure PCA"
    )

  ggsave(
    "pca_scatter_plot.png",
    p,
    width = 6,
    height = 5,
    dpi = 150
  )
}

# ============================================================
# ASSOCIATION PLOTS
# ============================================================

if (nrow(a) > 0) {

  a[
    ,
    neglogp :=
      -log10(
        pmax(
          p_value,
          .Machine$double.xmin
        )
      )
  ]

  # Manhattan plot

  man <- ggplot(
    a,
    aes(
      x = position,
      y = neglogp,
      color = factor(chromosome)
    )
  ) +
    geom_point(size = 0.8) +
    labs(
      x = "Genomic position",
      color = "Chromosome",
      y = "-log10(p-value)",
      title = "Association scan"
    )

  ggsave(
    "association_manhattan.png",
    man,
    width = 8,
    height = 4,
    dpi = 150
  )

  # QQ plot

  expected <- -log10(
    ppoints(
      nrow(a)
    )
  )

  observed <- -log10(
    sort(a$p_value)
  )

  qq <- ggplot(
    data.frame(
      expected = expected,
      observed = observed
    ),
    aes(
      expected,
      observed
    )
  ) +
    geom_point(size = 0.8) +
    geom_abline(
      slope = 1,
      intercept = 0
    ) +
    labs(
      x = "Expected -log10(p)",
      y = "Observed -log10(p)",
      title = "Association QQ plot"
    )

  ggsave(
    "association_qq.png",
    qq,
    width = 5,
    height = 5,
    dpi = 150
  )

} else {

  file.create(
    "association_manhattan.png"
  )

  file.create(
    "association_qq.png"
  )
}

# ============================================================
# ENSURE PCA PLOT EXISTS
# ============================================================

if (!file.exists("pca_scatter_plot.png")) {

  file.create(
    "pca_scatter_plot.png"
  )
}

# ============================================================
# FINAL SUMMARY
# ============================================================

cat(
  sprintf(
    paste0(
      "Association complete: ",
      "%d tested; ",
      "%d FDR-significant; ",
      "%d downstream candidates.\n"
    ),

    nrow(a),

    nrow(sig),

    nrow(cand)
  )
)