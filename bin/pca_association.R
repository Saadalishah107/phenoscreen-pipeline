#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(SNPRelate); library(data.table); library(ggplot2)
})
opt <- parse_args(OptionParser(option_list=list(
  make_option("--bfile_prefix", type="character"),
  make_option("--phenotype", type="character"),
  make_option("--phenotype_type", type="character", default="continuous"),
  make_option("--n_pcs", type="integer", default=10),
  make_option("--top_n", type="integer", default=20),
  make_option("--fdr", type="double", default=0.05)
)))

dir.create("pca_tmp", showWarnings=FALSE)
gds_file <- "pca_tmp/genotypes.gds"
snpgdsBED2GDS(paste0(opt$bfile_prefix,".bed"),
              paste0(opt$bfile_prefix,".fam"),
              paste0(opt$bfile_prefix,".bim"),
              gds_file, verbose=FALSE)
g <- snpgdsOpen(gds_file)
on.exit(snpgdsClose(g), add=TRUE)

sample_ids <- read.gdsn(index.gdsn(g,"sample.id"))
snp_id <- read.gdsn(index.gdsn(g,"snp.id"))

# Safely extract snp.rs.id with a fallback if SNPRelate drops the node
rsid_node <- index.gdsn(g, "snp.rs.id", silent = TRUE)
if (!is.null(rsid_node)) {
  snp_rsid <- read.gdsn(rsid_node)
} else {
  snp_rsid <- paste0("variant_", snp_id)
}

snp_chr <- read.gdsn(index.gdsn(g,"snp.chromosome"))
snp_pos <- read.gdsn(index.gdsn(g,"snp.position"))

n_pcs <- min(opt$n_pcs, max(1, length(sample_ids)-1), max(1, length(snp_id)-1))
pca <- snpgdsPCA(g, num.thread=2, eigen.cnt=n_pcs, verbose=FALSE)

# Cast as data.table to support data.table subsetting syntax below
pc <- as.data.table(pca$eigenvect[,seq_len(n_pcs),drop=FALSE])
colnames(pc) <- paste0("PC",seq_len(n_pcs))
pc$sample_id <- as.character(pca$sample.id)
fwrite(pc[,c("sample_id",paste0("PC",seq_len(n_pcs))),with=FALSE],
       "pca_eigenvectors.tsv", sep="\t")

vp <- data.frame(PC=seq_len(n_pcs),VarianceExplained=100*pca$varprop[seq_len(n_pcs)])
ggsave("pca_scree_plot.png", ggplot(vp,aes(PC,VarianceExplained))+
       geom_col()+labs(y="% variance explained",title="PCA scree plot"),
       width=6,height=4,dpi=150)

ph <- fread(opt$phenotype)
if (!all(c("sample_id","phenotype") %in% names(ph))) stop("Phenotype TSV needs sample_id and phenotype.")
ph[,sample_id:=as.character(sample_id)]
ph$phenotype <- if(opt$phenotype_type=="binary") as.factor(ph$phenotype) else as.numeric(ph$phenotype)
if(opt$phenotype_type=="binary" && length(unique(na.omit(ph$phenotype))) != 2) stop("Binary phenotype must contain exactly two observed classes.")
m <- merge(ph, pc, by="sample_id")
pcs <- paste0("PC",seq_len(min(5,n_pcs)))
terms <- paste(pcs,collapse=" + ")

geno <- snpgdsGetGeno(g, verbose=FALSE)
res <- vector("list",length(snp_id))
for(i in seq_along(snp_id)){
  m$genotype <- geno[,i][match(m$sample_id, pc$sample_id)]
  f <- if(opt$phenotype_type=="binary")
        as.formula(paste("phenotype ~ genotype +",terms))
        else as.formula(paste("phenotype ~ genotype +",terms))
  fit <- tryCatch(if(opt$phenotype_type=="binary") glm(f,data=m,family=binomial())
                  else lm(f,data=m), error=function(e) NULL)
  if(is.null(fit)) next
  co <- summary(fit)$coefficients
  if(!"genotype" %in% rownames(co)) next
  res[[i]] <- data.table(
    snp_id=snp_id[i], rsid=snp_rsid[i], chromosome=snp_chr[i], position=snp_pos[i],
    beta=co["genotype","Estimate"], se=co["genotype","Std. Error"],
    p_value=co["genotype",ncol(co)]
  )
}
a <- rbindlist(res,fill=TRUE)
a <- a[is.finite(p_value)]
a[,fdr:=p.adjust(p_value,method="BH")]
a[,effect_size:=beta]
setorder(a,p_value)
a[,selection_status:="not_selected"]
a[fdr < opt$fdr,selection_status:="significant_fdr"]
fwrite(a,"association_results.tsv",sep="\t")

sig <- a[selection_status=="significant_fdr"]
expl <- a[selection_status=="not_selected"][seq_len(min(opt$top_n,.N))]
if(nrow(sig)>0) cand <- sig else { cand <- expl; cand[,selection_status:="exploratory_top_n"] }
fwrite(cand,"candidate_loci.tsv",sep="\t")

if(nrow(pc)>=2 && "PC1" %in% names(pc) && "PC2" %in% names(pc)){
  p <- ggplot(pc,aes(PC1,PC2))+geom_point()+labs(title="Population-structure PCA")
  ggsave("pca_scatter_plot.png",p,width=6,height=5,dpi=150)
}
if(nrow(a)>0){
  a[,neglogp:=-log10(pmax(p_value,.Machine$double.xmin))]
  man <- ggplot(a,aes(x=position,y=neglogp,color=factor(chromosome)))+geom_point(size=.8)+
         labs(x="Genomic position",color="Chromosome",y="-log10(p-value)",title="Association scan")
  ggsave("association_manhattan.png",man,width=8,height=4,dpi=150)
  expected <- -log10(ppoints(nrow(a)))
  observed <- -log10(sort(a$p_value))
  qq <- ggplot(data.frame(expected=expected,observed=observed),aes(expected,observed))+
        geom_point(size=.8)+geom_abline(slope=1,intercept=0)+
        labs(x="Expected -log10(p)",y="Observed -log10(p)",title="Association QQ plot")
  ggsave("association_qq.png",qq,width=5,height=5,dpi=150)
} else {
  file.create("association_manhattan.png")
  file.create("association_qq.png")
}
if (!file.exists("pca_scatter_plot.png")) file.create("pca_scatter_plot.png")
cat(sprintf("Association complete: %d tested; %d FDR-significant; %d downstream candidates.\n",
            nrow(a),nrow(sig),nrow(cand)))