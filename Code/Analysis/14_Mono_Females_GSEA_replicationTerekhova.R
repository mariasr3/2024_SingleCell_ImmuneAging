library(AUCell);library(dplyr);library(GSEABase);library(stringr);library(fitdistrplus);library(Seurat);library(SingleCellExperiment); library(escape)

Csparse_validate = "CsparseMatrix_validate"

# Perform ennrichment of interferon and inflammatory genes of classical monocytes from Terekhova et al. 2023

#paths 
if (getwd()== "/Users/mariasopenar"){
  basepath <- "/Users/mariasopenar/cluster/"
  data_path <- "/Users/mariasopenar/cluster/Projects/scRNAseq/"
  
}else if(getwd()== "/home/mariasr"){
  basepath <- "/home/mariasr/cluster/"
  data_path <- "/home/mariasr/cluster/Projects/scRNAseq/"
}else{
  basepath <- "/gpfs/projects/bsc83/"
  data_path <- "/gpfs/projects/bsc83/Projects/scRNAseq/"
}
plots_path <- paste0(data_path, "msopena/02_OneK1K_Age/plots/")
robjects <- paste0(data_path, "msopena/02_OneK1K_Age/robjects/")
#dir.create(path_enrichScore, recursive = TRUE)

#functions
source(paste0(data_path, "/msopena/02_OneK1K_Age/scripts/functions.R"))

shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--sex"), action="store", default=NA, type='character',
              help="M, F"), 
  make_option(c("--gene_set"), action="store", default=NA, type='character',
              help="hallmark"))

opt = parse_args(OptionParser(option_list=option_list))

sex <- opt$sex
gene_set_name <- opt$gene_set
celltype <- "CD14_Mono"
method <- "AUCell"


# 0. Subset expression data for the classical monoyctes Terekhova----
# mtx <- readRDS(paste0(basepath, "Data/scRNAseq/Terekhova2023/GEX_HTO_processed/myeloid_cells/myeloid_cells_rna.rds"))
# md <- read.csv(paste0(basepath, "Data/scRNAseq/Terekhova2023/GEX_HTO_processed/myeloid_cells/myeloid_cells_metadata.csv"))
# 
# # get cell type we are interested in 
# md_class_mono <- md[md$Cluster_names == "Classical monocytes", ]
# dim(md_class_mono)
# rownames(md_class_mono) <- md_class_mono$X
# 
# #subset also matrix 
# mtx_class_mono <- mtx[, md_class_mono$X]
# dim(mtx_class_mono)
# #create Single Cell experiment object
# sce_class_mono <- SingleCellExperiment(list(counts=mtx_class_mono), colData=DataFrame(md_class_mono))
# 
# #save the subset 
# saveRDS(sce_class_mono, paste0(basepath, "Data/scRNAseq/Terekhova2023/GEX_HTO_processed/myeloid_cells/classical_monocytes_sce.rds"))


# 0. Subset expression data for CD14 Monocytes Oelen ----

so <- readRDS(paste0(data_path, "/aripol1/wijst-2020-hg19/v1/aging/00.so_split_by_celltype/v2/so.rds"))
mdata <- readRDS(paste0(data_path, "/aripol1/wijst-2020-hg19/v1/aging/00.so_split_by_celltype/v2/md.rds"))
identical(colnames(so), mdata$bare_barcode_lane)
sce_class_mono <- SingleCellExperiment(list(counts=so@assays$RNA@counts), colData=DataFrame(so@meta.data))


print(paste0("##### Computing Enrich GSEA using ",  method, " for Classical Monocytes #####"))

# 1. Get expression data 
print(paste0("1. Loading data and subseting cells for ", sex, "----------"))
#sce_class_mono <- readRDS(paste0(basepath, "Data/scRNAseq/Terekhova2023/GEX_HTO_processed/myeloid_cells/classical_monocytes_sce.rds"))
colnames_sex <- rownames(colData(sce_class_mono)[colData(sce_class_mono)$Gender == sex,])
md_class_mono_sex <- colData(sce_class_mono)[colData(sce_class_mono)$Gender == sex,]
sce_class_mono_sex <- sce_class_mono[, colnames_sex]

print("2. Loading senescence gene set ----------")
hallmark_ifng_response_genes.fn <- paste0(data_path, '/msopena/02_OneK1K_Age/robjects/03_Milo/03_EnrichmentISG/HALLMARK_INTERFERON_GAMMA_RESPONSE.v2024.1.Hs.gmt')
hallmark_ifna_response_genes.fn <-  paste0(data_path, '/msopena/02_OneK1K_Age/robjects/03_Milo/03_EnrichmentISG/HALLMARK_INTERFERON_ALPHA_RESPONSE.v2024.1.Hs.gmt')
hallmark_ifng_response_genes <- unlist(geneIds(getGmt(hallmark_ifng_response_genes.fn)))
hallmark_ifna_response_genes <- unlist(geneIds(getGmt(hallmark_ifna_response_genes.fn)))

## Inflammatory response
hallmark_inflammatory_response_genes.fn <- paste0(data_path, '/msopena/02_OneK1K_Age/robjects/03_Milo/03_EnrichmentISG/HALLMARK_INFLAMMATORY_RESPONSE.v2024.1.Hs.gmt')
hallmark_inflammatory_response_genes <- unlist(geneIds(getGmt(hallmark_inflammatory_response_genes.fn)))

## Select genes
hallmark_ifn_response_genes <- c(hallmark_ifng_response_genes, hallmark_ifna_response_genes)
hallmark_inflammatory_exclusive <- setdiff(hallmark_inflammatory_response_genes,hallmark_ifn_response_genes)
hallmark_ifn_exclusive <- setdiff(hallmark_ifn_response_genes, hallmark_inflammatory_response_genes)
gsea_hallmark.list <- list(hallmark_ifna_response_score = hallmark_ifna_response_genes,
                           hallmark_ifng_response_score = hallmark_ifng_response_genes,
                           hallmark_isg_score = hallmark_ifn_exclusive,
                           hallmark_inflammatory_score = hallmark_inflammatory_exclusive)


print("3.Running GSEA using AUCell ----------")

set.seed(42)
# bpparam <- MulticoreParam(progressbar=T, workers=8, log = T, stop.on.error = F)
# register(bpparam)

sce <- runEscape(sce_class_mono_sex, 
                 method = method,
                 gene.sets = gsea_hallmark.list, 
                 min.size = 0,
                 groups= 5000,
                 new.assay.name = paste0("escape.", method), BPPARAM=bpparam)



sce <- performNormalization(sce, 
                            assay = paste0("escape.", method), 
                            gene.sets = gsea_hallmark.list,  scale.factor = sce$nFeature_RNA)


scores.enrich <- t(altExp(sce)@assays@data[[paste0("escape.", method)]]) %>% as.data.frame()

mdata_scores.enrich <- merge(md_class_mono_sex, scores.enrich, by=0)

#save sce with the enrichments and the metadata 
#saveRDS(sce, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/03_EnrichmentsISG_replication/sce_classical_monocytes_Enrichment_Hallmarks_sex_",sex,".rds"  ))

saveRDS(mdata_scores.enrich, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/03_EnrichmentISG/Metadata_classical_monocytes_Enrichment_Hallmarks_sex_",s,".rds"   ))




mdata_scores.enrich_F <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/03_EnrichmentISG/Metadata_classical_monocytes_Enrichment_Hallmarks_sex_F.rds"  ))
mdata_scores.enrich_F$sex <- "Females"
mdata_scores.enrich_F <- mdata_scores.enrich_F[mdata_scores.enrich_F$cell_type == "CD14_Mono", ]
mdata_scores.enrich_M <-readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/03_EnrichmentISG/Metadata_classical_monocytes_Enrichment_Hallmarks_sex_M.rds"  ))
mdata_scores.enrich_M$sex <- "Males"
mdata_scores.enrich_M <- mdata_scores.enrich_M[mdata_scores.enrich_M$cell_type == "CD14_Mono", ]
mdata_scores.enrich <- rbind(mdata_scores.enrich_M, mdata_scores.enrich_F) %>% as.data.frame()
table(mdata_scores.enrich$sex)

mono_isg <- do.call(rbind, lapply(c("M", "F"), function(sex, celltype, gene) {
  colData(readRDS(paste0(data_path, "/aripol1/OneK1K_Age/MS_02_SubpopulationsByMarkers/" ,celltype, ".", sex, ".", gene, "_pos.SCE_Oelen.rds"))) %>% as.data.frame()
}, celltype = "CD14_Mono", gene = "ISG15"))

mdata_scores.enrich <- mdata_scores.enrich[mdata_scores.enrich$Row.names %in% mono_isg$bare_barcode_lane,]



library(glmmTMB)
model_enrichment <- function(hallmark){
  auc_df_m_term <- as.data.frame(mdata_scores.enrich[,c(hallmark, "Age","assignment", "sex", "batch")])
  colnames(auc_df_m_term)[1] <- "Score"
  auc_df_m_term$sex <- factor(auc_df_m_term$sex, levels=c("Males", "Females"))
  ecm_model <- glmmTMB(Score~sex+Age +(1 | assignment)+(1 | batch), data = auc_df_m_term)
  stats <- t(summary(ecm_model)$coefficients$cond["sexFemales", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]) %>% as.data.frame()
  colnames(stats) <- c("estimate", "std_error", "z_value", "p_value")
  rownames(stats) <- hallmark
return(stats)}

model_enrich <- lapply(c("hallmark_inflammatory_score", "hallmark_isg_score"), function(hallmark)model_enrichment(hallmark))
model_enrich_df <- do.call(rbind.data.frame, model_enrich)
model_enrich_df$fdr<- p.adjust(model_enrich_df$p_value)



saveRDS( model_enrich_df,paste0(, "/msopena/02_OneK1K_Age/robjects/03_Milo/03_EnrichmentISG/Metadata_EnrichmentSenescence_", method, "_Oelen_AllCells_DiffEnrich_.rds"))



ggplot(mdata_scores.enrich, aes(x=sex, y=hallmark_isg_score))+geom_boxplot()




