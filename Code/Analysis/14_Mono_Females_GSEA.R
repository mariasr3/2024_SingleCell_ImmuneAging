library(AUCell);library(dplyr);library(GSEABase);library(stringr);library(fitdistrplus);library(Seurat);library(SingleCellExperiment); library(escape)

Csparse_validate = "CsparseMatrix_validate"

# Get enrichment analysis for senescent genes using SenMayo from Saul et al. 2022 (https://www.nature.com/articles/s41467-022-32552-1)

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
  make_option(c("--direction"), action="store", default=NA, type='character',
              help="up, down"), 
  make_option(c("--signif"), action="store", default=NA, type='character',
              help="ns, ss"), 
  make_option(c("--gene_set"), action="store", default=NA, type='character',
              help="hallmark"))

opt = parse_args(OptionParser(option_list=option_list))

sex <- opt$sex
direction <- opt$direction
signif <- opt$signif
gene_set_name <- opt$gene_set

celltype <- "CD14_Mono"



method <- "AUCell"

metadata <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/metadata_processed.rds"))

print(paste0("##### Computing Enrich GSEA using ",  method, " for all cells #####"))


# 1. Get expression data 
print("1. Loading data ----------")

so <- readRDS(paste0(basepath, "Projects/scRNAseq/msopena/02_OneK1K_Age/robjects/03_Milo/01_Milo_NhoodMarkers/so_CD14_Mono_down_up_ns_Nhoods_Sex_M.rds"))

sce <- as.SingleCellExperiment(so)

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

sce <- runEscape(sce, 
                 method = method,
                 gene.sets = gsea_hallmark.list, 
                 min.size = 0,
                 groups= 5000,
                 new.assay.name = paste0("escape.", method), BPPARAM=bpparam)



sce <- performNormalization(sce, 
                            assay = paste0("escape.", method), 
                            gene.sets = gsea_hallmark.list,  scale.factor = sce$nFeature_RNA)


scores.enrich <- t(altExp(sce)@assays@data[[paste0("escape.", method)]]) %>% as.data.frame()

mdata_scores.enrich <- merge(metadata, scores.enrich, by=0)

#gene_set_prefix <-  strsplit(gene_set_name, "\\.")[[1]][1]

saveRDS(sce, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/03_EnrichmentISG/sce_CD14_Mono_Enrichment_Hallmarks_up_down_ns_Nhoods_SexM.rds"  ))

saveRDS(mdata_scores.enrich, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/03_EnrichmentISG/Metadata_CD14_Mono_Enrichment_Hallmarks_up_down_ns_Nhoods_SexM.rds"  ))


#Females 
# 1. Get expression data 
sex <- "F"
direction <- "up"
signif <- "ss"


print("1. Loading data ----------")
#sce <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/sce_data_objects/", cell_type, "_cell_type_sceraw.rds"))
so <- readRDS(paste0(basepath, "Projects/scRNAseq/msopena/02_OneK1K_Age/robjects/03_Milo/01_Milo_NhoodMarkers/so_CD14_Mono_up_ss_Nhoods_Sex_F.rds"))

sce <- as.SingleCellExperiment(so)

print("2. Loading senescence gene set ----------")
hallmark_ifng_response_genes.fn <- paste0(data_path, '/msopena/02_OneK1K_Age/robjects/03_Milo/03_EnrichmentISG/HALLMARK_INTERFERON_GAMMA_RESPONSE.v2024.1.Hs.gmt')
hallmark_ifna_response_genes.fn <-  paste0(data_path, '/msopena/02_OneK1K_Age/robjects/03_Milo/03_EnrichmentISG/HALLMARK_INTERFERON_ALPHA_RESPONSE.v2024.1.Hs.gmt')
hallmark_ifng_response_genes <- unlist(geneIds(getGmt(hallmark_ifng_response_genes.fn)))
hallmark_ifna_response_genes <- unlist(geneIds(getGmt(hallmark_ifna_response_genes.fn)))

## Inflammatory response
hallmark_inflammatory_response_genes.fn <- paste0(data_path, '/msopena/02_OneK1K_Age/robjects/03_Milo/03_EnrichmentISG/HALLMARK_INFLAMMATORY_RESPONSE.v2024.1.Hs.gmt')
hallmark_inflammatory_response_genes <- unlist(geneIds(getGmt(hallmark_inflammatory_response_genes.fn)))

## Select genes
hallmark_ifn_response_genes <- intersect(hallmark_ifng_response_genes, hallmark_ifna_response_genes)
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

sce <- runEscape(sce, 
                 method = method,
                 gene.sets = gsea_hallmark.list, 
                 min.size = 0,
                 groups= 5000,
                 new.assay.name = paste0("escape.", method), BPPARAM=bpparam)



sce <- performNormalization(sce, 
                            assay = paste0("escape.", method), 
                            gene.sets = gsea_hallmark.list,  scale.factor = sce$nFeature_RNA)


scores.enrich <- t(altExp(sce)@assays@data[[paste0("escape.", method)]]) %>% as.data.frame()

mdata_scores.enrich <- merge(metadata, scores.enrich, by=0)

#gene_set_prefix <-  strsplit(gene_set_name, "\\.")[[1]][1]
saveRDS(mdata_scores.enrich, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/03_EnrichmentISG/Metadata_CD14_Mono_Enrichment_Hallmarks_up_ss_Nhoods_SexF.rds"  ))
#saveRDS(sce, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/03_EnrichmentISG/sce_CD14_Mono_Enrichment_Hallmarks_up_ss_Nhoods_SexF.rds"  ))

mdata_scores.enrich_F <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/03_EnrichmentISG/Metadata_CD14_Mono_Enrichment_Hallmarks_up_ss_Nhoods_SexF.rds"  ))
mdata_scores.enrich_F$sex <- "Females"
mdata_scores.enrich_M <-readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/03_EnrichmentISG/Metadata_CD14_Mono_Enrichment_Hallmarks_up_down_ns_Nhoods_SexM.rds"  ))
mdata_scores.enrich_M$sex <- "Males"
mdata_scores.enrich <- rbind(mdata_scores.enrich_M, mdata_scores.enrich_F)


da_list_F <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Nhood-cell_id_df_Sex_F.rds" )) %>% dplyr::filter(cell_type_nhood== "CD14 Mono")
da_list_M <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Nhood-cell_id_df_Sex_M.rds" )) %>% dplyr::filter(cell_type_nhood== "CD14 Mono")

da_results_M <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DAA_results_cell_type_annotation__Sex_0.5_M.rds")) %>% dplyr::filter(cell_type=="CD14 Mono" & logFC < 0)
da_results_M$sex <- "Males"
da_results_F <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DAA_results_cell_type_annotation__Sex_0.5_F.rds")) %>% dplyr::filter(cell_type=="CD14 Mono" & SpatialFDR < 0.05)
da_results_F$sex <- "Females" 
da_results <- rbind(da_results_F, da_results_M)

da_list_F <- da_list_F[da_list_F$Nhood %in% da_results_F$Nhood,]
da_list_M <- da_list_M[da_list_M$Nhood %in% da_results_M$Nhood,]

mdata_scores.enrich <- mdata_scores.enrich[mdata_scores.enrich$Row.names %in% c(da_list_F$cell_id, da_list_M$cell_id),]

da_list <- rbind(da_list_F, da_list_M)

mdata_scores.enrich <- mdata_scores.enrich %>% left_join(da_list, by=c("Row.names"="cell_id"))


table(mdata_scores.enrich$sex)


library(glmmTMB)
model_enrichment <- function(hallmark){
  auc_df_m_term <- as.data.frame(mdata_scores.enrich[,c(hallmark, "Age", "assignment", "sex", "date", "Nhood")])
  colnames(auc_df_m_term)[1] <- "Score"
  auc_df_m_term$sex <- factor(auc_df_m_term$sex, levels=c("Males", "Females"))
  auc_df_m_term$date <- droplevels(auc_df_m_term$date)
  ecm_model <- glmmTMB(Score~sex+Age +(1|Nhood)+(1 | date)+(1 | assignment), data = auc_df_m_term)
  stats <- t(summary(ecm_model)$coefficients$cond["sexFemales", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]) %>% as.data.frame()
  colnames(stats) <- c("estimate", "std_error", "z_value", "p_value")
  rownames(stats) <- hallmark
return(stats)}

model_enrich <- lapply(c("hallmark_inflammatory_score", "hallmark_isg_score"), function(hallmark)model_enrichment(hallmark))
model_enrich_df <- do.call(rbind.data.frame, model_enrich)
model_enrich_df$fdr<- p.adjust(model_enrich_df$p_value)

saveRDS( model_enrich_df,paste0(path_senec, "Metadata_EnrichmentSenescence_", method, "_AllCells_DiffEnrich.rds"))








