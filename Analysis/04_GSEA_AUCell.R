
library(AUCell);library(dplyr);library(GSEABase);library(stringr);library(fitdistrplus);library(Seurat);library(SingleCellExperiment); #library(escape)

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
path_senec <- paste0(robjects, "/09_Senescence/")
#dir.create(path_enrichScore, recursive = TRUE)

#functions
source(paste0(data_path, "/msopena/02_OneK1K_Age/scripts/functions.R"))

shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--method"), action="store", default=NA, type='character',
              help="AUCell, UCell"),
  make_option(c("--geneset"), action="store", default=NA, type='character',
              help="Quiescence, Senescence"),
  make_option(c("--GeneSetFile"), action="store", default=NA, type='character',
              help="path to the gene set directory"))
opt = parse_args(OptionParser(option_list=option_list))

method <- opt$method 
GeneSet <- opt$geneset 
gmtFile <- opt$GeneSetFile
metadata <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/metadata_processed.rds"))

print(paste0("##### Computing Enrich GSEA using ",  method, " for all cells #####"))

# 1. Get expression data 
print("1. Loading data ----------")
#sce <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/sce_data_objects/", cell_type, "_cell_type_sceraw.rds"))
sce <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/sce_data_objects/AllCells_predicted.celltype.l1_sceraw.rds"))


print("2. Loading senescence gene set ----------")
print(paste0("File geneset: ", gmtFile))
#gmtFile <- "/home/mariasr/cluster/Projects/scRNAseq/msopena/02_OneK1K_Age/robjects/08_EnrichmentScores/GeneSets/SenMayo_GeneSet.gmt"
geneSets <- getGmt(gmtFile)

print("3.Running GSEA using UCell ----------")
set.seed(42)
bpparam <- MulticoreParam(progressbar=T, workers=8, log = T, stop.on.error = F)
register(bpparam)

sce <- runEscape(sce, 
                 method = method,
                 gene.sets = geneSets, 
                 min.size = 0,
                 groups= 5000,
                 new.assay.name = paste0("escape.", method), BPPARAM=bpparam)

sce <- performNormalization(sce, 
                            assay = paste0("escape.", method), 
                            gene.sets = geneSets,  scale.factor = sce$nFeature_RNA)


scores.enrich <- t(altExp(sce)@assays@data[[paste0("escape.", method)]]) %>% as.data.frame()

mdata_scores.enrich <- merge(metadata, scores.enrich, by=0)

saveRDS(sce, paste0(path_senec, "sce_Enrichment", GeneSet,"_",method,"_AllCells.rds"))

saveRDS(mdata_scores.enrich, paste0(path_senec, "Metadata_Enrichment",GeneSet,"_",method,"_AllCells.rds"))


mdata_scores.enrich <- readRDS(paste0(path_senec, "Metadata_Enrichment",GeneSet,"_",method,"_AllCells.rds"))
mdata_scores.enrich$celltype <- mdata_scores.enrich$cell_type
#mdata_scores.enrich <- mdata_scores.enrich[mdata_scores.enrich$celltype %in% order_cells$cell_type[1:19],]
mdata_scores.enrich<- reorder_cells(mdata_scores.enrich, reverse = T)
mdata_scores.enrich$Age_cat <- factor(mdata_scores.enrich$Age_cat, levels=c("O", "M","Y"))
ggplot(mdata_scores.enrich, aes(x=celltype, y=SAUL_SEN_MAYO, alpha=Age_cat))+geom_boxplot(outlier.shape = NA, fill=blue)+theme+facet_grid(celltype_l1~., scales = "free", space = "free")+coord_flip()



library(glmmTMB)
model_enrichment <- function(celltype, geneset){
  auc_df_m_term <- as.data.frame(mdata_scores.enrich[mdata_scores.enrich$celltype == celltype,c(geneset, "Age", "assignment", "Gender", "date")])
  colnames(auc_df_m_term)[1] <- "Score"
  ecm_model <- glmmTMB(Score~Age +Gender+(1 | date)+(1 | assignment), data = auc_df_m_term)
  stats <- t(summary(ecm_model)$coefficients$cond["Age", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]) %>% as.data.frame()
  colnames(stats) <- c("estimate", "std_error", "z_value", "p_value")
  rownames(stats) <- celltype
  return(stats)
}

if (GeneSet== "Quiescence"){
  model_enrich <- lapply(unique(as.character(mdata_scores.enrich$celltype)), function(celltype)model_enrichment(celltype, "SAUL_SEN_MAYO"))
}else{  model_enrich <- lapply(unique(as.character(mdata_scores.enrich$celltype)), function(celltype)model_enrichment(celltype, "QUIESCENCE_GENES"))
}
model_enrich_df <- do.call(rbind.data.frame, model_enrich)
saveRDS( model_enrich_df,paste0(path_senec, "DifferentialEnrichment_", GeneSet, "_", method, "_AllCells.rds"))
