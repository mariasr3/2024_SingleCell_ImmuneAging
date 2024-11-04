
# Explore SenCID results 

#packages
library(ggplot2); library(dplyr); library(RColorBrewer); library(UpSetR); library("AnnotationDbi"); library('org.Hs.eg.db');library(patchwork);
library(readxl);library(Seurat);library(plyr);library(ggsignif);library(clusterProfiler);library(tidyr);library(ggpubr);library(ggrepel);
library(scales)

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
plots_path <- paste0(data_path, "msopena/02_OneK1K_Age/plots/")
path_deg <- paste0(plots_path, "/01_DEG/")
dir.create(path_deg, recursive = TRUE)

#functions
source(paste0(data_path, "/msopena/02_OneK1K_Age/scripts/functions.R"))

#metadata
metadata <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/metadata_processed.rds"))
md_cols <- metadata[, c("Gender", "Age", "assignment", "date", "cell_type", "Age_cat")]
md_cols$cell_id <- rownames(md_cols)
order_cells_old<- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/order_cells.rds"))
order_cells<- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/new_order_cells.rds"))

celltype_l1 <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/celltypes_equivalence.rds"))
celltype_l1$cell_type <- factor(celltype_l1$cell_type, levels = order_cells$cell_type)


#read cell type list 
celltypes <- read.table(paste0(data_path, "msopena/02_OneK1K_Age/scripts/Analysis/cell_type.tab"))$V1

get_preditctions <- function(celltype){
  print(celltype)
  sen_path <- paste0(data_path, "/msopena/02_OneK1K_Age/robjects/14_SenCID/done/output/RecommendSID_Index_", celltype,"_cell_type_CountMatrix.txt")
  if(file.exists(sen_path)){
    # extract which model is the best prediction
    sen <- read.table(sen_path)
    sid <- names(table(sen$RecSID))[which.max(table(sen$RecSID))]
    
    # load the predictions
    predictions <- read.table(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/14_SenCID/done/output/predictions_", sid,"_",celltype,"_cell_type_CountMatrix.txt"))
    predictions$cell_id <- rownames(predictions)
    predictions <- merge(predictions, md_cols, by="cell_id")
    predictions$celltype <- celltype 
    return(predictions)
  }else{
    return(NULL)
  }
}


preds <- lapply(celltypes, function(cell) get_preditctions(cell))
predictions <- do.call(rbind.data.frame, preds)
predictions$celltype <- predictions$cell_type
predictions <- reorder_cells(predictions, neworder = T)
saveRDS(predictions, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/14_SenCID/Final_Matrix.rds"))

# model senescence score with age 

library(glmmTMB)
model_enrichment <- function(celltype){
  print(celltype)
  auc_df_m_term <- as.data.frame(predictions[predictions$celltype == celltype,c("SID_Score", "Age", "assignment", "Gender", "date")])
  colnames(auc_df_m_term)[1] <- "Score"
  ecm_model <- glmmTMB(Score~Age +Gender+(1 | date)+(1 | assignment), data = auc_df_m_term)
  stats <- t(summary(ecm_model)$coefficients$cond["Age", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]) %>% as.data.frame()
  colnames(stats) <- c("estimate", "std_error", "z_value", "p_value")
  rownames(stats) <- celltype
  return(stats)
}

stats <- lapply(unique(predictions$celltype), function(cell) model_enrichment(cell))
stats_df <- do.call(rbind.data.frame, stats)
stats_df$fdr <- p.adjust(stats_df$p_value, method = "fdr")
table(stats_df$fdr < 0.05)
stats_df$logpval <- -log10(stats_df$fdr)
stats_df$significance <- ifelse(stats_df$fdr < 0.05, "ss", "ns")
stats_df[is.na(stats_df$significance),]$significance <-"ns"
stats_df$celltype <- rownames(stats_df)

saveRDS(stats_df, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/14_SenCID/GLMM_model_results.rds"))
