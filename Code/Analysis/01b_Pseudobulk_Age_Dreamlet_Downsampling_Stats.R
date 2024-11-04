#!/usr/bin/env Rscript

# Paths
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
robjects_out <- paste0(data_path, "/msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/downsampling/")
plots_out <-  paste0(data_path, "/msopena/02_OneK1K_Age/plots/11_DEA_SexAge/analysis/")


# Plots results Differential Gene Expression Analysis

library(ggplot2); library(dplyr); library(RColorBrewer); library(UpSetR); library("AnnotationDbi"); library('org.Hs.eg.db');library(dplyr)
library(readxl);library(Seurat);library(plyr);library(ggsignif);library(clusterProfiler);library(tidyr);library(ggpubr);library(ggrepel);library(purrr)


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

#pal
# pal_direction <- c("up" = "#bc4749" , "down"= "#264653")
# computer <- "work"
# source(paste0(data_path, "/msopena/01_meta-analysis_Age_Sex/scripts/themes.R"))
pheno <- "Age"
cell_level <- "cell_type"

#data
metadata <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/metadata_processed.rds"))
order_cells<- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/order_cells.rds"))
#order_cells<- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/new_order_cells.rds"))
celltype_l1 <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/celltypes_equivalence.rds"))
celltype_l1$cell_type <- factor(celltype_l1$cell_type, levels = order_cells$cell_type)
mdata_donors <- metadata[!duplicated(metadata$assignment),]
saveRDS(mdata_donors, paste0(basepath, "Data/scRNAseq/Yazar2022/donor_metadata.rds"))



################# BEFORE DOWNSAMPLING ###############################

# 1. Number of samples after dreamlet filtering ######
major_cells <- order_cells$cell_type[1:13]
major_cells <- gsub(" ", "_", major_cells)

dea <- lapply(major_cells, function(cell_type) {
  print(cell_type)
  result_list <- lapply(c(1:20), function(n) {
    print(n)
    # Read the RDS file
    dea <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/downsampling/", cell_type, "_deaTopTable_downsampling_", n, ".rds"))$dea
    
    # Additional commands go here (if needed)
    dea <- dreamlet::topTable(dea, coef="Age", number = Inf)
    # n_it <- data.frame(length(dea[dea$adj.P.Val <0.05,]$ID))
    # colnames(n_it)<- cell_type
    # rownames(n_it) <- n
    return(dea)  # Return the result
  })
  return(result_list)  # Return the list of results
})

saveRDS(dea, paste0(data_path,"/msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/downsampling/AllCells_downsamplingTopTable.rds"))
