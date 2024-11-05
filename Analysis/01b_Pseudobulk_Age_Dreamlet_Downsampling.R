#!/usr/bin/env Rscript
library(ggplot2); library(dplyr); library(RColorBrewer); library(UpSetR); library("AnnotationDbi"); library('org.Hs.eg.db');library(dplyr)
library(readxl);library(Seurat);library(plyr);library(ggsignif);library(clusterProfiler);library(tidyr);library(ggpubr);library(ggrepel);library(purrr)


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

# Parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
make_option(c("--cell_type"), action="store", default=NA, type='character',
            help="Cell types in low resolution (predicted.celltype.l1) or high resolution (cell_type)"),
make_option(c("--idx"), action="store", default=NA, type='character',
            help="Number of the downsampling"))
opt = parse_args(OptionParser(option_list=option_list))

# Packages
shhh(library(Seurat))
shhh(library(MAST))
shhh(library(SingleCellExperiment))
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(ggplot2))
shhh(library(dreamlet))
shhh(library(zenith))
shhh(library(scater))
shhh(library(RColorBrewer))
Csparse_validate = "CsparseMatrix_validate"

################################## Set Variables and load Data ################################## 
cell_type <- opt$cell_type
n <- opt$idx

# Report
print(paste0('Cell type: ', cell_type))
print(paste0('R objects dir:', robjects_out))
print(paste0('Plots dir:', plots_out))

################################## Functions ################################## 

dreamlet.func <- function( ge_dge ){
  ### Defining the DEA formulas ###
  print('Defining the DEA formulas...')
  form <- as.formula(~Age+Gender+(1|date))
  #### Normalize and apply voom/voomWithDreamWeights ####
  print('Normalizing the pseudobulk-data...')
  system.time(res.proc <- processAssays(ge_dge, form, min.count=5))
  
  # View details of dropping samples
  details(res.proc)
  
  # Check nSamples and nGenes tested
  genes_all <- rownames(ge_dge)
  genes_tested <- rownames(as.data.frame(res.proc))
  genes_all.n <- nrow(ge_dge)
  genes_tested.n <- nrow(as.data.frame(res.proc))
  genes_tested.prop <- round(genes_tested.n/genes_all.n,3)
  samples_all <- colnames(ge_dge)
  samples_tested <- colnames(as.data.frame(res.proc))
  samples_all.n <- ncol(ge_dge)
  samples_tested.n <- ncol(as.data.frame(res.proc))
  samples_tested.prop <- round(samples_tested.n/samples_all.n,3)
  print(paste0('# Genes tested: ', genes_tested.n, ', out of ', genes_all.n, ' (', genes_tested.prop, ')'))
  print(paste0('# Samples tested: ', samples_tested.n, ', out of ', samples_all.n, ' (', samples_tested.prop, ')'))
  
  stats <- data.frame("Tested_genes" = genes_tested.n, "All_genes"= genes_all.n, "Tested_samples" = samples_tested.n, "All samples"= samples_all.n)
  rownames(stats) <- cell_type
  
  #save stats 
  stats_fn <- paste0(robjects_out,cell_type,'_dreamlet_stats.rds')
  saveRDS(stats, stats_fn)

  # Show voom plot for each cell clusters
  ## Here the mean-variance trend from voom is shown for each cell type. Cell types with sufficient number of cells and reads show a clear mean-variance trend. While in rare cell types like megakaryocytes, fewer genes have sufficient reads and the trend is less apparent.
  plotVoom.p <- plotVoom(res.proc)
  plotVoom.fn <- paste0(plots_out, cell_type, '_plotVoom.png')
  ggsave(plotVoom.fn, plotVoom.p)

  ### Differential expression ###
  ## Since the normalized expression data and metadata are stored within res.proc, only the regression formula remains to be specified.
  ## Here we only included the stimulus status, but analyses of larger datasets can include covariates and random effects.
  ## With formula ~ StimStatus, an intercept is fit and coefficient StimStatusstim log fold change between simulated and controls.
  ## Differential expression analysis within each assay, evaluated on the voom normalized data
  print('Running DEA...')
  system.time(res.dl <- dreamlet(res.proc, form))
  toptable <- variancePartition::topTable(res.dl[[gsub("_", " ", cell_type)]])
  ### Save outputs ###
  out <- list(dea = res.dl,
              topTable = toptable)
  out_fn <- paste0(robjects_out,cell_type,'_deaTopTable_downsampling_',n,'.rds')
  print(paste0('Saving dreamlet results: ',out_fn))
  saveRDS(out, out_fn)
  return(out)
}


################################## Analyses ################################## 

# 1. Read input data 
main_dir <- paste0(basepath, 'Data/scRNAseq/Yazar2022/sce_data_objects/')
in.fn <- paste0(main_dir, '/',cell_type, "_cell_type_sceraw.rds" )
print(paste0('Reading Seurat file in: ', in.fn))
system.time(sce <- readRDS(in.fn))

## Contrast order
Group_order.vec <- list(Age=c('Age'),
                        Gender = c('M','F'))
colData(sce)$Gender <- as.factor(colData(sce)$Gender)
colData(sce)[['Gender']] <- factor(colData(sce)[['Gender']],
                                   levels = Group_order.vec[['Gender']])

cnames <- c ('cell_type', 'assignment', 'date', 'Gender', 'Age')
colData(sce) <- colData(sce)[,cnames]

# Aggregate to pseudobulk
## Dreamlet, like muscat, performs analysis at the pseudobulk-level by summing raw counts across cells for a given sample and cell type. 
## aggregateToPseudoBulk is substantially faster for large on-disk datasets than muscat::aggregateData.
print('Performing pseudobulk...')
system.time(pb <- aggregateToPseudoBulk(sce,
                                        assay = "counts",     
                                        cluster_id = "cell_type", 
                                        sample_id = "assignment",
                                        verbose = FALSE))


# 2. Downsampling ---
# get samples that passed dreamlet filtering 
smp_df <- readRDS(paste0(data_path, "//msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/downsampling/", cell_type,'_samples_subsampling_n_', n,'.rds'))
smp <- rownames(smp_df)
pb_subset <- pb[,smp]

#pb_F <- pb_F[,sample(ncol(pb_F), ncol(pb_M))]

# 2. dreamlet
dreamlet.func(pb_subset)

