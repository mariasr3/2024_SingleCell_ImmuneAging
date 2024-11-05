#!/usr/bin/env Rscript

# setting working directory (cluster or local)
path_cluster <- '/gpfs/projects/bsc83/'
path_em <- '/home/aripol1/Desktop/bsc/'
path_opensuse <- '/home/bscuser/bsc/'

if(file.exists(path_cluster)){
  setwd(paste(path_cluster))
}else if(file.exists(path_em)){
  setwd(paste(path_em))
}else if(file.exists(path_opensuse)){
  setwd(paste(path_opensuse))
}

# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--cell_type"), action="store", default=NULL, type='character',
              help="One combination in MS_03_SubpopulationsDreamlet.celltypes.tab"),
  make_option(c("--covariates"), action="store", default="MS_03_SubpopulationsDreamlet.covariates.tab", type='character',
              help="Nhood markers"),
  make_option(c("--min_prop"), action="store", default=NULL type='double',
              help="In dreamlet::processAssays(), minimum proportion of retained samples with non-zero counts for a gene to be retained."),
  make_option(c("--vp_reduced"), action="store", default=FALSE, type='logical',
              help="Not estimate the effect of the batch factor in the VariancePartition analysis."),
  make_option(c("--in_dir"), action="store", default=NULL, type='character',
              help="Input directory."),
  make_option(c("--out_dir"), action="store", default=NULL, type='character',
              help="Output directory"))
opt = parse_args(OptionParser(option_list=option_list))

# Loading functions
main.dir <- 'Projects/scRNAseq/aripol1/OneK1K_Age/'
functions.fn <- paste0(main.dir, 'scripts/pseudobulkDEA_dreamlet_functions.R')
print(paste0('Loading functions from: ', functions.fn))
source(functions.fn)

################################## Set Variables and load Data ################################## 
# Testing
## Non-sex-stratified
opt$cell_type <- 'CD16_Mono' #CD14_Mono_sceraw.rds
opt$covariates <- 'Non_Sex_Stratified.covariates.tab'
opt$min_prop <- 0.2
opt$in_dir <- 'Data/scRNAseq/Yazar2022/sce_data_objects'
opt$out_dir <- 'Projects/scRNAseq/msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/01a_Pseudobulk_Age_Dreamlet'

## Sex-stratified
opt$cell_type <- 'CD16_Mono' #CD14_Mono_sceraw.rds
opt$covariates <- 'Sex_Stratified.covariates.tab'
opt$min_prop <- 0.4
opt$in_dir <- 'Data/scRNAseq/Yazar2022/sce_data_objects/F'
opt$out_dir <- 'Projects/scRNAseq/msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/01a_Pseudobulk_Age_Dreamlet/F'
opt$in_dir <- 'Data/scRNAseq/Yazar2022/sce_data_objects/M' #alternative
opt$out_dir <- 'Projects/scRNAseq/msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/01a_Pseudobulk_Age_Dreamlet/M' #alternative

## Downsampling
opt$cell_type <- 'CD16_Mono' #CD14_Mono_sceraw.rds
opt$covariates <- 'Non_Sex_Stratified.covariates.tab'
opt$min_prop <- 0.4 #check
opt$in_dir <- 'Data/scRNAseq/Yazar2022/sce_data_objects/Downsampling'
opt$out_dir <- 'Projects/scRNAseq/msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/01a_Pseudobulk_Age_Dreamlet/Downsampling'

## Subpopulations by markers
opt$cell_type <- 'B_naive.APOD_pos' #B_naive.APOD_pos_sceraw.rds
# opt$cell_type <- 'B_intermediate.APOD_pos' #B_intermediate.APOD_pos_sceraw.rds
# opt$cell_type <- 'CD14_Mono.ISG15_pos' #CD14_Mono.ISG15_pos_sceraw.rds
# opt$cell_type <- 'CD8TEM.GZMB_pos' #CD8TEM.GZMB_pos
opt$covariates <- 'Sex_Stratified.covariates.tab'
opt$min_prop <- 0.4
opt$in_dir <- 'Data/scRNAseq/Yazar2022/sce_data_objects/Subpopulations/F'
opt$out_dir <- 'Projects/scRNAseq/msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/01a_Pseudobulk_Age_Dreamlet/Subpopulations/F'
opt$in_dir <- 'Data/scRNAseq/Yazar2022/sce_data_objects/Subpopulations/M' #alternative
opt$out_dir <- 'Projects/scRNAseq/msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/01a_Pseudobulk_Age_Dreamlet/Subpopulations/M' #alternative

# Directories
in.dir <- paste0(main.dir, '/', opt$in_dir, '/')
out.dir <- paste0(main.dir, '/', opt$out_dir, '/', 
                  opt$cell_type, '/', 'min_prop_', as.character(opt$min_prop), '/')
if(opt$vp_reduced){out.dir <- paste0(out.dir, 'vp_reduced/')}
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

# Variables
covs.fn <- paste0(main.dir, 'scripts/', opt$covariates)
sce.fn <- paste0(in.dir, opt$cell_type, '.SCE.rds')
celltype <- str_split_fixed(opt$cell_type, '\\.', 3)[,1]
sex <- str_split_fixed(opt$cell_type, '\\.', 3)[,2]
gene <- str_split_fixed(opt$cell_type, '\\.', 3)[,3]

################################## Analyses ##################################
# Read covariates file
print(paste0('Reading covariates file in: ',covs.fn))
covs.df <- read.table(covs.fn, header = T)
covariates <- covs.df$covariate
phenotypes <- covs.df[covs.df$type=='fixed',]$covariate

# Read SCE object file
print(paste0('Reading SCE object file in: ',sce.fn))
system.time(sce <- readRDS(sce.fn))

# Save SCE object metadata
sce_md <- as.data.frame(colData(sce))
sce_md.fn <- paste0(out.dir, 'sce_metadata.rds')
print(paste0('Saving SCE object metadata file in: ',sce_md.fn))
saveRDS(sce_md, sce_md.fn)

# Get contrasts according to the phenotypes
contrast_coefName.list <- sapply(phenotypes, function(i) get_coefName(i, sce), simplify = FALSE)
contrast_coefName.list <- Filter(Negate(is.null), contrast_coefName.list)

# Print report
print(paste0('Cell type (manual): ', opt$cell_type))
print(paste0('Cell type: ', celltype))
print(paste0('Sex: ', sex))
print(paste0('Gene: ', gene))
print(paste0('Min prop: ', as.character(opt$min_prop)))
print(paste0('Covariates file: ', covs.fn))
print(paste0('nCells: ', ncol(sce)))
print(paste0('Input file: ', sce.fn))
print(paste0('Input directory: ', in.dir))
print(paste0('Output directory: ', out.dir))
cat('\n')

# Aggregate to pseudobulk: Dreamlet, like muscat, performs analysis at the pseudobulk-level by summing raw counts across cells for a given sample and cell type. 
# aggregateToPseudoBulk is substantially faster for large on-disk datasets than muscat::aggregateData.
print('Performing pseudobulk...')
colData(sce)$cell_type <- opt$cell_type
colData(sce)$cell_type <- as.factor(colData(sce)$cell_type)
system.time(pb <- aggregateToPseudoBulk(sce,
                                        assay = "counts",     
                                        cluster_id = "cell_type", 
                                        sample_id = "assignment",
                                        verbose = FALSE))
pb_fn <- paste0(out.dir, 'aggregateToPseudoBulk.rds')
print(paste0('Saving aggregateToPseudoBulk results: ', pb_fn))
saveRDS(pb, pb_fn)
pb_raw <- pb
cat('\n')

# DEA and VariancePartition with dreamlet 
# Testing
# pb_raw <- pb
# set.seed(123)
# genes <- sample(rownames(pb), 5000)
# pb <- pb[genes,]

print('Running dreamlet...')
system.time(dreamlet.res <- dreamlet.func(ge_dge = pb,
                                          covariates = covs.df, 
                                          min_prop = opt$min_prop,
                                          contrast_list = contrast_coefName.list, 
                                          vp_reduced = opt$vp_reduced,
                                          out_dir = out.dir))
cat('\n')

# Check results
print('### Checking results ###')
print('All results (top)...')
print(dreamlet.res$topTable)
cat('\n')
print('Significant (FDR) --> adj.P.Val < 0.05')
lapply(dreamlet.res$topTable, function(x) as.data.frame(x[x$adj.P.Val<=0.05,]))
cat('\n')
print('Significant (nominal) --> P.Value < 0.01')
lapply(dreamlet.res$topTable, function(x) as.data.frame(x[x$P.Value<=0.01,]))
