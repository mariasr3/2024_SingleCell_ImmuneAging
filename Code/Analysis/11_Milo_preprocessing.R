

# Run pre-processing for the single-cell differential abundance analysis with Milo 


#libraries
library(miloR); library(SingleCellExperiment); library(Seurat);library(scran); library(dplyr); library(patchwork);library(scater);library(ggplot2);library(BiocParallel)
Csparse_validate = "CsparseMatrix_validate"

#paths
basepath <- "/gpfs/projects/bsc83/"
data_path <- "/gpfs/projects/bsc83/Projects/scRNAseq/"

# arguments
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--num_rep"), action="store", default=NA, type='character',
              help="2-4"),
  make_option(c("--PerSex"), action="store", default=NA, type='boolean',
              help="TRUE/FALSE"),
  make_option(c("--sex"), action="store", default=NA, type='boolean',
              help="M, F, NA")
)
opt = parse_args(OptionParser(option_list=option_list))

n <- opt$num_rep
sex <- opt$sex

#functions
source(paste0(data_path, "/msopena/02_OneK1K_Age/scripts/Milo_functions.R"))

if(opt$PerSex){
  print(paste0('Pre-processing for DA testing with age for sex: ', sex))
}else{
  print(paste0('Pre-processing for DA testing with age for all data, subsampling: ', n))
}

# 1. Extract SingleCellExperiment object
print("0. Reading Seurat object ------")
if(opt$PerSex){
  sce <- readRDS(paste0(basepath, "/Data/scRNAseq/Yazar2022/sce_data_objects/AllCells_0.5_sce_harmony_",sex,".rds"))
}else{
  sce <- readRDS(paste0(basepath, "/Data/scRNAseq/Yazar2022/sce_data_objects/AllCells_0.25_predicted.celltype.l1_sce_harmony_",n,".rds"))
}

# Extract sample metadata to use for testing
print("1. Building MILO object ------")
traj_milo <- Milo(sce)
print("Milo object created ---")

# 2. Construct KNN graph ---
print("2. Building KNN graph ------")
traj_milo <- buildGraph(traj_milo, k = 100, d = 30, reduced.dim = "HARMONY")

# 3. Defining representative neighbourhoods
  #prop: the proportion of cells to randomly sample to start with (0.05 is sufficient for datasets with > 10k cells)
  # k: the k to use for KNN refinement (we recommend using the same k used for KNN graph building)
  # d: the number of reduced dimensions to use for KNN refinement (standard practices is d= 30)
# refined indicated whether you want to use the sampling refinement algorithm, or just pick cells at random (set as T)

print("3. Creating neighbourhoods ------")
traj_milo <- makeNhoods(traj_milo, prop = 0.05, k = 100, d=30, refined = TRUE, reduced_dims = "HARMONY", refinement_scheme = "graph" )

print("4. Counting cells per neighbourhood ---")
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), samples="assignment")

print("5. Building nhood graph ----")
traj_milo <- buildNhoodGraph(traj_milo, overlap=0)

print("6. Calculating per-neighbourhood gene expression ----")
traj_milo <- calcNhoodExpression(traj_milo)

if(opt$PerSex){
  saveRDS(traj_milo, paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/04_MiloObject_0.5_Preprocessed_",sex,".rds"))
}else{
  saveRDS(traj_milo, paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/04_MiloObject_Preprocessed_Subsampling_0.25_",n,".rds"))
  }


print("Milo pre-processing run and object saved")


