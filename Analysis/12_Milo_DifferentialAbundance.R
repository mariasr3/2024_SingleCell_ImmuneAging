
# Run single-cell differential abundance analysis with Milo 

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
  print(paste0('DA testing with age for sex: ', sex))
}else{
  print(paste0('DA testing with age for all data, subsampling: ', n))
}


# Read pre-processed milo object 
print("1. Read in milo object ----")
if(opt$PerSex){
  system.time(traj_milo <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/04_MiloObject_Preprocessed_Subsampling_0.25_",n,".rds")))
}else{
  system.time(traj_milo <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/04_MiloObject_Preprocessed_", sex,".rds")))
}
print("Object loaded")

# Create matrix design for the differential abundance analysis
print("2. Create matrix design  ----")
if(opt$PerSex){
  traj_design <- data.frame(colData(traj_milo))[,c("Age", "assignment", "date")]
  traj_design <- distinct(traj_design)
  traj_design$Age_binned <- ifelse(traj_design$Age<=40, 'Y', ifelse(traj_design$Age>=60, 'O', 'M'))
  rownames(traj_design) <- traj_design$assignment
  saveRDS(traj_design, paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/04_TrajDesign_Sex_",sex,".rds"))
}else{
  traj_design <- data.frame(colData(traj_milo))[,c("Age","Gender", "assignment", "date")]
  traj_design <- distinct(traj_design)
  #traj_design$Age_binned <- as.integer(paste0(substr(traj_design$Age, 1, 1), "0"))
  traj_design$Age_binned <- ifelse(traj_design$Age<=40, 'Y', ifelse(traj_design$Age>=60, 'O', 'M'))
  rownames(traj_design) <- traj_design$assignment
  saveRDS(traj_design, paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/04_TrajDesign_Subsampling_0.25_",n,".rds"))
}

# Run differential abundance analysis 
print("3.Differential abundance testing---")
set.seed(42)
bpparam <- MulticoreParam(progressbar=T, workers=8, log = T, stop.on.error = F)
register(bpparam)

# If we run the analysis per sex the variable "Gender" won't be included into the model 
if(opt$PerSex){
#test_NHoods_customed just adds a line of code to the funcion that saves the list with the GLMM results without buliding the final dataframe
system.time(da_results <- testNhoods_per_sex(traj_milo, design = ~ Age +(1|assignment) + (1|date), design.df = traj_design,sex=sex, fdr.weighting="graph-overlap",
                         glmm.solver="HE-NNLS", REML=TRUE, norm.method="TMM", BPPARAM = bpparam, fail.on.error=FALSE))
saveRDS(da_results,paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DifferentialAbundance_results_Sex_",sex,".rds"))
}else{
  system.time(da_results <- testNhoods_customed(traj_milo, design = ~ Age + Gender+(1|assignment) + (1|date), design.df = traj_design,n=n, fdr.weighting="graph-overlap",
                                                glmm.solver="HE-NNLS", REML=TRUE, norm.method="TMM", BPPARAM = bpparam, fail.on.error=FALSE))
  saveRDS(da_results,paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DifferentialAbundance_results_newPreprocessing_Subsampling_0.25_",n,".rds"))
  
}

print("4.Building output table --")
#Workararound when testNhoods does not work proprely. Basically we read the GLMM output and constructe the dataframe ourselves becasue the function was crashing 
if(opt$PerSex){
traj_milo <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/04_MiloObject_Preprocessed_Subsampling_0.25_Sex",sex,".rds"))
fit <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/Fit_cell_type_annotation_NewPreprocessing_Subsampling_0.25_Sex_",sex,".rds"))
}else{traj_milo <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/04_MiloObject_Preprocessed_Subsampling_0.25_",n,".rds"))
fit <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/Fit_cell_type_annotation_NewPreprocessing_Subsampling_0.25_",n,".rds"))}

n.nhoods <- length(fit)
half.n <- floor(n.nhoods * 0.5)
ret.beta <- 2
rand.levels <- 2
#names(rand.levels) <- c("date", "assignment")
res <- cbind.data.frame("logFC" = unlist(lapply(lapply(fit, `[[`, "FE"), function(traj_milo) traj_milo[ret.beta])),
                        "logCPM"=log2((rowMeans(nhoodCounts(traj_milo)/colSums2(nhoodCounts(traj_milo))))*1e6),
                        "SE"= unlist(lapply(lapply(fit, `[[`, "SE"), function(traj_milo) traj_milo[ret.beta])),
                        "tvalue" = unlist(lapply(lapply(fit, `[[`, "t"), function(traj_milo) traj_milo[ret.beta])),
                        "PValue" = unlist(lapply(lapply(fit, `[[`, "PVALS"), function(traj_milo) traj_milo[ret.beta])),
                        unlist(lapply(lapply(fit, `[[`, "Sigma"), function(traj_milo) traj_milo[rand.levels-1])),
                        unlist(lapply(lapply(fit, `[[`, "Sigma"), function(traj_milo) traj_milo[rand.levels])),
                        "Converged"=unlist(lapply(fit, `[[`, "converged")), "Dispersion" = unlist(lapply(fit, `[[`, "Dispersion")),
                        "Logliklihood"=unlist(lapply(fit, `[[`, "LOGLIHOOD")))

rownames(res) <- 1:length(fit)
colnames(res)[6:(6+length(rand.levels)-1)] <- paste(names(rand.levels), "variance")
res$Nhood <- as.numeric(rownames(res))
fdr.weighting <- "graph-overlap"
message("Performing spatial FDR correction with ", fdr.weighting[1], " weighting")
# res1 <- na.omit(res)
mod.spatialfdr <- graphSpatialFDR(x.nhoods=nhoods(traj_milo),
                                  graph=graph(traj_milo),
                                  weighting=fdr.weighting,
                                  k=traj_milo@.k,
                                  pvalues=res[order(res$Nhood), ]$PValue,
                                  indices=nhoodIndex(traj_milo),
                                  distances=nhoodDistances(traj_milo),
                                  reduced.dimensions=reducedDim(x, reduced.dim))
res$SpatialFDR[order(res$Nhood)] <- mod.spatialfdr
da_results <-  res

if(opt$PerSex){
 saveRDS(da_results,paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DifferentialAbundance_results_newPreprocessing_Subsampling_0.25_Sex_",sex,".rds"))
}else{ saveRDS(da_results,paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DifferentialAbundance_results_newPreprocessing_Subsampling_0.25_",n,".rds"))
}

print("5.Annotate cell types to DAA results ---")
system.time(da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "cell_type"))
da_results[is.na(da_results$SpatialFDR),]$SpatialFDR <- 1
if(opt$PerSex){
saveRDS(da_results,paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DAA_results_cell_type_annotation_NewPreprocessing_Subsampling_0.25_Sex_",sex,".rds"))
}else{
saveRDS(da_results,paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DAA_results_cell_type_annotation_NewPreprocessing_Subsampling_0.25_",n,".rds"))
}

print("6. Group nhoods fom the DAA results ---")
system.time(da_results <- groupNhoods(traj_milo, da_results, da.fdr = 0.1,  max.lfc.delta = 10))
saveRDS(da_results,paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DAA_results_cell_type_annotation_nhood_grouped_NewPreprocessing_Subsampling_0.25_Sex_",sex,"_log10.rds"))
saveRDS(da_results,paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DAA_results_cell_type_annotation_nhood_grouped_NewPreprocessing_Subsampling_0.25_",n,".rds"))

