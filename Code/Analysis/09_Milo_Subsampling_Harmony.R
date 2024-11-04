#Run Harmony
library(Seurat); library(SingleCellExperiment);library(dplyr)
Csparse_validate = "CsparseMatrix_validate"

# paths --
basepath <- "/gpfs/projects/bsc83/"

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
if(opt$PerSex){
  print(paste0('Subset and run harmony per sex for: ', sex))
}else{
  print(paste0('Subset and run harmony per all data, subsampling number: ', n))
}

if (opt$PerSex == T){
  # read sce object ----
  sce <- readRDS(paste0(basepath, "/Data/scRNAseq/Yazar2022/sce_data_objects/AllCells_",sex,".rds"))
  
}else{
  sce <- readRDS(paste0(basepath, "/Data/scRNAseq/Yazar2022/sce_data_objects/AllCells_cell_type_sceraw.rds"))
}
print("sce loaded")

#get the age classes to maintain the same distibution 
mdata <- colData(sce)
mdata$Age_cat <- ifelse(mdata$Age > 60, "O", ifelse(mdata$Age < 40, "Y", "M"))
mdata_donor <- mdata[!duplicated(mdata$assignment),]
age_class <- split(mdata_donor$assignment, mdata_donor$Age_cat)

print("Subsetting data to 25% ----")
extract_25_percent <- function(vec) {
  n <- length(vec)
  if(opt$PerSex==T){
    # if we do it per sex, the subsampling will be at 50% of each dataset (we had already split the data )
    sample_size <- ceiling(0.5 * n)
  }else{
    # if not we split by 25% 
    sample_size <- ceiling(0.25 * n)
  }
  sampled_elements <- sample(vec, sample_size)
  return(sampled_elements)
}

# Apply the function to each element of the list
samples_25 <- as.vector(unlist(lapply(age_class, extract_25_percent)))
cells_25 <- rownames(mdata[mdata$assignment %in% samples_25,])
sce_25 <- sce[,cells_25]
colData(sce_25) <- colData(sce_25)[cells_25,]

if(opt$PerSex==T){
saveRDS(sce_25, paste0(basepath, "/Data/scRNAseq/Yazar2022/sce_data_objects/AllCells_0.5_",sex,".rds"))
}else{
  saveRDS(sce_25, paste0(basepath, "/Data/scRNAseq/Yazar2022/sce_data_objects/AllCells_0.25_predicted.celltype.l1_sce_raw_",n,".rds"))
}
# 

# pre-process the new subsampled object with the standard Seuarat pipeline 
print('Object preprocessing ----')
so <- CreateSeuratObject(counts=counts(sce_25), meta.data=as.data.frame(colData(sce_25)))
so <- NormalizeData(so)
so <- FindVariableFeatures(so)
so <- ScaleData(so)
so <- RunPCA(so)
so <- FindNeighbors(so, dims = 1:50)
so <- FindClusters(so, resolution = 2)
so <- AddMetaData(so, as.data.frame(colData(sce_25)))

if(opt$PerSex==T){
saveRDS(so, paste0(basepath, "/Data/scRNAseq/Yazar2022/sce_data_objects/AllCells_0.5_so_preprocessed_",sex,".rds"))
}else{
  saveRDS(so, paste0(basepath, "/Data/scRNAseq/Yazar2022/sce_data_objects/AllCells_0.25_predicted.celltype.l1_so_preprocessed_",n,".rds"))
}

# Run batch correction method to remove batch effects 
library(harmony)
print('Running harmony----')
so <- RunHarmony(so, "date")
sce_25 <- as.SingleCellExperiment(so)


if(opt$PerSex==T){
  saveRDS(so, paste0(basepath, "/Data/scRNAseq/Yazar2022/sce_data_objects/AllCells_0.5_sce_harmony_",sex,".rds"))
}else{
  saveRDS(sce_25, paste0(basepath, "/Data/scRNAseq/Yazar2022/sce_data_objects/AllCells_0.25_predicted.celltype.l1_sce_harmony_",n,".rds"))
}
print("Data preprocessed and object saved")




