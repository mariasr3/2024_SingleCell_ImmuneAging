
# Run single-cell differential abundance analysis with Milo 

#libraries
library(miloR); library(SingleCellExperiment); library(Seurat);library(scran); library(dplyr); library(patchwork);library(scater);library(ggplot2);library(BiocParallel)
Csparse_validate = "CsparseMatrix_validate"

#paths
basepath <- "/gpfs/projects/bsc83/"
data_path <- "/gpfs/projects/bsc83/Projects/scRNAseq/"


shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--k"), action="store", default=NA, type='character',
              help="60-100"))
opt = parse_args(OptionParser(option_list=option_list))

sex <- opt$sex
k <- opt$k
print(paste0('DA testing with age for: ', sex))

source(paste0(data_path, "/msopena/02_OneK1K_Age/scripts/Milo_functions.R"))




# Extract SingleCellExperiment object
print("0. Reading Seurat object ------")
so <- readRDS(paste0(basepath, "/Data/scRNAseq/Yazar2022/sce_data_objects/AllCells_0.25_predicted.celltype.l1_sce_harmony_2.rds"))
sce <- as.SingleCellExperiment(so)

## Extract sample metadata to use for testing
print("1. Building MILO object ------")
#sce <- as.SingleCellExperiment(so)
traj_milo <- Milo(sce)
#saveRDS(traj_milo, paste0(data_path,  "/msopena/03_Menopause/robjects/03_Milo/01_MiloObject_Subsampling_0.25_",n,".rds"))
print("Milo object created ---")

# 2. Construct KNN graph ---
print("2. Building KNN graph ------")
#traj_milo <- readRDS( paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/02_MiloObject_Subsampling_0.25_",n,".rds"))
traj_milo <- buildGraph(traj_milo, k = k, d = 30, reduced.dim = "HARMONY")

# 3. Defining representative neighbourhoods
  #prop: the proportion of cells to randomly sample to start with (0.05 is sufficient for datasets with > 10k cells)
  # k: the k to use for KNN refinement (we recommend using the same k used for KNN graph building)
  # d: the number of reduced dimensions to use for KNN refinement (standard practices is d= 30)
  # refined indicated whether you want to use the sampling refinement algorithm, or just pick cells at random (set as T)

print("3. Creating neighbourhoods ------")
traj_milo <- makeNhoods(traj_milo, prop = 0.05, k = k, d=30, refined = TRUE, reduced_dims = "HARMONY", refinement_scheme = "graph" )
#saveRDS(traj_milo, paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/02_MiloObject_Nhoods_Subsampling_0.25_",n,".rds"))

print("4. Counting cells per neighbourhood ---")
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), samples="assignment")


#design DF 
traj_design <- data.frame(colData(traj_milo))[,c("Age", "assignment", "date")]
traj_design <- distinct(traj_design)
#traj_design$Age_binned <- as.integer(paste0(substr(traj_design$Age, 1, 1), "0"))
traj_design$Age_binned <- ifelse(traj_design$Age<=40, 'Y', ifelse(traj_design$Age>=60, 'O', 'M'))
rownames(traj_design) <- traj_design$assignment

#Check separation 
check.sep <- checkSeparation(traj_milo, design.df=traj_design, condition='Age_binned')


mat <- traj_milo@nhoods
nhoodSize <- colSums(mat != 0) %>% as.data.frame()
colnames(nhoodSize) <- "nh_size"

saveRDS(nhoodSize, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/00_Check_K/Nnhoods_AllData_", k,".rds" ))

saveRDS(check.sep, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/00_Check_K/checkSeparation_AllData_", k,".rds" ))

 
# plot nhood size all cells 
nhSize <-  do.call(rbind.data.frame, lapply(c(60, 70, 80, 90, 100), function(k){
 df <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/00_Check_K/Nnhoods_AllData_", k,".rds" ))
  df$k <- k
  return(df)}))

nhSize %>% group_by(k) %>% count()
nhSize$k <- as.factor(nhSize$k)
nh_size_p1 <- ggplot(nhSize, aes(x=nh_size, color=k, fill=k))+geom_density(alpha=0.6)+  
  geom_vline(data=medians, aes(xintercept=median_nh_size, color=k), size=0.5)+
  theme+scale_x_continuous(limits = c(0, 1300))+scale_color_manual(values = brewer.pal("Set1", n=5))+scale_fill_manual(values = brewer.pal("Set1", n=5))


# Plot distribution of nhood size 
nh_size_p2 <- ggplot(nhSize, aes(x=k, fill=k, y=nh_size))+geom_boxplot()+
  theme+scale_color_manual(values = brewer.pal("Set1", n=5))+scale_fill_manual(values = brewer.pal("Set1", n=5))



pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS7/FigS7_Nhood_size1_all.pdf"), width =3.60, height = 2.56 )
nh_size_p1
dev.off()

pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS7/FigS7_Nhood_size2_all.pdf"), width =3.60, height = 2.56 )
nh_size_p2
dev.off()



# check separation 
checkSep <- do.call(rbind.data.frame, lapply(c(60,80, 70, 90, 100), function(k){
  df <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/00_Check_K/checkSeparation_AllData_", k,".rds" )) %>% as.data.frame()
  df$k <- k
  return(df)}))

checkSep$k <- as.factor(checkSep$k)

ch_sep_p <- ggplot(checkSep,aes(x=.) )+geom_bar(stat="count", fill=alpha(blue, 0.6))+facet_grid(~k)+theme+ylab("# Nhoods")+xlab("perfect separation of nhoods")

pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS7/FigS7_checksep_all.pdf"), width =6.45, height = 2.56 )
ch_sep_p
dev.off()




