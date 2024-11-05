
# This script extract cells for significant nhoods per cell type and performs find markers ----

library(miloR); library(SingleCellExperiment); library(Seurat);library(ggplot2);library(scater);library(clusterProfiler)
Csparse_validate = "CsparseMatrix_validate"

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


# 1. Read milo object and DA results -------
print("1. Reading input data  --------")
milo <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/04_MiloObject_0.5_Preprocessed_",sex,".rds"))
da_results <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DAA_results_cell_type_annotation__Sex_0.5_",sex,".rds") )


# 2. Obtain a df of the cells per nhood with the cell type label -------
#extract cells per nhood 
sparse_mat <- milo@nhoods
colnames(sparse_mat) <- c(1:ncol(sparse_mat))
non_zero <- Matrix::summary(sparse_mat)

cell_nhood <- data.frame(
  cell_id = rownames(sparse_mat)[non_zero$i],
  Nhood = colnames(sparse_mat)[non_zero$j]
)

#extract cell_id from each cell 
md_milo <- data.frame(cell_id=milo$bare_barcode_lane, cell_type=milo$cell_type)

# add cell type label 
cell_type_nhood <- cell_nhood %>% left_join(md_milo, by="cell_id")

#add cell_type label of the nhood 
da_results_subset <- da_results[,c("Nhood", "logFC", "SpatialFDR", "cell_type")] %>% rename("cell_type_nhood"="cell_type")

#get all the dataframe we need for finding markers 
da_results_subset$Nhood <- as.character(da_results_subset$Nhood)
df_find_markers <- cell_type_nhood %>% left_join( da_results_subset, by="Nhood")
saveRDS(df_find_markers, paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Nhood-cell_id_df_Sex_", sex, ".rds"))


#split by cell type 
df_find_markers <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Nhood-cell_id_df_Sex_", sex, ".rds"))
df_find_markers$direction <- ifelse(df_find_markers$logFC < 0, "down", "up")
df_find_markers$signif <- ifelse(df_find_markers$SpatialFDR < 0.05,"ss", "ns")

# keep only cells that match the cell type with the cell type label assigned to the Nhood
df_find_markers <- df_find_markers[df_find_markers$cell_type ==df_find_markers$cell_type_nhood, ]
dim(df_find_markers)

# 3. Get a list for the cell id of each nhood and condition -------
result <- list()

# Split the dataframe by cell_type
split_by_celltype <- split(df_find_markers, df_find_markers$cell_type)

# Iterate over each cell_type
result <- lapply(names(split_by_celltype), function(ct) {
  celltype_data <- split_by_celltype[[ct]]
  
  # Split by direction
  split_by_direction <- split(celltype_data, celltype_data$direction)
  
  # Iterate over each direction 
  direction_lists <- lapply(names(split_by_direction), function(dir) {
    direction_data <- split_by_direction[[dir]]
    
    # Split by significance
    split_by_signif <- split(direction_data, direction_data$signif)
    
    # Create lists for ss and ns
    list(
      ss = split_by_signif$ss$cell_id,
      ns = split_by_signif$ns$cell_id
    )
  })
  
  # Assign direction lists to the corresponding direction in the final structure
  names(direction_lists) <- names(split_by_direction)
  
  return(direction_lists)
})

# Name the list by cell_type
names(result) <- names(split_by_celltype)

saveRDS(result,paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Cell_perNhood_list_",sex,".rds" ))


# 4.Obtain markers ---- 
#read the seurat object 
so <- readRDS(paste0(basepath, "/Data/scRNAseq/Yazar2022/sce_data_objects/AllCells_0.5_so_preprocessed_",sex,".rds"))

#read list with cell ids per Nhood 
da_list <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Cell_perNhood_list_",sex,".rds" ))

#subset to the cells from the nhoods label as particular cell type  
get_marker_genes <- function(celltype, direction){
  so_cell <- subset(so, subset =bare_barcode_lane  %in% unique(unlist(da_list[[celltype]])))
  #saveRDS(so_cell,paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/so_",celltype, "_Subsampling_",n,".rds" ) )
  print(celltype)
  # find the markers of nhoods for that cell type ---- 
  #we remove cells that are found in both nhoods ss and the rest 
  if(direction == "up"){
    cells_to_remove <- intersect(unique(da_list[[celltype]][["up"]][["ss"]]), unique(unlist(da_list[[celltype]][["down"]])))
    cells_1 <- setdiff(unique(da_list[[celltype]][["up"]][["ss"]]), cells_to_remove)
    cells_2 <- unique(setdiff(unique(unlist(da_list[[celltype]][["down"]])), cells_to_remove))
    
  }else{
    cells_to_remove <- intersect(unique(da_list[[celltype]][["down"]][["ss"]]), unique(unlist(da_list[[celltype]][["up"]])))
    cells_1 <- setdiff(unique(da_list[[celltype]][["down"]][["ss"]]), cells_to_remove)
    cells_2 <- unique(setdiff(unique(unlist(da_list[[celltype]][["up"]])), cells_to_remove))
  }
  
  # find highly variable genes to reduce the set of tested genes -----
  
  print("3.2. find HVG --")
  set.seed(101)
  so_cell <- FindVariableFeatures(so_cell, selection.method = "vst", nfeatures = 2000)
  hvgs <- head(VariableFeatures(so_cell), 2000)
  saveRDS(hvgs, paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/HVG_",direction, "_NonSignif_", celltype, "_Sex_",sex,".rds" ))
  
  markers <- FindMarkers(so_cell, ident.1 = cells_1, ident.2 = cells_2, features=hvgs, group.by= "assignment") 
  table(markers$p_val_adj < 0.05)
  saveRDS(markers,paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Markers_",direction, "_NonSignif_", celltype, "_Sex_",sex,".rds" ))
  
}
# cells that have nhoods significantly enriched 
if(sex == "M"){
  lapply(c("B naive", "B intermediate", "CD8 TEM", "CD4 TCM", "NK"), function(celltype) get_marker_genes(celltype, "up"))
  
}else{
  lapply(c("B naive", "B intermediate", "CD8 TEM", "CD4 TCM", "CD16 Mono", "CD14 Mono", "NK"), function(celltype) get_marker_genes(celltype, "up"))
  
}

# cells that have nhoods significantly depleted 
lapply(c("B naive", "B intermediate", "B memory", "CD8 TEM","CD8 Naive", "CD8 TCM", "MAIT", "CD4 Naive", "CD4 TCM","CD4 TEM"), function(celltype) get_marker_genes(celltype, "down"))



# 5.Perform enrichments of marker genes ---- 
enrichments_markers<-function(dir, celltype){ 
  print(celltype)
  if(dir=="up"){genes <-mk_list_up[[celltype]]}else{genes <-mk_df_ss_down[[celltype]]}
  universe <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/HVG_", dir,"_NonSignif_", celltype, "_Sex_",sex,".rds" ))
  go <- enrichGO(genes, OrgDb = "org.Hs.eg.db", universe = universe,ont = "BP", keyType = "SYMBOL")
  print(dotplot(go))
  go_df <- go@result
  return(go_df)
}


#enrichment up ---

if(sex == "M"){
  mk <-  lapply(c( "CD8 TEM", "B naive"), function(celltype) enrichments_markers( "up", celltype)) #only cell types that have enrichment
  names(mk) <-c( "CD8 TEM", "B naive")
}else{
  mk <- lapply(c("CD14 Mono", "B naive", "CD8 TEM"), function(celltype)enrichments_markers( "up", celltype))
names(mk) <-c("CD14 Mono", "B naive", "CD8 TEM")}

list_go_clean <- purrr::map(mk, ~ purrr::compact(.)) %>% purrr::keep(~length(.) != 0)
for(c in names(list_go_clean)){list_go_clean[[c]]$NhoodGroup <- c}
go_df <- do.call(rbind.data.frame, list_go_clean)
saveRDS(go_df,paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Markers_enrichments_up_NonSignif_Sex_",sex,".rds"))



# 6. reduction of enrichments for plotting ---- 
library(rrvgo);library(org.Hs.eg.db)

go_df_F <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Markers_enrichments_up_NonSignif_Sex_F.rds"))
go_df_F$sex <- "F"
go_df_M <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Markers_enrichments_up_NonSignif_Sex_M.rds"))
go_df_M$sex <- "M"
go_df <- rbind(go_df_F, go_df_M)
go_df_signif <- go_df[go_df$p.adjust < 0.05, ]
simMatrix <- calculateSimMatrix(go_df_signif$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
)


scores <- setNames(-log10(go_df_signif$qvalue), go_df_signif$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.9,
                                orgdb="org.Hs.eg.db")
saveRDS(reducedTerms, paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/01_Milo_NhoodMarkers/Markers_enrichments_up_NonSignif_Sex_",sex,".rds"))

