#!/usr/bin/env Rscript


# Compute cell type proportions using Oelen 2022 dataset --- 

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
  make_option(c("--celltype_sex"), action="store", default="Projects/scRNAseq/aripol1/OneK1K_Age/scripts/MS_02_SubpopulationsByMarkers.celltype_sex_Oelen.tab", type='character',
              help="Cell type and Sex"),
  make_option(c("--markers_dir"), action="store", default="Projects/scRNAseq/msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers", type='character',
              help="Nhood markers"),
  make_option(c("--markers_candidates"), action="store", default="Projects/scRNAseq/aripol1/OneK1K_Age/scripts/MS_02_SubpopulationsByMarkers.markers_candidates_Oelen.tab", type='character',
              help="Nhood markers candidates"),
  make_option(c("--so_dir"), action="store", default="Data/scRNAseq/Yazar2022/sce_data_objects", type='character',
              help="Input directory."),
  make_option(c("--out_dir"), action="store", default="Projects/scRNAseq/aripol1/OneK1K_Age/MS_02_SubpopulationsByMarkers/", type='character',
              help="Output directory"))
opt = parse_args(OptionParser(option_list=option_list))


# Packages
shhh(library(Seurat))
shhh(library(SeuratObject))
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(tidyr))
shhh(library(ggplot2))
shhh(library(ggpubr))
shhh(library(ggrepel))
shhh(library(viridis))
shhh(library(RColorBrewer))

################################## Set Variables and load Data ################################## 
# Output directory
out.dir <- paste0(opt$out_dir, '/')
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

################################## Functions ################################## 
# Check markers (optional)
# celltype <- celltypes[3]
# ct_sex_df = celltype_sex.df
# n = 10
check_markers <- function(celltype, ct_sex_df = celltype_sex.df, n = 10){
  print(celltype)
  df_i <- ct_sex_df[ct_sex_df$celltype==celltype,]
  sex <- df_i$sex
  direction <- df_i$direction
  print(sex)
  print(direction)
  
  # Markers
  hvgs_fn <- paste0(opt$markers_dir, '/Markers_',
                    direction, '_NonSignif_', celltype, '_Sex_', sex, '.rds')
  hvgs_df <- readRDS(hvgs_fn)
  markers_df <- hvgs_df[hvgs_df$p_val_adj<0.05 & hvgs_df$avg_log2FC>1,]
  # markers_df[rownames(markers_df)=='CD27',] #testing in B naive & B intermediate
  print(paste0('nMarkers =  ', nrow(markers_df)))
  
  ## sort markers
  markers.p_val_adj.df <- markers_df[order(markers_df$p_val_adj),]
  print('Sorted by: p_val_adj... ')
  head(markers.p_val_adj.df, n = n)
  cat('\n')
  
  print('Sorted by: avg_log2FC... ')
  markers.avg_log2FC.df <- markers_df[order(-markers_df$avg_log2FC),]
  head(markers.avg_log2FC.df, n = n)
  cat('\n')
  
  print('Sorted by: p_val_adj & avg_log2FC...')
  markers.p_val_adj.avg_log2FC.df <- markers_df[order(markers_df$p_val_adj, 
                                                      -markers_df$avg_log2FC),]
  head(markers.p_val_adj.avg_log2FC.df, n = n)
  cat('\n')
  
  print('Sorted by: p_val_adj & avg_log2FC & pct expression in both groups...')
  markers.p_val_adj.avg_log2FC.pcts.df <- markers_df[order(markers_df$p_val_adj,
                                                           -markers_df$avg_log2FC,
                                                           -markers_df$pct.1,
                                                           -markers_df$pct.2),]
  head(markers.p_val_adj.avg_log2FC.pcts.df, n = n)
  cat('\n')
  
  print('Sorted by: p_val_adj & pct expression in both groups...')
  markers.p_val_adj.pcts.df <- markers_df[order(markers_df$p_val_adj, 
                                                -markers_df$pct.1,
                                                -markers_df$pct.2),]
  head(markers.p_val_adj.pcts.df, n = n)
  cat('\n')
  
  print('Sorted by: avg_log2FC & pct expression in both groups...')
  markers.avg_log2FC.pcts.df <- markers_df[order(-markers_df$avg_log2FC, 
                                                 -markers_df$pct.1,
                                                 -markers_df$pct.2),]
  head(markers.avg_log2FC.pcts.df, n = n)
  cat('\n')
  
  # Check nCells with expression for the candidate genes
  markers.p_val_adj <- rownames(markers.p_val_adj.df)[1:n]
  markers.avg_log2FC <- rownames(markers.avg_log2FC.df)[1:n]
  
  cat('\n')
}

## Check markers candidates - Check in which pct of cells these genes are expressed in F/M conditions
# celltype <- celltypes[1]
# ct_sex_df = celltype_sex.df
# markers_list = markers.list
check_markers_candidates <- function(celltype, ct_sex_df = celltype_sex.df, markers_list = markers.list){
  print(celltype)
  df_i <- ct_sex_df[ct_sex_df$celltype==celltype,]
  sex <- df_i$sex
  direction <- df_i$direction
  print(sex)
  print(direction)
  
  # Markers
  hvgs_fn <- paste0(opt$markers_dir, '/Markers_',
                    direction, '_NonSignif_', celltype, '_Sex_', sex, '.rds')
  hvgs_df <- readRDS(hvgs_fn)
  markers_df <- hvgs_df[hvgs_df$p_val_adj<0.05 & hvgs_df$avg_log2FC>1,]
  # markers_df[rownames(markers_df)=='CD27',] #testing in B naive & B intermediate
  print(paste0('nMarkers =  ', nrow(markers_df)))
  
  # Check pct
  genes <- unlist(markers_list[[celltype]])
  markers.genes <- markers_df[rownames(markers_df)%in%genes,]
  print(markers.genes)
  return(markers.genes)
  
  cat('\n')
}

## Subset
# celltype <- celltypes.all[3]
# so_list <- so.list
# markers_list <- markers.list
subset_so <- function(celltype, so_list = so.list, markers_list = markers.list){
  celltypes <- unlist(str_split(celltype, '-'))
  celltypes.label <- unname(sapply(celltypes, function(x) gsub('_', ' ', x)))
  
  print(paste0('### ' , celltype, ' ###'))
  genes_list <- markers_list[[celltype]]
  
  gene <- genes_list[[1]]
  sex <- names(so_list)[1]
  lapply(genes_list, function(gene){
    print(paste0('# ', gene))
    sapply(names(so_list), function(sex){
      print(sex)
      so_sex.all <- so_list[[sex]]
      barcodes <- rownames(so_sex.all@meta.data[so_sex.all@meta.data$cell_type%in%celltypes.label,])
      so_sex <- so_sex.all[,barcodes]
      so_counts <- so_sex[['RNA']]$counts
      so_counts_genes <- so_counts[rownames(so_counts)%in%gene,]
      print('Cells with expression...')
      print(table(so_counts_genes>0))
      so_sex.pos <- so_sex[,so_counts_genes>0]
      sce_sex.pos <- as.SingleCellExperiment(so_sex.pos)
      so_sex.neg <- so_sex[,so_counts_genes==0]
      sce_sex.neg <- as.SingleCellExperiment(so_sex.neg)
      sce_filt <- list(positive = sce_sex.pos,
                       negative = sce_sex.neg)
      sce_sex.pos.fn <- paste0(out.dir,
                               celltype, '.', sex, '.',
                               gene, '_pos.SCE_Oelen.rds')
      sce_sex.neg.fn <- paste0(out.dir,
                               celltype, '.', sex, '.',
                               gene, '_neg.SCE_Oelen.rds')
      print(paste0('Saving SCE with expression of the gene in: ', sce_sex.pos.fn))
      system.time(saveRDS(sce_sex.pos, sce_sex.pos.fn))
      
      print(paste0('Saving SCE without expression of the gene in: ', sce_sex.neg.fn))
      system.time(saveRDS(sce_sex.neg, sce_sex.neg.fn))
      cat('\n')
      return(sce_filt)
    }, simplify = FALSE)
  })
}

# Now with the new subcelltype classification, recompute the proportions 
get_proportions <- function(df, donor_var, celltype_var){
  group_vars <- c(donor_var, celltype_var)
  df %>%
    dplyr::group_by_at(vars(donor_var)) %>%
    dplyr::summarise(n = n()) -> count.df
  count.df %>%
    dplyr::group_by_at(vars(donor_var)) %>%
    dplyr::mutate(freq = n / sum(n)) %>% as.data.frame() -> prop.df
  colnames(prop.df) <- c('donor', 'cell_type', 'n', 'freq')
  prop.df <- prop.df[order(prop.df$donor, -prop.df$freq),]
  prop.df %>%
    group_by(cell_type, donor) %>%
    summarise(n = n()) -> check_prop.df
  print(all(check_prop.df$n==1))
  return(prop.df)
}


################################## Analyses ################################## 
# Read Input files
## Seurat objects
so <- readRDS(paste0("Projects/scRNAseq/aripol1/wijst-2020-hg19/v1/aging/00.so_split_by_celltype/v2/so.rds"))
mdata <- readRDS(paste0("Projects/scRNAseq/aripol1/wijst-2020-hg19/v1/aging/00.so_split_by_celltype/v2/md.rds"))
identical(colnames(so), mdata$bare_barcode_lane)
#if true
so@meta.data$sex <- mdata$Gender
so@meta.data$cell_type <- gsub("_", " ", so@meta.data$cell_type)
# split by sex 
so.list <- list( subset(so, subset = sex == "M"), subset(so, subset = sex == "F"))
names(so.list) <- c("M", "F")

# Check markers (optional)
celltype_sex.df <- read.table(opt$celltype_sex, header=T)
celltype_sex.df$celltype <- gsub('_', ' ', celltype_sex.df$celltype)
celltypes <- celltype_sex.df$celltype
check_markers.list <- sapply(celltypes, function(i) check_markers(i), simplify = FALSE)

# Subset Seurat object by candidate markers
## Define candidates
markers.df <- read.table(opt$markers_candidates, header = T)
markers_tmp.list <- split(markers.df, markers.df$celltype)
markers.list <- lapply(markers_tmp.list, function(x) unlist(str_split(x$markers, ',')))
celltypes.all <- names(markers.list)

## Check in which pct of cells these genes are expressed in F/M conditions
check_markers.list <- sapply(celltypes, function(i) check_markers_candidates(i), simplify = FALSE)

## Subset
system.time(subset_so.list <- sapply(celltypes.all, function(i) subset_so(i), simplify = FALSE))


## Recompute the proportions 
out.dir <- opt$out_dir

mono_isg <- do.call(rbind, lapply(c("M", "F"), function(sex, celltype, gene) {
  colData(readRDS(paste0(out.dir, celltype, ".", sex, ".", gene, "_pos.SCE_Oelen.rds"))) %>% as.data.frame()
}, celltype = "CD14_Mono", gene = "ISG15"))
mono_ifi <- do.call(rbind, lapply(c("M", "F"), function(sex, celltype, gene) {
  colData(readRDS(paste0(out.dir, celltype, ".", sex, ".", gene, "_pos.SCE_Oelen.rds"))) %>% as.data.frame()
}, celltype = "CD14_Mono", gene = "IFI6"))


cd8_gzmb <- do.call(rbind, lapply(c("M", "F"), function(sex, celltype, gene) {
  colData(readRDS(paste0(out.dir, celltype, ".", sex, ".", gene, "_pos.SCE_Oelen.rds"))) %>% as.data.frame()
}, celltype = "CD8_TEM", gene = "GZMB"))
cd8_gzmh <- do.call(rbind, lapply(c("M", "F"), function(sex, celltype, gene) {
  colData(readRDS(paste0(out.dir, celltype, ".", sex, ".", gene, "_pos.SCE_Oelen.rds"))) %>% as.data.frame()
}, celltype = "CD8_TEM", gene = "GZMH"))



#add the new cell types to the metadata 
mdata$cell_type <- gsub("_", " ", mdata$cell_type)
#Cd14 Mono
mdata[mdata$cell_type == "CD14 Mono",]$cell_type <- "CD14 Mono ISG15-"
mdata[mdata$bare_barcode_lane %in%  intersect(mono_isg$bare_barcode_lane,mono_ifi$bare_barcode_lane ),]$cell_type <- "CD14 Mono ISG15+"
table(mdata$cell_type)
#Cd8 TEM 
mdata[mdata$cell_type == "CD8 TEM",]$cell_type <- "CD8 TEM GZMB-"
mdata[mdata$bare_barcode_lane %in%  intersect(cd8_gzmb$bare_barcode_lane, cd8_gzmh$bare_barcode_lane),]$cell_type <- "CD8 TEM GZMB+"
table(mdata$cell_type)
mdata.fn <- paste0(out.dir, 'mdata_new_celltypes_Oelen.rds')
saveRDS(mdata, mdata.fn)


#get proportions 
prop.df <- get_proportions(mdata, 'assignment', "cell_type")
prop.df_fn <- paste0(out.dir, 'proportions_all_cts_Oelen.rds')
print(paste0('Saving cell type proportions by donor sample (all cell types): ', prop.df_fn))
saveRDS(prop.df, prop.df_fn)

opt$cell_level <- "cell_type"
# Filter: We only considered cell types with >5 cells per donor in at least 5 donors
## which cell types do we consider
filtered_data <- prop.df %>%
  filter(n > 5) %>%
  group_by_at(vars(opt$cell_level)) %>%
  summarise(n_donors = n()) %>%
  arrange(desc(n_donors))
cts_in <- unique(filtered_data[filtered_data$n_donors>=4,][[opt$cell_level]])
setdiff(unique(mdata[[opt$cell_level]]), cts_in)

## filter cell metadata
cell_metadata_full.filt <- droplevels(mdata[mdata[[opt$cell_level]]%in%cts_in,])
prop_filt.df <- get_proportions(cell_metadata_full.filt, 'assignment', 'cell_type')
prop_filt.df_fn <- paste0(out.dir, 'proportions_Oelen.rds')
print(paste0('Saving cell type proportions by donor sample (all cell types): ', prop_filt.df_fn))
saveRDS(prop_filt.df, prop_filt.df_fn)
mdata_filt.fn <- paste0(out.dir, 'mdata_Oelen.rds')
saveRDS(cell_metadata_full.filt, mdata_filt.fn)





#test 
so_counts <- so[['RNA']]$counts
so_counts_genes <- so_counts[rownames(so_counts)%in%gene,]
print('Cells with expression...')
print(table(so_counts_genes>0))
so_sex.pos <- so_sex[,so_counts_genes>0]
sce_sex.pos <- as.SingleCellExperiment(so_sex.pos)
so_sex.neg <- so_sex[,so_counts_genes==0]
sce_sex.neg <- as.SingleCellExperiment(so_sex.neg)
sce_filt <- list(positive = sce_sex.pos,
                 negative = sce_sex.neg)


