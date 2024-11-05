
# Plots results Differential Gene Expression Analysis
library(ggplot2); library(dplyr); library(RColorBrewer); library(UpSetR); library("AnnotationDbi"); library('org.Hs.eg.db');library(ComplexHeatmap);library(dplyr)
library(readxl);library(Seurat);library(plyr);library(ggsignif);library(clusterProfiler);library(tidyr);library(ggpubr);library(ggrepel);library(purrr);library(formattable);library(ggside);
library(scales);library(patchwork);library(ggpubr);library(gridExtra);library(dendextend); library(openxlsx); library(patchwork)


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


#data
metadata <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/metadata_processed.rds"))
cells_to_keep <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/cells_to_keep.rds"))
order_cells<- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/new_order_cells.rds"))
celltype_l1 <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/celltypes_equivalence.rds"))
celltype_l1$cell_type <- factor(celltype_l1$cell_type, levels = order_cells$cell_type)


#themes and palettes 
pal_direction <- c("up" = "#bc4749" , "down"= "#264653")
computer <- "work"
source(paste0(data_path, "/msopena/01_meta-analysis_Age_Sex/scripts/themes.R"))
pheno <- "Age"
cell_level <- "cell_type"
df <- get_significant_degs(pheno, cell_level)

cell_types <- unique(df$celltype)
cell_types <- cell_types[!cell_types %in%c("Platelet")]
df <- df[df$celltype %in% cell_types,]
ncell <- metadata %>% dplyr::group_by(cell_type) %>% dplyr::count()
colnames(ncell) <- gsub("cell_type", "celltype", colnames(ncell))
colnames(ncell) <- gsub("n", "ncells", colnames(ncell))

#table of DEGs

all_genes <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/DEG_Age_cell_type_02.rds"))


# 1. Number of DEGs (Fig. 1A) -------------------------------------------------------------
# Extract the number of DEGs
degs <- get_significant_degs("Age", "cell_type")
degs <- reorder_cells(degs, reverse = T, neworder = T)
degs <- degs[degs$celltype != "Platelet",]

# Extract the numner of cells 
ncells <- metadata %>% dplyr::group_by(cell_type) %>% dplyr::count() %>% drop_na()
colnames(ncells)[1] <- "celltype"
ncells <- ncells[ncells$celltype %in% degs$celltype,] %>% as.data.frame()
ncells <- reorder_cells(ncells, reverse = T, neworder = T)
df <- df %>% left_join(ncells, by="celltype") 

# Compute the bias towards up/down regulation 
result_df <- data.frame(celltype = character(0), binom.test.p = numeric(0), prop_up = numeric(0))
for (cell in unique(degs$celltype)){ 
  print(cell)
  binom <- binom.test(length(degs[degs$direction == "down" & degs$celltype == cell,]$gene) , length(degs[degs$celltype == cell,]$gene))$p.value
  upregulated <- length(degs[degs$direction == "up" & degs$celltype == cell,]$gene)
  downregulated <- length(degs[degs$direction == "down" & degs$celltype == cell,]$gene)
  
  success_ratio = upregulated / (downregulated+upregulated)
  result_df <- rbind(result_df, data.frame(celltype = cell, binom.test.p = binom, prop_up = success_ratio))
}
result_df$p.val_dir <- ifelse(result_df$prop_up > 0.5, result_df$binom.test.p , -result_df$binom.test.p )
result_df$direction <- ifelse(result_df$prop_up > 0.5, "upregulation_bias" , "downregulation_bias" )

#saveRDS(result_df, paste0( data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/biomial_test_bias_02.rds"))


# Count the number of up/down DEGs per cell type 
df <- degs %>% dplyr::group_by(celltype, direction) %>%dplyr::summarise(n = n()) %>% dplyr::mutate(freq = n / sum(n)*100)
df$p.val <- unlist(lapply(df$celltype, function(cell) result_df[result_df$celltype == cell, ]$binom.test.p))
df$signif <- ifelse(df$p.val< 0.05, "*", "")
for (cell in unique(df$celltype)){
  if(length(df[df$celltype == cell,]$celltype) == 2  && df[df$celltype == cell & df$direction == "down",]$signif == "*"){
    df_cell <- df[df$celltype == cell,]
    if(df_cell[df_cell$direction == "up",]$n < df_cell[df_cell$direction == "down",]$n){
      df[df$celltype == cell & df$direction == "up", ]$signif <- ""
    } else{
      df[df$celltype == cell & df$direction == "down", ]$signif <- ""}}}

df <- reorder_cells(df, neworder = T)
df$celltype <- as.character(df$celltype)

# split the data 
up <- df[df$direction == "up", ]
up$signif_n <- paste0(up$n, " ", up$signif)
up <- reorder_cells(up, reverse = T, neworder = T)
down <- df[df$direction == "down", ]
down$signif_n <- paste0(down$signif, " ", down$n)
down <- reorder_cells(down, reverse = T, neworder = T)

# Plot 1- Number of cells
p_ncells <- ggplot(ncells, aes(y=celltype, x=n)) + geom_bar(stat="identity", fill=alpha(blue, 0.4))+  scale_x_continuous(labels = function(x) format(x / 1e6, scientific = FALSE)) +
  theme  + xlab("N cells (million)") + ylab("") +
  theme(axis.text.x=element_text(size=9), axis.text.y=element_text(size=11), axis.title=element_text(size=11),  strip.text=element_blank())+ facet_grid(celltype_l1~ ., scales="free", space="free")

# Plot 2- Number of DEGs
p_ndegs <- ggplot(degs, aes(y=celltype)) + geom_bar(stat="count", fill="#264653")+ 
  theme  + xlab("Number of DEGs") + ylab("")+ geom_text(stat='count', aes(label=..count..), size=3.5, position = position_stack(vjust = 1), vjust=0.5, hjust=-0.1)+
  theme(axis.text=element_text(size=12), strip.text=element_blank(),axis.title=element_text(size=11), axis.ticks.y=element_blank(), axis.text.y=element_blank())+
  facet_grid(celltype_l1~., scales = "free", space="free")+ scale_x_continuous(limits = c(0, 275000), breaks = c( 0,  1000, 2000)) 

# Plot 3 - Directionality DEGs
p_direction <- ggplot(df, aes(y=celltype, x=freq)) +
  # Upregulated genes
  geom_col(data = up, aes(y = celltype, x = freq, fill=direction), alpha=0.9, width=0.75) +
  # Downregulated genes
  geom_col(data = down, aes(y = celltype, x= -freq, fill=direction), width=0.75) +
  geom_vline(xintercept=0, linetype = "dashed", linewidth=0.3) +
  theme  + geom_text(data = up, aes(label=signif), size=5, position = position_stack(vjust = 1), vjust=0.7, hjust=-0.1) +  geom_text(data = down, aes(label=signif), position = position_stack(vjust = - 1), vjust=0.7, size=5, hjust=1) + 
  scale_fill_manual(values = c("up"=red, "down"=alpha(blue, 0.7)))+
  xlab("Percentage of DEGs") + ylab("") +scale_x_continuous(labels = abs)+ scale_y_discrete(limits = rev(levels(df$celltype)))+
  theme(axis.text=element_text(size=11), axis.text.y=element_blank(), axis.ticks.y=element_blank(),axis.title=element_text(size=10),strip.text=element_text(size=12), strip.background=element_rect(fill = "white", colour = "white"), strip.switch.pad.grid=unit(1, "pt"))+
  facet_grid(celltype_l1~ ., scales="free", space="free")+ scale_x_continuous(labels = abs, limits = c(-120, 120), breaks = c(-100, 0,  100)) 


p_ndegs_direction <- ggplot(df, aes(y=celltype, x=freq)) +
  # Upregulated genes
  geom_col(data = up, aes(y = celltype, x = n, fill=direction), alpha=0.9, width=0.75) +
  # Downregulated genes
  geom_col(data = down, aes(y = celltype, x= -n, fill=direction), width=0.75) +
  geom_vline(xintercept=0, linetype = "dashed", linewidth=0.3) +
  theme  + geom_text(data = up, aes(label=signif_n), hjust = -0.2, size = 4, position = position_dodge(width = 1)) +  geom_text(data = down, aes(label=signif_n),  vjust = 0.5, hjust = 1.5, size = 4, position = position_dodge(width = 1)) + 
  scale_fill_manual(values = c("up"=red, "down"=alpha(blue, 0.7)))+
  xlab("Number of DEGs") + ylab("") +scale_x_continuous(labels = abs)+ scale_y_discrete(limits = rev(levels(df$celltype)))+
  theme(axis.text=element_text(size=11),axis.title=element_text(size=12),strip.text=element_text(size=11), legend.position="top", strip.background=element_rect(fill = "white", colour = "white"), 
        strip.switch.pad.grid=unit(1, "pt"), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  facet_grid(celltype_l1~ ., scales="free", space="free")+ scale_x_continuous(labels = abs, limits = c(-2800, 2800), breaks = c(-2000,-1000, 0,  1000,2000)) 


# Figure 1A 
p_ncells+plot_spacer()+p_ndegs+plot_spacer()+p_direction+plot_layout(widths = c(8, -3.5, 22, -3.5, 17), guides = "collect")

Fig_1A <- p_ncells+plot_spacer()+p_ndegs_direction+plot_layout(widths = c(8, -3.5 , 22))


# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig1/Fig1A_numDEGs_02.pdf"), height =6, width= 5 )
# Fig_1A
# dev.off()



# 2. Sharing (Fig. 1B) ------------------------------------------------------------------
#  Cell type sharing
df <- degs  %>% dplyr::filter(!celltype %in% c( "cDC2", "NK_CD56bright", "gdT", "NK Proliferating", "Platelet")) %>% dplyr::group_by(gene) %>% dplyr::count() %>% as.data.frame()
df$x <- "DEGs"
df$label <- NA
#df[df$n > 10,]$label <- df[df$n > 10,"gene"]
df$concordant <- NA
non_concordant <- df[df$gene %in% degs[degs$direction == "up",]$gene & df$gene %in% degs[degs$direction == "down",]$gene ,]$gene
df[df$n > 1, ]$concordant <- ifelse(  df[df$n > 1, ]$gene %in% non_concordant, 0, 1)
df$count<- df$n
df$group <- ifelse(df$concordant == 1, "shared_concordant", ifelse(df$concordant == 0, "shared_non_concordant", "specific") )
df[is.na(df$group),]$group <- "specific"
df_perc <- df %>%   dplyr::group_by(group) %>% dplyr::summarise(cnt = n()) %>% mutate(freq = round(cnt / sum(cnt), 3))
df_perc$group <- factor(df_perc$group, levels= rev(c("shared_concordant", "shared_non_concordant",  "specific")))
df_perc$x <-"DEGs"
df_shared <-df[!is.na(df$concordant),] 
df_shared$concordant <- factor(df_shared$concordant, levels=c(1, 0))
df_shared$category <- NA
df_shared[df_shared$concordant == 1, ]$category <- "concordant"
df_shared[df_shared$concordant == 0, ]$category <- "non concordant"
#saveRDS(df_shared,paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/sharing_02.rds") )


p_sharing <- ggplot(df_perc, aes(x=x, y=freq, fill=group)) + geom_bar(stat="Identity")+theme+ xlab("")+ylab("") + scale_alpha_manual(values=rev(c(0.3, 0.7, 1)))+ylab("Frequency of DEGs")+
  geom_text(aes(label=cnt), position = position_stack(vjust=0.5), size=4, show.legend = F)+scale_alpha_manual(values = c(0.9, 0.6))+coord_flip()+scale_fill_manual(values=c("shared_non_concordant"= alpha(red, 0.9), "shared_concordant"= alpha(blue, 0.9), "specific"=alpha(blue, 0.5)))+
  theme(legend.position="top",axis.text.y=element_blank(), axis.text.x=element_text(size=11),axis.title=element_text(size=12), legend.title=element_blank(),axis.ticks.y=element_blank(), aspect.ratio=0.2, legend.key.size = unit(0.3, 'cm')) +
  guides(alpha = guide_legend(reverse=T))

p_concordance <- ggplot(df_shared, aes(x=as.factor(n), fill=category))+geom_bar(stat="count", position="dodge")+theme+coord_flip()+xlab("Number of cell types") +ylab("Number of DEGs")+
  scale_fill_manual(values=c("non concordant"= alpha("#bc4749", 0.8) , "concordant"= alpha("#264653", 0.8)))+
  theme(axis.title=element_text(size=13), legend.position="none", legend.key.size = unit(0.4, 'cm'), 
        legend.title=element_blank(), legend.text=element_text(size=12), aspect.ratio=1)


# Figure 1B 
Fig1B <- p_sharing / p_concordance 

# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig1/Fig1B_sharing_02.pdf"), height =6, width= 4 )
# Fig1B
# dev.off()



# 3. Example of highly shared gene (Fig 1C) ----------------------------------
save_plots_multiple <- function(cell_level, cells,gene){
  expr_list <- lapply(cells, function(cell_type) matrix_counts_gene(cell_level, cell_type, gene ))
  expr_gene <- do.call(rbind.data.frame, expr_list)
  degs <- get_significant_degs("Age", "cell_type")
  degs <- degs[degs$gene == gene, ]
  expr_gene$direction <- unlist(lapply(expr_gene$celltype, function(celltype) degs[degs$celltype==celltype, "direction"]))
  expr_gene <- reorder_cells(expr_gene, neworder= T)
  p <-  plot_pseudobulk_age(cell_level, NULL,gene,onecell = F, expr_gene = expr_gene )[[1]] + facet_wrap(~celltype, ncol = 3, scales="free_y")+
    theme(aspect.ratio=0.8, axis.text=element_text(size=11), strip.text=element_text(size=11))
  #Age cat all
  
  p2 <-  plot_pseudobulk_age(cell_level, NULL,gene,onecell = F, expr_gene = expr_gene )[[3]] + facet_wrap(~celltype, ncol = 4, scales="free_y")+
    theme(aspect.ratio=0.8, axis.text=element_text(size=11), strip.text=element_text(size=11))
  
  pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig1/Fig1B2M_sharing.pdf"),  height =5, width= 5 )
  print(p)
  dev.off()
  pdf(paste0(path_plots, "/", gene, "_Multiple_cells_Age_continous.pdf"),  height =7, width= 7 )
  print(p2)
  dev.off()
}


# Fig 1C 

cells <- unique(degs[degs$gene == "B2M",]$celltype)
cells <- gsub(" ", "_", cells)
save_plots_multiple("cell_type", cells,"B2M" )




# 4. Heatmap concordance (Fig. 1D) --------------------------------------
all_genes <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/DEG_Age_cell_type_02.rds"))
all_genes <- all_genes[!all_genes$celltype %in% c("Platelet", "Plasmablast", "NK_CD56bright", "NK Proliferating", "HSPC", "gdT", "CD4 CTL", "Eryth"),]
table(all_genes$celltype, all_genes$adj.P.Val<0.05)
tested_genes <- split(all_genes$gene, all_genes$celltype)
common_genes <- Reduce(intersect, tested_genes)
length(common_genes)

df <- degs[degs$gene %in% common_genes, ]
df <- df[!df$celltype %in% c("gdT", "Platelet", "NK_CD56bright", "cDC2") ,]
gene_list <- split(df$gene, df$celltype)
df_intersect <- crossprod(table(stack(gene_list)))

# Check % of concordant genes 
calculate_similarity_degs <- function(celltype1, celltype2){
  print(paste0( "celltype1: ", celltype1, " celltype2: ", celltype2))
  set1 <- df[df$celltype == celltype1, ]
  set2 <- df[df$celltype == celltype2, ]
  common_degs <- intersect(set1$gene, set2$gene)
  if(length(common_degs )!= 0){
    directionality_similarity <- sapply(common_degs, function(gene) {set1[set1$gene == gene, ]$direction ==set2[set2$gene == gene, ]$direction})
    jaccard_similarity <- sum(directionality_similarity) / length(common_degs) * 100
  }else{
    jaccard_similarity <- NA
  }
  
  return(jaccard_similarity)
}

cell_types <- unique(df$celltype)
similarity_matrix <- matrix(0, nrow = length(cell_types), ncol = length(cell_types))
rownames(similarity_matrix) <- cell_types
colnames(similarity_matrix) <- cell_types


for (i in 1:length(cell_types)) {
  for (j in 1:length(cell_types)) {
    similarity_matrix[i, j] <- calculate_similarity_degs(cell_types[i], cell_types[j])
  }
}

row_annotation <- celltype_l1[celltype_l1$cell_type %in% colnames(similarity_matrix),]
row_annotation <- row_annotation[match(colnames(similarity_matrix), row_annotation$cell_type),]
rownames(row_annotation) <-row_annotation$cell_type 
row_annotation <- row_annotation %>% dplyr::select(predicted.celltype.l1) %>% as.data.frame()
colnames(row_annotation) <- "cell_type_l1"  
row_annotation$immune_response <- ifelse(row_annotation$cell_type_l1 %in% c("NK", "Mono"), "innate", "adaptive")

annot_colors <-list(cell_type_l1= palette_majorcells,clusters=c("upregulation_bias"=red, "downregulation_bias"=blue),immune_response=c("innate"="#bea06bff", "adaptive"="#563a08ff") )
library(dendsort)

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...), type = "average"))
sort_hclust_cols <- function(...) as.hclust(dendsort(as.dendrogram(...), type = "complete"))

heatmap <- pheatmap::pheatmap(similarity_matrix,cluster_rows =sort_hclust(hclust(dist(similarity_matrix), method = "complete")), 
                              cluster_cols  =sort_hclust(hclust(dist(similarity_matrix), method = "single")))
rw <- heatmap$tree_col
clusters <- cutree(rw, 2) %>% as.data.frame() %>% dplyr::rename("clusters"=".")
clusters$clusters <- gsub(2, "upregulation_bias", clusters$clusters)
clusters$clusters <- gsub(1, "downregulation_bias", clusters$clusters)
clusters$clusters <- as.factor(clusters$clusters)

row_annotation <- merge(clusters, row_annotation, by="row.names")
rownames(row_annotation) <- row_annotation[,1]
row_annotation <- row_annotation[,-1]
mat_breaks <- seq(0, 100, by = 10)



heatmap_concordance <- pheatmap(similarity_matrix,display_numbers = round(similarity_matrix,digits = 1),number_color = "black",annotation_row = row_annotation,
                                show_colnames = T, annotation_colors = annot_colors,breaks = mat_breaks, annotation_names_col =F, 
                                border_color = "white",treeheight_row=F, treeheight_col=20, annotation_names_row = F,name = "% concordant DEGs",clustering_method = "median",
                                fontsize_number  =8,  fontsize_row = 12, fontsize_col = 12,color = colorRampPalette(rev(brewer.pal(n =8, name ="RdYlBu")))(10))


# Figure 1D
heatmap_concordance

# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig1/Fig1D_heatmap_concordance_02.pdf"), height =6, width= 7.5 )
# heatmap_concordance
# dev.off()



# 5. Heatmap sharing (Fig. 2A, 2B) ------------------------------------

# Read in sharing information --
sharing <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/sharing_02.rds"))
non_concordant <- sharing %>% dplyr::filter(category == "non concordant")

# Read in DEG results ---
all_genes <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/DEG_Age_cell_type_02.rds"))
all_genes$celltype <- factor(all_genes$celltype, levels=rev(order_cells$cell_type))

# Read in bias estimate -
bias <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/biomial_test_bias_02.rds"))
rownames(bias) <- bias$celltype

# Get siginificant DEGs
degs <-get_significant_degs("Age","cell_type")
cells <- unique(degs[!degs$celltype %in% c("NK Proliferating", "gdT", "NK_CD56bright", "Platelet", "HSPC", "cDC2"),]$celltype)
# Select cell types to plot 
df <- all_genes[all_genes$celltype %in% cells,] %>% as.data.frame()
df <- reorder_cells(df)

# Get the common tested genes across celltypes 
list <-  split(df$gene, df$celltype)
common_genes <- Reduce(intersect, list)

# Extract a matix of logFCs across celltypes 
df <- df[df$gene %in% common_genes,]
df$significance <- ifelse(df$adj.P.Val < 0.05, "·", " ")
logfc <- data.frame("gene"=common_genes)
for (celltype in levels(df$celltype)){
  cell <- df[df$celltype == celltype, c("gene", "logFC")]
  colnames(cell) <- gsub("logFC", paste0(celltype), colnames(cell))
  print(head(cell))
  print(dim(cell))
  logfc <- merge(logfc, cell, by="gene")
}
rownames(logfc) <- logfc$gene
logfc <- as.matrix(logfc[,-1])


# Get the matrix of siginificance values across celltypes 
significance <- data.frame("gene"=common_genes)
for (celltype in unique(df$celltype)){
  cell <- df[df$celltype == celltype, c("gene", "significance")]
  colnames(cell) <- gsub("significance", paste0(celltype), colnames(cell))
  print(head(cell))
  print(dim(cell))
  significance <- merge(significance, cell, by="gene")
}
rownames(significance) <- significance$gene

library(dendsort)
heatmap_logFC<- function(method,n, nonconcordant=T){
  # Extract non concordant genes DE in more than 6 cell types 
  logfc_na_shar <- logfc[rownames(logfc) %in% sharing[sharing$n > n,]$gene,]
  if (nonconcordant==T) {logfc_na_shar <- logfc_na_shar[rownames(logfc_na_shar) %in% non_concordant$gene,]}
  
  #Same for pvalues 
  pvals <- as.matrix(significance[,-1])
  pvals <- pvals[rownames(pvals) %in% sharing[sharing$n > n,]$gene,]
  if (nonconcordant==T) {pvals <- pvals[rownames(pvals) %in% non_concordant$gene,]}
  
  # Get column annotation 
  col_annotation <- celltype_l1[celltype_l1$cell_type %in% colnames(logfc_na_shar),]
  col_annotation <- col_annotation[match(colnames(logfc_na_shar), col_annotation$cell_type),]
  rownames(col_annotation) <-col_annotation$cell_type 
  col_annotation <- col_annotation %>% dplyr::select(predicted.celltype.l1) %>% as.data.frame()
  colnames(col_annotation) <- "cell_type_l1"  
  #annot_colors <-list(cell_type_l1= c("CD4 T"= "#00563f",  "CD8 T"= "#679267","other T"= "#c0d3ac", "NK"= "#E76F51", "B"= "#8580C2", "Mono" ="#d88d8b" ),clusters=c("1"=alpha("#264653", 0.8), "2"="#bc4749"))
  
  set.seed(3)
  row_clusters <- hclust(dist(logfc_na_shar), method = "average")
  col_clusters <- hclust(dist(t(logfc_na_shar)), method ="average")

  callback = function(hc, mat){
    sv = svd(t(mat))$v[,1]
    dend = reorder(as.dendrogram(hc), wts = sv, method = "average")
    as.hclust(dend)
  }

  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...), type = "average"))
  sort_hclust_cols <- function(...) as.hclust(dendsort(as.dendrogram(...), type = "ward.D2"))

   heatmap <- pheatmap::pheatmap(t(logfc_na_shar),  clustering_method = "complete")
  # # get gene clusters and change names 
  rw <- heatmap$tree_col
  clusters_row <- cutree(rw, 2) %>% as.data.frame() %>% dplyr::rename("cluster_genes"=".")
  clusters_row$cluster_genes <- as.factor(clusters_row$cluster_genes)
  clusters_row$cluster_genes <- gsub("1", "Ribosomal_function",clusters_row$cluster_genes)
  clusters_row$cluster_genes <- gsub("2", "Immunological_process",clusters_row$cluster_genes)

  # get cell bias annotation
  clusters_col <-bias[c("direction", "celltype")] %>% as.data.frame() %>% filter(celltype %in% colnames(logfc_na_shar))
  clusters_col <- clusters_col[match(colnames(logfc_na_shar), rownames(clusters_col)),] %>% as.data.frame()
  clusters_col <- clusters_col[,c("direction")] %>% as.data.frame() 
  colnames(clusters_col) <- "bias"
  rownames(clusters_col) <- colnames(logfc_na_shar)
  clusters_col$bias <- as.factor(clusters_col$bias)

  #colors 
  annot_colors <-list(cluster_genes= c("Ribosomal_function"="lightgrey", "Immunological_process"="#777777"),bias=c("downregulation_bias"=alpha("#264653", 0.8), "upregulation_bias"="#bc4749"))
  # annot_colors <-list(Genes_clust= c("Ribosomal_function"="lightgrey", "Immunological_process"="#777777"),Binomial_p.val=annotation_colors_for_row)
  
  # Heatmap
  mat_breaks <- seq(-0.01, 0.01, by = 0.001)
  heatmap_horizontal <- pheatmap::pheatmap(t(logfc_na_shar), display_numbers = t(pvals),
                                           method = "complete",
                                          cluster_rows = sort_hclust(hclust(dist(t(logfc_na_shar)), method = "ward.D")),
                                           annotation_col = clusters_row,
                                           annotation_row = clusters_col,
                                           annotation_names_row = FALSE,
                                           annotation_names_col = FALSE,
                                           treeheight_col = 10,treeheight_row = 0,
                                           fontsize_number = 12, fontsize_col = 3, fontsize_row = 10,
                                           breaks = mat_breaks, border_color = NA,
                                           color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(20),
                                           annotation_colors = annot_colors
  )
  return(list(heatmap_horizontal))
}


hmp <- heatmap_logFC("average", 7,nonconcordant = T)


# pdf(paste0(data_path, "/msopena/02_OneK1K_Age/figures/Fig2/Fig2A_Heatmap_highly_shared_genes.pdf"), width = 8.84, height=4.70 )
# hmp
# dev.off()


# Figure 2B- Enrichments clusters 
rw <- hmp[[1]]$tree_col
clusters_row <- cutree(rw, 2)
go_1 <- clusterProfiler::enrichGO(gene = names(clusters_row[clusters_row==2]), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL")@result %>% dplyr::filter( p.adjust < 0.05)
ggplot(go_1, aes(y=Description, x=as.factor(Count), alpha=p.adjust))+geom_point(size = 4)+theme+ylab("")+theme( axis.text.y=element_text(size=12), aspect.ratio=1)+xlab("# DEGs")+scale_fill_continuous(type = "gradient")
go_1$cluster <- "Cluster 2"


go_2 <- clusterProfiler::enrichGO(gene = names(clusters_row[clusters_row==1]), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", universe = all_genes$gene)@result %>% dplyr::filter( p.adjust < 0.05)
go_2_ordered <- go_2[order(go_2$Count, decreasing = TRUE), ]


# Create the bar plot
ggplot(go_2_ordered[1:5,], aes(y = Description, x = Count, fill=p.adjust)) +
  geom_bar(stat = "identity") +
  ylab("") +theme+scale_fill_gradient(high = "lightgrey", low = "#818589")+
  xlab("# DEGs") +
  theme(axis.text.y = element_text(size = 12), aspect.ratio=1)
go_2_ordered$cluster <- "Cluster 1"

go_all <- rbind(go_1[1:3,], go_2_ordered[1:5,])

Fig2B <- ggplot(go_all, aes(y = reorder(Description, Count), x = Count, color=cluster)) +
  geom_point(aes(size = p.adjust)) +scale_size_continuous(range = c(7, 4))+
  ylab("") +theme+scale_color_manual(values= c("Cluster 1"="lightgrey", "Cluster 2"="#777777"))+
  xlab("Number of DEGs") +facet_wrap(~cluster,scales="free_y", nrow=2)+
  theme(axis.text.y = element_text(size = 13), aspect.ratio=1, strip.text=element_text(size=14), axis.title=element_text(size=14), axis.text.x=element_text(size=12))

pdf(paste0(data_path, "/msopena/02_OneK1K_Age/figures/Fig2/Fig2B_EnrichmentsClustGenes_02.pdf"), width = 8.84, height=5.43 )
Fig2B
dev.off()



# 6.Cell-type specific DEGs (Fig 2C, D)--------------------------------------
degs <- get_significant_degs("Age", "cell_type")
df_all <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/DEG_Age_cell_type_02.rds"))
df_all <- df_all[df_all$celltype != "Platelet",]
universe_list <- split(df_all$gene, df_all$celltype)
#get common tested genes 
common_genes <- Reduce(intersect,universe_list)
length(common_genes)
#degs_common <- unique(degs[degs$gene %in% common_genes, ]$gene)

specific <-  degs %>% dplyr::group_by(gene) %>% dplyr::count() %>% dplyr::filter(n ==1)
specific_genes <-unique(specific$gene)
degs$specific <- ifelse(degs$gene %in% specific_genes, "specific", "shared")
table(degs$direction, degs$specific)


#Plot number of cell-type specific DEGs ---
specific_degs <- degs[degs$specific == "specific", ]
specific_degs <- reorder_cells(specific_degs, neworder = T, reverse = T)
specific_degs$direction <- factor(specific_degs$direction, levels=c("up", "down"))
Fig2C<- ggplot(specific_degs, aes(x=celltype)) +geom_bar(stat="count", alpha=0.9, aes(fill=direction))+theme+coord_flip()+facet_grid(celltype_l1~direction, scales="free_y", space="free")+
  geom_text(stat='count', aes(label=..count..), position=position_stack(vjust=0.9), color="black")+ylab("Number of cell-type specific DEGs")+xlab("")+
  scale_fill_manual(values=c("up"= red, "down"=alpha(colour = blue, alpha = 0.3)))+theme(legend.position="none")

# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig2/Fig2C_cell_type_spec_02.pdf"), height =5.76, width= 5.28 )
# Fig2C
# dev.off()


#6. Enrichments of cell-type specific DEGs (Fig 2D) --------

GO_enrichment_cell_spec <- function(celltype, ont="ALL"){
  print(paste0("split by cell type: ", celltype))
  signif_genes <- degs[degs$celltype == celltype ,]
  print(paste0("only specific genes"))
  signif_genes <- signif_genes[signif_genes$gene %in% specific_genes,]
  signif_genes_up <- signif_genes[signif_genes$direction == "up",]$gene
  signif_genes_down <- signif_genes[signif_genes$direction == "down",]$gene
  enrich <- enrichGO(signif_genes_up, universe = universe_list[[celltype]], OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont= ont)
  saveRDS(enrich, paste0(data_path, "/msopena//02_OneK1K_Age/robjects/02_Enrichments/cell_spec_erichments_Age_cell_type_", celltype, "_up_GO_specific_enrichments_0.2.rds"))
  enrich <- enrichGO(signif_genes_down, universe = universe_list[[celltype]], OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont= ont)
  saveRDS(enrich, paste0(data_path, "/msopena//02_OneK1K_Age/robjects/02_Enrichments/cell_spec_erichments_Age_cell_type_", celltype, "_down_GO_specific_enrichments_0.2.rds"))
  
}
lapply(unique(degs$celltype), function(celltype) GO_enrichment_cell_spec(celltype))


up <- lapply(unique(degs$celltype), function(celltype) {
  file_path <- paste0(data_path, "/msopena/02_OneK1K_Age/robjects/02_Enrichments/cell_spec_erichments_Age_cell_type_", celltype, "_up_GO_specific_enrichments.rds")
  tryCatch({
    file<- readRDS(file_path)@result
  }, error = function(e) {
    NULL
  })
  return(file)
})
names(up) <- unique(degs$celltype)
up <- up[ sapply(up, is.data.frame) ]
up<- up[sapply(up, function(x) dim(x)[1]) > 0]
for(n in names(up)){up[[n]]$celltype <- n}

up_df <- do.call(rbind.data.frame, up)
table(up_df$celltype)
up_df <- up_df[up_df$pvalue < 0.05,]
table(up_df$celltype)


down <- lapply(unique(degs$celltype), function(celltype) {
  file_path <- paste0(data_path, "/msopena/02_OneK1K_Age/robjects/02_Enrichments/cell_spec_erichments_Age_cell_type_", celltype, "_down_GO_specific_enrichments.rds")
  tryCatch({
    file<- readRDS(file_path)@result
    file$celltype <- celltype 
  }, error = function(e) {
    NULL
  })
  return(file)
})
names(down) <- unique(degs$celltype)

down <- down[ sapply(down, is.data.frame) ]
down<- down[sapply(down, function(x) dim(x)[1]) > 0]

for(n in names(down)){down[[n]]$celltype <- n}
down_df <- do.call(rbind.data.frame, down)
table(down_df$celltype)

# Reduce terms 
down_df$direction <- "down"
up_df$direction <- "up"
go_terms <- rbind(down_df, up_df)
saveRDS(go_terms, paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/02_Enrichments/cell_spec_erichments_Age_02.rds"))

#reduce terms for each GO family 

reduce_terms <- function(go_terms, ont){
  simMatrix <- calculateSimMatrix(go_terms$ID,
                                  orgdb="org.Hs.eg.db",
                                  ont=ont,
  )
  scores <- setNames(-log10(go_terms$qvalue), go_terms$ID)
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Hs.eg.db")
  saveRDS(reducedTerms, paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/02_Enrichments/cell_spec_erichments_Age_cell_type_",ont,"_reducedTerms_specific_02.rds"))
  return(reducedTerms)
}

# reducedTerms_BP <- reduce_terms(go_terms[go_terms$ONTOLOGY == "BP",], "BP")
# reducedTerms_MF <- reduce_terms(go_terms[go_terms$ONTOLOGY == "MF",], "MF")
# reducedTerms_CC <- reduce_terms(go_terms[go_terms$ONTOLOGY == "CC",], "CC")

reducedTerms_BP <-readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/02_Enrichments/cell_spec_erichments_Age_cell_type_BP_reducedTerms_specific_02.rds"))
reducedTerms_MF <-readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/02_Enrichments/cell_spec_erichments_Age_cell_type_MF_reducedTerms_specific_02.rds"))
reducedTerms_CC <-readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/02_Enrichments/cell_spec_erichments_Age_cell_type_CC_reducedTerms_specific_02.rds"))
reducedTerms <- rbind(reducedTerms_BP, reducedTerms_MF, reducedTerms_CC)

go_terms <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/02_Enrichments/cell_spec_erichments_Age_02.rds"))

plot_reduced_terms <- function(reducedTerms, go_df){
  go_df_signif_parentTerm <- merge(go_df, reducedTerms, by.x="Description", by.y="term")
  
  count_go <- go_df_signif_parentTerm %>% dplyr::group_by(celltype, parentTerm, direction ) %>% dplyr::count() %>% dplyr::arrange(n)
  total_terms <- go_df_signif_parentTerm %>% dplyr::group_by(celltype) %>% dplyr::count() %>% dplyr::rename("total" = "n")
  count_go_perc<- count_go %>% left_join(total_terms, by = "celltype")
  count_go_perc <- count_go_perc %>% dplyr::group_by(celltype,direction)%>% mutate(freq = round(n / total * 100))
  count_go_perc$celltype <- factor(count_go_perc$celltype)
  #count_go_perc <- reorder_cells(count_go_perc, neworder = T)
  count_go_perc$clust <- NA
  count_go_perc$clust <- ifelse(count_go_perc$celltype %in% c("MAIT", "CD8 Naive", "B memory", "B naive", "B intermediate", "CD4 Naive"), "Down_bias", "Up_bias")
  
  count_go_perc <- count_go_perc %>%
    mutate(parentTerm = factor(parentTerm, levels = parentTerm[order(celltype)]))
  count_go_perc$direction <- factor(count_go_perc$direction, levels=c("up", "down"))
  count_go_perc <- reorder_cells(count_go_perc, neworder = T)
  plot<- ggplot(count_go_perc, aes(x=celltype, y=parentTerm, color=direction))+geom_point( aes(size=n))+theme +scale_y_discrete(labels = label_wrap(30))+
    scale_size_continuous(name="# GO terms")+ xlab(" ")+ylab("Parent Terms")+ scale_color_manual(values= c("down"=alpha(blue, 0.8), "up"=alpha(red, 0.8)))+
    theme+ggtitle(paste0( "Cell-type specific DEGs"))+facet_grid(~direction, scales = "free_x")+
    theme(axis.text.x=element_text(angle=90, hjust = 0.95, vjust = 0.6), axis.text = element_text(size = 12),plot.title = element_text( face = "bold", hjust = 0.5, size=15) , strip.text=element_text(size=14))
return(list(plot,go_df_signif_parentTerm ))
  }


Fig2D <-plot_reduced_terms(reducedTerms, go_terms)[[1]]

go_df_signif_parentTerm <-plot_reduced_terms(reducedTerms, go_terms)[[2]]

# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig2/Fig2D_enrichments_allGO.pdf"), height =5.76, width= 6.85 )
# Fig2D
# dev.off()


#Save Supplementary Table 2
require(openxlsx)
df_tosave <- list("DEA"=all_genes, "DEG_Sharing"=df_shared, "Enrichment_ShareDEGs1"=go_1, "Enrichment_sharedDEGs2"=go_2, "Enrinchments_SpecificDEGs"=go_df_signif_parentTerm)
write.xlsx(df_tosave, paste0(data_path, '/msopena/02_OneK1K_Age/supplementary_tables/TableS2_DEA.xlsx'))


################ SUMPPLEMENTARY FIGURES ###################

#1. Downsampling results (Fig S1A) --------------
major_cells <- order_cells_old$cell_type[1:13]
major_cells <- gsub(" ", "_", major_cells)
stats <- do.call(rbind.data.frame, lapply(major_cells, function(cell) readRDS(paste0(data_path, "//msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/downsampling/", cell,'_dreamlet_stats_subsampling_1.rds'))))
stats$cell_type <- rownames(stats) 
stats <- stats %>% dplyr::select(cell_type, Tested_samples) 
stats$data <- "All_data"
#ndonors subsampling to MAIT (542 donors) ------
stats_downsampling <- data.frame(cell_type=major_cells, Tested_samples=452)
stats_downsampling$data <- "Downsampling"
stats <- rbind(stats, stats_downsampling)
stats$celltype <- gsub("_", " ", stats$cell_type)
stats <- reorder_cells(stats, neworder = T)

ndonors_plot <- ggplot(stats, aes(y=celltype, x=Tested_samples)) + geom_bar(stat="identity", aes(fill=data), position = "dodge", width = 0.8)+
  theme  + xlab("N donors") + ylab("") +scale_fill_manual(values=c("All_data"= alpha(blue, 0.9), "Downsampling"= alpha(blue, 0.5)))+
  scale_x_continuous(breaks = c(0, 500, 1000))+
  theme(axis.text.x=element_text(size=9), axis.text.y=element_text(size=11), axis.title=element_text(size=11),  legend.position="none",strip.text=element_blank())+ facet_grid(celltype_l1~ ., scales="free", space="free")


number_degs <-readRDS( paste0(data_path,"/msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/downsampling/AllCells_downsamplingTopTable_numbers.rds"))
major_cells <- gsub(" ", "_", major_cells)
names(number_degs) <- major_cells
for(cell in major_cells){number_degs[[cell]] <- do.call(rbind.data.frame, number_degs[[cell]])}
ndegs_df <- do.call(cbind.data.frame, number_degs) %>%
  pivot_longer(cols = everything(), names_to = "celltype", values_to = "value")
ndegs_df$celltype <- gsub("_", " ", ndegs_df$celltype)

mean_n <- ndegs_df %>% group_by(celltype) %>%  dplyr::summarize(value = mean(value, na.rm=TRUE))
ndegs_df$repetition <- "Downsampling_replicate"
mean_n$repetition <- "Downsampling_mean"
#mean_n$Tested_samples <- stats[stats$cell_type == "MAIT", "Tested_samples"]
mean_n$data <- "Downsampling_mean"
mean_n <- reorder_cells(mean_n,neworder = T )
mean_n <- mean_n%>% select(celltype, celltype_l1, data, value)


real_degs <-get_significant_degs("Age", "cell_type") %>% filter(celltype %in% gsub("_", " " , major_cells)) %>% dplyr::group_by(celltype) %>% dplyr::count()
real_degs$data <- "All_data"
#real_degs <-  real_degs %>% left_join(stats, by="celltype")
real_degs <- reorder_cells(real_degs, reverse = T)   
colnames(real_degs) <- gsub("n", "value", colnames(real_degs))

df <- rbind(real_degs, mean_n)
df <- reorder_cells(df, neworder = T)
ndegs_downsampling <- ggplot(df, aes(y=celltype, x=value)) +geom_bar(stat="identity",  aes(fill=data), position = "dodge")+
  xlab("Number of DEGs") + ylab("") +theme(axis.text=element_text(size=14), strip.text=element_text(size=14))+ theme+  scale_x_continuous(breaks = c(0, 1000, 2000))+
  facet_grid(celltype_l1~ ., scales="free", space="free")+theme(legend.position="top", legend.title=element_blank(), axis.text.y=element_blank())+
  scale_fill_manual(values=c(alpha(alpha=0.6, colour = "#264653" ),alpha(alpha=0.9, colour = "#264653" ) ))

library(patchwork)
FigS2A<- ndonors_plot+ndegs_downsampling+plot_layout(widths = c(8, 18))


pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS2/FigS2A_Downsampling.pdf"), heigh =6, width = 4.13 )
FigS2A
dev.off()

# 2. Number of concordant/non concodant (Fig S2B) ------------------
sharing <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/sharing.rds"))
sharing[sharing$n >= 6,]$n <- 6
sh <- sharing %>% dplyr::group_by(n, concordant) %>% dplyr::count() %>% dplyr::group_by(n) %>%
  dplyr::summarise(ratio = nn[concordant == 0] / nn[concordant == 1])

FigS2B <- ggplot(sh, aes(x=n, y=ratio))+geom_point(fill=blue)+geom_smooth(method = "lm", se = FALSE, color=blue)+  stat_cor(method = "spearman",   p.digits = 3,label.y=, size=3.5)+
ylab("non concordant / concordant")+xlab("Number of cell types")+theme+theme(legend.position="none", axis.title=element_text(size=11), axis.text=element_text(size=9.5), aspect.ratio=1)

pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS2/FigS2_CorrelationNSharing.pdf"), width =3, height = 3 )
FigS2B
dev.off()


# 3. Plot clustering of shared genes at several sharing thresholds (n=2:6) (Fig S2C-D) ---------------------
library(ggdendro);library(stringr);library(dendextend)

# Read in sharing information --
sharing <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/sharing_02.rds"))
non_concordant <- sharing %>% dplyr::filter(category == "non concordant")

# Read in DEG results ---
all_genes <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/DEG_Age_cell_type_02.rds"))
all_genes$celltype <- factor(all_genes$celltype, levels=rev(order_cells$cell_type))

# Get siginificant DEGs
degs <-get_significant_degs("Age","cell_type")
cells <- unique(degs[!degs$celltype %in% c("NK Proliferating", "gdT", "NK_CD56bright", "Platelet", "HSPC", "cDC2"),]$celltype)
# Select cell types to plot 
df <- all_genes[all_genes$celltype %in% cells,] %>% as.data.frame()
df <- reorder_cells(df)

# Get the common tested genes across celltypes 
list <-  split(df$gene, df$celltype)
common_genes <- Reduce(intersect, list)

# Extract a matix of logFCs across celltypes 
df <- df[df$gene %in% common_genes,]
df$significance <- ifelse(df$adj.P.Val < 0.05, "·", " ")
logfc <- data.frame("gene"=common_genes)
for (celltype in levels(df$celltype)){
  cell <- df[df$celltype == celltype, c("gene", "logFC")]
  colnames(cell) <- gsub("logFC", paste0(celltype), colnames(cell))
  print(head(cell))
  print(dim(cell))
  logfc <- merge(logfc, cell, by="gene")
}
rownames(logfc) <- logfc$gene
logfc <- as.matrix(logfc[,-1])


plot_dendrogram <- function(n, method, nonconcordant=T){
  if(n == 1){
    logfc_subset <- logfc
  }else{  logfc_subset <- logfc[rownames(logfc) %in% sharing[sharing$n >= n,]$gene,]
  if(nonconcordant == T){
    logfc_subset <- logfc_subset[rownames(logfc_subset) %in% non_concordant$gene,]
  }}
  clust <- hclust(dist(t(logfc_subset)), method)
  dhc <- dendro_data(clust)
  labels <- label(dhc)
  labels$celltype <- labels$label
  labels <- reorder_cells(labels)
  labels$point <- "."
  labels$label <- str_pad(labels$label, width = 14, side = "right")
  labels$label <- str_pad(labels$label, width = 17, side = "left")
  #labels$label <- ifelse(labels$label %in% c("     NK            ", "     Treg          ", "     B naive       " , "     MAIT          "), paste0("", labels$label), labels$label)
  d <- ggplot(data = labels,  aes(x = x, y = y)) +geom_segment(data = segment(dhc), aes(x = x, y = y, xend = xend, yend = yend)) +geom_text(aes(label = label, hjust=0),  size = 2) +
    geom_point( size = 2, aes(color=celltype_l1)) +
    coord_flip() +scale_y_reverse(expand = c(0.1, 0))+theme+theme(axis.text= element_blank(), axis.line=element_blank(), panel.border=element_blank(), axis.ticks=element_blank())+
    ylab("")+xlab("")+ggtitle(paste0("n=", n))  +theme(plot.title = element_text(hjust =0.5), legend.position="none")+scale_color_manual(values= c("CD8 T"= "#679267", "CD4 T"= "#00563f",  "other T"= "#c0d3ac", "NK"= "#E76F51", "B"= "#8580C2", "Mono" ="#d88d8b" ),)
  return(list(clust, d))
}


plots <- lapply(c(1, 2, 3, 4, 5, 6, 7,8,9), function(n) plot_dendrogram(n, "ave", T)[[2]])
FigS2D <- cowplot:: plot_grid(plotlist = plots, ncol = 3) 


pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS2/FigS2D_Dendrograms_04.pdf"), width =20, height = 9 )
FigS2D
dev.off()



#  Compare similarity of dendrograms at several n thresholds --
compare_dendrograms_n <- function(method, nonconcordant=T){
  #list <- lapply(c(1:10), function(n), plot_dendrogram(n, method , nonconcordant)[[1]])
  #n1 <- as.dendrogram(plot_dendrogram(1, method , nonconcordant)[[1]])
  n2 <- as.dendrogram(plot_dendrogram(2, method , nonconcordant)[[1]])
  n3 <- as.dendrogram(plot_dendrogram(3, method , nonconcordant)[[1]])
  n4 <- as.dendrogram(plot_dendrogram(4, method , nonconcordant)[[1]])
  n5 <- as.dendrogram(plot_dendrogram(5, method , nonconcordant)[[1]])
  n6 <- as.dendrogram(plot_dendrogram(6, method , nonconcordant)[[1]])
  n7 <- as.dendrogram(plot_dendrogram(7, method , nonconcordant)[[1]])
  n8 <- as.dendrogram(plot_dendrogram(8, method , nonconcordant)[[1]])
  n9 <- as.dendrogram(plot_dendrogram(9, method , nonconcordant)[[1]])
  n10 <- as.dendrogram(plot_dendrogram(9, method , nonconcordant)[[1]])
  dist_list <- dendlist("DEGs_2_cells" =n2, "DEGs_3_cells"= n3, "DEGs_4_cells"=n4, "DEGs_5_cells"=n5, "DEGs_6_cells"=n6, "DEGs_7_cells"=n7, "DEGs_8_cells"= n8, "DEGs_9_cells"=n9, "DEGs_10_cells"=n10)
  mat <- as.matrix(cor.dendlist(dist_list))
  mat[lower.tri(mat)] <- NA
  p <- ggcorrplot::ggcorrplot(mat,  method = "square", lab = T)+scale_fill_gradient2(low="white", high = "#A93F41",midpoint = 0.6 )+theme+
    theme(axis.text.x=element_text(angle=45, vjust = 0.5))+xlab("")+labs(x="", y="", fill="correlation")
  return(p)
}
FigS2C <- compare_dendrograms_n("average", T)

pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS2/FigS2C_CorrelationDendrograms_02.pdf"), width =6, height = 6 )
FigS2C
dev.off()


#4. Enrichments of shared-DEGs across n cell types (Figs S2E, S2F) --------------

get_enrichments_cluster_genes<-function(n,nclust=2,  equal=T){ 
  print(n)
  if(equal==T){logfc_na_shar <- logfc[rownames(logfc) %in% sharing[sharing$n == n,]$gene,]}else{
    logfc_na_shar <- logfc[rownames(logfc) %in% sharing[sharing$n > n,]$gene,]; logfc_na_shar <- logfc_na_shar[rownames(logfc_na_shar) %in% non_concordant$gene,]}
  set.seed(10)
  heatmap <- pheatmap::pheatmap(t(logfc_na_shar))
  rw <- heatmap$tree_col
  clst <- cutree(rw, nclust)
  
  return(list( names(clst[clst==1]),  names(clst[clst==2])))}

# # Two clusters 
# list <- lapply(c(2:8), function(n)barplot(clusterProfiler::enrichGO(gene = get_enrichments_cluster_genes(n, equal = F)[[1]], OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", universe = unique(all_genes$gene)), label_format=30, color="p.adjust")+theme+ggtitle(paste0("DEGs shared > ", n, " cell types")))
# ggarrange(plotlist=list)
# 
# list <- lapply(c(2:8), function(n)dotplot(clusterProfiler::enrichGO(gene = get_enrichments_cluster_genes(n, equal = F)[[2]], OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", universe = unique(all_genes$gene)), label_format=30, color="p.adjust")+theme+ggtitle(paste0("shared > ", n, " cell types (", length(get_enrichments_cluster_genes(n, equal = F)[[2]]), ")" )))
# ggarrange(plotlist=list)

library(rrvgo)
clst1 <- lapply(c(2:8), function(n)clusterProfiler::enrichGO( gene = get_enrichments_cluster_genes(n, equal = F)[[1]], OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", universe = unique(all_genes$gene))@result)

names(clst1) <- as.character(c(2:8))
for(n in names(clst1)){clst1[[n]]$n <- n}
clst1_df <- do.call(rbind.data.frame, clst1)
clst1_df <- clst1_df[clst1_df$p.adjust < 0.05,]

library(forcats)

clst1_df$Description <- fct_rev(fct_infreq(clst1_df$Description))

ggplot(clst1_df, aes(x=n, y=Description))+geom_point(color="#777777", aes(size=Count, alpha=p.adjust))+theme+scale_size(range = c(3, 6))+
  ylab("")+xlab("Number of cells DEGs")+ggtitle("Cluster 1- Translation process")+ scale_y_discrete(labels = label_wrap(60))+
  theme(axis.text.y=element_text(size=8), aspect.ratio=3, axis.title.x=element_text(size=11), plot.title = element_text(hjust = 0.5, size=12))

clst1_count <- clst1_df %>% group_by(Description) %>% count() %>% filter(n>2)

FigS2E <- ggplot(clst1_df[clst1_df$Description %in% clst1_count$Description,], aes(x=n, y=Description))+geom_point(color="#777777", aes(size=Count, alpha=p.adjust))+theme+scale_size(range = c(3, 6))+
  ylab("")+xlab("Number of cells DEGs")+ggtitle("Cluster 1- Translation process")+ scale_y_discrete(labels = label_wrap(50))+
  theme(axis.text.y=element_text(size=8), aspect.ratio=3, axis.title.x=element_text(size=11), plot.title = element_text(hjust = 0.5, size=12))




clst2 <- lapply(c(2:8), function(n)clusterProfiler::enrichGO(gene = get_enrichments_cluster_genes(n, equal = F)[[2]], OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", universe = unique(all_genes$gene))@result)

names(clst2) <- as.character(c(2:8))
for(n in names(clst2)){clst2[[n]]$n <- n}
clst2_df <- do.call(rbind.data.frame, clst2)
clst2_df <- clst2_df[clst2_df$p.adjust < 0.05, ]
library(forcats)

clst2_df$Description <- fct_rev(fct_infreq(clst2_df$Description))


ggplot(clst2_df, aes(x=n, y=Description))+geom_point(color="#777777", aes(size=Count, alpha=p.adjust))+theme+scale_size(range = c(3, 6))+
  ylab("")+xlab("Number of cells DEGs")+ggtitle("Cluster 2- Immunological proces")+ scale_y_discrete(labels = label_wrap(60))+
  theme(axis.text.y=element_text(size=8), aspect.ratio=3, axis.title.x=element_text(size=11))

clst2_count <- clst2_df %>% group_by(Description) %>% count() %>% filter(n>1)

FigS2F <- ggplot(clst2_df[clst2_df$Description %in% clst2_count$Description,], aes(x=n, y=Description))+geom_point(color="#777777", aes(size=Count, alpha=p.adjust))+theme+scale_size(range = c(3, 6))+
  ylab("")+xlab("Number of cells DEGs")+ggtitle("Cluster 2- Immunological proces")+ scale_y_discrete(labels = label_wrap(50))+
  theme(axis.text.y=element_text(size=8), aspect.ratio=3, axis.title.x=element_text(size=11), plot.title = element_text(hjust = 0.5, size=12))



pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS2/FigS2E_Enrichment1.pdf"), width =7.13, height = 6)
FigS2E
dev.off()

pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS2/FigS2D_Enrichment2.pdf"), width =5.85, height  = 3.70 )
FigS2F
dev.off()

FigS2E /FigS2F +plot_layout(heights =c(5, 2), widths = c(3, 5))




#5. Replication Oelen et al., 2020 (Figs S2G, S2E) ------------
# Dendrogram 
#read tested genes 
all_genes <- readRDS(paste0(data_path, "msopena//01_meta-analysis_Age_Sex/robjects/02_DEA/DEG_Age_cell_type_v2_0.2.rds"))
all_genes$gene <- all_genes$ID
all_genes <- all_genes[!all_genes$celltype %in% c("Platelet", "Plasmablast", "NK_CD56bright", "NK Proliferating", "HSPC", "gdT", "Eryth"),]
table(all_genes$celltype, all_genes$adj.P.Val<0.05)
all_genes$cell_type <- all_genes$celltype
tested_genes <- split(all_genes$gene, all_genes$celltype)
common_genes <- Reduce(intersect, tested_genes)

# get significant DEGs at p nominal 
degs <- all_genes[all_genes$P.Value < 0.05,]
degs <- degs[degs$celltype != "Platelet",]
degs$gene <- degs$ID


#sharing 
df <- degs  %>% dplyr::filter(!celltype %in% c(  "NK_CD56bright", "gdT", "Platelet")) %>% dplyr::group_by(gene) %>% dplyr::count() %>% as.data.frame()
df$x <- "DEGs"
df$label <- NA
df[df$n > 10,]$label <- df[df$n > 10,"gene"]
df$concordant <- NA
non_concordant <- df[df$gene %in% degs[degs$direction == "up",]$gene & df$gene %in% degs[degs$direction == "down",]$gene ,]$gene
df[df$n > 1, ]$concordant <- ifelse(  df[df$n > 1, ]$gene %in% non_concordant, 0, 1)
df$count<- df$n
df$group <- ifelse(df$count == 1, "Specific", "Shared" )
df_perc <- df %>%   dplyr::group_by(group) %>% dplyr::summarise(cnt = n()) %>% mutate(freq = round(cnt / sum(cnt), 3))
df_perc$group <- factor(df_perc$group, levels= rev(c("Specific", "Shared" )))
df_perc$x <-"DEGs"
df_shared <-df[!is.na(df$concordant),] 
df_shared$concordant <- factor(df_shared$concordant, levels=c(1, 0))

sharing <- df_shared
saveRDS(df_shared, paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/sharing_v2_02.rds"))


# sharing 

compare_v2_dendrogram <- function(all_genes, sharing,n){
  cells <- unique(degs[!degs$celltype %in% c("NK Proliferating", "NK CD56bright", "Platelet", "HSPC", "CD4 CTL", "dnT"),]$celltype)
  # Select cell types to plot 
  df <- all_genes[all_genes$celltype %in% cells,] %>% as.data.frame()
  #df <- reorder_cells(df)
  
  # Get the common tested genes across celltypes 
  list <-  split(df$gene, df$celltype)
  common_genes <- Reduce(intersect, list)
  
  # Extract a matix of logFCs across celltypes 
  df <- df[df$gene %in% common_genes,]
  logfc <- data.frame("gene"=common_genes)
  for (celltype in unique(df$celltype)){
    cell <- df[df$celltype == celltype, c("gene", "logFC")]
    colnames(cell) <- gsub("logFC", paste0(celltype), colnames(cell))
    print(head(cell))
    print(dim(cell))
    logfc <- merge(logfc, cell, by="gene")
  }
  rownames(logfc) <- logfc$gene
  logfc <- as.matrix(logfc[,-1])
  
  
  # Extract non concordant genes DE in more than 6 cell types 
  logfc_na_shar <- logfc[rownames(logfc) %in% sharing[sharing$n > n,]$gene,]
  logfc_na_shar <-  logfc_na_shar[,colnames(logfc_na_shar) !="pDC"]
  #logfc_na_shar <- logfc_na_shar[rownames(logfc_na_shar) %in% non_concordant$gene,]
  
  
  # Get column annotation 
  col_annotation <- celltype_l1[celltype_l1$cell_type %in% colnames(logfc_na_shar),]
  col_annotation <- col_annotation[match(colnames(logfc_na_shar), col_annotation$cell_type),]
  rownames(col_annotation) <-col_annotation$cell_type 
  col_annotation <- col_annotation %>% dplyr::select(predicted.celltype.l1) %>% as.data.frame()
  colnames(col_annotation) <- "cell_type_l1"  
  #annot_colors <-list(cell_type_l1= c("CD4 T"= "#00563f",  "CD8 T"= "#679267","other T"= "#c0d3ac", "NK"= "#E76F51", "B"= "#8580C2", "Mono" ="#d88d8b" ),clusters=c("1"=alpha("#264653", 0.8), "2"="#bc4749"))
  
  set.seed(3)
  #row_clusters <- hclust(dist(logfc_na_shar), method = "average")
  #col_clusters <- hclust(dist(t(logfc_na_shar)), method ="average")
  
  callback = function(hc, mat){
    sv = svd(t(mat))$v[,1]
    dend = reorder(as.dendrogram(hc), wts = sv, method = "average")
    as.hclust(dend)
  }
  
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...), type = "average"))
  sort_hclust_cols <- function(...) as.hclust(dendsort(as.dendrogram(...), type = "ward.D2"))
  
  heatmap <- pheatmap::pheatmap(t(logfc_na_shar),  clustering_method = "average")
  
  # get gene clusters and change names 
  rw <- heatmap$tree_col
  clusters_row <- cutree(rw, 2) %>% as.data.frame() %>% dplyr::rename("enrichments"=".")
  clusters_row$enrichments <- as.factor(clusters_row$enrichments)
  clusters_row$enrichments <- gsub("1", "Ribosomal_function",clusters_row$enrichments)
  clusters_row$enrichments <- gsub("2", "Immunological_process",clusters_row$enrichments)
  
  # get cell clusters and change names 
  
  annot_colors <-list(enrichments= c("Ribosomal_function"="lightgrey", "Immunological_process"="#777777"),bias=c("downregulation_bias"=alpha("#264653", 0.8), "upregulation_bias"="#bc4749"))
  # annot_colors <-list(Genes_clust= c("Ribosomal_function"="lightgrey", "Immunological_process"="#777777"),Binomial_p.val=annotation_colors_for_row)
  
  # Heatmap
  heatmap_horizontal <- pheatmap::pheatmap(t(logfc_na_shar),
                                           method = "average",
                                           cluster_rows = sort_hclust(hclust(dist(t(logfc_na_shar)), method = "average")),
                                           annotation_names_row = FALSE,
                                           annotation_names_col = FALSE,
                                           treeheight_col = 10,treeheight_row = 10,
                                           fontsize_number = 18, fontsize_col = 6, fontsize_row = 10,
                                           border_color = NA,
                                           color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(20),
                                           annotation_colors = annot_colors
  )
  
  clust <-   heatmap_horizontal$tree_row
  #clust <- hclust(dist(t(logfc_subset)), "ave")
  dhc <- dendro_data(clust)
  labels <- label(dhc)
  labels$celltype <- labels$label
  labels <- reorder_cells(labels)
  labels$point <- "."
  labels$label <- str_pad(labels$label, width = 14, side = "right")
  labels$label <- str_pad(labels$label, width = 17, side = "left")
  #labels$label <- ifelse(labels$label %in% c("     NK            ", "     Treg          ", "     B naive       " , "     MAIT          "), paste0("", labels$label), labels$label)
  d <- ggplot(data = labels,  aes(x = x, y = y)) +geom_segment(data = segment(dhc), aes(x = x, y = y, xend = xend, yend = yend)) +geom_text(aes(label = label, hjust=0),  size = 4) +
    geom_point( size = 3, aes(color=celltype_l1)) +
    coord_flip() +scale_y_reverse(expand = c(0.5, 0))+theme+theme(axis.text= element_blank(), axis.line=element_blank(), panel.border=element_blank(), axis.ticks=element_blank())+
    ylab("")+xlab("")+ggtitle(paste0("n=", n))  +theme(plot.title = element_text(hjust =0.5), legend.position="none")+scale_color_manual(values= c("CD8 T"= "#679267", "CD4 T"= "#00563f",  "other T"= "#c0d3ac", "NK"= "#E76F51", "B"= "#8580C2", "Mono" ="#d88d8b" ),)
  
  
  return(list(heatmap_horizontal, clust, d))
}

FigS2G <- compare_v2_dendrogram(all_genes, sharing,2)[[3]]

pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS2/FigS2G_DendrogramV2_04.pdf"), width =6, height = 3 )
FigS2G
dev.off()


# enrichment shared genes -------

sharing_onek <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/sharing.rds"))

# Read in DEG results ---
all_genes <- readRDS(paste0(data_path, "msopena//01_meta-analysis_Age_Sex/robjects/02_DEA/DEG_Age_cell_type_v2_0.2.rds"))
all_genes$gene <- all_genes$ID
all_genes <- all_genes[!all_genes$celltype %in% c("Platelet", "Plasmablast", "NK_CD56bright", "NK Proliferating", "HSPC", "Eryth"),]
tested_genes <- split(all_genes$gene, all_genes$celltype)
common_genes <- Reduce(intersect, tested_genes)
length(common_genes)
sharing <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/sharing_v2_02.rds"))

library(clusterProfiler)
oelen_enrich <- enrichGO(gene = sharing[sharing$n > 2,]$gene, universe=unique(all_genes$gene),OrgDb  = "org.Hs.eg.db", keyType = "SYMBOL", ont="BP")

# reduce terms 
library(rrvgo)
simMatrix <- calculateSimMatrix(oelen_enrich$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
)
scores <- setNames(-log10(oelen_enrich$qvalue), oelen_enrich$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.85,
                                orgdb="org.Hs.eg.db")



go_df_signif_parentTerm <- merge(oelen_enrich, reducedTerms, by.x="Description", by.y="term")
count_go_perc <- go_df_signif_parentTerm %>% dplyr::group_by(parentTerm) %>% dplyr::count() %>% dplyr::arrange(n)
total_terms <- nrow(oelen_enrich)
count_go_perc$total_terms<-  nrow(oelen_enrich)
count_go_perc$freq <- count_go_perc$n / count_go_perc$total_terms
#count_go_perc <- reorder_cells(count_go_perc, neworder = T)
count_go_perc$clust <- NA

FigS2H <- ggplot(count_go_perc, aes(x= freq, y=reorder(parentTerm, freq)))+geom_point( aes(size=n), color=blue)+theme +scale_y_discrete(labels = label_wrap(30))+
  scale_size_continuous(name="# GO terms")+ xlab(" ")+ylab("Parent Terms")+xlab("Percentage of terms")+
  theme+ggtitle(paste0( "Shared DEGs (n>2)"))+
  theme( axis.text = element_text(size = 11),plot.title = element_text( face = "bold", hjust = 0.5, size=12) , strip.text=element_text(size=14))


pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS2/FigS2H_EnrichmentsOelen_02.pdf"), width =5.28, height = 3.68 )
FigS2H
dev.off()





