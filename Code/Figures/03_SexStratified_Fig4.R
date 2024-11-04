

# Plots results of the Sex stratified analysis (Figure 5 and Supplemental Figure 5 )

library(ggplot2); library(dplyr); library(RColorBrewer);library(ggside);library(ggbeeswarm);library(ggpubr);library(ggrepel);library(tidyr); library(patchwork);library(scales)

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
path_prop <- paste0(plots_path, "/08_CellProps/")
dir.create(path_prop, recursive = TRUE)

#functions
source(paste0(data_path, "/msopena/02_OneK1K_Age/scripts/functions.R"))

#theme and palettes
computer <- "work"
source(paste0(data_path, "/msopena/01_meta-analysis_Age_Sex/scripts/themes.R"))


#data
metadata <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/metadata_processed.rds"))
order_cells<- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/new_order_cells.rds"))
celltype_l1 <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/celltypes_equivalence.rds"))
celltype_l1$cell_type <- factor(celltype_l1$cell_type, levels = order_cells$cell_type)
cells_to_keep <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/cells_to_keep.rds"))


# 0. Basic objects for the analysis-------------
deg_M <- get_significant_degs_sex("M",F) %>% tidyr::drop_na()
deg_F<- get_significant_degs_sex("F",F) %>% tidyr::drop_na()

# get overlap M and F
overlap <- lapply(unique(c(deg_F$cell_type, deg_M$cell_type)), function(cell) intersect(deg_M[deg_M$cell_type == cell,]$gene, deg_F[deg_F$cell_type == cell,]$gene ))
names(overlap) <- unique(c(deg_F$cell_type, deg_M$cell_type))

deg_F$overlap  <- NA
for(cell in unique(deg_F$celltype)){
  deg_F[deg_F$celltype == cell,]$overlap <- unlist(  lapply(deg_F[deg_F$celltype == cell,]$gene, function(gene) ifelse(gene %in% overlap[[cell]], "overlap", "only_Females")))
}
table(deg_F$overlap, deg_F$celltype)
table(deg_F$overlap)

deg_M$overlap  <- NA
for(cell in unique(deg_M$celltype)){
  deg_M[deg_M$celltype == cell,]$overlap <- unlist(  lapply(deg_M[deg_M$celltype == cell,]$gene, function(gene) ifelse(gene %in% overlap[[cell]], "overlap", "only_Males")))
}
table(deg_M$overlap, deg_M$celltype)
table(deg_M$overlap)


# Get overlap interaction 
only_F <- deg_F[deg_F$overlap== "only_Females",]
up_F <- only_F %>% dplyr::filter(direction == "up") %>% dplyr::group_by(gene) %>% dplyr::count()
down_F <- only_F %>% dplyr::filter(direction == "down") %>% dplyr::group_by(gene) %>% dplyr::count()
#saveRDS(only_F, paste0(data_path, "msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/DEG_genes_only_F.rds"))

only_M <- deg_M[deg_M$overlap== "only_Males",]
up_M <- only_M %>% dplyr::filter(direction == "up") %>% dplyr::group_by(gene) %>% dplyr::count()
down_M <- only_M %>% dplyr::filter(direction == "down") %>% dplyr::group_by(gene) %>% dplyr::count()
#saveRDS(only_M, paste0(data_path, "msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/DEG_genes_only_M.rds"))
overlap <- deg_M[deg_M$overlap == "overlap", ]


interaction <- readRDS( paste0(data_path, "msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/AllCells_Interaction_deaTopTable.rds"))
interaction$direction <- ifelse(interaction$logFC > 0, "up", "down")
interaction$celltype <- interaction$assay
interaction$cell_type<- gsub("_", " ", interaction$cell_type)
interaction <- interaction[interaction$P.Value < 0.01,]
interaction <- reorder_cells(interaction, reverse = T, neworder = T)
interaction$classification <- "not_DE"
interaction$classification <- unlist(lapply(unique(as.character(interaction$celltype)), 
                                            function(cell) ifelse(interaction[interaction$celltype == cell,]$gene %in% only_F[only_F$celltype == cell,]$gene, "DE_F", 
                                                                  ifelse(interaction[interaction$celltype == cell,]$gene %in% only_M[only_M$celltype == cell,]$gene, "DE_M",
                                                                         ifelse(interaction[interaction$celltype == cell,]$gene %in% overlap[overlap$celltype == cell,]$gene, "DE_both", "no_DE")))))

only_M$classification <- unlist(lapply(unique(as.character(only_M$celltype)), 
                                       function(cell) ifelse(only_M[only_M$celltype == cell,]$gene %in% interaction[interaction$celltype == cell & interaction$classification == "DE_M",]$gene, "interaction", "not_interaction")))

only_F$classification <- unlist(lapply(unique(as.character(only_F$celltype)), 
                                       function(cell) ifelse(only_F[only_F$celltype == cell,]$gene %in% interaction[interaction$celltype == cell & interaction$classification == "DE_F",]$gene, "interaction", "not_interaction")))

only_M_interaction <- only_M[only_M$classification == "interaction",]
only_F_interaction <- only_F[only_F$classification == "interaction", ]



#1. Plot DEA results per sex (Fig 4A)  ---------

plot_direction_DEG_sex <- function(balanced=F ){
  ncells <- metadata %>% dplyr::group_by(cell_type, Gender) %>% dplyr::count() %>% drop_na()
  colnames(ncells)[1] <- "celltype"
  #df <- df %>% left_join(ncells, by="celltype") 
  ncells$sex <- ncells$Gender
  df_F <- get_significant_degs_sex("M",F) %>% tidyr::drop_na()
  df_M <- get_significant_degs_sex("F",F) %>% tidyr::drop_na()
  df <- rbind(df_F, df_M)
  df <- df[!df$celltype %in% c("HSPC","Platelet"),]
  df$celltype <- df$cell_type
  ncells <- ncells[ncells$celltype %in% df$celltype,] %>% as.data.frame()
  ncells <- reorder_cells(ncells, reverse = T, neworder = T)
  #df <- df[!df$celltype %in% c("NK Proliferating" , "NK_CD56bright", "Platelet", "gdT" ),]
  df <- df %>% group_by(celltype, direction, sex) %>% dplyr::count()
  
  # Compute the bias towards up/down regulation 
  result_df <- data.frame(celltype = character(0), binom.test.p = numeric(0))
  for (sex in c("M", "F")){
    for (cell in unique(df[df$sex == sex,]$celltype)){
      print(cell)
      binom <- binom.test(length(degs[degs$direction == "down" & degs$celltype == cell & degs$sex == sex,]$gene) , length(degs[degs$celltype == cell &degs$sex == sex,]$gene))$p.value
      result_df <- rbind(result_df, data.frame(celltype = cell, binom.test.p = binom, sex= sex))
    }
  }
  
  df <- merge(df, result_df, by=c("celltype", "sex"))
  df$binom.test.fdr <- p.adjust(df$binom.test.p, method = "fdr")
  df$signif <- ifelse(df$binom.test.fdr < 0.05, "*", " ")
  #saveRDS(df, paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/biomial_test_bias_sex.rds"))
  
  for(sex in c("M", "F")){
    for (cell in unique(df[df$sex == sex,]$celltype)){
      if(length(df[df$celltype == cell & df$sex== sex,]$celltype) == 2  && df[df$celltype == cell & df$direction == "down"& df$sex== sex,]$signif == "*"){
        df_cell <- df[df$celltype == cell& df$sex== sex,,]
        if(df_cell[df_cell$direction == "up",]$n < df_cell[df_cell$direction == "down",]$n){
          df[df$celltype == cell & df$direction == "up" & df$sex== sex, ]$signif <- ""
        } else{
          df[df$celltype == cell & df$direction == "down"& df$sex== sex, ]$signif <- ""}}}
    
  }
  
  df$sex <- ifelse(df$sex=="M", "Males", "Females")
  up <- df[df$direction == "up", ]
  up <- reorder_cells(up, reverse = T, neworder = T)
  down <- df[df$direction == "down", ]
  down <- reorder_cells(down, reverse = T, neworder = T)
  up$signif_n <- paste0(up$n, " ", up$signif)
  down$signif_n <- paste0(down$signif, " ", down$n)
  p_ndegs<-
    ggplot(df, aes(x = celltype, y = n, fill = sex)) +
    geom_col(data = up, aes(x = celltype, y = n, fill = sex), alpha = 1, position = "dodge") +
    geom_col(data = down, aes(x = celltype, y = -n, fill = sex), alpha = 1, position = "dodge") +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_text(data = down, aes(label = n, x = celltype, y = -n), vjust = 0.5, hjust = 1, size = 3, position = position_dodge(width = 1)) +
    geom_text(data = up, aes(label = n, x = celltype, y = n), hjust = -0, size = 3, position = position_dodge(width = 1)) +
    geom_text(data = down, aes(label = signif, x = celltype, y = -500), vjust = 0.7, hjust = 2, size = 6, position = position_dodge(width = 1.5)) +
    geom_text(data = up, aes(label = signif, x = celltype, y = 500),vjust = 0.7, hjust = -1, size = 6, position = position_dodge(width = 1.5)) +
    ylab("Number of DEGs") +
    xlab("") +
    #scale_y_continuous(labels = abs, limits = c(-1000, 1000), breaks = c(-1000, -100, 0, 100, 1000), trans = pseudo_log_trans(base = 10)) +
    scale_y_continuous(labels = abs, limits = c(-2100, 1500), breaks = c(-2000, -1000, 0, 1000, 2000)) +
    scale_x_discrete(limits = rev(levels(df$celltype))) +
    facet_grid(celltype_l1 ~ sex, scales = "free", space = "free") +
    theme+
    theme(
      axis.text = element_text(size = 11),
      axis.text.x = element_text(size = 10),
      legend.position = "none",
      strip.background = element_blank(),
      legend.margin = margin(r = 10, l = 5, t = 5, b = 2),
      legend.key.size = unit(15, "pt"),
      legend.title=element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    coord_flip() +
    scale_fill_manual(values = c("Males" = red, "Females" = blue)) +
    guides(  color=guide_legend(override.aes = list(size=1)))
  p_ncells <- ggplot(ncells, aes(y=celltype, x=n)) + geom_bar(stat="identity", aes(fill=sex), position = "dodge", alpha=0.6, width = 0.8)+  scale_x_continuous(labels = function(x) format(x / 1e6, scientific = FALSE), limits = c(0,max(ncells$n)),breaks =c(0, 150e3) ) +
    theme  + xlab("N cells\n(thousands)") + ylab("") +scale_fill_manual(values=c("M"= red, "F"= blue))+
    theme(axis.text.x=element_text(size=9), axis.text.y=element_text(size=11), axis.title=element_text(size=11),  legend.position="none",strip.text=element_blank())+ facet_grid(celltype_l1~ ., scales="free", space="free")
  
  fig_5A<- p_ncells+plot_spacer()+p_ndegs+plot_layout(widths = c(5, -2.5, 22))
  
  
  return(fig_5A)
}

fig5A<- plot_direction_DEG_sex(F)

# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig5/Fig5A_DEA_results.pdf"), height =6.18, width= 5.42   )
# fig5A
# dev.off()


# 2. CODA results per sex  (Fig 4B) ---------

coda_M <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/12_CellProportions_Sex/CODA_resultsAge_M.rds"))
coda_M$sex <- "Males"
coda_F <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/12_CellProportions_Sex/CODA_resultsAge_F.rds"))
coda_F$sex <- "Females"
coda <- rbind(coda_M, coda_F)
coda$significance <- ifelse(coda$fdr< 0.05, "ss","ns")
coda$celltype <- gsub("_", " ", coda$celltype)
coda$celltype <- gsub(" CD56bright", "_CD56bright", coda$celltype)
#coda <- coda[coda$celltype %in% keep,]
md <- metadata%>% tidyr::drop_na("cell_type")
count_cells <- md %>% dplyr::group_by(cell_type) %>% dplyr::count() %>% as.data.frame()
count_cells <- count_cells[match(coda$celltype,count_cells[,1]),]
coda$n <- count_cells$n
coda <- reorder_cells(coda, reverse = T, neworder = T)
#coda <- coda[!coda$celltype_l1 %in% c("DC", "other"), ]
coda$direction <- ifelse(coda$estimate < 0, "down", "up")
coda <- coda[!coda$celltype_l1 %in% c( "other"), ]
coda <- coda[coda$celltype %in% unique(deg_F$celltype), ]
fig5B <- ggplot(coda, aes(x=estimate, y=celltype)) + 
  geom_point(aes(alpha=significance, fill=sex, color=sex), size=4) + xlab("Estimate")+ylab(NULL)+
  geom_pointrange(aes(xmin=conf.low, xmax=conf.high, alpha=significance, fill=sex, color=sex), fatten = .1) +
  geom_vline(xintercept=0, linetype = "dashed")  + scale_fill_manual(values=c("Males"=red, "Females"=blue))+scale_color_manual(values=c("Males"=red, "Females"=blue))+
  theme + theme(axis.text = element_text(size = 11), axis.title=element_text(size=12), legend.position="top", legend.title=element_blank(),legend.key.size = unit(0.5, "lines")) +
  scale_alpha_manual(values=c(0.4, 1)) +facet_grid(celltype_l1~sex, space="free", scale="free_y")+ scale_x_continuous(limits = c(-0.022, 0.025), breaks = seq(-0.02, 0.02, by = 0.02))


# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig5/Fig5B_CODA.pdf"), height =6.18, width= 4.25   )
# fig5B
# dev.off()


# 3. Correlation bias estimate and CoDA (Fig 4C) ------------

coda$sex2 <- ifelse(coda$sex == "Males", "M", "F")

get_prop_up <- function (sex) {
  degs <- get_significant_degs_sex(sex) %>% tidyr::drop_na()
  degs <- degs[!degs$celltype %in% c("gdT", "NK Proliferating", "Platelet"),]
  
  result_df <- data.frame(celltype = character(0), prop_up = numeric(0))
  for (cell in unique(degs$celltype)){ 
    print(cell)
    upregulated <- length(degs[degs$direction == "up" & degs$celltype == cell,]$gene)
    downregulated <- length(degs[degs$direction == "down" & degs$celltype == cell,]$gene)
    
    success_ratio = upregulated / (downregulated+upregulated)
    # Store the odds ratio for this cell type
    result_df <- rbind(result_df, data.frame(celltype = cell, prop_up = success_ratio))
  }
  
  bias_estimate <- merge(coda[coda$sex2==sex,], result_df, by="celltype")
  bias_estimate$significance <- ifelse(bias_estimate$fdr <0.05,"ss", "ns")
  return(bias_estimate)
}

bias_estimate_M <- get_prop_up("M")
bias_estimate_F <- get_prop_up("F")

bias_estimate_sex <- rbind(bias_estimate_M, bias_estimate_F)


binomial <- readRDS( paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/biomial_test_bias_sex.rds"))
binomial$sex2 <- binomial$sex

bias_estimate_sex <- merge(bias_estimate_sex, binomial, by=c("sex2", "celltype"))
bias_estimate_sex_M <- bias_estimate_sex[bias_estimate_sex$sex2=="M",]
bias_estimate_sex_F <- bias_estimate_sex[bias_estimate_sex$sex2=="F",]


fig5C <- ggplot(bias_estimate_sex ,aes(y=prop_up, x=estimate))+geom_point(aes(alpha=significance, size=-log10(binom.test.fdr),color=sex.x))+theme+  geom_smooth(method = "lm", se = FALSE, aes(color=sex.x))+ylab("Proportion up-regulated DEGs")+  
  scale_size(name = "-log10(binom.fdr)", breaks = c(10, 30, 50))+facet_grid(sex.x~.)+geom_text_repel(data=bias_estimate_sex_F[!duplicated(bias_estimate_sex_F$celltype),], aes( label=celltype))+geom_text_repel(data=bias_estimate_sex_M[!duplicated(bias_estimate_sex_M$celltype),], aes( label=celltype))+
  scale_alpha_manual(values=c("ss"=1, "ns"=0.3), name="signif")+xlab("CoDA estimate")+
  theme(aspect.ratio=0.95, legend.position="top", plot.title = element_text( size=15, hjust = 0.5) , legend.box = "vertical", legend.spacing = unit(0.000001, "cm"),  legend.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"), strip.text = element_text(size = 13))+
  scale_color_manual(values=c("Males"=red, "Females"=blue), guide="none")+stat_cor(method = "spearman", label.y = 1.15, label.x = -0.012,   p.digits = 1)+  guides(
    size = guide_legend(order = 1),   color = "none" ,
    alpha = guide_legend(order = 2), color = guide_legend(override.aes = list(size = 0.5))
  )


# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig5/Fig5C_CorrelationDEA_CODA.pdf"), height =6.18, width= 4.25   )
# fig5C
# dev.off()

# 4. Enrichment CD8 TEM in females (Fig 4D) ----------

all_genes_F <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/AllCells_F_deaTopTable.rds" ))
all_genes_F <- split(all_genes_F$gene, all_genes_F$celltype)
all_genes_M <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/AllCells_M_deaTopTable.rds" ))
all_genes_M <- split(all_genes_M$gene, all_genes_M$celltype)

library(rrvgo)
enrichment_GO<- function(gene_list, universe, celltype=NULL, dir=NULL){
  print(celltype)
  GO_BP_result <- enrichGO(  gene_list,universe = universe, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont= "BP", pAdjustMethod = "fdr")
  # GO_BP_result <- GO_BP_result@result
  # GO_BP_result$celltype <- celltype
  # GO_BP_result$direction <- dir
  go_cd8tem <- GO_BP_result@result
  go_cd8tem <- go_cd8tem[go_cd8tem$p.adjust <0.05,]
  go_cd8tem <- go_cd8tem[order(go_cd8tem$p.adjust, decreasing = F), ]
  
  simMatrix <- calculateSimMatrix(go_cd8tem$ID,
                                  orgdb="org.Hs.eg.db",
                                  ont="BP",
                                  method="Rel")
  
  scores <- setNames(-log10(go_cd8tem$qvalue), go_cd8tem$ID)
  go_cd8tem_red <- reduceSimMatrix(simMatrix,
                                   scores,
                                   threshold=0.999,
                                   orgdb="org.Hs.eg.db")
  
  go_cd8tem_red_count <- go_cd8tem_red  %>% group_by(parentTerm) %>%
    tally() %>%                          # Count occurrences
    mutate(percentage = (n / sum(n)) * 100)
  
  return(go_cd8tem_red_count)
  
}

go_cd8tem_red_count <- enrichment_GO(only_F[only_F$celltype == "CD8 TEM" & only_F$direction == "up",]$gene, universe = unique(all_genes_F$gene))
options(enrichplot.colours = rev(c(alpha(blue, 0.3),blue)))
clusterProfiler::dotplot(CD8_TEM, showCategory=8)+theme+theme(legend.key.size = unit(0.8, "lines"),  legend.text = element_text(size = 11), axis.text.y = element_text(lineheight=.7, size=11))+
  scale_y_discrete(labels = label_wrap(35))+ggtitle("CD8 TEM- Females")


Fig5D <- ggplot(go_cd8tem_red_count, aes(y = reorder(parentTerm, n), x = percentage)) +
  geom_point(  aes(size=n),color=alpha(blue, 0.9)) +
  ylab("") +theme+scale_fill_gradient(high = alpha(blue, 0.3), low = blue)+
  xlab("Percentage of terms") +scale_y_discrete(labels = label_wrap(40))+scale_size_continuous(name="# terms")+
  theme(axis.text.y = element_text(size = 10), aspect.ratio=2,plot.title = element_text( face = "bold") )+ggtitle("CD8 TEM - Females",    subtitle = "Up-regulated DEGs")

# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig5/Fig5D_EnrichmentsTEM.pdf"), height =3.64, width= 5.81   )
# Fig5D
# dev.off()

#NK females --- 
go_nk_red_count <- enrichment_GO(only_F[only_F$celltype == "NK" & only_F$direction == "up",]$gene, universe = unique(all_genes_F$gene))
FigS5N <- ggplot(go_nk_red_count, aes(y = reorder(parentTerm, n), x = percentage)) +
  geom_point(  aes(size=n),color=alpha(blue, 0.9)) +
  ylab("") +theme+scale_fill_gradient(high = alpha(blue, 0.3), low = blue)+
  xlab("Percentage of terms") +scale_y_discrete(labels = label_wrap(40))+scale_size_continuous(name="# terms")+
  theme(axis.text.y = element_text(size = 12), aspect.ratio=2,plot.title = element_text( face = "bold") )+ggtitle("NK - Females",    subtitle = "Up-regulated DEGs")


# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS5/FigS5N_EnrichmentsNK.pdf"), height =3.41, width= 5.27   )
# FigS5N
# dev.off()

# B memory males ---

go_b_red_count <- enrichment_GO(only_M[only_M$celltype == "B memory" & only_M$direction == "down",]$gene, universe = unique(all_genes_F$gene))
FigS5O <- ggplot(go_b_red_count, aes(y = reorder(parentTerm, n), x = percentage)) +
  geom_point(  aes(size=n),color=alpha(red, 0.9)) +
  ylab("") +theme+scale_fill_gradient(high = alpha(red, 0.3), low = red)+
  xlab("Percentage of terms") +scale_y_discrete(labels = label_wrap(40))+scale_size_continuous(name="# terms")+
  theme(axis.text.y = element_text(size = 11), aspect.ratio=2,plot.title = element_text( face = "bold") )+ggtitle("B memory - down",    subtitle = "Up-regulated DEGs")

# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS5/FigS5O_EnrichmentsBmem.pdf"), height =3.41, width= 5.27   )
# FigS5O
# dev.off()



# 5. Disease enrichment (Fig 4E) ---------------
all_genes_F <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/AllCells_F_deaTopTable.rds" ))
all_genes_F <- split(all_genes_F$gene, all_genes_F$celltype)
all_genes_M <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/AllCells_M_deaTopTable.rds" ))
all_genes_M <- split(all_genes_M$gene, all_genes_M$celltype)

library(DOSE);library(clusterProfiler);  library("AnnotationDbi");library(org.Hs.eg.db)

enrichment_diseases <- function(gene_list, universe, celltype, direction) {
  print(celltype)
  
  # Initialize empty data frames to store results
  DO_df <- data.frame()
  disgenet_df <- data.frame()
  
  if (length(gene_list) != 0){
    # Try executing each operation separately
    genes <- mapIds(org.Hs.eg.db, gene_list, 'ENTREZID', 'SYMBOL')
    universe_list <- mapIds(org.Hs.eg.db, universe, 'ENTREZID', 'SYMBOL')
    
    DO <- enrichDO(genes, ont="DO", pAdjustMethod = "fdr", universe = universe_list, minGSSize = 10, maxGSSize = 1000)
    #DO <- enrichDGN(genes, ont="DO", pAdjustMethod = "fdr", universe = universe_list, minGSSize = 10, maxGSSize = 500)
    
    if(!is.null(DO)){
      DO_df <- DO@result
      DO_df$direction <- direction 
      DO_df$celltype <- celltype 
      
      return(list(DO_df))
    }
  }
  
}
get_df_disease_erichments <- function(cell_types, sex, dir, all=T){
  print(sex)
  if(sex=="F"){
    DO_F_up <- lapply(cell_types, function(celltype)enrichment_diseases(only_F_interaction[only_F_interaction$celltype == celltype & only_F_interaction$direction == dir,]$gene,all_genes_F[[celltype]], celltype, dir)[[1]])
    DO_F_up_df <- do.call(rbind.data.frame, DO_F_up)
    DO_F_up_df <- DO_F_up_df[DO_F_up_df$p.adjust < 0.05, ]
    DO_F_up_df$sex <- sex
  }else{
    DO_F_up <- lapply(cell_types[cell_types != "Platelet"], function(celltype)enrichment_diseases(only_M[only_M$celltype == celltype & only_M$direction == dir,]$gene, all_genes_M[[celltype]], celltype, dir)[[1]])
    DO_F_up_df <- do.call(rbind.data.frame, DO_F_up)
    DO_F_up_df$sex <- sex
    DO_F_up_df <- DO_F_up_df[DO_F_up_df$p.adjust < 0.05, ]
  }
  
  return(DO_F_up_df)
}

# only sex-specific DEGs
DO_F_up_df <- get_df_disease_erichments(unique(only_F$celltype)[!unique(only_F$celltype) %in% c("CD16 Mono", "CD4 TEM", "CD4 TCM", "Treg", "CD4 CTL", "HSPC")],"F", "up" ) %>% drop_na()
#DO_F_down_df <- get_df_disease_erichments(unique(only_F$celltype)[!unique(only_F$celltype) %in% c("CD16 Mono", "CD4 CTL", "Treg", "Platelet"    ,"gdT",  "HSPC", "MAIT")],"F", "down" ) #no enrichments

#DO_M_up_df <- get_df_disease_erichments(unique(only_M$celltype),"M", "up" ) #no enrichments
DO_M_down_df <- get_df_disease_erichments(unique(only_M$celltype)[!unique(only_M$celltype) %in% c("CD14 Mono", "CD4 CTL", "B intermediate", "Platelet", "HSPC", "CD16 Mono")],"M", "down" )

DO_all <- rbind(DO_F_up_df, DO_M_down_df)
DO_all$sex <- gsub("F", "Females",DO_all$sex )
DO_all$sex <- gsub("M", "Males",DO_all$sex )
DO_all$direction <- gsub("up", "up-regulated",DO_all$direction )
DO_all$direction <- gsub("down", "down-regulated",DO_all$direction )
DO_all$direction  <- factor(DO_all$direction , levels=c("up-regulated", "down-regulated"))
Fig5H <- ggplot(DO_all[DO_all$p.adjust < 0.05,], aes(x=celltype, y=Description))+geom_point(aes(size=Count, alpha=p.adjust, color=sex))+
  theme+ylab("")+xlab("")+scale_alpha_continuous(range = c(1, 0.5))+scale_y_discrete(labels = label_wrap(40))+scale_color_manual(values=c("Females"= blue, "Males"= red ))+
  scale_y_discrete(limits=rev)+theme( axis.text.y=element_text(size=11), plot.title = element_text( face = "bold",))+facet_grid(~direction, scales="free_x")




# &. GWAS overlap (Fig 4F) ---------------
gwas_allgenes <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/15_OverlapGWAS/FisherResults/FisherResults_AllDEGs.rds"))
gwas_signif <- gwas_allgenes[gwas_allgenes$AdjustedPValue < 0.05, ]
gwas_c <- gwas_signif %>% group_by(sex) %>% dplyr::count()
gwas_signif$celltype <- gsub("_", " ",gwas_signif$CellType )
Fig5I <-ggplot(gwas_signif, aes(x=celltype, y=Trait))+geom_point(aes(size=OR, alpha=AdjustedPValue), color=blue)+theme+facet_wrap(~sex, scales="free_x")+xlab("")+ylab("GWAS trait")+scale_y_discrete(labels = label_wrap(30))+
  theme+scale_color_manual(values=c("Males"=red, "Females"=blue))+scale_alpha_continuous(range = c(1, 0.5), name = "FDR")+theme(strip.text=element_text(size=13))


Fig5HF <- Fig5H + Fig5I+plot_layout(widths = c(3, 1.5))
# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig5/Fig5HF_DiseaseEnrichments.pdf"), height =3.35, width= 10.45   )
# Fig5HF
# dev.off()
gwas_SexSpec<- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/15_OverlapGWAS/FisherResults/FisherResults_SexSpecificDEGs.rds"))


#save supplementary table 14

sheets <- list("DEA_females" = deg_F, "DEA_males"=deg_M, "CODA_Females"=coda_F, "CODA_Males"= coda_M, "DEA_interaction"= interaction, "CODA_interaction"= coda, "CD8TEM_F_Enrichment"=go_cd8tem_red_count, "NK_F_Enrichment"=go_nk_red_count,"Bmem_M_Enrichment"=go_b_red_count,  "Disease_Enrichment"= DO_all, "GWAS_overlap_allDEGs"= gwas_allgenes, "GWAS_overlap_SexSpec_DEGs"= gwas_SexSpec)
write.xlsx(sheets, paste0(data_path, '/msopena/02_OneK1K_Age/supplementary_tables/TableS4_Sex_Stratified_DEA_CODA.xlsx'))






# 1. Heatmap expression by sex (Fig S4A, B) ------------

# Read in DEG results 
all_genes <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/AllCells_", sex, "_deaTopTable.rds"))
all_genes$celltype <- factor(all_genes$celltype, levels=rev(order_cells$cell_type))
bias <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/biomial_test_bias.rds"))
rownames(bias) <- bias$celltype
bias[bias$celltype == "CD8 Naive", ]$binom.test.p  <- .Machine$double.xmin
bias$log10_biom.p <- -log10(bias$binom.test.p)
bias$log10_biom.p_sign <- NA
bias[bias$direction == "up",]$log10_biom.p_sign <- bias[bias$direction == "up",]$log10_biom.p
bias[bias$direction == "down",]$log10_biom.p_sign <- -(bias[bias$direction == "down",]$log10_biom.p)
annot_bias <- bias[,c("celltype", "log10_biom.p_sign")]


plot_heatmap_sex <- function(degs){
  sharing <- degs %>% group_by(gene) %>% count()
  
  cells <- unique(degs[!degs$celltype %in% c("NK Proliferating", "gdT", "NK_CD56bright", "Platelet", "CD4 CTL", "HSPC"),]$celltype)
  # Select cell types to plot 
  df <- all_genes[all_genes$celltype %in% cells,] %>% as.data.frame()
  df <- reorder_cells(df)
  df <- df[!df$celltype %in% c("CD4 CTL"),]
  
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
  


  # Extract non concordant genes DE in more than 6 cell types 
  logfc_na_shar <- logfc[rownames(logfc) %in% sharing[sharing$n > 5,]$gene,]
  #if (nonconcordant==T) {logfc_na_shar <- logfc_na_shar[rownames(logfc_na_shar) %in% non_concordant$gene,]}
  
  #Same for pvalues 
  pvals <- as.matrix(significance[,-1])
  pvals <- pvals[rownames(pvals) %in% rownames(logfc_na_shar),]
  #if (nonconcordant==T) {pvals <- pvals[rownames(pvals) %in% non_concordant$gene,]}
  
  
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
  # get gene clusters and change names 
  rw <- heatmap$tree_col
  clusters_row <- cutree(rw, 2) %>% as.data.frame() %>% dplyr::rename("cluster_genes"=".")
  clusters_row$cluster_genes <- as.factor(clusters_row$cluster_genes)
  clusters_row$cluster_genes <- gsub("1", "Ribosomal_function",clusters_row$cluster_genes)
  clusters_row$cluster_genes <- gsub("2", "Immunological_process",clusters_row$cluster_genes)
  
  # get cell clusters and change names 
  cl <- heatmap$tree_row
  clusters_col <- cutree(cl, 2) %>% as.data.frame() %>% dplyr::rename("bias"=".")
  clusters_col$bias <- "upregulation_bias"
  clusters_col[c("CD4 Naive", "CD8 Naive", "B memory", "B naive", "B intermediate", "MAIT"),] <- "downregulation_bias"
  clusters_col$bias <- as.factor(clusters_col$bias)
  clusters_col$bias <- gsub("2", "downregulation_bias",clusters_col$bias)
  clusters_col$bias <- gsub("1", "upregulation_bias",clusters_col$bias)
  annot_bias <- annot_bias[rownames(clusters_col),]
  annot_bias_cont <- annot_bias[,-1] %>% as.data.frame()%>% dplyr::rename("Binomial_p.val"=".")
  rownames(annot_bias_cont) <- rownames(annot_bias)
  
  
  #colors 
  mat_breaks <- seq(-0.01, 0.01, by = 0.001)
  min_value <- min(annot_bias_cont, na.rm = TRUE)
  max_value <- max(annot_bias_cont, na.rm = TRUE)
  
  # Create a color function centered at 0 using colorRamp2 from the circlize package
  breaks <- seq(min_value, max_value, length.out = 10)
  continuous_colors_func <- colorRamp2(c(min_value, 0, max_value), c(blue, "white", red))
  
  # Generate the color mapping for annotation_row
  annotation_colors_for_row <- continuous_colors_func(seq(min_value, max_value, length.out = 10))
  
  annot_colors <-list(cluster_genes= c("Ribosomal_function"="lightgrey", "Immunological_process"="#777777"),bias=c("downregulation_bias"=alpha("#264653", 0.8), "upregulation_bias"="#bc4749"))
  # annot_colors <-list(Genes_clust= c("Ribosomal_function"="lightgrey", "Immunological_process"="#777777"),Binomial_p.val=annotation_colors_for_row)
  
  # Heatmap
  heatmap_horizontal <- pheatmap::pheatmap(t(logfc_na_shar), display_numbers = t(pvals),
                                           method = "complete",
                                          cluster_rows = sort_hclust(hclust(dist(t(logfc_na_shar)), method = "complete")),
                                           annotation_col = clusters_row,
                                           annotation_row = clusters_col,
                                           annotation_names_row = FALSE,
                                           annotation_names_col = FALSE,
                                           treeheight_col = 10,treeheight_row = 10,
                                           fontsize_number = 12, fontsize_col = 6, fontsize_row = 10,
                                           breaks = mat_breaks, border_color = NA,
                                           color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(20),
                                           annotation_colors = annot_colors
  )
  
  return(list(heatmap_horizontal))
}
#Males 
degs_M <- get_significant_degs_sex("M", F)
figS5A <- plot_heatmap_sex(degs_M)[[1]]

degs_F<- get_significant_degs_sex("F", F) %>% drop_na()
figS5B <- plot_heatmap_sex(degs_F)[[1]]



# 2. Compare heatmaps with non-stratified analysis (Fig S4C, D) ----------

get_info_dendrogram <- function(sex){
all_genes <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/AllCells_", sex, "_deaTopTable.rds" ))
all_genes$celltype <- factor(all_genes$celltype, levels=rev(order_cells$cell_type))
all_genes$direction <- ifelse(all_genes$logFC > 0, "up", "down")
# concordance <- lapply(unique(all_genes$gene), function(gene) ifelse(length(unique(all_genes[all_genes$gene ==gene,]$direction)) == 1, "concordant", "non_concordant"))
# names(concordance) <- unique(all_genes$gene)
# concordance_df <- data.frame(gene = names(concordance), concordance = unlist(concordance), row.names = NULL)
concordance_df <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/Concordance_", sex, "_unbalanced.rds" ))

# Get siginificant DEGs
degs <-get_significant_degs_sex(sex, F) %>% tidyr::drop_na()

# sharing 
sharing <- degs %>% group_by(gene) %>% dplyr::count()
concordance_df <- as.data.frame(concordance_df[concordance_df$gene %in% sharing$gene, ])
non_concordant <- concordance_df[concordance_df$concordance == "non_concordant",]

cells <- unique(degs[!degs$celltype %in% c("NK Proliferating", "gdT", "NK_CD56bright", "Platelet"),]$celltype)
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
return(list(sharing, non_concordant, logfc))
}

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
  d <- ggplot(data = labels,  aes(x = x, y = y)) +geom_segment(data = segment(dhc), aes(x = x, y = y, xend = xend, yend = yend)) +geom_text(aes(label = label, hjust=0),  size = 4) +
    geom_point( size = 3, aes(color=celltype_l1)) +
    coord_flip() +scale_y_reverse(expand = c(0.5, 0))+theme+theme(axis.text= element_blank(), axis.line=element_blank(), panel.border=element_blank(), axis.ticks=element_blank())+
    ylab("")+xlab("")+ggtitle(paste0("n=", n))  +theme(plot.title = element_text(hjust =0.5), legend.position="none")+scale_color_manual(values= c("CD8 T"= "#679267", "CD4 T"= "#00563f",  "other T"= "#c0d3ac", "NK"= "#E76F51", "B"= "#8580C2", "Mono" ="#d88d8b" ),)
  return(list(clust, d))
}



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
  if(sex=="F"){n9 <- as.dendrogram(plot_dendrogram(9, method , nonconcordant)[[1]])
  dist_list <- dendlist("DEGs_2_cells" =n2, "DEGs_3_cells"= n3, "DEGs_4_cells"=n4, "DEGs_5_cells"=n5, "DEGs_6_cells"=n6, "DEGs_7_cells"=n7, "DEGs_8_cells"= n8, "DEGs_9_cells"=n9)}else{
    dist_list <- dendlist("DEGs_2_cells" =n2, "DEGs_3_cells"= n3, "DEGs_4_cells"=n4, "DEGs_5_cells"=n5, "DEGs_6_cells"=n6, "DEGs_7_cells"=n7, "DEGs_8_cells"= n8)
  }
  mat <- as.matrix(cor.dendlist(dist_list))
  mat[lower.tri(mat)] <- NA
  p <- ggcorrplot::ggcorrplot(mat,  method = "square", lab = T)+scale_fill_gradient2(low="white", high = "#A93F41",midpoint = 0.6 )+theme+
    theme(axis.text.x=element_text(angle=45, vjust = 0.5))+xlab("")+labs(x="", y="", fill="correlation")
  return(p)
}


#females
sharing <- get_info_dendrogram("F")[[1]]
non_concordant <- get_info_dendrogram("F")[[2]]
logfc <- get_info_dendrogram("F")[[3]]


figS5c <- compare_dendrograms_n("average", T)

#males
sharing <- get_info_dendrogram("M")[[1]]
non_concordant <- get_info_dendrogram("M")[[2]]
logfc <- get_info_dendrogram("M")[[3]]

figS5d <-compare_dendrograms_n("average", T)

# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS5/FigS5D_CorrHeatmaps_F.pdf"), width =5.59, height = 3.72 )
# figS5c
# dev.off()
# 
# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS5/FigS5D_CorrHeatmaps_M.pdf"), width =5.59, height = 3.72 )
# figS5d
# dev.off()
# 3. Overlap with all the data (Fig S4E ) --------

all_genes_F <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/AllCells_F_deaTopTable.rds" ))
all_genes_F <- split(all_genes_F$gene, all_genes_F$celltype)
all_genes_M <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/AllCells_M_deaTopTable.rds" ))
all_genes_M <- split(all_genes_M$gene, all_genes_M$celltype)
all_genes <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/DEG_Age_cell_type.rds"))
all_genes <- split(all_genes$gene, all_genes$celltype)

deg_all <- get_significant_degs("Age", "cell_type")
deg_M <- get_significant_degs_sex("M") %>% drop_na()
deg_F <- get_significant_degs_sex("F") %>% drop_na()
deg_M$overlap <-NA
deg_M[deg_M$direction == "up",]$overlap <- unlist(lapply(unique(deg_M[deg_M$direction == "up",]$cell_type), function(cell) 
  ifelse(deg_M[deg_M$direction == "up" & deg_M$cell_type == cell,]$gene %in% deg_all[deg_all$direction == "up" & deg_all$celltype == cell,]$gene, 
         "overlap", "not_overlap")))
deg_M[deg_M$direction == "down",]$overlap <- unlist(lapply(unique(deg_M[deg_M$direction == "down",]$cell_type), function(cell) 
  ifelse(deg_M[deg_M$direction == "down" & deg_M$cell_type == cell,]$gene %in% deg_all[deg_all$direction == "down" & deg_all$celltype == cell,]$gene, 
         "overlap", "not_overlap")))
deg_F$overlap <-NA
deg_F[deg_F$direction == "up",]$overlap <- unlist(lapply(unique(deg_F[deg_F$direction == "up",]$cell_type), function(cell) 
  ifelse(deg_F[deg_F$direction == "up" & deg_F$cell_type == cell,]$gene %in% deg_all[deg_all$direction == "up" & deg_all$celltype == cell,]$gene, 
         "overlap", "not_overlap")))
deg_F[deg_F$direction == "down",]$overlap <- unlist(lapply(unique(deg_F[deg_F$direction == "down",]$cell_type), function(cell) 
  ifelse(deg_F[deg_F$direction == "down" & deg_F$cell_type == cell,]$gene %in% deg_all[deg_all$direction == "down" & deg_all$celltype == cell,]$gene, 
         "overlap", "not_overlap")))
deg_sex <- rbind(deg_M, deg_F)
deg_sex <- deg_sex[!deg_sex$celltype %in% c("Platelet", "HSPC"),]

df_M <- deg_sex %>% dplyr::group_by( overlap, celltype, sex) %>% dplyr::count()
df_M <- reorder_cells(df_M, neworder = T, reverse = T)
df_M$sex <- gsub("M", "Males", df_M$sex)
df_M$sex <- gsub("F", "Females", df_M$sex)
FIGS5H <- ggplot(df_M, aes(x=celltype, y=n, alpha=overlap, fill=sex)) +
  geom_col( width=0.75) + scale_alpha_manual(values=c(0.5, 1))+
  coord_flip() + theme   + # Add this line
  ylab("Number of DEGs") + xlab("") +scale_y_continuous()+ scale_y_continuous(labels = abs, limits = c(0, 2100), breaks = c(0, 1000, 2000)) +facet_grid(celltype_l1~sex, scales="free", space="free")+
  theme(axis.text=element_text(size=11), strip.text=element_text(size=11), legend.position="top", axis.title=element_text(size=12), legend.text=element_text(size=11), legend.title=element_text(size=11))+ 
  guides(alpha = guide_legend(title = "overlap all data")) + scale_fill_manual(values=c("Females"=blue, "Males"=red), guide = "none")

# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS5/FigS5G_OverlapAllData.pdf"), width =4.22, height = 5.66 )
# FIGS5H
# dev.off()


# 4. Overlap between males and females (Fig S4F) ------
all_degs <- rbind(deg_F, deg_M)
all_degs <- reorder_cells(all_degs, neworder = T, reverse = T)

all_degs$concordance <- NA
for(cell in unique( all_degs[all_degs$overlap == "overlap",]$celltype)){
  all_degs[all_degs$overlap == "overlap" & all_degs$celltype == cell,]$concordance <- unlist(  all_degs[all_degs$overlap == "overlap" & all_degs$sex == "M"  & all_degs$celltype == cell ,]$direction ==  
                                                                                     all_degs[all_degs$overlap == "overlap" & all_degs$sex == "F" & all_degs$celltype == cell,]$direction)
}


all_degs[all_degs$overlap == "overlap",]$overlap <- ifelse(all_degs[all_degs$overlap == "overlap",]$concordance == T , "overlap_concordant", "overlap_non_concordant")

FIGS5I<-ggplot(all_degs[!all_degs$celltype %in% c("Platelet", "HSPC"),], aes(x = celltype, fill = overlap)) +
  geom_bar(position = "stack", stat = "count") +
  scale_fill_manual(values = c("only_Females" =blue, "only_Males" = red, "overlap_concordant" = "#777","overlap_non_concordant"="grey" )) +
  facet_grid(celltype_l1 ~ ., scales = "free", space = "free") +
  xlab("") +theme+  theme(legend.title=element_blank(), legend.position="top",legend.key.size=) +  # Ensure your theme function is properly called
  coord_flip() + scale_y_continuous(labels = abs, limits = c(0, 4000), breaks = c(0, 2000, 4000)) +
  ylab("Frequency of genes")+   guides(  color=guide_legend(override.aes = list(size=1)))


# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS5/FigS5H_OverlapM_F.pdf"), width =3.34, height = 5.66 )
# FIGS5I
# dev.off()


# 5. Correlate cellular proportions with age between M and F and with all the data (Fig S4G) -------
coda_all <- readRDS(paste0(data_path, "aripol1/wijst-2020-hg19/v1/aging/25.cell_type_proportions.OneK1K_age_corrected/cell_type/CLR/phenotype_stats.rds"))$Age
coda_all$celltype <- gsub("_", " ", coda_all$celltype )
colnames(coda_all) <- gsub("estimate", "estimate_all", colnames(coda_all))
coda_M$estimate_M <- coda_M$estimate
coda_F$estimate_F <- coda_F$estimate

cM <- merge(coda_all, coda_M, by="celltype")
cF <- merge(coda_all, coda_F, by="celltype")
c <- merge(cM, cF, by="celltype")

c_sex_vs_all <- rbind(cF[, -which(names(cF) == "estimate_F")], cM[, -which(names(cM) == "estimate_M")])
c_sex_vs_all$significance <- ifelse(c_sex_vs_all$fdr.y < 0.05, "ss", "ns")
c_sex_vs_all <-reorder_cells(c_sex_vs_all)
c_sex_vs_all <- c_sex_vs_all[c_sex_vs_all$celltype_l1 != "other",]
c_sex_vs_all <- c_sex_vs_all[c_sex_vs_all$celltype %in% bias_estimate_M$celltype,]
figS5e <- ggplot(c_sex_vs_all, aes(x=estimate_all, y=estimate))+geom_point(aes(color=sex, alpha=significance))+theme+  geom_smooth(method = "lm", se = FALSE, aes(color=sex))+theme+
  ylab("CoDA estimate sex-stratified")+xlab("CoDA estimate all")+  stat_cor(method = "spearman",  p.digits = 1 ,aes(color = sex))+
  scale_alpha_manual(values=c("ss"=1, "ns"=0.5))+scale_color_manual(values=c(blue, red))+theme(aspect.ratio=1)#+geom_text_repel(aes( label=celltype), position="identity")
# 
# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS5/FigS5E_Correlation_CoDA_all.pdf"), width =5.59, height = 3.72 )
# figS5e
# dev.off()





# 6. Binomial FDR in males and in females (Fig S4H) -------
bias <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/biomial_test_bias_sex.rds"))
bias_F <- bias[bias$sex == "F",]
bias_F <- bias_F[!duplicated(bias_F$celltype),]
bias_M <- bias[bias$sex == "M",]
bias_M <- bias_M[!duplicated(bias_M$celltype),]

bias_cor <- merge(bias_M, bias_F, by="celltype")
bias_cor[bias_cor$celltype == "CD8 Naive", ]$binom.test.p.x <- .Machine$double.xmin
bias_cor[bias_cor$celltype == "CD8 Naive", ]$binom.test.p.y <- .Machine$double.xmin
bias_cor <- reorder_cells(bias_cor, neworder = T)
FigS5F <- ggplot(bias_cor, aes(y=-log10(binom.test.fdr.x), x=-log10(binom.test.fdr.y)))+geom_point(size=2.5, color=blue)+theme+ylab("-log10(fdr binomial males)")+xlab("-log10(fdr binomial females)")+
  theme(aspect.ratio=1, legend.position="right")+scale_x_continuous()+
  geom_text_repel( data=bias_cor, position = "identity",aes( label=celltype))

# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS5/FigS5F_Correlation_Binomials.pdf"), width =5.59, height = 3.72 )
# FigS5F
# dev.off()



# 7. Correlate coda estimates M and F (Fig S4I) ---------

coda_M <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/12_CellProportions_Sex/CODA_resultsAge_M.rds"))
coda_M$sex <- "Males"
coda_F <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/12_CellProportions_Sex/CODA_resultsAge_F.rds"))
coda_F$sex <- "Females"

coda_estimate <- merge(coda_F[,c("sex", "celltype", "estimate", "fdr", "p.value")], coda_M[,c("sex", "celltype", "estimate", "fdr", "p.value")], by="celltype")
coda_estimate <- coda_estimate[!coda_estimate$celltype %in% c("CD4 CTL", "NK Proliferating", "CD4 Proliferating"),]
coda_estimate <- reorder_cells(coda_estimate)
coda_estimate <- coda_estimate[!coda_estimate$celltype_l1 %in% c("other", "DC"),]

FigS5F <- ggplot(coda_estimate, aes(y=-log10(fdr.x), x=-log10(fdr.y)))+geom_point(size=2.5, color=blue)+theme+ylab("-log10(p.value CoDA females)")+xlab("-log10(p.value CoDA males)")+
  theme(aspect.ratio=1, legend.position="right")+scale_color_manual(values = palette_allcells)+scale_x_continuous()+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+geom_vline(xintercept = -log10(0.05), linetype="dashed")+
  geom_text_repel( data=coda_estimate, position = "identity",aes( label=celltype))

# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS5/FigS5H_Correlation_CoDA.pdf"), width =5.59, height = 3.72 )
# FigS5F
# dev.off()



# 8. Overlap with the interaction ------------

interaction$classification  <- factor(interaction$classification , levels=c("no_DE", "DE_both", "DE_M", "DE_F"))
interaction <- interaction[interaction$cell_type != "Plasmablast", ]
FIGS5J <- ggplot(interaction, aes(x = celltype, fill = classification)) +
  geom_bar(stat = "count") +
  facet_grid(celltype_l1 ~ ., scales = "free", space = "free") +
  theme+
  theme(
    axis.text.x = element_text(size = 13),
    legend.position = "top",
    legend.title = element_blank(),
    legend.margin = margin(r = 70, l = 5, t = 5, b = 2),
    legend.key.size = unit(0.5, "cm"),
  ) +
  scale_fill_manual(
    values = c(
      "DE_F" = blue, 
      "DE_M" = red, 
      "no_DE" = "darkgrey", 
      "DE_both" = "lightgrey"
    )
  ) + coord_flip() + xlab("") + ylab("Number of age-DEGs")


# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS5/FigS5H_OverlapInteraction.pdf"), width =3.34, height = 5.66 )
# FIGS5J
# dev.off()


coda<- readRDS(paste0(data_path, "/aripol1/wijst-2020-hg19/v1/aging/25.cell_type_proportions.OneK1K_age_corrected/cell_type/CLR/interaction/phenotype_stats.rds"))
coda <- coda$`Age:SexF`
coda$celltype <- gsub("_", " ", coda$celltype)
coda$celltype <- gsub(" CD56bright", "_CD56bright", coda$celltype)
coda <- coda[coda$celltype %in% cells_to_keep,]
#recompute FDR with only the cells we selected 
coda$fdr <- p.adjust(coda$p.value, method = "BH")
coda$significance <- ifelse(coda$p.value< 0.01, "ss","ns")
coda <- reorder_cells(coda, reverse = T, neworder = T)
#coda <- coda[!coda$celltype_l1 %in% c("DC", "other"), ]
coda$direction <- ifelse(coda$estimate < 0, "down", "up")


FIGS5K <- ggplot(coda, aes(x=estimate, y=celltype)) + 
  geom_point(aes(alpha=significance,color=direction), size=4) + xlab("CoDA estimate")+ylab(NULL)+
  geom_pointrange(aes(xmin=conf.low, xmax=conf.high, alpha=significance, color=direction), fatten = .1) +
  geom_vline(xintercept=0, linetype = "dashed")  + 
  theme + theme( axis.text.x = element_text(size = 12),axis.title=element_text(size=12), legend.position="top", legend.title=element_blank(),legend.key.size = unit(0.5, "lines") ) +
  scale_fill_manual(values = pal_direction)+scale_color_manual(values=pal_direction)+
  scale_alpha_manual(values=alpha_vec) +facet_grid(celltype_l1~., space="free", scale="free")+ scale_x_continuous(limits = c(-0.015, 0.015), breaks = seq(-0.01, 0.01, by = 0.01))

# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS5/FigS5K_CoDA_Interaction.pdf"), width =2.9, height = 5.66 )
# FIGS5K
# dev.off()






##### Supplemenatry figure S6 ####

#1. Number of enriched/depleted nhoods per cell type and sex ----
da_results_M <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DAA_results_cell_type_annotation__Sex_0.5_M.rds"))
da_results_M$sex <- "Males"
da_results_F <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DAA_results_cell_type_annotation__Sex_0.5_F.rds"))
da_results_F$sex <- "Females" 
da_results <- rbind(da_results_F, da_results_M)

da_results <- da_results[!is.na(da_results$SpatialFDR),]
da_results$celltype <- da_results$cell_type
da_results$celltype <- gsub(" CD56bright", "_CD56bright", da_results$celltype)

da_results <- da_results[da_results$celltype %in% cells_to_keep,]
da_results$direction <- ifelse(da_results$logFC < 0, "depleted", "enriched")
da_results$significance <- ifelse(da_results$SpatialFDR < 0.1, "ss", "ns")
da_results <- reorder_cells(da_results, reverse = T, neworder = T)
# da_results$direction <- ifelse(da_results$logFC < 0, "down", "up")
da_results[da_results$SpatialFDR > 0.05,]$direction <- NA
da_results <- da_results %>%
  arrange(!is.na(direction))


da_M <- da_results[da_results$sex =="Males",  ] 
da_count_M <-da_M %>% group_by(direction, celltype) %>%   summarise(count = n() ) %>% ungroup() %>% group_by(celltype) %>%
   mutate( prop = count / sum(count) )

da_count_M$sex <- "Males"
da_count_F <-  da_results[da_results$sex =="Females",  ] %>% group_by(direction, celltype) %>%   summarise(count = n() ) %>% ungroup() %>% group_by(celltype) %>%
  mutate( prop = count / sum(count) )
da_count_F$sex <- "Females"

da_count <- rbind(da_count_F, da_count_M)
da_count <- reorder_cells(da_count, neworder = T)
da_count <-da_count[!is.na(da_count$direction),]

p_perc_nhoods <- ggplot(da_count, aes(y=celltype, x=prop , fill=sex)) + geom_bar( stat="identity", position="dodge")+scale_fill_manual(values= c("Males"=red, "Females"=blue))+
  theme+facet_grid(celltype_l1~direction, scales="free_y", space="free")+xlab("Proportion of Nhoods")+scale_x_continuous(breaks = c(0, 0.5, 1))+theme(axis.text=element_text(size=10), axis.text.y = element_blank(), legend.position="top")+
geom_text(aes(label=count), position =  position_dodge(width = 0.9), size=3)+ylab("")



nnhoods <- da_results %>% group_by(celltype, sex) %>% count()
nnhoods <- reorder_cells(nnhoods, neworder = T)
p_nnhoods <- ggplot(nnhoods, aes(y=celltype, x=n)) + geom_bar(stat="identity", aes(fill=sex), position = "dodge", alpha=0.6, width = 0.8)+  scale_x_continuous(labels = function(x) format(x / 1e3, scientific = FALSE), limits = c(0,max(ncells$n)),breaks =c(0, 2000) ) +
  theme  + xlab("N Nhoods") + ylab("") +scale_fill_manual(values=c("Males"= red, "Females"= blue))+
  theme(axis.text.x=element_text(size=9), axis.text.y=element_text(size=11), axis.title=element_text(size=11),  legend.position="none",strip.text=element_blank())+ facet_grid(celltype_l1~ ., scales="free", space="free")


library(patchwork)
FigS6A <- p_nnhoods +p_perc_nhoods+plot_layout(width=c(1, 5))

#2. GZMB+ expression markers females ----
library(Seurat)
so_CD8TEM <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/01_Milo_NhoodMarkers/so_CD8TEM_upNhoods_Sex_",sex,".rds"))

ft_GZMH <- FeaturePlot(so_CD8TEM, features = "GZMH" , pt.size = 0.2, cols = c("lightgrey", red) )&theme&xlab("")&ylab("")&
  theme(aspect.ratio=1, axis.title=element_text(size=10),plot.title = element_text(size = 15, hjust = 0.5, face="bold.italic"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.position="none" )
ft_FGFBP2 <- FeaturePlot(so_CD8TEM, features = "FGFBP2" , pt.size = 0.2, cols = c("lightgrey", red) )&theme&xlab("")&ylab("")&
  theme(aspect.ratio=1, axis.title=element_text(size=10),plot.title = element_text(size = 15, hjust = 0.5, face="bold.italic"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.position="left" )
ft_GZMB <- FeaturePlot(so_CD8TEM, features = "GZMB" , pt.size = 0.2, cols = c("lightgrey", red) )&theme&xlab("")&ylab("")&
  theme(aspect.ratio=1,axis.title=element_text(size=10),plot.title = element_text(size = 15, hjust = 0.5, face="bold.italic"),
        panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.position="none" )

Fig4C_top <- ft_GZMH+ft_FGFBP2+ft_GZMB+ plot_layout(nrow = 3)

# pdf(paste0(data_path, "/msopena/02_OneK1K_Age/figures/FigS6/FigS6C_GeneMarkers_UMAP.pdf"), width = 5.5, height = 5.5)
# Fig4C_top
# dev.off()


#3.2. barplot 
expr <- so_CD8TEM@assays$RNA$counts
cell_nhood_df <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Nhood-cell_id_df_Sex_",sex,".rds"))
cells_to_keep <- unique(cell_nhood_df[cell_nhood_df$SpatialFDR < 0.05,]$cell_type_nhood)
cells_to_keep <- cells_to_keep[!cells_to_keep %in% c("CD4 CTL", "dnT", "HSPC")]
cell_nhood_df <- cell_nhood_df %>% tidyr::drop_na()
cell_nhood_df <- cell_nhood_df[cell_nhood_df$cell_type_nhood %in% c(cells_to_keep),]
cell_nhood_df$ss <- ifelse(cell_nhood_df$SpatialFDR < 0.05, "ss", "ns")
cell_nhood_df$direction <- ifelse(cell_nhood_df$logFC > 0, "enriched", "depleted")
cell_nhood_df$color <- "ns"
cell_nhood_df[cell_nhood_df$ss == "ss",  ]$color <- cell_nhood_df[cell_nhood_df$ss == "ss",  ]$direction

# extract expression matrix of this subsampling --- 

plot_example_nhood <- function( gene, cd8tem=F, cells=NULL){
  expr_gene <- expr[gene,] %>% as.data.frame()
  colnames(expr_gene) <- "Expression"
  expr_gene$cell_id <- rownames(expr_gene)  
  expr_gene <- expr_gene %>% left_join(cell_nhood_df, by="cell_id")
  expr_nhood <- expr_gene %>% tidyr::drop_na() %>% dplyr::rename(celltype =cell_type_nhood) %>% group_by(Nhood, celltype, color ) %>% 
    dplyr::summarize(across(Expression, mean, na.rm = TRUE)) 
  
  df <- reorder_cells(expr_nhood, neworder = T)
  if(cd8tem == F){
    return (ggplot(df, aes(x=celltype, y=Expression))+geom_boxplot(aes(fill=color), outlier.shape = NA)+theme+
              scale_fill_manual(values=c("enriched"=red, "depleted"=blue, "ns"="grey"), )+xlab("")+labs(fill="Nhoods")+
              ggtitle(paste0(gene))+theme(plot.title = element_text(size = 15, face = "bold.italic", hjust = 0.5), axis.text.x=element_text(size=10, angle=90,  hjust = 0.95, vjust = 0.6)))
  }else{
    ggplot(df[df$celltype %in% cells, ], aes(x=celltype, y=Expression))+geom_boxplot(aes(fill=color), outlier.shape = NA)+theme+
      scale_fill_manual(values=c("enriched"=red, "depleted"=blue, "ns"="grey"), )+xlab("")+labs(fill="Nhoods")+
      ggtitle(paste0(gene))+theme(plot.title = element_text(size = 15, face = "bold.italic", hjust = 0.5), axis.text.x=element_text(size=12))
  }}

GZMH_all <- plot_example_nhood("GZMH")+theme(axis.text.x=element_blank())
#GZMH <- plot_example_nhood("GZMH", cd8tem=T, cells=c("CD8 TEM", "NK"))

FGFBP2_all <- plot_example_nhood("FGFBP2")+theme(axis.text.x=element_blank())
#FGFBP2 <- plot_example_nhood("FGFBP2",  cd8tem=T, cells=c("CD8 TEM", "NK"))

GZMB_all <- plot_example_nhood("GZMB")
#GZMB <- plot_example_nhood("GZMB",  cd8tem=T, cells=c("CD8 TEM", "NK"))

Fig4C_down <- GZMH_all +plot_spacer()+FGFBP2_all+ylab("")+plot_spacer()+GZMB_all+ylab("")+plot_layout( ncol = 1)


# pdf(paste0(data_path, "/msopena/02_OneK1K_Age/figures/FigS6/FigS6C_GeneMarkers_Boxplot.pdf"), width = 4.24, height = 8.22)
# Fig4C_down
# dev.off()


# 3. B cell subpopulation markers ----



