

# Plots results of the Compositional Anlysis (Figure 3 and Supplemenatal Figure 3 )

library(ggplot2); library(dplyr); library(RColorBrewer);library(ggside);library(ggbeeswarm);library(ggpubr);library(ggrepel)


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



#data
metadata <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/metadata_processed.rds"))
order_cells<- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/order_cells.rds"))
cells_to_keep <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/cells_to_keep.rds"))
#order_cells<- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/new_order_cells.rds"))
celltype_l1 <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/celltypes_equivalence.rds"))
celltype_l1$cell_type <- factor(celltype_l1$cell_type, levels = order_cells$cell_type)
mdata_donors <- metadata[!duplicated(metadata$assignment),]
#saveRDS(mdata_donors, paste0(basepath, "Data/scRNAseq/Yazar2022/donor_metadata.rds"))


#themes and palettes 
alpha_vec <- c(0.4, 1)
computer <- "work"
source(paste0(data_path, "/msopena/01_meta-analysis_Age_Sex/scripts/themes.R"))

#1. Plot CoDA results -------

coda <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/coda_results.rds"))$Age
coda$celltype <- gsub("_", " ", coda$celltype)
deg <- get_significant_degs("Age", "cell_type")
coda <- coda[coda$celltype %in% unique(degs$celltype),]
coda$fdr <- p.adjust(coda$p.value, method = "fdr")
coda$celltype <- gsub(" CD56bright", "_CD56bright", coda$celltype)


plot_proportions <- function(cell_level){
  coda <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/coda_results.rds"))$Age
  coda$celltype <- gsub("_", " ", coda$celltype)
  coda$celltype <- gsub(" CD56bright", "_CD56bright", coda$celltype)
  coda <- coda[coda$celltype %in% unique(degs$celltype),]
  #recompute FDR with only the cells we selected 
  coda$fdr <- p.adjust(coda$p.value, method = "fdr")
  coda$significance <- ifelse(coda$fdr< 0.05, "ss","ns")
  md <- metadata%>% tidyr::drop_na(paste0(cell_level))
  count_cells <- md %>% dplyr::group_by(md[[cell_level]]) %>% dplyr::count() %>% as.data.frame()
  count_cells <- count_cells[match(coda$celltype,count_cells[,1]),]
  coda$n <- count_cells$n
  coda <- reorder_cells(coda, reverse = T, neworder = T)
  #coda <- coda[!coda$celltype_l1 %in% c("DC", "other"), ]
  coda$direction <- ifelse(coda$estimate < 0, "down", "up")
  
   
  p <- ggplot(coda, aes(x=estimate, y=celltype)) + 
    geom_point(aes(alpha=significance,color=direction), size=4) + xlab("CoDA estimate")+ylab(NULL)+
    geom_pointrange(aes(xmin=conf.low, xmax=conf.high, alpha=significance, color=direction), fatten = .1) +
    geom_vline(xintercept=0, linetype = "dashed")  + 
    theme + theme(axis.text.y = element_blank(), axis.text.x = element_text(size = 12),axis.title=element_text(size=12), legend.position="top", legend.title=element_blank(),legend.key.size = unit(0.5, "lines") ) +
    scale_fill_manual(values = pal_direction)+scale_color_manual(values=pal_direction)+
    scale_alpha_manual(values=alpha_vec) +facet_grid(celltype_l1~., space="free", scale="free")+ scale_x_continuous(limits = c(-0.015, 0.015), breaks = seq(-0.01, 0.01, by = 0.01))
  # + geom_ysidecol(aes(x=n), color="#4347534C")  +
  # scale_ysidex_continuous(guide = guide_axis(angle = 0), minor_breaks = NULL, n.breaks = 2)+ ggside( y.pos = "right")
  return(p)
}


coda_plot <- plot_proportions("cell_type") + theme( axis.text = element_text(size = 13))


# Extract the numner of cells 
ncells <- metadata %>% dplyr::group_by(cell_type) %>% dplyr::count() %>% tidyr::drop_na()
colnames(ncells)[1] <- "celltype"
ncells <- ncells[ncells$celltype %in% unique(degs$celltype),] %>% as.data.frame()
ncells <- reorder_cells(ncells, reverse = T, neworder = T)


p_ncells <- ggplot(ncells, aes(y=celltype, x=n)) + geom_bar(stat="identity", fill=alpha(blue, 0.4))+  scale_x_continuous(labels = function(x) format(x / 1e6, scientific = FALSE)) +
  theme  + xlab("N cells (million)") + ylab("") +
  theme(axis.text.x=element_text(size=9), axis.text.y=element_text(size=11), axis.title=element_text(size=11),  strip.text=element_blank())+ facet_grid(celltype_l1~ ., scales="free", space="free")

library(patchwork)
fig2A <- p_ncells+coda_plot+plot_layout(widths = c(5, 15))
  
# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig3/Fig3A_CODA_allcells.pdf"), height =6.68, width= 4.49 )
# fig2A
# dev.off()

#2. Plot correlation coda and expression bias -------

coda <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/coda_results.rds"))$Age
biomial_pval <-readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/biomial_test_bias_02.rds"))

#coda <- readRDS(paste0(data_path, "aripol1/wijst-2020-hg19/v1/aging/25.cell_type_proportions.OneK1K_age_corrected/cell_type/CLR/coda_results.rds"))$Age
coda$celltype <- gsub("_", " ", coda$celltype)
coda$direction <- ifelse(coda$estimate < 0, "down", "up")
degs <- get_significant_degs("Age", "cell_type")
degs <- degs[!degs$celltype %in% c("gdT", "NK Proliferating", "Platelet"),]


biomial_pval$binom.test.fdr <- p.adjust(biomial_pval$binom.test, method = "fdr")
bias_estimate <- merge(coda, biomial_pval, by="celltype")

bias_estimate$significance <- ifelse(bias_estimate$fdr <0.05,"ss", "ns")
bias_estimate[bias_estimate$celltype == "NK",]$prop_up <- 0.89
#bias_estimate <- bias_estimate[bias_estimate$p.val < 0.05,]

bias_estimate$p.val_dir <- ifelse(bias_estimate$prop_up > 0.5, bias_estimate$binom.test.p , -bias_estimate$binom.test.p )


fig2b <- ggplot(bias_estimate, aes(y=prop_up, x=estimate))+geom_point( color=blue, aes(alpha=significance, size=-log10( binom.test.fdr)))+theme+  
  geom_smooth(method = "lm", se = FALSE, color=blue)+ylab("Proportion up-regulated DEGs")+  scale_size_continuous(range = c(0.5, 3))+
  geom_text_repel( data=bias_estimate, aes(x=estimate, y=prop_up, label=celltype), size=3.1)+geom_vline(xintercept = 0, color="#777777" , linetype="dashed", size=0.5)+
  scale_alpha_manual(values=c("ss"=1, "ns"=0.5))+xlab("CoDA estimate")+theme(legend.position="right", axis.title=element_text(size=11), legend.title=element_text(size=10))+
  stat_cor(method = "spearman", label.y = 1, label.x = -0.0085,   p.digits = 1, size=3.4)+labs(alpha= "CoDA singnificance", size="-log10(binom.FDR)")


# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig3/Fig3B_Correlation_allcells.pdf"), width =5.59, height = 3.72 )
# fig2b
# dev.off()


#3. correlate coda estimate with cell proliferation  -------

# Seurat 
pal_cell_cycle <- c("G2M" = alpha("#264653", 1) , "S"= alpha("#264653", 0.7), "G1/G0"="#bc4749")
degs <- get_significant_degs("Age", "cell_type")
cc_composition <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/05_CellCycle/coda_results_cellcycle_seurat_allcells.rds"))
cc_composition$celltype <- gsub("_", " ", cc_composition$celltype)
coda <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/coda_results.rds"))$Age
coda$celltype <- gsub("_", " ", coda$celltype)
cc_composition <- cc_composition %>% full_join(coda, by="celltype")
cc_composition <- cc_composition[!is.na(cc_composition$cc_phase),]
cc_composition<- cc_composition[cc_composition$celltype %in%  unique(degs$celltype), ]
#cc_composition <- cc_composition[!cc_composition$celltype %in%  c("Platelet", "gdT", "NK Proliferating"), ]
cc_composition$cc_phase <- factor(cc_composition$cc_phase, levels=c("G1/G0", "S", "G2M"))
cc_composition <- cc_composition[cc_composition$celltype %in% unique(degs$celltype),]

fig3c <- ggplot(cc_composition[cc_composition$celltype != "gdT", ], aes(y=estimate.x, x=estimate.y))+geom_point(aes(color=cc_phase))+theme+geom_smooth(method = "lm", se = FALSE, aes(color=cc_phase))+  stat_cor(method = "spearman",   p.digits = 3,label.y=0.007, size=3.5)+
  scale_color_manual(values = pal_cell_cycle)+ylab("CoDA estimate cell cycle")+xlab("CoDA estimate cell proportion")+theme(legend.position="none", axis.title=element_text(size=11), axis.text=element_text(size=9.5))+facet_wrap(~cc_phase, nrow = 1)+
  geom_text_repel( data=cc_composition[cc_composition$celltype != "gdT", ], aes(label=celltype), size=3)

# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig3/Fig3C_CorrelationCCylce_allcells.pdf"), width =7.43, height = 3 )
# fig3c
# dev.off()



# Tricycle 

cc_composition <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/05_CellCycle/coda_results_cellcycle_tricycle_allcells.rds"))
coda <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/coda_results.rds"))$Age
coda$celltype <- gsub("_", " ", coda$celltype)
cc_composition <- cc_composition %>% full_join(coda, by="celltype")
cc_composition<- cc_composition[cc_composition$celltype %in%  unique(degs$celltype), ]
cc_composition <- cc_composition[!cc_composition$celltype %in%  c("Platelet", "gdT", "NK Proliferating"), ]
cc_composition$cc_phase <- factor(cc_composition$cc_phase, levels=c("G1/G0", "S", "G2M"))


figS3C <- ggplot(cc_composition[cc_composition$celltype %in%  unique(degs$celltype), ], aes(y=estimate.x, x=estimate.y))+geom_point(aes(color=cc_phase))+theme+geom_smooth(method = "lm", se = FALSE, aes(color=cc_phase))+  stat_cor(method = "spearman",   p.digits = 3,label.y=0.007)+
  scale_color_manual(values = pal_cell_cycle)+ylab("CoDA estimate cell cycle")+xlab("CoDA estimate cell proportion")+theme(legend.position="none")+
  geom_text_repel( data=cc_composition[cc_composition$celltype %in%  unique(degs$celltype), ], aes(label=celltype), size=3.5)+facet_wrap(~cc_phase, nrow = 3)

# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS3/FigS3C_CorrelationCCylce_tricycle_allcells.pdf"), height =7.43, width = 3 )
# figS3C
# dev.off()




############ SUPPLEMENTARY ##################


# 4. CoDA results from cell cycle scoring ------

# Seurat 
tidy_mod.df <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/05_CellCycle/coda_results_cellcycle_seurat_allcells.rds"))
tidy_mod.df$cc_phase <- factor(tidy_mod.df$cc_phase, levels=c("G1/G0", "S", "G2M"))
tidy_mod.df$celltype <- gsub("NK CD56bright", "NK_CD56bright", tidy_mod.df$celltype)
tidy_mod.df <- reorder_cells(tidy_mod.df, reverse = F, neworder = T)

figS3A <- ggplot(tidy_mod.df, aes(x=estimate, y=cc_phase)) + 
  geom_point(aes(alpha=significance,color=direction), size=3) + xlab("Estimate")+ylab(NULL)+
  geom_pointrange(aes(xmin=conf.low, xmax=conf.high, alpha=significance, color=direction), fatten = .00) +
  geom_vline(xintercept=0, linetype = "dashed")  + 
  theme + theme(axis.text = element_text(size = 9), axis.title=element_text(size=12), legend.position="right") +
  scale_fill_manual(values = c("up"=red, "down"=blue))+scale_color_manual(values=c("pos"=red, "neg"=blue))+ theme( strip.placement = "outside", # Place facet labels outside x axis labels.
                                                                                                                   strip.background = element_blank(), # Make facet label background white.
                                                                                                                   legend.position = "top",
                                                                                                                   strip.text.y = element_text(angle = 0, hjust = -0))+
  scale_alpha_manual(values=c(0.5, 0.9)) +facet_grid(celltype~., space="free", scale="free")+ scale_x_continuous(limits = c(-0.015, 0.015), breaks = seq(-0.01, 0.01, by = 0.01))

# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS3/figS3A_CoDA_CC_Seurat_allcells.pdf"), height =7, width = 3 )
# figS3A
# dev.off()

seurat_cc <- tidy_mod.df



#Tricycle 
tidy_mod.df <-  readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/05_CellCycle/coda_results_cellcycle_tricycle_allcells.rds"))
tidy_mod.df <- reorder_cells(tidy_mod.df, reverse = F, neworder = T)
tidy_mod.df$cc_phase <- factor(tidy_mod.df$cc_phase, levels=c("G1/G0", "S", "G2M"))
figS3B <- ggplot(tidy_mod.df, aes(x=estimate, y=cc_phase)) + 
  geom_point(aes(alpha=significance,color=direction), size=3) + xlab("Estimate")+ylab(NULL)+
  geom_pointrange(aes(xmin=conf.low, xmax=conf.high, alpha=significance, color=direction), fatten = .00) +
  geom_vline(xintercept=0, linetype = "dashed")  + 
  theme + theme(axis.text = element_text(size = 9), axis.title=element_text(size=12), legend.position="right") +
  scale_fill_manual(values = c("up"=red, "down"=blue))+scale_color_manual(values=c("pos"=red, "neg"=blue))+ theme( strip.placement = "outside", # Place facet labels outside x axis labels.
                                                                                                                   strip.background = element_blank(), # Make facet label background white.
                                                                                                                   legend.position = "top",
                                                                                                                   strip.text.y = element_text(angle = 0, hjust = -0))+
  scale_alpha_manual(values=c(0.5, 0.9)) +facet_grid(celltype~., space="free", scale="free")+ scale_x_continuous(limits = c(-0.015, 0.015), breaks = seq(-0.01, 0.01, by = 0.01))

# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS3/figS3B_CoDA_CC_Tricycle_allcells.pdf"), height =7, width = 3 )
# figS3B
# dev.off()

tricyle <- tidy_mod.df

# 5. GSEA enrichments ------
# Quiescence 

model_enrich_df <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects//09b_Quiescence/DifferentialEnrichment_Quiescence_AUCell_AllCells.rds"))
model_enrich_df$celltype <- rownames(model_enrich_df)
#model_enrich_df <- model_enrich_df[model_enrich_df$celltype %in% unique(deg$celltype),]
model_enrich_df$celltype <- rownames(model_enrich_df)
model_enrich_df$p.adjust <- p.adjust(model_enrich_df$p_value, method = "fdr")

model_enrich_df$logpval <- -log10(model_enrich_df$p.adjust)
model_enrich_df$significance <- ifelse(model_enrich_df$p.adjust < 0.05, "ss", "ns")
model_enrich_df <- model_enrich_df[!is.na(model_enrich_df$significance),]
figS3d <- ggplot(model_enrich_df, aes(x=estimate, y=logpval))+geom_point(aes(color=significance))+theme+ylab("-log10(p.adjust)")+ggtitle("AUCell - Quiescence")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+ scale_color_manual(values=c("ns"="lightgrey", ss=alpha(blue, 0.8)))+xlab("DE estimate")+
  geom_text(data=model_enrich_df[model_enrich_df$significance == "ss",], aes(x=estimate+0.00001, y=logpval+0.1, label=celltype))

# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS3/FigS3D_DEQuiescence.pdf"), width =5, height = 3.72 )
# figS3d
# dev.off()

# compositional anlaysis 
coda <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/coda_results.rds"))$Age
coda$celltype <- gsub("_", " ", coda$celltype)
coda$celltype <- gsub(" CD56bright", "_CD56bright", coda$celltype)

corr <- merge(coda, model_enrich_df, by = "celltype")
corr <- reorder_cells(corr)
corr <- corr[corr$celltype  %in%  unique(deg$celltype),]

figS3g <- ggplot(corr, aes(y=estimate.y, x=estimate.x))+geom_point(size=2, color=blue)+theme+  geom_smooth(method = "lm", se = FALSE, color=blue)+ylab("DE estimate")+  ggtitle("Quiescence")+
  geom_text_repel( data=corr, aes(x=estimate.x, y=estimate.y, label=celltype))+geom_vline(xintercept = 0, color="#777777" , linetype="dashed", size=0.5)+
  scale_alpha_manual(values=c("ss"=1, "ns"=0.5))+xlab("CoDA estimate")+theme(aspect.ratio=1, legend.position="top")+
  stat_cor(method = "spearman",   p.digits = 1)


pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS3/FigS3G_Correlation_Quiescence.pdf"), width =5, height = 3.72 )
figS3g
dev.off()

stats_quiesc <- model_enrich_df

# Senescence - SenMayo

model_enrich_df <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects//09_Senescence/DifferentialEnrichment_Senescence_AUCell_AllCells.rds"))
model_enrich_df$celltype <- rownames(model_enrich_df)
model_enrich_df <- model_enrich_df[model_enrich_df$celltype %in% unique(deg$celltype),]
model_enrich_df$p.adjust <- p.adjust(model_enrich_df$p_value, method = "fdr")

model_enrich_df$logpval <- -log10(model_enrich_df$p.adjust)
model_enrich_df$significance <- ifelse(model_enrich_df$p.adjust < 0.05, "ss", "ns")
model_enrich_df <- model_enrich_df[!is.na(model_enrich_df$significance),]
figS3e <- ggplot(model_enrich_df, aes(x=estimate, y=logpval))+geom_point(aes(color=significance))+theme+ylab("-log10(p.value)")+ggtitle("AUCell - Senescence ")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+ scale_color_manual(values=c("ns"="lightgrey", ss=alpha(blue, 0.8)))+xlab("DE estimate")+
  geom_text_repel(data=model_enrich_df[model_enrich_df$significance == "ss",], aes(x=estimate+0.00001, y=logpval+0.1, label=celltype))

model_enrich_df[model_enrich_df$estimate <0,]
# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS3/FigS3E_DESenMayo.pdf"), width =5, height = 3.72 )
# figS3e
# dev.off()

#compositional analysis results 
coda <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/coda_results.rds"))$Age
coda$celltype <- gsub("_", " ", coda$celltype)
coda$celltype <- gsub(" CD56bright", "_CD56bright", coda$celltype)

corr <- merge(coda, model_enrich_df, by = "celltype")
corr <- reorder_cells(corr)
corr <- corr[corr$celltype  %in%  unique(deg$celltype),]


figS3h <- ggplot(corr, aes(y=estimate.y, x=estimate.x))+geom_point(size=2, color=blue)+theme+  geom_smooth(method = "lm", se = FALSE, color=blue)+ylab("DE estimate")+  ggtitle("AUCell - Senescence")+
  geom_text_repel( data=corr, aes(x=estimate.x, y=estimate.y, label=celltype))+geom_vline(xintercept = 0, color="#777777" , linetype="dashed", size=0.5)+
  scale_alpha_manual(values=c("ss"=1, "ns"=0.5))+xlab("CoDA estimate")+theme(aspect.ratio=1, legend.position="top")+
  stat_cor(method = "spearman",   p.digits = 1)
pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS3/FigS3H_Correlation_SenMayo.pdf"), width =5, height = 3.72 )
figS3h
dev.off()

stats_df_senmayo <- model_enrich_df


# Senescence - SenCID 
stats_df <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/14_SenCID/GLMM_ model_results.rds"))
stats_df <- stats_df[stats_df$celltype %in% unique(degs$celltype),]
stats_df$fdr <- p.adjust(stats_df$p_value, method = "fdr")
stats_df$logpval <- -log10(stats_df$fdr)
figS3f <- ggplot(stats_df, aes(x=estimate, y=logpval))+geom_point(aes(color=significance))+theme+ylab("-log10(p.adjusted)")+xlab("DE estimate")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+ scale_color_manual(values=c("ns"="lightgrey", ss=alpha(blue, 0.8)))+ggtitle("SenCID - Senescence")+
  geom_text_repel(data=stats_df[stats_df$significance == "ss",], aes(x=estimate+0.00001, y=logpval+0.1, label=celltype))

# 
# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS3/FigS3F_DE_SenCID.pdf"), width =5, height = 3.72 )
# figS3f
# dev.off()


coda <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/coda_results.rds"))$Age
coda$celltype <- gsub("_", " ", coda$celltype)
coda$celltype <- gsub(" CD56bright", "_CD56bright", coda$celltype)

corr <- merge(coda, stats_df, by = "celltype")
corr <- reorder_cells(corr)
corr <- corr[corr$celltype %in% unique(degs$celltype),]

figS3i <- ggplot(corr, aes(y=estimate.y, x=estimate.x))+geom_point(size=2, color=blue)+theme+  geom_smooth(method = "lm", se = FALSE, color=blue)+ylab("Estimate SenCID")+  
  geom_text_repel( data=corr, aes(x=estimate.x, y=estimate.y, label=celltype))+geom_vline(xintercept = 0, color="#777777" , linetype="dashed", size=0.5)+
  scale_alpha_manual(values=c("ss"=1, "ns"=0.5))+xlab("CoDA estimate")+theme(aspect.ratio=1, legend.position="top")+ggtitle("SenCID - Senescence")+
  stat_cor(method = "spearman",   p.digits = 1)

pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS3/FigS3I_Correlation_SenCID.pdf"), width =5, height = 3.72 )
figS3i
dev.off()

stats_df_sencid <- stats_df

# correlation DEA and CoDA Oelen 2022 
degs_oelen_all <- readRDS(paste0(data_path, "msopena//01_meta-analysis_Age_Sex/robjects/02_DEA/DEG_Age_cell_type_v2_0.2.rds"))
#modify some names for plotting 
degs_oelen_all$celltype <- gsub("_", " ", degs_oelen_all$celltype)
degs_oelen_all$direction <- ifelse(degs_oelen_all$logFC > 0, "up", "down")
degs_oelen_all$gene <- degs_oelen_all$ID
degs_oelen <- degs_oelen_all[degs_oelen_all$P.Value < 0.01, ]
result_df_oelen <- data.frame(celltype = character(0), binom.test.p = numeric(0), prop_up = numeric(0))
for (cell in unique(degs_oelen$celltype)){ 
  print(cell)
  binom <- binom.test(length(degs_oelen[degs_oelen$direction == "down" & degs_oelen$celltype == cell,]$gene) , length(degs_oelen[degs_oelen$celltype == cell,]$gene))$p.value
  upregulated <- length(degs_oelen[degs_oelen$direction == "up" & degs_oelen$celltype == cell,]$gene)
  downregulated <- length(degs_oelen[degs_oelen$direction == "down" & degs_oelen$celltype == cell,]$gene)
  
  success_ratio = upregulated / (downregulated+upregulated)
  result_df_oelen <- rbind(result_df_oelen, data.frame(celltype = cell, binom.test.p = binom, prop_up = success_ratio))
}
result_df_oelen$p.val_dir <- ifelse(result_df_oelen$prop_up > 0.5, result_df_oelen$binom.test.p , -result_df_oelen$binom.test.p )
result_df_oelen$direction <- ifelse(result_df_oelen$prop_up > 0.5, "upregulation_bias" , "downregulation_bias" )
result_df_oelen$binom.test.fdr <- p.adjust(result_df_oelen$binom.test.p, method = "fdr")

coda_oelen <- readRDS(paste0(data_path, "/aripol1/wijst-2020-hg19/v1/aging/25.cell_type_proportions.CMV/all_data/v2/cell_type/CLR/phenotype_stats.rds"))$Age
coda_oelen$celltype <- gsub("_", " ", rownames(coda_oelen))

coda_dea_oelen <- coda_oelen %>% left_join(result_df_oelen, by="celltype")
coda_dea_oelen <- coda_dea_oelen[coda_dea_oelen$celltype %in% unique(deg$celltype),]
coda_dea_oelen$significance <- ifelse(coda_dea_oelen$fdr < 0.05, "ss", "ns")
cor.test(coda_dea_oelen$estimate, coda_dea_oelen$prop_up)

ggplot(coda_dea_oelen, aes(y=prop_up, x=estimate))+geom_point( color=blue)+theme+  
  geom_smooth(method = "lm", se = FALSE, color=blue)+ylab("Proportion up-regulated DEGs")+  scale_size_continuous(range = c(0.5, 3))+
  geom_text_repel( data=coda_dea_oelen, aes(x=estimate, y=prop_up, label=celltype), size=3.1)+geom_vline(xintercept = 0, color="#777777" , linetype="dashed", size=0.5)+
  scale_alpha_manual(values=c("ss"=1, "ns"=0.5))+xlab("CoDA estimate")+theme(legend.position="right", axis.title=element_text(size=11), legend.title=element_text(size=10))+
  stat_cor(method = "spearman", label.y = 1, label.x = -0.05,   p.digits = 1, size=3.4)+labs(alpha= "CoDA singnificance", size="-log10(binom.FDR)")




#save supplementary table 2 
coda_oelen <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/coda_results.rds"))$Age
coda_oelen$celltype <- gsub("_", " ", coda$celltype)
deg <- get_significant_degs("Age", "cell_type")
coda <- coda[coda$celltype %in% unique(degs$celltype),]
coda$fdr <- p.adjust(coda$p.value, method = "fdr")
coda$celltype <- gsub(" CD56bright", "_CD56bright", coda$celltype)

library(openxlsx)
sheets<- list("Compositional Analysis"= coda, "DEA Oelen 2022"=degs_oelen, "CoDA Oelen 2022" = coda_oelen, "CellCycle_Seurat" = seurat_cc, "CellCycle_Tricycle"=tricyle,  "Quiescence"=stats_quiesc, "Senescence_SenMayo" = stats_df_senmayo, "Senescence_SenCID"= stats_df_sencid)
write.xlsx(sheets, paste0(data_path, '/msopena/02_OneK1K_Age/supplementary_tables/TableS3_Compositional_Analysis.xlsx'))



