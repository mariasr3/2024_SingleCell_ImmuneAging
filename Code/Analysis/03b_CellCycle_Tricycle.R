
library(tricycle);library(dplyr);library(GSEABase);library(stringr);library(Seurat);library(SingleCellExperiment);

Csparse_validate = "CsparseMatrix_validate"

# Compute cell cyclce score using Tricycle 

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
robjects <- paste0(data_path, "msopena/02_OneK1K_Age/robjects/")
path_senec <- paste0(robjects, "/08_EnrichmentScores/")
#dir.create(path_enrichScore, recursive = TRUE)

#functions
source(paste0(data_path, "/msopena/02_OneK1K_Age/scripts/functions.R"))

#metadata <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/metadata_processed.rds"))
order_cells_old<- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/order_cells.rds"))
order_cells<- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/new_order_cells.rds"))
celltype_l1 <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/celltypes_equivalence.rds"))
celltype_l1$cell_type <- factor(celltype_l1$cell_type, levels = order_cells$cell_type)




metadata <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/metadata_processed.rds"))


# 1. Get expression data 
print("1. Loading data ----------")

sce <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/sce_data_objects/AllCells_predicted.celltype.l1_sceraw.rds"))
sce <- project_cycle_space(sce, gname.type="SYMBOL", species="human")
sce <- estimate_cycle_position(sce)
mdata <- colData(sce)
saveRDS( mdata,paste0(data_path, "/msopena/02_OneK1K_Age/robjects/08_EnrichmentScores/Metadata_CellCycle_tricycle.rds"))

mdata_cc <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/08_EnrichmentScores/Metadata_CellCycle_tricycle.rds"))

mdata_cc$tricyclePosition <- as.numeric(mdata_cc$tricyclePosition)
mdata_cc$cell_cycle <- NA
mdata_cc[mdata_cc$tricyclePosition >= 0.5*pi & mdata_cc$tricyclePosition < pi,]$cell_cycle <- "S"
mdata_cc[mdata_cc$tricyclePosition >=pi & mdata_cc$tricyclePosition <1.75*pi, ]$cell_cycle <- "G2M"
mdata_cc[mdata_cc$tricyclePosition <= 0.25*pi  | mdata_cc$tricyclePosition >= 1.75*pi,]$cell_cycle <- "G1/G0"  


# run CoDA

compute_proportions <- function(celltype){
  print(celltype)
  celltype <- gsub("NK CD56bright", "NK_CD56bright", celltype)
  mdata_cc <- mdata_cc[!is.na(mdata_cc$cell_cycle),]
  mdata <- mdata_cc[mdata_cc$cell_type == celltype, ]
  prop <- as.data.frame(prop.table(table(mdata$assignment, mdata$cell_cycle), margin = 1))
  colnames(prop) <- c("assignment", "CC_Phase", "freq")
  prop <- prop %>% left_join(mdata_cc[!duplicated(mdata_cc$assignment),], by="assignment") %>% dplyr::select(c("assignment", "CC_Phase", "freq", "Age", "Gender", "date"))
  prop$celltype <- celltype
  return(prop)
}

prop_list <- lapply(unique(mdata_cc$cell_type), function(celltype)compute_proportions(celltype))
props.df_all <- do.call(rbind.data.frame, prop_list)
saveRDS(props.df_all, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/05_CellCycle/proportions_cellcycle_tricycle.rds"))

props.df_all <-readRDS( paste0(data_path, "/msopena/02_OneK1K_Age/robjects/05_CellCycle/proportions_cellcycle_tricycle.rds"))


model_CC_propotions <- function(celltype) {
  celltype <- gsub("NK CD56bright", "NK_CD56bright", celltype)
  print(celltype)
  props.df <- props.df_all[props.df_all$celltype == celltype, ]
  Geom_mean <- function(x){
    exp(mean(log(x)))
  }
  CLR <- function(D){
    log2(D / Geom_mean(D))
  }
  props.df[is.na(props.df$cell_type),]
  df <- props.df
  clr_func <- function(df){
    # Long to wide --> rows = cell_type, cols = donor, value = freq
    DF_proportions <- reshape2::dcast(df, assignment ~ CC_Phase, value.var = "freq") %>% as.data.frame()
    DF_proportions[is.na(DF_proportions)] <- 0
    rownames(DF_proportions) <- DF_proportions[,1]
    DF_proportions <- DF_proportions[,-1]
    DF_proportions <- t(DF_proportions)
    DF_proportions[is.na(DF_proportions)] <- 0
    
    # Add pseudocount
    pseudocount <- 1/3*min(DF_proportions[!DF_proportions==0])
    
    # CLR
    DF_proportions %>% 
      apply(2, function(x){CLR(x+pseudocount)}) -> DF_proportions.clr
    
    # # Check
    # apply(DF_proportions, 1, summary)
    # apply(DF_proportions.clr, 1, summary)
    
    return(DF_proportions.clr)
  }
  
  ## Apply main function
  props_clr.df <- clr_func(props.df)
  
  lm_by_ct <- function(cc_phase, df = props_clr.df, md = metadata[!duplicated(metadata$assignment),]){
    print(cc_phase)
    vec_i <- df[rownames(df)==cc_phase, ]
    df_i <- as.data.frame(vec_i)
    colnames(df_i) <- 'freq'
    df_i$assignment <- rownames(df_i)
    df_i <- merge(df_i, md, by = 'assignment')
    rownames(df_i) <- df_i$assignment
    df_i <- df_i[,-1]
    
    
    # Formula
    fmla <- paste(c("Age", "Gender","(1|date)"),collapse='+')
    fmla <- paste0('freq ~ ',fmla)
    form <- as.formula(fmla)
    print(paste0('Fitting lmer: ',fmla))
    
    # fit model
    print('lmer...')
    mod <-  lmerTest::lmer(form, data = df_i)
    # tidy model
    tidy_mod <- broom.mixed::tidy(mod, conf.int = TRUE, effects = "fixed")
    
    
    # tidy to dataframe
    cnames <- c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high")
    tidy_mod <- tidy_mod[,which(colnames(tidy_mod)%in%cnames)]
    tidy_mod.df <- as.data.frame(tidy_mod)
    tidy_mod.df$cc_phase <- cc_phase
    return(tidy_mod.df)
  }
  
  ## Apply function
  tidy_mod.list <- sapply(c( "G2M", "S","G1/G0" ), function(i) lm_by_ct(i), simplify = FALSE)
  tidy_mod.df <- do.call("rbind", tidy_mod.list)
  tidy_mod.by_phe <- split(tidy_mod.df, tidy_mod.df$term)
  tidy_mod.Age <- tidy_mod.by_phe$Age
  tidy_mod.Age$fdr <- p.adjust(tidy_mod.Age$p.value, 'fdr')
  tidy_mod.Age <- tidy_mod.Age[order(tidy_mod.Age$fdr),]
  table(tidy_mod.Age$fdr < 0.05)
  
  
  # List to DF
  
  tidy_mod.Age$direction <- ifelse(tidy_mod.Age$estimate>0, 'pos', 'neg')
  tidy_mod.Age$celltype <- celltype 
  return(tidy_mod.Age)}

tidy_mod.list <- sapply(order_cells$cell_type[!order_cells$cell_type %in% c("CD4 Proliferating", "Doublet", "cDC1", "CD8 Proliferating", "NK Proliferating")], function(celltype) model_CC_propotions(celltype), simplify = FALSE)
tidy_mod.df <- do.call(rbind.data.frame, tidy_mod.list)
for(cc_phase in unique(tidy_mod.df$cc_phase)){
  tidy_mod.df[tidy_mod.df$cc_phase == cc_phase, ]$fdr <-  p.adjust(tidy_mod.df[tidy_mod.df$cc_phase == cc_phase,]$p.value, 'fdr')}
#tidy_mod.df <- tidy_mod.df[tidy_mod.df$celltype %in%cell_types, ]
tidy_mod.df <- reorder_cells(tidy_mod.df, neworder = T, reverse = F)
tidy_mod.df$cc_phase <- factor(tidy_mod.df$cc_phase, levels=rev(c("G1/G0", "S", "G2M")))
tidy_mod.df$significance <- ifelse(tidy_mod.df$fdr < 0.05, "ss", "ns")

saveRDS(tidy_mod.df, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/05_CellCycle/coda_results_cellcycle_tricycle.rds"))





tidy_mod.df <- tidy_mod.df[tidy_mod.df$celltype %in% degs$celltype,]
ggplot(tidy_mod.df[!tidy_mod.df$celltype %in% c("Platelet","NK_CD56bright", "NK Proliferating", "gdT", "HSPC"), ], aes(x=estimate, y=cc_phase)) + 
  geom_point(aes(alpha=significance,color=direction), size=3) + xlab("Estimate")+ylab(NULL)+
  geom_pointrange(aes(xmin=conf.low, xmax=conf.high, alpha=significance, color=direction), fatten = .00) +
  geom_vline(xintercept=0, linetype = "dashed")  + 
  theme + theme(axis.text = element_text(size = 10), axis.title=element_text(size=12), legend.position="right") +
  scale_fill_manual(values = c("up"=red, "down"=blue))+scale_color_manual(values=c("pos"=red, "neg"=blue))+ theme( strip.placement = "outside", # Place facet labels outside x axis labels.
                                                                                                                   strip.background = element_blank(), # Make facet label background white.
                                                                                                                   legend.position = "none",
                                                                                                                   strip.text.y = element_text(angle = 0, hjust = -0))+
  scale_alpha_manual(values=c(0.5, 0.9)) +facet_grid(celltype~., space="free", scale="free")+ scale_x_continuous(limits = c(-0.015, 0.015), breaks = seq(-0.01, 0.01, by = 0.01))


cc_composition <-tidy_mod.df
coda <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/phenotype_stats.rds"))$Age
coda$celltype <- gsub("_", " ", coda$celltype)
cc_composition <- cc_composition %>% full_join(coda, by="celltype")
cc_composition <- cc_composition[!is.na(cc_composition$cc_phase),]
cc_composition<- cc_composition[cc_composition$celltype %in%  unique(degs$celltype), ]
cc_composition <- cc_composition[!cc_composition$celltype %in%  c("Platelet", "gdT", "NK Proliferating", "HSPC"), ]









