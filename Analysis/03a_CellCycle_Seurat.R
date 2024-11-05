#!/usr/bin/env Rscript

# setting working directory (cluster or local)
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


# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--celltype"), action="store", default=NA, type='character',
              help="cell type in l2"))
opt = parse_args(OptionParser(option_list=option_list))
celltype <- opt$celltype

print(paste0("--------- Running Cell Cycle Annotation ---------"))
print(paste0("Celltype: ", celltype))
# Packages
shhh(library(Seurat))
shhh(library(SeuratObject))
shhh(library(data.table))
shhh(library(plyr))
shhh(library(reshape2))
shhh(library(dplyr))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(ggplot2))
shhh(library(ggrepel))
shhh(library(ggpubr))
shhh(library(RColorBrewer))
shhh(library(org.Hs.eg.db))


################################## Set Variables and load Data ################################## 
# Main directory
main_dir <- paste0(data_path, "/aripol1/wijst-2020-hg19/v1/aging/")

# Output directory
plotsdir <- paste0(data_path,  "msopena/02_OneK1K_Age/robjects/05_CellCycle/")
if(!dir.exists(plotsdir)){dir.create(plotsdir, recursive = T)}
robjectsdir <- paste0(data_path,  "msopena/02_OneK1K_Age/plots/05_CellCycle/")
if(!dir.exists(robjectsdir)){dir.create(robjectsdir, recursive = T)}

# assay_var <- assay_var.vec[1]
pbmc <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/sce_data_objects/", celltype, "_cell_type_so_harmony.rds"))
print(" Seurat object loaded ----")

################################## Functions ##################################
source(paste0(basepath, '/Projects/GTEx_v8/aripol1/utils/Rscripts/global_functions.R'))
source(paste0(data_path, "/msopena/02_OneK1K_Age/scripts/functions.R"))


CellCycleScoring_1 <- function(object.cc, s.features, g2m.features, ctrl = NULL, set.ident = FALSE) {
  name <- "Cell.Cycle"
  features <- list(S.Score = s.features, G2M.Score = g2m.features)
  if (is.null(x = ctrl)) {
    ctrl <- min(vapply(X = features, FUN = length, FUN.VALUE = numeric(length = 1)))
  }
 # object.cc <- AddModuleScore(object = object, features = features, 
  #                            name = name, ctrl = ctrl, ...)
  cc.columns <- grep(pattern = name, x = colnames(x = object.cc[[]]), 
                     value = TRUE)
  cc.scores <- object.cc[[cc.columns]]
 # rm(object.cc)
  CheckGC()
  assignments <- apply(X = cc.scores, MARGIN = 1, FUN = function(scores, 
                                                                 first = "S", second = "G2M", null = "G1") {
    if (all(scores < 0)) {
      return(null)
    }
    else {
      if (length(which(x = scores == max(scores))) > 1) {
        return("Undecided")
      }
      else {
        return(c(first, second)[which(x = scores == max(scores))])
      }
    }
  })
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments), 
                     by = 0)
  colnames(x = cc.scores) <- c("rownames", "S.Score", "G2M.Score", 
                               "Phase")
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c("S.Score", "G2M.Score", "Phase")]
  object.cc[[colnames(x = cc.scores)]] <- cc.scores
  if (set.ident) {
    object.cc[["old.ident"]] <- Idents(object = object.cc)
    Idents(object = object.cc) <- "Phase"
  }
  return(object.cc)
}

################################## Analysis ##################################

# Cell-Cycle Scoring and Regression: https://hbctraining.github.io/scRNA-seq/lessons/cell_cycle_scoring.html
## Read a list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can segregate this list into markers of G2/M phase and markers of S phase
cc_fn <- paste0(main_dir, 'seurat/cell_cycle_vignette_files/Homo_sapiens.csv')
cell_cycle_genes.df <- read.csv(cc_fn)

gene_ens <- unique(cell_cycle_genes.df$geneID)
gene_symbols <- mapIds(org.Hs.eg.db, 
                   keys = gene_ens, keytype = "ENSEMBL", column="SYMBOL")
gene_symbols <- gene_symbols[!is.na(gene_symbols)]
gene_symbols.df <- data.frame(geneID = names(gene_symbols),
                              gene_id = unname(gene_symbols))
cell_cycle_markers <- merge(cell_cycle_genes.df, gene_symbols.df, by = 'geneID')
colnames(cell_cycle_markers)[colnames(cell_cycle_markers)=='geneID'] <- 'ENSEMBL'
colnames(cell_cycle_markers)[colnames(cell_cycle_markers)=='gene_id'] <- 'SYMBOL'

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S" & SYMBOL %in% rownames(pbmc) )%>%
  pull("SYMBOL") 

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M" & SYMBOL %in% rownames(pbmc) ) %>%
  pull("SYMBOL")

# Assign to Seurat
features <- list('S.Score' = s_genes, 'G2M.Score' = g2m_genes)
nbin_var <- 24
nbin_vars <- nbin_var:1
for(i in nbin_vars){
  print(i)
  pbmc_check <- try(AddModuleScore(object = pbmc, 
                                   features = features,
                                   nbin = i,
                                   name = 'Cell.Cycle'))
  if(class(pbmc_check)!='try-error'){break}
}
pbmc <- pbmc_check
pbmc_rna <- CellCycleScoring_1(pbmc, s.features = s_genes, g2m.features = g2m_genes)

# Visualize the PCA, grouping by cell cycle phase
umap.Phase <- DimPlot(pbmc_rna,
                      reduction = "umap",
                      group.by= "Phase")
umap.Phase_fn <- paste0(plotsdir, celltype,".UMAP.Phase.png")
ggsave(umap.Phase_fn, umap.Phase, width = 5, height = 4)
print("UMAP saved")

umap.Phase <- DimPlot(pbmc_rna,
                      reduction = "harmony",
                      group.by= "Phase")
umap.Phase_fn <- paste0(plotsdir,celltype,".celltype.Phase.png")
ggsave(umap.Phase_fn, umap.Phase, width = 5, height = 4)
print("UMAP saved")


# Proportion of Phase - all genes
md <- pbmc_rna@meta.data %>% dplyr::select( S.Score, G2M.Score, Phase)
md$celltype <- celltype
saveRDS(md, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/05_CellCycle/", celltype,"_CellCycleInfo.rds" ))
ct_phase <- prop.table(table(md$celltype, md$Phase), margin = 1)
saveRDS(ct_phase, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/05_CellCycle/", celltype,"_CellCycleInfo_table.rds" ))



#merge CellCycle info 
cells_to_keep <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/cells_to_keep.rds"))

cc_mdata <- do.call(rbind.data.frame, lapply(gsub(" ", "_", cells_to_keep), function(celltype) 
                      readRDS( paste0(data_path, "/msopena/02_OneK1K_Age/robjects/05_CellCycle/", celltype,"_CellCycleInfo.rds" ))))

saveRDS(cc_mdata, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/05_CellCycle/All_Cells_CellCycleInfo.rds" ))



# perform compositional analysis with the cell types with DEGs
mdata_all <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/metadata_processed.rds"))

# get proportions 
compute_proportions <- function(celltype){
  print(celltype)
  mdata_cc <- readRDS( paste0(data_path, "/msopena/02_OneK1K_Age/robjects/05_CellCycle/", celltype,"_CellCycleInfo.rds" ))
  #merge with assignment info 
  mdata <- merge(mdata_cc, mdata_all, by="row.names")
  prop <- as.data.frame(prop.table(table(mdata$assignment, mdata$Phase.y), margin = 1))
  colnames(prop) <- c("assignment", "CC_Phase", "freq")
  prop <- prop %>% left_join(mdata_all[!duplicated(mdata_all$assignment),], by="assignment") %>% dplyr::select(c("assignment", "CC_Phase", "freq", "Age", "Gender", "date"))
  prop$celltype <- celltype
  return(prop)
}

degs <- get_significant_degs("Age", "cell_type")

prop_list <- lapply(gsub(" ", "_", unique(degs$celltype)), function(celltype)compute_proportions(celltype))
props.df_all <- do.call(rbind.data.frame, prop_list)
saveRDS(props.df_all, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/05_CellCycle/proportions_cellcycle_seurat.rds"))

props.df_all <-readRDS( paste0(data_path, "/msopena/02_OneK1K_Age/robjects/05_CellCycle/proportions_cellcycle_seurat.rds"))
props.df_all$CC_Phase <- as.character(props.df_all$CC_Phase)
props.df_all[props.df_all$CC_Phase =="G1",]$CC_Phase <- "G1/G0"

model_CC_propotions <- function(celltype) {
  celltype <- gsub(" ", "_", celltype)
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

tidy_mod.list <- sapply(unique(degs$celltype), function(celltype) model_CC_propotions(celltype), simplify = FALSE)
tidy_mod.df <- do.call(rbind.data.frame, tidy_mod.list)
for(cc_phase in unique(tidy_mod.df$cc_phase)){
  tidy_mod.df[tidy_mod.df$cc_phase == cc_phase, ]$fdr <-  p.adjust(tidy_mod.df[tidy_mod.df$cc_phase == cc_phase,]$p.value, 'fdr')}
tidy_mod.df$celltype <- gsub("_", " ", tidy_mod.df$celltype)
tidy_mod.df$cc_phase <- factor(tidy_mod.df$cc_phase, levels=rev(c("G1/G0", "S", "G2M")))
tidy_mod.df$significance <- ifelse(tidy_mod.df$fdr < 0.05, "ss", "ns")


saveRDS(tidy_mod.df, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/05_CellCycle/coda_results_cellcycle_seurat_allcells.rds"))












