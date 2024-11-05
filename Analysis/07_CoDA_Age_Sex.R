#!/usr/bin/env Rscript

# Paths
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

robjects_out <- paste0(data_path, "/msopena/03_Menopause/robjects/01_DEA/")
plots_out <-  paste0(data_path, "/msopena/03_Menopause/plots/")

# Parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
shhh(library(Seurat))
shhh(library(SeuratDisk))
shhh(library(SingleCellExperiment))
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(ggplot2))
shhh(library(dreamlet))
shhh(library(zenith))
shhh(library(scater))
shhh(library(RColorBrewer))
Csparse_validate = "CsparseMatrix_validate"

shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--sex"), action="store", default=NA, type='character',
              help="F or M"))
opt = parse_args(OptionParser(option_list=option_list))
sex <- opt$sex

# read in files 
metadata <-readRDS(paste0(basepath,  "Data/scRNAseq/Yazar2022/16.Azimuth_celltype_proportions/cell_type/metadata.rds"))
metadata[metadata$sex == 1,]$sex <- "M"
metadata[metadata$sex == 2,]$sex <- "F"
metadata$cell_type <- metadata$predicted.celltype.l2
metadata$donor <- metadata$individual
mdata_donor <- metadata[!duplicated(metadata$donor), ]
colnames(mdata_donor) <- gsub("assignment", "donor", colnames(mdata_donor))
#order_cells<- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/new_order_cells.rds"))
props.df <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/12_CellProportions_Sex/proportions_", sex, ".rds" ))

# ## Calculate proportions
# calculate_props <- function(sex){
#   metadata[metadata$sex == sex,]%>%
#     group_by_at(c("individual", "cell_type")) %>%
#     summarise(n = n()) -> count.df
#   count.df %>%
#     group_by(individual) %>%
#     mutate(freq = n / sum(n)) %>% as.data.frame() -> prop.df
#   colnames(prop.df) <- c('donor', 'cell_type', 'n', 'freq')
#   prop.df$dataset <- 'OneK1K'
#   prop.df %>%
#     group_by(cell_type, donor) %>%
#     summarise(n = n()) -> check_prop.df
#   saveRDS(prop.df, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/12_CellProportions_Sex/proportions_", sex, ".rds" ))
# 
# }
# lapply(c("M", "F"), function(sex) calculate_props(sex))
# Add metadata to the proportions

# colnames(mdata_donor) <- gsub("assignment", "donor", colnames(mdata_donor))
# props.df <- proportions %>% left_join(mdata_donor, by="donor") %>% select(c(colnames(proportions), "Age", "date"))
# saveRDS(props.df, (paste0(data_path, "/msopena/02_OneK1K_Age/robjects/12_CellProportions_Sex/proportions_", sex, ".rds" )))



# Summarise info by celltype
props.df %>%
  group_by(cell_type, .drop = FALSE) %>%
  summarise(n = sum(n)) -> count.df
count.df %>%
  mutate(freq = n / sum(n)) %>% 
  arrange(desc(freq)) %>% as.data.frame() -> celltype.stats.df
celltype.stats.df$celltype.label <- paste0(celltype.stats.df$cell_type,
                                           ' (n=', celltype.stats.df$n, ';freq=', as.character(round(celltype.stats.df$freq,2)), ')')


# Transform data
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
  DF_proportions <- reshape2::dcast(df, donor ~ cell_type, value.var = "freq") %>% as.data.frame()
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


## Function
# cell_type <- celltypes[1]
# covs_df <- covs.df
df <- props_clr.df
md = mdata_donor
# interaction = opt$interaction
lm_by_ct <- function(cell_type, df = props_clr.df, md = mdata_donor){
  print(cell_type)
  vec_i <- df[rownames(df)==cell_type, ]
  df_i <- as.data.frame(vec_i)
  colnames(df_i) <- 'freq'
  df_i$donor <- rownames(df_i)
  df_i <- merge(df_i, md, by = 'donor')
  rownames(df_i) <- df_i$donor
  df_i <- df_i[,-1]
  
  
  # Formula
  fmla <- paste(c("age", "(1|pool)"),collapse='+')
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
  tidy_mod.df$celltype <- cell_type
  return(tidy_mod.df)
}

## Apply function
tidy_mod.list <- sapply(unique(props.df$cell_type[!is.na(props.df$cell_type)]), function(i) lm_by_ct(i), simplify = FALSE)
tidy_mod.df <- do.call("rbind", tidy_mod.list)

# Split by phenotype and compute FDR
tidy_mod.by_phe <- split(tidy_mod.df, tidy_mod.df$term)
tidy_mod.Age <- tidy_mod.by_phe$age
tidy_mod.Age$fdr <- p.adjust(tidy_mod.Age$p.value, 'fdr')
tidy_mod.Age <- tidy_mod.Age[order(tidy_mod.Age$fdr),]
table(tidy_mod.Age$fdr < 0.05)


# List to DF

tidy_mod.Age$direction <- ifelse(tidy_mod.Age$estimate>0, 'pos', 'neg')
tidy_mod.Age <- merge(tidy_mod.Age, celltype.stats.df[,c('cell_type','celltype.label')], by.x = 'celltype', by.y="cell_type")
celltype.label <- unique(celltype.stats.df$celltype.label)

# Check results
phe_stats.list.fn <- paste0(data_path, "/msopena/02_OneK1K_Age/robjects/12_CellProportions_Sex/CODA_resultsAge_", sex, ".rds")
saveRDS(tidy_mod.Age, phe_stats.list.fn)


