


#functions --------------

#library("AnnotationDbi"); library('org.Hs.eg.db');library(readxl)

#variables
red <-  "#bc4749" 
blue <-  "#264653"
orange <- "#E18335"
purple <- "#855799"


#extract significant DEGs - meta 
get_significant_degs<- function(pheno, cell_level, consortium = F){
  #read datasets for each pheno and cell_level 
  deg <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/DEG_", pheno, "_", cell_level, ".rds"))
  if(consortium){ deg <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/DEG_", pheno, "_", cell_level, ".rds"))}
  #modify some names for plotting 
  deg$direction <- ifelse(deg$logFC > 0, "up", "down")
  #Select significant genes 
  deg_signif <- deg[deg$adj.P.Val < 0.05 ,]
  return(deg_signif)
}

get_significant_degs_metanalysis <- function(pheno, cell_level){
  #read datasets for each pheno and cell_level 
  deg <- readRDS(paste0(data_path, "msopena//01_meta-analysis_Age_Sex/robjects/02_DEA/DEG_pseudobulk_", pheno, "_", cell_level, ".rds"))
  #modify some names for plotting 
  deg$celltype <- gsub("_", " ", deg$celltype)
  deg <- deg[deg$dataset == "meta", ]
  deg$direction <- ifelse(deg$logFC > 0, "up", "down")
  deg <- deg[deg$fdr < 0.05,]
  return(deg)
}

get_significant_degs_MAST<- function(pheno, cell_level){
  mast <- readRDS(paste0(data_path, "/aripol1/wijst-2020-hg19/v1/aging/17.dea_MAST_glmer_split.by_cell_type.DEGs.OneK1Kagecorrected_coefs/", cell_level, "/l1_0_l2_0/", pheno, "/df_fcHurdle.rds"))
  all_genes_mast <- mast$primerid
  mast <- mast[mast$chemistry=='OneK1K_age_corrected',]
  mast <- split(mast, mast$cell_type)
  mast<- lapply(mast, function(x){
    x$fdr_re <- p.adjust(x[[2]], 'fdr')
    x$fdr_re.ss <- ifelse(x$fdr_re<=0.05, 'ss', 'ns')
    return(x)
  })
  mast <- do.call(rbind.data.frame, mast)
  mast <- mast[mast$fdr_re < 0.05,]
  mast$celltype <- gsub("_", " ", mast$cell_type)
  mast$gene <- mast$primerid
  return(mast)
  
}

get_significant_degs_sex<- function(sex, balanced=T){
  #read datasets for each pheno and cell_level 
  if(balanced){
    deg  <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/AllCells_", sex, "_deaTopTable.rds" ))
  }else{
    deg <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/AllCells_", sex, "_unbalanced_deaTopTable.rds" ))
    
  }
  #modify some names for plotting 
  deg$direction <- ifelse(deg$logFC > 0, "up", "down")
  #Select significant genes 
  deg_signif <- deg[deg$adj.P.Val < 0.05 ,]
  return(deg_signif)
}


#annotate chromosoms Sex genes
dea_sex_annotate <- function(degs.df, gencode_v26.fn =paste0(basepath, '/Data/gene_annotation/gencode/release_26/gencode.v26.gene_id.chr_relabel.rds'), oliva_s3.fn = paste0(basepath, '/Projects/GTEx_v8/aripol1/degs_literature/Sex/01.degs_literature/Oliva2020/Table_S3.xlsx')){
  # DEGs symbols to ensembl
  gene_symbols <- unique(degs.df$gene)
  gene_ens <- mapIds(org.Hs.eg.db, keys = gene_symbols, keytype = "SYMBOL", column="ENSEMBL")
  gene_ens.vec <- gene_ens[!is.na(gene_ens)]
  degs.df$gene_ens <- gene_ens.vec[degs.df$gene]
  
  # Read summarised annotation file (chromosome info): Gencode release 19 --> GRCh37 (hg37), used in the seurat object we're using for the DEA
  gencode_v26.df <- readRDS(gencode_v26.fn)
  
  # Read sex-genes annotation from Oliva 2020
  col_types.vec <- c(rep("text",3), rep("numeric", 52))
  oliva_s3.dt <- read_excel(oliva_s3.fn, sheet = "GTEx v8 X-linked sex-biased gen", col_types = col_types.vec)
  oliva_s3.df <- as.data.frame(oliva_s3.dt)
  colnames(oliva_s3.df)[11] <- 'Reported_Escapee'
  oliva_s3.df$gene_ens <- gsub('.[0-9]+$', '', oliva_s3.df$ENSEMBL_gene_id)
  
  # Check and Merge --> Oliva_2020 and Gencode_v19_GRCh37
  table(oliva_s3.df$gene_ens%in%gencode_v26.df$gene_ens)
  annotation.df <- merge(gencode_v26.df, oliva_s3.df, by = 'gene_ens', all = TRUE)
  nrow(annotation.df)==nrow(gencode_v26.df)
  cnames <- c('gene_id', 'ENSEMBL_gene_id', 'gene_id', 'gene_ens', 'HUGO_gene_id', 'Genetype', 'transcript_type', 'chr', 'chr_broad', 'chr_spec', 'Reported_Escapee')
  anno.df <- annotation.df[,colnames(annotation.df)%in%cnames]
  table(!is.na(anno.df$Reported_Escapee), anno.df$chr_spec)
  
  # Add new columns for chr, chr_broad and chr_spec + Reported_Escapee
  cnames_chr <- colnames(anno.df)[grep('chr',colnames(anno.df))]
  anno.out <- anno.df
  i <- cnames_chr[1]
  for(i in cnames_chr){
    print(i)
    cname <- paste0(i, '.escapee')
    anno.out[[cname]] <- ifelse(!is.na(anno.out$Reported_Escapee) & anno.out$Reported_Escapee==1, paste0(anno.out[[i]], '.escapee'), anno.out[[i]])
  }
  anno_genes.df <- anno.out[,-which(colnames(anno.out)%in%c('HUGO_gene_id', 'gene_id'))]
  par_genes <- anno_genes.df[anno_genes.df$chr=='PAR',]$gene_ens
  df_pre_par <- anno_genes.df[anno_genes.df$gene_ens%in%par_genes,]
  df_not_par <- anno_genes.df[!anno_genes.df$gene_ens%in%par_genes,]
  
  # PAR genes
  df_par <- data.frame()
  for(i in par_genes){
    df_i <- df_pre_par[df_pre_par$gene_ens==i,]
    if(nrow(df_i)>=2){
      print('PAR gene...')
      df_i <- df_i[df_i$chr=='PAR',]
      print(df_i$gene_ens)
    }
    df_par <- rbind(df_par, df_i)
  }
  anno_genes_clean.df <- rbind(df_not_par, df_par)
  
  # Merge
  degs_anno.df <- merge(degs.df, anno_genes_clean.df, by = 'gene_ens')
  kept_genes <- round(length(unique(degs_anno.df$gene))/length(unique(degs.df$gene)),3)
  print(paste0('Mapped genes: prop = ', kept_genes, ', n = ', length(unique(degs_anno.df$gene))))
  return(degs_anno.df)
}


#plot beesward DA analysis
plot_beeswarm <- function(da.res, group.by = NULL, alpha = 0.05, subset.nhoods = NULL) {
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(da.res)) {
      stop(group.by, " is not a column in da.res. Have you forgot to run annotateNhoods(x, da.res, ", 
           group.by, ")?")
    }
    if (is.numeric(da.res[, group.by])) {
    }
    da.res <- mutate(da.res, group_by = da.res[, group.by])
  }
  else {
    da.res <- mutate(da.res, group_by = "g1")
  }
  if (!is.factor(da.res[, "group_by"])) {
    message("Converting group_by to factor...")
    da.res <- mutate(da.res, group_by = factor(group_by, 
                                               levels = unique(group_by)))
  }
  if (!is.null(subset.nhoods)) {
    da.res <- da.res[subset.nhoods, ]
  }
  beeswarm_pos <- ggplot_build(da.res %>% mutate(is_signif = ifelse(SpatialFDR < 
                                                                      alpha, 1, 0)) %>% arrange(group_by) %>% ggplot(aes(group_by, 
                                                                                                                         logFC)) + geom_quasirandom(  bandwidth = 0.1))
  pos_x <- beeswarm_pos$data[[1]]$x
  pos_y <- beeswarm_pos$data[[1]]$y
  n_groups <- levels(da.res$group_by) %>% length()
  da.res_mutated <- da.res %>% mutate(is_signif = ifelse(SpatialFDR < alpha, 
                                       1, 0)) %>% mutate(logFC_color = ifelse(is_signif == 1, 
                                                                              logFC, NA)) %>% arrange(group_by) %>% mutate(Nhood = factor(Nhood, 
                                                                                                                                          levels = unique(Nhood))) %>% mutate(pos_x = pos_x, pos_y = pos_y)
    
    
    ggplot(da.res_mutated,aes(pos_x, pos_y, color = logFC_color)) + scale_color_gradient(low="#4b8aa4", high="#c15557", na.value = "lightgrey") + 
    guides(color = "none") + xlab(group.by) + ylab("Log Fold Change") + 
    scale_x_continuous(breaks = seq(1, n_groups), labels = setNames(levels(da.res$group_by), 
                                                                    seq(1, n_groups))) + geom_point(size=0.7) + coord_flip() + 
    theme + theme(strip.text.y = element_text(angle = 0),aspect.ratio=1.5)
    
    
    da.res$direction <- ifelse(da.res$logFC < 0, "down", "up")
    da.res[da.res$SpatialFDR > 0.05,]$logFC <- NA
    
    ggplot(da.res_mutated,aes(x=celltype, y=logFC, color=direction))+  geom_quasirandom()+scale_color_gradient(low="#4b8aa4", high="#c15557", na.value = "lightgrey")
    
}


#reoreder cells
reorder_cells <- function(df, celltype_df= celltype_l1, order=order_cells, reverse=F, neworder=F){
  if(neworder == T){
    order <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/new_order_cells_dual.rds"))
  }
  df$celltype <- as.character(df$celltype)
  df$celltype_l1 <- unlist(lapply(df$celltype, function(celltype) celltype_df[celltype_df$cell_type==celltype, "predicted.celltype.l1"]))
  df$celltype_l1 <- factor(df$celltype_l1, levels =order$predicted.celltype.l1[order$predicted.celltype.l1 %in% unique(df$celltype_l1)])
  if (reverse){df$celltype <- factor(df$celltype, levels =rev(order$cell_type[order$cell_type %in% unique(df$celltype)]))}
  else{  df$celltype <- factor(df$celltype, levels =order$cell_type[order$cell_type %in% unique(df$celltype)])}
  return(df)
}


matrix_counts_gene <- function( cell_level, cell_type,gene, pheno="Age"){
  print(cell_type)
  expression <- readRDS(paste0(basepath, "/Data/scRNAseq/Yazar2022/sce_data_objects/", cell_type, "_", cell_level, "_pseudobulk_counts.rds"))
  expr <- expression[gene,]
  if (cell_type == "NK_CD56bright"){  expr_gene <- as.data.frame(t(as.matrix(expr@assays@data[[cell_type]])))}
  else{  expr_gene <- as.data.frame(t(as.matrix(expr@assays@data[[gsub("_", " ", cell_type)]])))}
  expr_gene$Age_cat_all <-ifelse(expr$Age < 40, "Y", ifelse( expr$Age > 60, "O", "M"))
  expr_gene$Age_cat <- ifelse(expr$Age < 40, "Y", "O")
  expr_gene$Age <- expr$Age
  expr_gene$Sex <- expr$Gender
  if (cell_type == "NK_CD56bright"){expr_gene$celltype <- cell_type}
  else{expr_gene$celltype <- gsub("_", " ", cell_type)} 
  return(expr_gene)
}
perform_counts_dge <- function(exprs.data, test.model, gene.offset=gene.offset,
                                model.contrasts=NULL, n.coef=NULL){
  
  i.dge <- DGEList(counts=exprs.data,
                   lib.size=colSums(exprs.data))
  
  if(isTRUE(gene.offset)){
    n.gene <- apply(exprs.data, 2, function(X) sum(X > 0))
      test.model <- cbind(test.model, n.gene)
      colnames(test.model) <- c(colnames(test.model)[seq_len(ncol(test.model)-1)], "NGenes")
    }
  
  i.dge <- estimateDisp(i.dge, test.model)
  i.fit <- glmQLFit(i.dge, test.model, robust=TRUE)
  
  if(!is.null(model.contrasts)){
    mod.constrast <- makeContrasts(contrasts=model.contrasts, levels=test.model)
    i.res <- as.data.frame(topTags(glmQLFTest(i.fit, contrast=mod.constrast),
                                   sort.by='none', n=Inf))
  } else{
    if(is.null(n.coef)){
      n.coef <- ncol(test.model)
    }
    i.res <- as.data.frame(topTags(glmQLFTest(i.fit, coef=n.coef), sort.by='none', n=Inf))
  }
  return(i.res)
}

# total_n_gens <- function(pheno, cell_level){
#   DVG <- readRDS(paste0(data_path, "msopena//01_meta-analysis_Age_Sex/DEG_", pheno, "_", cell_level, ".rds"))
#   DVG <- DVG[DVG$dataset == "meta",]
#   count_degs <- DVG %>% group_by(celltype) %>% count()
#   count_degs$celltype <- gsub( "_", " ", count_degs$celltype)
#   colnames(count_degs) <- gsub( "n", "ntotal", colnames(count_degs))
#   return(count_degs)
# }
# 
# 
# 
# #extract significan DEGs - pseudobulk

extract_signif_DEG_pseudobulk <- function(pheno, cell_level, method="holm"){
  #read datasets for each pheno and cell_level
  DVG <- readRDS(paste0(data_path, "msopena//01_meta-analysis_Age_Sex/robjects/02_DEA/DEG_pseudobulk_", pheno, "_", cell_level, ".rds"))
  #modify some names for plotting
  DVG$celltype <- gsub("_", " ", DVG$celltype)
  DVG$dataset <- gsub("v2", "data 1", DVG$dataset )
  DVG$dataset <- gsub("v3", "data 2", DVG$dataset )
  DVG$dataset <- gsub("pilot3", "data 3", DVG$dataset )
  #Select significant genes
  if(method == "fdr"){
    DVG_signif <- DVG[DVG$fdr < 0.05,]
  }else{
    DVG_signif <- DVG[DVG$holm < 0.05,]
  }
  DVG_signif$direction <- ifelse(DVG_signif$logFC < 0, "down", "up")
  return(DVG_signif)
}


plot_example_nhood <- function( gene, cd8tem=F, cells=NULL){
  expr_gene <- expr[gene,] %>% as.data.frame()
  colnames(expr_gene) <- "Expression"
  expr_gene$cell_id <- rownames(expr_gene)  
  expr_gene <- expr_gene %>% left_join(cell_nhood_df, by="cell_id")
  expr_nhood <- expr_gene %>% tidyr::drop_na() %>% dplyr::rename(celltype =cell_type_nhood) %>% group_by(Nhood, celltype, color ) %>% 
    dplyr::summarize(across(Expression, mean, na.rm = TRUE)) 
  expr_nhood <- expr_nhood[expr_nhood$color != "ns", ]
  df <- reorder_cells(expr_nhood, neworder = T)
  if(cd8tem == F){
    return (ggplot(df, aes(x=celltype, y=Expression))+geom_boxplot(aes(fill=color), outlier.shape = NA)+theme+
              scale_fill_manual(values=c("enriched"=red, "depleted"=blue, "ns"="grey"), )+xlab("")+labs(fill="Nhoods")+
              ggtitle(paste0(gene))+theme(plot.title = element_text(size = 20, face = "bold.italic", hjust = 0.5), axis.text.x=element_text(size=10, angle=90,  hjust = 0.95, vjust = 0.6)))
  }else{
    ggplot(df[df$celltype %in% cells, ], aes(x=celltype, y=Expression))+geom_boxplot(aes(fill=color), outlier.shape = NA)+theme+
      scale_fill_manual(values=c("enriched"=red, "depleted"=blue, "ns"="grey"), )+xlab("")+labs(fill="Nhoods")+
      ggtitle(paste0(gene))+theme(plot.title = element_text(size = 20, face = "bold.italic", hjust = 0.5),  axis.text.x=element_text(size=10, angle=90,  hjust = 0.95, vjust = 0.6))
    
  }
}


# total_n_gens_pseudo <- function(pheno, cell_level){
#   DVG <- readRDS(paste0(data_path, "msopena//01_meta-analysis_Age_Sex/DEG_pseudobulk_", pheno, "_", cell_level, ".rds"))
#   DVG <- DVG[DVG$dataset == "meta",]
#   count_degs <- DVG %>% group_by(celltype) %>% count()
#   count_degs$celltype <- gsub( "_", " ", count_degs$celltype)
#   colnames(count_degs) <- gsub( "n", "ntotal", colnames(count_degs))
#   return(count_degs)
# }
# 
# 
# 
# 
# 
# #extract significant DEGs - reanalysis 
# 
# extract_signif_DEG_reanalysis <- function(reanalysis, cell_level){
#   comb_S <- reanalysis[reanalysis$fdr < 0.5, ]
#   comb_S$gene <- comb_S$primerid
#   colnames(comb_S) <- gsub("_", "", colnames(comb_S))
#   comb_S$celltype <- gsub("_", " ", comb_S$celltype)
#   comb_S$direction <- gsub("up", "increase", comb_S$direction)
#   comb_S$direction <- gsub("down", "decrease", comb_S$direction)
#   return(comb_S)
# }
# 

# extract_signif_DVG <- function(pheno, cell_level, adj="holm"){
#   #read datasets for each pheno and cell_level 
#   DVG <- readRDS(paste0(data_path, "msopena/01_meta-analysis_Age_Sex/DVG_", pheno, "_", cell_level, ".rds"))
#   #modify some names for plotting 
#   DVG$celltype <- gsub("_", " ", DVG$celltype)
#   DVG$dataset <- gsub("v2", "data 1", DVG$dataset )
#   DVG$dataset <- gsub("v3", "data 2", DVG$dataset )
#   DVG$dataset <- gsub("pilot3", "data 3", DVG$dataset )
#   DVG$direction <- ifelse(DVG$CCV_diff > 0, "increase", "decrease")
#   #Select significant genes 
#   if(adj=="fdr"){
#     DVG_signif <- DVG[DVG$fdr < 0.05 ,]
#   }else{
#     DVG_signif <- DVG[DVG$holm < 0.05 ,]
#   }
#   
#   #DVG_signif <- DVG_signif[!is.na(DVG_signif$celltype),]
#   #DVG_signif <- inner_join( DVG_signif,ncells, by=c('celltype'='celltype', 'dataset'='dataset'))
#   return(DVG_signif)
# }
# 
# total_n_gens_DVG <- function(pheno, cell_level){
#   DVG <- readRDS(paste0(data_path, "msopena/01_meta-analysis_Age_Sex/DVG_", pheno, "_", cell_level, ".rds"))
#   DVG <- DVG[DVG$dataset == "meta",]
#   count_degs <- DVG %>% group_by(celltype) %>% count()
#   count_degs$celltype <- gsub( "_", " ", count_degs$celltype)
#   colnames(count_degs) <- gsub( "n", "ntotal", colnames(count_degs))
#   return(count_degs)
# }
# 
# extract_signif_genes<- function(DVG, cell_level){
#   #read datasets for each pheno and cell_level 
#   ncells <- readRDS(paste0(data_path, "msopena/01_meta-analysis_Age_Sex/NCells_", cell_level, ".rds"))
#   colnames(ncells)[3]<- "ncells"
#   #modify some names for plotting 
#   DVG$celltype <- gsub("_", " ", DVG$celltype)
#   DVG$dataset <- gsub("v2", "data 1", DVG$dataset )
#   DVG$dataset <- gsub("v3", "data 2", DVG$dataset )
#   DVG$dataset <- gsub("pilot3", "data 3", DVG$dataset )
#   #Select significant genes 
#   DVG_signif <- DVG[DVG$fdr < 0.05,]
#   DVG_signif <- inner_join( DVG_signif,ncells, by=c('celltype'='celltype', 'dataset'='dataset'))
#   return(DVG_signif)
# }
# 
# 
# #################### OLD
# meta_vs_reanalysis_pseudo <- function(cell_level,fc_age_rean=0, fc_sex_rean=0, fc_age_meta= 0, fc_sex_meta=0){
#   degs_age <- readRDS(paste0(data_path, "/msopena/01_meta-analysis_Age_Sex/robjects/02_DEA/pseudobulk_Age_", cell_level, ".rds"))
#   meta_age <- degs_age[, c("celltype", "fdr", "logFC")] %>% filter(fdr < 0.05, abs(logFC) > fc_age_meta)
#   meta_age$pheno <- "Age"
#   meta_age$analysis <- "meta"
#   degs_sex <- readRDS(paste0(data_path, "/msopena/01_meta-analysis_Age_Sex/robjects/02_DEA/pseudobulk_Gender_", cell_level, ".rds"))
#   meta_sex <- degs_sex[, c("celltype", "fdr", "logFC")]  %>% filter(fdr < 0.05, abs(logFC) > fc_sex_meta)
#   meta_sex$pheno <- "Sex"
#   meta_sex$analysis <- "meta"
#   reanalysis_sex <- readRDS(paste0(data_path, "/msopena/01_meta-analysis_Age_Sex/robjects/02_DEA/pseudobulk_reanalysis_Gender_", cell_level, ".rds"))
#   colnames(reanalysis_sex) <- gsub("adj.P.Val", "fdr",   colnames(reanalysis_sex))
#   rean_sex <- reanalysis_sex[, c("celltype", "fdr", "logFC")]  %>% filter(fdr < 0.05, abs(logFC) > fc_sex_rean)
#   rean_sex$pheno <- "Sex"
#   rean_sex$analysis <- "reanalysis"
#   reanalysis_age <- readRDS(paste0(data_path, "/msopena/01_meta-analysis_Age_Sex/robjects/02_DEA/pseudobulk_reanalysis_Age_", cell_level, ".rds"))
#   colnames(reanalysis_age) <- gsub("adj.P.Val", "fdr",   colnames(reanalysis_age))
#   rean_age <- reanalysis_age[, c("celltype", "fdr", "logFC")] %>% filter(fdr < 0.05, abs(logFC) > fc_age_rean)
#   rean_age$pheno <- "Age"
#   rean_age$analysis <- "reanalysis"
#   all <- rbind(meta_age,meta_sex, rean_sex, rean_age )
#   all<-   all %>% group_by(pheno, analysis, celltype) %>% count()
#   all$celltype <- gsub("_", " ", all$celltype)
#   all$color <- paste0(all$analysis, "_", all$pheno)
#   order <- order_cells[[cell_level]]
#   order <- gsub("_", " ", order)
#   all$cell_type <- factor(all$celltype, levels = order)
#   pal <- c("meta_Sex"="#46644df7",     "meta_Age"= "#76404A"  ,   "reanalysis_Sex"="#5c896662", "reanalysis_Age"="#BFA5A8")
#   all$pheno <- factor(all$pheno, levels= c("Sex", "Age"))
#   p <- ggplot(all, aes(x=cell_type, y=n, fill=pheno, alpha=analysis))+geom_col(position= position_dodge(preserve="single"))+ theme+ylab("Number of DEGs") +xlab("")+scale_fill_manual(values=c("Sex"="#46644df7", "Age"= "#76404A" ))+
#     facet_wrap(~ pheno, ncol = 1, scales="free_y", strip.position = "right") + theme(strip.background = element_rect(fill=NA, linewidth = 0), strip.text = element_blank())+ 
#     theme(legend.position = "top")+scale_alpha_discrete(range = c(0.3, 1)) +theme(axis.text.x = element_text(angle=90, hjust = 0.95, vjust = 0.6))
#   return(p)}
# 
# 
