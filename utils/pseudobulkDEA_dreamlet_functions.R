#!/usr/bin/env Rscript

############################### Load R packages ################################
print('Loading R packages...')
shhh <- suppressPackageStartupMessages
shhh(library(Seurat))
shhh(library(MAST))
shhh(library(SingleCellExperiment))
shhh(library(dreamlet))
shhh(library(zenith))
shhh(library(scater))
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(ggplot2))
shhh(library(RColorBrewer))

######################## Functions used in pseudobulkDEA_limmadream.R ##################
# 1. Get contrasts
get_coefName <- function(phe, sce){
  md <- colData(sce)
  contrast_var <- phe
  if(is.factor(md[[phe]])){
    contrast_var.levels <- levels(md[[phe]])
    if(length(contrast_var.levels)>1){
      contrast_var <- paste0(phe, contrast_var.levels[[2]])
    }else{
      print(paste0(phe, ' only has 1 level: ', contrast_var.levels, '. NOT consider this phenotype.'))
      contrast_var <- NULL
    }
  }
  return(contrast_var)
}

# 2. Define formula (VP or DEA)
define_form <- function(gt, df, vp){
  # forms
  print(gt)
  fixed_var_dea <- df[df$type=='fixed',]$covariate
  random_var_dea <- df[df$type=='random',]$covariate
  if(vp){
    if(gt=='VP'){
      print('Not considering the batch effect in the VariancePartition...')
      random_var_dea.idx <- which(df$covariate==random_var_dea)
      df <- df[-random_var_dea.idx,]
    }
  }
  model_vars <- df$covariate
  
  if(gt=='VP'){
    df_to_mod <- df[df$type=='fixed' & df$class=='factor',]
    if(nrow(df_to_mod)>0){
      df[df$type=='fixed' & df$class=='factor',]$type <- 'random'
      }
  }
  
  fixed_var <- df[df$type=='fixed',]$covariate
  fixed.fmla <- paste(fixed_var,collapse='+')
  random_var <- df[df$type=='random',]$covariate
  random.fmla <- NULL
  if(length(random_var)>0){
    random.fmla <- paste(paste0('(1|',random_var,')'),collapse='+')
  }
  form_vars <- paste(c(fixed.fmla,random.fmla), collapse='+')
  form_vars <- paste0('~',form_vars)
  print(paste0('Fitting lmer: ',form_vars))
  form <- as.formula(form_vars)
  
  # specificy colors
  sex.hex <- brewer.pal(9, 'Greens')[7]
  age.hex <- brewer.pal(9, 'Blues')[7]
  random.hex <- brewer.pal(9, 'Greys')[7]
  residuals.hex <- brewer.pal(9, 'Greys')[3]
  cols_vars <- c(sex.hex, age.hex, random.hex, residuals.hex)
  names(cols_vars) <- c(fixed_var_dea[c(1,2)], random_var_dea[1], 'Residuals')
  if(length(random_var_dea)>1 | length(fixed_var_dea)>1){
    if(length(random_var_dea)>1){
      random_added.hex <- brewer.pal(9, 'Greys')[5]
      names(random_added.hex) <- random_var_dea[2]
      cols_vars <- c(cols_vars, random_added.hex)
    }
    if(length(fixed_var_dea)>1){
      fixed_added.hex <- brewer.pal(9, 'Oranges')[5]
      names(fixed_added.hex) <- fixed_var_dea[3]
      cols_vars <- c(cols_vars, fixed_added.hex)
    }
  }
  model_vars_in <- c(model_vars, 'Residuals')
  cols_vars <- cols_vars[names(cols_vars)%in%model_vars_in]
  cols_vars <- cols_vars[match(model_vars_in, names(cols_vars))]
  
  # output
  out <- list(form = form,
              cols = cols_vars)
  
  return(out)
}

# 3. DEA extract and plots
extract_plots <- function(i, dea_res, contrast_var, vp_res, cols, o_dir){
  # Extract results (topTable --> DEGs)
  ### Each entry in res.dl stores a model fit by dream(), and results can be extracted using topTable() as in limma by specifying the coefficient of interest. 
  ### The results shows the gene name, log fold change, average expression, t-statistic, p-value, FDR (i.e. adj.P.Val).
  genes <- rownames(dea_res[[1]]$residuals)
  topTable.res <- topTable(dea_res, 
                           coef = contrast_var,
                           number = length(genes))
  degs <- topTable.res[topTable.res$adj.P.Val<=0.05,]$ID
  
  ### DEA plots (for all genes) ###
  # Volcano plots
  ## The volcano plot can indicate the strength of the differential expression signal with each cell type. Red points indicate FDR < 0.05.
  plotVolcano.p <- plotVolcano(dea_res, coef = contrast_var)
  plotVolcano.fn <- paste0(o_dir, 'plotVolcano.png')
  print(paste0('Saving plotVolcano in: ', plotVolcano.fn))
  ggsave(plotVolcano.fn, plotVolcano.p)
  
  ### DEA and VP plots (only if DEGs) ###
  if(length(degs)>0){
    ### DEA plots ###
    # Gene-level heatmap
    ## For each cell type and specified gene, show z-statistic from dreamlet analysis. 
    ## Grey indicates that insufficient reads were observed to include the gene in the analysis.
    plotGeneHeatmap.p <- plotGeneHeatmap(dea_res, coef=contrast_var, genes=degs)
    plotGeneHeatmap.fn <- paste0(o_dir, 'plotGeneHeatmap.png')
    print(paste0('Saving plotGeneHeatmap in: ', plotGeneHeatmap.fn))
    ggsave(plotGeneHeatmap.fn, plotGeneHeatmap.p)
    
    ## Forest plot
    ## A forest plot shows the log fold change and standard error of a given gene across all cell types. The color indicates the FDR.
    # os_dir <- paste0(o_dir, '/plotForest/')
    # if(!dir.exists(os_dir)){dir.create(os_dir, recursive = T)}
    # plotForest.save <- lapply(degs, function(i){
    #   plotForest.p <- plotForest(dea_res, coef = contrast_var, gene = i)
    #   plotForest.fn <- paste0(os_dir, i, '.png')
    #   print(paste0('Saving plotForest in: ', plotForest.fn))
    #   ggsave(plotForest.fn, plotForest.p)
    #   return(NULL)
    # })
    
    ### VP plots ###
    # Pick only DEGs in the VP results
    vp_res.degs <- vp_res[vp_res$gene%in%degs,]
    cnames <- colnames(vp_res.degs)
    vp_res.degs <- vp_res.degs[,match(cnames, colnames(vp_res.degs))]
    vp_res.degs <- vp_res.degs[match(degs, vp_res.degs$gene),]
    
    # plotPercentBars --> some genes show differences when estimating only Gender/Age or Gender/Age/date
    plotPercentBars.p <- plotPercentBars(vp_res.degs, cols)
    plotPercentBars.fn <- paste0(o_dir, 'plotPercentBars.png')
    print(paste0('Saving plotPercentBars in: ', plotPercentBars.fn))
    ggsave(plotPercentBars.fn, plotPercentBars.p)
    
    # plotVarPart
    plotVarPart.p <- plotVarPart(vp_res.degs, cols, label.angle=60) 
    plotVarPart.fn <- paste0(o_dir, 'plotVarPart.png')
    print(paste0('Saving plotVarPart in: ', plotVarPart.fn))
    ggsave(plotVarPart.fn, plotVarPart.p)
  }
  cat('\n')
  res <- topTable.res
  return(res)
}

# 4. DEA + VP extract and plots (by phenotype)
extract_plots_by_phe <- function(phe, dea_res, vp_res, c_list, cols, o_dir){
  print(phe)
  
  # create output dir
  out_sdir <- paste0(o_dir, '/', phe, '/')
  if(!dir.exists(out_sdir)){dir.create(out_sdir, recursive = T)}
  
  # pick contrast variable
  contrast_coefName <- c_list[[phe]]
  
  # extract and plots
  res <- extract_plots(i = phe,
                       dea_res = dea_res,
                       contrast_var = contrast_coefName,
                       vp_res = vp_res,
                       cols = cols,
                       o_dir = out_sdir)
  return(res)
}

# 5. dreamlet
dreamlet.func <- function(ge_dge, covariates, min_prop, contrast_list, vp_reduced, out_dir, gene_test = c('VP','DEA')){
  ### Defining the VP/DEA formulas ###
  print('Defining the VP/DEA formulas...')
  gene_test.forms <- sapply(gene_test, function(i) define_form(i, covariates, vp_reduced), simplify = FALSE)
  
  #### Normalize and apply voom/voomWithDreamWeights ####
  # Run processAssays()
  form <- gene_test.forms$DEA$form
  print('Normalizing the pseudobulk-data...')
  system.time(res.proc <- processAssays(ge_dge, form, 
                                        min.cells = 5, 
                                        min.count = 5, 
                                        min.samples = 4, 
                                        min.prop = min_prop))
  
  # View details of dropping samples
  details(res.proc)
  
  # Check nSamples and nGenes tested
  genes_all <- rownames(ge_dge)
  genes_tested <- rownames(as.data.frame(res.proc))
  genes_all.n <- nrow(ge_dge)
  genes_tested.n <- nrow(as.data.frame(res.proc))
  genes_tested.prop <- round(genes_tested.n/genes_all.n,3)
  samples_all <- colnames(ge_dge)
  samples_tested <- colnames(as.data.frame(res.proc))
  samples_all.n <- ncol(ge_dge)
  samples_tested.n <- ncol(as.data.frame(res.proc))
  samples_tested.prop <- round(samples_tested.n/samples_all.n,3)
  print(paste0('# Genes tested: ', genes_tested.n, ', out of ', genes_all.n, ' (', genes_tested.prop, ')'))
  print(paste0('# Samples tested: ', samples_tested.n, ', out of ', samples_all.n, ' (', samples_tested.prop, ')'))
  
  # Show voom plot for each cell clusters
  ## Here the mean-variance trend from voom is shown for each cell type. Cell types with sufficient number of cells and reads show a clear mean-variance trend. While in rare cell types like megakaryocytes, fewer genes have sufficient reads and the trend is less apparent.
  plotVoom.p <- plotVoom(res.proc)
  plotVoom.fn <- paste0(out_dir, 'plotVoom.png')
  ggsave(plotVoom.fn, plotVoom.p)
  
  ### Differential expression ###
  ## Since the normalized expression data and metadata are stored within res.proc, only the regression formula remains to be specified.
  ## Here we only included the stimulus status, but analyses of larger datasets can include covariates and random effects.
  ## With formula ~ StimStatus, an intercept is fit and coefficient StimStatusstim log fold change between simulated and controls.
  ## Differential expression analysis within each assay, evaluated on the voom normalized data
  print('Running DEA...')
  system.time(res.dl <- dreamlet(res.proc, form))
  
  ### Variance partitioning ###
  ## The variancePartition package uses linear and linear mixed models to quanify the contribution of multiple sources of expression variation at the gene-level.
  ## For each gene it fits a linear (mixed) model and evalutes the fraction of expression variation explained by each variable.
  ## Variance fractions can be visualized at the gene-level for each cell type using a bar plot, or genome-wide using a violin plot.
  # Mymic DEA model: https://github.com/GabrielHoffman/dreamlet/issues/4#issuecomment-1507767030
  # vp.lst = fitVarPart(res.proc, form) # not working --> 2 alternatives (use option 1):
  #### 1. Use ~(1|Sex)+Age+(1|Batch). variancePartition works best when categorical variables are modeled as a random effects. It't not an issue that this formula isn't identical the the differential expression formula.
  #### 2. We can try to regress out the Batch variable, and do fitVarPart() on the residuals --> You could do that, but you'd have to use variancePartition::fitExtractVarPartModel() directly.
  print('Running VariancePartition...')
  system.time(vp.lst <- fitVarPart(res.proc, gene_test.forms$VP$form))
  cols_vars <- gene_test.forms$VP$cols
  
  ### Extract results and plots (DEA and VP) ###
  extract_plots_by_phe.res <- sapply(names(contrast_list),
                                     function(i) extract_plots_by_phe(phe = i,
                                                                      dea_res = res.dl,
                                                                      vp_res = vp.lst,
                                                                      c_list = contrast_list,
                                                                      cols = cols_vars,
                                                                      o_dir = out_dir), simplify = FALSE)
  
  ### Save outputs ###
  res <- list(processed = res.proc,
              dea = res.dl,
              vp = vp.lst,
              topTable = extract_plots_by_phe.res)
  res_fn <- paste0(out_dir, 'dea_vp_topTable.rds')
  print(paste0('Saving dreamlet results: ', res_fn))
  saveRDS(res, res_fn)
  
  return(res)
}

