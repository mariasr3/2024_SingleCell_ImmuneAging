# Load necessary libraries
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

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
#dir.create(path_enrichScore, recursive = TRUE)

#functions
source(paste0(data_path, "/msopena/02_OneK1K_Age/scripts/functions.R"))


# Step 1: Read the pre-existing Disease-Gene Dictionary from RDS File
read_disease_gene_dict <- function(rds_file) {
  disease_gene_dict <- readRDS(rds_file)
  return(disease_gene_dict)
}

# Step 2: Separate Files Based on Cell Type
split_by_celltype <- function(input_file, sex, output_dir) {
  data <- readRDS(input_file)
  cell_types <- unique(data$celltype)
  
  for (cell_type2 in cell_types) {
    cell_data <- data %>% filter(celltype == cell_type2)
    output_file <- file.path(output_dir, paste0(sex, "_", cell_type2, ".rds"))
    saveRDS(cell_data, output_file)
  }
}

# Step 3: Perform Fisher's Exact Test, Print, and Collect Significant Results
perform_fisher_test_and_collect <- function(cell_file, disease_gene_dict, sex_spec_degs=F, sex=NA) {
  data <- readRDS(cell_file)
  results <- data.frame()
  
  for (trait in names(disease_gene_dict)) {
    if(sex_spec_degs == T){
      cell <-  str_extract(cell_file, "(?<=M_|F_)[^\\.rds]+")
      sex_spec <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/DEG_genes_only_", sex, ".rds"))
      sex_spec_genes <- sex_spec %>% filter( celltype == cell) %>% pull(gene)
      de_genes <- data %>% filter(gene %in% sex_spec_genes) %>% pull(gene)
      non_de_genes <- data %>% filter(!gene %in% de_genes) %>% pull(gene)
    }else{
      de_genes <- data %>% filter(adj.P.Val < 0.05) %>% pull(gene)
      non_de_genes <- data %>% filter(adj.P.Val >= 0.05) %>% pull(gene)
    }

    in_gwas_de <- sum(de_genes %in% disease_gene_dict[[trait]])
    not_in_gwas_de <- sum(!de_genes %in% disease_gene_dict[[trait]])
    in_gwas_non_de <- sum(non_de_genes %in% disease_gene_dict[[trait]])
    not_in_gwas_non_de <- sum(!non_de_genes %in% disease_gene_dict[[trait]])
    
    fisher_matrix <- matrix(c(in_gwas_de, not_in_gwas_de, in_gwas_non_de, not_in_gwas_non_de),
                            nrow = 2, byrow = TRUE)
    
    fisher_test <- fisher.test(fisher_matrix)
    
    results <- rbind(results, data.frame(
      CellType = sub("^[FM]_", "", basename(cell_file)),  # Remove 'F_' or 'M_' prefix
      Trait = trait,
      PValue = fisher_test$p.value,
      OR = fisher_test$estimate
    ))
  }
  
  # Adjust p-values using FDR
  results <- results %>%
    mutate(AdjustedPValue = p.adjust(PValue, method = "fdr"))
  
  return(results)
}

# Example usage:

# Set file paths
disease_gene_dict_file <- paste0(data_path, "msopena/02_OneK1K_Age/robjects/15_OverlapGWAS/ordered_gene_dictionary.rds")  # RDS file containing the disease-gene dictionary
all_f_file <- paste0(data_path, "msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/AllCells_F_deaTopTable.rds")  # RDS file for females
all_m_file <- paste0(data_path, "msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/AllCells_M_deaTopTable.rds")   # RDS file for males
output_dir_f <- paste0(data_path, "msopena/02_OneK1K_Age/robjects/15_OverlapGWAS/F/")
output_dir_m <- paste0(data_path, "msopena/02_OneK1K_Age/robjects/15_OverlapGWAS/M/")
fisher_output_file <- paste0(data_path, "msopena/02_OneK1K_Age/robjects/15_OverlapGWAS/FisherResults/")  # Output file for all significant results

# Ensure that the output directory exists
output_dir <- dirname(fisher_output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read the pre-existing disease-gene dictionary from the RDS file
disease_gene_dict <- read_disease_gene_dict(disease_gene_dict_file)

# Split files by cell type for both females and males
split_by_celltype(all_f_file, "F", output_dir_f)
split_by_celltype(all_m_file, "M", output_dir_m)

# a. Perform fisher test with all age-DEGs --
# Collect significant results for females
cell_files_f <- list.files(output_dir_f, full.names = TRUE)
results_f <- bind_rows(lapply(cell_files_f, perform_fisher_test_and_collect, disease_gene_dict))
results_f$sex <- "Females"
# Collect significant results for males
cell_files_m <- list.files(output_dir_m, full.names = TRUE)
results_m <- bind_rows(lapply(cell_files_m, perform_fisher_test_and_collect, disease_gene_dict))
results_m$sex <- "Males"
# Combine all significant results
all_results <- bind_rows(results_f, results_m)
# Save all significant results to a single CSV file
saveRDS(all_results, paste0(fisher_output_file, "FisherResults_AllDEGs.rds"))

# b. Perform Fisher test with sex-specific age_DEGs --
cell_files_f <- list.files(output_dir_f, full.names = TRUE)
results_f <- bind_rows(lapply(cell_files_f, perform_fisher_test_and_collect, disease_gene_dict, T, "F"))
results_f$sex <- "Females"
# Collect significant results for males
cell_files_m <- list.files(output_dir_m, full.names = TRUE)
results_m <- bind_rows(lapply(cell_files_m, perform_fisher_test_and_collect, disease_gene_dict,T, "M"))
results_m$sex <- "Males"
# Combine all significant results
all_results <- bind_rows(results_f, results_m)
# Save all significant results to a single CSV file
saveRDS(all_results, paste0(fisher_output_file, "FisherResults_SexSpecificDEGs.rds"))
# The `all_significant_results` variable now contains all significant findings across all cell types and sexes.
