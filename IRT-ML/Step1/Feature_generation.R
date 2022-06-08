# install.packages("h2o", type="source", repos=(c("http://h2o-release.s3.amazonaws.com/h2o/latest_stable_R")))
# Library
# packages
library(biomaRt)
library(RMariaDB)
library(easier)
library(tidyverse)
library(tximport)
library(tximportData)
library(purrr)

feature <- function(cancer_type){
  con_RAW_COUNT <- DBI::dbConnect(drv = MariaDB(), host = "192.168.0.91", port = 3306, user = "root", password = "sempre813!", 
                                  dbname = "TCGA_Gene_expression_IlluminaHiSeq_RNASeqV2_count")
  con_TPM <- DBI::dbConnect(drv = MariaDB(), host = "192.168.0.91", port = 3306, user = "root", password = "sempre813!", 
                            dbname = "TCGA_Gene_expression_IlluminaHiSeq_RNASeqV2_TPM")
  
  if(cancer_type == "PANCAN"){
    all_list <- intersect(con_RAW_COUNT %>% dbListTables(), con_TPM %>% dbListTables())
    pan_list_count <- lapply(X = all_list, FUN = function(value){
      tbl(con_RAW_COUNT, value) %>% collect() %>% 
        dplyr::select(Gene, contains("-01A")) %>% 
        return()
    })
    pan_list_tpm <- lapply(X = all_list, FUN = function(value){
      tbl(con_TPM, value) %>% collect() %>%
        dplyr::select(Gene, contains("-01A")) %>% 
        return()
    })
    
    # tpm <- pan_list_tpm %>% purrr::reduce(., inner_join, by = "Gene")
    tpm <- pan_list_tpm %>% bind_cols()
    gene_name <- pan_list_tpm[[1]] %>% dplyr::select(Gene)
    tpm <- tpm %>% select_if(is.numeric) %>% 
      bind_cols(gene_name, .) %>% 
      distinct(Gene, .keep_all = TRUE)
    
    tpm_m <- tpm %>% dplyr::select(-Gene)
    rownames(tpm_m) <- tpm$Gene
    
    counts <- pan_list_count %>% bind_cols()
    gene_name <- pan_list_count[[1]] %>% dplyr::select(Gene)
    counts <- counts %>% select_if(is.numeric) %>% 
      bind_cols(gene_name, .) %>% 
      distinct(Gene, .keep_all = TRUE)
    
    counts_m <- counts %>% dplyr::select(-Gene) %>% as.matrix()
    rownames(counts_m) <- counts$Gene
    
  } else {
    tpm <- tbl(con_TPM, cancer_type) %>% collect() %>% 
      distinct(Gene, .keep_all = TRUE)
    
    tpm_m <- tpm %>% dplyr::select(-Gene) 
    rownames(tpm_m) <- tpm$Gene
    
    counts <- tbl(con_RAW_COUNT, cancer_type) %>% collect() %>% 
      distinct(Gene, .keep_all = TRUE)
    counts_m <- counts %>% dplyr::select(-Gene) %>% as.matrix()
    rownames(counts_m) <- counts$Gene
  }
  
  # ----------------------------------------- #
  # compute system-based derived signatures
  # ----------------------------------------- #
  
  # Computation of cell fractions
  cell_fractions <- compute_cell_fractions(RNA.tpm=tpm_m)
  
  # normalized
  scale_model <- cell_fractions %>% as_tibble() %>% preProcess(method = "range")
  cell_fractions_scale <- predict(scale_model, cell_fractions) %>% as_tibble() %>% 
    rename_with(~paste0("Cell_Fractions.", .))
  
  # Computation of pathway activity
  pathways_activity <- compute_pathways_scores(RNA.counts=counts_m, remove.genes.ICB_proxies=TRUE)
  
  # normalized
  scale_model <- pathways_activity[[1]] %>% as_tibble() %>% preProcess(method = "range")
  pathways_activity_scale <- predict(scale_model, pathways_activity[[1]]) %>% as_tibble() %>% 
    rename_with(~paste0("Pathways_activity.", .))
  
  # Computation of TF activity
  TF_activity <- compute_TF_activity(RNA.tpm=tpm_m, remove.genes.ICB_proxies=FALSE)
  
  # normalized
  scale_model <- TF_activity[[1]] %>% as_tibble() %>% preProcess(method = "range")
  TF_activity_scale <- predict(scale_model, TF_activity[[1]]) %>% as_tibble() %>% 
    rename_with(~paste0("TF_activity.", .))
  
  # Computation of LR pairs weights
  lrpairs_weights <- compute_LR_pairs(RNA.tpm=tpm_m, remove.genes.ICB_proxies=FALSE, compute.cytokines.pairs=FALSE, cancertype="pancan") 
  
  # normalized
  scale_model <- lrpairs_weights[[1]] %>% as_tibble() %>% preProcess(method = "range")
  lrpairs_weights_scale <- predict(scale_model, lrpairs_weights[[1]]) %>% as_tibble() %>% 
    rename_with(~paste0("LR_Pairs_weights.", .))
  
  # Computation of Cell-Cell scores
  ccpairsgrouped_scores <- compute_CC_pairs_grouped(lrpairs=lrpairs_weights$LRpairs, cancertype="pancan")
  
  # normalized
  scale_model <- ccpairsgrouped_scores[[1]] %>% as_tibble() %>% preProcess(method = "range")
  ccpairsgrouped_scores_scale <- predict(scale_model, ccpairsgrouped_scores[[1]]) %>% as_tibble() %>% 
    rename_with(~paste0("CC_Pairs_grouped.", .))
  
  
  # ----------------------------------------- #
  # compute scores of immune response
  # ----------------------------------------- #
  
  tasks <- c("CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS") # Ock_IS was derived from publication
  tmp_file_path <- c("temp/")
  tasks_values <- compute_gold_standards(RNA.tpm=tpm_m, list_gold_standards=tasks, cancertype=cancer_type, output_file_path=tmp_file_path)
  
  # returns a list, convert to matrix
  # Into matrix
  immune_response <- do.call(cbind, lapply(tasks, function(X){
    immune_response <- t(tasks_values[[X]])
    return(immune_response)
  }))
  
  # Assess correlation between chemokines and the other correlated tasks
  tasks_cormat <- cor(immune_response)
  cor_sign <- sign(tasks_cormat[,"chemokines"])
  cor_sign <- cor_sign[names(cor_sign) != "chemokines"]
  if (all(cor_sign == -1)){
    immune_response[,"chemokines"] <- -immune_response[,"chemokines"]
  }
  
  # normalized
  scale_model <- immune_response %>% as_tibble() %>% preProcess(method = "range")
  immune_response_scale <- predict(scale_model, immune_response) %>% as_tibble() %>% 
    rename_with(~paste0("Immue_response_Y.", .))
  
  bind_cols(
    immune_response_scale,
    cell_fractions_scale,
    pathways_activity_scale,
    TF_activity_scale,
    lrpairs_weights_scale,
    ccpairsgrouped_scores_scale
  ) %>% return()
}

feature_non_scale <- function(cancer_type){
  con_RAW_COUNT <- DBI::dbConnect(drv = MariaDB(), host = "192.168.0.91", port = 3306, user = "root", password = "sempre813!", 
                                  dbname = "TCGA_Gene_expression_IlluminaHiSeq_RNASeqV2_count")
  con_TPM <- DBI::dbConnect(drv = MariaDB(), host = "192.168.0.91", port = 3306, user = "root", password = "sempre813!", 
                            dbname = "TCGA_Gene_expression_IlluminaHiSeq_RNASeqV2_TPM")
  
  if(cancer_type == "PANCAN"){
    all_list <- intersect(con_RAW_COUNT %>% dbListTables(), con_TPM %>% dbListTables())
    pan_list_count <- lapply(X = all_list, FUN = function(value){
      tbl(con_RAW_COUNT, value) %>% collect() %>% 
        dplyr::select(Gene, contains("-01A")) %>% 
        return()
    })
    pan_list_tpm <- lapply(X = all_list, FUN = function(value){
      tbl(con_TPM, value) %>% collect() %>%
        dplyr::select(Gene, contains("-01A")) %>% 
        return()
    })
    
    # tpm <- pan_list_tpm %>% purrr::reduce(., inner_join, by = "Gene")
    tpm <- pan_list_tpm %>% bind_cols()
    gene_name <- pan_list_tpm[[1]] %>% dplyr::select(Gene)
    tpm <- tpm %>% select_if(is.numeric) %>% 
      bind_cols(gene_name, .) %>% 
      distinct(Gene, .keep_all = TRUE)
    
    tpm_m <- tpm %>% dplyr::select(-Gene)
    rownames(tpm_m) <- tpm$Gene
    
    counts <- pan_list_count %>% bind_cols()
    gene_name <- pan_list_count[[1]] %>% dplyr::select(Gene)
    counts <- counts %>% select_if(is.numeric) %>% 
      bind_cols(gene_name, .) %>% 
      distinct(Gene, .keep_all = TRUE)
    
    counts_m <- counts %>% dplyr::select(-Gene) %>% as.matrix()
    rownames(counts_m) <- counts$Gene
    
  } else {
    tpm <- tbl(con_TPM, cancer_type) %>% collect() %>% 
      dplyr::select(Gene, contains("-01A")) %>% 
      distinct(Gene, .keep_all = TRUE)
    
    tpm_m <- tpm %>% dplyr::select(-Gene) %>% as.matrix()
    rownames(tpm_m) <- tpm$Gene
    
    counts <- tbl(con_RAW_COUNT, cancer_type) %>% collect() %>% 
      distinct(Gene, .keep_all = TRUE)
    counts_m <- counts %>% dplyr::select(-Gene) %>% as.matrix()
    rownames(counts_m) <- counts$Gene
  }
  
  # ----------------------------------------- #
  # compute system-based derived signatures
  # ----------------------------------------- #
  
  # Computation of cell fractions
  feature_list <- list()
  cell_fractions <- compute_cell_fractions(RNA_tpm=tpm_m)
  feature_list[[2]] <- cell_fractions %>% as_tibble() %>% 
    rename_with(~paste0("Cell_Fractions.", .)) %>% 
    bind_cols(rownames(cell_fractions) %>% tibble(Sample_barcode = .), .)
  
  # Computation of pathway activity
  pathways_activity <- compute_pathway_activity(RNA_counts=counts_m, remove_sig_genes_immune_response =TRUE)
  feature_list[[3]] <- pathways_activity %>% as_tibble() %>% 
    rename_with(~paste0("Pathways_activity.", .)) %>% 
    bind_cols(rownames(pathways_activity) %>% tibble(Sample_barcode = .), .)
  
  # Computation of TF activity
  TF_activity <- compute_TF_activity(RNA_tpm=tpm_m)
  feature_list[[4]] <- TF_activity %>% as_tibble() %>% 
    rename_with(~paste0("TF_activity.", .)) %>% 
    bind_cols(rownames(TF_activity) %>% tibble(Sample_barcode = .), .)
  
  # Computation of LR pairs weights
  lrpairs_weights <- compute_LR_pairs(RNA_tpm=tpm_m, cancer_type="pancan") 
  feature_list[[5]] <- lrpairs_weights %>% as_tibble() %>% 
    rename_with(~paste0("LR_Pairs_weights.", .)) %>% 
    bind_cols(rownames(lrpairs_weights) %>% tibble(Sample_barcode = .), .)
  
  # Computation of Cell-Cell scores
  ccpairs_scores <- compute_CC_pairs(lrpairs=lrpairs_weights, cancer_type="pancan")
  feature_list[[6]] <- ccpairs_scores %>% as_tibble() %>% 
    rename_with(~paste0("CC_Pairs_scores.", .)) %>% 
    bind_cols(rownames(ccpairs_scores) %>% tibble(Sample_barcode = .), .)
  
  # ----------------------------------------- #
  # compute scores of immune response
  # ----------------------------------------- #
  
  tasks <- c("CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS") # Ock_IS was derived from publication
  tmp_file_path <- c("temp/")
  tasks_values <- compute_scores_immune_response(RNA_tpm=tpm_m, selected_scores=tasks)
  feature_list[[1]] <- tasks_values %>% as_tibble() %>% 
    rename_with(~paste0("Gold_Standard.", .)) %>% 
    bind_cols(rownames(tasks_values) %>% tibble(Sample_barcode = .), .)
  # returns a list, convert to matrix
  # Into matrix
  # immune_response <- do.call(cbind, lapply(tasks, function(X){
  #   immune_response <- t(tasks_values[[X]])
  #   return(immune_response)
  # }))
  
  # Assess correlation between chemokines and the other correlated tasks
  # tasks_cormat <- cor(immune_response)
  # cor_sign <- sign(tasks_cormat[,"chemokines"])
  # cor_sign <- cor_sign[names(cor_sign) != "chemokines"]
  # if (all(cor_sign == -1)){
  #   immune_response[,"chemokines"] <- -immune_response[,"chemokines"]
  # }
  
  # normalized
  # feature_list[[1]] <- immune_response %>% as_tibble() %>% 
  #   rename_with(~paste0("Gold_Standard.", .)) %>% 
  #   bind_cols(rownames(immune_response) %>% tibble(Sample_barcode = .), .)
  
  feature_list %>% purrr::reduce(inner_join, by = "Sample_barcode") %>% return()
}

## RUN ----
# type <- read_delim(file = "train/traning_cancer_type.txt", delim = "\t", col_names = F) %>% pull(1)

feature <- feature_non_scale(cancer_type = "LUSC") 

feature$Sample_barcode <- str_replace(feature$Sample_barcode, "01A+-[:alnum:]+-[:alnum:]+-[:digit:]+", "01")
write_delim(feature,file = "/workspace/rawdata/autoencoder/LUSC/LUSC.txt", delim = "\t")

{
  # # non-scale
  # solid_cancer_scale <- lapply(X = type, FUN = function(t){
  #   pr <- tibble(Project = t)
  #   bind_cols(pr, feature(cancer_type = t)) %>%
  #     return()
  # })
  # # Response variable continous to category with quantile[4] -- Q3(75%)
  # load("~/Immune-response-wmbio/immune_solid_feature_non_scale.RData")
  # lapply(X = solid_cancer, FUN = function(DF){
  #   DF_temp <- DF %>% mutate(
  #     Immue_response_Y.CYT = ifelse(quantile(Immue_response_Y.CYT)[4] > Immue_response_Y.CYT, 0, 1),
  #     Immue_response_Y.Roh_IS = ifelse(quantile(Immue_response_Y.Roh_IS)[4] > Immue_response_Y.Roh_IS, 0, 1),
  #     Immue_response_Y.chemokines = ifelse(quantile(Immue_response_Y.chemokines)[4] > Immue_response_Y.chemokines, 0, 1),
  #     Immue_response_Y.Davoli_IS = ifelse(quantile(Immue_response_Y.Davoli_IS)[4] > Immue_response_Y.Davoli_IS, 0, 1),
  #     Immue_response_Y.IFNy = ifelse(quantile(Immue_response_Y.IFNy)[4] > Immue_response_Y.IFNy, 0, 1),
  #     Immue_response_Y.Ayers_expIS = ifelse(quantile(Immue_response_Y.Ayers_expIS)[4] > Immue_response_Y.Ayers_expIS, 0, 1),
  #     Immue_response_Y.Tcell_inflamed = ifelse(quantile(Immue_response_Y.Tcell_inflamed)[4] > Immue_response_Y.Tcell_inflamed, 0, 1),
  #     Immue_response_Y.RIR = ifelse(quantile(Immue_response_Y.RIR)[4] > Immue_response_Y.RIR, 0, 1),
  #     Immue_response_Y.TLS = ifelse(quantile(Immue_response_Y.TLS)[4] > Immue_response_Y.TLS, 0, 1)
  # 
  #   )
  # 
  #   pr_name <- DF_temp$Project %>% unique()
  #   if(length(pr_name) > 1)
  #     pr_name <- "SOLID_CANCER"
  # 
  #   write_delim(DF_temp, file = paste0(getwd(), "/", pr_name,".txt"), delim = "\t")
  # })
  }