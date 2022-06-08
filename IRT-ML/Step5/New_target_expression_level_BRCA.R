library(RMariaDB)
library(tidyverse)

# BRCA의 경우
load("/workspace/dummy/BRCA/TCGA_Gene_expression_TPM_BRCA.RData")
df_brca <- result_tpm

gene <- read.table("/workspace/report/Target ID/IRT_ML/STAD/Gene_cluster_irtml_STAD_final.txt", header = T)
gene_list <- gene[[2]]

df_list <- list()
for (i in gene_list){
  df <- df_brca %>%
    dplyr::select(sample, contains("-01")) %>% 
    dplyr::filter(sample == i) %>%
    t() %>% as.data.frame() %>% 
    rownames_to_column(., "Sample_barcode")
  
  df <- df[-1,]
  colnames(df) <- c("Sample_barcode", i)
  
  df_list[[i]] <- df
}

total_df <- purrr::reduce(df_list, full_join, by = "Sample_barcode")

exp2 <- function(x){
  (2^x)-0.001
}
total_df[,gene_list]  <- lapply(total_df[,gene_list], as.numeric) 
total_df[,gene_list] <- lapply(total_df[,gene_list], exp2)
fun <- function(s){
  ifelse(s < 0, 0, s)
}
total_df[,gene_list] <- lapply(total_df[,gene_list], fun)

total_df1 <- total_df %>% column_to_rownames(., var = "Sample_barcode") %>% 
  t() %>% as.data.frame() %>% rownames_to_column(., var = "ensembl") %>% 
  full_join(gene, ., by = "ensembl")

df_final <- total_df1[,-2] %>% column_to_rownames(., var = "Gene") %>% t() %>% 
  as.data.frame() %>% rownames_to_column(., "Sample_barcode")
df_final$Project <- c("BRCA")
write_delim(df_final,file = "TCGA_gene_expression_irtml_target_STAD_final_BRCA.txt", delim = "\t")
