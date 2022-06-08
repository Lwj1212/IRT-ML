library(RMariaDB)
library(tidyverse)
library(data.table)

a <- fread("/workspace/Machine_Learning/rawdata/DEA_eachcancer/Solid_tumor_pr.txt")
a
cancer_list <- a$TCGA
cancer_list

gene <- read.table("/workspace/report/Target ID/IRT_ML/STAD/Gene_cluster_irtml_STAD.txt", header = T)

con <- DBI::dbConnect(drv = MariaDB(), host = "192.168.0.91", port = 3306, 
                      user = "root", 
                      password = "sempre813!",
                      dbname = "TCGA_Gene_expression_TOIL_RSEM_TPM")
gene_list <- gene[[2]]
# gene_list <- "ENSG00000147255"
dbListTables(con)
cancer_type <- cancer_list
cancer_type <- c("STAD")

df_list <- list()
for (i in gene_list){
  df <- tbl(con, cancer_type) %>% collect() %>% 
    dplyr::select(sample, contains(c("-01"))) %>% 
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
df_final$Project <- cancer_type
write_delim(df_final, file = "TCGA_gene_expression_irtml_target_STAD.txt", delim = "\t")


