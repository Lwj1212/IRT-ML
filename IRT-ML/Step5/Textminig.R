library(RMariaDB)
library(tidyverse)
library(data.table)


setwd("/workspace/report/Target ID/IRT_ML/LUAD")
gene <- read.table("/workspace/report/Target ID/IRT_ML/LUAD/Gene_cluster_irtml_LUAD.txt", header = T)

con <- DBI::dbConnect(drv = MariaDB(), host = "192.168.0.91", port = 3306, 
                      user = "root", 
                      password = "sempre813!",
                      dbname = "Textmining")
gene_list <- gene[[1]]
# gene_list <- "ENSG00000147255"
dbListTables(con)
# cancer_type <- cancer_list
cancer_type <- c("LUAD")

df <- tbl(con, cancer_type) %>% collect() %>% 
    dplyr::filter(gene %in% gene_list) %>% 
  rename("Gene" = "gene")

gene_table <- gene_list %>% as.data.frame() %>% rename ("." = "Gene") %>%
  full_join(., df, by = "Gene")

write_delim(df, file = "Textmining_irtml_target_LUAD_eachcancer.txt", delim = "\t")
