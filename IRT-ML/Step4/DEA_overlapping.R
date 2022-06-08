setwd("/workspace/Machine_Learning/Pipeline_1/LUSC/overlapping")


library(tidyverse)

src_dir <- c("/workspace/Machine_Learning/Pipeline_1/LUSC/DEA")    
src_file <- list.files(src_dir, pattern = "DEA_result_final_FDR0.05_logFC\\+\\-1_*", full.names = T, recursive = T)
src_file_cnt <- length(src_file)

#########그룹 데이터 하나의 파일로 통합

df_list <- list()
for(i in 1:src_file_cnt) { # write.table one by one automatically, using loop program
  tmp <- read.table(
    paste(src_file[i], sep=""),
    sep="\t",
    header=T,
    stringsAsFactors = F)
    
  
  # colnames(tmp) <- c("Sample_barcode", paste0("km_tsne_", i))
  df_list[[i]] <- tmp
}


save(df_list, file = "df_list1.RData")
load("df_list1.RData")
temp_dea <- purrr::reduce(df_list, full_join, by = "X")

temp_dea_final_100 <- temp_dea[apply(is.na(temp_dea), 1, sum) <= as.integer(ncol(temp_dea) * 0),]

kingene <- read.table("/workspace/Machine_Learning/Pipeline_1/Pancan/DEA_overlapping_result/kinase_gene.txt", header=T)
colnames(kingene) <- c("X")

temp_dea_100_kin <- inner_join(temp_dea_final_100, kingene, by = "X")
final_100 <- dplyr::select(temp_dea_final_100, "X")

dea_2463 <- read.table("/workspace/Machine_Learning/Pipeline_1/LUSC/DEA/km_tsne_2463/DEA_result_final_FDR0.05_logFC+-1_km_tsne_2463.txt", header = T) %>%  
  rownames_to_column(., var = "X") %>%  dplyr::select(X,log2FoldChange, padj)
colnames(dea_2463) <- c("X", "log2FC_tsne2463", "padj_tsne2463")

dea_4794 <- read.table("/workspace/Machine_Learning/Pipeline_1/LUSC/DEA/km_tsne_4794/DEA_result_final_FDR0.05_logFC+-1_km_tsne_4794.txt", header = T) %>% 
  rownames_to_column(., var = "X") %>%  dplyr::select(X,log2FoldChange, padj)
colnames(dea_4794) <- c("X", "log2FC_tsne4794", "padj_tsne4794")

dea_1104 <- read.table("/workspace/Machine_Learning/Pipeline_1/LUSC/DEA/km_tsne_1104/DEA_result_final_FDR0.05_logFC+-1_km_tsne_1104.txt", header = T) %>% 
  rownames_to_column(., var = "X") %>%  dplyr::select(X,log2FoldChange, padj)
colnames(dea_1104) <- c("X", "log2FC_tsne1104", "padj_tsne1104")

dea_3963 <- read.table("/workspace/Machine_Learning/Pipeline_1/LUSC/DEA/km_tsne_3963/DEA_result_final_FDR0.05_logFC+-1_km_tsne_3963.txt", header = T) %>% 
  rownames_to_column(., var = "X") %>%  dplyr::select(X,log2FoldChange, padj)
colnames(dea_3963) <- c("X", "log2FC_tsne3963", "padj_tsne3963")

dea_4938 <- read.table("/workspace/Machine_Learning/Pipeline_1/LUSC/DEA/km_tsne_4938/DEA_result_final_FDR0.05_logFC+-1_km_tsne_4938.txt", header = T) %>% 
  rownames_to_column(., var = "X") %>%  dplyr::select(X,log2FoldChange, padj)
colnames(dea_4938) <- c("X", "log2FC_tsne4938", "padj_tsne4938")

dea_4693 <- read.table("/workspace/Machine_Learning/Pipeline_1/LUSC/DEA/km_tsne_4693/DEA_result_final_FDR0.05_logFC+-1_km_tsne_4693.txt", header = T) %>% 
  rownames_to_column(., var = "X") %>%  dplyr::select(X,log2FoldChange, padj)
colnames(dea_4693) <- c("X", "log2FC_tsne4693", "padj_tsne4693")

dea_4852 <- read.table("/workspace/Machine_Learning/Pipeline_1/LUSC/DEA/km_tsne_4852/DEA_result_final_FDR0.05_logFC+-1_km_tsne_4852.txt", header = T) %>% 
  rownames_to_column(., var = "X") %>%  dplyr::select(X,log2FoldChange, padj)
colnames(dea_4852) <- c("X", "log2FC_tsne4852", "padj_tsne4852")

dea_1219 <- read.table("/workspace/Machine_Learning/Pipeline_1/LUSC/DEA/km_tsne_1219/DEA_result_final_FDR0.05_logFC+-1_km_tsne_1219.txt", header = T) %>% 
  rownames_to_column(., var = "X") %>%  dplyr::select(X,log2FoldChange, padj)
colnames(dea_1219) <- c("X", "log2FC_tsne1219", "padj_tsne1219")

dea_2980 <- read.table("/workspace/Machine_Learning/Pipeline_1/LUSC/DEA/km_tsne_2980/DEA_result_final_FDR0.05_logFC+-1_km_tsne_2980.txt", header = T) %>% 
  rownames_to_column(., var = "X") %>%  dplyr::select(X,log2FoldChange, padj)
colnames(dea_2980) <- c("X", "log2FC_tsne2980", "padj_tsne2980")

dea_4980 <- read.table("/workspace/Machine_Learning/Pipeline_1/LUSC/DEA/km_tsne_4980/DEA_result_final_FDR0.05_logFC+-1_km_tsne_4980.txt", header = T) %>% 
  rownames_to_column(., var = "X") %>%  dplyr::select(X,log2FoldChange, padj)
colnames(dea_4980) <- c("X", "log2FC_tsne4980", "padj_tsne4980")

dea_list <- list(dea_2463,dea_4794,dea_1104,dea_3963,dea_4938,dea_4693, dea_4852, dea_1219, dea_2980, dea_4980)
dea_top10 <- purrr::reduce(dea_list, full_join, by = "X")

dea_final <- inner_join(final_100, dea_top10, by = "X")


write.table(dea_final,"/workspace/Machine_Learning/Pipeline_1/LUSC/overlapping/Gene_final_overlapping_100_top10.txt",sep ='\t', row.names = TRUE,col.names=NA,quote = F)



