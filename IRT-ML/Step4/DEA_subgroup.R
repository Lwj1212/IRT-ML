if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("survival")
install.packages("survminer")
install.packages("tidyverse")
install.packages("data.table")
install.packages("ggplot2")
install.packages("curl")
BiocManager::install("Matrix")
BiocManager::install("DESeq2")
BiocManager::install("biomaRt", force = TRUE)
BiocManager::install("TCGAbiolinks")
BiocManager::install("EDASeq", force = TRUE)
BiocManager::install('Rsamtools')
BiocManager::install('Rhtslib')
BiocManager::install('RCurl')
BiocManager::install('maftools')
BiocManager::install('UsingR')

library(UsingR)
library(maftools)
library(survival)
library(survminer)
library(data.table)
library(ggplot2)
library(DESeq2)
library(biomaRt)
library(TCGAbiolinks)
library(EDASeq)
library(Matrix)
library(BiocParallel)
library(tidyverse)
library(curl)

#########transpose 함수

transpose_df <- function(df) {
  t_df <- data.table::transpose(df)
  colnames(t_df) <- rownames(df)
  rownames(t_df) <- colnames(df)
  t_df <- t_df %>%
    tibble::rownames_to_column(.data = .) %>%
    tibble::as_tibble(.)
  
  col_name <- t_df[1, ] %>% as.character()
  t_df <- t_df[-1, ]
  colnames(t_df) <- col_name
  
  # chr to numeric
  t_df %>% mutate_at(vars(-names(.)[1]), as.numeric) %>%
    return()
}

setwd("/workspace/Machine_Learning/Pipeline_1/LUSC")       # Base 폴더 지정

########데이터 불러오기

deseq_count <- read.table("/workspace/Machine_Learning/Pipeline_1/LUSC/DEA_count/DEA_count.txt", row.names = 1, header = T, sep = '\t')
count_keep <- rowSums(deseq_count) >= 10
colnames(deseq_count) <- gsub("\\.","\\-",colnames(deseq_count))
deseq_countfilt <- deseq_count[count_keep,]

###########Official gene symbol 치환
httr::set_config(httr::config(ssl_verifypeer = FALSE))
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

gene_name <- deseq_countfilt %>% 
  rownames() %>% lapply(X = ., FUN = function(value){
    strsplit(value, split = "\\.") %>% unlist() %>% .[1] %>% return()
  }) %>% unlist() %>% 
  getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
        filters = "ensembl_gene_id", values = .,
        mart = mart) %>% as_tibble()

gene_name_filt <- gene_name %>% filter(hgnc_symbol != '' & hgnc_symbol != 'LINC01238' & hgnc_symbol != 'POLR2J4')  # duplicate gene 제거
ensembl_name <- deseq_countfilt %>% 
  rownames_to_column('rn') %>% pull(1) %>% 
  lapply(X = ., FUN = function(value){
    strsplit(value, split = "\\.") %>% unlist() %>% .[1] %>% as_tibble() %>% return()
  }) %>% bind_rows()

deseq_countfilt_ensembl <- bind_cols(ensembl_name, deseq_countfilt) %>% 
  dplyr::rename(., "ensembl_gene_id" = "value")

deseq_countfilt_og <- inner_join(gene_name_filt, deseq_countfilt_ensembl, by = 'ensembl_gene_id')


########### zero gene remove

deseq_gene_final <- deseq_countfilt_og[apply(deseq_countfilt_og != 0, 1, sum) >= as.integer(ncol(deseq_countfilt_og) * 0.8), ] %>%
  subset(,select=-c(ensembl_gene_id)) %>% column_to_rownames(.,'hgnc_symbol')


# row_name <- deseq_gene_final[,1] %>% as.character()
# deseq_gene_final <- deseq_gene_final[,-1]
# rownames(deseq_gene_final) <- row_name
src_dir <- c("/workspace/Machine_Learning/Pipeline_1/LUSC/DEA_model") 
src_file <- list.files(src_dir)
src_file_cnt <- length(src_file)
df_list <- list()


for(i in 1:src_file_cnt) { # write.table one by one automatiically, using loop program
  deseq_group <- read.table(
    paste(src_dir, "/", src_file[i], sep=""),
    sep="\t",
    header=T,
    stringsAsFactors = F)
  
  tmp_tsne <- names(deseq_group)[[2]]
  
  colnames(deseq_group) <- c("Sample_barcode", "tsne")
  deseq_group$tsne <- as.factor(deseq_group$tsne)

  ########### data 잘라서 돌려보기 위한 작업
  
  # temp1 <- deseq_count[1:1000,]
  # temp3 <- temp1 %>% mutate_all(function(x){x+1})
  # colnames(temp1) <- gsub("\\.","\\-",colnames(temp1))
  # temp2 <- data.frame(deseq_group)%>% dplyr::slice(1:1000)
  # temp2$tsne1 <- as.factor(temp2$tsne1)
  
  
  ############# DESeq2
  
  dds <- DESeqDataSetFromMatrix(deseq_gene_final, deseq_group, ~ tsne)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds <- DESeq(dds, test="Wald", parallel = TRUE, BPPARAM=MulticoreParam(4))
  res <- results(dds, alpha = 0.05)
  resOrdered <- res[order(res$pvalue),]
  sum(res$padj < 0.05, na.rm=TRUE)
  resSig <- subset(resOrdered, padj < 0.05)
  resSig <- data.frame(resSig)
  resfinal <- subset(resSig, abs(log2FoldChange) > 1)   # !!! group=1일 때 OS 감소 => log2FC > 1 gene / OS 증가 => log2FC < -1 gene
  sum(resfinal$padj < 0.05, na.rm=TRUE)
  summary(resfinal)
  
  dir.create(paste0("/workspace/Machine_Learning/Pipeline_1/LUSC/DEA/", tmp_tsne))
  
  png(paste0("/workspace/Machine_Learning/Pipeline_1/LUSC/DEA/", tmp_tsne,"/DEA_MAplot_",tmp_tsne,".png"), width = 700)
  plotMA(res, ylim=c(-15,15), main = paste0("DEA MA-plot in ",tmp_tsne)) ### 꼭 저장할 것!!
  dev.off()
  
  write.table(resfinal,file = paste0("/workspace/Machine_Learning/Pipeline_1/LUSC/DEA/", tmp_tsne,"/DEA_result_final_FDR0.05_logFC+-1_", tmp_tsne,".txt"),sep = '\t', row.names = TRUE,col.names=NA,quote = F)
  write.table(res,file = paste0("/workspace/Machine_Learning/Pipeline_1/LUSC/DEA/", tmp_tsne,"/DEA_result_final_", tmp_tsne, ".txt"),sep = '\t', row.names = TRUE,col.names=NA,quote = F)
  
  
}





