a <- fread("/workspace/Machine_Learning/rawdata/DEA_eachcancer/Solid_tumor_pr.txt")
a
df_list <- a$TCGA
df_list



run_deseq_normal <- function(pr_name, rdata_path){
  register(MulticoreParam(20))
  suppressMessages({
    if((!file.exists(paste0(rdata_path,"/", pr_name, "_normal.RData"))) | 
       (!file.exists(paste0(rdata_path,"/", pr_name, "_RnaseqSE_normal.RData")))){
      query <- GDCquery(project = paste0("TCGA-", pr_name), 
                        data.category = "Gene expression",
                        data.type = "Gene expression quantification",
                        experimental.strategy = "RNA-Seq",
                        platform = "Illumina HiSeq",
                        file.type = "results",
                        sample.type = c("Primary Tumor", "Solid Tissue Normal"), 
                        legacy = TRUE)
      
      GDCdownload(query)
      RnaseqSE <- GDCprepare(query)
      
      save(RnaseqSE, file = paste0(rdata_path,"/", pr_name, "_RnaseqSE_normal.RData"))
      
      Rnaseq_CorOutliers <- assay(RnaseqSE) # to matrix
      
      # normalization of genes, # quantile filter of genes
      dataNorm <- TCGAanalyze_Normalization(tabDF = Rnaseq_CorOutliers, geneInfo =  geneInfo)
      dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                        method = "quantile", 
                                        qnt.cut =  0.25)
      
      save(dataFilt, file = paste0(rdata_path, pr_name, "_normal.RData"))
    } else {
      load(paste0(rdata_path,"/", pr_name, "_RnaseqSE_normal.RData"))
      load(paste0(rdata_path,"/", pr_name, "_normal.RData"))
    }
    
    
    # selection of normal samples "NT"
    samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), typesample = c("NT")) %>% 
      as_tibble() %>% 
      mutate(group = 0) %>% 
      dplyr::rename(sample = value)
    
    # selection of tumor samples "TP"
    samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), typesample = c("TP")) %>% 
      as_tibble() %>% 
      mutate(group = 1) %>% 
      dplyr::rename(sample = value)
    
    metadata <- bind_rows(samplesNT, samplesTP) %>% 
      mutate(group = ifelse(group == 0, "NT", "TP"))
    metadata$group <- factor(metadata$group, levels = c("NT", "TP"))
    
    tcga_se <- DESeqDataSetFromMatrix(countData = dataFilt, colData = metadata, design = ~ group)
    tcga_deseq <- DESeq(tcga_se, parallel = TRUE)
    
    tcga_deseq_result <- results(tcga_deseq, contrast=c("group", "TP", "NT"))
    tcga_deseq_result_tidy <- results(tcga_deseq, tidy = TRUE, contrast=c("group", "TP", "NT"))
    
  #   # volcano plot
  #   p <- EnhancedVolcano(tcga_deseq_result,
  #                        lab = rownames(tcga_deseq_result),
  #                        x = 'log2FoldChange',
  #                        y = 'padj',
  #                        title = 'Primary Tumor versus Solid Tissue Normal',
  #                        pCutoff = 0.05,
  #                        FCcutoff = 1.5,
  #                        pointSize = 3.0,
  #                        labSize = 6.0)
  #   
  #   ggsave(plot = p, filename = paste0(deg_path, pr_name, "_DESEQ2_normal_volcano.png"), height = 8, width = 12, dpi = 70)    
  })
  
  return(tcga_deseq_result_tidy)
}

# df <- lapply(X = df_list, FUN = function(pr_name){
#   temp <- run_deseq_normal(pr_name,"/workspace/rawdata/DEA_eachcancer/gdc")
# })

  
LUSC <-  run_deseq_normal("LUSC","/workspace/Machine_Learning/rawdata/DEA_eachcancer/") 

write.table(LUSC, "/workspace/Machine_Learning/rawdata/DEA_eachcancer/LUSC_DEA.txt", quote = F, row.names = F, col.names = T, sep = '\t')
#파일에서 "row" 열 이름을 "X"로 직접 변경
LUSC <- read.table("/workspace/Machine_Learning/rawdata/DEA_eachcancer/LUSC_DEA.txt", header = T)
tmp <- read.table("/workspace/Machine_Learning/Pipeline_1/LUSC/overlapping/Gene_final_overlapping_100_top10.txt", header = T)
LUSC <- LUSC %>% dplyr::select(., c(1,3,7))
final <- dplyr::left_join(tmp, LUSC, by = "X")

write.table(final, "/workspace/Machine_Learning/Pipeline_1/LUSC/overlapping/LUSC_final_overlapping_DEA_LUSC.txt", quote = F, row.names = F, col.names = T, sep = '\t')

final %>% filter(., padj < 0.1) %>% filter(., log2FoldChange > 0) %>% na.omit()
