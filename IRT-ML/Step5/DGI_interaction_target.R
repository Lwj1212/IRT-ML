# BiocManager::install("httr", force = T)


library(httr)
library(tidyverse)

# DGIdb
dgidb_interaction <- function(gene_name){
  base_url <- "https://dgidb.org"
  request_url <- paste0(base_url, "/api/v2/interactions.json?")
  result_list <- list()
  
  # chunk id
  id_chunk <- split(gene_name, ceiling(seq_along(gene_name)/200))
  
  for(index in 1:length(id_chunk)){
    # print(index)
    payload <-  list(genes = paste0(id_chunk[[index]], collapse = ","),
                     fda_approved_drug="true")
    
    # output
    dgidb_result <- POST(request_url, body = payload, encode = "form") %>%  
      httr::content(encoding = "UTF-8") 
    
    result_list[[index]] <- lapply(X = dgidb_result$matchedTerms, FUN = function(dgidb_element){
      gene_category <- dgidb_element$geneCategories %>% 
        sapply(X = ., FUN = function(value) {value$name}) %>% 
        paste0(collapse = ",")
      
      interaction <- dgidb_element$interactions %>% 
        sapply(X = ., FUN = function(value){
          drug_name <- value$drugName
          score <- value$score
          types <- value$interactionTypes %>% unlist() %>% paste0(collapse = "&")
          
          paste0(c(drug_name, score, types), collapse = ";") %>% 
            as_tibble() %>% 
            return()
          
          # return(drug_name)  
        }) %>% unlist() %>% 
        paste0(., collapse = "&")
      
      tibble(
        gene = dgidb_element$geneName,
        DGI_GENE_CATEGORY = gene_category, 
        `DGI(DRUG_NAME;SCORE;TYPE)` = interaction,
        DGI_COUNT = length(dgidb_element$interactions)
      )  %>% return()
      
    }) %>% bind_rows() 
  }
  
  result_list %>% bind_rows() %>% return()
}

gene <- read.table("/workspace/report/Target ID/IRT_ML/PAAD/Gene_cluster_irtml_PAAD.txt", header = T)
gene_list <- gene[[1]]

a <- dgidb_interaction(gene_list)
