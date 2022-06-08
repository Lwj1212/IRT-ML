# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# install.packages("devtools")
# install.packages("survival")
# install.packages("survminer")
# install.packages("tidyverse")
# install.packages("data.table")
# install.packages("ggplot2")
# 
# devtools::install_github("r-lib/conflicted")
# BiocManager::install("Matrix")
# BiocManager::install("DESeq2")
# BiocManager::install("biomaRt", force = TRUE)
# BiocManager::install(version = "3.14","XVector")
# BiocManager::install("TCGAbiolinks")
# BiocManager::install("EDASeq", force = TRUE)
# BiocManager::install('Rsamtools', force = TRUE)
# BiocManager::install('Rhtslib', force = TRUE)
# BiocManager::install('RCurl', force = TRUE)
# BiocManager::install('maftools')
# BiocManager::install('UsingR')
# BiocManager::install('survcomp')


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
library(survcomp)
library(tidyverse)

########### transpose 함수

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

########### Grouping data 가져오기

setwd("/workspace/Machine_Learning/Pipeline_1/LUSC")       # Base 폴더 지정

########## 그룹데이터 모아놓은 폴더 입력 (폴더에는 pipeline에서 뽑은 그룹데이터만 포함할 것!!!)

src_dir <- c("/workspace/Machine_Learning/Pipeline_1/LUSC/AutoEncoder/group")    
src_file <- list.files(src_dir)
src_file_cnt <- length(src_file)

#########그룹 데이터 하나의 파일로 통합

df_list <- list()
for(i in 1:src_file_cnt) { # write.table one by one automatiically, using loop program
  tmp <- read.table(
    paste(src_dir, "/", src_file[i], sep=""),
    sep="\t",
    header=T,
    stringsAsFactors = F) %>% 
    dplyr::select(-X)
  
  colnames(tmp) <- c("Sample_barcode", paste0("km_tsne_", i))
  df_list[[i]] <- tmp
}

temp_res <- purrr::reduce(df_list, inner_join, by = "Sample_barcode")
temp_res <- temp_res %>% filter(Sample_barcode != "TCGA-21-1076-01" & Sample_barcode != "TCGA-6-156-01")  # overlap sample 제거


sil <- read.table("/workspace/Machine_Learning/Pipeline_1/LUSC/AutoEncoder/group_silhouette.csv",header = T, sep = ',')
anv <- read.table("/workspace/Machine_Learning/Pipeline_1/LUSC/AutoEncoder/group_anova.csv",header = T, sep = ',')
rf <- read.table("/workspace/Machine_Learning/Pipeline_1/LUSC/AutoEncoder/group_rf.csv",header = T, sep = ',')

metric <- sil %>% inner_join(anv, by = "Sample_barcode") %>% dplyr::inner_join(rf, by = "Sample_barcode")

temp_val <- metric %>% cbind(.,transpose_df(temp_res)) %>% dplyr::select(-5)

temp_res_raw <- temp_val %>% dplyr::select(-c(2,3,4)) %>% transpose_df()

temp_res_filt <- temp_val %>% dplyr::filter(., Silhouette > 0.5 & Silhouette != 1 & svm_anova > 0.9 & svm_anova != 1 &
                                              svm_rf > 0.9 & svm_rf != 1) %>%
  dplyr::select(-c(2,3,4)) %>% transpose_df()



write.table(temp_res_raw, "group_alldata.txt", row.names = F, col.names = T, quote = F, sep = '\t') # 그룹 통합데이터 저장
write.table(temp_res_filt, "group_filtered.txt", row.names = F, col.names = T, quote = F, sep = '\t') # 그룹 filtering 데이터 저장




######### Raw count data 가져오기 (for DESeq2)

temp_count <- read_delim("/workspace/Machine_Learning/rawdata/tcga_gene_expected_count.txt", delim = '\t')
temp_count <- temp_count %>% mutate_if(is.numeric, function(x) {
  ifelse(x < -9.9, 0, 2^x-1) %>% round(.,digits = 0) %>% 
    return()}
)
temp_count <- transpose_df(temp_count)
colnames(temp_count)[1] <- "Sample_barcode"

############## Survival data 가져오기

temp_surv <- read.delim("/workspace/Machine_Learning/rawdata/Survival_SupplementalTable_S1_20171025_xena_sp.txt")
temp_surv <- dplyr::select(temp_surv, c(sample,OS,OS.time,DSS,DSS.time,DFI,DFI.time,PFI,PFI.time))
names(temp_surv)[names(temp_surv) == "sample"] <- c("Sample_barcode")


########## inner_join (model selection을 위해)

temp_OS <- inner_join(temp_res_filt,temp_surv)
temp <- inner_join(temp_OS,temp_count)
temp_sample <- temp %>% dplyr::select(c("Sample_barcode"))

############## Model selection
coldata <- colnames(temp_res_filt[,2:906])
coldata <- data.frame(coldata)

col_names <- temp %>% colnames()
col_names <- col_names[str_detect(col_names, "km")]

df_log_rank <- lapply(X = col_names, function(col_name){
  surv_diff <- survdiff(
    formula = as.formula(paste0("Surv(time = OS.time, event = OS) ~", col_name)), 
    data = temp)
  # surv_diff_result <- summary(surv_diff)
  chisq_p_value = surv_diff$chisq
  tmp_obs = surv_diff$obs
  tmp <- data.frame(chisq_p_value)
}) %>% bind_rows()

df_log_rank <- bind_cols(coldata,df_log_rank)

df_log_rank <- df_log_rank %>%  arrange(desc(chisq_p_value))

############ Model selection data 저장 (Survdiff 에 대한 chisq 표기)

write.table(df_log_rank, "group_survrank.txt", row.names = F, col.names = T, quote = F, sep = '\t')


df_log_rank_top <- df_log_rank[1:905,]
df_sil_top <- df_log_rank %>% inner_join(sil, by = c("coldata" = "Sample_barcode")) %>% arrange(desc(Silhouette)) %>% 
  dplyr::select(c(1,3))

write.table(df_sil_top, "group_silhouette.txt", row.names = F, col.names = T, quote = F, sep = '\t')

df_sil_top500 <- df_sil_top[1:3238,]

############## Selected model 변수로 재저장

temp_resa <- dplyr::select(temp_res_filt, c("Sample_barcode", df_log_rank_top$coldata))
# temp_res <- rename(temp_res, "km_tsne" = "km_tsne_9")   # km_tsne_10 자리에는 chisq 가장 큰 model 선택

############ Survival validation

temp_col <- colnames(temp_resa) %>% as.data.frame()
colnames(temp_col) <- c("tsne")
temp_col <- temp_col %>% filter(tsne!=c("Sample_barcode"))




df_km_list <- list()
df_km_p_list <- list()
df_km_concordance <- list()

for (i in temp_col$tsne){
 
  
  temp_dummy <- temp_resa %>% dplyr::select(c("Sample_barcode", i)) %>% as.data.frame()

  
  temp_OSa <- inner_join(temp_dummy,temp_surv) %>% dplyr::select(c("Sample_barcode",i,"OS","OS.time")) %>% as.data.frame()
  # temp_OSa[,2] <- as.factor(temp_OSa[,2])
  temp_tsne <- colnames(temp_OSa[,2])
  tempa <- inner_join(temp_OSa, temp_count)
  temp_km <- inner_join(temp_OSa,temp_sample)

  temp_km_name <- names(temp_km)[[2]]
# km_object <- Surv(event = temp_km$OS, time = as.numeric(temp_km$OS.time))

km_os <- survfit(as.formula(paste0("Surv(event = OS, time = OS.time) ~", i)), data=temp_km)

km_os_summary <- summary(km_os)$table %>% as.data.frame()

group_reverse <- TRUE

if(km_os_summary[1,5] > km_os_summary[2,5]){
  group_reverse <- FALSE
}

if(group_reverse){
  temp_km <- temp_km %>% mutate_at(c(i), ~ (ifelse(. == 0 , 1, 0)))
  }



km_os <- survfit(as.formula(paste0("Surv(event = OS, time = OS.time) ~", i)), data=temp_km)

km_os_summary <- summary(km_os)$table %>% as.data.frame()
km_p_value <- surv_pvalue(km_os)[,2] %>% as.data.frame()
rownames(km_p_value) <- i
colnames(km_p_value) <- c("pvalue")
fit <- coxph(as.formula(paste0("Surv(event = OS, time = OS.time) ~", i)), data=temp_km)
a <- concordance(fit)
km_concordance <- a$concordance %>% as.data.frame()
rownames(km_concordance) <- i
colnames(km_concordance) <- c("CI")

df_km_list[[i]] <- km_os_summary
df_km_p_list[[i]] <- km_p_value
df_km_concordance[[i]] <- km_concordance



survp <- ggsurvplot(km_os,
                    conf.int = TRUE, # 신뢰구간 표현 여부
                    risk.table = FALSE, # 테이블 표시 여부
                    risk.table.height = 0.4, # 테이블 높이 설정
                    ggtheme = theme_bw(), # 데이터 테마 설정
                    # palette = c("#2E9FDF"), 
                    font.x = c(10), # x축 제목 크기 설정
                    font.y = c(10), # y축 제목 크기 설정
                    font.tickslab = c(10), # 축 값 크기 설정
                    surv.median.line = "hv", # 50% 생존지점 표시
                    break.time.by = 1000, xlim = c(0,5000),
                    pval = T,
                    pval.coord = c(0.2,0.2)
)

ggexport(filename = paste0("/workspace/Machine_Learning/Pipeline_1/LUSC/OS_validation/survplot_",i,".png"), plot = survp$plot, width = 700)




############## inner join (DEA를 위한 데이터 변환을 위해)

temp_survfilt <- tempa[,c(1,3:4)]
temp_countfilt <- tempa[,-c(2:4)]
temp_resfilt <- tempa[,c(1:2)]
temp_countfilt <- transpose_df(temp_countfilt)
temp_resfilt <- as.tibble(temp_resfilt)

######### DESeq2 format 맞춰주기

namea = temp_countfilt[,1]
deseq_count <- data.frame(temp_countfilt, row.names = namea$Sample_barcode)
deseq_count <- deseq_count %>% dplyr::select(-Sample_barcode)
colnames(deseq_count) <- gsub("\\.","\\-",colnames(deseq_count))
nameb = temp_resfilt[,1]
deseq_group <- data.frame(temp_resfilt, row.names = nameb$Sample_barcode)
deseq_group <- deseq_group %>% dplyr::select(i)
all(rownames(deseq_group) == colnames(deseq_count))
write.table(deseq_group,file = paste0("/workspace/Machine_Learning/Pipeline_1/LUSC/DEA_model/DEA_group_",i,".txt"),sep = '\t', row.names = TRUE,col.names=NA)
}

df_total_ci <- df_km_concordance %>% bind_rows() %>% 
  write.table(., file = "/workspace/Machine_Learning/Pipeline_1/LUSC/Total_model_OS_validation_CI.txt", sep = '\t',row.names = T, col.names = T)

df_total_p <- df_km_p_list %>% bind_rows() %>%
  write.table(., file = "/workspace/Machine_Learning/Pipeline_1/LUSC/Total_OS_validation_pval.txt", sep = '\t',row.names = T, col.names = T)


df_filtered_ci <- df_km_concordance %>% bind_rows() %>%  filter(.,CI >= 0.7) %>% 
  write.table(., file = "/workspace/Machine_Learning/Pipeline_1/LUSC/Significant_model_OS_validation_CI.txt", sep = '\t',row.names = T, col.names = T)

df_filtered_p <- df_km_p_list %>% bind_rows() %>% filter(.,pvalue < 0.05) %>% 
  write.table(., file = "/workspace/Machine_Learning/Pipeline_1/LUSC/Significant_model_OS_validation_pval.txt", sep = '\t',row.names = T, col.names = T)

df_total_km <- df_km_list %>% bind_rows() %>% 
  write.table(., file = "/workspace/Machine_Learning/Pipeline_1/LUSC/Significant_model_OS_validation_summary.txt", sep = '\t',row.names = T, col.names = T)

df_total_sil <- df_sil_top %>% bind_rows() %>% 
  write.table(., file = "/workspace/Machine_Learning/Pipeline_1/LUSC/Total_model_silhouette_score.txt", sep = '\t',row.names = T, col.names = T)

write.table(deseq_count,"/workspace/Machine_Learning/Pipeline_1/LUSC/DEA_count/DEA_count.txt",sep = '\t', row.names = TRUE,col.names=NA)
