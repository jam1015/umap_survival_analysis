#creating files to normalize genomic data
rm(list = ls())
setwd("C:\\Users/JAM526/pancan_fpkmuq/met_analysis/scripts/")
library("GenomicFeatures")
library("biomaRt")
library("tidyverse")
library("tseries")
library("magrittr")
exp_dec <- function(x) {relu((2^x) - 0.001)}


transpose_df <- function(df) {
gene <- df$gene
numeric <-  within(df, rm("gene") )
samples <- colnames(numeric)
rm(list = deparse(substitute(df)),pos = 1)
gc()
rm(list = quo_name(quo(df)))
gc()
nr <- nrow(numeric)
nc <- ncol(numeric)
numeric <- c(numeric)

mat <- NULL
for(j in 1:nc){
 mat <- rbind(mat,numeric[1:nr])
 numeric <- tail(numeric,-nr)  

 if(mod(j,100) == 0){print(j)}

 }







# 
 colnames(numeric) <- gene

 mat <- as.tibble(mat)
 mat <- add_column(mat,sample = samples,.before = 1)
   return(mat)
}



data_dir <- "C:/Users/JAM526/pancan_fpkmuq/met_analysis/data"

#doing lee
lee_fpkm <- readRDS(file.path(data_dir,"Lee","lee_fpkm.rds"))#"C:/Users/JAM526/pancan_fpkmuq/met_analysis/data/Lee/lee_fpkm.rds")
lee_gene <- rownames(lee_fpkm)
lee_fpkm <- as.tibble(lee_fpkm)
lee_fpkm <- add_column(lee_fpkm,gene = lee_gene,.before = 1)


lee_meta <- read_csv(file.path(data_dir,"Lee","lee_case_ids.csv"))
lee_meta <- dplyr::rename(lee_meta,"sample" = "samples",  "sample_type" = "sample_type.samples" ,  "biopsy_tissue" = "Site" ,  "dataset" = "Provider"     )
lee_meta <- (lee_meta %>% mutate(biopsy_tissue = replace(biopsy_tissue, is.na(biopsy_tissue),"Breast")) )
lee_meta <- mutate(lee_meta,tissue = "Breast")
lee_meta <- mutate(lee_meta,disease = "BRCA")


#doing met500

       met500_fpkm <- read_delim(file.path(data_dir,"met500","met500_fpkm.txt"),delim = "\t")
       split_genes <- strsplit(x = met500_fpkm$gene,split = ".",fixed = TRUE)
  ens_gene_id_base <- unlist(lapply(split_genes,"[[",1))
ens_gene_id_suffix <- unlist(lapply(split_genes,"[[",2))
  met500_fpkm$gene <- ens_gene_id_base


met500_meta <- read_delim(file.path(data_dir,"met500","met500_meta.tsv"),delim = "\t")
met500_meta <- met500_meta[,2:ncol(met500_meta)] #deleting extraneous "samples" variable
met500_meta <- dplyr::rename(met500_meta,  "sample"="samples",  "disease"="cohort","patient_id"="sample_source","run_id"="run.id")
met500_meta <-  (met500_meta %>% mutate(sample_type = ifelse(biopsy_tissue == tissue,"Primary Tumor","Metastatic"))) #adding a sample type variable
met500_meta <-   mutate(met500_meta,dataset = "met500") #what dateset
met500_meta <- within(met500_meta, rm("run_id") )
met500_meta[is.na(met500_meta$sample_type),"sample_type"] = "Metastatic"

#tcga


#this was the preprocessing that took a long time
#tcga_fpkm <- read_delim(file.path(data_dir,"gdc_tcga","tcga_RSEM_gene_fpkm"),delim = "\t") %>% 
# mutate_if(.predicate = is.double,.fun = (function(x) relu({(2^x) - 0.001 })))  %>% dplyr::rename( gene  = sample)



#split_genes <- strsplit(x = tcga_fpkm$gene,split = ".",fixed = TRUE)
#ens_gene_id_base <- unlist(lapply(split_genes,"[[",1))
#ens_gene_id_suffix <- unlist(lapply(split_genes,"[[",2))
#tcga_fpkm$gene <- ens_gene_id_base
#tcga_fpkm_preprocess is here
tcga_fpkm <- readRDS(file.path(data_dir,"gdc_tcga","tcga_fpkm_preprocess.rds"))


tcga_meta_1 <- read_delim(file.path(data_dir,"gdc_tcga","tissue_patientid.tsv"),delim = "\t")
tcga_meta_2 <- read_delim(file.path(data_dir,"gdc_tcga","sample_cancertype_sampletype.tsv"),delim = "\t")
tcga_meta <- full_join(tcga_meta_1,tcga_meta_2, by = "sample")
tcga_meta <- dplyr::select(tcga_meta,sample,"cancer type abbreviation","_primary_site","_PATIENT", sample_type ) %>%
  dplyr::rename(disease = "cancer type abbreviation"  , tissue = "_primary_site" ,patient_id = "_PATIENT")

tcga_meta$biopsy_tissue <- NA
tcga_meta$dataset <- "TCGA"


tcga_meta$biopsy_tissue <- ifelse( tcga_meta$sample_type == "Primary Tumor", tcga_meta$tissue, tcga_meta$biopsy_tissue)




  #joining
joined_meta <- bind_rows(met500_meta,lee_meta,tcga_meta)  #putting the two datasets together
joined_meta <- mutate(joined_meta,tissue = tolower(tissue))
joined_meta <- mutate(joined_meta,biopsy_tissue = tolower(biopsy_tissue))
#joined_meta <- mutate(joined_meta,biopsy_tissue = replace(biopsy_tissue,which(biopsy_tissue == "bone"),"bone_marrow"))


joined_fpkm_pre <- full_join(full_join(met500_fpkm,lee_fpkm,by = "gene"),tcga_fpkm,by = "gene")
rm("tcga_fpkm")
gc()
joined_fpkm_tbl <- transpose_df(joined_fpkm_pre)

gc()


joined_meta <- joined_meta[joined_meta$sample %in% joined_fpkm_tbl$sample,]
joined_fpkm_tbl <- arrange(joined_fpkm_tbl,sample)
joined_meta <- arrange(joined_meta,sample)

saveRDS(joined_fpkm_tbl,file.path(data_dir,"merged","joined_fpkm_tbl.rds"))
saveRDS(joined_meta,file.path(data_dir,"merged","joined_meta.rds"))
rm("joined_meta")


sample_code <- as.vector(select_if(.tbl = joined_fpkm_tbl, .predicate = is.character))  
sample_code <- sample_code$sample

joined_fpkm_tbl <- select_if(.tbl = joined_fpkm_tbl, .predicate = is.double)  
joined_fpkm_tbl <- as.matrix(joined_fpkm_tbl) #now it is actually a matrix
rownames(joined_fpkm_tbl) <- sample_code
saveRDS(joined_fpkm_tbl,file.path(data_dir,"merged","joined_fpkm_mat.rds"))




# lee_fpkm_mat <- data.matrix(lee_fpkm[,2:dim(lee_fpkm)[2]])
# rownames(lee_fpkm_mat) <- lee_fpkm$gene
# saveRDS(lee_fpkm_mat,file.path(data_dir,"Lee","lee_fpkm.rds"))
# saveRDS(lee_meta,file.path(data_dir,"Lee","lee_meta.csv"))
# 
# met500_fpkm_mat <- data.matrix(met500_fpkm[,2:dim(met500_fpkm)[2]])
# rownames(met500_fpkm_mat) <- met500_fpkm$gene
# saveRDS(met500_fpkm_mat,file.path(data_dir,"met500","met500_fpkm.rds"))
# saveRDS(met500_meta,file.path(data_dir,"met500","met500_meta.rds"))
# 
# 
# 
# tcga_fpkm_mat <- data.matrix(tcga_fpkm[,2:dim(tcga_fpkm)[2]])
# rownames(tcga_fpkm_mat) <- tcga_fpkm$gene
# saveRDS(tcga_fpkm_mat,file.path(data_dir,"gdc_tcga","tcga_fpkm.rds"))
# saveRDS(tcga_meta,file.path(data_dir,"gdc_tcga","tcga_meta.rds"))


#ditch the conversion to tpm here
# met500_fpkm_numeric_indx <- unlist(lapply(met500_fpkm,is.numeric))
# met500_fpkm_numeric <- met500_fpkm[,met500_fpkm_numeric_indx]
# lib_depths <- colSums(data.matrix(met500_fpkm_numeric))
# met500_tpm <- met500_fpkm
# met500_tpm[,met500_fpkm_numeric_indx] <- sweep(met500_fpkm_numeric,2,lib_depths,`/`)*100000
# met500_tpm$gene <-   sapply(strsplit(met500_tpm$gene,'.',fixed = TRUE),function(x) x[1])
#  
# 
# 
# 
# joined_tpm <- inner_join(met500_tpm, lee_quasi_tpm,by = "gene")
# joined_tpm_numeric_indx <- unlist(lapply(joined_tpm,is.numeric))
# joined_tpm_gene_names <- joined_tpm$gene
# 
# joined_tpm_numeric <- data.matrix(joined_tpm[,joined_tpm_numeric_indx])
# rownames(joined_tpm_numeric) <- joined_tpm_gene_names
# joined_tpm_numeric <- t(joined_tpm_numeric)
# joined_tpm_tibble <- as_tibble(joined_tpm_numeric,rownames = "sample")
# joined_tpm_tibble <- joined_tpm_tibble[match(joined_meta$sample,joined_tpm_tibble$sample),]



