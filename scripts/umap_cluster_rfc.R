rm(list = ls())
graphics.off()
library("org.Hs.eg.db")
library("GenomicFeatures")
library("dplyr")
library("biomaRt")
library(rsample)      # data splitting 
library(randomForest) # basic implementation
library(ranger)       # a faster implementation of randomForest
library(caret)        # an aggregator package for performing many machine learning models
library(h2o)          # an extremely fast java-based platform


in_directory <- c("C:/Users/JAM526/pancan_fpkmuq/met_analysis/output/immune_pancan_tsne_has_joined/")
wanted_files <- in_directory %>% dir() %>% grepl("data_accumulation.rds",.,fixed = TRUE) %>% "["(dir(in_directory),.)




for (f in wanted_files){
 fp <- file.path(in_directory,f)
 all_data <- readRDS(fp)
wanted_df <- all_data$joined_umap
 all_data %>% names() %>% print()
 entrez_to_ensembl   <- mapIds(org.Hs.eg.db,colnames(wanted_df),"ENSEMBL","SYMBOL") 
 genes <- entrez_to_ensembl[!is.na(entrez_to_ensembl)] %>% names()
 rfc_df <- wanted_df[,c(genes,"Cluster")]

 
 OOB_RMSE <- vector(mode = "numeric", length = 100)
 
 for(i in seq_along(OOB_RMSE)) {
   
   optimal_ranger <- ranger(
     formula         = Cluster ~ ., 
     data            = rfc_df, 
     num.trees       = 500,
     mtry            = 24,
     min.node.size   = 5,
     sample.fraction = .8,
     importance      = 'impurity'
   )
   
   OOB_RMSE[i] <- sqrt(optimal_ranger$prediction.error)
 }
 
 hist(OOB_RMSE, breaks = 20)
 
 
 optimal_ranger$variable.importance %>% 
   tidy() %>%
   dplyr::arrange(desc(x)) %>%
   dplyr::top_n(25) %>%
   ggplot(aes(reorder(names, x), x)) +
   geom_col() +
   coord_flip() +
   ggtitle("Top 25 important variables")
 
 
 
 }