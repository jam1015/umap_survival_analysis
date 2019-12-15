#converts to tpm from fpkm
#includes some normalization that I do in advance in othe rscrip
# add option to not display labels on output plot
#add option to loop through order layers
rm(list = ls())
setwd("C:\\Users/JAM526/pancan_fpkmuq/met_analysis/scripts/")
graphics.off()
library("Matrix.utils")
library("pracma")
library("dplyr")
library("stringr")
library("umap")
library("scales")
library("ggplot2")
library("ggthemes")
library("dbscan")
library("MASS")
library("RVAideMemoire")
library("made4")
library("plotly")
library("Rtsne")
library("randomForest")
library("org.Hs.eg.db")
library("GenomicFeatures")
library("dplyr")
library("biomaRt")
library("randomcoloR")
library("magrittr")
library("rlang")
#options

#how should we add noise



#don't include these things in the title or the legend 


base_dir <-  "C:/Users/JAM526/pancan_fpkmuq/met_analysis"
data_dir <- "data/merged"
fpkm_file <- "joined_fpkm_mat.rds"





#make variable to sort rows of DF for plotting


#reading the expression data

  quasi_tpm <- readRDS(file.path(base_dir,data_dir,fpkm_file))





  nrows <- size(quasi_tpm)[1]
  ncols <- size(quasi_tpm)[2]
  for (r in 1:nrows){
    wanted_row <- quasi_tpm[r,]
    
    rowsum <- sum(wanted_row)
    rm(list = "wanted_row")
    
    
    
    for (c in 1:ncols){
      to_modify <-   quasi_tpm[r,c]
      new_val <- (to_modify/rowsum)*10^6
      quasi_tpm[r,c] <- to_modify

      #if(mod(c,10000)==0){print(paste("col is",c))}
      
    }
    #if(mod(r,10000)==0){print(paste("col is",r))}
  
  }
  saveRDS(quasi_tpm,file.path(base_dir,data_dir,"joined_tpm_mat.rds"))
  print("saved")
  
