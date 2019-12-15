#plots things in a quadrant structure and getw their survival
#main correlation matrix plotting script.
rm_list <- ls()
rm_list <- rm_list[!(rm_list %in% c("expression_data","tpm_loaded"))]
rm(list = rm_list)
#reading the expression data
setwd("C:\\Users/JAM526/pancan_fpkmuq/met_analysis/scripts/")
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
library("TCGAbiolinks")
library("survminer")
library("survival")
library("ComplexHeatmap")
library("reshape2")
library("Hmisc")
library("EMT")
#options
add_grey <- TRUE
chi_sq_var <- "tissue"
color_selection <- quo(Region)
shape_selection <- quo(sample_type)
print_p_threshold <- 0.01
use_tpm <- TRUE
standardize <-  FALSE #FALSE#TRUE   #summing to one
spherize <- FALSE
draw_lines <- FALSE  #draw lines between matched samples
filter_dataset <- TRUE   #are we jsut seleecting a few datasets
dataset_filt <- "TCGA"
filter_sample_type <- TRUE

filter_disease <- TRUE

symbol_size <- 2
run_stats <- TRUE
subset_for_stats <- FALSE
ntile <- .4
save_image <- TRUE
axis_numbers <- TRUE
show_p_vals <- TRUE
#should we display every run? 
disease_filt_list <- c("BLCA")
show_one_iteration  <- FALSE
only_dup <- TRUE
div_increment <- 0.01
log_transform <- TRUE
divide_by_second <- TRUE
cluster_hm_list <- c(TRUE,FALSE)

sample_type_filt <- c("Primary Tumor","Solid Tissue Normal")
#ov
#hnsc
#lung squamous


#don't include these things in the title or the legend 
exclude_title_legend <- c("patient_id","sample","tissue","sample_tissue","biopsy_tissue")
colors_file = "C:/Users/JAM526/pancan_fpkmuq/met_analysis/data/color_palettes/pa200_set1_chroma23_lightness23.RDS"

base_dir <-  "C:/Users/JAM526/pancan_fpkmuq/met_analysis"
data_dir <- "data/merged"
data_dir_1 <- "data"
tpm_file <- "joined_tpm_mat.rds"
fpkm_file <- "joined_fpkm_mat.rds"
output_dir <- "output"
pathway_lists_dir <- "pathway_lists"

#the file names of the pathways we are using
pathway_files <- c(
  "cell_cycle.csv",
  "chol.csv",
  "fabo.csv",
  "glycolysis.csv",
  "hippo.csv",
  "mrp.csv",
  "mito_cyto_rp.csv",
  "myc.csv",
  "notch.csv",
  "pent_phos.csv",
  "pi3k.csv",
  "purine.csv",
  "pyrimidine.csv",
  "rp.csv",
  "tca.csv",
  "tgfb.csv",
  "tp53.csv",
  "wnt.csv",
  "nfe2l2keap1.csv",
  "nfe2l2_targets.csv",
  "nfe2l2.csv")


pathway_names <- c(
  "Cell Cycle",
  "Cholesterol Biosynthesis",
  "Fatty Acid Beta-Oxidation",
  "Glycolysis",
  "Hippo",
  "Mitochondrial Ribosomal Proteins",
  "Mitochondrial and Cytoplasmic Ribosomal Proteins",
  "Myc",
  "Notch",
  "Pentose Phosphate",
  "PI3-Kinase",
  "purine.csv",
  "Pyrimidine Biosynthesis",
  "Cytoplasmic Ribosomal Proteins",
  "TCA Cycle",
  "TGF-Beta",
  "Tp53",
  "Wnt",
  "NFE2L2 and KEAP1",
  "NFE2L2 Targets",
  "NFE2L2")
for(disease_filt in disease_filt_list){
  
  #################
  
  names(pathway_files) <- pathway_names
  pathway <- c("NFE2L2", "NFE2L2 Targets")
  pathway <- pathway_files[pathway]
  
  
  #WHAT ORDER DO THINGS APPEAER IN CANCER
  descendant <-     list( "Solid Tissue Normal" = "Primary Tumor",  "Primary Tumor" =   "Metastatic",   "Metastatic" = "None" )
  
  #GIVING BETTER NAMES TO PLOT THINGS BY
  plot_names <- list("sample" = "Sample" ,      "disease"  = "Disease" ,    "tissue" = "Tissue" ,"biopsy_tissue" = "Biopsy Tissue", "patient_id" = "Patient ID",    "sample_type" = "Sample Type",   "dataset" = "Dataset" )      
  
  
  
  #make variable to sort rows of DF for plotting
  source("functions/sample_order.r")
  source("functions/all_dup.r")
  #exponentiates and decrements
  source("functions/exp_dec.r")
  #collapses contingency table to be yes/no for one category 
  source("functions/contingency_collapse.r")
  #takes a patient, data frame with sample information, and a plot; outputs plot with added line between paired samples
  source("functions/plot_arrows.r")
  #function to load genetic data
  source("functions/load_gen_data.R")
  
  
  
  #loading the mapping between different gene names
  load(file.path(base_dir,data_dir_1,"genes_entrez_ensembl.rda"))
  
  
  expression_data <- load_gen_data()        
  
  #loading the metadata about each sample
  all_metadata <- readRDS(file.path(base_dir,data_dir,"joined_meta.rds")) #recall that a single primary can have multiple mets
  surv_data    <- readRDS(file.path(base_dir,data_dir,"gdc_survival.rds"))
  surv_data    <- surv_data[,c("sample","OS","OS.time")]
  all_metadata <- inner_join(all_metadata,surv_data,by = "sample")
  ##########################################looping through the pathways
  #if(show_one_iteration){ graphics.off()}
  metadata <- dplyr::filter(all_metadata,!is.na("disease"))
  if(filter_dataset){ metadata <- dplyr::filter(metadata,dataset %in% dataset_filt)  }
  if(filter_sample_type){  metadata <- dplyr::filter(metadata,sample_type %in% sample_type_filt)}
  if(filter_disease){  metadata <- dplyr::filter(metadata,disease %in% disease_filt)}
  
  
  #good if we're working on paired samples
  #doing a gene name conversion
  
  #reading the actual pathway
  tx1 <- levels(read.csv(file.path(base_dir,data_dir_1,pathway_lists_dir,pathway[1]),header = F)[,])
  tx2 <- levels(read.csv(file.path(base_dir,data_dir_1,pathway_lists_dir,pathway[2]),header = F)[,])
  tx <- c(tx1,tx2)
  
  
  entrez_to_ensembl   <- mapIds(org.Hs.eg.db,tx,"ENSEMBL","SYMBOL") #actual gene names
  entrez_to_gene_name <- names(entrez_to_ensembl) 
  #getting the expression of the genes in the pathway
  wanted_tpm <- expression_data[,entrez_to_ensembl,drop = FALSE]#matrix(expression_data[,entrez_to_ensembl1],dimnames = list(names(expression_data[,entrez_to_ensembl1]),entrez_to_ensembl1))
  #labeling with hgnc gene name  
  colnames(wanted_tpm) <- entrez_to_gene_name
  #filtering to values we want
  wanted_tpm <- wanted_tpm[metadata$sample,,drop = FALSE]
  #removing NA
  wanted_tpm <- wanted_tpm[, colSums(is.na(wanted_tpm)) == 0]
  #recall that a single primary can have multiple mets
  #matching up rownames of metadata and expression
  metadata_end <- dim(metadata)[2]
  #making the expression data into a tibble that can do it right 
  tpm_df <- dplyr::as_tibble(wanted_tpm)
  tpm_df$sample <- rownames(wanted_tpm)
  #merging metadata and expression data
  meta_tpm_joined <- merge(metadata,tpm_df,by = "sample")
  #which indices of the matrix have expression data
  expression_indices <- (metadata_end+1):NCOL(meta_tpm_joined)
  #deleting samples with no metadata
  has_gen <-  apply(meta_tpm_joined[,tx], 1, function(x) !any(is.na(x)))
  meta_tpm_joined <- meta_tpm_joined[has_gen,]
  scatter_matrix <- data.matrix(meta_tpm_joined[,tx,drop = FALSE])
  rownames(scatter_matrix) <- meta_tpm_joined$sample
  
  
  
  #linearly normalizing to one
  if(standardize){
    scatter_matrix <- sweep(scatter_matrix,1,rowSums(scatter_matrix),"/")
    scatter_matrix[is.na(scatter_matrix)] <- 0
  }
  
  
  #centering mean and projecting onto unit hypersphere
  if(spherize){
    centroid <- colMeans(scatter_matrix)
    scatter_matrix <- sweep(scatter_matrix,2,centroid,"-")
    radius <- sqrt(rowSums(scatter_matrix^2))
    scatter_matrix <- sweep(scatter_matrix,1,radius,"/")
  }
  meta_tpm_joined[,tx] <- scatter_matrix 
  meta_tpm_joined <- dplyr::filter(meta_tpm_joined,(sample_type == sample_type_filt[1] | sample_type == sample_type_filt[2]))
  if(only_dup){
    meta_tpm_joined <- tidyr::unite(meta_tpm_joined,sample_type_patient_code,sample_type,patient_id,remove = FALSE)
    meta_tpm_joined <- dplyr::filter(meta_tpm_joined,all_dup(meta_tpm_joined$patient_id))
    meta_tpm_joined <- dplyr::filter(meta_tpm_joined,!duplicated(sample_type_patient_code))
  }
  
  sample_one <- filter(meta_tpm_joined,sample_type == sample_type_filt[1]) %>% arrange(patient_id) %>% dplyr::filter(!duplicated(patient_id))
  sample_two <- filter(meta_tpm_joined,sample_type == sample_type_filt[2]) %>% arrange(patient_id) %>% dplyr::filter(!duplicated(patient_id))
  if(only_dup){
    sample_one <- semi_join(x = sample_one,y = sample_two, by = "patient_id")
    sample_two <- semi_join(x = sample_two,y = sample_one, by = "patient_id")
  }
  scatter_matrix_1 <- data.matrix(sample_one[,tx])+ div_increment
  scatter_matrix_2 <- data.matrix(sample_two[,tx])+ div_increment
  mean_of_second <- colMeans(scatter_matrix_2)
  if(divide_by_second){ 
    if(only_dup){
      scatter_matrix <- (scatter_matrix_1)/(scatter_matrix_2)
    } else{
      scatter_matrix <- sweep(x = scatter_matrix_1,MARGIN = 2,STATS = mean_of_second,FUN ="/" )}
  } else{
    scatter_matrix <- scatter_matrix_1 
  } 
  
  scatter_matrix_non_log <- scatter_matrix
  if(log_transform){scatter_matrix %<>% log()}
  
  
  
  sample_one[,tx] <- scatter_matrix
  actual_scatter_matrix <- as_tibble(sample_one)
  actual_scatter_matrix <- dplyr::select(actual_scatter_matrix,-sample)
  
  
  
  
  #making title by what is unique
  pway_name_string <- names(pathway)
  is_one <- apply(actual_scatter_matrix,2,function(x) {length(unique(x))} ) == 1
  title <-(actual_scatter_matrix[,is_one])
  title %<>% dplyr::select_if(~!all(is.na(.)))
  title <- dplyr::select(title, - !!((exclude_title_legend[exclude_title_legend %in% colnames(title)])))
  
  
  title_string <-  apply( title[1,] , 1 , paste , collapse = " " ) #pasting first line of title df with spaces
  
  
  #Did we apply these normalizations to data? if so include in title
  if(spherize) { title_string <- paste(title_string, "Spherized")}  
  if(standardize) { title_string <- paste(title_string, "Standardized")}
  if(divide_by_second){title_string <- paste(title_string, "Sample Type Normalized")}
  if(log_transform){title_string <- paste(title_string, "Log TxFormed")}
  if(only_dup){title_string <- paste(title_string,"Only Matched")}
  
  actual_scatter_matrix %<>% mutate_if(is.character,tools::toTitleCase)
  
  
  ### Making first plot
  to_plot <- colnames(actual_scatter_matrix[,tx])
  
  vectx1 <- scatter_matrix[,tx1,drop = FALSE] %>% rowMeans() #first pathway
  vectx2 <- scatter_matrix[,tx2,drop = FALSE] %>% rowMeans() #second pathway
  

  
  actual_scatter_matrix$meantx1 <- vectx1
  actual_scatter_matrix$meantx2 <- vectx2
  
  
  #### Making the color scales
  contrasting_colors <- tolower(readRDS(file = colors_file )) #reading from the colors file (a list of color hex values)
  if(add_grey){
    contrasting_colors[1] <- "#cccfcd"
  }
  color_column_string <- quo_name(color_selection)
  color_column <- actual_scatter_matrix[[color_column_string]] %>% factor() %>% levels()
  wanted_colors <- contrasting_colors[1:length(color_column)]
  names(wanted_colors) <- color_column
  col_scale <- scale_colour_manual(values = wanted_colors ,aesthetics = c("colour","fill")) 
  fill_scale <-  scale_fill_manual(values = wanted_colors ) 
  
  
  
  #this is the scatter 
  out_plot <- ggplot(actual_scatter_matrix,aes(x=meantx1,y = meantx2))+
    geom_point(size = symbol_size,aes(shape= !! shape_selection)) + 
    labs(title = title_string, x = names(pathway)[1], y = names(pathway)[2],shape = plot_names[[quo_name(shape_selection)]])+  col_scale +
    geom_smooth(mapping = aes(x=meantx1,y =meantx2),method = "lm",color = "red")+stat_cor(method = "pearson")
  print(out_plot)
  
  
  #saving image of scatter
  if(save_image){ #saving image
    fp_string <- file.path(base_dir,output_dir,paste(title_string,"_scatter",sep = ""))
    fp_string <- gsub(" ","_", fp_string,fixed = TRUE)                  
    fp_string <- paste(fp_string,".jpg",sep = "")   
    ggsave(fp_string,device = "jpeg")
  }
  
  for(cluster_hm in cluster_hm_list){
    #making correlation heatmap
    if(cluster_hm){title_string <- paste(title_string,"Clustered")} else {title_string <- gsub(" Clustered","",title_string) }
    cor <- rcorr(data.matrix(actual_scatter_matrix[,tx]))  #getting rid of the averages
    r_all <- cor$r
    r <- r_all
    r_p <- cor$P   
    r[r_p > 0.05]  <- NA
    r_for_test <- r[upper.tri(r)]
    nas <- sum(is.na(r_for_test))
    positive <- sum(r_for_test[!is.na(r_for_test)] >0)
    negative <- sum(r_for_test[!is.na(r_for_test)]<=0)
    
    denom <- length(r_for_test)
    expect <- (positive+negative)/2
    probs <- c(nas,expect,expect)/denom
    
    binom <- RVAideMemoire::multinomial.test(c(nas,positive,negative),probs)
    
    col<- colorRampPalette(c("blue", "white", "red"))(20)
    hm <- Heatmap(cor$r,na_col = "black",name = "Correlation",
                  cell_fun = function(j, i, x, y, width, height, fill= "black"){
                    if(is.na(r[i,j])){ grid.rect(x = x,y = y,width = width, height = height,gp = gpar(col ="black", fill ="black"),just = "centre")
                    }
                  },
                  cluster_rows = cluster_hm,cluster_columns = cluster_hm,column_title = title_string)
    
    
    draw(hm)
    fp_string <- file.path(base_dir,output_dir,paste(title_string,sep = ""))
    fp_string <- gsub(" ","_", fp_string,fixed = TRUE)                  
    fp_string <- paste(fp_string,".png",sep = "")   
    png(filename = fp_string,width = 3000,height = 3000,res = 300,pointsize = 9)
    draw(hm)
    dev.off()
    
    
    
    
  }
  fp_string <- file.path(base_dir,output_dir,"binomial_test")
  fp_string <- gsub(" ","_", fp_string,fixed = TRUE) 
  fp_string <- paste(fp_string,"_binom_test.txt",sep = "")  
  fp_string %<>% gsub("_Clustered","",.)
  
  
  
  
  
  
  txt <- paste(title_string,"P value",as.character(binom$p.value),"Probability given number of NAS",positive/(positive+negative),sep = "\n")
  write(txt,file= fp_string,append = TRUE)
  
  
  #printing out the t test for nfe2l2 
  tt <- t.test(x = actual_scatter_matrix$meantx1)
  fp_string <- file.path(base_dir,output_dir)
  fp_string <- gsub(" ","_", fp_string,fixed = TRUE) 
  fp_string <- paste(fp_string,"/ttest_1.txt",sep = "")  
  fp_string %<>% gsub("_Clustered","",.)
  
  
  if(log_transform){meen <- as.character(exp(tt$estimate))[1]} else{meen <- as.character((tt$estimate))[1]}
  txt <- paste("----------------------------",title_string,tx1[1],dim(scatter_matrix)[1],dim(scatter_matrix)[2],"P value",as.character(tt$p.value),"mean",meen,sep = "\n")
  write(txt,file = fp_string,append = TRUE)
  
  
  #printing out the t test for nfe2l2 targets
  tt <- t.test(x = actual_scatter_matrix$meantx2)
  fp_string <- file.path(base_dir,output_dir)
  fp_string <- gsub(" ","_", fp_string,fixed = TRUE) 
  fp_string <- paste(fp_string,"/ttest_2.txt",sep = "")  
  fp_string %<>% gsub("_Clustered","",.)
  print(binom$p.value)
  print(binom$conf.int)
  
  if(log_transform){meen <- as.character(exp(tt$estimate))[1]} else{meen <- as.character((tt$estimate))[1]}
  
  txt <- paste("----------------------------",title_string,tx2[1],dim(scatter_matrix)[1],dim(scatter_matrix)[2],"P value",as.character(tt$p.value),"mean",meen,sep = "\n")
  write(txt,file = fp_string,append = TRUE)
  
  
  
  
}
