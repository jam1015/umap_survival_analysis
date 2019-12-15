
#this is broken right now.  it only works for tcga data.  other data needs to be appropriately wrangeled.  if downloaded from github file pathws should be changeed for the local computer.


#add loop for output formats
dummy  <- runif(5)
rm_list <- ls()
rm_list <- rm_list[!(rm_list %in% c("expression_data","tpm_loaded"))]
rm(list = rm_list)
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
library(survival)
library(survminer)
library(survMisc)
library(rsample)      # data splitting 
library(randomForest) # basic implementation
library(ranger)       # a faster implementation of randomForest
library(caret)        # an aggregator package for performing many machine learning models
library(h2o)     

library("grid")
library("ggplotify")
library("gridExtra")
#options
chi_sq_var <- "tissue"
color_selection <- quo(sample_type)
color_selection_2 <- quo(Cluster)
shape_selection <- NULL#quo(disease)
print_p_threshold <- 0.01
use_tpm <- TRUE
standardize <-  TRUE #FALSE#TRUE   #summing to one
spherize <- TRUE
draw_lines <- FALSE  #draw lines between matched samples
filter_dataset <- TRUE   #are we jsut seleecting a few datasets
dataset_filt <- "TCGA"
filter_sample_type <- TRUE
#sample_type_filt <- "Primary Blood Derived Cancer - Peripheral Blood"
sample_type_filt <- "Primary Tumor"
filter_disease <- TRUE

symbol_size <- 2
num_trees <- 500

subset_for_stats <- FALSE
stats_filter_string <- "Primary Tumor" #string for selecting for subtype
save_image <- TRUE
axis_numbers <- FALSE
show_p_vals <- TRUE
#should we display every run? 
show_one_iteration  <- FALSE

tpm_file <- "joined_tpm_mat.rds"
fpkm_file <- "joined_fpkm_mat.rds"

fisher <- TRUE
font_size <- 14

save_data <- TRUE

#how should we add noise
jitter_cluster <- FALSE
jitter_plot <- FALSE
jitter_factor <- NULL
jitter_amount <- .005


#umap parameters
custom.config <- umap.defaults
custom.config$n_epochs <- 2000
custom.config$min_dist <- 10^(-95)
custom.config$n_neighbors <- 10
custom.config$verbose <- TRUE
hdbscan_min_pts <- 25
run_stats <- TRUE
survival_analysis <- TRUE
contingency_analysis <- FALSE
survival_zero <- FALSE
run_rfc <- TRUE
used_pathway_files <- c("tcell_checkpoint.csv")
disease_filt <- "KIRP"
out_fig_width <- 6.5
out_fig_height <- 6.5
out_fig_units <- "in"
filter_noise <- TRUE



#don't include these things in the title or the legend 
exclude_title_legend <- c("patient_id","sample","tissue","sample_tissue","biopsy_tissue")
colors_file <- "C:/Users/JAM526/pancan_fpkmuq/met_analysis/data/color_palettes/pa200_set1_chroma23_lightness23.RDS"
colors_file_2 <- "C:/Users/JAM526/pancan_fpkmuq/met_analysis/data/color_palettes/new_pastel.RDS" 

base_dir <-  "C:/Users/JAM526/pancan_fpkmuq/met_analysis"
data_dir <- "data/merged"
data_dir_1 <- "data"
output_dir <- "output"
pathway_lists_dir <- "pathway_lists"

#the file names of the pathways we are using
pathway_files <- c(
  "cell_cycle.csv",#
  "chol.csv",#
  "fabo.csv",
  "glycolysis.csv",#
  "hippo.csv",
  "mrp.csv",#
  "mito_cyto_rp.csv",#
  "myc.csv",#
  "notch.csv",#
  "pent_phos.csv",#
  "pi3k.csv",#
  "purine.csv",#
  "pyrimidine.csv",#
  "rp.csv",#
  "tca.csv",#
  "tgfb.csv",#
  "tp53.csv",#
  "wnt.csv",#
  "ifng.csv",#
  "nfkb_tnf.csv",#
  "tcell_checkpoint.csv")#

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
  "Purine Biosynthesis",
  "Pyrimidine Biosynthesis",
  "Cytoplasmic Ribosomal Proteins",
  "TCA Cycle",
  "TGF-Beta",
  "Tp53",
  "Wnt",
  "IFN Gamma",
  "NF-kB TNF",
  "T-Cell and Checkpoint")

save_list <- list()  #making an empty list of things to save

names(pathway_files) <- pathway_names

pathway_files <- pathway_files[pathway_files %in% used_pathway_files]


#WHAT ORDER DO THINGS APPEAER IN CANCER
descendant <-     list( "Solid Tissue Normal" = "Primary Tumor",  "Primary Tumor" =   "Metastatic",   "Metastatic" = "None" )

#GIVING BETTER NAMES TO PLOT THINGS BY
plot_names <- list("sample" = "Sample" ,      "disease"  = "Disease" ,    "tissue" = "Tissue" ,"biopsy_tissue" = "Biopsy Tissue", "patient_id" = "Patient ID",    "sample_type" = "Sample Type",   "dataset" = "Dataset" )      


ggcolor_hue_minus <- function(n) {
  hues = seq(15, 375, length = n + 1)
  full <- hcl(h = hues, l = 65, c = 100)[1:n]
  minus_one <-full[2:numel(full)]
  return(minus_one)
}
#make variable to sort rows of DF for plotting
source("functions/sample_order.r")
source("functions/all_dup.r")
#exponentiates and decrements
source("functions/exp_dec.r")
#collapses contingency table to be yes/no for one category 
source("functions/contingency_collapse.r")
#takes a patient, data frame with sample information, and a plot; outputs plot with added line between paired samples
source("functions/plot_arrows.r")
source("functions/chisq_res_p.r")

#function to load genetic data
source("functions/load_gen_data.R")


expression_data <- load_gen_data()

#loading the mapping between different gene names
load(file.path(base_dir,data_dir_1,"genes_entrez_ensembl.rda"))

#loading the metadata about each sample
all_metadata <-  readRDS(file.path(base_dir,data_dir,"joined_meta.rds")) #recall that a single primary can have multiple mets
surv_data <- readRDS(file.path(base_dir,data_dir,"gdc_survival.rds"))
surv_data <- surv_data[,c("sample","OS","OS.time")]
all_metadata <- inner_join(all_metadata,surv_data,by = "sample")


##########################################looping through the pathways
for(pathway in pathway_files){
  
  if(show_one_iteration){ graphics.off()}
  
  metadata <- dplyr::filter(all_metadata,!is.na("disease"))
  
  
  if(filter_dataset){ metadata <- dplyr::filter(metadata,dataset %in% dataset_filt)  }
  if(filter_disease){  metadata <- dplyr::filter(metadata,disease %in% disease_filt)}
  if(filter_sample_type){  metadata <- dplyr::filter(metadata,sample_type %in% sample_type_filt)}
  
  
  
  #doing a gene name conversion
  genes_entrez_ensembl_filtered_ordered = genes_entrez_ensembl[match(colnames(expression_data),genes_entrez_ensembl$ensembl_gene_id),]
  
  #reading the actual pathway
  tx <- read.csv(file.path(base_dir,data_dir_1,pathway_lists_dir,pathway),header = F)$V1
  transcript_selector <- match(tx,genes_entrez_ensembl_filtered_ordered$hgnc_symbol)
  #Selecting based on the entrez name rather than logical vector to be safe. deleting NA's from this. 
  transcript_selector <- transcript_selector[!is.na(transcript_selector)] 
  entrez_to_ensembl   <- mapIds(org.Hs.eg.db,as.character(tx),"ENSEMBL","SYMBOL") #actual gene names
  entrez_to_gene_name <- names(entrez_to_ensembl)                  #genes_entrez_ensembl_filtered_ordered[transcript_selector,   "hgnc_symbol"]  #these are the actual gene names
  #getting the expression of the genes in the pathway
  wanted_tpm <- expression_data[,entrez_to_ensembl]
  #labeling with hgnc gene bane
  colnames(wanted_tpm) <- entrez_to_gene_name
  #filtering to values we want
  wanted_tpm <- wanted_tpm[metadata$sample,]
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
  has_gen <-  apply(meta_tpm_joined[,expression_indices], 1, function(x) !any(is.na(x)))
  meta_tpm_joined <- meta_tpm_joined[has_gen,]
  #what we actually run umap on 
  umap_matrix <- data.matrix(meta_tpm_joined[,expression_indices])
  #linearly normalizing to one
  if(standardize){
    umap_matrix <- sweep(umap_matrix,1,rowSums(umap_matrix),"/")
    umap_matrix[is.na(umap_matrix)] <- 0
  }
  
  
  #centering mean and projecting onto unit hypersphere
  if(spherize){
    centroid <- colMeans(umap_matrix)
    umap_matrix <- sweep(umap_matrix,2,centroid,"-")
    radius <- sqrt(rowSums(umap_matrix^2))
    umap_matrix <- sweep(umap_matrix,1,radius,"/")
  }
  rownames(umap_matrix) <- meta_tpm_joined$sample
  #running umap
  embedding <- umap(umap_matrix,config = custom.config, method = "umap-learn",verbose = TRUE)
  layout_df <- data.frame(embedding$layout)
  
  #joining umap layout with meta_tpm joined.  important to not reorder meta_tpm_joined before running umap
  joined_umap <- bind_cols(meta_tpm_joined,layout_df)  
  
  #RUNNING CLUSTERING AND APPLYTING JITTER IF NECESSARY
  to_clust <- cbind(joined_umap$X1,joined_umap$X2)
  colnames(to_clust) <- c("X1","X2")
  if(jitter_cluster){
    to_clust %<>% jitter(amount = jitter_amount,factor = jitter_factor)
  } else{ }
  
  joined_umap$Cluster <- factor(hdbscan(to_clust,minPts = hdbscan_min_pts)$cluster)  #running hdbscan
  
  
  #making title
  pway_name_string <- names(pathway_files)[pathway_files %in% pathway] 
  is_one <- (apply(joined_umap,2,function(x) {length(unique(x))} ) == 1) & !(colnames(joined_umap) %in% levels(tx) )
  
  
  
  title <-(joined_umap[,is_one])
  title %<>% dplyr::select_if(~!all(is.na(.)))
  title <- dplyr::select(title, - !!((exclude_title_legend[exclude_title_legend %in% colnames(title)])))
  
  
  title_string <- paste("UMAP", 
                        pway_name_string, 
                        apply( title[1,] , 1 , paste , collapse = " " )) #pasting first line of title df with spaces
  
  
  #Did we apply these normalizations to data? if so include in title
  if(spherize) { title_string <- paste(title_string, "Spherized")}  
  if(standardize) { title_string <- paste(title_string, "Standardized")}
  
  
  
  joined_umap <- as_tibble(joined_umap)
  joined_umap %<>% mutate_if(is.character,tools::toTitleCase)
  ##### APPLYING JITTER TO PLOT IF WE SPECIFIED
  if(jitter_plot){joined_umap$X1 <- to_clust[,"X1"]
  joined_umap$X2 <- to_clust[,"X2"]
  title_string <- paste(title_string, "Jittered")}else{}
  
  
  caption_string <- paste0("Parameters: n-Neighbors = ",custom.config$n_neighbors ,", Min. Dist. = ",custom.config$min_dist,", n-Epochs = ",custom.config$n_epochs)
  #### Making the color scales
  contrasting_colors <- tolower(readRDS(file = colors_file )) #reading from the colors file (a list of color hex values)
  color_column_string <- quo_name(color_selection)
  color_column <- joined_umap[[color_column_string]] %>% factor() %>% levels()
  wanted_colors <- contrasting_colors[1:length(color_column)]
  names(wanted_colors) <- color_column
  col_scale <- scale_colour_manual(values = wanted_colors ,aesthetics = c("colour","fill")) 
  
  
  
  contrasting_colors_2 <- tolower(readRDS(file = colors_file_2 )) #reading from the colors file (a list of color hex values)
 names(contrasting_colors_2) <- (1:length(contrasting_colors_2)) -1
   color_column_string_2 <- quo_name(color_selection_2)
  color_column_2 <- joined_umap[[color_column_string_2]] %>% factor() %>% levels()
  wanted_colors_2 <- contrasting_colors_2[1:length(color_column_2)]
  names(wanted_colors_2) <- color_column_2
  col_scale_2 <- scale_colour_manual(values = contrasting_colors_2 ,aesthetics = c("colour","fill")) 
  
  
  
  ###   SORTING THE PLOT DATAFRAME SO THAT THE METASTASES ARE ON TOP 
  joined_umap$to_sort <- lapply(FUN = sample_order, X = joined_umap$sample_type) %>% unlist()
  joined_umap %<>% arrange(to_sort)
  
  
  ##################### 
  #Makign first plot
  out_plot <- ggplot(joined_umap,aes(x=X1,y=X2))+
    geom_point(size = symbol_size,aes(color= !! color_selection, shape = !! shape_selection)) + 
    labs(title = title_string, y="UMAP 2", x="UMAP 1",caption = caption_string,color = plot_names[[quo_name(color_selection)]],shape = plot_names[[quo_name(shape_selection)]])+  col_scale+  theme(aspect.ratio=1,text = element_text(size=font_size))+theme(plot.title = element_text(hjust = 0.5))
  
  if(!axis_numbers){
    out_plot <- out_plot + theme(axis.text = element_blank(),axis.ticks = element_blank()) 
  }
  umap_scatter <- out_plot
  
  
  first_out_plot <- out_plot
  save_list <- c(save_list,list(first_scatter = first_out_plot))
  
  #THIS DRAWS THE LINES BY PUTTING THINGS ON THE OUTPUT PLOT 
  if(draw_lines){
    
    
    for (d in 1:dim(joined_umap)[1]){ #looping through output and adding arrows
      st_in <- joined_umap$sample[d]
      out_plot <- plot_arrows(st_in,joined_umap,out_plot)
    }
    #priting the version with lines
    umap_scatter_lines <- out_plot
    
    first_out_plot_lines <- out_plot
    save_list <- c(save_list,list(first_scatter_lines = first_out_plot_lines))
    
  } else{}
  
  
  
  if(filter_noise){ joined_umap_scatter <- filter(joined_umap,Cluster != 0)} else { joined_umap_scatter <- joined_umap}
  
  
  
  # PLOTTING THE UMAP OUTPUT COLORED BY CLUSTER
  title_string2 <- paste("HDBSCAN of" ,disease_filt,pway_name_string)
  
  # shape = plot_names[[quo_name(shape_selection)]]
  out_plot <- ggplot(joined_umap_scatter,aes(x=X1,y=X2)) + 
    geom_point(size = symbol_size,aes(color=Cluster,shape = !!shape_selection)) + 
    labs(title = title_string2, y="UMAP 2", x="UMAP 1",caption = caption_string,shape = plot_names[[quo_name(shape_selection)]])+ 
    theme(aspect.ratio=1,text = element_text(size=font_size)) + theme(plot.title = element_text(hjust = 0.5)) + col_scale_2
  
  
  
  # 
  # if(!any(joined_umap_scatter$Cluster == 0)){
  #   pal = contrasting_colors_2[2:(numel(levels(joined_umap_scatter$Cluster))+1)]
  #   col_scale_base <- scale_colour_manual(values = pal ,aesthetics = c("colour","fill")) 
  #   out_plot <- out_plot + col_scale_base
  # } else{     pal = contrasting_colors_2[1:(numel(levels(joined_umap_scatter$Cluster)))]
  # col_scale_base <- scale_colour_manual(values = pal ,aesthetics = c("colour","fill")) 
  # out_plot <- out_plot + col_scale_base}
  # 
  
  
  
  
  if(!axis_numbers){
    out_plot <- out_plot + theme(axis.text = element_blank(),axis.ticks = element_blank()) 
  }
  
  umap_scatter_cluster <- out_plot;
  save_list <- c(save_list,list(umap_scatter_cluster = umap_scatter_cluster))
  
  ######STARTING THE PART WHERE WE RUN STATS
  if(run_stats){
    #the subset of tumors
    
    if(subset_for_stats){    
      subset_umap <- joined_umap[(joined_umap$sample_type) %in% stats_filter_string,]
    } else{ subset_umap <- joined_umap}
    
    
    
    #a contingency table of cluster vs. met site
    tbl <- table(as.character(subset_umap$Cluster),unlist(lapply(subset_umap[[chi_sq_var]],as.character)))
    
    #converting the above table to matrix
    tbl_mat <- as.matrix(tbl)
    
    #running chi square on the table
    cs <- chisq.test(tbl, correct = TRUE)
    save_list <- c(save_list,list(global_chi_square = cs))
    
    #running post-hoc fisher test on the frequency table 
    if(contingency_analysis){
      
      if(fisher) {
        post_hoc <- matrix(NA,dim(tbl_mat)[1],dim(tbl_mat)[2])
        rownames(post_hoc) <- rownames(tbl_mat)
        colnames(post_hoc) <- colnames(tbl_mat)
        
        
        for (e in 1:dim(tbl_mat)[1]){
          for (f in 1:dim(tbl_mat)[2]){
            collapsed <- contingency_collapse(tbl_mat,e,f)
            
            ft_out <- fisher.test(collapsed)
            post_hoc[e,f] <- ft_out$p.value
            
            
          }
        }
        post_hoc_fisher <- post_hoc
        save_list <- c(save_list,list(post_hoc_fisher = post_hoc_fisher))
      } else {
        post_hoc_res <- cs$stdres
        rownames(post_hoc_res) <- rownames(tbl_mat)
        colnames(post_hoc_res) <- colnames(tbl_mat)
        post_hoc <- pchisq((post_hoc_res)^2,df = 1,lower.tail = FALSE)
        post_hoc_res
        save_list <- c(save_list,list(post_hoc_residuals = post_hoc_res))
      }
      
      
      
      
      #############################!!!!!!!!!!!!!!!!!!!!!!!
      freq_tab <- dplyr::count(subset_umap,(!!sym(chi_sq_var)),Cluster)
      
      #APPENDING THE POST-HOC P VALUE TO DATA FRAME
      #############################!!!!!!!!!!!!!!!!!!!!!!!
      freq_tab$fisher_exact_p <- post_hoc[cbind( as.character(freq_tab$Cluster),as.character(freq_tab[[chi_sq_var]]))]
      fish_sig <- freq_tab$fisher_exact_p < print_p_threshold
      
      freq_tab %<>% mutate(fisher_exact_p = paste("P =",formatC(.data$fisher_exact_p, format = "e", digits = 2)))
      freq_tab$fisher_exact_p[!fish_sig] <- ""
      freq_tab %<>% mutate_if(is.character,tools::toTitleCase)
      save_list <- c(save_list,list(frequency_table = freq_tab))
      
      #MAKING THE bar PLOT
      title_string3 <- paste0(" Cluster Composition of ",title_string2)
      cs_format <- formatC(cs$p.value, format = "e", digits = 2)
      caption_string <- paste0("Chi-sq P = ",cs_format,"; P < ",as.character(print_p_threshold) ," displayed above")
      
      cluster_comp_bp <- ggplot(freq_tab,aes(x = Cluster,y = n,label = fisher_exact_p))+fill_scale+
        labs(title = c(title_string3),x = "Cluster",y = "Frequency",caption = caption_string, fill =  plot_names[[quo_name(color_selection)]] )+ 
        geom_bar(stat = "identity")  +aes(fill =  !!sym(color_column_string))
      
      
      
      save_list <- c(save_list,list(bar_plot = cluster_comp_bp))
      
      
      if(show_p_vals){
        cluster_comp_bp <- cluster_comp_bp + geom_text(size = 3, position = position_stack(vjust = 0.5))
        cluster_comp_bp_pvals <- cluster_comp_bp
        save_list <- c(save_list,list(bar_plot_pvals = cluster_comp_bp_pvals))
        
        
        
      }
    }
    
    
    
    if(survival_analysis){
      
      
      names(contrasting_colors_2) <- paste0("Cluster","=",names(contrasting_colors_2))
      
      
      title_string_surv <- paste0("Survival of ",title_string," ", 
                                  pway_name_string," ", 
                                  apply( title[1,] , 1 , paste , collapse = " " )) 
      
      if(!survival_zero){
        subset_umap <- dplyr::filter(subset_umap,Cluster != "0")
        subset_umap$Cluster <- droplevels(subset_umap$Cluster)
        }
      
      
      color_selector <- levels(subset_umap$Cluster) %>% strtoi()  %>% sort() 
      color_selector <- color_selector+1
      contrasting_colors_2 <- contrasting_colors_2[color_selector]
     
      
      
       fit <- survfit(Surv(OS.time,OS) ~ Cluster,data = subset_umap)
      save_list <- c(save_list,list(full_survival_fit = fit))
      
      surv_plot_obj <- ggsurvplot(fit,data = subset_umap, pval = TRUE,title = title_string_surv,ggtheme = theme_gray(),ylim = c(0,1),xlim = c(0,max(subset_umap$OS.time,na.rm = TRUE)),palette = contrasting_colors_2,risk.table = TRUE,xlab = "Time (Days)",font.x = font_size,font.y =font_size,font.xtickslab = font_size,font.ytickslab = font_size,font.main = font_size,font.legend = font_size,censor.shape  = 124)
      
      save_list <- c(save_list,list(survival_plot_object = surv_plot_obj))
      
      
      surv_plot_obj$plot <- surv_plot_obj$plot + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))+ theme(aspect.ratio=1)+theme(plot.title = element_text(hjust = 0.5))
      surv_plot_obj
      surv_plot <- surv_plot_obj$plot
      pw <- pairwise_survdiff(Surv(OS.time,OS) ~ Cluster,data = subset_umap,p.adjust.method = "none")
      pw_p <- formatC(pw$p.value, format = "e", digits = 2)
      
      save_list <- c(save_list,list(pairwise_surv = pw))
      save_list <- c(save_list,list(subset_umap = subset_umap))
      
      pw_tab <- tableGrob(pw_p,theme = ttheme_minimal())
      
      g1 <-  arrangeGrob(grobs  = list(surv_plot_obj$plot,
                                       surv_plot_obj$table,
                                       textGrob("Pairwise P-Values:"),
                                       pw_tab),heights = c(10,3,3))
      
      lay<- rbind(1,1,1,1,1,1,1,2,3,3,3)
      
      
      g1 <-  grid.arrange(grobs =  list(surv_plot_obj$plot,
                                        textGrob("Pairwise P-Values:"),
                                        pw_tab),layout_matrix = lay)
      
      
      
      
      
      
      
      
      if(run_rfc) {
        gen_df <- as_tibble(umap_matrix, rownames = "sample")
        info_df <- subset_umap[,c("sample","Cluster")]
        rfc_df <- left_join(gen_df,info_df,by = "sample",keep = FALSE) %>% dplyr::select(-sample) %>% dplyr::filter(Cluster != 0)
        
        hyper_grid <- expand.grid(
          mtry       = seq(20, 30, by = 2),
          node_size  = seq(3, 9, by = 2),
          sampe_size = c(.55, .632, .70, .80),
          OOB_RMSE   = 0
        )
        
        for(i in 1:nrow(hyper_grid)) {
          
          
          
          
          
          
          # train model
          model <- ranger(
            formula         = Cluster ~ ., 
            data            = rfc_df, 
            num.trees       = 500,
            mtry            = hyper_grid$mtry[i],
            min.node.size   = hyper_grid$node_size[i],
            sample.fraction = hyper_grid$sampe_size[i],
            seed            = 42069
          )
          
          # add OOB error to grid
          hyper_grid$OOB_RMSE[i] <- sqrt(model$prediction.error)
        }
        
        RMSE <- hyper_grid$OOB_RMSE
        min_error <- RMSE == min(RMSE)
        params <- hyper_grid[min_error,][1,]
        
        
        OOB_RMSE <- vector(mode = "numeric", length = 100)
        
        for(i in seq_along(OOB_RMSE)) {
          
          rfc_out <- ranger(
            formula         = Cluster ~ ., 
            data            = rfc_df, 
            num.trees       = num_trees,
            mtry            = params$mtry,
            min.node.size   = params$node_size,
            sample.fraction = params$sampe_size,
            importance      = 'impurity'
          )
          
          OOB_RMSE[i] <- sqrt(rfc_out$prediction.error)
        }
        
        
        
        
        rfc_plot_data <-     rfc_out$variable.importance %>% 
          tidy() %>%
          dplyr::arrange(desc(x)) %>%
          dplyr::top_n(15) 
        
        colnames(rfc_plot_data) <- c("Gene","Importance")
        rfc_plot_data %<>% arrange(-Importance)
        caption_string <- paste0("Parameters: n-trees = ",num_trees ,", Mtry  = ",params$mtry  ,", Min Node Size = ",params$node_size, ", Sample Fraction = ",params$sampe_size )
        rfc_plot <-   ggplot(rfc_plot_data,aes(x = reorder(Gene, Importance), y = Importance)) +
          geom_col() +
          coord_flip()  + labs(title = paste(title_string,"Random Forest Gene Importance"),y = "Predictive Importance", x = "Gene",
                               caption = caption_string)+theme(aspect.ratio=1,text = element_text(size=font_size))+theme(plot.title = element_text(hjust = 0.5))
        
        
        save_list <- c(save_list,list(rfc_plot = rfc_plot))
        save_list <- c(save_list,list(rfc_model = rfc_out))
        
      }

      
    }
    
    
    if(save_image){
      plots <- quos(  umap_scatter,
                      umap_scatter_lines,
                      umap_scatter_cluster,
                      surv_plot,
                      cluster_comp_bp,
                      rfc_plot)
      for (w in plots){
        var_name <- w %>% quo_name()
        if(var_name %>%  exists()) {
          
          
          
          print(eval_tidy(w))
          fp_string <- file.path(base_dir,output_dir,paste(title_string,"_",var_name,sep = ""))            
          fp_string <- gsub(" ","_", fp_string,fixed = TRUE)
          fp_string <- paste(fp_string,".png",sep = "")
          ggsave(fp_string,device = "png", height = out_fig_height,units = out_fig_units)
          
          fp_string <- file.path(base_dir,output_dir,paste(title_string,"_",var_name,sep = ""))            
          fp_string <- gsub(" ","_", fp_string,fixed = TRUE)
          fp_string <- paste(fp_string,".svg",sep = "")
          ggsave(fp_string,device = "svg", height = out_fig_height,units = out_fig_units)
          
          
        }
        
      }
      
      
    }
    saved_once <- TRUE
    
    
    if(save_data){ fp_string <- file.path(base_dir,output_dir,paste(title_string,"_data_accumulation",sep = ""))
    fp_string <- gsub(" ","_", fp_string,fixed = TRUE)
    fp_string <- paste(fp_string,".rds",sep = "")
    saveRDS(object = save_list,file = fp_string)} 
    
  }
  
}