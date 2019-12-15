#includes some normalization that I do in advance in othe rscrip
# add option to not display labels on output plot
rm_list <- ls()
rm_list <- rm_list[!(rm_list %in% "quasi_tpm")]
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

#to do: get rid of use of quo below, and use !! sym(.) instead
#options
add_white <- TRUE
chi_sq_var <- "hpv_status_by_p16_testing"
color_selection <- quo(hpv_status_by_p16_testing)
shape_selection <- quo(sample_type)
print_p_threshold <- 0.01
tpm_convert <- FALSE
standardize <-  TRUE #FALSE#TRUE   #summing to one
spherize <- FALSE
draw_lines <- FALSE  #draw lines between matched samples
filter_dataset <- TRUE   #are we just seleecting a few datasets
dataset_filt <- "TCGA"
filter_sample_type <- TRUE
sample_type_filt <- "Primary Tumor"
filter_disease <- TRUE
disease_filt <- "HNSC"
symbol_size <- 2
run_stats <- TRUE
subset_for_stats <- TRUE
filter_exp <- quo(sample_type == "Primary Tumor" & hpv_status_by_p16_testing != "Missing" & Cluster != 0 )

save_image <- FALSE
axis_numbers <- TRUE
show_p_vals <- TRUE
#should we display every run? 
show_one_iteration  <- FALSE
hdbscan_min_pts <- 20
remove_missing <- FALSE


#how should we add noise
jitter_cluster <- FALSE
jitter_plot <- FALSE
jitter_factor <- NULL
jitter_amount <- 1


#umap parameters
custom.config <- umap.defaults
custom.config$n_epochs <- 1000
custom.config$min_dist <- .000000001
custom.config$n_neighbors <- 10
custom.config$verbose <- TRUE
umap_method <- "umap-learn"

#don't include these things in the title or the legend 
exclude_title_legend <- c("patient_id","sample","tissue","sample_tissue","biopsy_tissue")
colors_file = "C:/Users/JAM526/pancan_fpkmuq/met_analysis/data/color_palettes/pa200_set1_chroma23_lightness23.RDS"

base_dir <-  "C:/Users/JAM526/pancan_fpkmuq/met_analysis"
data_dir <- "data/merged"
data_dir_1 <- "data"
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
  "wnt.csv")

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
  "Wnt")

names(pathway_files) <- pathway_names
used_pathway_files <-  'tp53.csv'
pathway_files <- pathway_files[pathway_files %in% used_pathway_files]

#WHAT ORDER DO THINGS APPEAER IN CANCER
descendant <-     list( "Solid Tissue Normal" = "Primary Tumor",  "Primary Tumor" =   "Metastatic",   "Metastatic" = "None" )

#GIVING BETTER NAMES TO PLOT THINGS BY
plot_names <- list("sample" = "Sample" ,      "disease"  = "Disease" ,    "tissue" = "Tissue" ,"biopsy_tissue" = "Biopsy Tissue", "patient_id" = "Patient ID",    "sample_type" = "Sample Type",   "dataset" = "Dataset" )      

#make variable to sort rows of DF for plotting
source("sample_order.r")
source("all_dup.r")
#exponentiates and decrements
source("exp_dec.r")
#collapses contingency table to be yes/no for one category 
source("contingency_collapse.r")
#takes a patient, data frame with sample information, and a plot; outputs plot with added line between paired samples
source("plot_arrows.r")
#reading the expression data
if(!exists("quasi_tpm")){
quasi_tpm <- readRDS(file.path(base_dir,data_dir,"joined_fpkm_mat.rds"))
}

if(tpm_convert){quasi_tpm <- sweep(quasi_tpm,1,rowSums(quasi_tpm),'/')*10^6}
#loading the mapping between different gene names
load(file.path(base_dir,data_dir_1,"genes_entrez_ensembl.rda"))

#loading the metadata about each sample
all_metadata <-  readRDS(file.path(base_dir,data_dir,"joined_meta.rds")) #recall that a single primary can have multiple mets
hnsc_hpv <- readr::read_tsv("C:/Users/JAM526/pancan_fpkmuq/met_analysis/data/hnsc_hpv.txt")
all_metadata <-  left_join(all_metadata,hnsc_hpv,by = "sample")
##########################################looping through the pathways
for(pathway in pathway_files){
  
  if(show_one_iteration){ graphics.off()}
  
  metadata <- dplyr::filter(all_metadata,!is.na("disease"))
  
  
  if(filter_dataset){ metadata <- dplyr::filter(metadata,dataset %in% dataset_filt)  }
  if(filter_sample_type){  metadata <- dplyr::filter(metadata,sample_type %in% sample_type_filt)}
  if(filter_disease){  metadata <- dplyr::filter(metadata,disease %in% disease_filt)}
 
  if(remove_missing){
   metadata <- dplyr::filter(metadata,!(hpv_status_by_p16_testing %in% "Missing"))
  }
  
  #doing a gene name conversion
  genes_entrez_ensembl_filtered_ordered = genes_entrez_ensembl[match(colnames(quasi_tpm),genes_entrez_ensembl$ensembl_gene_id),]
  
  #reading the actual pathway
  tx <- read.csv(file.path(base_dir,data_dir_1,pathway_lists_dir,pathway),header = F)$V1
  transcript_selector <- match(tx,genes_entrez_ensembl_filtered_ordered$hgnc_symbol)
  #Selecting based on the entrez name rather than logical vector to be safe. deleting NA's from this. 
  transcript_selector <- transcript_selector[!is.na(transcript_selector)] 
  entrez_to_ensembl   <- mapIds(org.Hs.eg.db,as.character(tx),"ENSEMBL","SYMBOL") #actual gene names
  entrez_to_gene_name <- names(entrez_to_ensembl)                  #genes_entrez_ensembl_filtered_ordered[transcript_selector,   "hgnc_symbol"]  #these are the actual gene names
  #getting the expression of the genes in the pathway
  wanted_tpm <- quasi_tpm[,entrez_to_ensembl]
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
  
  #running umap
  embedding <- umap(umap_matrix,config = custom.config, method = umap_method,verbose = TRUE)
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
  is_one <- apply(joined_umap,2,function(x) {length(unique(x))} ) == 1
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
  if(add_white){
  contrasting_colors[1] <- "#cccfcd"
  }
  color_column_string <- quo_name(color_selection)
  color_column <- joined_umap[[color_column_string]] %>% factor() %>% levels()
  wanted_colors <- contrasting_colors[1:length(color_column)]
  names(wanted_colors) <- color_column
col_scale <- scale_colour_manual(values = wanted_colors ,aesthetics = c("colour","fill")) 
fill_scale <- scale_fill_manual(values = wanted_colors ) 
  

###   SORTING THE PLOT DATAFRAME SO THAT THE METASTASES ARE ON TOP 
joined_umap$to_sort <- lapply(FUN = sample_order, X = joined_umap$sample_type) %>% unlist()
joined_umap %<>% arrange(to_sort)


### Makign first plot
out_plot <- ggplot(joined_umap,aes(x=X1,y=X2))+
    geom_point(size = symbol_size,aes(color= !! color_selection,shape= !! shape_selection)) + 
    labs(title = title_string, y="UMAP 2", x="UMAP 1",caption = caption_string,color = plot_names[[quo_name(color_selection)]],shape = plot_names[[quo_name(shape_selection)]])+  col_scale
    
  
if(!axis_numbers){
  out_plot <- out_plot + theme(axis.text = element_blank(),axis.ticks = element_blank()) 
}
  print(out_plot)
  
  if(save_image){ #saving image
    fp_string <- file.path(base_dir,output_dir,paste(title_string,sep = ""))
    fp_string <- gsub(" ","_", fp_string,fixed = TRUE)                  
    fp_string <- paste(fp_string,".jpg",sep = "")   
    ggsave(fp_string,device = "jpeg")
  }

  
  #THIS DRAWS THE LINES BY PUTTING THINGS ON THE OUTPUT PLOT 
  if(draw_lines){
    
    
    for (d in 1:dim(joined_umap)[1]){ #looping through output and adding arrows
      st_in <- joined_umap$sample[d]
      out_plot <- plot_arrows(st_in,joined_umap,out_plot)
    }

    print( out_plot)  #priting the version with lines
    
    if(save_image){
      fp_string <- file.path(base_dir,output_dir,paste(title_string,"Lines",sep = ""))
      fp_string <- gsub(" ","_", fp_string,fixed = TRUE)                  
      fp_string <- paste(fp_string,".jpg",sep = "")   
      ggsave(fp_string,device = "jpeg")
    }
  } else{}
  
  
  
  
  
  
  # PLOTTING THE UMAP OUTPUT COLORED BY CLUSTER
  title_string2 <- paste0("HBDSCAN of " ,pway_name_string)
  
                                                  # shape = plot_names[[quo_name(shape_selection)]]
  out_plot <- ggplot(joined_umap,aes(x=X1,y=X2)) + 
    geom_point(size = symbol_size,aes(color=Cluster,shape = !!shape_selection)) + 
    labs(title = title_string2, y="UMAP 2", x="UMAP 1",caption = caption_string,shape = plot_names[[quo_name(shape_selection)]])
    
    if(!axis_numbers){
    out_plot <- out_plot + theme(axis.text = element_blank(),axis.ticks = element_blank()) 
    }
  print(out_plot)
  if(save_image){
    fp_string <- file.path(base_dir,output_dir,paste(title_string,"_HDBSCAN",sep = ""))
    fp_string <- gsub(" ","_", fp_string,fixed = TRUE)
    fp_string <- paste(fp_string,".jpg",sep = "")
    ggsave(fp_string,device = "jpeg")
  }
  ######STARTING THE PART WHERE WE RUN STATS
  if(run_stats){
  #the subset of tumors
  
  if(subset_for_stats){    
  subset_umap <-  filter(joined_umap,!! filter_exp)
  } else{ subset_umap <- joined_umap}
  #################
  #a contingency table of cluster vs. met site
    
   tbl <- table(as.character(subset_umap$Cluster),unlist(lapply(subset_umap[[chi_sq_var]],as.character)))
  #converting the above table to matrix
  tbl_mat <- as.matrix(tbl)
  
  
 
  #running chi square on the table
  cs <- chisq.test(tbl)
  
  
  #running post-hoc fisher test on the frequency table 
  post_hoc_res <- cs$stdres
  rownames(post_hoc_res) <- rownames(tbl_mat)
  colnames(post_hoc_res) <- colnames(tbl_mat)
  post_hoc <- pchisq((post_hoc_res)^2,df = 1,lower.tail = FALSE)



  
  
#this next line is redundant
  freq_tab <- dplyr::count(subset_umap,(!!sym(chi_sq_var)),Cluster)

    #APPENDING THE POST-HOC P VALUE TO DATA FRAME
  #############################!!!!!!!!!!!!!!!!!!!!!!!
  freq_tab$residual_p <- post_hoc[cbind( as.character(freq_tab$Cluster),as.character(freq_tab[[chi_sq_var]]))]
  fish_sig <- freq_tab$residual_p < print_p_threshold
  
  freq_tab %<>% mutate(residual_p = paste("P =",formatC(.data$residual_p, format = "e", digits = 2)))
  freq_tab$residual_p[!fish_sig] <- ""
  freq_tab %<>% mutate_if(is.character,tools::toTitleCase)
  freq_tab %<>% filter( Cluster != 0)
  freq_tab %<>% filter( !!sym(chi_sq_var) != "Missing")
  #MAKING THE bar PLOT
  title_string3 <- paste0(" Cluster Composition of ",title_string2)
  cs_format <- formatC(cs$p.value, format = "e", digits = 2)
  caption_string <- paste0("Chi-sq P = ",cs_format,"; P < ",as.character(print_p_threshold) ," displayed above")
 
   bp <- ggplot(freq_tab,aes(x = Cluster,y = n,label = residual_p))+fill_scale+
     labs(title = c(title_string3),x = "Cluster",y = "Frequency",caption = caption_string, fill =  plot_names[[quo_name(color_selection)]] )+ 
     geom_bar(stat = "identity")  +aes(fill =  !!sym(color_column_string))
 
  print(bp)
   if(save_image){
    fp_string <- file.path(base_dir,output_dir,paste(title_string,"_cluster_comp",sep = ""))
    fp_string <- gsub(" ","_", fp_string,fixed = TRUE)
    fp_string <- paste(fp_string,".jpg",sep = "")
    ggsave(fp_string,device = "jpeg")
   }
  
  if(show_p_vals){
    bp <- bp + geom_text(size = 3, position = position_stack(vjust = 0.5))
    print(bp)
   
    
     if(save_image){
      fp_string <- file.path(base_dir,output_dir,paste(title_string,"_cluster_comp_pval",sep = ""))
      fp_string <- gsub(" ","_", fp_string,fixed = TRUE)
      fp_string <- paste(fp_string,".jpg",sep = "")
      ggsave(fp_string,device = "jpeg")
    }
    
     }
  
  
}
  
}