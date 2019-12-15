#includes some normalization that I do in advance in othe rscrip
# add option to not display labels on output plot
#add option to loop through order layers
#clustering removed. Will add in later. 


rm_list <- ls()
rm_list <- rm_list[!(rm_list %in% c("expression_data","tpm_loaded"))]
rm(list = rm_list)
#reading the expression data


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
library("TCGAbiolinks")
library("survminer")
library("survival")

ggcolor_hue_minus <- function(n) {
  hues = seq(15, 375, length = n + 1)
  full <- hcl(h = hues, l = 65, c = 100)[1:n]
  minus_one <-full[2:numel(full)]
  return(minus_one)
}
#options
add_white <- TRUE
chi_sq_var <- "tissue"
color_selection <- quo(Expression)
shape_selection <- quo(sample_type)
print_p_threshold <- 0.01
use_tpm <- TRUE
standardize <-  FALSE #FALSE#TRUE   #summing to one
spherize <- FALSE
draw_lines <- FALSE  #draw lines between matched samples
filter_dataset <- TRUE   #are we jsut seleecting a few datasets
dataset_filt <- "TCGA"
filter_sample_type <- TRUE
sample_type_filt <- "Primary Tumor"
filter_disease <- TRUE
disease_filts <- c("OV","LUSC","ESCA","HNSC","SARC","LIHC","BLCA","LUAD")
symbol_size <- 2
run_stats <- TRUE
subset_for_stats <- FALSE
ntile <- .25
save_image <- TRUE
axis_numbers <- TRUE
show_p_vals <- TRUE
#should we display every run? 
show_one_iteration  <- FALSE
font_size  <- 14
filter_exp <- quo(sample_type == "Primary Tumor" & hpv_status != "Missing" & Cluster != 0 ) #string for selecting for subtype
file_type <- ".svg"
projects <- getGDCprojects()$project_id
targets_to_find <- c("TCGA","TARGET") 

#getting project variables
log_vecs <- sapply(targets_to_find,str_detect,string = projects)
projects <- projects[rowSums(log_vecs)>0]
log_vecs <-  sapply(disease_filts,str_detect,string = projects)
projects <- projects[rowSums(log_vecs) >0]

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
  "pyriIntermediate Levelsine.csv",
  "rp.csv",
  "tca.csv",
  "tgfb.csv",
  "tp53.csv",
  "wnt.csv",
  "nfe2l2keap1.csv",
  "nfe2l2_targets.csv",
  "nfe2l2.csv",
  "nfe2l2_and_targets.csv")


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
  "PyriIntermediate Levelsine Biosynthesis",
  "Cytoplasmic Ribosomal Proteins",
  "TCA Cycle",
  "TGF-Beta",
  "Tp53",
  "Wnt",
  "NFE2L2 and KEAP1",
  "NFE2L2 Targets",
  "NFE2L2",
  "NFE2L2 and Targets")



names(pathway_files) <- pathway_names
y_pathway_files <- c("nfe2l2_targets.csv")
x_pathway_files <- c("nfe2l2.csv")

y_pathway_files <- pathway_files[pathway_files %in% y_pathway_files]
x_pathway_files <- pathway_files[pathway_files %in% x_pathway_files]

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
all_metadata <-  readRDS(file.path(base_dir,data_dir,"joined_meta.rds")) #recall that a single primary can have multiple mets
surv_data <- readRDS(file.path(base_dir,data_dir,"gdc_survival.rds"))
surv_data <- surv_data[,c("sample","OS","OS.time")]
all_metadata <- inner_join(all_metadata,surv_data,by = "sample")
##########################################looping through the pathways
for(disease_filt in disease_filts) {
  for(j in 1:length(y_pathway_files)){
    if(show_one_iteration){ graphics.off()}
    
    metadata <- dplyr::filter(all_metadata,!is.na("disease"))
    
    
    if(filter_dataset){ metadata <- dplyr::filter(metadata,dataset %in% dataset_filt)  }
    if(filter_sample_type){  metadata <- dplyr::filter(metadata,sample_type %in% sample_type_filt)}
    if(filter_disease){  metadata <- dplyr::filter(metadata,disease %in% disease_filt)}
    
    
    
    #getting survival data directly from tcga
    
    
    
    #doing a gene name conversion
    genes_entrez_ensembl_filtered_ordered = genes_entrez_ensembl[match(colnames(expression_data),genes_entrez_ensembl$ensembl_gene_id),]
    
    #reading the actual pathway
    tx_y <- read.csv(file.path(base_dir,data_dir_1,pathway_lists_dir,y_pathway_files[j]),header = F)$V1 %>% levels()
    
    tx_x <- read.csv(file.path(base_dir,data_dir_1,pathway_lists_dir,x_pathway_files[j]),header = F)$V1 %>% levels()
    
    tx <- c(tx_y,tx_x)
    ### OLD WAY OF GETTING GENE NAMES: DEPRECATED
    # transcript_selector <- match(tx,genes_entrez_ensembl_filtered_ordered$hgnc_symbol)
    #Selecting based on the entrez name rather than logical vector to be safe. deleting NA's from this. 
    #transcript_selector <- transcript_selector[!is.na(transcript_selector)] 
    
    
    entrez_to_ensembl   <- mapIds(org.Hs.eg.db,as.character(tx),"ENSEMBL","SYMBOL") #actual gene names
    entrez_to_gene_name <- names(entrez_to_ensembl)                  #genes_entrez_ensembl_filtered_ordered[transcript_selector,   "hgnc_symbol"]  #these are the actual gene names
    #getting the expression of the genes in the pathway
    wanted_tpm <- expression_data[,entrez_to_ensembl]
    #labeling with hgnc gene name
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
    has_gen <-  apply(meta_tpm_joined[,tx], 1, function(x) !any(is.na(x)))
    meta_tpm_joined <- meta_tpm_joined[has_gen,]
    #what we actually run umap on 
    scatter_matrix <- data.matrix(meta_tpm_joined[,tx])
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
    
    
    
    #calling it joined scatter. 
    meta_tpm_joined <- meta_tpm_joined
    meta_tpm_joined[,tx] <- scatter_matrix 
    
    
    #making title
    pway_name_string_y <- names(y_pathway_files)[y_pathway_files %in% y_pathway_files[j]]
    pway_name_string_x <- names(x_pathway_files)[x_pathway_files %in% x_pathway_files[j]]
    is_one <- apply(meta_tpm_joined,2,function(x) {length(unique(x))} ) == 1
    title <-(meta_tpm_joined[,is_one])
    title %<>% dplyr::select_if(~!all(is.na(.)))
    title <- dplyr::select(title, - !!((exclude_title_legend[exclude_title_legend %in% colnames(title)])))
    
    
    title_string <- paste("Scatter of ", 
                          paste(pway_name_string_y,"vs.",pway_name_string_x), 
                          apply( title[1,] , 1 , paste , collapse = " " )) #pasting first line of title df with spaces
    
    
    #Did we apply these normalizations to data? if so include in title
    if(spherize) { title_string <- paste(title_string, "Spherized")}  
    if(standardize) { title_string <- paste(title_string, "Standardized")}
    
    
    
    
    
    meta_tpm_joined %<>% mutate_if(is.character,tools::toTitleCase)
    
    
    
    ### Makign first plot
    to_plot <- colnames(meta_tpm_joined[,tx])
    
    vec_x <- rowMeans(meta_tpm_joined[,tx_x,drop = FALSE])
    vec_y <- rowMeans(meta_tpm_joined[,tx_y,drop = FALSE])
    
    
    
    up_ntile_x_val <-   quantile(vec_x,1-ntile)
    down_ntile_x_val <-  quantile(vec_x,ntile)
    up_ntile_y_val <-  quantile(vec_y,1-ntile) 
    down_ntile_y_val <-  quantile(vec_y,ntile)
    
    
    up_ntile_x <- vec_x  >  up_ntile_x_val
    down_ntile_x <- vec_x <  down_ntile_x_val 
    up_ntile_y <- vec_y  > up_ntile_y_val
    down_ntile_y <- vec_y <  down_ntile_y_val
    
    
    meta_tpm_joined$Expression     <- "Intermediate Levels"
    meta_tpm_joined$Expression[up_ntile_y] <- "Elevated Levels"
    meta_tpm_joined$Expression[down_ntile_y] <- "Decreased Levels"
    meta_tpm_joined$vec_x <- vec_x
    meta_tpm_joined$vec_y <- vec_y
    #  meta_tpm_joined <- meta_tpm_joined[meta_tpm_joined$Expression != "Intermediate Levels",]
    
    #### Making the color scales
    contrasting_colors <- c("#0000ff","#ff0000")#tolower(readRDS(file = colors_file )) #reading from the colors file (a list of color hex values)
    if(add_white){
      contrasting_colors <- c("#cccfcd",contrasting_colors)
    }
    color_column_string <- quo_name(color_selection)
    color_column <- meta_tpm_joined[[color_column_string]] %>% factor() %>% levels()
    
    if(any(color_column == "Intermediate Levels")){
      color_column <- color_column[color_column != "Intermediate Levels"]
      color_column <- c("Intermediate Levels",color_column)
    }
    
    
    
    wanted_colors <- contrasting_colors[1:length(color_column)]
    names(wanted_colors) <- color_column
    col_scale <- scale_colour_manual(values = wanted_colors ,aesthetics = c("colour","fill")) 
    fill_scale <- scale_fill_manual(values = wanted_colors ) 
    
    #### making scatter plot
    out_plot <- ggplot(meta_tpm_joined,aes_string(x=vec_x,y = vec_y))+
      geom_point(size = symbol_size,aes(color=!! color_selection,shape= !! shape_selection)) + 
      labs(title = title_string, x= names(x_pathway_files[j]), y= names(y_pathway_files[j]),shape = plot_names[[quo_name(shape_selection)]],caption = paste("N-tile cut is",as.character(ntile)))+  col_scale+
      geom_hline(yintercept = down_ntile_y_val )+geom_hline(yintercept =up_ntile_y_val)+theme(aspect.ratio=1)  #+geom_vline(xintercept=up_ntile_x_val )+geom_vline(xintercept = down_ntile_x_val)+
    
   out_plot <- ggpar(out_plot,font.x = font_size,font.y =font_size,font.xtickslab = font_size,font.ytickslab = font_size,font.main = font_size,font.legend = font_size)
    
    
    print(out_plot)
    
    if(save_image){ #saving image
      fp_string <- file.path(base_dir,output_dir,paste(title_string,sep = ""))
      fp_string <- gsub(" ","_", fp_string,fixed = TRUE)                  
      fp_string <- paste(fp_string,".svg",sep = "")   
      ggsave(fp_string,device = "svg")
    }
    
    
    #making survival plot
    title_string_surv <- paste("Survival of", 
                               pway_name_string_y, 
                               apply( title[1,] , 1 , paste , collapse = " " )) #pasting first line of title df with spaces
    
    names(wanted_colors) <- paste0(quo_name(color_selection),"=",names(wanted_colors))
    
    meta_tpm_joined_surv <- filter(meta_tpm_joined,Expression != "Intermediate Levels")
    fit <- survfit(Surv(OS.time,OS) ~ Expression,data = meta_tpm_joined_surv)
    sp <- ggsurvplot(fit,
                     data = meta_tpm_joined,palette = c("blue","red"), 
                     pval = TRUE,title = title_string_surv,
                     ggtheme = theme_gray(),
                     ylim = c(0,1),
                     xlim = c(0,max(meta_tpm_joined$OS.time,na.rm = TRUE)),font.x = font_size,font.y =font_size,font.xtickslab = font_size,font.ytickslab = font_size,font.main = font_size,font.legend = font_size,censor.shape  = 124,xlab = "Time (Days)") + 
      labs(caption = paste("N-tile cut is",as.character(ntile))) 
    sp$plot <- sp$plot + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))+ theme(aspect.ratio=1)
    print(sp)
    
    
    if(save_image){
      fp_string <- file.path(base_dir,output_dir,paste(title_string,"_surv_comp",sep = ""))
      fp_string <- gsub(" ","_", fp_string,fixed = TRUE)
      fp_string <- paste(fp_string,".svg",sep = "")
      ggsave(fp_string,device = "svg")
    }
    
    
    
    ## making survival plot with risk table 
    title_string_surv <- paste("Survival of", 
                               pway_name_string_y, 
                               apply( title[1,] , 1 , paste , collapse = " " )) #pasting first line of title df with spaces
    
    names(wanted_colors) <- paste0(quo_name(color_selection),"=",names(wanted_colors))
    
    meta_tpm_joined_surv <- filter(meta_tpm_joined,Expression != "Intermediate Levels")
    fit <- survfit(Surv(OS.time,OS) ~ Expression,data = meta_tpm_joined_surv)
    sp <- ggsurvplot(fit,risk.table = TRUE,
                     data = meta_tpm_joined,palette = c("blue","red"), 
                     pval = TRUE,title = title_string_surv,
                     ggtheme = theme_gray(),
                     ylim = c(0,1),
                     xlim = c(0,max(meta_tpm_joined$OS.time,na.rm = TRUE)),font.x = font_size,font.y =font_size,font.xtickslab = font_size,font.ytickslab = font_size,font.main = font_size,font.legend = font_size,censor.shape  = 124,xlab = "Time (Days)") + 
      labs(caption = paste("N-tile cut is",as.character(ntile))) 
    sp$plot <- sp$plot + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
    print(sp)
    
    
    if(save_image){
      fp_string <- file.path(base_dir,output_dir,paste(title_string,"_surv_comp_rt",sep = ""))
      fp_string <- gsub(" ","_", fp_string,fixed = TRUE)
      fp_string <- paste(fp_string,".svg",sep = "")
      ggsave(fp_string,device = "svg")
    }
    
    
    
    
    
    
    ## making pairwise surv table
    
    pw <- pairwise_survdiff(Surv(OS.time,OS) ~ Expression,data = meta_tpm_joined_surv,)
    if(save_image){
      
      
      fp_string <- file.path(base_dir,output_dir,paste(title_string,sep = ""))
      fp_string <- gsub(" ","_", fp_string,fixed = TRUE)                  
      fp_string <- paste(fp_string,".csv",sep = "")   
      write.csv(x = data.frame(pw$p.value),file = fp_string)
      
    }
    
    
    
  }
}
