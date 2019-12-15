#includes some normalization that I do in advance in othe rscript

rm(list = ls())
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


 fn <- c(
 "cell_cycle.csv",
 "chol.csv",
 "fabo.csv",
 "glycolysis.csv",
 "hippo.csv",
 "mrp.csv",
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
 "wnt.csv"
 )



for (file_name in fn) {
  graphics.off()
  print(file_name)
#file_name <- ".csv"
set.seed(33)
base_dir <-  "C:/Users/JAM526/pancan_fpkmuq/met_analysis"
data_dir_name <- "met500"
data_dir <- paste("data/",data_dir_name, sep = "")
pathway_lists_dir <- "pathway_lists"


allDup <- function (value){
  duplicated(value) | duplicated(value, fromLast = TRUE)
}

exp_dec <- function(x) {
  (2^x)-1
}

contingency_collapse <- function(table,r,c){
  out_mat <- matrix(NA,2,2)
  out_mat[1,1]  <- table[r,c]
  out_mat[1,2]  <- sum(table[r,])-table[r,c]
  out_mat[2,1]  <- sum(table[,c])-table[r,c]
  out_mat[2,2]  <- sum(table) -(out_mat[1,1] + out_mat[1,2] +out_mat[2,1])
  return(out_mat)
}

plot_arrows <- function(st_in,mat,plot_in){
  plot_out <- plot_in
  ind <- match(rownames(mat),st_in,nomatch=0) > 0
  
  start_point <- unlist(mat[ind,c("X1","X2")])
  
  
  sample_type <- unlist((mat$sample_type.samples)[ind]);
  descendant <-     list( "Solid Tissue Normal" = "Primary Tumor",  "Primary Tumor" =   "Metastatic",   "Metastatic" = "None" )
  selected_descendant <- unlist(descendant[sample_type]);
  
  if (!strcmp(selected_descendant,"None")){
    patient_code <- mat$patient_id[ind]
    
    sub_mat_selector <- which(mat$patient_id %in% patient_code)
    
    # print(length(sub_mat_selector))
    
    if (length(sub_mat_selector)>1){
      sub_mat <- mat[sub_mat_selector,]
      # print(dim(sub_mat))
      sub_mat <- sub_mat[which(sub_mat$sample_type.samples %in% selected_descendant),]
      #  print(dim(sub_mat))
      
      
      
      rep_num <- dim(sub_mat)[1]
      
      x_end <- (sub_mat$X1)
      y_end <- (sub_mat$X2)
      
      x_start = repmat(start_point[1],rep_num,1)
      y_start = repmat(start_point[2],rep_num,1)
      
      plot_mat <- data.frame(x1 = x_start,x2= x_end,y1=y_start,y2 = y_end)
      plot_out <- plot_out + geom_segment(aes(x = x_start,y = y_start,
                                              xend = x_end ,yend =y_end),data = plot_mat,show.legend = FALSE)
    }
  }
  return(plot_out)
}



#some of these lines are not necessary because of the new command 'mapIds'
#load(file.path(base_dir,data_dir,"quasi_tpm.rda"))



quasi_tpm <- readRDS(file.path(base_dir,data_dir,"fpkm.rds"))
#quasi_tpm <- sweep(quasi_tpm,1,rowSums(quasi_tpm),'/')*10^6
quasi_tpm <- sweep(quasi_tpm,2,colSums(quasi_tpm),'/')*10^6
  
  
load(file.path(base_dir,data_dir,"genes_entrez_ensembl.rda"))
genes_entrez_ensembl_filtered_ordered = genes_entrez_ensembl[match(rownames(quasi_tpm),genes_entrez_ensembl$ensembl_gene_id),]

#####################################################
# EDIT THIS LINE TO MAKE GENE NAME CHARACTER VECTOR #
#####################################################

tx <- read.csv(file.path(base_dir,data_dir,pathway_lists_dir,file_name),header = F)$V1
transcript_selector <- match(tx,genes_entrez_ensembl_filtered_ordered$hgnc_symbol)

transcript_selector <- transcript_selector[!is.na(transcript_selector)] 
#Selecting based on the entrez name rather than logical vector to be safe

entrez_to_ensembl   <- mapIds(org.Hs.eg.db,as.character(tx),"ENSEMBL","SYMBOL") #actual gene names
entrez_to_gene_name <- names(entrez_to_ensembl)                  #genes_entrez_ensembl_filtered_ordered[transcript_selector,   "hgnc_symbol"]  #these are the actual gene names

wanted_tpm <- quasi_tpm[entrez_to_ensembl,]
rownames(wanted_tpm) <- entrez_to_gene_name


wanted_tpm_normalized <- sweep(wanted_tpm,2,colSums(wanted_tpm),"/")
wanted_tpm_normalized <- t(wanted_tpm_normalized)

metadata <- read.table(file.path(base_dir,data_dir,"case_ids.csv"),row.names = 1, header = TRUE, sep = ",",as.is = TRUE) #recall that a single primary can have multiple mets

metadata <- metadata[rownames(wanted_tpm_normalized),] #matching up rownames of metadata and expression
metadata_end <- dim(metadata)[2]
tpm_df <- data.frame(wanted_tpm_normalized)

meta_tpm_joined <- transform(merge(metadata,tpm_df,by=0),row.names = Row.names, Row.names = NULL)

has_tcga <-   grepl('TCGA',row.names(meta_tpm_joined),fixed=TRUE)
#getting rid of has tcga
meta_tpm_joined <- meta_tpm_joined[!has_tcga,]
#is_primary <- grepl('Primary',meta_tpm_joined[,"sample_type.samples"],fixed=TRUE)
#meta_tpm_joined <- meta_tpm_joined[!is_primary,]

expression_indices <- (metadata_end+1):NCOL(meta_tpm_joined)
umap_matrix <- data.matrix(meta_tpm_joined[,expression_indices])
umap_matrix_sq <- umap_matrix^2
umap_matrix_sq_sum <- rowSums(umap_matrix^2)
umap_matrix_sq_sum_rep <- replicate(dim(umap_matrix)[2],umap_matrix_sq_sum)
umap_matrix <- sqrt(umap_matrix/umap_matrix_sq_sum_rep)
#umap_matrix_norm <-  

heatplot_matrix <- umap_matrix

#running umap
custom.config <- umap.defaults
custom.config$n_epochs <- 1000
custom.config$min_dist <-0.001

######################################writing numeric_tensorboard file
write.table(umap_matrix, file=paste("tb",data_dir_name,file_name,sep = "_"), quote=FALSE, sep='\t', col.names = FALSE)
embedding <- umap(umap_matrix,config = custom.config, method ="naive")
layout_df <- data.frame(embedding$layout)
#M220S
#trimming down the metadata to what we actually ran
metadata <- meta_tpm_joined[rownames(layout_df),]

joined_umap <-transform( merge(layout_df,metadata,by= 0),row.names = Row.names, Row.names = NULL)

Legend <- paste(joined_umap[,"sample_type.samples"],joined_umap[,"Provider"],joined_umap[,"Site"])


##########################writing meta tensorboard file
write.table(cbind(Legend,merge(layout_df,metadata,by= 0)),file = paste("tb",data_dir_name,file_name,"meta",sep = "_"),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)


# 
# Legend <- trimws(Legend)
# Legend <- factor(Legend)
# 
# joined_umap$Cluster <- factor(hdbscan(cbind(joined_umap$X1,joined_umap$X2),minPts = 5)$cluster)
# joined_umap$Legend <- Legend
# 
# 
# 
# 
# #MAKE HEATMAP
# # heatplot(umap_matrix,"scale" = "row",classvec = joined_umap$Site)# classvec2 = joined_umap$sample_type.samples )
# # unlist(lapply(joined_umap$Legend, as.character)) %in% "Primary Tumor Lee"
# # joined_umap_primary_lee <- joined_umap[,]
# 
# out_plot <- ggplot(data.frame(joined_umap),aes(x=X1,y=X2))+geom_point(size =3,aes(color=Legend,shape=sample_type.samples)) + 
#   labs(title="BRCA UMAP", y="UMAP 2", x="UMAP 1")+theme_classic()+ theme(axis.text = element_blank(),axis.ticks = element_blank())+ scale_colour_brewer(palette = "Set1") 
# 
# for (d in 1:dim(joined_umap)[1]){
#   st_in <- (rownames(joined_umap))[d]
#   out_plot <- plot_arrows(st_in,joined_umap,out_plot)
# }
# out_plot
# 
# out_plot_2 <- ggplot(data.frame(joined_umap),aes(x=X1,y=X2))+geom_point(size =3,aes(color=Cluster,shape=sample_type.samples)) + 
#   labs(title="BRCA UMAP", y="UMAP 2", x="UMAP 1")+theme_classic()+ theme(axis.text = element_blank(),axis.ticks = element_blank()) 
# 
# #plotting the arrows on that 
# for (d in 1:dim(joined_umap)[1]){
#   st_in <- (rownames(joined_umap))[d]
#   out_plot_2 <- plot_arrows(st_in,joined_umap,out_plot_2)
# }
# out_plot_2
# 
# #dividing joined umap into met and primary
# met_umap <- joined_umap[(joined_umap$sample_type.samples) %in% "Metastatic",]
# met_umap_numeric <- as.matrix(met_umap[11:84,])
# 
# prim_umap <- joined_umap[!((joined_umap$sample_type.samples) %in% "Metastatic"),]
# 
# #heatplot(as.matrix(met_umap[7:86]),"scale" = "column",classvec = met_umap$Site)# classvec2 = joined_umap$sample_type.samples )
# 
# #running chi square
# tbl <- table(as.character(met_umap$Cluster),unlist(lapply(met_umap$Legend,as.character)))
# tbl_mat <- as.matrix(tbl)
# cs <- chisq.test(tbl)
# cs
# post_hoc <- matrix(NA,dim(tbl_mat)[1],dim(tbl_mat)[2])
# rownames(post_hoc) <- rownames(tbl_mat)
# colnames(post_hoc) <- colnames(tbl_mat)
# for (e in 1:dim(tbl_mat)[1]){
#   for (f in 1:dim(tbl_mat)[2]){
#     
#     collapsed <- contingency_collapse(tbl_mat,e,f)
#     ft_out <- fisher.test(collapsed)
#     
#     post_hoc[e,f] <- ft_out$p.value
#   }
# }
# 
# #making chi squared bar graph
# bp <- ggplot(met_umap,aes(Cluster,fill=Legend))
# bp <- bp +geom_bar()
# bp +  scale_fill_brewer(palette="Set1")
}
