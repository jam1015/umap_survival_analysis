rm(list = ls())
graphics.off()
library("shiny")
library("tidyverse")
library("reshape2")
library("org.Hs.eg.db")
library("CORElearn")
library("e1071")
library("pbapply")
library("EnhancedVolcano")

#loading data
joined_fpkm                                     <- readRDS("./data/joined_fpkm.rds")
joined_meta                                    <- readRDS("./data/joined_meta.rds")

#filtering
joined_meta <- joined_meta[joined_meta$dataset=="met500",]
joined_meta <- joined_meta[joined_meta$disease=="BRCA",]
joined_fpkm <- joined_fpkm[joined_fpkm$sample %in% joined_meta$sample,]
has_na <- unlist(lapply(lapply(joined_fpkm,is.na),any))
joined_fpkm <- joined_fpkm[,!has_na]

#selecting numeric
joined_fpkm_numeric_indx                        <- unlist(lapply(joined_fpkm,is.numeric))
joined_fpkm_numeric                             <- joined_fpkm[,joined_fpkm_numeric_indx]


inc  <-     function(x)  {x <- x+1}
meanNorm <- function(x) {x/mean(x)}
sumNorm <- function(x){x/sum(x)}

txform <- "inc_meanNorm_log2_var"
switch(txform,
                var = {         txform_fun <-     function(x){var(x)}},
               inc_var  ={      txform_fun <-     function(x){x <- var(inc(x)) }},   
       inc_meanNorm_var = {     txform_fun <- function(x){x <- var(meanNorm(inc(x)))}},
       inc_meanNorm_log2_var = {txform_fun <- function(x){x <- var( log2(meanNorm(inc(x))))}}
       
       )

joined_fpkm_numeric_txform               <- pbapply(joined_fpkm_numeric      ,FUN = txform_fun,  MARGIN = 2)
joined_fpkm_numeric_txform_rank          <-    rank(joined_fpkm_numeric_txform,ties.method = "first") 

threshold <- 5000
pass_ind                                   <-  which(joined_fpkm_numeric_txform_rank >=  (length(joined_fpkm_numeric_txform_rank)  -threshold))

#############  #Using the pass indices to get tissue-wise means for each gene
tpm_pass <- bind_cols(joined_fpkm[,!joined_fpkm_numeric_indx],joined_fpkm[,pass_ind])  #adding on sample variable
joined_meta_tpm_pass                                <- full_join(joined_meta,tpm_pass,by = "sample")  
joined_meta_tpm_pass_numeric_ind                    <- unlist(lapply(joined_meta_tpm_pass,is.numeric))#the metadata have no numeric variables
joined_meta_tpm_pass_numeric <- joined_meta_tpm_pass[,joined_meta_tpm_pass_numeric_ind]
joined_meta_tpm_pass_numeric_names              <- colnames(joined_meta_tpm_pass_numeric) #gettign the response variables, no p values yet
joined_meta_tpm_pass_numeric_inc <- pbapply(joined_meta_tpm_pass_numeric      ,FUN = inc,  MARGIN = 2)
pass_inc_vals_tissue_mean <- as_tibble(aggregate(joined_meta_tpm_pass_numeric_inc, 
                        list(by = joined_meta_tpm_pass$biopsy_tissue), 
                        FUN = mean ))               #these are the actual tissue means

pass_inc_vals_tissue_mean_num_ind <-  unlist(lapply(pass_inc_vals_tissue_mean,is.numeric))  #getting the indices numeric vals of the tissuewise means (getting rid of 'by' variable)
pass_inc_vals_tissue_mean_num     <-  pass_inc_vals_tissue_mean[,pass_inc_vals_tissue_mean_num_ind] 
############## getting mean across ALL tissues, dividing by mean, normalizing by log, 
joined_meta_tpm_pass_numeric_inc_mean <- colMeans(joined_meta_tpm_pass_numeric_inc) #the mean across all tissues
pass_inc_vals_tissue_mean_num_norm_log2 <- log2(sweep(x=pass_inc_vals_tissue_mean_num,MARGIN=2,STATS=joined_meta_tpm_pass_numeric_inc_mean,FUN = "/")) #normalizing each by the total mean
#adding back on non num (the by column)
pass_inc_vals_tissue_mean_norm_log2 <- bind_cols(pass_inc_vals_tissue_mean[!pass_inc_vals_tissue_mean_num_ind],pass_inc_vals_tissue_mean_num_norm_log2)



#running anova
used_df_num <- log(joined_meta_tpm_pass_numeric[,joined_meta_tpm_pass_numeric_names  ]+1)

g <- factor(joined_meta_tpm_pass$biopsy_tissue)
ano_fun <- function(x) {pairwise.wilcox.test(x ,g)}
anova_base <- pbapply(X =used_df_num ,FUN = ano_fun,MARGIN = 2)
get_p <- function(x){x$p.value}
p_vals_base <- unlist(lapply(anova_base,get_p))


#deciding which will be the cutoff
p_vals_base_rank <- rank(p_vals_base,ties.method = "first")
p_threshold <- round(threshold*0.1)
#########################################
p_vals_wanted <- p_vals_base_rank[p_vals_base_rank <= p_threshold]
genes_p_filtered <- names(p_vals_wanted)
genes_p_filtered <- unlist(lapply(strsplit(genes_p_filtered,".",fixed = TRUE),FUN = function(x) x[1]))
#getting the mean values from original tpm
p_tibble <- bind_cols(as_tibble("p"),as_tibble(t(p_vals_base))) #making a row of p-values
colnames(p_tibble) <- colnames(pass_inc_vals_tissue_mean_norm_log2)   
pass_inc_vals_tissue_mean_p <- bind_rows(pass_inc_vals_tissue_mean_norm_log2,p_tibble) #putting the p values with the log2 fold changes
pass_inc_vals_tissue_mean_p_numeric <- pass_inc_vals_tissue_mean_p[,unlist(lapply(pass_inc_vals_tissue_mean_p,is.numeric))] #getting the numeric of the above

numeric_mat <- t(data.matrix(pass_inc_vals_tissue_mean_p_numeric,rownames.force = TRUE)) #making anumeric matrix from above and transposing the numeric matrix
colnames(numeric_mat) <- pass_inc_vals_tissue_mean_p$by  #making the "by" the new column names, have to do this because we transpose above
numeric_mat <- numeric_mat[!is.na(numeric_mat[,"p"]),] #getting rid of nan p values 
numeric_mat_df <- data.frame(numeric_mat)     #converting back to dataframe

numeric_mat_df$gene <- rownames(numeric_mat_df) #adding a gene column
all_cols <- colnames(numeric_mat_df)
gather_key <- all_cols[all_cols != "p"  &all_cols != "gene" ]

numeric_mat_df_gathered <- gather(numeric_mat_df,key = "biopsy_tissue",value = "expression",gather_key)   #converting data to tall form 
numeric_mat_df_gathered_united <- unite(numeric_mat_df_gathered,col = "gene_tissue",c("gene","biopsy_tissue"),sep = "_") #binding gene and 
EnhancedVolcano(numeric_mat_df_gathered_united,
                lab = rownames(numeric_mat_df_gathered_united),
                x = 'expression',
                y = 'p',
               # xlim = c(-6, 6),
                title = 'N061011 versus N61311',
                pCutoff = 10e-12,
                FCcutoff = 1.5,
                transcriptPointSize = 1.5,
                transcriptLabSize = 3.0,
                colAlpha = 1,
                cutoffLineType = 'blank',
                cutoffLineCol = 'black',
                cutoffLineWidth = 0.8,
               # hline = c(10e-12, 10e-36, 10e-60,10e-75 ,10e-84),
                #hlineCol = c('grey0', 'grey25','grey50','grey75'),
                hlineType = 'longdash',
                hlineWidth = 0.8,
                gridlines.major = FALSE,
                gridlines.minor = FALSE)  #making stem and leaf plot
  
##############################################


###### running UMAP
#rownames(joined_meta_tpm_pass)<- joined_meta_tpm_pass$sample

#getting rid of has tcga

#is_primary <- grepl('Primary',meta_tpm_joined[,"sample_type.samples"],fixed=TRUE)
#meta_tpm_joined <- meta_tpm_joined[!is_primary,]

expression_indices <- unlist(lapply(joined_meta_tpm_pass, is.numeric))
wanted_meta <- joined_meta_tpm_pass[,!expression_indices]
umap_matrix <- data.matrix(joined_meta_tpm_pass[,expression_indices])
rownames(umap_matrix)<- wanted_meta$sample
umap_matrix <- umap_matrix[,genes_p_filtered]


heatplot_matrix <- umap_matrix

#running umap
custom.config <- umap.defaults
custom.config$n_epochs <- 1000
custom.config$min_dist <-0.001

embedding <- umap(umap_matrix,config = custom.config, method ="naive")
layout_df <- data.frame(embedding$layout)
#M220S
#trimming down the metadata to what we actually ran
metadata <- joined_meta_tpm_pass[match(rownames(layout_df),wanted_meta$sample),]

joined_umap <- merge(layout_df,metadata,by.x = 0,by.y = "sample")

Legend <- paste(joined_umap[,"sample_type"],joined_umap[,"dataset"],joined_umap[,"biopsy_tissue"])
Legend <- trimws(Legend)
Legend <- factor(Legend)

joined_umap$Cluster <- factor(hdbscan(cbind(joined_umap$X1,joined_umap$X2),minPts = 5)$cluster)
joined_umap$Legend  <- Legend


#MAKE HEATMAP
#heatplot(umap_matrix,"scale" = "row",classvec = joined_umap$Site)# classvec2 = joined_umap$sample_type.samples )


# unlist(lapply(joined_umap$Legend, as.character)) %in% "Primary Tumor Lee"
# 
# joined_umap_primary_lee <- joined_umap[,]

out_plot <- ggplot(data.frame(joined_umap),aes(x=X1,y=X2))+geom_point(size =3,aes(color=Legend,shape=sample_type)) + 
  labs(title="BRCA UMAP", y="UMAP 2", x="UMAP 1")+theme_classic()+ theme(axis.text = element_blank(),axis.ticks = element_blank())+ scale_colour_brewer(palette = "Paired") 

# for (d in 1:dim(joined_umap)[1]){
#   st_in <- (rownames(joined_umap))[d]
#   out_plot <- plot_arrows(st_in,joined_umap,out_plot)
# }
out_plot

 
