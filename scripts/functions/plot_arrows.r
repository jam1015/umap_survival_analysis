#takes a patient, data frame with sample information, and a plot; outputs plot with added line between paired samples
plot_arrows <- function(st_in,mat,plot_in){
  plot_out <- plot_in
  ind <- match(mat$sample,st_in,nomatch=0) > 0
  
  start_point <- unlist(mat[ind,c("X1","X2")])
  
  
  sample_type <- unlist((mat$sample_type)[ind]);
  descendant <-     list( "Solid Tissue Normal" = "Primary Tumor",  "Primary Tumor" =   "Metastatic",   "Metastatic" = "None" )
  selected_descendant <- unlist(descendant[sample_type]);
  
  if (!strcmp(selected_descendant,"None")){
    patient_code <- mat$patient_id[ind]
    
    sub_mat_selector <- which(mat$patient_id %in% patient_code)
    
    # print(length(sub_mat_selector))
    
    if (length(sub_mat_selector)>1){
      sub_mat <- mat[sub_mat_selector,]
      sub_mat <- sub_mat[which(sub_mat$sample_type %in% selected_descendant),]
      
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