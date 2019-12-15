sample_order <- function(in_str){
  if(strcmp(in_str,"Primary Tumor")){
    return(1) } else if(strcmp(in_str,"Metastatic")){
      return(2) }else if(strcmp(in_str,"Solid Tissue Normal")){
        return(3) }else {
          return(NA) }
  
}
