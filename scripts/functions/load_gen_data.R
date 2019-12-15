#the idea is to load the data we want only if it is not already loaded
#is a function but used like a script
load_gen_data <- function(){
if(!exists("expression_data")){
  
  
  if(use_tpm){ 
    expression_data <- readRDS(file.path(base_dir,data_dir,tpm_file))
    tpm_loaded <<- TRUE
  } else {
    expression_data <- readRDS(file.path(base_dir,data_dir,fpkm_file))
    tpm_loaded <<- FALSE 
  } 
  
  
} else {
  
  if(use_tpm){
    if(tpm_loaded){
    }else{
      expression_data <- readRDS(file.path(base_dir,data_dir,tpm_file))
      tpm_loaded <<- TRUE
    }
  } else {
    if(tpm_loaded){
      expression_data <- readRDS(file.path(base_dir,data_dir,fpkm_file))
      tpm_loaded <<- FALSE
    }else{
    }
    
  }  
}
return(expression_data)
}