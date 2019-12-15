#reads a color palette and extends it by permuting it

library(combinat)
basic_vector <- data.matrix(read.csv("../data/color_palettes/pa200_set1_chroma23_lightness23.txt",header = FALSE))

dims <- dim(basic_vector)

perms <- permn(1:dims[2])
length_perms <- length(perms)
template_vector <- matrix(data=NA,nrow=dims[1]*length_perms,ncol=dims[2])
for (j in 1:length_perms ){
 
  
   cat_mat <- basic_vector[,perms[[j]]]
  start_ind <- ((j-1)*dims[1])+1
  end_ind <-   j*dims[1]
  
  print(start_ind)
  print(end_ind)
  template_vector[start_ind:end_ind,] <- cat_mat
  }

template_vector <- as.character(as.hexmode(round(template_vector *255)))
output_colors <- matrix(template_vector,byrow = TRUE,ncol = dims[2])
output_colors <- apply(X = output_colors, MARGIN = 1, FUN = paste0, collapse = "")
output_colors <- paste0("#",output_colors)
saveRDS(output_colors)

