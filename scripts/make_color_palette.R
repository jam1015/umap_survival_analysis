#makes color palette from a premade set
library(combinat)
rm(list = ls())
in_string = "../data/color_palettes/new_pastel"
basic_vector <- data.matrix(read.csv(paste0(in_string,".txt"),header = FALSE))

dims <- dim(basic_vector)




basic_vector <- as.character(as.hexmode(basic_vector))
output_colors <- matrix(basic_vector,byrow = FALSE,ncol = dims[2])
output_colors <- apply(X = output_colors, MARGIN = 1, FUN = paste0, collapse = "")
output_colors <- paste0("#",output_colors)
saveRDS(output_colors,file = paste0(in_string,".RDS"))

