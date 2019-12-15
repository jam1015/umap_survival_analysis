#main script to analyze hpv+/- cesc tumors
rm_list <- ls()
rm_list <- rm_list[!(rm_list %in% "quasi_tpm")]
rm(list = rm_list)
setwd("C:\\Users/JAM526/pancan_fpkmuq/met_analysis/scripts/")
graphics.off()