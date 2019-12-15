#collapses contingency table to be yes/no for one category 
contingency_collapse <- function(table,r,c){
  out_mat <- matrix(NA,2,2)
  out_mat[1,1]  <- table[r,c]
  out_mat[1,2]  <- sum(table[r,])-table[r,c]
  out_mat[2,1]  <- sum(table[,c])-table[r,c]
  out_mat[2,2]  <- sum(table) -(out_mat[1,1] + out_mat[1,2] +out_mat[2,1])
  return(out_mat)
}
