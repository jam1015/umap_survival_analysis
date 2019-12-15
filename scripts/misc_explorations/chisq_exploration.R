mat <-matrix(c(27,19,26,9,17,31,41,40,4,8,13,28),nrow = 3,byrow = TRUE)
rownames(mat) <- c("science","social","art")
colnames(mat) <- c("low","mid","upmid","upper")
cs <- chisq.test(mat)
cs_res <- cs$stdres



asr_test_1 <- function(mat){
  siz <- size(mat)
  nrows <- siz[1]
  ncols <- siz[2]
  lenth <- nrows*ncols
  res_mat <- matrix(rep(NA,prod(siz)),nrow = nrows)
  N <- sum(mat)

  
  
  for (r in 1:nrows){
    for (c in 1:ncols){
      
      O <- mat[r,c]
      nA <- sum(mat[r,])
      nB <- sum(mat[,c])
      E <- (nA*nB)/N
      factor_3 <- (1-(nA/N))
      factor_4 <- (1-(nB/N))
      sqrt_arg <- (factor_3*factor_4*E)
      asr <- (O-E)/sqrt(sqrt_arg)
      res_mat[r,c] <- asr 
    
    }
  }
  return(res_mat)
}

asr_test_2 <- function(mat){
  siz <- size(mat)
  nrows <- siz[1]
  ncols <- siz[2]
  lenth <- nrows*ncols
  res_mat <- matrix(rep(NA,prod(siz)),nrow = nrows)
  N <- sum(mat)
  
  
  
  for (r in 1:nrows){
    for (c in 1:ncols){
      
      O <- mat[r,c]
      nA <- sum(mat[r,])
      nB <- sum(mat[,c])
      E <- (nA*nB)/N
      factor_3 <- (1-(nA/N))
      factor_4 <- (1-(nB/N))
      sqrt_arg <- (nA*nB*factor_3*factor_4)/N
      asr <- (O-E)/sqrt(sqrt_arg)
      res_mat[r,c] <- asr 
      
    }
  }
  return(res_mat)
}

custom1 <- asr_test(mat)
custom2 <- asr_test(mat)