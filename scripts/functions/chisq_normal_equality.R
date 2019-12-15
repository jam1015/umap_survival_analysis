x <-1.96
(pnorm(abs(x),lower.tail = FALSE))*2
  pchisq(x^2,df = 1,lower.tail = FALSE)