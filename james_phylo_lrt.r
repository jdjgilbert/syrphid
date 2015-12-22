
### write a function to perform a LRT on a pgls

lrt <- function(m1,m2) {
  LR=2*(logLik(m1)-logLik(m2))
  DF=(attr(logLik(m1),'df')-attr(logLik(m2),'df'))
  P=1-pchisq(LR, DF)
  list(LR=as.vector(LR), DF=DF, P=as.vector(P))
}

