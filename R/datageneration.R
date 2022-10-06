datageneration<-function(trt,probfunc){
  X=runif(1)
  ptmp=pfunc(trt,X,probfunc)
  Y=rbinom(1,1,ptmp)
  return(c(round(X,4),round(Y)))
}
