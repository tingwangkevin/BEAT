PPoSBEAT<-function(m1,m2,r1,r2,n1,n2,eta,delta){#Calculate P(P(p1post-p2post>eta|p1,p2,n1,n2)>=delta)
  ptmp=seq(0,1,by=0.01)
  PPoS1=0
  for (rf1 in 0:n1){
    for (rf2 in 0:n2){
      #print(c(rf1,rf2))
      PPoS1=PPoS1+mean(dbinom(rf1,n1,ptmp)*dbeta(ptmp,1+r1,1+m1-r1))*
                  mean(dbinom(rf2,n2,ptmp)*dbeta(ptmp,1+r2,1+m2-r2))*
                  mean(rbeta(20000, 1+r1+rf1, 1+m1+n1-r1-rf1)>rbeta(20000, 1+r2+rf2, 1+m2+n2-r2-rf2))
    }
  }
  return(PPoS1)
}