PPoSF<-function(part,delta,ScreenSample,Sample,tau,tautmp,chat,nk,k){
  library(RBesT)
  options(warn=-1)
  if (part==0){#Calculate PPoF
    sortedtau=unique(sort(c(tau[1:k],0,tautmp)))
    ordertautmp=which(sortedtau==tautmp)
    pij=matrix(0,nrow=2,ncol=ordertautmp-1)
    nij=pij; mij=pij;
    for (i in 1:(ordertautmp-1)){
      for (j in 1:2){# 1 for placebo 2 for trt
        mij[j,i]=length(which(ScreenSample[,3]==(j-1) & ScreenSample[,1]>sortedtau[i] & ScreenSample[,1]<=sortedtau[i+1]))
        enrolled=which(Sample[,3]==(j-1) & Sample[,1]>sortedtau[i] & Sample[,1]<=sortedtau[i+1])
        nij[j,i]=length(enrolled)
        if (length(Sample[enrolled,2])>0){
          pij[j,i]=sum(Sample[enrolled,2])/length(Sample[enrolled,2])
        }else{
          pij[j,i]=0
        }
      }
    }
    pCk=sum(pij[1,]*mij[1,]/sum(mij[1,]))
    pTk=sum(pij[2,]*mij[2,]/sum(mij[2,]))
    oCk=sum(nij[1,]); oTk=sum(nij[2,])
    wCk=round(oCk*pCk); wTk=round(oTk*pTk)
    ToobsC=floor(sum(nk[(k+1):length(nk)])*sum(mij[1,])/length(which(ScreenSample[,3]==0)))
    ToobsT=floor(sum(nk[(k+1):length(nk)])*sum(mij[2,])/length(which(ScreenSample[,3]==1)))
    PPoSFf=pos2S(mixbeta(c(1, 1+wCk,1+oCk-wCk)), mixbeta(c(1, 1+wTk,1+oTk-wTk)), ToobsC, ToobsT,  decision2S(c(delta) , c(0), FALSE))
    if (is.na(PPoSFf(mixbeta(c(1, 1+wCk,1+oCk-wCk)), mixbeta(c(1, 1+wTk,1+oTk-wTk))))){
      #print("run PPOSBEAT")
      PPoSFv=PPoSBEAT(oCk,oTk,wCk,wTk,ToobsC,ToobsT,0,delta)
    }else{
      PPoSFv=PPoSFf(mixbeta(c(1, 1+wCk,1+oCk-wCk)), mixbeta(c(1, 1+wTk,1+oTk-wTk)))
    }
  }else{
    sortedtau=unique(sort(c(tau[1:k],chat,1),decreasing = TRUE))
    orderchat=which(sortedtau==chat)
    pij=matrix(0,nrow=2,ncol=orderchat-1)
    nij=pij; mij=pij;
    for (i in 1:(orderchat-1)){
      for (j in 1:2){# 1 for placebo 2 for trt
        mij[j,i]=length(which(ScreenSample[,3]==(j-1) & ScreenSample[,1]<=sortedtau[i] & ScreenSample[,1]>sortedtau[i+1]))
        enrolled=which(Sample[,3]==(j-1) & Sample[,1]<=sortedtau[i] & Sample[,1]>sortedtau[i+1])
        nij[j,i]=length(enrolled)
        if (length(Sample[enrolled,2])>0){
          pij[j,i]=sum(Sample[enrolled,2])/length(Sample[enrolled,2])
        }else{
          pij[j,i]=0
        }
      }
    }
    pCk=sum(pij[1,]*mij[1,]/sum(mij[1,]))
    pTk=sum(pij[2,]*mij[2,]/sum(mij[2,]))
    oCk=sum(nij[1,]); oTk=sum(nij[2,])
    wCk=round(oCk*pCk); wTk=round(oTk*pTk)
    ToobsC=floor(sum(nk[(k+1):length(nk)])*sum(mij[1,])/length(which(ScreenSample[,3]==0)))
    ToobsT=floor(sum(nk[(k+1):length(nk)])*sum(mij[2,])/length(which(ScreenSample[,3]==1)))
    PPoSFf=pos2S(mixbeta(c(1, 1+wTk,1+oTk-wTk)), mixbeta(c(1, 1+wCk,1+oCk-wCk)), ToobsT, ToobsC, decision2S(c(delta) , c(0), FALSE))
    if (is.na(PPoSFf(mixbeta(c(1, 1+wTk,1+oTk-wTk)), mixbeta(c(1, 1+wCk,1+oCk-wCk))))){
      #print("run PPOSBEAT")
      PPoSFv=PPoSBEAT(oTk,oCk,wTk,wCk,ToobsT,ToobsC,0,delta)
    }else{
      PPoSFv=PPoSFf(mixbeta(c(1, 1+wTk,1+oTk-wTk)), mixbeta(c(1, 1+wCk,1+oCk-wCk)))
    }
  }
  return(PPoSFv)
}