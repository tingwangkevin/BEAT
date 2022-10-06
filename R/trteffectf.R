trteffectf<-function(tau,chat,ScreenSample,Sample,ctrue){
  sortedtau=unique(sort(c(tau,0,1,chat)))
  totalcutn=length(sortedtau)
  orderchat=which(sortedtau==chat)
  pij=matrix(0,nrow=2,ncol=totalcutn-1)
  nij=pij; mij=pij;
  for (i in 1:(totalcutn-1)){
    for (j in 1:2){# 1 for placebo 2 for trt
      mij[j,i]=length(which(ScreenSample[,3]==(j-1) & ScreenSample[,1]>sortedtau[i] & ScreenSample[,1]<=sortedtau[i+1]))
      enrolled=which(Sample[,3]==(j-1) & Sample[,1]>sortedtau[i] & Sample[,1]<=sortedtau[i+1])
      nij[j,i]=length(enrolled)
      pij[j,i]=sum(Sample[enrolled,2])/length(Sample[enrolled,2])
    }
  }
  
  pCpos=sum(pij[1,orderchat:(totalcutn-1)]*mij[1,orderchat:(totalcutn-1)]/sum(mij[1,orderchat:(totalcutn-1)]))
  pTpos=sum(pij[2,orderchat:(totalcutn-1)]*mij[2,orderchat:(totalcutn-1)]/sum(mij[2,orderchat:(totalcutn-1)]))
  pCneg=sum(pij[1,1:(orderchat-1)]*mij[1,1:(orderchat-1)]/sum(mij[1,1:(orderchat-1)]))
  pTneg=sum(pij[2,1:(orderchat-1)]*mij[2,1:(orderchat-1)]/sum(mij[2,1:(orderchat-1)]))
  pCall=sum(pij[1,]*mij[1,]/sum(mij[1,]))
  pTall=sum(pij[2,]*mij[2,]/sum(mij[2,]))
  
  VpCpos=sum(pij[1,orderchat:(totalcutn-1)]*(1-pij[1,orderchat:(totalcutn-1)])/nij[1,orderchat:(totalcutn-1)]*
            (mij[1,orderchat:(totalcutn-1)]/sum(mij[1,orderchat:(totalcutn-1)]))^2)
  VpTpos=sum(pij[2,orderchat:(totalcutn-1)]*(1-pij[2,orderchat:(totalcutn-1)])/nij[2,orderchat:(totalcutn-1)]*
            (mij[2,orderchat:(totalcutn-1)]/sum(mij[2,orderchat:(totalcutn-1)]))^2)
  VpCneg=sum(pij[1,1:(orderchat-1)]*(1-pij[1,1:(orderchat-1)])/nij[1,1:(orderchat-1)]*
            (mij[1,1:(orderchat-1)]/sum(mij[1,1:(orderchat-1)]))^2)
  VpTneg=sum(pij[2,1:(orderchat-1)]*(1-pij[2,1:(orderchat-1)])/nij[2,1:(orderchat-1)]*
            (mij[2,1:(orderchat-1)]/sum(mij[2,1:(orderchat-1)]))^2)
  VpCall=sum(pij[1,]*(1-pij[1,])/nij[1,]*(mij[1,]/sum(mij[1,]))^2)
  VpTall=sum(pij[2,]*(1-pij[2,])/nij[2,]*(mij[2,]/sum(mij[2,]))^2)
  
  rejpos=ifelse((pTpos-pCpos)/sqrt(VpCpos+VpTpos)>qnorm(0.95),1,0)
  rejneg=ifelse((pTneg-pCneg)/sqrt(VpCneg+VpTneg)<qnorm(0.05),1,0)
  rejall=ifelse((pTall-pCall)/sqrt(VpCall+VpTall)>qnorm(0.95),1,0)
  
  truetrteffects=truetrteffect(chat,ctrue)
  poschatin=ifelse(truetrteffects[1]>=(pTpos-pCpos)-1.96*sqrt(VpCpos+VpTpos) & 
                   truetrteffects[1]<=(pTpos-pCpos)+1.96*sqrt(VpCpos+VpTpos), 1, 0)
  posctruein=ifelse(truetrteffects[2]>=(pTpos-pCpos)-1.96*sqrt(VpCpos+VpTpos) & 
                    truetrteffects[2]<=(pTpos-pCpos)+1.96*sqrt(VpCpos+VpTpos), 1, 0)
  negchatin=ifelse(truetrteffects[3]>=(pTneg-pCneg)-1.96*sqrt(VpCneg+VpTneg) & 
                   truetrteffects[3]<=(pTneg-pCneg)+1.96*sqrt(VpCneg+VpTneg), 1, 0)
  negctruein=ifelse(truetrteffects[4]>=(pTneg-pCneg)-1.96*sqrt(VpCneg+VpTneg) & 
                    truetrteffects[4]<=(pTneg-pCneg)+1.96*sqrt(VpCneg+VpTneg), 1, 0)
  overallin=ifelse(truetrteffects[5]>=(pTall-pCall)-1.96*sqrt(VpCall+VpTall) & 
                   truetrteffects[5]<=(pTall-pCall)+1.96*sqrt(VpCall+VpTall), 1, 0)
  
  return(c(rejpos,(pTpos-pCpos),sqrt(VpCpos+VpTpos),poschatin, posctruein,
           rejneg,(pTneg-pCneg),sqrt(VpCneg+VpTneg),negchatin, negctruein,
           rejall,(pTall-pCall),sqrt(VpCall+VpTall),overallin))
}