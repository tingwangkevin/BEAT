markerestf<-function(Sample){
  mylogit <- glm(response ~ marker + trt + marker*trt, data = Sample, family = "binomial")
  gamma2=summary(mylogit)$coefficient[3,1]
  gamma3=summary(mylogit)$coefficient[4,1]
  xhat=-gamma2/gamma3
  chat=min(max(0,xhat),1)
  vxhat=t(c(0,0,-1/gamma3,gamma2/gamma3^2))%*%vcov(mylogit)%*%t(t(c(0,0,-1/gamma3,gamma2/gamma3^2)))
  mu=xhat; sigma=sqrt(vxhat)
  Xsample<-rnorm(100000, mean=mu, sd=sigma)
  X01sample<-Xsample[which(Xsample>0 & Xsample<1)]
  alpha1=(0-mu)/sigma; alpha2=(1-mu)/sigma
  if (xhat>0 & xhat<1){
    vchat=xhat^2*pnorm(0,mean=mu,sd=sigma)+ mean((X01sample-xhat)^2) + (1-xhat)^2*pnorm(1,mean=mu,sd=sigma,lower.tail = FALSE)
          #(sigma^2+mu^2-sigma^2*(alpha2*dnorm(alpha2)-alpha1*dnorm(alpha1))/(pnorm(alpha2)-pnorm(alpha1))-
          #2*mu*sigma*(dnorm(alpha2)-dnorm(alpha1))/(pnorm(alpha2)-pnorm(alpha1)))-
          #2*xhat*(mu-sigma*(dnorm(alpha2)-dnorm(alpha1))/(pnorm(alpha2)-pnorm(alpha1)))+
          #xhat^2*(pnorm(1,mean=mu,sd=sigma)-pnorm(0,mean=mu,sd=sigma))+
  }else if (xhat<=0){
    vchat=mean(X01sample^2)+ pnorm(1,mean=mu,sd=sigma,lower.tail = FALSE)
          #(sigma^2+mu^2-sigma*(alpha2*dnorm(alpha2)-alpha1*dnorm(alpha1))/(pnorm(alpha2)-pnorm(alpha1))-
          #2*mu*sigma*(dnorm(alpha2)-dnorm(alpha1))/(pnorm(alpha2)-pnorm(alpha1)))+
  }else{
    vchat=pnorm(0,mean=mu,sd=sigma)+mean((X01sample-1)^2)
          #(sigma^2+mu^2-sigma*(alpha2*dnorm(alpha2)-alpha1*dnorm(alpha1))/(pnorm(alpha2)-pnorm(alpha1))-
          #2*mu*sigma*(dnorm(alpha2)-dnorm(alpha1))/(pnorm(alpha2)-pnorm(alpha1)))-
          #2*(mu-sigma*(dnorm(alpha2)-dnorm(alpha1))/(pnorm(alpha2)-pnorm(alpha1)))+
          #(pnorm(1,mean=mu,sd=sigma)-pnorm(0,mean=mu,sd=sigma))
  }
  return(c(chat,sqrt(vchat)))
}
