#' BEAT design
#'
#' Conduct BEAT design and obtain all important design characteristics
#' @param K Number of blocks
#' @param Kstar After Kstar blocks' enrollment, if estimate of biomarker is still 1, stop for futility.
#' @param t User-specified parameter that is used to flexibly define the enrollment cutoff.
#' @param nk A vector of length K with number of enrollment per arm per block
#' @param truebio True biomarker value
#' @param deltaeta A vector containing four values: deltapos, etapos, deltaneg, etaneg for PPoS and PPoF
#' @param probfunc A vector of length 2 response probablity functions for control and treatment arms
#' @return (1) biomarker estimator and estimated standard deviation;
#' (2) stop indicator, screened number;
#' (3) enrolled and excluded positive number;
#' (4) treatment effect for positive, negative and overall (rejection probability, estimator, estimated standard deviation, 95\% cover probability based on estimated / true biomarker threshold)
#' @examples
#' K=4; Kstar=3; t=0.5; nk=c(111,111,110,110); truebio=0.6; deltaeta=c(0.5,0.4,0.7,0.7);
#' probfunc=c('exp(-2.2)/(1+exp(-2.2))', 'exp(-4+3*x)/(1+exp(-4+3*x))')
#' BEAT(K, Kstar, t, nk, truebio, deltaeta, probfunc)
#'
#' @export
BEAT<-function(K,Kstar,t,nk,truebio,deltaeta,probfunc){

  deltapos=deltaeta[1]
  etapos=deltaeta[2]
  deltaneg=deltaeta[3]
  etaneg=deltaeta[4]

  tau=c(0,rep(NA,K-1))
  ScreenSampleP=c()
  nScreenP=0
  SampleP=c()
  nP=0
  ScreenSampleT=c()
  nScreenT=0
  SampleT=c()
  nT=0
  Enrollpos=0
  Excludpos=0
  k=1
  Stopind=0

  while (k<=K & Stopind==0){
    print(k)
    #Enroll qualified patients
    while (nP<sum(nk[1:k])){
      datatmp<-datageneration(0,probfunc)
      ScreenSampleP=rbind(ScreenSampleP,datatmp)
      nScreenP=nScreenP+1
      if (datatmp[1]>=tau[k]){
        SampleP=rbind(SampleP,datatmp)
        nP=nP+1
        if (datatmp[1]>truebio){
          Enrollpos=Enrollpos+1
        }
      }else{
        if (datatmp[1]>truebio){
          Excludpos=Excludpos+1
        }
      }
    }
    while (nT<sum(nk[1:k])){
      datatmp<-datageneration(1,probfunc)
      ScreenSampleT=rbind(ScreenSampleT,datatmp)
      nScreenT=nScreenT+1
      if (datatmp[1]>=tau[k]){
        SampleT=rbind(SampleT,datatmp)
        nT=nT+1
        if (datatmp[1]>truebio){
          Enrollpos=Enrollpos+1
        }
      }else{
        if (datatmp[1]>truebio){
          Excludpos=Excludpos+1
        }
      }
    }

    #Estimate biomarker
    ScreenSample<-rbind(cbind(ScreenSampleP,rep(0,nScreenP)),cbind(ScreenSampleT,rep(1,nScreenT)))
    ScreenSample<-as.data.frame(ScreenSample)
    colnames(ScreenSample)<-c("marker","response","trt")
    rownames(ScreenSample)<-c()
    Sample<-rbind(cbind(SampleP,rep(0,nP)),cbind(SampleT,rep(1,nT)))
    Sample<-as.data.frame(Sample)
    colnames(Sample)<-c("marker","response","trt")
    rownames(Sample)<-c()
    markerest<-markerestf(Sample)
    print(markerest)

    estmean=markerest[1]; estsd=markerest[2]

    if (k<K){
      tautmp=max(0,estmean-t*estsd)

      #Check negative and revise tau if needed
      if (tautmp>0){
        PPoFn<-PPoSF(0,deltaneg,ScreenSample,Sample,tau,tautmp,estmean,nk,k)
        print(PPoFn)
        if (PPoFn<=etaneg){
          tautmp=0
        }
      }

      tau[k+1]=tautmp

      #Check positive to decide next step
      if(estmean==1 & k<Kstar){
        k=k+1
      }else if(estmean==1 & k>=Kstar){
        Stopind=ifelse(k<K,1,0)
        k=K+1
      }else{
        PPoSp<-PPoSF(1,deltapos,ScreenSample,Sample,tau,tautmp,estmean,nk,k)
        print(PPoSp)
        if (PPoSp<=etapos){
          Stopind=ifelse(k<K,1,0)
          k=K+1
        }else{
          k=k+1
        }
      }
    }else{
      k=k+1
    }
  }
  #Final estimate
  if (Stopind==0){
    trteffect=trteffectf(tau,estmean,ScreenSample,Sample,truebio)
    return(c(estmean,estsd,Stopind,dim(ScreenSample)[1],Enrollpos,Excludpos,trteffect))
  }else{
    return(c(NA,NA,Stopind,dim(ScreenSample)[1],NA,NA,rep(NA,14)))
  }
}
