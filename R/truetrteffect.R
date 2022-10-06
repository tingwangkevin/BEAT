truetrteffect<-function(chat,ctrue){
  xtmp=seq(chat,1,by=0.001)
  poschat=mean(pfunc(1,xtmp,probfunc))-mean(pfunc(0,xtmp,probfunc))

  xtmp=seq(ctrue,1,by=0.001)
  posctrue=mean(pfunc(1,xtmp,probfunc))-mean(pfunc(0,xtmp,probfunc))

  xtmp=seq(0,chat,by=0.001)
  negchat=mean(pfunc(1,xtmp,probfunc))-mean(pfunc(0,xtmp,probfunc))

  xtmp=seq(0,ctrue,by=0.001)
  negctrue=mean(pfunc(1,xtmp,probfunc))-mean(pfunc(0,xtmp,probfunc))

  xtmp=seq(0,1,by=0.001)
  overall=mean(pfunc(1,xtmp,probfunc))-mean(pfunc(0,xtmp,probfunc))

  return(c(poschat,posctrue,negchat,negctrue,overall))
}
