pfunc<-function(trt,x,probfunc){
  if (trt==0){
    probv=eval(parse(text = probfunc[1]))
  }else{
    probv=eval(parse(text = probfunc[2]))
  }
  return(probv)
}
