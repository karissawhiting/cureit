coxsplit=function(y, nfolds){
  N=nrow(y)
  tem=data.frame(y, i=seq(N), foldid=0)
  tem=tem[order(y[, "time"], y[, "status"]), ]
  n1=sum(y[, "status"]);n2=N-n1
  
  tem$foldid[tem[, "status"]==1]=sample(rep(seq(nfolds), length=n1))
  tem$foldid[tem[, "status"]==0]=sample(rep(seq(nfolds), length=n2))
  
  foldid=tem$foldid[order(tem$i)]
  return(foldid)
}

coxsplitb=function(y, b, nfolds){
  N=nrow(y)
  tem=data.frame(y, b, i=seq(N), foldid=0)
  tem=tem[order(y[, "time"], y[, "status"]), ]

  for (d in 0:1){
    
    for (batch in 1:max(b)){
      
      tem$foldid[tem[, "status"]==d & tem[,"b"] == batch] = sample(rep(seq(nfolds), length=sum(tem[, "status"]==d & tem[,"b"] == batch)))
      
    }
    
  }
  
  foldid=tem$foldid[order(tem$i)]
  return(foldid)
  
}