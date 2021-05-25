chooselambda=function(object,d){
  cvm=object$cvm
  cvsd=object$cvsd
  limit=min(cvm+d*cvsd)
  maxless <- max(cvm[cvm <= limit])
  index=which(cvm == maxless)
  object$lambda[index]
}

