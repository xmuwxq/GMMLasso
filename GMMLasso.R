

newGMMLasso=function(y=y,Gen=Gen,Fun_K=Fun_Klinear,index=index){
  
  test_index <- index  ##test index
  train_index <- c(1:length(y))[-index]  ##train index
  
  Pheno <- y[train_index]  ##training Phenotypes
  
  GenTrain=list()
  
  for (j in 1:length(Gen)) {
    
    GenTrain[[j]]=Gen[[j]][train_index,]   ##training Genomic matrix
    
  }
  
  list_K=lapply(GenTrain,Fun_K) ##get kernel function
  
  numK=length(list_K)
  list_K[[numK+1]]=diag(length(train_index))
  
  n=length(Pheno);   
  V=diag(n)
  E=eigen(V)$vectors
  A=E                       ##get A matrix
  
  
  lasso_Pheno = c(t(A)%*%Pheno%*%t(Pheno)%*%A)   ##training Phenotypes used in GMMLasso
  T=NULL
  for (j in 1:length(list_K)) {
    T=cbind(T,c(t(A)%*%list_K[[j]]%*%A))    ##training Genomic matrix used in GMMLasso
  }
  
  p.fac = c(rep(1,length(Gen)),0)
  fit_cv <- cv.glmnet(T, lasso_Pheno, alpha=1,lower.limits =0,penalty.factor = p.fac)     ## fit the model
  
  lambda=ownse(fit_cv,0.1)      ##search the optimal lambda
  
  fit_best <- glmnet(T, lasso_Pheno, alpha = 1, lambda = lambda, lower.limits =0, penalty.factor = p.fac)       ## fit the model
  beta= t(as.matrix(fit_best$beta))   #get the parameter estimations 
  mu=coef(fit_best)[1]                #get the constant term
  
  
  
  ALL_list_K=lapply(Gen,Fun_K)
  
  Sigma=matrix(0,length(y),length(y))                
  
  for (j in 1:length(Gen)) {
    Sigma=Sigma + fit_best$beta[j] * ALL_list_K[[j]]
  }
  
  Sigma=Sigma+fit_best$beta[length(list_K)]*diag(length(y))           ##get the variance matrix
  
  
  Sigma11=Sigma[test_index,test_index]
  Sigma12=Sigma[test_index,train_index]
  Sigma21=Sigma[train_index,test_index]
  Sigma22=Sigma[train_index,train_index]
  
  
  
  Pheno_pred=mu+Sigma12 %*% MASS::ginv(Sigma22) %*%(Pheno-mu)
  
  combineout=cbind(pred=Pheno_pred,true=y[test_index])
  
  out=list(out=combineout,beta=beta)
  
  out
}






