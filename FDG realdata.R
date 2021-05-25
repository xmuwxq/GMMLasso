library(rlist)
library(PhenotypeSimulator)
library(glmnet)
library(CompQuadForm)
library(MASS)
library(Matrix)
library(rrBLUP)
library(matrixcalc)
library(genio)

source("./impute.R")
source("./LinearKernel.R")
source("./chooselambda.R")



###get the phenotype data
phe <- read.csv("./phe.csv")
head(phe)

###get the subjects ID, which means there are a total of 808 subjects
sampleindex <- read.csv("./sampleindex.csv")
head(sampleindex)
dim(sampleindex)
index=sampleindex$index
head(index)

###extract the phenotype data of 808 subjects
newphe=phe[which(phe$PTID %in% index),]
head(newphe)
dim(newphe)
###extract the phenotype data of 808 subjects under baseline
newphe=newphe[newphe$VISCODE=="bl",]
head(newphe)
dim(newphe)

###only save the outcomes and ID for those 808 subjects
fdg=newphe[,c("PTID","FDG")]
av45=newphe[,c("PTID","AV45")]
head(fdg)
head(av45)
dim(fdg)
dim(av45)

###delete the missing values
fdg=fdg[!is.na(fdg$FDG),]
av45=av45[!is.na(av45$AV45),]
head(fdg)
head(av45)
dim(fdg)
dim(av45)


######################################################################################################
#data preprocess######################################################################################

#construct a function to transpose the matrix
transpose<- function(x) {
  tx=t(x)
  tx
}
#construct a function to remove the SNP which has no variants from each gene
remove=function(mydf){
  
  mydf=as.data.frame(mydf) 
  mydf=Filter(function(x) length(unique(x))>1, mydf)
  mydf=as.matrix(mydf)
}

Total=list.load("./Total.rds") #load the original gene data
length(Total)
dim(Total[[1]])
dim(Total[[2]])
dim(Total[[95]])

Total =lapply(Total, transpose)#transpose every gene matrix data
Total =lapply(Total, impute)#impute the NAs in every gene matrix data
Total =lapply(Total, remove)#remove the SNPs which has no variants from every gene matrix data
Total =lapply(Total, standardiseGenotypes)#scale every gene matrix data

#get the 95 gene name
genename=readxl::read_excel("./gene name.xlsx")
head(genename)
dim(genename)
genename=genename$genename



#####################################################################################################
#real data application###############################################################################

Gen=Total 
y=fdg$FDG

totalnum=1:length(y)
test_index <- sample(length(y),100) #randomly abstract 100 sample ID of test data
train_index <- totalnum[-test_index] #get the left sample ID of train data

#get the outcomes of training data
Pheno <- y[train_index]
head(Pheno)

##get the genotype matrix of training data
GenTrain=list()
for (j in 1:length(Total)) {
  GenTrain[[j]]=Gen[[j]][train_index,]
  
}


list_Gen=GenTrain 
list_K=lapply(list_Gen,Fun_Klinear) #get the kernel matrix for each gene

numK=length(list_K)
list_K[[numK+1]]=diag(length(train_index))#add a kernel matrix of the error item into the list


# get A  
n=length(Pheno); 
V=diag(n)
E=eigen(V)$vectors
A=E


lasso_Pheno = c(t(A)%*%Pheno%*%t(Pheno)%*%A)
T=NULL
for (j in 1:length(list_K)) {
  T=cbind(T,c(t(A)%*%list_K[[j]]%*%A))
}

# fit the model
p.fac = c(rep(1,length(Total)),0)
fit_cv <- cv.glmnet(T, lasso_Pheno, alpha=1,lower.limits =0,penalty.factor = p.fac)

lambda=chooselambda(fit_cv,0.1)
fit_best <- glmnet(T, lasso_Pheno, alpha = 1, lambda = lambda, lower.limits =0, penalty.factor = p.fac)
coef=fit_best$beta
rownames(coef)=c(genename,"error item")
coef
mu=coef(fit_best)[1]

#########################################################################################################################
#get the variance matrix


ALL_list_Gen=Gen  
ALL_list_K=lapply(ALL_list_Gen,Fun_Klinear)

Sigma=matrix(0,length(y),length(y))
for (j in 1:length(Total)) {
  Sigma=Sigma + fit_best$beta[j] * ALL_list_K[[j]]
}

Sigma=Sigma+fit_best$beta[length(Total)+1]*diag(length(y))

#########################################################################################################################

Sigma11=Sigma[test_index,test_index]
Sigma12=Sigma[test_index,train_index]
Sigma21=Sigma[train_index,test_index]
Sigma22=Sigma[train_index,train_index]


Pheno_pred=mu+Sigma12 %*% MASS::ginv(Sigma22) %*%(Pheno-mu)

combineout=cbind(Pheno_pred,y[test_index])
head(combineout)
res=Pheno_pred-y[test_index]
MSE=mean(res^2)
COR=cor(Pheno_pred,y[test_index])



