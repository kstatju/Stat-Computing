require('quadprog')


##################################################################
#####defining the functions

marr_kernel <- function(x1,x2,gamma){
    K<- prod((1-((x1-x2)/gamma)^2)*exp(-1/2*((x1-x2)/gamma)^2))
    return(K)
}

mor_kernel <- function(x1,x2,gamma){   
    K<- prod((cos(1.75*(x1-x2)/gamma))*exp(-1/2*((x1-x2)/gamma)^2))
    return(K)
}

morrbf_kernel <- function(x1,x2,gamma){  
    k1 = prod((cos(1.75*(x1-x2)/gamma[1]))*exp(-1/2*((x1-x2)/gamma[1])^2))
    K<- exp(-gamma[2]*(2-2*k1))
    return(K)
  }

matrix_kernel <- function(x1,x2,gamma){   
  K<- prod((cos(1.75*(x1-x2)/gamma))*exp(-1/2*((x1-x2)/gamma)^2))
  return(K)
}
  
  
  
  
###########################################################################

svmtrain <- function(X,Y,C=Inf, gamma=1.5,esp=1e-10){

N<-length(Y)
Dm<-matrix(0,N,N)
X<-as.matrix(X);Y<-as.vector(Y)

for(i in 1:N){
for(j in 1:N){
Dm[i,j]<-Y[i]*Y[j]*morrbf_kernel(X[i,],X[j,],gamma)
}
}
Dm<-Dm+diag(N)*1e-12 # adding a very small number to the diag, some trick
dv<-t(rep(1,N))
meq<-1
Am<-cbind(matrix(Y,N),diag(N))
bv<-rep(0,1+N) # the 1 is for the sum(alpha)==0, others for each alpha_i >= 0
if(C!=Inf){
# an upper bound is given
Am<-cbind(Am,-1*diag(N))
bv<-c(cbind(matrix(bv,1),matrix(rep(-C,N),1)))
}
alpha_org<-solve.QP(Dm,dv,Am,meq=meq,bvec=bv)$solution
indx<-which(alpha_org>esp,arr.ind=TRUE)
alpha<-alpha_org[indx]
nSV<-length(indx)
if(nSV==0){
throw("QP is not able to give a solution for these data points")
}
Xv<-X[indx,]
Yv<-Y[indx]
Yv<-as.vector(Yv)
## choose one of the support vector to compute b. for safety reason,


#--- Averaging over all support vectors ---

sum0 = 0
for (i in 1 : length(alpha)){
	sum = 0
	for (j in 1 : length(alpha)){
		temp = alpha[j] * Yv[j] * morrbf_kernel(Xv[j,],Xv[i,],gamma)
		sum = sum + temp
	}
	sum0 = sum0 + Yv[i] - sum
}
b = sum0 / length(alpha)

list(alpha=alpha, alpha_org = alpha_org, b=b, nSV=nSV, Xv=Xv, Yv=Yv, gamma=gamma)
}

### Predict the class of an object X


svmpredict <- function(x,X,Y,C0=Inf, gamma0=1.5,esp0=1e-10){
model = svmtrain(X,Y,C=C0, gamma=gamma0,esp=esp0)
alpha<-model$alpha
b<-model$b
Yv<-model$Yv
Xv<-model$Xv
nSV<-model$nSV
gamma<-model$gamma
sum = 0
for (i in 1 : length(alpha)){
	temp = alpha[i] * Yv[i] * morrbf_kernel(Xv[i,],x,gamma)
	sum = sum + temp
}
hat = sum + b
result<-sign(hat)
return(result)
}

svmtrainError <- function(X,Y,C0=Inf, gamma0=1.5,esp0=1e-10){
error = 0
for (i in 1 : length(Y)){
	Ypred = svmpredict(X[i, ],X,Y,C0=C0, gamma0=gamma0,esp0=esp0)
	if (Ypred != Y[i]) error = error + 1
}
return(error)
}

setwd('D:/UCF/6106/Exam/1')
A = read.table("D:/UCF/6106/Exam/1/pb2.txt", head = FALSE)
X = A[, 2 : 5]
Y0 = A[, 1]
Y = c()
for (i in 1 : length(Y0)){
	if (Y0[i] == 1) Y[i] = 1
	if (Y0[i] == 2) Y[i] = -1
}


z = c(18, 17, 33, 26)

gamma = c(7,1.5)
#---support vector:
res = svmtrain(X,Y,C=Inf, gamma = gamma)
res$alpha_org
res$nSV
#---prediction:
svmpredict(z, X, Y, gamma0 = gamma)
#---train errors:
svmtrainError(X,Y,C=Inf, gamma0 = gamma)



####d
##create 25 matrices
X = list()
for (j in 1:13){
  x2 = rnorm(4, 0, 1)
  for (i in 1:4){
    x1 =  rnorm(4, 2, 1)
    x2= cbind(x2, x1)
  }
  x2=unname(as.matrix(x2), force = TRUE)
  X[j]=list(x2)
}
for (j in 14:25){
  x2 = rnorm(4, 2, 1)
  for (i in 1:4){
    x1 =  rnorm(4, 0, 1)
    x2= cbind(x2, x1)
  }
  x2=unname(as.matrix(x2), force = TRUE)
  X[j]=list(x2)
}
Y = c(rep(1,13), rep(-1, 12))

## Defining the Matrix Gaussian kernel
rbf_kernel <- function(x1,x2,gamma){
  #x1, x2 two matrices
  tr = sum(diag(t(x1 - x2) %*% (x1 - x2)))
  K<-exp(-(1/gamma^2)*tr)
  return(K)
}

## Defining the Matrix kernel
morlet_kernel <- function(x1,x2,gamma){  
  ## Defining the Morlet-RBF wavelet kernel
  K = prod((cos(1.75*(x1-x2)/gamma))*exp(-1/2*((x1-x2)/gamma)^2))
  return(K)
  
}

svmtrain <- function(X,Y,C=Inf, gamma=1.5,esp=1e-10){
  #matirx form, X is a list containing all the explanatory matrices
  N<-length(Y)
  Dm<-matrix(0,N,N)
  Y<-as.vector(Y)
  
  for(i in 1:N){
    for(j in 1:N){
      Dm[i,j]<-Y[i]*Y[j]*morlet_kernel(X[[i]],X[[j]],gamma)
    }
  }
  Dm<-Dm+diag(N)*1e-12 # adding a very small number to the diag, some trick
  dv<-t(rep(1,N))
  meq<-1
  Am<-cbind(matrix(Y,N),diag(N))
  bv<-rep(0,1+N) # the 1 is for the sum(alpha)==0, others for each alpha_i >= 0
  if(C!=Inf){
    # an upper bound is given
    Am<-cbind(Am,-1*diag(N))
    bv<-c(cbind(matrix(bv,1),matrix(rep(-C,N),1)))
  }
  alpha_org<-solve.QP(Dm,dv,Am,meq=meq,bvec=bv)$solution
  indx<-which(alpha_org>esp & alpha_org < C, arr.ind=TRUE)
  alpha<-alpha_org[indx]
  nSV<-length(indx)
  if(nSV==0){
    throw("QP is not able to give a solution for these data points")
  }
  Xv<-X[indx]
  Yv<-Y[indx]
  
  #--- Averaging over all support vectors ---
  
  sum0 = 0
  for (i in 1 : length(alpha)){
    sum = 0
    for (j in 1 : length(alpha)){
      temp = alpha[j] * Yv[j] * morlet_kernel(Xv[[j]],Xv[[i]],gamma)
      sum = sum + temp
    }
    sum0 = sum0 + Yv[i] - sum
  }
  b = sum0 / length(alpha)
  
  list(alpha_org=alpha_org, b=b)
}

### Predict the class of an object X


svmpredict <- function(x,X,Y,C0=Inf, gamma0=1.5,esp0=1e-10){
  #matrix form
  model = svmtrain(X,Y,C=C0, gamma=gamma0,esp=esp0)
  alpha_org<-model$alpha_org
  b<-model$b
  
  sum = 0
  for (i in 1 : length(Y)){
    temp = alpha_org[i] * Y[i] * morlet_kernel(X[[i]],x,gamma0)
    sum = sum + temp
  }
  hat = sum + b
  result<-sign(hat)
  return(result)
}

svmtrain(X, Y)
svmpredict(X[[1]],X,Y)














