require('quadprog')
X = as.list(numeric())
for (j in 1:15){
  x2 = rnorm(5, 2, 1)
for (i in 1:4){
  x1 =  rnorm(5, 2, 1)
  x2= cbind(x2, x1)
}
  x2=unname(as.matrix(x2), force = TRUE)
  X[j]=list(x2)
}
for (j in 16:25){
  x2 = rnorm(5, 0, 1)
  for (i in 1:4){
    x1 =  rnorm(5, 0, 1)
    x2= cbind(x2, x1)
  }
  x2=unname(as.matrix(x2), force = TRUE)
  X[j]=list(x2)
}


Y = c(rep(1,15), rep(-1,10))

## Defining the Gaussian kernel
rbf_kernel <- function(x1,x2,gamma){
  x1 = as.matrix(x1)
  x2 = as.matrix(x2)
  return(exp(-(1/gamma^2)*sum(diag(t(x1 - x2) %*% (x1 - x2)))))
}





kcalculator <- function(X, ker_par){
  X=as.matrix(X)
  N<-dim(X)[1]
  K<-matrix(0,N,N)
  for(i in 1:N){
    for(j in 1:N){
      K[i,j]<-rbf_kernel(X[i,][[1]],X[j,][[1]],ker_par)
    }
  }
  return(K)
}



bcalculator <- function(Y, X, alpha, ker_par){
  N<-length(Y)
  K = kcalculator(X, ker_par)
  w01=rowSums((alpha*Y)*K)
  w0 = mean(Y-w01)
  
}



svmtrain <- function(X, Y, C=Inf, ker_par =1.5, esp=1e-2){
  N<-length(Y)
  X<-as.matrix(X)
  Y<-as.vector(Y)
  
  K = kcalculator(X, ker_par)
  Dm = (Y %*% t(Y))*K
  Dm<-Dm+diag(N)*1e-8 # adding a very small number to the diag, some trick
  dv<-t(rep(1,N))
  meq<-1
  Am<-cbind(matrix(Y,N),diag(N))
  bv<-rep(0,1+N) # the 1 is for the sum(alpha)==0, others for each alpha_i >= 0
  if(C!=Inf){
    # an upper bound is given
    Am<-cbind(Am,-1*diag(N))
    bv<-c(cbind(matrix(bv,1),matrix(rep(-C,N),1)))
  }
  alpha_org <- solve.QP(Dm,dv,Am,bvec=bv, meq=meq)$solution
  indx<-which(alpha_org>esp,arr.ind=TRUE)
  alpha<-alpha_org[indx]
  nSV<-length(indx)
  if(length(indx)==0){
    throw("QP is not able to give a solution for these data points")
  }
  Xv<-X[indx,]
  Yv<-as.vector(Y[indx])
  # choose one of the support vector to compute b. for safety reason,
  # select the one with max alpha
  
  b = bcalculator(Yv, Xv, alpha, ker_par)
  
  return(list(alpha=alpha, b=b, nSV=nSV, Xv=Xv, Yv=Yv, ker_par=ker_par))
}

### Predict the class of an object X

model1 = svmtrain(X, Y, ker_par = 34)
model1

svmpredict <- function(x,model){
  x = as.matrix(x)
  kernel = model$kernel
  ker_per = model$ker_par
  alpha<-model$alpha
  b<-model$b
  Yv<-model$Yv
  Xv<-as.matrix(model$Xv)
  ker_par<-model$ker_par
  # wstar<-model$wstar
  result = as.vector(rep(0,dim(x)[1]))
  for (k in 1:dim(x)[1]){
    sum = 0
      for (i in 1 : length(alpha)){
        sum = sum + alpha[i] * Yv[i] * rbf_kernel(Xv[i,][[1]],x[k,][[1]],ker_per)
      }
       result[k]<-sign(sum + b)

  }
  margin = sum(alpha)
  return(list(result = result, margin = sqrt(1/margin), primal = margin))
}



svmpredict(model1$Xv, model1)
svmpredict(X, model1)
