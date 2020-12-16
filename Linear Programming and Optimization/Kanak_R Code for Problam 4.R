#######################################################################
#
#       SVM
#
#######################################################################



require('quadprog')
data <- read.table("C:/Users/ka746940/Desktop/UCF/STA 6106 - Statistical Computing/Assignments/Final/pb22.txt")
#data <- read.table("D:/UCF/STA 6106 Statistical Computing/Assignments/Midterm/pb2.txt")
data[,1][data[,1]==2] <- -1
X = data[,2:5]
mmean = colMeans(X)
cvar = diag(var(X))
Y = data[,1]
X1 = as.data.frame(scale(X))


kernel_cal <- function(x1,x2,ker_par, kernel = "GAUSSIAN"){
  if(toupper(kernel)=="GAUSSIAN"){     
    ## Defining the Gaussian kernel
    K<-exp(-(1/(ker_par^2))*t(x1-x2)%*%(x1-x2))
    return(K)
    
  } else if(toupper(kernel)=="POLYNOMIAL"){ 
    ## Defining the Polinomial kernel
    K<- (t(as.matrix(x1)) %*% as.matrix(x2) + ker_par[1])^ker_par[2]
    return(K)
    
  } else if(toupper(kernel)=="MARR"){  
    ## Defining the  Marr wavelet kernel
    K<- prod((1-((x1-x2)/ker_par)^2)*exp(-1/2*((x1-x2)/ker_par)^2))
    return(K)
    
  } else if(toupper(kernel)=="MORLET"){   
    ## Defining the  Morlet wavelet kernel
    K<- prod((cos(1.75*(x1-x2)/ker_par))*exp(-1/2*((x1-x2)/ker_par)^2))
    return(K)
    
  } else if(toupper(kernel)=="MORLET-RBF"){   
    ## Defining the Morlet-RBF wavelet kernel
    k1 = prod((cos(1.75*(x1-x2)/ker_par[1]))*exp(-1/2*((x1-x2)/ker_par[1])^2))
    K<- exp(-ker_par[2]*(2-2*k1))
    return(K)
    
  }
}





kcalculator <- function(X, kernel, ker_par){
  X=as.matrix(X)
  N<-dim(X)[1]
  K<-matrix(0,N,N)
  for(i in 1:N){
    for(j in 1:N){
      K[i,j]<-kernel_cal(X[i,],X[j,],ker_par, kernel = kernel)
    }
  }
  return(K)
}


bcalculator <- function(Y, X, alpha, kernel, ker_par){
  N<-length(Y)
  K = kcalculator(X, kernel, ker_par)
  w01=rowSums((alpha*Y)*K)
  w0 = mean(Y-w01)
  
}



svmtrain <- function(X, Y, C=Inf, kernel = "Gaussian", ker_par =1.5, esp=1e-2){
  if(!(toupper(kernel) %in% c("GAUSSIAN","POLYNOMIAL", "MARR", "MORLET", "MORLET-RBF"))){
    cat("There is no Kernel named ", kernel)
    return()
  }
  N<-length(Y)
  X<-as.matrix(X)
  Y<-as.vector(Y)
  
  K = kcalculator(X, kernel, ker_par)
  Dm = (Y %*% t(Y))*K
  Dm<-Dm+diag(N)*1e-5 # adding a very small number to the diag, some trick
  dv<-t(rep(1,N))
  meq<-1
  Am<-cbind(matrix(Y,N),diag(N))
  bv<-rep(0,1+N) # the 1 is for the sum(alpha)==0, others for each alpha_i >= 0
  if(C!=Inf){
    # an upper bound is given
    Am<-cbind(Am,-1*diag(N))
    bv<-c(cbind(matrix(bv,1),matrix(rep(-C,N),1)))
  }
  alpha_org<-solve.QP(Dm,dv,Am,bvec=bv, meq=meq)$solution
  indx<-which(alpha_org>esp,arr.ind=TRUE)
  alpha<-alpha_org[indx]
  nSV<-length(indx)
  if(length(indx)==0){
    stop("QP is not able to give a solution for these data points")
  }
  Xv<-X[indx,]
  Yv<-as.vector(Y[indx])
  w<-unname(t(Xv)%*%(alpha*Yv), force = TRUE)
  # choose one of the support vector to compute b. for safety reason,
  # select the one with max alpha
  
  b = bcalculator(Yv, Xv, alpha, kernel, ker_par)
  
  return(list(alpha=alpha, wstar=w, b=b, nSV=nSV, Xv=Xv, Yv=Yv, kernel = kernel ,ker_par=ker_par))
}

### Predict the class of an object X



svmpredict <- function(x,model){
  x = as.matrix(x)
  kernel = model$kernel
  ker_per = model$ker_par
  alpha<-model$alpha
  b<-model$b
  Yv<-model$Yv
  Xv<-model$Xv
  ker_par<-model$ker_par
  # wstar<-model$wstar
  result = as.vector(rep(0,dim(x)[1]))
  for (k in 1:dim(x)[1]){
    sum = 0
    for (i in 1 : length(alpha)){
      s1 = alpha[i] * Yv[i] * kernel_cal(Xv[i,],x[k,],ker_per, kernel)
      sum = sum + s1
    }
    result[k]<-sign(sum + b)
    
    
  }
  
  return(result)
}



model =svmtrain(X, Y, kernel = "gaussian", ker_par = 2.5)
model
model1 = svmtrain(X1, Y, kernel = "Polynomial", ker_par = c(23,2))
model1
model2 = svmtrain(X, Y, kernel = "marr", ker_par = 5)
model2
model3 = svmtrain(X, Y, kernel = "morlet", ker_par = 10)
model3
model4 = svmtrain(X, Y, kernel = "MORLET-RBF", ker_par = c(10,2))
model4
p = svmpredict(X, model)
cat("Training Erro", (length(Y)-sum(Y==p))/length(Y)*100, "%")
p1 = svmpredict(X1, model1)
cat("Training Erro", (length(Y)-sum(Y==p1))/length(Y)*100, "%")
p2 = svmpredict(X, model2)
cat("Training Erro", (length(Y)-sum(Y==p2))/length(Y)*100, "%")
p3 = svmpredict(X, model3)
cat("Training Erro", (length(Y)-sum(Y==p3))/length(Y)*100, "%")
p4 = svmpredict(X, model4)
cat("Training Erro", (length(Y)-sum(Y==p4))/length(Y)*100, "%")


a = matrix(c(18,17,33,26), 1,4,byrow = TRUE)
a1 = (a-mmean)/cvar

pr = svmpredict(a, model)
cat("Prediction", pr)

pr1 = svmpredict(a1, model1)
cat("Prediction", pr1)

pr2 = svmpredict(a, model2)
cat("Prediction", pr2)

pr3 = svmpredict(a, model3)
cat("Prediction", pr3)

pr4 = svmpredict(a, model4)
cat("Prediction", pr4)


# part b problem 2
kk = c("GAUSSIAN", "MARR", "MORLET")
require(caTools)
set.seed(101) 
sample = sample.split(Y, SplitRatio = .75)
X_tr = subset(X1, sample == TRUE)
X_ts = subset(X1, sample == FALSE)
Y_tr = subset(Y, sample == TRUE)
Y_ts = subset(Y, sample == FALSE)             

cseq = seq(2, 500, length = 10)
kerseq = seq(.1, 200, length = 100)
combbi = expand.grid(cseq, kerseq)
kkmin = list(kernel = character(), NSVM = numeric(), C = numeric(), Sigma = numeric(), Tr_e = numeric(), Ts_e = numeric())
for (j in 1:length(kk)){
  min_k = numeric()
  t_e = numeric()
  ts_e = numeric()
  for (i in 1:dim(combbi)[1]){
    min_k1 = svmtrain(X=X_tr, Y=Y_tr, C=combbi[i,1], kernel = kk[j], ker_par =combbi[i,2])
    min_k[i] = min_k1$nSV
    p = svmpredict(X_tr, min_k1)
    t_e[i] = (length(Y_tr)-sum(Y_tr==p))/length(Y)*100
    t_p = svmpredict(X_ts, min_k1)
    ts_e[i] = (length(Y_ts)-sum(Y_ts==t_p))/length(Y)*100
  }
  pp = which(ts_e == min(ts_e))
  pp1 = which(min_k == min(min_k[pp]))
  pp2 = combbi[pp1,][which(combbi[pp1,1] == min(combbi[pp1,1])),]
  kkmin$kernel[j] = kk[j]
  kkmin$NSVM[j] = min_k[pp1[1]]
  kkmin$C[j] = pp2[1,1]
  kkmin$Sigma[j] = pp2[1,2]
  kkmin$Ts_e = ts_e[pp[1]]
  kkmin$Tr_e = t_e[pp[1]]
}

kkmin

for (i in 1:3){
  trn = svmtrain(X=X1, Y=Y, C=kkmin$C[i], kernel = kkmin$kernel[i], ker_par =kkmin$Sigma[i])
  p = svmpredict(X1, trn)
  cat("Training Erro", (length(Y)-sum(Y==p))/length(Y)*100, "%")
  a = matrix(c(18,17,33,26), 1,4,byrow = TRUE)
  a1 = scale(a, center = mmean, scale = cvar)
  cat("Prediction", svmpredict(a1, trn))
  
}





#######################################################################
#
#       SMM
#
#######################################################################

require('quadprog')
set.seed(0)
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



kernel_cal <- function(x1,x2,ker_par, kernel = "GAUSSIAN"){
  if(toupper(kernel)=="GAUSSIAN"){     
    ## Defining the Gaussian kernel
    return(exp(-(1/ker_par^2)*sum(diag(t(x1 - x2) %*% (x1 - x2)))))
    
  } else if(toupper(kernel)=="MORLET"){   
    ## Defining the Morlet-RBF wavelet kernel
    K = prod((cos(1.75*(x1-x2)/ker_par))*exp(-1/2*((x1-x2)/ker_par)^2))
    return(K)
    
  }
}


## Defining the Gaussian kernel





kcalculator <- function(X, ker_par, kernel){
  X=as.matrix(X)
  N<-dim(X)[1]
  K<-matrix(0,N,N)
  for(i in 1:N){
    for(j in 1:N){
      K[i,j]<-kernel_cal(X[i,][[1]],X[j,][[1]],ker_par=ker_par, kernel=kernel)
    }
  }
  return(K)
}



bcalculator <- function(Y, X, alpha, ker_par, kernel){
  N<-length(Y)
  K = kcalculator(X, ker_par, kernel)
  w01=rowSums((alpha*Y)*K)
  w0 = mean(Y-w01)
  
}



smmtrain <- function(X, Y, C=Inf, ker_par =1.5, kernel = "GAUSSIAN", esp=1e-2){
  N<-length(Y)
  X<-as.matrix(X)
  Y<-as.vector(Y)
  
  K = kcalculator(X, ker_par, kernel)
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
  
  b = bcalculator(Yv, Xv, alpha, ker_par, kernel)
  
  return(list(alpha=alpha, b=b, nSV=nSV, Xv=Xv, Yv=Yv, ker_par=ker_par, kernel = kernel))
}

### Predict the class of an object X

model1 = smmtrain(X, Y, ker_par = 34, kernel = "MORLET")
model1

smmpredict <- function(x,model){
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
      sum = sum + alpha[i] * Yv[i] * kernel_cal(Xv[i,][[1]],x[k,][[1]],ker_per, kernel)
    }
    result[k]<-sign(sum + b)
    
  }
  margin = sum(alpha)
  return(list(pred = result, margin = sqrt(1/margin), primal = margin, kernel = kernel))
}



p = smmpredict(X, model1)
cat("Training Erro", (length(Y)-sum(Y==p$pred))/length(Y)*100, "%")


# part b problem 2
kk = c("GAUSSIAN", "MORLET")
require(caTools)
set.seed(101) 
sample = sample.split(Y, SplitRatio = .75)
X_tr = subset(X, sample == TRUE)
X_ts = subset(X, sample == FALSE)
Y_tr = subset(Y, sample == TRUE)
Y_ts = subset(Y, sample == FALSE)             

cseq = seq(2, 500, length = 5)
kerseq = seq(.1, 200, length = 5)
combbi = expand.grid(cseq, kerseq)
kkmin = list(kernel = character(), NSVM = numeric(), C = numeric(), Sigma = numeric(), Tr_e = numeric(), Ts_e = numeric())
for (j in 1:length(kk)){
  min_k = numeric()
  t_e = numeric()
  ts_e = numeric()
  for (i in 1:dim(combbi)[1]){
    min_k1 = smmtrain(X=X_tr, Y=Y_tr, C=combbi[i,1], kernel = kk[j], ker_par =combbi[i,2])
    min_k[i] = min_k1$nSV
    p = smmpredict(X_tr, min_k1)
    t_e[i] = (length(Y_tr)-sum(Y_tr==p$pred))/length(Y)*100
    t_p = smmpredict(X_ts, min_k1)
    ts_e[i] = (length(Y_ts)-sum(Y_ts==t_p$pred))/length(Y)*100
  }
  pp = which(ts_e == min(ts_e))
  pp1 = which(min_k == min(min_k[pp]))
  pp2 = combbi[pp1,][which(combbi[pp1,1] == min(combbi[pp1,1])),]
  kkmin$kernel[j] = kk[j]
  kkmin$NSVM[j] = min_k[pp1[1]]
  kkmin$C[j] = pp2[1,1]
  kkmin$Sigma[j] = pp2[1,2]
  kkmin$Ts_e = ts_e[pp[1]]
  kkmin$Tr_e = t_e[pp[1]]
}

kkmin

for (i in 1:2){
  trn = smmtrain(X=X, Y=Y, C=kkmin$C[i], kernel = kkmin$kernel[i], ker_par =kkmin$Sigma[i])
  p = smmpredict(X, trn)
  cat("Training Erro", (length(Y)-sum(Y==p$pred))/length(Y)*100, "%")
  a = matrix(c(18,17,33,26), 1,4,byrow = TRUE)
  a1 = as.matrix(list(X[[1]]))
  cat("Prediction", smmpredict(a1, trn)$pred)
  
}