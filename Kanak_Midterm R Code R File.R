#########################################################################
#     Problem 1
#########################################################################

rnormalBX <- function(n=1, mean = 0, var = 1, seed = NULL){
  
  if(!is.null(seed))  set.seed(seed)  # set seed for the random number
  
  # generate two uniform and then calculate r and theta and finally generate normal random number
  
  u = runif(n, 0, 1)
  v = runif(n, 0, 1)
  r = sqrt(-2*log(v))
  theta = 2*pi*u
  norm = cbind(r*cos(theta), r*sin(theta))
  if(n != 1) {norm = (norm - colMeans(norm))/sqrt(diag(var(norm)))}
  
  # transform normal random number with specified mean and variance 
  
  if (mean != 0 || var != 1){
    norm = sqrt(var)*norm + mean
    
  }
  return(unname(norm, force = TRUE))
}

# # calling rnormalBX function to generate normal random number using Box-Muller algorithm
# x = rnormalBX(n = 1000, mean = 0, var = 1, seed = 123)
# summary(x)


rcauchydist <- function(n=1,  seed = NULL){
  x = rnormalBX(n = n, mean = 0, var = 1, seed = seed)
  return(x[,1]/x[,2])
}

x = rcauchydist(1000, seed = 4578)
#Part 1 Independence of two normal random variable generated using Box-Muller algorithm



# Density Curve for the normal sample. Function densityCurve has two parameters, data = data set for the density curve, dist = distribution.


densityCurve <- function(data , nbreaks = NULL){
  
  if(is.null(nbreaks)){ nbreaks = pretty(data)}
  title = "Histogram with Cauchy Curve"
  his  <- hist(data, breaks = nbreaks, freq = FALSE, col="red", xlab="x", main=title, xlim = c(-30,30))
  xfit <- his$mids 
  yfit <- dcauchy(xfit, location = 0, scale =  1)
  lines(xfit, yfit, col="blue", lwd=2)
  
}


# Calling densityCurve function to create density curve.
densityCurve(x, nbreaks = 500)

# problem 1 part 2

cauchyMeanForEachSample <- function(samSize = 1:100, seed = NULL){
  if(!is.null(seed))  set.seed(seed)
  return(apply(as.array(samSize), 1, function(x) mean(rcauchydist(x))))
}

seq = seq.int(1,10^4, 50)
cmean = cauchyMeanForEachSample(seq, seed = 245)
plot(seq, cmean)



# problem 1 part 3

cauchyMeanSim <- function(samSize = 100, distSize = 1000 , seed = NULL){
  if(!is.null(seed))  set.seed(seed)
  meandist =  replicate(distSize, mean(rcauchydist(samSize)))
  
  return(meandist)
}

cmean = cauchyMeanSim(40, distSize = 2000, seed = 365)
densityCurve(cmean, nbreaks = 800)



cmean = cauchyMeanSim(100, distSize = 2000, seed = 365)
densityCurve(cmean, nbreaks = 1300)


cmean = cauchyMeanSim(500, distSize = 2000, seed = 4125)
densityCurve(cmean, nbreaks = 1500)




#########################################################################
#     Problem 2
#########################################################################


require('quadprog')
data <- read.table("C:/Users/ka746940/Desktop/UCF/STA 6106 - Statistical Computing/Assignments/Midterm/pb2.txt")
#data <- read.table("D:/UCF/STA 6106 Statistical Computing/Assignments/Midterm/pb2.txt")
data[,1][data[,1]==2] <- -1
X = data[,2:5]
mmean = colMeans(X)
cvar = diag(var(X))
Y = data[,1]
X1 = scale(X)
## Defining the Gaussian kernel
rbf_kernel <- function(x1,x2,ker_par){
  x1 = as.matrix(x1)
  x2 = as.matrix(x2)
  K<-exp(-(1/(ker_par^2))*t(x1-x2)%*%(x1-x2))
  return(K)
}


poli_kernel <- function(x1, x2, c, d){
  K<- (t(as.matrix(x1)) %*% as.matrix(x2) + c)^d
  return(K)
}


kcalculator <- function(X, kernel, ker_par){
  X=as.matrix(X)
  N<-dim(X)[1]
  K<-matrix(0,N,N)
  if (toupper(kernel)== "GAUSSIAN"){
    for(i in 1:N){
      for(j in 1:N){
        K[i,j]<-rbf_kernel(X[i,],X[j,],ker_par)
      }
    }
  }
  if (toupper(kernel)== "POLYNOMIAL"){
    for(k in 1:N){
      for(l in 1:N){
        K[k,l]<-poli_kernel(X[k,],X[l,],ker_par[1], ker_par[2])
      }
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
  N<-length(Y)
  X<-as.matrix(X)
  Y<-as.vector(Y)
  
  K = kcalculator(X, kernel, ker_par)
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
  alpha_org<-solve.QP(Dm,dv,Am,bvec=bv, meq=meq)$solution
  indx<-which(alpha_org>esp,arr.ind=TRUE)
  alpha<-alpha_org[indx]
  nSV<-length(indx)
  if(length(indx)==0){
    throw("QP is not able to give a solution for these data points")
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
    if (toupper(kernel)== "GAUSSIAN"){
      for (i in 1 : length(alpha)){
        s1 = alpha[i] * Yv[i] * rbf_kernel(Xv[i,],x[k,],ker_per)
        sum = sum + s1
      }
      result[k]<-sign(sum + b)
    }
    
    if (toupper(kernel)== "POLYNOMIAL"){
      for (i in 1 : length(alpha)){
        s1 = alpha[i] * Yv[i] * poli_kernel(Xv[i,],x[k,],ker_per[1], ker_per[2])
        sum = sum + s1
      }
      result[k]<-sign(sum + b)
    }
  }
  
  return(result)
}


model1 = svmtrain(X1, Y, kernel = "Polynomial", ker_par = c(23,2))
model1
model =svmtrain(X, Y, kernel = "Gaussian", ker_par = 1.5)
model
p1 = svmpredict(X1, model1)
cat("Training Erro", (length(Y)-sum(Y==p1))/length(Y)*100, "%")
p2 = svmpredict(X, model)
cat("Training Erro", (length(Y)-sum(Y==p2))/length(Y)*100, "%")
a = matrix(c(18,17,33,26), 1,4,byrow = TRUE)
a1 = (a-mmean)/cvar
p1 = svmpredict(a1, model1)
cat("Prediction", p1)
p2 = svmpredict(a, model)
cat("Prediction", p2)


# part b problem 2

sigm = seq(.5, 100, length = 100)
itsig = matrix(0, length(sigm), 3)

for (i in 1:length(sigm)){
  fit = svmtrain(X, Y, kernel = "Gaussian", ker_par = sigm[i])
  pred = svmpredict(X, fit)
  error = (length(Y)-sum(Y==pred))/length(Y)*100
  nSV = fit$nSV
  itsig[i,] = c(error, nSV, sigm[i])
}


plot(itsig[,3], itsig[,2], xlab = "Sigma", ylab = "Number of Support Vectors")
hist(itsig[,2])



#########################################################################
#     Problem 3
#########################################################################


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

