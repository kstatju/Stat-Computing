require('quadprog')
data <- read.table("D:/UCF/STA 6106 Statistical Computing/Assignments/Final/pb2.txt")
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



model =svmtrain(X, Y, kernel = "gaussian", ker_par = 1.5)
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

sigm = seq(.5, 500, length = 1000)
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
