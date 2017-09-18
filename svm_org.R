require('quadprog')
data <- read.table("D:/UCF/STA 6106 Statistical Computing/Assignments/Midterm/pb2.txt")
data[,1][data[,1]==2] <- -1
X = data[,2:5]
Y = data[,1]

## Defining the Gaussian kernel
rbf_kernel <- function(x1,x2,gamma){
K<-exp(-(1/gamma^2)*t(x1-x2)%*%(x1-x2))
return(K)
}

svmtrain <- function(X,Y,C=Inf, gamma=1.5,esp=1e-10){

N<-length(Y)
Dm<-matrix(0,N,N)
X<-as.matrix(X);Y<-as.vector(Y)

for(i in 1:N){
for(j in 1:N){
Dm[i,j]<-Y[i]*Y[j]*rbf_kernel(X[i,],X[j,],gamma)
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
# choose one of the support vector to compute b. for safety reason,
# select the one with max alpha
idv<-which.max(alpha)
w<-t(Xv)%*%(alpha*Yv)
b<-Yv[idv]-rbf_kernel(w,Xv[idv,],gamma)
list(alpha=alpha, wstar=w, b=b, nSV=nSV, Xv=Xv, Yv=Yv, gamma=gamma)
}

### Predict the class of an object X



svmpredict <- function(x,model){
alpha<-model$alpha
b<-model$b
Yv<-model$Yv
Xv<-model$Xv
nSV<-model$nSV
gamma<-model$gamma
wstar<-model$wstar
result<-sign(rbf_kernel(wstar,x,gamma)+b)
return(result)
}

model1 =svmtrain(X, Y)

svmpredict(X, model)

