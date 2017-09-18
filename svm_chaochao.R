require('quadprog')
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
## choose one of the support vector to compute b. for safety reason,
## select the one with max alpha
#idv<-which.max(alpha)
#w<-t(Xv)%*%(alpha*Yv)
#b<-Yv[idv]-rbf_kernel(w,Xv[idv,],gamma)

#--- Averaging over all support vectors ---

sum0 = 0
for (i in 1 : length(alpha)){
	sum = 0
	for (j in 1 : length(alpha)){
		temp = alpha[j] * Yv[j] * rbf_kernel(Xv[j,],Xv[i,],gamma)
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
	temp = alpha[i] * Yv[i] * rbf_kernel(Xv[i,],x,gamma)
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


A = read.table("pb2.txt", head = FALSE)
X = A[, 2 : 5]
Y0 = A[, 1]
Y = c()
for (i in 1 : length(Y0)){
	if (Y0[i] == 1) Y[i] = 1
	if (Y0[i] == 2) Y[i] = -1
}


z = c(18, 17, 33, 26)

gamma = 1.5
#---support vector:
res = svmtrain(X,Y,C=Inf, gamma = gamma)
res$alpha_org
res$nSV
#---prediction:
svmpredict(z, X, Y, gamma0 = gamma)
#---train errors:
svmtrainError(X,Y,C=Inf, gamma0 = gamma)
















#---polynomial kernel---

pol_kernel <- function(x1,x2,gamma){
K<-(1 + t(x1)%*%(x2))^gamma
return(K)
}


svmPoltrain <- function(X,Y,C=Inf, gamma=1.5,esp=1e-10){

N<-length(Y)
Dm<-matrix(0,N,N)
X<-as.matrix(X);Y<-as.vector(Y)

for(i in 1:N){
for(j in 1:N){
Dm[i,j]<-Y[i]*Y[j]*pol_kernel(X[i,],X[j,],gamma)
}
}
Dm<-Dm+diag(N)*1e-1 # adding a very small number to the diag, some trick
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
## select the one with max alpha
#idv<-which.max(alpha)
#w<-t(Xv)%*%(alpha*Yv)
#b<-Yv[idv]-pol_kernel(w,Xv[idv,],gamma)

#--- Averaging over all support vectors ---

sum0 = 0
for (i in 1 : length(alpha)){
	sum = 0
	for (j in 1 : length(alpha)){
		temp = alpha[j] * Yv[j] * pol_kernel(Xv[j,],Xv[i,],gamma)
		sum = sum + temp
	}
	sum0 = sum0 + Yv[i] - sum
}
b = sum0 / length(alpha)

list(alpha=alpha, alpha_org = alpha_org, b=b, nSV=nSV, Xv=Xv, Yv=Yv, gamma=gamma)
}

### Predict the class of an object X


svmPolpredict <- function(x,X,Y,C0=Inf, gamma0=1.5,esp0=1e-10){
model = svmPoltrain(X,Y,C=C0, gamma=gamma0,esp=esp0)
alpha<-model$alpha
b<-model$b
Yv<-model$Yv
Xv<-model$Xv
nSV<-model$nSV
gamma<-model$gamma
sum = 0
for (i in 1 : length(alpha)){
	temp = alpha[i] * Yv[i] * pol_kernel(Xv[i,],x,gamma)
	sum = sum + temp
}
hat = sum + b
result<-sign(hat)
return(result)
}


gamma = 2
#---support vector:
res = svmPoltrain(X,Y,C=Inf, gamma = gamma)
res$alpha_org
res$nSV
#---prediction:
svmPolpredict(z, X, Y, gamma0 = gamma)



