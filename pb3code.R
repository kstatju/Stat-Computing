require('quadprog')
## Defining the Matrix Gaussian kernel
rbf_kernel <- function(x1,x2,gamma){
x1 = as.matrix(x1)
x2 = as.matrix(x2)
return(exp(-(1/gamma^2)*sum(diag(t(x1 - x2) %*% (x1 - x2)))))
}

svmtrain <- function(X,Y,C=Inf, gamma=1.5,esp=1e-10){
#matirx form, X is a list containing all the explanatory matrices
N<-length(Y)
Dm<-matrix(0,N,N)
Y<-as.vector(Y)

for(i in 1:N){
for(j in 1:N){
Dm[i,j]<-Y[i]*Y[j]*rbf_kernel(X[[i]],X[[j]],gamma)
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
		temp = alpha[j] * Yv[j] * rbf_kernel(Xv[[j]],Xv[[i]],gamma)
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
	temp = alpha_org[i] * Y[i] * rbf_kernel(X[[i]],x,gamma0)
	sum = sum + temp
}
hat = sum + b
result<-sign(hat)
return(result)
}

