####Problem 1
#compute in the Gelman.Rubin function ###
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi)      #row means
  B <- n * var(psi.means)         #between variance est.
  psi.w <- apply(psi, 1, "var")   #within variances
  W <- mean(psi.w)                #within est.
  v.hat <- W*(n-1)/n + (B/n)      #upper variance est.
  r.hat <- v.hat / W              #G-R statistic
  return(r.hat)
}

# evaluate the Rayleigh density
f<- function(x,sigma){
  if(any(x<0)) return(0)
  stopifnot(sigma>0)
  return((x/sigma^2)*exp(-x^2/(2*sigma^2)))
}

# written as function for M-H sampler

normal.chain=function(x0,n)
{
  set.seed(1234)
  r=rep(0,n)
  x=x0
  for(k in 1:n)
  {
    u=runif(1)
    y=rnorm(1,x,.1)
    if(u<alpha(x,y))
    {
      x=y
    }
    else
    {
      x=x0
    }
    r[k]=x
  }
  return(r)
}
p1=0.7
alpha = function(x,y) {
  # Acceptance probability calculation
  return( min( 1, (p1*dnorm(y,7,.5)+(1-p1)*dnorm(y,10,.5))/(p1*dnorm(x,7,.5)+(1-p1)*dnorm(x,10,.5))))
}

n <- 10000    #length of chains
sigma <- .2 #parameter of proposal distribution
k <- 3 #number of chains to generate
b <- 1000 #burn-in length

#choose initial values as 1,2,3,4
x0 <- c(0, 7, 15)

#generate the chains
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k){
  X[i, ] <- normal.chain(x0[i],sigma)
}

#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi[,1:n]))

#plot psi for the four chains
par(mfrow=c(2,2))
for (i in 1:k)
  plot(psi[i, (b+1):n], type="l",
       xlab=i, ylab=bquote(psi))
par(mfrow=c(1,1)) #restore default

#plot the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)




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
