######STA 6106-Final Exam-An Sun###
#####Problem 1#####










#####Problem 2#####
install.packages("mcsm")
library(mcsm)
data(challenger)
attach(challenger)
dim(challenger)
challenger

y<-challenger$oring;y
x0<-rep(1,23);x0
x1<-challenger$temp;x1
f3<-function(beta){
  
  grad1<- mean(y-1/(1+1/exp(beta[1]*x0+beta[2]*x1)))
  grad2<- mean(y*x1-x1/(1+1/exp(beta[1]*x0+beta[2]*x1)))
  gradient<-c(grad1,grad2)
  
  H11<- mean(-exp(beta[1]*x0+beta[2]*x1)/(1+exp(beta[1]*x0+beta[2]*x1))^2)
  H12<- mean(-x1*exp(beta[1]*x0+beta[2]*x1)/(1+exp(beta[1]*x0+beta[2]*x1))^2)
  H22<- mean(-x1^2*exp(beta[1]*x0+beta[2]*x1)/(1+exp(beta[1]*x0+beta[2]*x1))^2)
  Hessian<-matrix(c(H11,H12,H12,H22),2,2)
  return(list(gradient,Hessian))
}

#Part a solution
#IRLS
z<-cbind(x0,x1);head(z)
dim(z) #n*2 matrix
n<-dim(z)[1]
y<-matrix(y,n,1) #n*1 matrix
dim(y)

library(boot)
ftn7<-function(beta){ 
  
  pi<-inv.logit(z%*%beta) #inverse logit function to get the pi
  
  # create the diagonal matrix w
  w<-diag(n)
  for (i in 1:n){
    w[i,i]<-pi[i]*(1-pi[i])
  }
  pi<-matrix(pi,n,1)	#pi is n*1 matrix
  grad<-t(z)%*%(y-pi)  #gradient is z*1 matrix 
  Hessian<--t(z)%*%w%*%z #Hessian is z*z matrix
  return(list(grad,Hessian))
}

newtonglm <- function(f3, beta0, tol = 1e-6, n.max = 100) {
  # Newton's method for optimization, starting at beta0
  beta0<-matrix(beta0,2,1)  
  beta.new <- beta0
  f3.x <- f3(beta.new)
  iter <- 0
  
  while ((max(abs(f3.x[[1]]))>tol) & (iter < n.max)) {
    beta.new <- beta.new - solve(f3.x[[2]])%*%f3.x[[1]]
    f3.x <- f3(beta.new)
    iter <- iter + 1
  }
  if (iter == n.max) {
    cat('newton failed to converge\n')
  } else {
    return(list(iter, beta.new))
  }
}

beta0<-c(0,0)
cat(beta0,'-->',newtonglm(ftn7,beta0)[[2]],',number of iterations:', newtonglm(ftn7,beta0)[[1]], '\n')
####check####
glm( y ~ x1, family = binomial(logit))


#Part b
beta<-c(15.0429,-.2321627)
x<-c(1,31)
pi<-1/(1+exp(-sum(beta*x)));pi
e.x<-23*pi;e.x












#####Problem 3#######
###Part 1.a####
x <- 1                       # Initialize x
y <- 1                       # Initialize y
x.vector <- c(x,y)           # take vector of x,y
epsilon <- .00001            # testing min
stepsize <- .01              # stepsize tolerance 
for (i in 1:200) {
  print("Current step:")
  print(i)
  g <- 4*x*y + (x+y^2)^2                                      # Function g to be minimized
  g.grad <- c(2*(x + y*(y +2)), 4*(x*y + x + y^3))            # gradient of g
  g.hessian.x <- c(2, 4*(y+1))
  g.hessian.y <- c( 4*(y+1), 4*(x + 3*y^2))
  g.hessian <- matrix(c(t(g.hessian.x), t(g.hessian.y)), nrow=2, ncol=2)   # Hessian of g
  
  if (g.grad[1]<epsilon & g.grad[2]<epsilon)                  # End  loop if gradient of g is within epsilon of 0 when evaluated at current x
  {
    print("Function is minimized at:")
    print(x.vector)
    print("Minimum of function:")
    print(g)
    break
  }
  else
  delta <- -(solve(g.hessian)%*%g.grad)
  x.vector.0 <- x.vector
  x.vector <- x.vector + delta
  g0 <- g
  x <- x.vector[1]
  y <- x.vector[2]
  g <- 4*x*y + (x+y^2)^2
  while (g> g0 + stepsize)
  { print("Delta not sufficiently small. Halving delta.")
    delta <- delta/2
    x.vector <- x.vector.0 + delta
    x <- x.vector[1]
    y <- x.vector[2]
    g <- 4*x*y + (x+y^2)^2
  }
  x <- x.vector[1]
  y <- x.vector[2]
  print("value of g at current x")
  print(g)
}

######Part 1.b#####
#Steepest descent method with golden section method
gsection <- function(ftn, x.l, x.r, x.m, tol = 1e-9) {
  # applies the golden-section algorithm to maximise ftn
  # we assume that ftn is a function of a single variable
  # and that x.l < x.m < x.r and ftn(x.l), ftn(x.r) <= ftn(x.m)
  # the algorithm iteratively refines x.l, x.r, and x.m and terminates
  # when x.r - x.l <= tol, then returns x.m
  
  # golden ratio plus one
  gr1 <- 1 + (1 + sqrt(5))/2
  
  # successively refine x.l, x.r, and x.m
  f.l <- ftn(x.l)
  f.r <- ftn(x.r)
  f.m <- ftn(x.m)
  while ((x.r - x.l) > tol) {
    if ((x.r - x.m) > (x.m - x.l)) {
      y <- x.m + (x.r - x.m)/gr1
      f.y <- ftn(y)
      if (f.y >= f.m) {
        x.l <- x.m
        f.l <- f.m
        x.m <- y
        f.m <- f.y
      } else {
        x.r <- y
        f.r <- f.y
      }
    } else {
      y <- x.m - (x.m - x.l)/gr1
      f.y <- ftn(y)
      if (f.y >= f.m) {
        x.r <- x.m
        f.r <- f.m
        x.m <- y
        f.m <- f.y
      } else {
        x.l <- y
        f.l <- f.y
      }
    }
  }
  return(x.m)
}

#line search
line.search <- function(f, x, y, tol = 1e-9, a.max = 2^5) {
  # f is a real function that takes a vector of length d
  # x and y are vectors of length d
  # line.search uses gsection to find a >= 0 such that
  #   g(a) = f(x + a*y) has a local maximum at a,
  #   within a tolerance of tol
  # if no local max is found then we use 0 or a.max for a
  # the value returned is x + a*y
  
  if (sum(abs(y)) == 0) return(x) # g(a) constant
  
  g <- function(a) return(f(x + a*y))
  
  # find a triple a.l < a.m < a.r such that
  # g(a.l) <= g(a.m) and g(a.m) >= g(a.r)
  # a.l
  a.l <- 0
  g.l <- g(a.l)
  # a.m
  a.m <- 1
  g.m <- g(a.m)
  while ((g.m < g.l) & (a.m > tol)) {
    a.m <- a.m/2
    g.m <- g(a.m)
  }
  # if a suitable a.m was not found then use 0 for a
  if ((a.m <= tol) & (g.m < g.l)) return(x)
  # a.r
  a.r <- 2*a.m
  g.r <- g(a.r)
  while ((g.m < g.r) & (a.r < a.max)) {
    a.m <- a.r
    g.m <- g.r
    a.r <- 2*a.m
    g.r <- g(a.r)
  }
  # if a suitable a.r was not found then use a.max for a
  if ((a.r >= a.max) & (g.m < g.r)) return(x + a.max*y)
  
  # apply golden-section algorithm to g to find a
  a <- gsection(g, a.l, a.r, a.m)
  return(x + a*y)
}


ascent <- function(f, grad.f, x0, tol = 1e-9, n.max = 100) {
  # steepest ascent algorithm
  # find a local max of f starting at x0
  # function grad.f is the gradient of f
  
  x <- x0
  x.old <- x
  x <- line.search(f, x, grad.f(x))
  n <- 1
  while ((f(x) - f(x.old) > tol) & (n < n.max)) {
    x.old <- x
    x <- line.search(f, x, grad.f(x))
    n <- n + 1
  }
  return(x)
}

f<-function(x){
  f<- -1*(4*x[1]*x[2]+(x[1]+x[2]^2)^2)
  return (f)
}

grad.f<-function(x){
  grad1<- -1*(4*x[2]+2*(x[1]+x[2]^2))
  grad2<- -1*(4*x[1]+4*x[2]*(x[1]+x[2]^2))
  gradient<-c(grad1,grad2)
  return (gradient)
}

x0<-c(1,1)
result<-ascent(f,grad.f,x0)
result
g<--1*f(result)
g

###Part 2.a####
##Drawing a simple graph to see where the max value exists
x<-seq(.1,10,0.01)
y<-log(x)/(1+x)
plot(x,y,type = "l")
#Part a Newton's method
ftn3<-function(x){
    fx<-log(x)/(1+x)
    dfx1<-(x-x*log(x)+1)/(x*(1+x)^2)
    dfx2<-(-3*x^2+2*x^2*log(x)-4*x-1)/(x^2*(1+x)^3)
   return(c(fx,dfx1,dfx2))
}
Newton <- function(f3, x0, tol = 1e-6, n.max = 100) {
##Newton's method for optimization, starting at x0
##f3 is a function that given x returns the vector
##(g(x), g'(x), g''(x)), for some g
     x <- x0
     f3.x <- f3(x)
     n <- 0
     while ((abs(f3.x[2]) > tol) & (n < n.max)) {
       x <- x - f3.x[2]/f3.x[3]
      f3.x <- f3(x)
       n <- n + 1
       }
    if (n == n.max) {
      cat('newton failed to converge\n')
       } else {
         return(c(n,x))
         }
}
result<-Newton(ftn3,3.0,1e-6);result
ftn3(result[2])[1]

###Part 2.b####
#Secant method
ftn4<-function(x){
    fx<-log(x[2])/(1+x[2])
    dfx1<-numeric(2)
    dfx1[1]<-(x[1]-x[1]*log(x[1])+1)/(x[1]*(1+x[1])^2)
    dfx<-dfx1[2]<-(x[2]-x[2]*log(x[2])+1)/(x[2]*(1+x[2])^2)
    pseudo.dfx2<-(dfx1[2]-dfx1[1])/(x[2]-x[1])
    return(c(fx,dfx,pseudo.dfx2))
}
Secant <- function(f3, x0, tol = 1e-6, n.max = 100) {
# Secant method for optimization, starting at x0
# f3 is a function that given x returns the vector
# (f(x), f'(x), pseudo.f''((x)), for some f  
    x <- x0
    f3.x <- f3(x)
    n <- 0
    while ((abs(f3.x[2]) > tol) & (n < n.max)) {
        temp<-x[2]
        x[2] <- x[2] - f3.x[2]/f3.x[3]
        x[1] <-temp
        f3.x <- f3(x)
        n <- n + 1
        }
    if (n == n.max) {
        cat('secant failed to converge\n')
      } else {
    return(c(n,x))
        }
}
result<-Secant(ftn4,c(3,4),1e-6);result
ftn4(result[2:3])[1]


####Problem 4#####













  