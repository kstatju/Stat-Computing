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