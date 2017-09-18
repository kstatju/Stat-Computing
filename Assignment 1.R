library(distr)

#######################################################################
#  Problem 1
# Function to generate random number from a discrete distribution. The function has four parameters. n = sample size, xrange  =  range of random variables, prob = probability for the random variable, seed = Random seed

rDiscrete <- function(n=1, xrange = c(0,1), prob = c(0.5,0.5), seed = NULL){
  if(!is.null(seed))  set.seed(seed)  # set seed for the random number
  
  # checking length of xrange and prob and length should be equal
  
  if (length(xrange) != length(prob)){
    return(message("The dimenstion of both xrange and 
                   prob should be equal"))
  }
  
  # checking for probability range. It should be between 0 and 1.
  
  if ((any(prob >1) || any(prob < 0)) ){
    return(message("Probability should be between 0 and 1"))
  }
  
  # checking sum of probability. Sum of probability should be 1.
  
  if(sum(prob) != 1){
    return(message("Sum of Probability must be 1"))
  }
  
  prob = cumsum(prob)		# Cumulative probability 
  
  uni = runif(n, min = 0, max = 1)	# uniform random number of size n 
  rnd = NULL				
  
  # assigning random number for each xrange to rand
  
  for (i in 1:length(xrange)){
    if (i != 1){
      # if uni lies between a & b choose x
      rnd[which(uni > prob[i-1] & uni < prob[i])] <- xrange[i]                
    }else {
      rnd[which(uni > 0 & uni < prob[i])] <- xrange[i]
    }
  }
  
  rm(uni)
  return(rnd)  
  
  }


# calling function rDiscrete to generate random number form a probability distribution.

random_num = rDiscrete(n = 100, xrange = c(0,1,3), 
                       prob = c(.3,.2,.5), seed = 2321)
random_num


# Probability for random variable calculated from the sample.

table(random_num)/length(random_num)

###################################################################################



#############################################################################
# Problem 2
# Part 2 Box-Muller Algorithm for Normal random variate. Function rnormalBX generates 2 sets of random normal variate using Box-Muller Algorithm. The function contains four parameter, n = sample size, mean = mean of normal distribution, var = variance of normal distribution, seed = random seed.

rnormalBX <- function(n=1, mean = 0, var = 1, seed = NULL){
  
  if(!is.null(seed))  set.seed(seed)  # set seed for the random number
  
  # checking for n whether it is numeric or greater than 1.
  
  if (n < 1 || !is.numeric(n)){
    return(message("n must be greater than or equal to 1 and 
                   numeric value"))
  }
  
  # checking for var whether it is less than equal 0 or numeric and length greater than 1.
  
  if (var <= 0 || !is.numeric(var) || length(var) > 1){
    return(message("var must be greater than or equal to 1 and
                   length 1"))
  }
  
  # checking for mean whether it is numeric and length greater than 1.
  
  if (length(mean) >1 || !is.numeric(mean)){
    return(message("mean must be length 1 and numeric value"))
  }
  
  # generate two uniform and then calculate r and theta and finally generate normal random number
  
  u = runif(n, 0, 1)
  v = runif(n, 0, 1)
  r = sqrt(-2*log10(v))
  theta = 2*pi*u
  norm = cbind(r*cos(theta), r*sin(theta))
  norm = (norm - colMeans(norm))/sqrt(diag(var(norm)))
  rm(u,v, r, theta)
  
  # transform normal random number with specified mean and variance 
  
  if (mean != 0 || var != 1){
    norm = sqrt(var)*norm + mean
    
  }
  return(unname(norm, force = TRUE))
  }

# calling rnormalBX function to generate normal random number using Box-Muller algorithm
x = rnormalBX(n = 1000, mean = 52, var = 10, seed = 123)
summary(x)


#Part 1 Independence of two normal random variable generated using Box-Muller algorithm

# Density Curve for the normal sample. Function densityCurve has two parameters, data = data set for the density curve, dist = distribution.

densityCurve <- function(data , dist = "Normal"){
  
  # checking and converting data into data frame.
  
  if (!is.data.frame(data)){  data = as.data.frame(data)  }
  
  # dividing graphical outlet based on dimension of data
  
  par(mfrow = c(dim(data)[2],1))
  
  # creating histogram with normal density curve for each variable
  
  for (i in 1:dim(data)[2]){
    his  <- hist(data[,i], col="red", xlab="x", main="Histogram with Normal Curve") 
    xfit <- seq(min(data[,i]),max(data[,i]),length=nrow(data)) 
    yfit <- dnorm(xfit, mean = mean(data[,i]), sd = sd(data[,i])) 
    yfit <- yfit*diff(his$mids[1:2])*length(data[,i]) 
    lines(xfit, yfit, col="blue", lwd=2)
  }
  par(mfrow = c(1,1))
}

# Calling densityCurve function to create density curve.
densityCurve(x)



# Independence Test for the normal random variables. Function indepTest has one parameter, data = data set. The function will output correlation test, QQ Plot, Kolmogorov-Smirnov test, Wilcoxon Signed Rank Test for checking whether the sample comes from normal and they are independent. 

indepTest <- function(data){
  
  #checking for data dimension
  
  if(dim(data)[2] >2){
    return(message("Data contains more than two variables"))
  }
  
  #correlation test independence 
  
  cortest = cor.test(data[,1], data[,2])
  
  # QQ-plot and KS test for each variable for normality check
  
  par(mfrow = c(dim(data)[2],1))
  kstest = NULL
  
  for (i in 1:dim(data)[2]){
    qqnorm(data[,i])
    kstest[[i]] = ks.test(data[,i], "pnorm")
  }
  
  #Wilcoxon Signed Rank Test for Normal population
  
  wilcos = wilcox.test(data)
  
  # creating output report
  
  if (cortest$p.value > 0.5 & kstest[[1]][2] < 0.05 & 
      kstest[[2]][2] < 0.05 & wilcos$p.value < 0.05){
    comm = "Correlation test indicates two samples are independent, Kolmogorov-Smirnov \n and Wilcoxon Signed Rank Test indicted that the two sample come from Normal distribution"
  } else {
    
    comm = paste("p-value for all tests-","\n Correlation Test", cortest$p.value[1], "\n Kolmogorov-Smirnov Test", kstest[[1]][2][[1]],
                 'and', kstest[[2]][2][[1]], "\n Wilcoxon Signed Rank Test",
                 wilcos$p.value[1])
  }
  message(comm)
  return(c(CorTest = cortest, KStest = kstest, WilCos = wilcos, 
           Comment = comm))
}

#calling function for Independence and Normality check 

indepTest(x)



###################################################################################



###################################################################################
# Problem 3
####################################################################################
# Problem 3 



# Importance Sampling from Poisson distribution
# function rpoinsImportance contains five parameters, n = sample size, lambda = Poission parameter of the target distribution,
# imp.func.pram = Parameter value of support distribution, m = simulation times for getting average of the importance sampling function
# seed = random seed.

rpoinsImportance <- function(n = 1, lambda = 1, imp.func.pram = 1, m = 100, seed = NULL){
  
  if(!is.null(seed))  set.seed(seed)  # set seed for the random number
  
  # checking for n whether it is numeric or greater than 1.
  if (n < 1 || !is.numeric(n)){
    return(message("n must be greater than or equal to 1 and numeric value"))
  }
  
  # replicate rpoinsImportance function n times to generate n size random sample.
  
  if (n != 1){
    rnd = replicate(n, rpoinsImportance(lambda = lambda, imp.func.pram = imp.func.pram, m = m))
    return(rnd)
  }
  
  # generate random number of size m from Poission distribution
  
  fx = rpois(m, lambda = imp.func.pram)
  
  # generate random sample from our required Poisson distribution using importance sampling method
  
  rnd = round(mean(fx*dpois(fx, lambda = lambda)/dpois(fx, lambda = imp.func.pram)))
  
  return(rnd)
  
}



x = rpoinsImportance(n = 10000, lambda = 4.5, imp.func.pram = 2.4653)
tab = table(x)
barplot(tab)

mean(x)
var(x)



# Function for calculating Likelyhood Ratio Test. function likelytest contains 

likelytest <-  function(sample = NULL, dist.with.pram = 1, null.hypo = 1, sample_size = 1000, 
                        sim_p.value = FALSE, replicat.times = 10000, imp.sampling = FALSE, 
                        imp.func.pram = 1, seed = NULL, CI = TRUE){
  if(!is.null(seed))  set.seed(seed)  # set seed for the random number
  
  if (sim_p.value){
    if (!is.null(sample)) dist.with.pram = mean(sample)
    
    likelyhoodlist = replicate(replicat.times, likelytest(sample = NULL, dist.with.pram = null.hypo, 
                                                          null.hypo = null.hypo, sample_size = sample_size, 
                                                          imp.sampling = imp.sampling, imp.func.pram = imp.func.pram,
                                                          CI = FALSE))
  }
  
  if (is.null(sample)){
    if (!imp.sampling) { rand = rpois; pand = dpois    
    rand_var = rand(n = sample_size, lambda = dist.with.pram)
    } else if (imp.sampling) { rand = rpoinsImportance; pand = dpois    
    rand_var = rand(n = sample_size, lambda = dist.with.pram, imp.func.pram = imp.func.pram, m = 100, seed = seed)
    } 
  } else rand_var = sample
  
  if (CI & !imp.sampling){
    lambda.dist = replicate(replicat.times, mean(rand(n = sample_size, lambda = dist.with.pram)))
  } else if (CI & imp.sampling){ 
    lambda.dist = replicate(replicat.times, mean(rand(n = sample_size, lambda = dist.with.pram, 
                                                      imp.func.pram = imp.func.pram, m = 100, seed = seed)))
  }
  alter_hypo = mean(rand_var)
  
  likelyhood = prod(pand(x = rand_var, lambda = null.hypo))/prod(pand(rand_var, lambda = alter_hypo))
  
  
  if (sim_p.value) {
    
    ucl = quantile(lambda.dist, probs = c(0.025, .975))
    return(c(lambda.hat = alter_hypo, CI = c(LCL = ucl[1], UCL = ucl[2]), 
             test.statistics = likelyhood,  
             p.value = 1 - sum(likelyhoodlist["test.statistics",] > likelyhood)/length(likelyhoodlist["test.statistics",])))
    #return(likelyhoodlist)
    
  }else {
    return(c(lambda.hat = alter_hypo, test.statistics = likelyhood))
  }
}

kk = likelytest(dist.with.pram = 2.5, null.hypo = 2, replicat.times = 1000, sample_size = 10, 
                sim_p.value = TRUE)
kk

kk = likelytest(dist.with.pram = 2.5, null.hypo = 2, replicat.times = 1000, sample_size = 10, 
                sim_p.value = TRUE, imp.sampling = TRUE, imp.func.pram = 2.46)
kk






###################################################################################
# Problem 4
# Normal Random number generation using Accept-Reject Algorithm. Function rnormAR contains seven parameters, n = sample size, M.constant = value of M, dist = Distribution , tar.parameter = parameter of normal distribution, support.dist = Support Distribution Name, sub.parameter = Parameters of support distribution, seed = Random Seed. Function will output random number list.

rnormAR <-  function(n = 1, M.constant = 1, dist = "Normal", 
                     tar.parameter = c(0,1), support.dist = 
                       "Double Exponential", sub.parameter = NULL, 
                     seed =NULL){
  
  if(!is.null(seed))  set.seed(seed)  # set seed for the random number

  
  # checking for n whether it is numeric or greater than 1.
  
  if (n < 1 || !is.numeric(n)){
    return(message("n must be greater than or equal to 1 and numeric value"))
  }
  
  # checking for M whether it is numeric or lest thank 1
  
  if (M.constant < 1 || !is.numeric(n)){
    return(message("M.constant must be greater than or equal to 1 and numeric value"))
  }
  
  if (dist  == "Normal"){
    
    # checking for var whether it is less than equal 0 or numeric and length greater than 1.
    
    if (tar.parameter[2] < 1 || !is.numeric(tar.parameter[2])){
      return(message("var must be greater than or equal to 1 and
                     length 1"))
    }
    
    # checking for mean whether it is numeric and length greater than 1.
    
    if (length(tar.parameter[1]) >1 || !is.numeric(tar.parameter[1])){
      return(message("mean must be length 1 and numeric value"))
    }
    
    # setting parameter for Double Exponential distribution
    
    if (is.null(sub.parameter)){
      de = DExp(1)
    }else {
      de = DExp(sub.parameter)
    }
    sam.size = 0
    sample = NULL
    
    # Generating sample from normal distribution using f(x)/(M*g(x)) of size n
    
    while (sam.size < n){
      uni = runif((ceiling((n-sam.size)*1.5)), 0,1)
      rde = r(de)(ceiling((n-sam.size)*1.5))    
      gx = dnorm(rde)/(M.constant*d(de)(rde))
      sample = rbind(sample, rde[which(uni < gx)])
      sam.size = length(sample)
    }
    rm(uni, rde, gx, sam.size, de)
    sample = sample[1:n]
    if (tar.parameter[1] != 0 || tar.parameter[2] != 1){
      sample = sqrt(tar.parameter[2])*sample + tar.parameter[1]
      
    }
    }
  return(sample)
}

# calling function and checking density curve.

x = rnormAR(10000, M.constant = sqrt(2/pi)*exp(1/2), tar.parameter = c(0,1))
densityCurve(x)



###################################################################################


# Problem 5


f3 <- function(x1, x2){
  if (((x1-1)^2+(x2-1)^2 <= 1) & ((x1-1)^2+(x2+1)^2 <= 1)){
    return(x1^2+x2^2)
    
  }else {return(NA)}
}

f1 <- function(x1, x2){
  if ((x1-1)^2+(x2-1)^2 <= 1){
    return(x1^2+x2^2)
    
  }else {return(NA)}
}

f2 <- function(x1, x2){
  if ( (x1-1)^2+(x2+1)^2 <= 1){
    return(x1^2+x2^2)
    
  }else {return(NA)}
}


f <- function(x1, x2){
  return(x1^2+x2^2)
  
}

x1r = seq(-2, 2, by = 0.1)
x2r = seq(-2, 2, by = 0.1)


for (i in 1:length(x1r)){
  y = NULL
  for (j in 1:length(x2r)){
    x = f(x1r[i], x2r[j])
    y = rbind(y, x)
  }
  if (i == 1) {z = y; dim(z)
  } else {dim (z); z = cbind(z, y)}
}

for (i in 1:length(x1r)){
  y = NULL
  for (j in 1:length(x2r)){
    x = f1(x1r[i], x2r[j])
    y = rbind(y, x)
  }
  if (i == 1) {z1 = y
  } else {z1 = cbind(z1, y)}
}
for (i in 1:length(x1r)){
  y = NULL
  for (j in 1:length(x2r)){
    x = f2(x1r[i], x2r[j])
    y = rbind(y, x)
  }
  if (i == 1) {z2 = y;
  } else { z2 = cbind(z2, y)}
}

for (i in 1:length(x1r)){
  y = NULL
  for (j in 1:length(x2r)){
    x = f3(x1r[i], x2r[j])
    y = rbind(y, x)
  }
  if (i == 1) {z3 = y;
  } else { z3 = cbind(z3, y)}
}

contour(x2r, x1r, as.matrix(z))
contour(x2r, x1r, as.matrix(z1), add = TRUE, col = 2)
contour(x2r, x1r, as.matrix(z2), add = TRUE, col = 3)
contour(x2r, x1r, as.matrix(z3), add = TRUE, col = 4)




####################################################################################
# Problem 6 


# Random number for Exponential - Gamma Mixture density function
## Random number generation from Exponential-Gamma Mixture distribution. Function rExpGammaMix contains
## three parameter n = sample size, gamma.parm = parameter for gamma distribution, and seed = random seed.


rEGMix <- function(n = 1, gamma.pram = c(r = 1, b = 1), seed = NULL){
  
  if(!is.null(seed))  set.seed(seed)  # set seed for the random number
  
  # checking for n whether it is numeric or greater than 1.
  
  if (n < 1 || !is.numeric(n)){
    return(message("n must be greater than or equal to 1 and numeric value"))
  }
  
  # checking Gamma Parameter whether it is numeric or greater than 1.
  
  if (any(gamma.pram <= 0) || !is.numeric(gamma.pram) || length(gamma.pram) > 2){
    return(message("gamma.pram must be greater than 0 and length 2"))
  }
  
  # Generate Gamma random variable of size n with given parameters
  
  gma = rgamma(n, shape = gamma.pram[[1]], scale = gamma.pram[[2]])
  
  # Using Gamma random variable as the parameter for exponential distribution generate Exponential Random variable.
  
  return(rexp(n, rate = gma))
  
}

# Calling funtion to generate Exponential-Gamma Mixture random Variable.

x = rEGMix(10000, gamma.pram = c(4,2))
hist(x)



# alternative way to generate random sample from Exponential-Gamma Mixture Distribution.

# Pareto probability density function


dpereto <- function(x, shape = 1, scale = 1, func = "Lomex", seed = NULL){
  
  if(!is.null(seed))  set.seed(seed)  # set seed for the random number
  
  if (n < 1 || !is.numeric(n)){
    return(message("n must be greater than or equal to 1 and numeric value"))
  }
  if (shape < 1 || !is.numeric(shape) || length(shape) > 1){
    return(message("shape parameter must be greater than 0 and length 1"))
  }
  if (scale < 1 || !is.numeric(scale) || length(scale) > 1){
    return(message("scale parameter must be greater than 0 and length 1"))
  }
  
  
  return(shape/scale*(1/(x/scale+1)^(shape+1)))
  
}


# Random number for Pareto probability density function

rpereto <- function(n = 1, shape = 1, scale = 1, seed = NULL){
  
  if(!is.null(seed))  set.seed(seed)  # set seed for the random number
  
  if (n < 1 || !is.numeric(n)){
    return(message("n must be greater than or equal to 1 and numeric value"))
  }
  if (shape <= 0 || !is.numeric(shape) || length(shape) > 1){
    return(message("shape parameter must be greater than 0 and length 1"))
  }
  if (scale <= 0 || !is.numeric(scale) || length(scale) > 1){
    return(message("scale parameter must be greater than 0 and length 1"))
  }
  
  uni = runif(n, min = 0, max = 1)
  return(scale*(1-uni)^(-1/shape)-scale)
  
}

x = rpereto(10000, 4,2)
y = dpereto(x, 3, 4)
hist(x)



