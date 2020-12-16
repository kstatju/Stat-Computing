#############################################
#### 1a
x <- c(28,-44,29,30,26,27,22,23,33,16,24,40,21,31,34,-2,25,19)
N=length(x)
B=1000
theta=mean(x)
sigma2=var(x)
alpha = 0.05
## statistic T

hi = list()
for (i in 1:B){
  y = sample(x,N,replace = T)
  t=(y-mean(y))^2/var(y)
  hi[i]= unlist(quantile(t, 1-alpha))
}
h =mean(unlist(hi))
h
if((38-theta)^2/sigma2 > h){
  print("Outlier")}else {print("Not Outlier")}

  