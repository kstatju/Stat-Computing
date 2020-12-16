# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 21:08:43 2016

@author: Kanak
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pyplot as pltbar
from mpl_toolkits.mplot3d import Axes3D
from math import exp
from scipy.stats import chi2
from scipy.stats import binom
from scipy.special import gamma


########################################################################
#
#                Problem 1
#
########################################################################

burn = 500
a = 5
b = 10
n = 25

def ncr(n, r):
    r = min(r, n-r)
    if r == 0: return 1
    numer = 1
    for i in range(n, n-r, -1):
        numer *= i
    denom = 1
    for i in range(1, r+1):
        denom *= i
    return numer//denom
    
    
def prob_func(data, a, b, n):
    length = len(x)
    for i in range(0, length):
        data[i].append(ncr(n, data[i][0])*(data[i][1]**(data[i][0]+a-1)
        )*((1-data[i][1])**(n-data[i][0]+b-1)))
        
    return(data)    
    
def gib_sample(a, b, n, sample_size = 1, burn = 500):
    chain = sample_size + burn
    x = np.random.binomial(n, 0.5, size=None)
    y = np.random.beta(x + a, n - x +b, size=None)
    sam = []
    sam.append([x, y])
    
    for _ in range(1, chain):
        x = np.random.binomial(n, y, size=None)
        y = np.random.beta(x+a, n-x+b, size=None)
        sam.append([x, y])
    return(sam[burn:chain])
    
    
def marginalpdf(x, n, a, b):
    fx = ncr(n, x)*gamma(a+b)*gamma(x+a)*gamma(n-x+b)/(gamma(a)*gamma(b)*gamma(a+b+n))
    return fx   

x = gib_sample(sample_size = 10000, a = a, b = b, n = n)
    
xp = pd.DataFrame(prob_func(x, a, b, n))
xp.columns = ["x", "y", "p"]

fig = plt.figure(figsize=(12,10))
ax = fig.add_subplot(111, projection='3d')

z = xp["p"]
x = xp["x"]
y = xp["y"]

ax.scatter(x, y, z, c='r', marker='o')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Probability')

plt.show()
plt.clf()
plt.cla()
plt.close()

xran = np.array(range(0, n))
bp = []
for i in xran:
    bp.append(np.mean(binom.pmf(i, n, xp["y"])))

fig, axs=pltbar.subplots(1,2, figsize=(10, 6), sharex='col', sharey='row')
width = 1
axs[0].bar(xran, bp, width, color="blue")

bpm = []
for i in xran:
    bpm.append(marginalpdf(i, n, a, b))
    
axs[1].bar(xran, bpm, width, color="blue")
pltbar.show()
pltbar.clf()
pltbar.cla()
pltbar.close()

meanbp = np.sum(xran*bp)
varbp = np.sum(xran**2*bp) - meanbp**2

meanbpm = np.sum(xran*bpm)
varbpm = np.sum(xran**2*bpm) - meanbpm**2

print("Simulated Mean: ", meanbp, ", Simulated Var: ", varbp, 
      ", Marginal Dist Mean: ", meanbpm, ", Marginal Dist Var: ", varbpm, )

########################################################################
#
#                Problem 2
#
########################################################################
def rlpdf(x, sigma):
    if x < 0 or sigma <= 0: return
    return((x/sigma**2)*exp(-x**2/(2*sigma**2)))
    
def rand_rayleigh(n, sigma, x1 = None):
    if not x1 == None:
        x = []
        x.append(x1)
    else:
        x = []
        x.append(np.random.chisquare(1))
    k = 0
    u = np.random.uniform(size = n)
    
    for i in range(1, n):
        xt = x[i-1]
        y = np.random.chisquare(xt)
        num = rlpdf(y, sigma)*chi2.pdf(xt, y)
        den = rlpdf(xt, sigma)*chi2.pdf(y, xt)
        if (u[i] <= num / den):
            x.append(y)
        else:
            x.append(xt)
            k += 1
    return (np.array(x))
    
def gelman_rubin(x):
    if not isinstance(x, np.ndarray):
        try:
            x = np.array(x)
        except:
            print ("x is not possible to convert as array")
            return    
    n = x.shape[1]
    x_means = np.mean(x, axis = 1)
    b = n * np.var(x_means)
    x_w = np.var(x, axis = 1)
    W = np.mean(x_w)
    v_hat = W* (n-1) / n + (b/n)
    r_hat = np.sqrt(v_hat / W)
    return(r_hat)
        


n = 15000
b = 500
sigma = 20
x0 = [0.5, 5, 10, 20, 50, 100, 300, 500]
indexlow = 0
indexup = 2000

x = rand_rayleigh(n, sigma)

plt.plot(range(indexlow, indexup), x[indexlow:indexup])


  
  
x = []
for i in range(0, len(x0)):
    x.append(rand_rayleigh(n, sigma, x0[i]))
        
x = np.array(x)

cumsumx = np.cumsum(x, axis = 1)
ncol = np.array(range(1, x.shape[1]+1))
cumsumx = cumsumx / ncol

   
fig ,axs=plt.subplots(len(x0),1, figsize=(10, 50))
for i in range(0, len(x0)):
    axs[i].plot(range(0, n-b), cumsumx[i, b:n])
    
plt.show()
plt.clf()
plt.cla()
plt.close()
rhat = []
xrange = np.array(range(b, n))
for i in xrange:
    rhat.append(gelman_rubin(cumsumx[:,0:i]))
    
rhat = np.array(rhat)    
h = 1.2
lst = min(xrange[rhat<1.2])

plt.plot(xrange,rhat)
plt.axvline(x=lst, c="red", linewidth=0.5,zorder=0)
plt.axhline(y=h, c="red", linewidth=0.5,zorder=0)
plt.show()

print("Burning numbers", lst)




n = 15000
b = 10
sigma = [10, 20, 30]
x0 = [5, 10, 20, 50, 100]


listn = []   
for j in range(0, len(sigma)):
       
    x = []
    for i in range(0, len(x0)):
        x.append(rand_rayleigh(n, sigma[j], x0[i]))
            
    x = np.array(x)
    
    cumsumx = np.cumsum(x, axis = 1)
    ncol = np.array(range(1, x.shape[1]+1))
    cumsumx = cumsumx / ncol

    rhat = []
    xrange = np.array(range(b, n))
    for i in xrange:
        rhat.append(gelman_rubin(cumsumx[:,0:i]))
        
    rhat = np.array(rhat)    
    h = 1.2
    lst = min(xrange[rhat<1.2])
    listn.append(lst)
    
    plt.plot(xrange,rhat)
    plt.axvline(x=lst, c="red", linewidth=0.5,zorder=0)
    plt.axhline(y=h, c="red", linewidth=0.5,zorder=0)
    plt.show()
    
    print("Burning numbers", lst)
    
print(listn)