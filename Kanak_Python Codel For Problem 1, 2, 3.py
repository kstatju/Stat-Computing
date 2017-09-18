# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 21:08:43 2016

@author: Kanak
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm



########################################################################
#                Problem 1
# Question a
########################################################################
def norm_mix(data, mu, sigma, p):
    a = p*norm(mu[0], sigma[0]).pdf(data)+(1-p)*norm(mu[1], sigma[1]).pdf(data)     
    return(a)
    
def rand_norm_mix(n, mu, sigma, p, pro_par, x1 = None):
    if not x1 == None:
        x = []
        x.append(x1)
    else:
        x = []
        x.append(norm(0, pro_par).rvs(1)[0])
    u = np.random.uniform(size = n)
    
    for i in range(1, n):
        xt = x[i-1]
        y = norm(xt, pro_par).rvs(1)[0]
        num = norm_mix(y, mu, sigma, p)*norm(y, pro_par).pdf(xt)
        den = norm_mix(xt, mu, sigma, p)*norm(xt, pro_par).pdf(y)
        if (u[i] <= num / den):
            x.append(y)
        else:
            x.append(xt)
    return (np.array(x))
    
    
n = 10000
mu = [7, 10]
sigma = [.5, .5]
p = 0.7
pro_par = 0.01
x0 = [0,7,15]
b=100
indexlow = b
indexup = n

x = []
for i in range(0, len(x0)):
    x.append(rand_norm_mix(n, mu, sigma, p, pro_par, x0[i]))
        
x = np.array(x)
#fig ,axs=plt.subplots(len(x0),1, figsize=(10, 50))
for i in range(0,len(x0)):
    plt.plot(range(indexlow, indexup), x[i,indexlow:indexup])
    plt.show()
    plt.clf()
    plt.cla()
    plt.close()



#fig ,axs=plt.subplots(len(x0),1, figsize=(10, 50))
for i in range(0,len(x0)):
    prob = norm_mix(x[i,indexlow:indexup], mu, sigma, p)
    myHist = plt.hist(x[i, indexlow:indexup], 40, normed=True)
    xx = np.arange(min(x[i, indexlow:indexup])-2, max(x[i,indexlow:indexup])+2,0.001)
    prob = norm_mix(xx, mu, sigma, p)
    h = plt.plot(xx, prob, lw=2)
    plt.show()
    plt.clf()
    plt.cla()
    plt.close()



########################################################################
#                Problem 1
# Question b
########################################################################


n = 10000
mu = [7, 10]
sigma = [.5, .5]
p = 0.7
pro_par = 1
x0 = [0,7,15]
b=1000
indexlow = b
indexup = n

x = []
for i in range(0, len(x0)):
    x.append(rand_norm_mix(n, mu, sigma, p, pro_par, x0[i]))
        
x = np.array(x)
#fig ,axs=plt.subplots(len(x0),1, figsize=(10, 50))
for i in range(0,len(x0)):
    plt.plot(range(indexlow, indexup), x[i,indexlow:indexup])
    plt.show()
    plt.clf()
    plt.cla()
    plt.close()

#fig ,axs=plt.subplots(len(x0),1, figsize=(10, 50))
for i in range(0,len(x0)):
    prob = norm_mix(x[i,indexlow:indexup], mu, sigma, p)
    myHist = plt.hist(x[i, indexlow:indexup], 40, normed=True)
    xx = np.arange(min(x[i, indexlow:indexup])-2, max(x[i,indexlow:indexup])+2,0.001)
    prob = norm_mix(xx, mu, sigma, p)
    h = plt.plot(xx, prob, lw=2)
    plt.show()
    plt.clf()
    plt.cla()
    plt.close()

########################################################################
#                Problem 1
# Question c
########################################################################

 
  
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
    

n = 10000
mu = [7, 10]
sigma = [.5, .5]
p = 0.7
pro_par = 1
x0 = [0,7,15]
b=100
indexlow = b
indexup = n

x = []
for i in range(0, len(x0)):
    x.append(rand_norm_mix(n, mu, sigma, p, pro_par, x0[i]))
    
x = np.array(x)        
cumsumx = np.cumsum(x, axis = 1)
ncol = np.array(range(1, x.shape[1]+1))
cumsumx = cumsumx / ncol

   
#fig ,axs=plt.subplots(len(x0),1, figsize=(10, 50))
for i in range(0, len(x0)):
    plt.plot(range(0, n-b), cumsumx[i, b:n])
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
condi = xrange[rhat<1.2]
if (sum(condi) == 0):
    lst = n
else:
    lst = min(condi)
plt.plot(xrange,rhat)
plt.axvline(x=lst, c="red", linewidth=0.5,zorder=0)
plt.axhline(y=h, c="red", linewidth=0.5,zorder=0)
plt.show()
plt.clf()
plt.cla()
plt.close()


n = 10000
mu = [7, 10]
sigma = [.5, .5]
p = 0.7
pro_par = [0.1, .5, 1, 2, 5]
x0 = [0,7,15]
b=100
indexlow = b
indexup = n

listn = []   
for j in range(0, len(pro_par)):
       
    x = []
    for i in range(0, len(x0)):
        x.append(rand_norm_mix(n, mu, sigma, p, pro_par[j], x0[i]))
            
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
    condi = xrange[rhat<1.2]
    if (sum(condi) == 0):
        lst = n
    else:
        lst = min(condi)
    listn.append(lst)
    
    plt.plot(xrange,rhat)
    plt.axvline(x=lst, c="red", linewidth=0.5,zorder=0)
    plt.axhline(y=h, c="red", linewidth=0.5,zorder=0)
    plt.show()
    plt.clf()
    plt.cla()
    plt.close()
    
    print("Burning numbers", lst)
    
print(listn)






########################################################################
#                Problem 2
#
########################################################################




# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 20:17:09 2016

@author: ka746940
"""

import numpy as np
import pandas as pd
from numpy.linalg import inv
import matplotlib.pyplot as plt

dt = pd.read_csv("C:/Users/ka746940/Desktop/UCF/STA 6106 - Statistical Computing/Assignments/Final/Final_chao/challenge.csv");
dt = np.array(dt)

def probability(theta, x):
    prob = (1.0/(1.0+np.exp(-np.dot(x, theta))))
    return prob


#Get Log Likelihood of Dataset
def log_likelihood(Y,p):
    loglikelihood = Y*np.log(p+1e-24) + (1-Y)*np.log(1-p+1e-24)
    return -1*loglikelihood.sum()


def gradient_hessian(X,error, W):
    """Calculate gradient and Hessian for the current set of weights and features."""
    gred = (np.dot(error,X))
    hessi = np.dot(X.T,np.dot(W,X))
    return {"gradient": gred, "hessian": hessi}


def newton_raphson(Y, X, theta0 = None, tol = 1e-6, maxite = 200):
    if not isinstance(X, np.ndarray) or not isinstance(Y, np.ndarray):
        try:
            X = np.array(X)
            Y = np.array(Y)
        except:
            print ("X is not possible to convert as array")
            return
            
    ncol = X.shape[1]
    
    if (theta0 == None):
        theta = np.repeat(0, ncol)
    else:
        theta = theta0
    
    itter = 0
    conv = 999999999

    while (itter < maxite  and conv> tol):
        itter = itter + 1
        p = probability(theta, X)
        W = np.diag(p)
        error = Y - p
        gre_heis = gradient_hessian(X, error, W)
        update = np.dot(inv(gre_heis["hessian"]), gre_heis["gradient"])
        theta = theta + update
        conv = max(np.absolute(update))
    
    if (itter == maxite):
        print("Failed to converge")
        return
    else:
        p = probability(theta, X)
        likeli =   log_likelihood(Y, p)
        error = Y - p
        W = np.diag(p)
        gr_hes = gradient_hessian(X, error, W)
        return({"prem": theta, "gradient": gr_hes["gradient"], \
        "Hessian": gr_hes["hessian"], "error": error, "iteration": itter, \
        "likelihood": likeli, "probability": p})


def pred_prob(x, model):
    if not isinstance(x, np.ndarray):
        try:
            x = np.array(x)
        except:
            print ("x is not possible to convert as array")
            return  
    prem = model["prem"]
    pred = probability(prem,x)
    return({"pred": pred})


def feature_mat(dt):
    X = []    
    for i in range(len(dt)):
        X.append([1,dt[i,1]])
    Y = dt[:,0]
    return(np.array(X), np.array(Y))

X, Y = feature_mat(dt)

model = newton_raphson(Y, X, maxite = 2000)

pred = pred_prob([1, 31], model)

SSE = np.sum((Y - pred["pred"])**2)

print(model)
print(pred)
print(SSE)



plt.scatter(X[:,1], Y, label='Scatter Plot')
plt.plot(X[:,1], model["probability"], 'k', label='Fitted Line', lw = 2, color = 'r')
for i in range(len(Y)):
    if(i == 0):
        plt.plot([X[i,1], X[i,1]], [model["probability"][i], Y[i]], 'k', label = "Error", color = 'g', lw = 2)
    else:
        plt.plot([X[i,1], X[i,1]], [model["probability"][i], Y[i]], 'k', color = 'g', lw = 2)
legend = plt.legend(loc='best')
plt.show()
plt.clf()
plt.cla()
plt.close()




########################################################################
#                Problem 3
#
########################################################################





# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 20:17:09 2016

@author: ka746940
"""

import numpy as np
from numpy.linalg import inv


###########################################################
#   Problem 3 Part 1
#       Question 1
###########################################################
def gxy(x):
    return(4*x[0]*x[1]+(x[0]+x[1]**2)**2)
    


def gradient_hessian(X):
    """Calculate gradient and Hessian"""
    gred1 = 4*X[1]+2*(X[0]+X[1]**2)
    gred2 = 4*X[0]+4*X[1]*(X[0]+X[1]**2)
    gred = np.array([gred1, gred2])
    hessi11 = 2.0
    hessi12 = 4*X[1]+4
    hessi21 = 4*X[1]+4
    hessi22 = 4*X[0]+12*X[1]**2
    hessi = np.array([[hessi11, hessi12], [hessi21, hessi22]])
    return {"gradient": gred, "hessian": hessi}


def newton_raphson(func, f_gred, X = None,  tol = 1e-6, maxite = 200):
    if (X == None):
        X = np.repeat(0, 2)
        
    if not isinstance(X, np.ndarray):
        try:
            X = np.array(X, dtype='f')
        except:
            print ("X is not possible to convert as array")
            return
    
    update = 999999
    conv = update
    itter = 0

    while (itter < maxite  and conv> tol):
        itter = itter + 1
        gre_heis = f_gred(X)
        update = np.dot(inv(gre_heis["hessian"]), gre_heis["gradient"])
        X = X - update
        conv = max(np.absolute(update))
    
    if (itter == maxite):
        print("Failed to converge")
        return
    else:
        g = func(X)
        gr_hes = f_gred(X)
        return({"prem": X, "gradient": gr_hes["gradient"], \
        "Hessian": gr_hes["hessian"], "iteration": itter, "objective": g})

newton_raphson(gxy, gradient_hessian, [1,1])


#       Question 2
###########################################################


def gradient(X):
    """Calculate gradient"""
    gred1 = 4*X[1]+2*(X[0]+X[1]**2)
    gred2 = 4*X[0]+4*X[1]*(X[0]+X[1]**2)
    gred = np.array([gred1, gred2])
    return gred




def steep_des_GS(func, f_gred, x, a, b, tol = 1e-6, maxite = 200):
    def update(X, gred, alpha):
        upd = X - alpha * gred
        return(upd)
    if not isinstance(x, np.ndarray):
        try:
            x = np.array(x, dtype='f')
        except:
            print ("X is not possible to convert as array")
            return
    bb = b
    aa = a
    upd = (a+b)/2*f_gred(x)
    tau = (np.sqrt(5)-1.0)/2
    itter = 0
    while (np.sqrt(np.dot(upd, upd)) > tol and itter < maxite):
        a = aa
        b = bb
        gred = f_gred(x)
        x1 = a + (1.0-tau)*(b-a)
        x2 = a + tau * (b-a)
        f1 = func(update(x, gred, x1))
        f2 = func(update(x, gred, x2))
        while((b-a) > tol):
            if (f1 > f2):        
                a = x1
                x1 = x2
                f1 = f2
                x2 = a + tau * (b-a)
                f2 = func(update(x, gred, x2))
            else:
                b = x2
                x2 = x1
                f2 = f1
                x1 = a + (1.0-tau) * (b-a)
                f1 = func(update(x, gred, x1))
        upd = ((a+b)/2.0)*gred
        x = x - upd
        itter += 1

    if (itter == maxite):
        print("Failed to converge")
        return
    else:
        g = func(x)
        gr_hes = f_gred(x)
        return({"prem": x, "gradient": gr_hes, "iteration": itter, "objective": g})

steep_des_GS(gxy, gradient, np.array([1,1]), 0, 10)

'''###########################################################'''



###########################################################
#   Problem 3 Part 2
#       Question 1
###########################################################


def gxy(x):
    return(np.log(x)/(1+x))
    

def gradient_hessian(X):
    """Calculate gradient and Hessian"""
    gred = np.array(((1+X)/X-np.log(X))/((1+X)**2))
    hessi = np.array(((-X**(-2) - X**(-1)) * (1 + X)**2 - 2 * (1 + X) * ((1 + X) / X - np.log(X))) / ((1 + X)**4))
    return {"gradient": gred, "hessian": hessi}


def newton_raphson(func, f_gred, X = None,  tol = 1e-6, maxite = 200):
    if (X == None):
        X = np.repeat(0, 2)
        
    if not isinstance(X, np.ndarray):
        try:
            X = np.array(X, dtype='f')
        except:
            print ("X is not possible to convert as array")
            return
    
    update = 999999
    conv = update
    itter = 0

    while (itter < maxite  and conv> tol):
        itter = itter + 1
        gre_heis = f_gred(X)
        update = np.dot((1./gre_heis["hessian"]), gre_heis["gradient"])
        X = X - update
        conv = np.absolute(update)
    
    if (itter == maxite):
        print("Failed to converge")
        return
    else:
        g = func(X)
        gr_hes = f_gred(X)
        return({"prem": X, "gradient": gr_hes["gradient"], \
        "Hessian": gr_hes["hessian"], "iteration": itter, "objective": -g})

newton_raphson(gxy, gradient_hessian, [1])



###########################################################
#   Problem 3 Part 2
#       Question 2
###########################################################

def gradient(X):
    """Calculate gradient."""
    gred = np.array(((1 + X) / X - np.log(X)) / ((1 + X)**2))
    return (gred)


def secant_opt(func, f_gred, X,  tol = 1e-6, maxite = 200):
        
    if not isinstance(X, np.ndarray):
        try:
            X = np.array(X, dtype='f')
        except:
            print ("X is not possible to convert as array")
            return

    update = 999999999
    itter = 0
    conv = update

    while (itter < maxite  and conv> tol):
        itter = itter + 1
        gre_heis0 = f_gred(X[0])
        gre_heis1 = f_gred(X[1])
        hessi = ((X[1] - X[0])/(gre_heis1-gre_heis0))
        update = hessi*gre_heis1
        X[0] = X[1]
        X[1] = X[0] - update
        conv = np.absolute(update)
    
    if (itter == maxite):
        print("Failed to converge")
        return
    else:
        g = func(X[1])
        gre_heis = f_gred(X[1])
        return({"prem": X[1], "gradient": gre_heis, \
        "iteration": itter, "objective": -g})

secant_opt(gxy, gradient, [1, 1.1])


