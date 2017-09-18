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