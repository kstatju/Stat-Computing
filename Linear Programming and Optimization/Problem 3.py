# alternative way

# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 23:58:23 2016

@author: Kanak
"""
import numpy as np
import pandas as pd
x = [28, -44, 29, 30, 26, 27, 22, 23, 33, 16, 24, 40, 21, 31, 34, -2, 25, 19]

def tstat(x, theta = None, sigsq = None):
    if theta is None:
        theta = x.mean()
    if sigsq is None:
        sigsq = x.var()
    
    return (((x-theta)**2)/sigsq)

def baggging(x, size = None, alpha = 0.05, boottimes = 100):
    if not isinstance(x, pd.DataFrame):
        try:
            x = pd.DataFrame(x)
        except:
            print ("x is not possible to convert as Data Frame")
            return
    if size == None:
        size = x.shape[0]
    hi = []
    for _ in range(boottimes):
        y = x.sample(n = size, replace = True)
        m = y.mean().values.flatten()
        sd = y.var().values.flatten()
        t = tstat(y, theta = m, sigsq = sd)
        hi.append(t.quantile(q = 1-alpha).values.flatten())
    h = pd.DataFrame(hi).mean()
    
    return {"h": h.values.flatten(),
            "Mean": x.mean().values.flatten(),
            "Variance": x.var().values.flatten()}

def predictoutlier(model, x):
    par = model['h']
    m = model["Mean"]
    sd = model["Variance"]
    t = np.array(tstat(x, m, sd))
    if np.ndarray.min(t) < 0 or np.ndarray.max(t) > np.ndarray.min(par):
        out = "Outlier"
        com = ("The observation {} with t-value {} is an OUTLIER because"+
              " given value does not fall between 0 and {}").format(x, t, par)
        print(com)
    else:
        out = "Not Outlier"
        com = ("The observation {} with t-value {} is NOT OUTLIER" +
              " because given value falls between 0 and {}").format(x, t, par)
        print(com)
    return {"h": par,
            "Mean": m,
            "Variance": sd,
            "t": t, 
            "Decision": out, 
            "Comment": com}        
                 
           
model = baggging(x, boottimes = 9999)
print(model)

pred = predictoutlier(model, 38)
print(pred)
