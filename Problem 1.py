# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 14:00:27 2016

@author: Kanak
"""

##################################################################
#
# Problem 1                       Part 1
#
##################################################################


import pandas as pd
import numpy as np
import math

df = pd.read_csv("D:/UCF/STA 6106 Statistical Computing/Assignments/Midterm 2/Exam2_pb1.txt", header=None, names=["a"])

def boot_bias_sde(x, boottimes = 10, method = 'Mean', size = None):
    if not isinstance(x, pd.DataFrame):
        print ("x is not a Data Frame")
        return
    if size == None:
        size = x.shape[0]
    if method.upper() == 'MEAN':
        theta = x.mean()
        thetab = pd.DataFrame([(x.sample(n = size, replace = True)).mean() \
                               for _ in range(boottimes)])
    if method.upper() == 'MEDIAN':
        theta = x.median()
        thetab = pd.DataFrame([(x.sample(n = size, replace = True)).median() \
                               for _ in range(boottimes)])
    botse = thetab.std()
    bias = (thetab - theta).mean()
    return {"Parameter":  theta.values.flatten(), 
            "Bias": bias.values.flatten(),
            "Standard Error": botse.values.flatten()}
                               
    
a = boot_bias_sde(df, boottimes = 9999, method = 'median', size = None)

print (a)

                            
##################################################################
#
# Problem 1                       Part 2
#
##################################################################

def jackknife(x, method = 'Mean'):
    if not isinstance(x, pd.DataFrame):
        print ("x is not a Data Frame")
        return
    n = x.shape[0]
    if method.upper() == 'MEAN':
        theta = x.mean()
        thetab = pd.DataFrame([(x.drop(i)).mean() for i in range(n)])
    if method.upper() == 'MEDIAN':
        theta = x.median()
        thetab = pd.DataFrame([(x.drop(i)).median() for i in range(n)])
    
    jkse = ((n-1)/math.sqrt(n))*thetab.std()
    bias = (n-1)*(thetab.mean() - theta)
    

    return {"Parameter":  theta.values.flatten(), 
            "Bias": bias.values.flatten(),
            "Standard Error": jkse.values.flatten()}
    
jk = jackknife(df, method = "mean")

print(jk)

##################################################################
#
# Problem 1                      Part 3
#
##################################################################


def boot_ci(x, alpha = 0.05, boottimes = 500, r = 100, method = 'Mean'):
    if not isinstance(x, pd.DataFrame):
        try:
            x = pd.DataFrame(x)
        except:
            print ("x is not possible to convert as Data Frame")
            return
    size = x.shape[0]

    def boot_se(y, r, method, size):
        if method.upper() == 'MEAN':
            thetab = pd.DataFrame([(y.sample(n = size, replace = True)).mean() \
                               for _ in range(r)])
        if method.upper() == 'MEDIAN':
            thetab = pd.DataFrame([(y.sample(n = size, replace = True)).median() \
                                   for _ in range(r)])
        return thetab.std().values.flatten()
        
    se = []
    thetab = []
    if method.upper() == 'MEAN':
        theta = x.mean().values.flatten()
        for _ in range(boottimes):
            y =x.sample(n = size, replace = True)
            se.append(boot_se(y, r, method, size))
            thetab.append(y.mean().values.flatten())    

    if method.upper() == 'MEDIAN':
        theta = x.median().values.flatten()
        for _ in range(boottimes):
            y =x.sample(n = size, replace = True)
            se.append(boot_se(y, r, method, size))
            thetab.append(y.median()) 
            

    t = pd.DataFrame((thetab - np.mean(thetab))/ se)
    thetab = pd.DataFrame(thetab)
    sd = thetab.std().values.flatten()
    qinterval = thetab.quantile(q = (alpha/2, 1-alpha/2))
    tqt = t.quantile(q = (alpha/2, 1-alpha/2)).abs().values.flatten()
    bci = theta + (tqt * [-1, 1] * sd)
    bias = (thetab - theta).mean()
    ci = {"Parameter":  theta,
          "Bootstrap CI": bci, 
          "Percentile Interval": qinterval.values.flatten(),
          "Bias": bias.values.flatten(),
          "Standard Error": sd}
    return ci
    

ci = boot_ci(df, alpha = 0.05, boottimes = 500, method = 'median')

print(ci)