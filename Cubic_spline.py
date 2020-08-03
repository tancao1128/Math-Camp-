import numpy as np
from math import *

def clamped_cubic_spline(mydata):
    n = len(mydata)
    a = np.zeros(n,float)
    b = np.zeros(n,float)
    c = np.zeros(n,float)
    d = np.zeros(n,float)
    fpo = mydata[0][2]
    fpn = mydata[n-1][2]
    h = np.zeros(n,float)
    l = np.zeros(n,float)
    mu = np.zeros(n,float)
    z = np.zeros(n,float)
    alpha = np.zeros(n,float)
    for i in range(n):
        a[i] = mydata[i][1]
    for i in range(n-1):
        h[i] = mydata[i+1][0] - mydata[i][0]
        
    alpha[0] = 3*(a[1]-a[0])/h[0] - 3*fpo
    alpha[n-1] = 3*fpn - 3*(a[n-1]-a[n-2])/h[n-2]
    for i in range(1,n-1):
        alpha[i] = (3/h[i])*(a[i+1]-a[i])-(3/h[i-1])*(a[i]-a[i-1])
    l[0] = 2*h[0]
    mu[0] = 0.5
    z[0] = alpha[0]/l[0]
    for i in range(1,n-1):
        l[i] = 2*(mydata[i+1][0]-mydata[i-1][0])-h[i-1]*mu[i-1]
        mu[i] = h[i]/l[i]
        z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i]
    
    l[n-1] = h[n-2]*(2-mu[n-2])
    z[n-1] = (alpha[n-1]-h[n-2]*z[n-2])/l[n-1]
    c[n-1] = z[n-1]
    j = n-2
    while True:
        c[j] = z[j] - mu[j]*c[j+1]
        b[j] = (a[j+1] - a[j])/h[j] - h[j]*(c[j+1]+2*c[j])/3
        d[j] = (c[j+1]-c[j])/(3*h[j])
        j -=1
        if j == -1:
            break 
    return a, b, c, d