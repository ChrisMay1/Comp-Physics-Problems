#!/usr/bin/python2.7
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import time

# (hbar^2/2m) = 1

rmax = 8
rmin = 0.5
dr = 1e-1
hbar2m = 1

def j_l(x,l):
    j = [0]*(l+2)
    j[0] = np.sin(x)/x
    j[1] = (np.sin(x)/(x*x)) - (np.cos(x)/x)
    if l == 0:
        return j[0]
    elif l == 1:
        return j[1]
    elif (l > 1):
        for n in range(1,l):
            j[n+1] = ((2*n+1)/x)*j[n] - j[n-1]
        return j[l]

def n_l(x,l):
    j = [0]*(l+2)
    j[0] = -np.cos(x)/x
    j[1] = (-np.cos(x)/(x*x)) - (np.sin(x)/x)
    if l == 0:
        return j[0]
    elif l == 1:
        return j[1]
    elif (l > 1):
        for m in range(1,l):
            j[m+1] = ((2*m+1)/x)*j[m] - j[m-1]
        return j[l]

def p_l(x,l):
    p = [0]*(l+2)
    p[0], p[1] = 1,x
    if l == 0:
        return p[0]
    elif l == 1:
        return p[1]
    elif (l > 1):
        for n in range(1,l):
            p[n+1] = ((2*n+1)*x*p[n] - n*p[n-1])/(n+1)
        return p[l]

def delta_l(r1,r2,u1,u2,E,l):
    k = np.sqrt(E)
    K = (r1*u2)/(r2*u1)
    d_l = np.arctan((K*j_l(k*r1,l)-j_l(k*r2,l))/(K*n_l(k*r1,l)-n_l(k*r2,l)))
    return d_l

def F(l,r,E):
    F = (hbar2m*(l*l+l)/(r*r)) - E + V(r)
    return F

def V(r):
    rho = 3.57
    eps = (5.9*6.12)/(3.57*3.57)
    pr = rho/r
    pr2 = pr*pr
    pr6 = pr2*pr2*pr2
    pr12 = pr6*pr6
    V = eps*(pr12-2*pr6)
    return V

def wave(rmin,rmax,dr,E,l):
    half_wave = np.pi/np.sqrt(E)
    r = np.arange(rmin,rmax+half_wave,dr)
    w,w[1] = np.zeros(len(r)), dr**(l+1)
    norm = w[1]*w[1]*dr
    for n in range(1,len(w)-1):
        f = F(l,r[n],E)
        u = w[n]/(1-(f*dr*dr/12))
        w[n+1] = 2*w[n] - w[n-1] + dr*dr*f*u
        norm += w[n+1]*w[n+1]*dr
    w = w/np.sqrt(norm)
    u1,r1 = w[int((rmax-rmin)/dr)+1],r[int((rmax-rmin)/dr)+1]
    u2,r2 = w[len(w)-1],r[len(w)-1]
    return r[0:int((rmax-rmin)/dr)],w[0:int((rmax-rmin)/dr)],u1,u2,r1,r2

de = 0.01
E = np.arange(0.1,3.5,de)
sigma = np.zeros(len(E))
for e in range(len(E)):
    if e % int(len(E)/10) == 0:
        print('Completion ratio: %.2f'%(float(e/len(E))))
    for l in range(10):
        r,w,u1,u2,r1,r2, = wave(rmin,rmax,dr,E[e],l)
        dl = delta_l(r1,r2,u1,u2,E[e],l)
        sigma[e] += np.sin(dl)*np.sin(dl)*(2*l+1)
    sigma[e] = sigma[e]*4*np.pi/(E[e]*3.57*3.57)

r,w,u1,u2,r1,r2 = wave(rmin,rmax,dr,max(E),l)
plt.figure()
plt.plot(E,sigma)
plt.xlabel('Energy meV')
plt.ylabel(r'Total Cross Section ($\rho^2$)')
plt.title('Total Cross Section as function of energy for LJ potential')

plt.figure()
plt.plot(r,w,'k')
plt.plot(r[1:],V(r[1:]),'--r')
plt.ylim([min(V(r[1:])),max(w)+1])
plt.xlim([0,10])
plt.title('Wave Function Energy %.3g'%(max(E)))
plt.legend(["u(r)","V(r)"])
plt.xlabel('r')
plt.ylabel('W(r)')
plt.show()
