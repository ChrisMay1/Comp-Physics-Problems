#!/usr/bin/python2.7
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import time

t1=time.time()
def duffInt(x0,x1,dt,pts):
    tjh=time.time()
    m = 1
    a = 1/2
    b = 1/4
    F = 2.0
    omega = 2.4
    gamma = 0.1
    poin = {'pos':[],'vel':[]}
    period = 2*np.pi/(omega)

    tf = pts*(2*np.pi/omega)
    t = np.arange(0,tf,dt)
    x = np.zeros(len(t))
    v = np.zeros(len(t))
    x[0] = x0; x[1] = x1; v[1] = (x1-x0)/(2*dt)
    wtf = int(len(t)/pts)
    print('Initial Cond.: %f, Periods: %ld, dt: %e'%(x0,pts,dt))
    k = 1
    for n in range(1,len(t)-1):
        x[n+1]=2*x[n]-x[n-1]+((dt**2)/m)*(F*np.cos(omega*t[n])+2*a*x[n]-4*b*(x[n]**3)-gamma*v[n])
        v[n+1] = (x[n+1]-x[n])/(2*dt)
        if abs(t[n] - k*period) < dt:
            poin['pos'].append(x[n])
            poin['vel'].append(v[n])
            k += 1
        if n % int(0.1*pts*wtf) == 0:
            tjh = time.time()
            print'Percent Complete: ',(round(100*n/len(t),2)),' Time Elapsed: ',(round((tjh-t1),5))
            print('Time Remaining: %g'%float((len(t)-n)/(6000*(tjh-t1))))
    return t,x,v,poin

x0 = [0.5,0.5001]
dt = 1e-2
periods = 25000

t,x,v,poin = duffInt(x0[0],x0[1],dt,periods)
t2 = time.time()
print ('Time to compute strange attractor: %g\n\n'%(t2-t1))

def dimension(n):
    exist = {'x':[],'y':[]}
    asdf = np.linspace(-3,3,2*(2**n)+1)
    for ix in range(len(asdf)-1):
        for iy in range(len(asdf)-1):
            xmin = asdf[ix]; xmax = asdf[ix+1]
            ymin = asdf[iy]; ymax = asdf[iy+1]
            for k in range(len(poin['pos'])):
                xi = poin['pos'][k]; vi = poin['vel'][k]
                if (xi >= xmin and xi <= xmax) and (vi >= ymin and vi <= ymax):
                    exist['x'].append(ix)
                    exist['y'].append(iy)
                    break
    return exist
# The dimension calculation takes forever and isn't accurate.
##N = np.zeros(7)
##B = np.zeros(7)
##tdim = time.time()
##for l in range(1,8):
##    tl = time.time()
##    print('Starting l: %d/7, Dimension Calc. Time Elapsed: %g'%(l,(tl-tdim))) 
##    N[l-1] = np.log(len(dimension(l)['y']))
##    B[l-1] = np.log(6/l)
##t4 = time.time()
##print ('Time to compute dimension: %g\n\n'%(t4-t2))
##der = 0
##for d in range(6):
##    der += ((N[d+1]-N[d])/(B[d+1]-B[d]))/6

plt.figure()
plt.plot(poin['pos'],poin['vel'],'.k',markersize = 0.75)
plt.title('Strange Attractor')
plt.xlabel('x')
plt.ylabel('p')
plt.xlim([-3,3])
plt.ylim([-3,3])
plt.grid(True)

##plt.figure()
##plt.plot(B[:],N[:])
##plt.xlabel('log(b)')
##plt.ylabel('log(N(b))')
##plt.grid(True)
##plt.title('Fractal Dimension: %g'%der)
t3 = time.time()
print('Time to plot %g:'%(t3-t2))
t5 = time.time()
print('Total simulation time: %g'%(t5-t1))
plt.show()
