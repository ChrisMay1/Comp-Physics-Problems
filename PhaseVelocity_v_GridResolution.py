#!/usr/bin/python2.7
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import math
## This script is meant to be a solution to problem 2.3 in Comput. Elect.

def vphase(R,cdt):
    k = R*np.arccos(1+(1/(cdt*cdt))*(np.cos(2*np.pi*cdt/R)-1))
    v = 2*np.pi/k
    return v,k

def vgroup(R,cdt):
    k = vphase(R,cdt)[1]
    vg = cdt*np.sin(k/R)/np.sin(2*np.pi*cdt/R)
    return vg

cdt = [0.99,0.9,0.5,0.1]
R = np.arange(3,20,1)

plt.figure()
for i in cdt:
    plt.plot(R,vphase(R,i)[0])
plt.xlabel('Grid Space Resolution (# of grid points)')
plt.ylabel('Numerical Phase Velocity')
plt.title('Dependence of Numerical Phase Velocity on Spatial Resolution')
plt.legend([cdt[0],cdt[1],cdt[2],cdt[3]])

plt.figure()
for i in cdt:
    plt.plot(R,vgroup(R,i))
plt.xlabel('Grid Space Resolution (# of grid points)')
plt.ylabel('Numerical Group Velocity')
plt.title('Dependence of Numerical Group Velocity on Spatial Resolution')
plt.legend([cdt[0],cdt[1],cdt[2],cdt[3]])
plt.show()



