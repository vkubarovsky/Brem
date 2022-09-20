#!/usr/local/bin/python3
import numpy as np
from matplotlib import pyplot as plt
import math
import os

#
#  H. W. Koch and J. W. Motz, Rev. Mod. Phys.31, 920 (1959).
#
r0=2.82E-13 
def M0(k,E0,Z):
    E=E0-k
    return 1./( math.pow(k/(2*E0*E),2)+math.pow((math.pow(Z,1./3.)/111.),2))

def b(k,E0,Z):
    E=E0-k
    return 2*E0*E*math.pow(Z,1./3.)/(111.*k)

def f(Z):
    a=Z/137.
    return a**2 * ( 1./(1+a**2) + 0.20206 - 0.0369*a**2 + 0.0083*a**4 - 0.002*a**6)

def X0(Z):
#    return 1. / ( 1./716.408 * (5.31-f(1.)+6.144))
    return 63.05      # g cm-2

def k_sig5(k,E0,Z):
    E=E0-k
    y=E/E0
    if(k<E0):
        xs = 1.0E27*(4*Z**2*r0**2)/137. * \
        (1+y**2-(2/3)*y) * (math.log(2*E0*E/k)-1/2)
#    (math.log(M0(k,E0,Z))+1-2/b(k,E0,Z)*math.atan(b(k,E0,Z))) + \
#    y * 0. *(2/b(k,E0,Z)**2*math.log(1+b(k,E0,Z)**2) + \
#    4 * (2-b(k,E0,Z)**2)/(3*b(k,E0,Z)**3)*math.atan(b(k,E0,Z)) - \
#    8/(3*b(k,E0,Z)**2) + 2/9) \
#    )
        return xs
    else:
        return 0.0

def k_sig1(k,E0,Z):
    E=E0-k
    y=E/E0
    xs = 1.0E27*(2*Z**2*r0**2)/137. * \
    (
    (1+y**2-(2/3)*y) * \
    (math.log(M0(k,E0,Z))+1-2/b(k,E0,Z)*math.atan(b(k,E0,Z))) + \
    y * (2/b(k,E0,Z)**2*math.log(1+b(k,E0,Z)**2) + \
    4 * (2-b(k,E0,Z)**2)/(3*b(k,E0,Z)**3)*math.atan(b(k,E0,Z)) - \
    8/(3*b(k,E0,Z)**2) + 2/9) \
    )
    if(k<E0 and xs>0):
        return xs
    else:
        return 0.0
    

def k_sig2(k,E0,Z):
    E=E0-k
    y=E/E0
    xs =  1.0E27 * (4*Z*Z*r0*r0)/137. * (( 1+y**2-(2/3)*y) * math.log(183./math.pow(Z,1/3))+(1/9)*y)
    if(k<E0 and xs>0):
        return xs
    else:
        return 0.0

def k_sig3(k,E0,Z):
    E=E0-k
    y=k/E0
    xs =  1.0E27 * 1/63.05/6.022E23 * \
    ( (4./3.)-(4./3.)*y+y**2)
    if(k<E0 and xs>0):
        return xs
    else:
        return 0.0

def k_sig4(k,E0,Z):
    E=E0-k
    y=k/E0
    xs =  1.0E27 * 4*r0**2/137. * \
    (( 4/3 - 4/3*y + y**2) * (Z**2*(5.31-f(1.)+Z*6.144)) + 1/9*(1-y)*(Z**2+Z)) 
    if(k<E0 and xs>0):
        return xs
    else:
        return 0.0

Z=1.
E_0=10.6
x = np.linspace(E_0*0.01, E_0*1.1, 1000)
e = x   * 1000/0.51
E0= E_0 * 1000/0.51
y1=x*0.; y2=x*0.; y3=x*0.; y4=x*0.; y5=x*0.; y6=x*0.

for i in range(np.size(x)):
    y1[i] = k_sig1(e[i],E0,Z)
    y2[i] = k_sig2(e[i],E0,Z)
    y3[i] = k_sig3(x[i],E_0,Z)
    y4[i] = k_sig4(x[i],E_0,Z)
    y5[i] = k_sig5(e[i],E0,Z)

#    z1[i] = k_sig3(e[i],E0,   Z)
#    z2[i] = k_sig2(e[i],E0/10,Z)

#    print(x[i],y1[i],y2[i])
plt.figure (num=0,dpi=120)
plt.title  ("Bremsstrahlung spectrum")
plt.xlabel ("k - Photon energy /GeV")
plt.ylabel (r'$k\cdot \frac{d\sigma}{dk},mb$')
#plt.xscale("log")

plt.plot(x,y1)
#plt.plot(x,y2)
#plt.plot(x,y3)
#plt.plot(x,y4)
plt.plot(x,y5)

#plt.plot(x,z1)
#plt.plot(x,z2)

#plt.show()
plt.ylim(bottom=0)
#os.system("rm -f brem.pdf")
plt.savefig('brem.pdf')
os.system("open -a preview brem.pdf")



    

