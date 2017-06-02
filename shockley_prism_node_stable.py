#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import shockley_pair as shpair
import prism_pair as prpair

def bisect(x,y):
    """ return the root of the data (x,y) using the bisection method"""
    # sort the data according to y

    for i in range(len(y)-1):
        if np.sign(y[i+1]) != np.sign(y[i]):
            return ( (x[i]+x[i+1])/2.0, (y[i]+y[i+1])/2.0 )
    return None

npnts = 1000
beta = np.array(range(0,npnts))/float(npnts-1)*np.pi
E_shock = [0]*npnts
E_prism = [0]*npnts

# shpair and prpair have different energy units due to the prefector
#  and due to the norm of b
# New energy unit is E0 = mu*a^2/PI
# b_prism^2 = a^2/4, and b_shock^2 = a^2/3
# correct for E0: E0_sh = mub^2/8/pi, E0_pr = mub^2/4/pi
#  E0_sh = mua^2/24/pi and E0_pr = mua^2/16/pi
gamma_B = 198e-3/(131e9*3.232e-10)*np.pi
gamma_P = 135e-3/(131e9*3.232e-10)*np.pi

nu = 0.333333

for i in range(npnts):
    E_shock[i] = shpair.Es(nu,beta[i]) + shpair.Eif(nu,beta[i],gamma_B)
    E_prism[i] = prpair.Es(nu,0) + prpair.Eif(nu,0,gamma_P)

E_shock = np.array(E_shock)*24
E_prism = np.array(E_prism)*16

dE_shock_dB = np.diff(E_shock)/np.diff(beta) #derv E_D^B wrt beta
E_shock = E_shock[:len(E_shock)-1] #len(derivative) = N-1
E_prism = E_prism[:len(E_prism)-1] #len(derivative) = N-1
beta = beta[:len(beta)-1]

F_node = np.cos(beta)*E_shock - np.sin(beta)*dE_shock_dB-E_prism
root = bisect(beta,F_node)

fig1 = plt.figure()
plt.xlabel(r'$\theta$'+ ' (rad)')
plt.ylabel(r'$\Sigma$'+'Force '+r'$\left(\frac{\mu a^2}{\pi}\right)$')

ax = fig1.gca()

#ax.plot(beta,E_shock,label="basal disloc")
#ax.plot(beta,E_prism,label="prism disloc")
ax.plot(beta,F_node,label="prism-basal node")
ax.text(root[0]+0.1,root[1]+0.1,"(%1.4f, %1.2f)"%(root[0],root[1]))
ax.axhline(y=root[1],linewidth=0.5,color='k')
ax.axvline(x=root[0],linewidth=0.5,color='k')
ax.legend(loc=0)
plt.show()
