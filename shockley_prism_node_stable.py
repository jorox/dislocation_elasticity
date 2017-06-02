#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import blt_energy as blt

def bisect(x,y):
    """ return the root of the data (x,y) using the bisection method"""
    # sort the data according to y

    for i in range(len(y)-1):
        if np.sign(y[i+1]) != np.sign(y[i]):
            return ( (x[i]+x[i+1])/2.0, (y[i]+y[i+1])/2.0 )
    return None
##################################################################################
# Common parameters
mu =34.0e9 #Pa
nu = 0.3970589126
a = 3.232e-10 #A

# Basal - Prism parameters - SFE and angle between b_perfect and b_partials
sfe_basal = 198e-3 #J/m^2
sfe_prism = 135e-3
alpha_basal = np.pi/6.0 #rad
alpha_prism = 0
b_basal = a/np.sqrt(3.0)
b_prism = a/2.0


npnts = 1000
beta = np.array(range(0,npnts))/float(npnts-1)*np.pi
E_shock = [0]*npnts
E_prism = [0]*npnts


efac = mu*a*a/np.pi
print("efac %1.5e"%(efac))
for i in range(npnts):
    E_prism[i] = blt.E_blt(b_prism,b_prism,beta[i]+alpha_prism,beta[i]-alpha_prism,
                           sfe_prism, mu, nu, 1000e-10,1e-10)/efac
    E_shock[i] = blt.E_blt(b_basal,b_basal,beta[i]+alpha_basal,beta[i]-alpha_basal,
                           sfe_basal, mu, nu, 1000e-10,1e-10)/efac

E_shock = np.array(E_shock)
E_prism = np.array(E_prism)

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
