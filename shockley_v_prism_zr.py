#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import shockley_pair as shpair
import prism_pair as prpair

npnts = 200
beta = np.array(range(0,npnts))/float(npnts-1)*np.pi/2.0
E_shock = [0]*npnts
E_prism = [0]*npnts

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

fig1 = plt.figure()
plt.xlabel(r'$\theta$'+ ' (rad)')
plt.ylabel('Elastic energy '+r'$\left(\frac{\mu a^2}{\pi}\right)$')

ax = fig1.gca()

ax.plot(beta,E_shock,label="basal disloc")
ax.plot(beta,E_prism,label="prism disloc")

ax.legend(loc=0)
plt.show()
