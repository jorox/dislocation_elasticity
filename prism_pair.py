#!/usr/bin/python

"""
=================================
Elastic energy of a Prismatic Pair
=================================

This script plots the elastic energy as a function of the angle of
a mixed dislocation dissociated into two partials.
The Burgers of the partials are not rotated wrt the original Burgers 

An example application is a dislocation in the prismatic plane of zirconium

Parameters:
------------
E0 = mu*b^2/4/PI, take as the unit of energy [J/m]
nu = Poisson's ratio, unitless
gamma = Stacking-fault energy give in units of = E0/b 

we take R = 1000d and r0 = d/2
which makes ln(R/d) ~= to 3*2.3 = 6.9
and makes   ln(R/r0) = ln(R/2/d) = 3*2.3-ln(2) = 6.2
""" 

import numpy as np
import matplotlib.pyplot as plt
from basic_units import radians, degrees, cos

def d_e(nu,beta,gamma):
    """ returns the equilibrium separation between Shockley pair for 
    units = b
    """
    return 1/(1-nu)/gamma*(1+(2-nu)*np.cos(2*beta))

def Es(nu, beta):
    """ returns the self energies for the pair of Shockley
    """
    return 2/(1-nu)*(1-nu*np.cos(beta)*np.cos(beta))*6.2

def Eif(nu,beta,gamma):
    """ interaction energy + stacking-fault energy
    """
    return d_e(nu,beta,gamma)*gamma*(6.9+1)

def main():
    npnts = 200
    beta = np.array(range(0,npnts))/float(npnts-1)*np.pi/2.0
    E = [0]*npnts

    gamma = 135e-3/(131e9*3.232e-10/2)*4*np.pi
    nu = 0.33333

    for i in range(npnts):
        E[i] = Es(nu,0) + Eif(nu,0,gamma)

    fig1 = plt.figure()
    plt.xlabel(r'$\beta$'+ ' (rad)')
    plt.ylabel('Elastic energy '+r'$\left(\frac{\mu b^2}{8\pi}\right)$')

    ax = fig1.gca()
    ax.plot(beta,E,label="Dissociated Basal")

    plt.show()

if __name__=="__main__":
    main()
