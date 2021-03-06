#!/usr/bin/python

"""
=================================
Elastic energy of a Shockley pair
=================================

This script plots the elastic energy as a function of the angle of
a mixed dislocation dissociated into two Shockley partials.
For example, a dislocation in the basal plane of zirconium

Parameters:
------------
E0 = mu*b^2/8/PI, take as te unit of energy [J/m]
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
    return (2-nu)/(1-nu)/gamma*(1-2*nu*np.cos(2*beta)/(2-nu))

def Es(nu, beta):
    """ returns the self energies for the pair of Shockley
    """
    return (4-nu-2*nu*np.cos(beta)*np.cos(beta))/(1-nu)*6.2

def Eif(nu,beta,gamma):
    """ interaction energy + stacking-fault energy
    """
    d = d_e(nu,beta,gamma)
    return d*gamma*(6.2-np.log(d)+1)

def main():
    # application to zirconium
    npnts = 200
    beta = np.array(range(0,npnts))/float(npnts-1)*np.pi/2.0
    E = [0]*npnts
    nu = 0.33333
    mu = 34e9 #Pa
    b = 3.232e-10 #A
    
    gamma = 198e-3/(mu*b*np.sqrt(3)/3)*8*np.pi
    

    for i in range(npnts):
        E[i] = Es(nu,beta[i]) + Eif(nu,beta[i],gamma)

    fig1 = plt.figure()
    plt.xlabel(r'$\beta$'+ ' (rad)')
    plt.ylabel('Elastic energy '+r'$\left(\frac{\mu b^2}{8\pi}\right)$')

    ax = fig1.gca()
    ax.plot(beta,E,label="Dissociated Basal")

    plt.show()
