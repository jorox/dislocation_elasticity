#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
"""
Calculate the energy of a dislocation using Bacon line-tension model
"""

def E_lt(b,theta,mu,nu,R,r0):
    """ calculate energy per unit length using line tension model 
    """
    #print(b,theta,mu,nu,R,r0)
    Elt =  mu*b*b/4./np.pi/(1-nu)*(1-nu*np.cos(theta)*np.cos(theta))*np.log(R/r0)
    #print("Elt = %1.6e"%(Elt))
    return Elt

def E_blt(b1,b2,theta1,theta2,sfenergy,mu,nu,R,r0):
    """ calculate the energy per unit length using bacon line tension model
    """
    E0 = E_lt(b1,theta1,mu,nu,R,r0) + E_lt(b2,theta2,mu,nu,R,r0)
    F0 = mu*b1*b2/2/np.pi/(1-nu)
    gamma1 = np.cos(theta1)
    gamma2 = np.cos(theta2)
    gamma3 = np.cos(theta1-theta2)
    gammad = F0 * (gamma3 - nu * gamma1 * gamma2)

    
    d = gammad/sfenergy;
    E1 =  gammad*(1.0+np.log(R/d))
    
    efac = mu*3.232*3.232e-20/np.pi
    print("      E0            E1            F0            g1            g2             g3             gd              d       ")
    print("------------- ------------- ------------- ------------- ------------- --------------- ---------------- -------------")

    print("%1.6e, %1.6e, %1.6e, %1.6e, %1.6e, %1.6e, %1.6e,    %1.6e"%
          (E0/efac,E1/efac, F0, gamma1, gamma2, gamma3, gammad,d/3.232e-10))
    
    return E0+E1

def main():

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

    # x-coordinate 
    npnts = 2
    beta = [0., np.pi/2.]
    #npnts = 200
    #beta = np.array(range(0,npnts))/float(npnts-1)*np.pi/2.0

    # y-coordinate
    E_basal = [0]*npnts
    E_prism = [0]*npnts

    efac = mu*a*a/np.pi
    print("efac %1.5e"%(efac))
    for i in range(npnts):
        E_prism[i] = E_blt(b_prism,b_prism,beta[i]+alpha_prism,beta[i]-alpha_prism,
                           sfe_prism, mu, nu, 1000e-10,1e-10)/efac
        E_basal[i] = E_blt(b_basal,b_basal,beta[i]+alpha_basal,beta[i]-alpha_basal,
                           sfe_basal, mu, nu, 1000e-10,1e-10)/efac
       


    E_basal = np.array(E_basal)
    E_prism = np.array(E_prism)

    print("Prism screw/edge: %1.4e %1.4e" % (E_prism[0],E_prism[-1]))
    print("Basal screw/edge: %1.4e %1.4e" % (E_basal[0],E_basal[-1]))

    fig1 = plt.figure()
    plt.xlabel(r'$\theta$'+ ' (rad)')
    plt.ylabel('Elastic energy '+r'$\left(\frac{\mu a^2}{\pi}\right)$')

    ax = fig1.gca()

    ax.plot(beta,E_basal,label="basal disloc")
    ax.plot(beta,E_prism,label="prism disloc")

    ax.legend(loc=0)
    plt.show()

if __name__=="__main__":
    main()
