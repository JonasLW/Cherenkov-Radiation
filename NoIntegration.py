# -*- coding: utf-8 -*-
"""
Created on Wed May 09 14:31:51 2018

@author: Jonas
"""
from __future__ import division
import pylab as py
import Functions as f
import Plotting as p
py.close('all')

# Script for plotting functions which do not involve a numerical integral over omega

# Values of parameters:
# To plot over a variable (or two), set its (their) value(s) to a tuple: ( min, max ) 
# To make a sperical plot, set sphere to True.
res = 1000
sphere = False
n = 1.3                         # Refraction coefficient of medium
m = 5.11*10**5                  # Particle mass
v = 0.9                         # Particle velocity
w = 6.59                        # Photon frequency (eV)
chi = ( 0, py.pi )              # Spinor polar angle
phi = 0                         # Spinor azimuthal angle
a = py.pi/4                     # Polarization basis parameter (alpha)

if sphere is True:
    res = 300
    chi = ( 0, py.pi )
    phi = ( 0, 2*py.pi )
    v = 1.000001/n     # 
    w = 0.9*f.MaxW(n,m,v)       # Example values
    
n, m, v, w, chi, phi, a, proj, axes = p.VarProd( res, variables=(n,m,v,w,chi,phi,a) )

ch = f.ChAngle( n, m, v, w ) 
conc = f.Concurrence( n, m, v, w, chi, phi, a )
conc2 = f.Concurrence( n, m, v, w, chi, py.pi, a )

if sphere is True:
    k = 1.1*p.Delta( ch-axes[0], 0.02 )*p.Delta( axes[1], 0.02 )*py.amax( z )
    if type(z) is not tuple:
        z = ( z, )
    z = z + ( k, )
    
labels = ( '$C(\chi)$', '$\\theta_C$',
           '$|\mathcal{M}_{1,2-}|^2$', '$|\mathcal{M}_{1,2+}|^2$'  )
fig, ax = p.MultPlot( (chi,), conc, title='',
                    labels=labels, xlabel='$\chi$', ylabel='',
                    labelsize=25,
                    projection='polar', sphere=sphere )
fig, ax = p.MultPlot( (-chi,), conc2, title='',
                    labels=(), xlabel='$\chi$', ylabel='',
                    labelsize=25,
                    projection='polar', sphere=sphere,
                    fig=fig, ax=ax )
ax.legend( bbox_to_anchor=(1.10, 0.1), prop={'size':17} ) # bbox_to_anchor=(1.17, 1.17) # for manual location
ax.tick_params( labelsize=13.5 )
#fig.savefig('ConcW16V077.png', bbox_inches='tight',dpi=600)
fig.show()


