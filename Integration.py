# -*- coding: utf-8 -*-
"""
Created on Thu May 03 10:37:17 2018

@author: Jonas
"""
from __future__ import division
import pylab as py
import Functions as f
import Plotting as p
py.close('all')

# Script for plotting functions which involve a numerical integral over omega

# Values of parameters
# To plot over a variable (or two), set its (their) value(s) to a tuple: ( min, max ) 
# To make a sperical plot, set sphere to True.
res = 1000
sphere = False
n = 1.3                         # Refraction coefficient of medium
m = 5.11*10**5                  # Particle mass
v = 0.9                     # Particle velocity
w = 65.9                        # Photon frequency (eV)
wlims = ( 0, w )                # Integration limits for photon frequency
chi = py.pi/4                   # Spinor polar angle
phi = ( 0, 2*py.pi)                     # Spinor azimuthal angle
a = py.pi/4                     # Polarization basis parameter (alpha) 
    
n, m, v, w, chi, phi, a, proj, axes = p.VarProd( res, variables=(n,m,v,w,chi,phi,a) )

v = f.MaxV(n,m,v,w,chi,phi,a)
print v
c = 299792458
r = f.IntegrateOmega( f.PlusMinDiff,n,m,v,wlims,chi,phi,a )*2*py.pi/(v*c)
#r2 = f.IntegrateOmega( f.PlusMinDiff,1.03,m,v,wlims,chi,phi,a )
#r3 = f.IntegrateOmega( f.RateTot,n,m,v,wlims,chi,phi,a )
#r4 = f.IntegrateOmega( f.RateTot,1.03,m,v,wlims,chi,phi,a )
#R = f.IntegrateOmega( f.RateTot,n,m,v,wlims,chi,phi,a )*fac/(v*c)
r1 = (r >= 0)*r
r2 = (r < 0)*(-r)
ch = f.ChAngle( n, m, v, w )

z = ( r1, r2 )

labels = ( '$\Delta N_{1,+}(\phi)$', '$-\Delta N_{1,+}(\phi)$', '$\chi=\\frac{3}{8}\pi$', '$\chi=\\frac{1}{2}\pi$'  )
fig, ax = p.MultPlot( axes, z, title='',
                    labels=labels, xlabel='$\phi$', ylabel='', labelsize=25,
                    projection='polar', sphere=sphere )
ax.legend( bbox_to_anchor=(1.17, 1.17), prop={'size':17} ) # bbox_to_anchor=(1.17, 1.17) # for manual location
ax.tick_params( labelsize=13 )
#ax.set_ylim(-0.005,None)
#fig.savefig('PhotonNumPolChi14W17V86.png', bbox_inches='tight',dpi=600)
fig.show()












    
    