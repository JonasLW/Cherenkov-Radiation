# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 10:58:13 2018

@author: Jonas
"""
from __future__ import division
import pylab as py    
from scipy import integrate
from scipy import optimize

# NOTE: These functions employ natural units. c = hbar = 1. Give masses and frequencies
#       in units of eV.

# Functions for basic quantities related to the physical system
def Gamma( v ):
    """Returns the Lorentz factor as a function of velocity."""
    return 1/py.sqrt(1-v**2)

def CosChCl( n, v ):
    """Returns classical value of cosine of the Cherenkov angle."""
    result = 1/(n*v)
    if py.ndim(result) > 0:
        result[ result>1 ] = py.nan
        result[ result<0 ] = py.nan
    elif result > 1 or result < 0:
        result = py.nan
    return result
    
def CosCh( n, m, v, w ):
    """Returns cosine of Cherenkov angle as function of refraction coefficient, mass, velocity, and frequency."""
    f = 2*Gamma(v)*m    # For readability
    result = ( 1 + (n**2-1)*w/f )/(n*v)
    if py.ndim(result) > 0: 
        result[ result>1 ] = py.nan
        result[ result<0 ] = py.nan
    elif result > 1 or result < 0:
        result = py.nan
    return result
    
def ChAngleCl( n, v ):
    """Returns classical Cherenkov angle as function of refraction coefficient, mass, velocity, and frequency."""
    return py.arccos( CosChCl(n,v) )  # Angle is always positive

def ChAngle( n, m, v, w ):
    """Returns Cherenkov angle as function of refraction coefficient, mass, velocity, and frequency."""
    return py.arccos( CosCh(n,m,v,w) )  # Angle is always positive
    
def SinCh( n, m, v, w ):
    """Returns sine of Cherenkov angle as function of refraction coefficient, mass, velocity, and frequency."""
    return py.sin( ChAngle(n,m,v,w) )
    

# Amplitudes summed over final state spin
# Expressions do not include S_k, as this cancels in calculating physical quantities
def F( n, m, v, w, chi, phi ):
    """Function containing the spinor angle dependence of the squared amplitudes."""
    return py.sin(chi)*py.cos(phi)*SinCh(n,m,v,w)/Gamma(v) - py.cos(chi)*( n*v-CosCh(n,m,v,w) )

def TermOne( n, m, v, w ):
    """Recurring term in amplitudes."""
    return 0.5*( v*SinCh(n,m,v,w) )**2
    
def TermTwo( n, m, v, w ):
    """Recurring term in amplitudes."""
    return 0.5*w/( Gamma(v)*m )
    
def AmpMinSq( n, m, v, w, chi, phi, a ):
    """Returns squared transition amplitude of left-elliptical light, as function of spinor angles phi and chi, and basis a."""
    A = TermOne( n, m, v, w )
    B = TermTwo( n, m , v, w )
    return A + (n**2-1)*B**2 - py.cos(2*a)*A + py.sin(2*a)*B*F(n,m,v,w,chi,phi) 
    
def AmpPlusSq( n, m, v, w, chi, phi, a ):
    """Returns squared transition amplitude of right-elliptical light, as function of spinor angles phi and chi, and basis a."""
    A = TermOne( n, m, v, w )
    B = TermTwo( n, m, v, w )
    return A + (n**2-1)*B**2 + py.cos(2*a)*A - py.sin(2*a)*B*F(n,m,v,w,chi,phi)
    
def AmpTotSq( n, m, v, w ):
    """Returns squared transition amplitude summed over polarization states."""
    A = TermOne( n, m, v, w )
    B = TermTwo( n, m, v, w )
    return 2*A + 2*(n**2-1)*B**2  
    

# Helper functions used in calculation of individual amplitudes
def ACrossPlusThree( n, m, v, w, chi, phi, a ):
    """Three-component of cross product ['a'X'epsilon'].
    
    Note: the vector a in the cross-product is not the polarization parameter a."""
    g = Gamma(v)
    E = g*m
    p = g*m*v
    sinchi = py.sin(chi)
    coschi = py.cos(chi)
    sinphi = py.sin(phi)
    cosphi = py.cos(phi)
    sint = py.sin( ChAngle(n,m,v,w) )
    cost = CosCh(n,m,v,w)
    sina = py.sin(a)
    cosa = py.cos(a)
    one = p*( sina*cosphi*sinchi - 1j*cost*sinphi*sinchi*cosa )
    two = n*(E+m)*( sina*sint*coschi - sina*sinchi*cosphi*cost + 1j*sinchi*sinphi*cosa )
    return -1j*w*( one + two )
    
def Eps13( n, m, v, w, chi, phi ):
    """Three-component of epsilon_1."""
    sint = py.sin( ChAngle(n,m,v,w) )
    return sint*py.cos(chi) - CosCh(n,m,v,w)*py.sin(chi)*py.cos(phi)

def Eps23( chi, phi ):
    """Three-component of epsilon_2."""
    return -py.sin(chi)*py.sin(phi)
    
def K3( n, m, v, w, chi, phi ):
    """Three-component of k."""
    sint = py.sin( ChAngle(n,m,v,w) )
    return CosCh(n,m,v,w)*py.cos(chi) + sint*py.sin(chi)*py.cos(phi)
    
    
# Individual amplitudes   
def AmpOnePlus( n, m, v, w, chi, phi, a ):
    """Returns amplitude of the ( 1 -> 1,+ ) process."""
    E = Gamma(v)*m
    p = Gamma(v)*m*v
    B = 2*E+2*m-w
    dot = p*py.cos(a)*SinCh(n,m,v,w)
    cross = ACrossPlusThree(n,m,v,w,chi,phi,a)
    amp = dot*B + 1j*cross
    fac = 2*E*py.sqrt( (E+m)*(E+m-w) )  # Factor (M_k)/(S_k). Aligns the prefactor with that of squared amplitudes
    return -amp/fac    
    
def AmpTwoPlus( n, m, v, w, chi, phi, a ):
    """Returns amplitude of the ( 1 -> 2,+ ) process."""
    E = Gamma(v)*m
    p = Gamma(v)*m*v
    A = p*CosCh(n,m,v,w) - n*(E+m)
    B = p*SinCh(n,m,v,w)
    e13 = Eps13( n, m, v, w, chi, phi )
    e23 = Eps23( chi, phi )
    k3 = K3( n, m, v, w, chi, phi )    
    re = -w*py.sin(a)*e23*( A*e13 - B*k3 )/py.sqrt( 1 - e23**2 )
    im = -w*( A*py.cos(a)*(1-e23**2) + A*py.sin(a)*k3 + B*py.sin(a)*e13 )/py.sqrt(1-e23**2)
    fac = 2*E*py.sqrt( (E+m)*(E+m-w) )  # Factor (M_k)/(S_k). Aligns the prefactor with that of squared amplitudes
    return (re + 1j*im)/fac

def AmpOneMin( n, m, v, w, chi, phi, a ):
    """Returns amplitude of the ( 1 -> 1,- ) process."""
    return 1j*AmpOnePlus( n, m, v, w, chi, phi, a+0.5*py.pi )
    
def AmpTwoMin( n, m, v, w, chi, phi, a ):
    """Returns amplitude of the ( 1 -> 2,- ) process."""
    return 1j*AmpTwoPlus( n, m, v, w, chi, phi, a+0.5*py.pi )
    
def Concurrence( n, m, v, w, chi, phi, a ):
    """Returns two times the absolute value of the determinant of the matrix of scattering amplitudes.
    
    Used as a measure of entanglement."""
    om = AmpOneMin( n, m, v, w, chi, phi, a )
    op = AmpOnePlus( n, m, v, w, chi, phi, a )
    tm = AmpTwoMin( n, m, v, w, chi, phi, a )
    tp = AmpTwoPlus( n, m, v, w, chi, phi, a )
    return 2*py.absolute(op*tm - om*tp)/AmpTotSq(n,m,v,w)


# Squared individual amplitudes
def AmpOnePlusSq( n, m, v, w, chi, phi, a ):
    """Returns squared amplitude of the ( 1 -> 1,+ ) process."""
    E = Gamma(v)*m
    p = Gamma(v)*m*v
    B = 2*E+2*m-w
    dot = p*py.cos(a)*SinCh(n,m,v,w)
    cross = ACrossPlusThree(n,m,v,w,chi,phi,a)
    amp = dot*B + 1j*cross
    fac = 4*E**2*(E+m)*(E+m-w)
    return py.absolute(amp)**2/fac

def AmpTwoPlusSq( n, m, v, w, chi, phi, a ):
    """Returns squared amplitude of the ( 1 -> 2,+ ) process."""
    return AmpPlusSq(n,m,v,w,chi,phi,a)-AmpOnePlusSq(n,m,v,w,chi,phi,a)

def AmpOneMinSq( n, m, v, w, chi, phi, a ):
    """Returns squared amplitude of the ( 1 -> 1,+ ) process."""
    return AmpOnePlusSq( n, m, v, w, chi, phi, a+0.5*py.pi )

def AmpTwoMinSq( n, m, v, w, chi, phi, a ):
    """Returns squared amplitude of the ( 1 -> 2,- ) process."""
    return AmpMinSq(n,m,v,w,chi,phi,a)-AmpOneMinSq(n,m,v,w,chi,phi,a)  
    

# Physical quantities resulting from amplitudes. Power and Rates.
# NOTE: ALl quantities are densities in w and phi.
#       Integrate over w and phi for total quantities
def PowFac( v, w ): 
    """Converts a squared transition amplitude into Power density, in SI-units [J/rad].
    
    Contains a factor 1/2pi as there is in general no integral over phi."""
    c = 299792458
    q = 1.6*10**(-19)
    h = 1.1*10**(-34)
    W = w*1.6*10**(-19)/h           # Converting unit from eV to 1/s
    return 10**(-7)*q**2*W*c/(v*2*py.pi)
    
def PowMin( n, m, v, w, chi, phi, a ):
    return PowFac(v,w)*AmpMinSq(n,m,v,w,chi,phi,a)
    
def PowPlus( n, m, v, w, chi, phi, a ):
    return PowFac(v,w)*AmpPlusSq(n,m,v,w,chi,phi,a)
    
def PowTot( n, m, v, w, chi, phi, a ):
    return PowFac(v,w)*AmpTotSq(n,m,v,w)
    
def PowTotQM( n, m, v, w, chi, phi, a ):
    c = 299792458
    q = 1.6*10**(-19)
    h = 1.1*10**(-34)
    W = w*1.6*10**(-19)/h           # Converting unit from eV to 1/s
    return 10**(-7)*q**2*W*v*c*(1-1/(n**2*v**2))/(2*py.pi)


def RateFac( v ):
    """Converts a squared transition amplitude into a transition rate density, in SI-units [1/rad].
    
    Integrate transition rate density over w and phi for total transition rate of process."""
    c = 299792458
    q = 1.6*10**(-19)
    h = 1.1*10**(-34)
    return 10**(-7)*q**2*c/(h*v*2*py.pi)

def RateMin( n, m, v, w, chi, phi, a ):
    return RateFac(v)*AmpMinSq(n,m,v,w,chi,phi,a)
    
def RatePlus( n, m, v, w, chi, phi, a ):
    return RateFac(v)*AmpPlusSq(n,m,v,w,chi,phi,a)
    
def RateTot( n, m, v, w, chi, phi, a ):
    return RateFac(v)*AmpTotSq(n,m,v,w)
    
def RateTotQM( n, m, v, w, chi, phi, a ):
    E = w*1.6*10**(-19)             # Converting unit from eV to J
    return PowTotQM(n,m,v,w,chi,phi,a)/E  

def PlusMinDiff( n, m, v, w, chi, phi, a ):
    """Calculates the difference between rates in the plus/minus modes""" 
    return ( RatePlus(n,m,v,w,chi,phi,a)- RateMin(n,m,v,w,chi,phi,a) )


# Specific values of parameters
def MaxV( n, m, v, w, chi, phi, a ):
    """Returns velocity at which the difference between right-handed and left-handed photons is maximal."""
    D = lambda V: -py.absolute( PlusMinDiff(n,m,V,w,chi,0,a)/V )
    F = lambda V: -py.absolute( PlusMinDiff(n,m,V,w,chi,py.pi,a)/V )
    optD = optimize.minimize_scalar( D, method='bounded',bounds=(1/n,1) )
    optF = optimize.minimize_scalar( F, method='bounded',bounds=(1/n,1) )
    if optD.fun < optF.fun:
        return float(optD.x)
    else:
        return float(optF.x)
        
def MinChi( n, m, v, w ):
    """Returns angle, chi, which minimizes left-handed polarization amplitude."""
    den = Gamma(v)*( n*v-CosCh(n,m,v,w) )
    return py.arctan( SinCh(n,m,v,w)/den )  # chi always positive, <= pi
    
def MinAlpha( n, m, v, w ):
    """Returns polarization basis parameter, a, which minimizes left-handed polarization amplitude at chi=MinChi."""
    A = 0.5*( v*SinCh(n,m,v,w) )**2
    B = 0.5*w/( Gamma(v)*m )
    den = A + (n**2-1)*B**2
    return 0.5*py.arccos( A/den )   # alpha always positive, <= pi/2
    
def MaxW( n, m, v ):
    """Returns the maximum photon frequency permitted by energy conservation."""
    num = ( n*v-1 )*2*m*Gamma(v)
    den = n**2-1
    return num/den


# Utility functions for integration
def QuadArr( Func, a, b, args=()):
    """Integrates Func from a to b for args. One element of args may be array. 
    
    Outputs array of integrated values.
    args is tuple of floats and/or ints and/or numpy arrays.    
    Currently only one element of args may be array.
    The purpose of this function is to calculate an array of values where each element requires an integral."""
    l = 1
    for arg in args:
        if type(arg) is not float and type(arg) is not int:
            l = len(arg)
            break
    argarr = py.transpose( py.array([ py.ones(l)*arg for arg in args ]) )
    return py.array([ integrate.quad(Func,a,b,args=tuple(ars))[0] for ars in argarr ])
    
def IntegrateOmega( Func, n, m, v, wlims, chi, phi, a ):
    """Takes a function Func and integrates it over w.
    
    Func must take the appropriate arguments.
    wlims may be array used otherwise for plotting. Will then integrate over these values.
    Can not handle wlims input as meshgrid.
    Output in SI-units, however performs integral in natural units.
    Assumes Func output in SI-units."""
    h = 1.1*10**(-34)       
    fac = 1.6*10**(-19)/h   
    R = lambda W,N,M,V,CHI,PHI,A : py.nan_to_num( Func(N,M,V,W,CHI,PHI,A) )
    return fac*QuadArr( R, wlims[0], wlims[-1], args=(n,m,v,chi,phi,a) )  

def Delta( x, tolerance ):
    """Approximates delta function for plotting: Equals one when x is within tolerance of zero."""
    y = py.ones_like(x)
    y[ x > tolerance ] = 0
    y[ x <= -tolerance ] = 0
    return y

 