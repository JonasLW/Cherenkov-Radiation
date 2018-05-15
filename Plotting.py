# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 14:35:33 2018

@author: Jonas
"""
from __future__ import division
import pylab as py
import itertools
from mpl_toolkits.mplot3d import Axes3D

def VarProd( res, variables=() ):
    """Outputs desired values or ranges of values for variables to input in functions.
    
    If an element in variables is a tuple, it will be exchanged for a numpy array spanning its values.
    If two elements are tuples, they will be replaced by meshgrids.
    No more than two tuple elements allowed in input.
    Outputs also the necessary projection variable for a subsequent plot.
    Outputs also copies of the variables to be used as axes in a tuple."""
    vals = [ ]
    projection = [ None ]
    axes = [ ]
    counter = 0
    pos = [ ]
    for i, var in enumerate(variables):
        if type(var) is tuple:
            counter += 1
            pos.append(i)
            var = py.linspace( var[0], var[1], res )
            if counter == 2:
                vals[pos[0]], var = py.meshgrid( vals[pos[0]], var )
                projection = ['3d']
        vals.append(var)
    for i in pos:
        axes.append( vals[i] )
    axes = [tuple(axes)]
    return vals + projection + axes
    
py.ioff()
def SpherePlot( theta, phi, data, 
                title='', xlabel='x', ylabel='y', zlabel='z', color='b',
                fig=py.figure(), ax=None, pos=111 ):
    """Plots data in terms of polar angles THETA and PHI.
    
    x must be a tuple of meshgrids:
    Polar angle, THETA, as first element of x, 
    azimuthal angle, PHI, as second element of x.
    Note: Sphere plots require higher resolution than mesh plots. 300 is good."""
    X = py.sin(theta)*py.cos(phi)*data
    Y = py.sin(theta)*py.sin(phi)*data
    Z = py.cos(theta)*data
    if ax is None:
        ax = fig.add_subplot( pos, projection='3d' )
        ax.ticklabel_format(style='sci',scilimits=(-1,2),axis='both')      #
        ax.xaxis.major.formatter._useMathText = True                       # 
        ax.yaxis.major.formatter._useMathText = True                       # 
        ax.zaxis.major.formatter._useMathText = True                       # Sets scientific notation
        ax.plot_surface( X, Y, Z, color=color )
        oldlim = 0
    else:
        ax.plot_surface( X, Y, Z, color=color )
        oldlim = ax.get_xlim()[1]
    lim = py.amax(data)
    lim = max(lim,oldlim)
    ax.set_xlim(-lim,lim)
    ax.set_ylim(-lim,lim)
    ax.set_zlim(-lim,lim)
    return fig, ax   
    
def MultPlot( x, data, fig=py.figure(), ax=None, pos=111, projection=None, sphere=False,
              title='', labels=(), xlabel='x', ylabel='y', zlabel='z', labelsize=17,
              linewidth = 1.5, styles = ('solid','solid','solid','solid','solid'),
              colors = ('blue','red','cyan', 'orange', 'yellow', 'green' ) ):
    """Plots multiple data sets on one set of axes (x).
    
    Breaks if dfferent lengths between x and elements of data.
    Can also handle inputs as meshgrids for 3d plotting (must set projection to '3d').
    Can plot in spherical coordinates as well (set sphere to True)."""
    if ax is None:
        ax = fig.add_subplot( pos, projection=projection )
        if projection == 'polar':
            #ax.set_theta_zero_location("N")
            #ax.set_rlabel_position(-30)
            ax.ticklabel_format(style='sci',scilimits=(-2,2),axis='y')
            ax.yaxis.set_offset_position('right')
        else:
            ax.ticklabel_format(style='sci',scilimits=(-2,2),axis='both')   #
        ax.xaxis.major.formatter._useMathText = True                        # 
        ax.yaxis.major.formatter._useMathText = True                        # Sets scientific notation
    if ax.get_title() == '':
        ax.set_title( title, loc='center' )
    if type(data) is not tuple:
        data = ( data, )
    for dat, lab, clr, stl in itertools.izip_longest( data, labels, colors, styles ) :
        if dat is None:
            break
        elif sphere is True:
            fig, ax = SpherePlot( x[0], x[1], dat, fig=fig, ax=ax, color=clr )
            ax.zaxis.major.formatter._useMathText = True
            ax.set_zlabel( zlabel, size=labelsize )
        elif projection == '3d':
            Axes3D.plot_wireframe( ax, x[0], x[1], dat, label=lab, linewidth=linewidth, color=clr )
            ax.zaxis.major.formatter._useMathText = True
            ax.set_zlabel( zlabel, size=labelsize )
        else:
            ax.plot( x[0], dat, label=lab, linewidth=linewidth, color=clr, ls=stl )
    ax.set_xlabel( xlabel, size=labelsize )
    ax.set_ylabel( ylabel, size=labelsize )
    return fig, ax
py.ion()


