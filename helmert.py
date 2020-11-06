"""Calculate of parametres of Helmert transformation from point pairs"""

from math import sqrt, atan2, sin, cos
from statistics import median
import numpy as np
from numpy.linalg import inv

def trp_helmert(ptl, trp):
    """Calculate Helmert transformation of ptl [(x,y,z), ...]
       by parameters trp (dx, dy, a, b)"""
    return [(trp[0]+pt[0]*trp[2]+pt[1]*trp[3],
             trp[1]-pt[0]*trp[3]+pt[1]*trp[2]) for pt in ptl]

def helmert_lsq(ptl1, ptl2):
    """Calculates parameters from 1 to 2 coordsys, by point lists
       (ptl1 and ptl2 [(x,y,z), ...]), less squeres method"""
    npt=len(ptl1)
    if len(ptl2)!=npt:
        raise ValueError('Different list length (ptl1 and ptl2)')
    prex=(0.0, 0.0, 1.0, 0.0)
    a=np.zeros([npt*2,4], dtype=np.double)
    p=np.identity(npt*2, dtype=np.double)
    l=np.zeros([npt*2,1], dtype=np.double)
    for i in range(npt):
        a[2*i  ][0] = 1.0
        a[2*i  ][1] = 0.0
        a[2*i  ][2] = ptl1[i][0]
        a[2*i  ][3] = ptl1[i][1]
        a[2*i+1][0] = 0.0
        a[2*i+1][1] = 1.0
        a[2*i+1][2] = ptl1[i][1]
        a[2*i+1][3] =-ptl1[i][0]
        l[2*i  ][0] =(prex[0]+prex[2]*ptl1[i][0]+prex[3]*ptl1[i][1]-ptl2[i][0])
        l[2*i+1][0] =(prex[1]-prex[3]*ptl1[i][0]+prex[2]*ptl1[i][1]-ptl2[i][1])
    x=-inv(a.transpose()@p@a)@(a.transpose()@p@l)
    return (prex[0]+x[0][0], prex[1]+x[1][0], prex[2]+x[2][0], prex[3]+x[3][0])
    

def wmed(wml):
    """Wighted median from wml [(weight, vaule), ...]"""
    wmls=sorted(wml, key=lambda element:element[1])
    wmsh=sum(map(lambda element:element[0], wmls))/2.0
    wmsi=0.0
    for wme in wmls:
        wmsi+=wme[0]
        if wmsi>wmsh:
            return wme[1]
    
def helmert_rob(ptl1, ptl2):
    """Calculates parameters from 1 to 2 coordsys, by point lists
       (ptl1 and ptl2 [(x,y,z), ...]), robust method"""
    npt=len(ptl1)
    if len(ptl2)!=npt:
        raise ValueError('Different list length (ptl1 and ptl2)')
    dml=[]
    al=[]
    bl=[]
    for ipt in range(npt):
        for jpt in range(ipt):
            c1=complex(ptl1[jpt][0]-ptl1[ipt][0], ptl1[jpt][1]-ptl1[ipt][1])
            c2=complex(ptl2[jpt][0]-ptl2[ipt][0], ptl2[jpt][1]-ptl2[ipt][1])
            cr=(c2/c1).conjugate()
            dml.append( 0.5*(sqrt(c1.real**2+c1.imag**2)+sqrt(c2.real**2+c2.imag**2)) )
            al.append(cr.real)
            bl.append(cr.imag)            
    dmed=median(dml)
    wl=[min(dm, 1.5*dmed) for dm in dml]
    a=wmed(zip(wl,al))
    b=wmed(zip(wl,bl))
    dx=median([pp[1][0]-a*pp[0][0]-b*pp[0][1] for pp in zip(ptl1, ptl2)])
    dy=median([pp[1][1]+b*pp[0][0]-a*pp[0][1] for pp in zip(ptl1, ptl2)])
    return (dx, dy, a, b)
    
    
