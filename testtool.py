"""Helper functions for testing transformation methods"""

import random
from math import pi, sin, cos, sqrt

def trcoordlist(ptl1, ptl2, trfunc, olpt=()):
    """Compare the points of ptl1 transformed by trfunc,
and the points of ptl2. Calculate the squared mean with
and without outliers, whose indices are in the olpt.
The result is printes to stdout"""
    trp=trfunc(ptl1, ptl2)
    print('The parameters of the transformation: ',trp)
    sumn=0.0
    sumnall=0.0
    for pp in zip(ptl1, ptl2, [(i in olpt) for i in range(len(ptl1))]):
        x1=pp[0][0]
        y1=pp[0][1]
        x2=pp[1][0]
        y2=pp[1][1]
        x2r=trp[0]+trp[2]*x1+trp[3]*y1
        y2r=trp[1]-trp[3]*x1+trp[2]*y1
        dx=x2-x2r
        dy=y2-y2r
        dr=sqrt(dx**2+dy**2)
        if not pp[2]:
            sumn+=dr**2
        sumnall+=dr**2
        print(f'{x1:9.2f}   {y1:9.2f}   {x2:9.2f}   {y2:9.2f}   {x2r:9.2f}   \
{y2r:9.2f}   {dx:9.2f}   {dy:9.2f}   {dr:9.2f}')
    print(' squared mean: ',sqrt(sumn/(2*(len(ptl1)-len(olpt)))))
    print('with outliers: ',sqrt(sumnall/(2*len(ptl1))))

def randtrp(offs=(-1000.0,-1000.0,1000.0,1000.0),
            scaleiv=(0.999, 1.001), rotiv=(-pi, pi)):
    """Create a random Helmert transformation"""
    scale=random.uniform(scaleiv[0], scaleiv[1])
    rot=random.uniform(rotiv[0], rotiv[1])
    a=scale*cos(rot)
    b=scale*sin(rot)
    dx=random.uniform(offs[0], offs[2])
    dy=random.uniform(offs[1], offs[3])
    return (dx, dy, a, b)

def randptp(ptn, trp, stderr=0.07071, pt1win=(0.0,0.0,1000.0,1000.0),
            oln=0, oldst=(1.0,10.0)):
    """Create a list of ptn random point in pt1win area (ptl1),
the transformed point by a Helmert transformation of trp (ptl2),
and the random modified points by stderr normal variate error (ptl2e).
oln points will be outlier, whose indices are in olpti.
The function returns (ptl1, ptl2, ptl2e, olpti).
"""
    ptl1=[(random.uniform(pt1win[0], pt1win[2]),
           random.uniform(pt1win[1], pt1win[3])) for i in range(ptn)]
    ptl2e=[(trp[0]+trp[2]*pt[0]+trp[3]*pt[1],
            trp[0]-trp[3]*pt[0]+trp[2]*pt[1]) for pt in ptl1]
    ptl2=[(pt[0]+random.normalvariate(0,stderr),
           pt[1]+random.normalvariate(0,stderr)) for pt in ptl2e]    
    olpti=random.sample(range(ptn),oln)
    for oi in olpti:
        ptx,pty=ptl2[oi]
        ptx+=random.uniform(10.0,50.0)*random.choice([-1,1])
        pty+=random.uniform(10.0,50.0)*random.choice([-1,1])
        ptl2[oi]=(ptx, pty)
    return (ptl1, ptl2, ptl2e, olpti)

def trpstat(ptl1, ptl2, olpt=[]):
    """Calculate the squared mean and the mean of the absolute values
of the difference between the points of ptl1 and ptl2 lists.
The calculation skipts points, whose indices are in the olpt.
The function returns (squared mean, absolute mean) in a tuple."""
    nn=len(ptl1)
    n=nn-len(olpt)
    sumr2=0.0
    sumra=0.0
    for ptp in zip(ptl1, ptl2, range(nn)):
        if not ptp[2] in olpt:
            dx=ptp[1][0]-ptp[0][0]
            dy=ptp[1][1]-ptp[0][1]
            dr2=dx**2+dy**2
            sumr2+=dr2
            sumra+=sqrt(dr2)
    return (sqrt(sumr2/n), sumra/n)
