import numpy as np
from Useful import hartree2wn

def Potentials(r,R,A,b,ra,rb,d,em0a,em1a,em2a,em3a,ep0a,ep1a,ep2a,ep3a):
    return A*np.exp(-b*(R-(r-ra)*d))-\
        (em0a+em1a*(r-ra)+em2a/2*(r-ra)**2+em3a/6*(r-ra)**3)*(R-rb)**-3*\
        (1+(ep0a+ep1a*(r-ra)+ep2a/2*(r-ra)**2+ep3a/6*(r-ra)**3)*(R-rb)**-3)

def c0R(R,A,b,ra,rb,d,em0a,em1a,em2a,em3a,ep0a,ep1a,ep2a,ep3a):
    return A*np.exp(-b*(R))-(em0a)*(R-rb)**-3*(1+(ep0a)*(R-rb)**-3)
def c0dR(R,A,b,ra,rb,d,em0a,em1a,em2a,em3a,ep0a,ep1a,ep2a,ep3a):
    return -A*b * np.exp(-b * (R)) + 3*(em0a) * (R - rb) ** -4 * (1 + 2*(ep0a) * (R - rb) ** -3)
def c0dR2(R,A,b,ra,rb,d,em0a,em1a,em2a,em3a,ep0a,ep1a,ep2a,ep3a):
    return A*b**2 * np.exp(-b * (R)) - 12*(em0a) * (R - rb) ** -5 * (1 + 7/2*(ep0a) * (R - rb) ** -3)
def c0dR3(R,A,b,ra,rb,d,em0a,em1a,em2a,em3a,ep0a,ep1a,ep2a,ep3a):
    return -A*b**3 * np.exp(-b * (R)) + 60*(em0a) * (R - rb) ** -6 * (1 + 28/5*(ep0a) * (R - rb) ** -3)
def c0dR4(R,A,b,ra,rb,d,em0a,em1a,em2a,em3a,ep0a,ep1a,ep2a,ep3a):
    return +A*b**4 * np.exp(-b * (R)) - 360*(em0a) * (R - rb) ** -7 * (1 + 6*7*8*9/360*(ep0a) * (R - rb) ** -3)


def c1R(R,A,b,ra,rb,d,em0a,em1a,em2a,em3a,ep0a,ep1a,ep2a,ep3a):
    return A*b*d*np.exp(-b*(R))-(em1a)*(R-rb)**-3*(1+(ep0a+ep1a*em0a/em1a)*(R-rb)**-3)
def c1dR(R,A,b,ra,rb,d,em0a,em1a,em2a,em3a,ep0a,ep1a,ep2a,ep3a):
    return -A*b**2*d*np.exp(-b*(R))+3*(em1a)*(R-rb)**-4*(1+2*(ep0a+ep1a*em0a/em1a)*(R-rb)**-3)
def c1dR2(R,A,b,ra,rb,d,em0a,em1a,em2a,em3a,ep0a,ep1a,ep2a,ep3a):
    return A*b**3*d*np.exp(-b*(R))-12*(em1a)*(R-rb)**-5*(1+6*7/12*(ep0a+ep1a*em0a/em1a)*(R-rb)**-3)
def c1dR3(R,A,b,ra,rb,d,em0a,em1a,em2a,em3a,ep0a,ep1a,ep2a,ep3a):
    return -A*b**4*d*np.exp(-b*(R))+60*(em1a)*(R-rb)**-6*(1+6*7*8/60*(ep0a+ep1a*em0a/em1a)*(R-rb)**-3)
def c1dR4(R,A,b,ra,rb,d,em0a,em1a,em2a,em3a,ep0a,ep1a,ep2a,ep3a):
    return +A*b**5*d*np.exp(-b*(R))-360*(em1a)*(R-rb)**-7*(1+6*7*8*9/360*(ep0a+ep1a*em0a/em1a)*(R-rb)**-3)
def c1dR5(R,A,b,ra,rb,d,em0a,em1a,em2a,em3a,ep0a,ep1a,ep2a,ep3a):
    return -A*b**6*d*np.exp(-b*(R))+360*7*(em1a)*(R-rb)**-8*(1+6*8*9*10/360*(ep0a+ep1a*em0a/em1a)*(R-rb)**-3)

def c2R(R,A,b,ra,rb,d,em0a,em1a,em2a,em3a,ep0a,ep1a,ep2a,ep3a):
    return A*b**2*d**2*np.exp(-b*(R))-(em2a)*(R-rb)**-3*(1+(ep0a+2*ep1a*em1a/em2a+ep2a*em0a/em2a)*(R-rb)**-3)
def c2dR(R,A,b,ra,rb,d,em0a,em1a,em2a,em3a,ep0a,ep1a,ep2a,ep3a):
    return -A*b**3*d**2*np.exp(-b*(R))+3*em2a*(R-rb)**-4*(1+2*(ep0a+2*ep1a*em1a/em2a+ep2a*em0a/em2a)*(R-rb)**-3)
def c2dR2(R,A,b,ra,rb,d,em0a,em1a,em2a,em3a,ep0a,ep1a,ep2a,ep3a):
    return A*b**4*d**2*np.exp(-b*(R))-12*(em2a)*(R-rb)**-5*(1+6*7/12*(ep0a+2*ep1a*em1a/em2a+ep2a*em0a/em2a)*(R-rb)**-3)
def c2dR3(R,A,b,ra,rb,d,em0a,em1a,em2a,em3a,ep0a,ep1a,ep2a,ep3a):
    return -A*b**5*d**2*np.exp(-b*(R))+60*(em2a)*(R-rb)**-6*(1+6*7*8/60*(ep0a+2*ep1a*em1a/em2a+ep2a*em0a/em2a)*(R-rb)**-3)
def c2dR4(R,A,b,ra,rb,d,em0a,em1a,em2a,em3a,ep0a,ep1a,ep2a,ep3a):
    return +A*b**6*d**2*np.exp(-b*(R))-360*(em2a)*(R-rb)**-7*(1+6*7*8*9/360*(ep0a+2*ep1a*em1a/em2a+ep2a*em0a/em2a)*(R-rb)**-3)

def c3R(R,A,b,ra,rb,d,em0a,em1a,em2a,em3a,ep0a,ep1a,ep2a,ep3a):
    return A*b**3*d**3*np.exp(-b*(R))-(em3a)*(R-rb)**-3*(1+(ep0a)*(R-rb)**-3)

def Evsr(drs,Rvsr,a,d,ps):
    z=np.exp(-a*drs)
    return -d*(z-1)**2-c0R(Rvsr,*ps)-c1R(Rvsr,*ps)*drs-c2R(Rvsr,*ps)*drs**2/2-c3R(Rvsr,*ps)*drs**3/6

def wvsr(drs,Rvsr,a,d,ps):
    ka=2*d*a*a
    wa=np.sqrt(ka/1830)*hartree2wn
    z=np.exp(-a*drs)
    return wa/2*(1-(2*z*z-z)-c2R(Rvsr,*ps)/ka-c3R(Rvsr,*ps)/ka*drs)