import numpy as np
from Useful import hartree2wn,Ang2Bohr

def Potentials(r,R,A,b,ra,rb,d,em0a,dm1,ep0a,ep1a,ep2a):
    return A*np.exp(-b*(R-(r-ra)*d))-\
        (em0a+em0a/dm1*(r-ra))*(R-rb)**-3*\
        (1+(ep0a+ep1a*(r-ra)+ep2a/2*(r-ra)**2)*(R-rb)**-3)

def c0R(R,A,b,ra,rb,d,em0a,dm1,ep0a,ep1a,ep2a):
    return A*np.exp(-b*(R))-em0a*(R-rb)**-3*(1+ep0a*(R-rb)**-3)
def c0dR(R,A,b,ra,rb,d,em0a,dm1,ep0a,ep1a,ep2a):
    return -A*b * np.exp(-b * (R)) + 3*em0a * (R - rb) ** -4 * (1 + 2*ep0a * (R - rb) ** -3)
def c0dR2(R,A,b,ra,rb,d,em0a,dm1,ep0a,ep1a,ep2a):
    return A*b**2 * np.exp(-b * (R)) - 12*em0a * (R - rb) ** -5 * (1 + 7/2*ep0a * (R - rb) ** -3)
def c0dR3(R,A,b,ra,rb,d,em0a,dm1,ep0a,ep1a,ep2a):
    return -A*b**3 * np.exp(-b * (R)) + 60*em0a * (R - rb) ** -6 * (1 + 28/5*ep0a * (R - rb) ** -3)
def c0dR4(R,A,b,ra,rb,d,em0a,dm1,ep0a,ep1a,ep2a):
    return +A*b**4 * np.exp(-b * (R)) - 360*em0a * (R - rb) ** -7 * (1 + 6*7*8*9/360*ep0a* (R - rb) ** -3)
def c0Relec(R,A,b,ra,rb,d,em0a,dm1,ep0a,ep1a,ep2a):
    return -A*b*np.exp(-b*R)
def c0Rdipole(R,A,b,ra,rb,d,em0a,dm1,ep0a,ep1a,ep2a):
    return 3*em0a*(R-rb)**-4
def c0Rpolar(R,A,b,ra,rb,d,em0a,dm1,ep0a,ep1a,ep2a):
    return 6*em0a*ep1a*(R-rb)**-7

def c1R(R,A,b,ra,rb,d,em0a,dm1,ep0a,ep1a,ep2a):
    return A*b*d*np.exp(-b*(R))-(em0a/dm1)*(R-rb)**-3*(1+(ep0a+ep1a*dm1)*(R-rb)**-3)
def c1dR(R,A,b,ra,rb,d,em0a,dm1,ep0a,ep1a,ep2a):
    return -A*b**2*d*np.exp(-b*(R))+3*(em0a/dm1)*(R-rb)**-4*(1+2*(ep0a+ep1a*dm1)*(R-rb)**-3)
def c1dR2(R,A,b,ra,rb,d,em0a,dm1,ep0a,ep1a,ep2a):
    return A*b**3*d*np.exp(-b*(R))-12*(em0a/dm1)*(R-rb)**-5*(1+6*7/12*(ep0a+ep1a*dm1)*(R-rb)**-3)
def c1dR3(R,A,b,ra,rb,d,em0a,dm1,ep0a,ep1a,ep2a):
    return -A*b**4*d*np.exp(-b*(R))+60*(em0a/dm1)*(R-rb)**-6*(1+6*7*8/60*(ep0a+ep1a*dm1)*(R-rb)**-3)
def c1dR4(R,A,b,ra,rb,d,em0a,dm1,ep0a,ep1a,ep2a):
    return +A*b**5*d*np.exp(-b*(R))-360*(em0a/dm1)*(R-rb)**-7*(1+6*7*8*9/360*(ep0a+ep1a*dm1)*(R-rb)**-3)
def c1dR5(R,A,b,ra,rb,d,em0a,dm1,ep0a,ep1a,ep2a):
    return -A*b**6*d*np.exp(-b*(R))+360*7*(em0a/dm1)*(R-rb)**-8*(1+6*8*9*10/360*(ep0a+ep1a*dm1)*(R-rb)**-3)

def c2R(R,A,b,ra,rb,d,em0a,dm1,ep0a,ep1a,ep2a):
    return A*b**2*d**2*np.exp(-b*(R))-(2*ep1a*em0a/dm1+ep2a*em0a)*(R-rb)**-6
def c2dR(R,A,b,ra,rb,d,em0a,dm1,ep0a,ep1a,ep2a):
    return -A*b**3*d**2*np.exp(-b*(R))+6*(2*ep1a*em0a/dm1+ep2a*em0a)*(R-rb)**-7
def c2dR2(R,A,b,ra,rb,d,em0a,dm1,ep0a,ep1a,ep2a):
    return A*b**4*d**2*np.exp(-b*(R))-6*7*(2*ep1a*em0a/dm1+ep2a*em0a)*(R-rb)**-8
def c2dR3(R,A,b,ra,rb,d,em0a,dm1,ep0a,ep1a,ep2a):
    return -A*b**5*d**2*np.exp(-b*(R))+6*7*8*(2*ep1a*em0a/dm1+ep2a*em0a)*(R-rb)**-9
def c2dR4(R,A,b,ra,rb,d,em0a,dm1,ep0a,ep1a,ep2a):
    return +A*b**6*d**2*np.exp(-b*(R))-6*7*8*9*(2*ep1a*em0a/dm1+ep2a*em0a)*(R-rb)**-10



def Rr(R,ka,a,ps):
    return (-ka-c2R(R,*ps)+np.sqrt((+ka+c2R(R,*ps))**2+6*a*ka*c1R(R,*ps)))/(-3*a*ka)

def Rdr(R,ka,a,ps):
    r=Rr(R,ka,a,ps)
    z=np.exp(-a*r)
    return -(ka*(2*z**2-z)+c2R(R,*ps))/(c1dR(R,*ps)+c2dR(R,*ps)*r)

def Rdr2(R,ka,a,ps):
    r=Rr(R,ka,a,ps)
    dr=Rdr(R,ka,a,ps)
    z=np.exp(-a*r)
    return (a*ka*(4*z**2-z)-dr*(2*c2dR(R,*ps)+dr*(c1dR2(R,*ps)+c2dR2(R,*ps)*r)))/(c1dR(R,*ps)+c2dR(R,*ps)*r)

def Rdr3(R,ka,a,ps):
    r=Rr(R,ka,a,ps)
    dr=Rdr(R,ka,a,ps)
    dr2=Rdr2(R,ka,a,ps)
    z=np.exp(-a*r)
    return (-a**2*ka*(8*z**2-z)-c1dR3(R,*ps)*dr**3-3*c1dR2(R,*ps)*dr2*dr)/c1dR(R,*ps)
def Rdr4(R,ka,a,ps):
    r=Rr(R,ka,a,ps)
    dr=Rdr(R,ka,a,ps)
    z=np.exp(-a*r)
    return (a**3*ka*(16*z**2-z)-c1dR4(R,*ps)*Rdr(R,r,ka,a,ps)**4-6*c1dR3(R,*ps)*Rdr2(R,r,ka,a,ps)*Rdr(R,r,ka,a,ps)**2
            -c1dR2(R,*ps)*(4*Rdr3(R,r,ka,a,ps)*Rdr(R,r,ka,a,ps)+3*Rdr2(R,r,ka,a,ps)**2))/c1dR(R,*ps)

def c0dr(R,ka,a,ps):
    return c0dR(R,*ps)*Rdr(R,ka,a,ps)
def c1dr(R,ka,a,ps):
    return c1dR(R,*ps)*Rdr(R,ka,a,ps)
def c2dr(R,ka,a,ps):
    return c2dR(R,*ps)*Rdr(R,ka,a,ps)
def c0dr2(R,ka,a,ps):
    return c0dR2(R,*ps)*Rdr(R,ka,a,ps)**2+c0dR(R,*ps)*Rdr2(R,ka,a,ps)

def c0dr3(R,ka,a,ps):
    return c0dR3(R,*ps)*Rdr(R,ka,a,ps)**3+c0dR2(R,*ps)*Rdr2(R,ka,a,ps)*Rdr(R,ka,a,ps)+c0dR(R,*ps)*Rdr3(R,ka,a,ps)

def Evsr(drs,Rvsr,a,d,ps):
    z=np.exp(-a*drs)
    return -d*(z-1)**2-c0R(Rvsr,*ps)-c1R(Rvsr,*ps)*drs-c2R(Rvsr,*ps)*drs**2/2

def wvsr(drs,Rvsr,a,d,ps):
    ka=2*d*a*a
    wa=np.sqrt(ka/1830/Ang2Bohr**2)*hartree2wn
    z=np.exp(-a*drs)
    return wa/2*(1-(2*z*z-z)-c2R(Rvsr,*ps)/ka)#-c3R(Rvsr,*ps)/ka*drs)