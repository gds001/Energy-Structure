import numpy as np
from numpy import cos,sin, pi

def norm2(list): return np.sqrt(list[0]*list[0]+list[1]*list[1])
def norm22(list): return norm2(list)*norm2(list)
def redmass(m1,m2): return m1*m2/(m1+m2)
h=6.62607004E-34 # Js
c=2.99792458E8 # m/s
k=1.38064852E-23 # J/K
Cal2J = 4.184 # J/cal
J2Hartree=2.294E17 #Hatree/J
Ang2Bohr=1.88973 #Bohr/Ang
mDynepA_to_EhpA2 = 1/15.569141*Ang2Bohr**2
Deg2Rad=np.pi/180 #radian/degree
amu2kg=1.66054E-27 # kg/amu
amu2me=1822.89 # me/amu
s2au=2.418884236E-17 # s
Na=6.02214E23 #molecules/mol
hbar=h/2/pi #Js
e = 1.602176634E-19 #C
ao = 5.29177210903E-11 #m
Rh=1.09677583E7 #1/m
eo=8.85418781762039E-12 #C^2/Nm^2 =
c=2.998E8 #m/s
po=101325 #Pa
J2Ev=6.624E18 # eV/J
hbarev=hbar*J2Ev
hartree2wn=219474.63
oscstrength2kmmol=5.33049e6*4.70165e-7##?????? Check if it's right
atm2pascal=101325
aP2wn=s2au/c/100
#amu
mH=1.00784
mD=2.014
mT=3.0160492
mB=10.811
mC=12.0107
mN=14.0067
mO=15.999
mF=18.998403
mCl=35.453
mBr=79.904

def R2(yt,yr):
    sst=0
    m=np.mean(yt)
    for i in yt: sst+=(i-m)**2
    ssr=0
    for i in range(len(yt)):
        ssr+=(yt[i]-yr[i])**2
    return 1-ssr/sst

def chi2(w,e,f):
    num=0
    for i in range(len(e)):
        num+=w[i]*(e[i]-f[i])**2
    return num/sum(w)

def RMSD(e,f):
    result=(e-f)**2
    return np.sqrt(np.sum(result,where=~np.isnan(result))/len(result[~np.isnan(result)]))
def RMSD2(e,f1,f2):
    result=0
    for i in range(len(e)):
        if np.isnan(f1[i]): continue
        result+=np.amin([(e[i]-f1[i])**2,(e[i]-f2[i])**2])
    return np.sqrt(result/len(e))

def MAE(e,f):
    result=np.abs(e-f)
    return np.sum(result,where=~np.isnan(result))/len(result[~np.isnan(result)])
def MAE2(e,f1,f2):
    result=0
    for i in range(len(e)):
        if np.isnan(f1[i]): continue
        result+=np.amin([np.abs(e[i]-f1[i]),np.abs(e[i]-f2[i])])
    return result/len(e)

def wRMSD(e,f,ws=None):
    #if ws==None: ws=np.ones(len(e))
    result=(e-f)**2*ws
    denom=np.sum(ws,where=~np.isnan(f))
    return np.sqrt(np.sum(result,where=~np.isnan(f))/denom)
def wRMSD2(e,f1,f2,ws):
    result=0
    for i in range(len(e)):
        if np.isnan(f1[i]): continue
        result+=np.amin([(e[i]-f1[i])**2,(e[i]-f2[i])**2])*ws[i]
    return np.sqrt(result/np.sum(ws))

def wMAE(e,f,ws):
    result=np.abs(e-f)*ws
    denom=np.sum(ws,where=~np.isnan(f))
    return np.sum(result,where=~np.isnan(f))/denom

def wMAE2(e,f1,f2,ws):
    result=0
    for i in range(len(e)):
        if np.isnan(f1[i]): continue
        result+=np.amin([np.abs(e[i]-f1[i]),np.abs(e[i]-f2[i])])*ws[i]
    return result/np.sum(ws)

def m1ton(n):
    if n%2==0: return 1
    return -1

def vec2mat(vec):
    n=int(np.sqrt(len(vec)))
    return np.reshape(vec,(n,n))
def mat2vec(mat):
    n=len(mat[0])
    return np.reshape(mat,(n*n,1))

def CD1(xs,ys):
    dx=xs[1]-xs[0]
    for i in range(2,len(xs)):
        temp=xs[i]-xs[i-1]
        if np.abs(temp-dx)>1e-12:
            print("Non-even Grid", dx,temp)
            return None
    if len(xs)==3:
        return (ys[2]/2-ys[0]/2)/dx
    if len(xs)==5:
        return (ys[0]/12-ys[1]*2/3+ys[3]*2/3-ys[4]/12)/dx
    if len(xs)==7:
        return (ys[6]/60-ys[5]*3/20+ys[4]*3/4-ys[2]*3/4+ys[1]*3/20-ys[0]/60)/dx

def CD2(xs,ys):
    dx=xs[1]-xs[0]
    for i in range(2,len(xs)):
        temp=xs[i]-xs[i-1]
        if np.abs(temp-dx)>1e-12:
            print("Non-even Grid", dx,temp)
            return None
    if len(xs)==3:
        return (ys[2]-2*ys[1]+ys[0])/dx**2
    if len(xs)==5:
        return (-ys[0]/12+ys[1]*4/3-ys[2]*5/2+ys[3]*4/3-ys[4]/12)/dx**2
    if len(xs)==7:
        return (ys[0]/90-ys[1]*3/20+ys[2]*3/2-ys[3]*49/18+ys[4]*3/2-ys[5]*3/20+ys[6]/90)/dx**2

def CD3(xs,ys):
    dx=xs[1]-xs[0]
    for i in range(2,len(xs)):
        temp=xs[i]-xs[i-1]
        if np.abs(temp-dx)>1e-12:
            print("Non-even Grid", dx,temp)
            return None
    if len(xs)==5:
        return (-ys[0]/2+ys[1]-ys[3]+ys[4]/2)/dx**3
    if len(xs)==7:
        return (ys[0]/8-ys[1]+ys[2]*13/8-ys[4]*13/8+ys[5]-ys[6]/8)/dx**3

def CD4(xs,ys):
    dx=xs[1]-xs[0]
    for i in range(2,len(xs)):
        temp=xs[i]-xs[i-1]
        if np.abs(temp-dx)>1e-12:
            print("Non-even Grid", dx,temp)
            return None
    if len(xs)==5:
        return (ys[0]-4*ys[1]+6*ys[2]-4*ys[3]+ys[4])/dx**4
    if len(xs)==7:
        return (-ys[0]/6+2*ys[1]-13/2*ys[2]+28/3*ys[3]-13/2*ys[4]+2*ys[5]-ys[6]/6)/dx**4

def chain1(dydu,dudx):
    return dydu*dudx
def chain2(dydu,d2ydu2,dudx,d2udx2):
    return d2ydu2*dudx**2+dydu*d2udx2
def chain3(dydu,d2ydu2,d3ydu3,dudx,d2udx2,d3udx3):
    return d3ydu3*dudx**3+3*d2ydu2*dudx*d2udx2+dydu*d3udx3
def chain4(dydu,d2ydu2,d3ydu3,d4ydu4,dudx,d2udx2,d3udx3,d4udx4):
    return d4ydu4*dudx**4+6*d3ydu3*dudx**2*d2udx2+d2ydu2*(4*dudx*d3udx3+3*d2udx2**2)+dydu*d4udx4

if __name__ == "__main__":
    print(redmass(mH,mO))