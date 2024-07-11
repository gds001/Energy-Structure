import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import matplotlib as mpl
import pandas as pd
import json
import scipy.interpolate as inter
from Useful import CD1,CD2,CD3
import sys
from SemiRigidScans_Simpler_Relax_rb_fix_dpr.ExpansionAndDerivatives import Potentials
print("hf_hf")
np.set_printoptions(suppress=True,precision=4,linewidth=1000)
np.seterr(all='ignore', divide='ignore', over='ignore', under='ignore', invalid='ignore')
data=np.genfromtxt("hf_hf_linear_mp2_avtz_rahs_rab.csv",delimiter=',',dtype=np.double).T
don=np.genfromtxt("../hf/hf_mp2_avtz_scan.csv",delimiter=',',dtype=np.double).T
acc=np.genfromtxt("../hf/hf_mp2_avtz_scan.csv",delimiter=',',dtype=np.double).T


def c0f(R,A,b,ra,rb,g,d,em0a,em1a,em2a,n):
    return A*np.exp(-b*(R-ra*d))-(em0a)*(R-ra*g-rb)**(-n)
def c1f(R,A,b,ra,rb,g,d,em0a,em1a,em2a,n):
    return A*b*d*np.exp(-b*(R-ra*d))\
        -(em0a)*(R-ra*g-rb)**(-n-1)*n*g-(em1a)*(R-ra*g-rb)**(-n)
def c2f(R,A,b,ra,rb,g,d,em0a,em1a,em2a,n):
    return A*b*b*d*d*np.exp(-b*(R-ra*d))\
        -(em0a)*(R-ra*g-rb)**(-n-2)*n*(n+1)*g*g\
        -(em1a)*(R-ra*g-rb)**(-n-1)*n*g*2\
        -(em2a)*(R-ra*g-rb)**(-n)
def reR(R,a,da,A,b,ra,rb,g,d,em0a,em1a,em2a,n):
    return -1/a*np.log(1/2+np.sqrt(1/4+A*b*d/(2*a*da)*np.exp(-b*(R-d*ra))-em1a/(2*a*da*(R-g*ra-rb)**n)-n*g*em0a/(2*a*da*(R-g*ra-rb)**(n+1))))

def MSD(params,x1,x2,y,func):
    return np.sqrt(np.sum((func(x1,x2,*params)-y)**2)/len(y))
n=np.argmin(don[1])
ra=don[0][n]
print(ra)
#re=dimer2d[0][n]
msl=[n-3,n-2,n-1,n,n+1,n+2,n+3]
kA=CD2(don[0][msl],don[1][msl])
d3=CD3(don[0][msl],don[1][msl])
print(kA,d3)

a=d3/kA/-3
da=kA/2/a**2
print(a,da)
plt.scatter(don[0],don[1]-np.amin(don[1]))
plt.plot(don[0],da*(np.exp(-a*(don[0]-ra))-1)**2)
plt.close()

data[3]-=np.amin(don[1])+np.amin(acc[1])
m=np.argmin(data[3])
rs=data[0][data[1]==data[1][m]]
Rs=data[1][data[0]==data[0][m]]
ers=data[3][data[0]==data[0][m]]
data=data[:,data[1,:]>Rs[np.argmax(ers<0)]]
data=data[:,data[0,:]>rs[1]]
data=data[:,data[0,:]<rs[-2]]
m=np.argmin(data[3])
rs=data[0][data[1]==data[1][m]]
Rs=data[1][data[0]==data[0][m]]
ers=data[3][data[0]==data[0][m]]
es=data[3]-da*(np.exp(-a*(data[0]-ra))-1)**2
Es=np.reshape(es,(len(Rs),len(rs)))
s=1e-12
while True:
    try:
        file=open("Fits_HF_HF_Polar.json",'r')
        dic = json.loads(file.read())
    except:
        print("Starting Fresh")
        basic = {'x': [212.43, 3.4, ra, 0, 0.42, 0.11, 0.6593, 2.3289, 44.66], 'fun': 1}
        dic = {
            'fitdda': basic,
            'fitddda': basic,
        }
    ni = 10_000
    print(ni, "hf_hf")
    t = dic['fitdda']['fun']
    temp = np.array(dic['fitdda']['x'], dtype=np.double)
    # if fitdd.fun < t:
    #    t = fitdd.fun
    #    temp = fitdd.x
    # if fitda.fun < t:
    #    t = fitda.fun
    #    temp = fitda.x
    fitdda = opt.basinhopping(MSD, x0=temp,
                              minimizer_kwargs={'args': (data[0], data[1], es, Potentials),
                                                'bounds': ((0, np.inf), (0, 10),
                                                           (ra - s, ra + s), (-1,1),
                                                           (0, 1),
                                                           (0, np.inf), (0.6035 - s, 0.6035 + s),
                                                           (0 - s, 0 + s),
                                                           (0, np.inf),
                                                           ),
                                                'method': "Nelder-Mead"
                                                }, niter=ni)
    print('fitdda', fitdda.fun * 627.5, repr(fitdda.x))

    t = dic['fitddda']['fun']
    temp = np.array(dic['fitddda']['x'], dtype=np.double)
    # if fitddd.fun < t:
    #    t = fitddd.fun
    #    temp = fitddd.x
    # if fitdda.fun < t:
    #    t = fitdda.fun
    #    temp = fitdda.x
    fitddda = opt.basinhopping(MSD, x0=temp,
                               minimizer_kwargs={'args': (data[0], data[1], es, Potentials),
                                                 'bounds': ((0, np.inf), (0, 10),
                                                            (ra - s, ra + s), (-1,1),
                                                            (0, 1),
                                                            (0, np.inf), (0.6035 - s, 0.6035 + s),
                                                            (2.2783 - s, 2.2783 + s),
                                                            (0, np.inf),
                                                            ),
                                                 # 'method': "Nelder-Mead"
                                                 }, niter=ni)
    print('fitddda', fitddda.fun * 627.5, repr(fitddda.x))

    print(
        dic['fitdda']['fun'] - fitdda.fun,
        dic['fitddda']['fun'] - fitddda.fun, )
    improve = np.abs(
        np.abs(dic['fitdda']['fun'] - fitdda.fun) + \
        np.abs(dic['fitddda']['fun'] - fitddda.fun))
    if improve < 1e-12:
        print("Converged")

    for i in fitdda:
        try:

            fitdda[i] = fitdda[i].tolist()
            fitddda[i] = fitddda[i].tolist()
        except:
            continue
    dict = {
        'fitdda': {'fun': fitdda.fun, 'x': fitdda.x},
        'fitddda': {'fun': fitddda.fun, 'x': fitddda.x}, }
    str = json.dumps(dict, indent=4)
    file=open("Fits_HF_HF_Polar.json",'w')
    file.write(str)
    file.close()
    if improve<1e-12: break
sys.exit()
