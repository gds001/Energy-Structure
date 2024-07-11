import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import matplotlib as mpl
import pandas as pd
import json
import scipy.interpolate as inter
from Useful import CD1,CD2,CD3
from SemiRigidScans_Simpler_Relax_rb.ExpansionAndDerivatives import Potentials
import sys
print("h2o_nh3")
np.set_printoptions(suppress=True,precision=4,linewidth=1000)
np.seterr(all='ignore', divide='ignore', over='ignore', under='ignore', invalid='ignore')
data=np.genfromtxt("h2o_nh3_linear_mp2_avtz_rahs_rab.csv",delimiter=',',dtype=np.double).T
don=np.genfromtxt("../h2o/h2o_mp2_avtz_scan.csv",delimiter=',',dtype=np.double).T
acc=np.genfromtxt("../nh3/nh3_mp2_avtz_scan.csv",delimiter=',',dtype=np.double).T

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
        file=open("Fits_H2O_NH3_Polar.json",'r')
        dic = json.loads(file.read())
    except:
        print("Starting Fresh")
        basic = {'x': [300, 3, ra, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 'fun': 1}
        dic = {
            'fitda': basic,
            'fitdda': basic,
            'fitdaa': basic,
            'fitddda': basic,
            'fitddaa': basic,
            'fitdaaa': basic,
            'fitdddda': basic,
            'fitdddaa': basic,
            'fitddaaa': basic,
            'fitdaaaa': basic,
            'fitddddaa': basic,
            'fitdddaaa': basic,
            'fitddaaaa': basic,
            'fitddddaaa': basic,
            'fitdddaaaa': basic,
            'fitddddaaaa': basic,
        }
    ni = 1_000
    print(ni, "h2o_nh3")
    # if fitd.fun < dic['fitda']['fun']:
    #    temp = fitd.x
    # else:
    temp = np.array(dic['fitda']['x'], dtype=np.double)
    fitda = opt.basinhopping(MSD, x0=temp,
                             minimizer_kwargs={'args': (data[0], data[1], es, Potentials),
                                               'bounds': ((0, np.inf), (0, 10),
                                                          (ra - s, ra + s), (-np.inf, np.inf),
                                                          (0, 1),
                                                          (0, np.inf), (0 - s, 0 + s), (0 - s, 0 + s), (0 - s, 0 + s),
                                                          (0, np.inf), (0 - s, 0 + s), (0 - s, 0 + s), (0 - s, 0 + s),
                                                          ),
                                               'method': "Nelder-Mead"
                                               }, niter=ni)
    print('fitda', fitda.fun * 627.5, repr(fitda.x))

    t = dic['fitdda']['fun']
    temp = np.array(dic['fitdda']['x'], dtype=np.double)
    if fitda.fun < t:
        t = fitda.fun
        temp = fitda.x
    fitdda = opt.basinhopping(MSD, x0=temp,
                              minimizer_kwargs={'args': (data[0], data[1], es, Potentials),
                                                'bounds': ((0, np.inf), (0, 10),
                                                           (ra - s, ra + s), (-np.inf, np.inf),
                                                           (0, 1),
                                                           (0, np.inf), (-np.inf, np.inf), (0 - s, 0 + s),
                                                           (0 - s, 0 + s),
                                                           (0, np.inf), (0 - s, 0 + s), (0 - s, 0 + s),
                                                           (0 - s, 0 + s),
                                                           ),
                                                'method': "Nelder-Mead"
                                                }, niter=ni)
    print('fitdda', fitdda.fun * 627.5, repr(fitdda.x))
    if fitda.fun < dic['fitdaa']['fun']:
        temp = fitda.x
    else:
        temp = np.array(dic['fitdaa']['x'], dtype=np.double)
    fitdaa = opt.basinhopping(MSD, x0=temp,
                              minimizer_kwargs={'args': (data[0], data[1], es, Potentials),
                                                'bounds': ((0, np.inf), (0, 10),
                                                           (ra - s, ra + s), (-np.inf, np.inf),
                                                           (0, 1),
                                                           (0, np.inf), (0 - s, 0 + s), (0 - s, 0 + s),
                                                           (0 - s, 0 + s),
                                                           (0, np.inf), (-np.inf, np.inf), (0 - s, 0 + s),
                                                           (0 - s, 0 + s),
                                                           ),
                                                # 'method': "Nelder-Mead"
                                                }, niter=ni)
    print('fitdaa', fitdaa.fun * 627.5, repr(fitdaa.x))
    t = dic['fitddda']['fun']
    temp = np.array(dic['fitddda']['x'], dtype=np.double)
    if fitdda.fun < t:
        t = fitdda.fun
        temp = fitdda.x
    fitddda = opt.basinhopping(MSD, x0=temp,
                               minimizer_kwargs={'args': (data[0], data[1], es, Potentials),
                                                 'bounds': ((0, np.inf), (0, 10),
                                                            (ra - s, ra + s), (-np.inf, np.inf),
                                                            (0, 1),
                                                            (0, np.inf), (-np.inf, np.inf), (-np.inf, np.inf),
                                                            (0 - s, 0 + s),
                                                            (0, np.inf), (0 - s, 0 + s), (0 - s, 0 + s),
                                                            (0 - s, 0 + s),
                                                            ),
                                                 # 'method': "Nelder-Mead"
                                                 }, niter=ni)
    print('fitddda', fitddda.fun * 627.5, repr(fitddda.x))
    t = dic['fitddaa']['fun']
    temp = np.array(dic['fitddaa']['x'], dtype=np.double)
    if fitdda.fun < t:
        t = fitdda.fun
        temp = fitdda.x
    if fitdaa.fun < t:
        t = fitdaa.fun
        temp = fitdaa.x
    fitddaa = opt.basinhopping(MSD, x0=temp,
                               minimizer_kwargs={'args': (data[0], data[1], es, Potentials),
                                                 'bounds': ((0, np.inf), (0, 10),
                                                            (ra - s, ra + s), (-np.inf, np.inf),
                                                            (0, 1),
                                                            (0, np.inf), (-np.inf, np.inf), (0 - s, 0 + s),
                                                            (0 - s, 0 + s),
                                                            (0, np.inf), (-np.inf, np.inf), (0 - s, 0 + s),
                                                            (0 - s, 0 + s),
                                                            ),
                                                 # 'method': "Nelder-Mead"
                                                 }, niter=ni)
    print('fitddaa', fitddaa.fun * 627.5, repr(fitddaa.x))
    if fitdaa.fun < dic['fitdaaa']['fun']:
        temp = fitdaa.x
    else:
        temp = np.array(dic['fitdaaa']['x'], dtype=np.double)
    fitdaaa = opt.basinhopping(MSD, x0=temp,
                               minimizer_kwargs={'args': (data[0], data[1], es, Potentials),
                                                 'bounds': ((0, np.inf), (0, 10),
                                                            (ra - s, ra + s), (-np.inf, np.inf),
                                                            (0, 1),
                                                            (0, np.inf), (0 - s, 0 + s), (0 - s, 0 + s), (0 - s, 0 + s),
                                                            (0, np.inf), (-np.inf, np.inf), (-np.inf, np.inf),
                                                            (0 - s, 0 + s),
                                                            ),
                                                 'method': "Nelder-Mead"
                                                 }, niter=ni)
    print('fitdaaa', fitdaaa.fun * 627.5, repr(fitdaaa.x))
    t = dic['fitdddda']['fun']
    temp = np.array(dic['fitdddda']['x'], dtype=np.double)
    if fitddda.fun < t:
        t = fitddda.fun
        temp = fitddda.x
    fitdddda = opt.basinhopping(MSD, x0=temp,
                                minimizer_kwargs={'args': (data[0], data[1], es, Potentials),
                                                  'bounds': ((0, np.inf), (0, 10),
                                                             (ra - s, ra + s), (-np.inf, np.inf),
                                                             (0, 1),
                                                             (0, np.inf), (-np.inf, np.inf), (-np.inf, np.inf),
                                                             (-np.inf, np.inf),
                                                             (0, np.inf), (0 - s, 0 + s), (0 - s, 0 + s),
                                                             (0 - s, 0 + s),
                                                             ),
                                                  # 'method': "Nelder-Mead"
                                                  }, niter=ni)
    print('fitdddda', fitdddda.fun * 627.5, repr(fitdddda.x))
    t = dic['fitdddaa']['fun']
    temp = np.array(dic['fitdddaa']['x'], dtype=np.double)
    if fitddda.fun < t:
        t = fitddda.fun
        temp = fitddda.x
    if fitddaa.fun < t:
        t = fitddaa.fun
        temp = fitddaa.x
    fitdddaa = opt.basinhopping(MSD, x0=temp,
                                minimizer_kwargs={'args': (data[0], data[1], es, Potentials),
                                                  'bounds': ((0, np.inf), (0, 10),
                                                             (ra - s, ra + s), (-np.inf, np.inf),
                                                             (0, 1),
                                                             (0, np.inf), (-np.inf, np.inf), (-np.inf, np.inf),
                                                             (0 - s, 0 + s),
                                                             (0, np.inf), (-np.inf, np.inf), (0 - s, 0 + s),
                                                             (0 - s, 0 + s),
                                                             ),
                                                  # 'method': "Nelder-Mead"
                                                  }, niter=ni)
    print('fitdddaa', fitdddaa.fun * 627.5, repr(fitdddaa.x))
    t = dic['fitddaaa']['fun']
    temp = np.array(dic['fitddaaa']['x'], dtype=np.double)
    if fitddaa.fun < t:
        t = fitddaa.fun
        temp = fitddaa.x
    if fitdaaa.fun < t:
        t = fitdaaa.fun
        temp = fitdaaa.x
    fitddaaa = opt.basinhopping(MSD, x0=temp,
                                minimizer_kwargs={'args': (data[0], data[1], es, Potentials),
                                                  'bounds': ((0, np.inf), (0, 10),
                                                             (ra - s, ra + s), (-np.inf, np.inf),
                                                             (0, 1),
                                                             (0, np.inf), (-np.inf, np.inf), (0 - s, 0 + s),
                                                             (0 - s, 0 + s),
                                                             (0, np.inf), (-np.inf, np.inf), (-np.inf, np.inf),
                                                             (0 - s, 0 + s),
                                                             ),
                                                  # 'method': "Nelder-Mead"
                                                  }, niter=ni)
    print('fitdddaa', fitdddaa.fun * 627.5, repr(fitdddaa.x))
    if fitdaaa.fun < dic['fitdaaaa']['fun']:
        temp = fitdaaa.x
    else:
        temp = np.array(dic['fitdaaaa']['x'], dtype=np.double)
    fitdaaaa = opt.basinhopping(MSD, x0=temp,
                                minimizer_kwargs={'args': (data[0], data[1], es, Potentials),
                                                  'bounds': ((0, np.inf), (0, 10),
                                                             (ra - s, ra + s), (-np.inf, np.inf),
                                                             (0, 1),
                                                             (0, np.inf), (0 - s, 0 + s), (0 - s, 0 + s),
                                                             (0 - s, 0 + s),
                                                             (0, np.inf), (-np.inf, np.inf), (-np.inf, np.inf),
                                                             (-np.inf, np.inf),
                                                             ),
                                                  'method': "Nelder-Mead"
                                                  }, niter=ni)
    print('fitdaaaa', fitdaaaa.fun * 627.5, repr(fitdaaaa.x))

    t = dic['fitddddaa']['fun']
    temp = np.array(dic['fitddddaa']['x'], dtype=np.double)
    if fitdddda.fun < t:
        t = fitdddda.fun
        temp = fitdddda.x
    if fitdddaa.fun < t:
        t = fitdddaa.fun
        temp = fitdddaa.x
    fitddddaa = opt.basinhopping(MSD, x0=temp,
                                 minimizer_kwargs={'args': (data[0], data[1], es, Potentials),
                                                   'bounds': ((0, np.inf), (0, 10),
                                                              (ra - s, ra + s), (-np.inf, np.inf),
                                                              (0, 1),
                                                              (0, np.inf), (-np.inf, np.inf), (-np.inf, np.inf),
                                                              (-np.inf, np.inf),
                                                              (0, np.inf), (-np.inf, np.inf), (0 - s, 0 + s),
                                                              (0 - s, 0 + s),
                                                              ),
                                                   # 'method': "Nelder-Mead"
                                                   }, niter=ni)
    print('fitddddaa', fitddddaa.fun * 627.5, repr(fitddddaa.x))

    t = dic['fitdddaaa']['fun']
    temp = np.array(dic['fitdddaaa']['x'], dtype=np.double)
    if fitdddaa.fun < t:
        t = fitdddaa.fun
        temp = fitdddaa.x
    if fitddaaa.fun < t:
        t = fitddaaa.fun
        temp = fitddaaa.x
    fitdddaaa = opt.basinhopping(MSD, x0=temp,
                                 minimizer_kwargs={'args': (data[0], data[1], es, Potentials),
                                                   'bounds': ((0, np.inf), (0, 10),
                                                              (ra - s, ra + s), (-np.inf, np.inf),
                                                              (0, 1),
                                                              (0, np.inf), (-np.inf, np.inf), (-np.inf, np.inf),
                                                              (0 - s, 0 + s),
                                                              (0, np.inf), (-np.inf, np.inf), (-np.inf, np.inf),
                                                              (0 - s, 0 + s),
                                                              ),
                                                   # 'method': "Nelder-Mead"
                                                   }, niter=ni)
    print('fitdddaaa', fitdddaaa.fun * 627.5, repr(fitdddaaa.x))

    t = dic['fitddaaaa']['fun']
    temp = np.array(dic['fitddaaaa']['x'], dtype=np.double)
    if fitddaaa.fun < t:
        t = fitddaaa.fun
        temp = fitddaaa.x
    if fitdaaaa.fun < t:
        t = fitdaaaa.fun
        temp = fitdaaaa.x
    fitddaaaa = opt.basinhopping(MSD, x0=temp,
                                 minimizer_kwargs={'args': (data[0], data[1], es, Potentials),
                                                   'bounds': ((0, np.inf), (0, 10),
                                                              (ra - s, ra + s), (-np.inf, np.inf),
                                                              (0, 1),
                                                              (0, np.inf), (-np.inf, np.inf), (0 - s, 0 + s),
                                                              (0 - s, 0 + s),
                                                              (0, np.inf), (-np.inf, np.inf), (-np.inf, np.inf),
                                                              (-np.inf, np.inf),
                                                              ),
                                                   # 'method': "Nelder-Mead"
                                                   }, niter=ni)
    print('fitdddaaa', fitdddaaa.fun * 627.5, repr(fitdddaaa.x))

    t = dic['fitddddaaa']['fun']
    temp = np.array(dic['fitddddaaa']['x'], dtype=np.double)
    if fitddddaa.fun < t:
        t = fitddddaa.fun
        temp = fitddddaa.x
    if fitdddaaa.fun < t:
        t = fitdddaaa.fun
        temp = fitdddaaa.x
    fitddddaaa = opt.basinhopping(MSD, x0=temp,
                                  minimizer_kwargs={'args': (data[0], data[1], es, Potentials),
                                                    'bounds': ((0, np.inf), (0, 10),
                                                               (ra - s, ra + s), (-np.inf, np.inf),
                                                               (0, 1),
                                                               (0, np.inf), (-np.inf, np.inf), (-np.inf, np.inf),
                                                               (-np.inf, np.inf),
                                                               (0, np.inf), (-np.inf, np.inf), (-np.inf, np.inf),
                                                               (0 - s, 0 + s),
                                                               ),
                                                    # 'method': "Nelder-Mead"
                                                    }, niter=ni)
    print('fitdddddaa', fitddddaaa.fun * 627.5, repr(fitddddaaa.x))

    t = dic['fitdddaaaa']['fun']
    temp = np.array(dic['fitdddaaaa']['x'], dtype=np.double)
    if fitdddaaa.fun < t:
        t = fitdddaaa.fun
        temp = fitdddaaa.x
    if fitddaaaa.fun < t:
        t = fitddaaaa.fun
        temp = fitddaaaa.x
    fitdddaaaa = opt.basinhopping(MSD, x0=temp,
                                  minimizer_kwargs={'args': (data[0], data[1], es, Potentials),
                                                    'bounds': ((0, np.inf), (0, 10),
                                                               (ra - s, ra + s), (-np.inf, np.inf),
                                                               (0, 1),
                                                               (0, np.inf), (-np.inf, np.inf), (-np.inf, np.inf),
                                                               (0 - s, 0 + s),
                                                               (0, np.inf), (-np.inf, np.inf), (-np.inf, np.inf),
                                                               (-np.inf, np.inf),
                                                               ),
                                                    # 'method': "Nelder-Mead"
                                                    }, niter=ni)
    print('fitdddaaaa', fitdddaaaa.fun * 627.5, repr(fitdddaaaa.x))

    t = dic['fitddddaaaa']['fun']
    temp = np.array(dic['fitddddaaaa']['x'], dtype=np.double)
    if fitddddaaa.fun < t:
        t = fitddddaaa.fun
        temp = fitddddaaa.x
    if fitdddaaaa.fun < t:
        t = fitdddaaaa.fun
        temp = fitdddaaaa.x
    fitddddaaaa = opt.basinhopping(MSD, x0=temp,
                                   minimizer_kwargs={'args': (data[0], data[1], es, Potentials),
                                                     'bounds': ((0, np.inf), (0, 10),
                                                                (ra - s, ra + s), (-np.inf, np.inf),
                                                                (0, 1),
                                                                (0, np.inf), (-np.inf, np.inf), (-np.inf, np.inf),
                                                                (-np.inf, np.inf),
                                                                (0, np.inf), (-np.inf, np.inf), (-np.inf, np.inf),
                                                                (-np.inf, np.inf),
                                                                ),
                                                     # 'method': "Nelder-Mead"
                                                     }, niter=ni)
    print('fitddddaaaa', fitddddaaaa.fun * 627.5, repr(fitddddaaaa.x))

    print(
        dic['fitda']['fun'] - fitda.fun, '\n',
        dic['fitdda']['fun'] - fitdda.fun,
        dic['fitdaa']['fun'] - fitdaa.fun, '\n',
        dic['fitddda']['fun'] - fitddda.fun,
        dic['fitddaa']['fun'] - fitddaa.fun,
        dic['fitdaaa']['fun'] - fitdaaa.fun, '\n',
        dic['fitdddda']['fun'] - fitdddda.fun,
        dic['fitdddaa']['fun'] - fitdddaa.fun,
        dic['fitddaaa']['fun'] - fitddaaa.fun,
        dic['fitdaaaa']['fun'] - fitdaaaa.fun, '\n',
        dic['fitddddaa']['fun'] - fitddddaa.fun,
        dic['fitdddaaa']['fun'] - fitdddaaa.fun,
        dic['fitddaaaa']['fun'] - fitddaaaa.fun, '\n',
        dic['fitddddaaa']['fun'] - fitddddaaa.fun,
        dic['fitdddaaaa']['fun'] - fitdddaaaa.fun, '\n',
        dic['fitddddaaaa']['fun'] - fitddddaaaa.fun, )
    improve = np.abs(dic['fitda']['fun'] - fitda.fun) + \
              np.abs(dic['fitdda']['fun'] - fitdda.fun) + \
              np.abs(dic['fitdaa']['fun'] - fitdaa.fun) + \
              np.abs(dic['fitddda']['fun'] - fitddda.fun) + \
              np.abs(dic['fitddaa']['fun'] - fitddaa.fun) + \
              np.abs(dic['fitdaaa']['fun'] - fitdaaa.fun) + \
              np.abs(dic['fitdddda']['fun'] - fitdddda.fun) + \
              np.abs(dic['fitdddaa']['fun'] - fitdddaa.fun) + \
              np.abs(dic['fitddaaa']['fun'] - fitddaaa.fun) + \
              np.abs(dic['fitdaaaa']['fun'] - fitdaaaa.fun) + \
              np.abs(dic['fitddddaa']['fun'] - fitddddaa.fun) + \
              np.abs(dic['fitdddaaa']['fun'] - fitdddaaa.fun) + \
              np.abs(dic['fitddaaaa']['fun'] - fitddaaaa.fun) + \
              np.abs(dic['fitddddaaa']['fun'] - fitddddaaa.fun) + \
              np.abs(dic['fitdddaaaa']['fun'] - fitdddaaaa.fun) + \
              np.abs(dic['fitddddaaaa']['fun'] - fitddddaaaa.fun)
    if improve < 1e-12:
        print("Converged")

    for i in fitda:
        try:
            fitda[i] = fitda[i].tolist()
            fitdda[i] = fitdda[i].tolist()
            fitdaa[i] = fitdaa[i].tolist()
            fitddda[i] = fitddda[i].tolist()
            fitddaa[i] = fitddaa[i].tolist()
            fitdaaa[i] = fitdaaa[i].tolist()
            fitdddda[i] = fitdddda[i].tolist()
            fitdddaa[i] = fitdddaa[i].tolist()
            fitddaaa[i] = fitddaaa[i].tolist()
            fitdaaaa[i] = fitdaaaa[i].tolist()
            fitddddaa[i] = fitddddaa[i].tolist()
            fitdddaaa[i] = fitdddaaa[i].tolist()
            fitddaaaa[i] = fitddaaaa[i].tolist()
            fitddddaaa[i] = fitddddaaa[i].tolist()
            fitdddaaaa[i] = fitdddaaaa[i].tolist()
            fitddddaaaa[i] = fitddddaaaa[i].tolist()
        except:
            continue
    dict = {'fitda': {'fun': fitda.fun, 'x': fitda.x},
            'fitdda': {'fun': fitdda.fun, 'x': fitdda.x},
            'fitdaa': {'fun': fitdaa.fun, 'x': fitdaa.x},
            'fitddda': {'fun': fitddda.fun, 'x': fitddda.x},
            'fitddaa': {'fun': fitddaa.fun, 'x': fitddaa.x},
            'fitdaaa': {'fun': fitdaaa.fun, 'x': fitdaaa.x},
            'fitdddda': {'fun': fitdddda.fun, 'x': fitdddda.x},
            'fitdddaa': {'fun': fitdddaa.fun, 'x': fitdddaa.x},
            'fitddaaa': {'fun': fitddaaa.fun, 'x': fitddaaa.x},
            'fitdaaaa': {'fun': fitdaaaa.fun, 'x': fitdaaaa.x},
            'fitddddaa': {'fun': fitddddaa.fun, 'x': fitddddaa.x},
            'fitdddaaa': {'fun': fitdddaaa.fun, 'x': fitdddaaa.x},
            'fitddaaaa': {'fun': fitddaaaa.fun, 'x': fitddaaaa.x},
            'fitddddaaa': {'fun': fitddddaaa.fun, 'x': fitddddaaa.x},
            'fitdddaaaa': {'fun': fitdddaaaa.fun, 'x': fitdddaaaa.x},
            'fitddddaaaa': {'fun': fitddddaaaa.fun, 'x': fitddddaaaa.x}}
    str = json.dumps(dict, indent=4)
    file=open("Fits_H2O_NH3_Polar.json",'w')
    file.write(str)
    file.close()
    if improve<1e-12: break
sys.exit()
