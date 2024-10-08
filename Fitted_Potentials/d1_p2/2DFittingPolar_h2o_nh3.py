import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import matplotlib as mpl
import pandas as pd
import json
import scipy.interpolate as inter
from Useful import CD1,CD2,CD3,hartree2wn
from A_Figures_for_Paper_1.ExpansionAndDerivatives import Potentials
import sys
print("h2o_nh3")
np.set_printoptions(suppress=True,precision=4,linewidth=1000)
np.seterr(all='ignore', divide='ignore', over='ignore', under='ignore', invalid='ignore')
data=np.genfromtxt("h2o_nh3_linear_mp2_avtz_rahs_rab.csv",delimiter=',',dtype=np.double).T
don=np.genfromtxt("../h2o/h2o_mp2_avtz_scan.csv",delimiter=',',dtype=np.double).T
acc=np.genfromtxt("../nh3/nh3_mp2_avtz_scan.csv",delimiter=',',dtype=np.double).T

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
file=open("../h2o_h2o/Fits_H2O_H2O_Polar.json",'r')
dic=json.loads(file.read())
ps=dic['fitddaaa']['x']
print(ps)
dm1d2 =ps[6]
while True:
    try:
        file=open("Fits_H2O_NH3_Polar.json",'r')
        dic=json.loads(file.read())
    except:
        print("Starting Fresh")
        dic={
             'fitddaaa':{'x':[
            222.76610874638135,
            3.445084061619073,
            0.961300000001,
            0.22086139739776348,
            0.4235819221082713,
            0.11382995548193271,
            dm1d2,
            40.34174744450823,
            0.7685673098210626,
            211.26255784032185],'fun':1}
         }
    ni=100_000
    print(ni,"h2o_nh3")
    t = dic['fitddaaa']['fun']
    temp = np.array(dic['fitddaaa']['x'], dtype=np.double)
    fitddaaa = opt.basinhopping(MSD, x0=temp,
                                 minimizer_kwargs={'args': (data[0], data[1], es, Potentials),
                                                   'bounds': ((0, np.inf), (0, 10),
                                                              (ra - s, ra + s), (-np.inf, np.inf),
                                                              (0, 1),
                                                              (0, np.inf), (dm1d2-s,dm1d2+s),
                                                              (0, np.inf), (-np.inf,np.inf), (-np.inf,np.inf)
                                                              ),
                                                    #'method': "Nelder-Mead"
                                                   }, niter=ni)
    print('fitddaaa', fitddaaa.fun * hartree2wn, repr(fitddaaa.x))

    print(dic['fitddaaa']['fun'] - fitddaaa.fun)
    improve=np.abs(dic['fitddaaa']['fun']-fitddaaa.fun)

    if improve<1e-12:
        print("Converged")

    for i in fitddaaa:
        try:
            fitddaaa[i]=fitddaaa[i].tolist()
        except:
            continue
    dict={'fitddaaa':{'fun':fitddaaa.fun,'x':fitddaaa.x}}
    str=json.dumps(dict,indent=4)
    file=open("Fits_H2O_NH3_Polar.json",'w')
    file.write(str)
    file.close()
    if improve<1e-12: break
sys.exit()
