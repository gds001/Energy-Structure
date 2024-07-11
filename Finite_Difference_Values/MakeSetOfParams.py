import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import scipy.optimize as opt
import scipy.interpolate as inter
from Useful import CD1, CD2, CD3, RMSD
plot=False
def extrapSurf(rs,ra,alpha,da,cs):
    rs=rs-ra
    ys=np.exp(-alpha*rs)
    return da*(ys-1)**2+cs[0]+cs[1]*rs+cs[2]/2*rs**2

import json
import matplotlib as mpl

show=True
plt.rcParams.update({'font.size':14})
np.set_printoptions(suppress=True,precision=10)

donor="h2o"
acceptor="h2o"
bs="avtz"
linear="_linear"
l=""
#if linear=="_linear": l="_l"
if donor=='hf' and acceptor=='nh3':
    linear=""
    l=""
monomer=np.genfromtxt("Outputs/monomers/{}/{}_mp2_{}_scan.csv".format(donor,donor,bs),delimiter=',',dtype=np.double).T
acc=np.genfromtxt("Outputs/monomers/{}/{}_mp2_{}_scan.csv".format(acceptor,acceptor,bs),delimiter=',',dtype=np.double).T
dimer2d=np.genfromtxt("Outputs/Dimers/{}_{}{}/{}_{}{}_mp2_{}_rahs_rab.csv".format(donor,acceptor,l,donor,acceptor,linear,bs),delimiter=',',dtype=np.double).T

rabs=np.unique(dimer2d[1])

Ed=np.amin(monomer[1])
Ea=np.amin(acc[1])
Ehb=np.amin(dimer2d[3])-Ed-Ea
m=np.argmin(monomer[1])
n=np.argmin(dimer2d[3])
monomer[1]-=Ed
dimer2d[3]-=-Ehb+np.amin(dimer2d[3])
rA=monomer[0][m]
#re=dimer2d[0][n]
Re=dimer2d[1][np.argmin(dimer2d[3])]
msl=[m-3,m-2,m-1,m,m+1,m+2,m+3]
print(monomer[0])
print(monomer[0][msl],monomer[1][msl])
kA=CD2(monomer[0][msl],monomer[1][msl])
d3=CD3(monomer[0][msl],monomer[1][msl])
print(kA,d3)

alpha=d3/kA/-3
dA=kA/2/alpha**2
#dre=(re-rA)
print(kA,dA)
try:
    file=open("./DerivParams_All.json",'r')
    string=file.read()
    file.close()
    data=json.loads(string)
except: data={}

for i in range(len(rabs)):
    Rab=rabs[i]
    #if Rab<Re: continue
    print(dimer2d,rabs[i])
    dimer=dimer2d[:,dimer2d[1]==rabs[i]]
    m=1000000
    m=np.where(dimer[0]==rA)[0][0]
    print(dimer)
    #if np.amin(dimer[3])>0: continue
    #dimsurf=inter.CubicSpline(dimer[0],dimer[3])
    #rs=np.linspace(rA-.1,rA+.2,1_000_000)
    #es= dimsurf(rs)
    print(dimer[0])
    mst = [m - 3, m - 2, m - 1, m, m + 1, m + 2, m + 3]
    mst = [m - 1, m, m + 1]
    ms = [m - 6, m - 4, m - 2, m, m + 2, m + 4, m + 6]
    msl = [m - 9, m - 6, m - 3, m, m + 3, m + 6, m + 9]
    drs=dimer[0][mst[1:]]-dimer[0][mst[:-1]]
    print(drs)
    for i in range(1,len(drs)):
        if drs[i-1]-drs[i]>1e-4:
            print("Wrong Minimum", dimer[0][mst], drs)
            input()
    print(dimer[0][mst],dimer[0][ms],dimer[0][msl],)
    print(dimer[3][ms])
    Ezdr=dimer[3][m]
    c0=Ezdr
    s=1e-1
    print(CD1(dimer[0][mst],dimer[3][mst]))
    c1=CD1(dimer[0][mst],dimer[3][mst])
    print(dimer[0][mst],dimer[3][mst])
    c2=CD2(dimer[0][mst],dimer[3][mst])-kA
    #c3=CD3(dimer[0][mst],dimer[3][mst])-d3
    c3=0
    #print(c1,c2,c3)

    def fdv(x,a,d,rA,c1,c2,c3):
        y=np.exp(-a*(x-rA))
        return 2*a*d*(y-y**2)+c1+c2*(x-rA)+1/2*c3*(x-rA)**2
    re=scipy.optimize.fsolve(fdv,x0=1.0,args=(alpha,dA,rA,c1,c2,c3))[0]
    print(re)
    khb=kA*(2*np.exp(-2*alpha*(re-rA))-np.exp(-alpha*(re-rA)))+c2+c3*(re-rA)
    print(kA, khb)
    print("{},{},{},{},"
      "{},{},{},{},"
          "{},{},"
      "{},{},{},{}".format(rA,kA,dA,alpha,
                           c0,c1,c2,c3,
                           Re,Rab,
                           RMSD(dimer[3][2:-2],extrapSurf(dimer[0][2:-2],rA,alpha,dA,[c0,c1,c2]))*627.5,
                           re,0,-np.amin(dimer[3])*627.5))
    rs = dimer[0][dimer[3] <= np.amin(dimer[3])+1]-rA
    data["{}_{}_{}_{}".format(donor,acceptor,bs,Rab)]={'donor':donor,
                                                       'acceptor':acceptor,
                                                       'bs':bs,
                                                       'ra':rA,
                                                       'ka':kA,
                                                       'da':dA,
                                                       'alpha':alpha,
                                                       'c0':c0,
                                                       'c1':c1,
                                                       'c2':c2,
                                                       're':re,
                                                       'R':Rab,
                                                       'Re':Re,
                                                       'khb':khb,
                                                       'Ee':-np.amin(dimer[3])*627.5,
                                                       'RMSD':RMSD(dimer[3][2:-2],extrapSurf(dimer[0][2:-2],rA,alpha,dA,[c0,c1,c2]))*627.5}
    #rs=np.linspace(-0.15,0.2,1000)
    rs=np.linspace(-.1,.2,1_000_000)
    mon=dA*(np.exp(-alpha*(rs))-1)**2
    plt.scatter(monomer[0],monomer[1],c='tab:orange',marker='.')
    plt.scatter(dimer[0],dimer[3],c='tab:blue',marker='.')
    plt.plot(rs+rA,mon,c='tab:orange')
    plt.plot(rs+rA,c0+c1*rs+0.5*c2*rs**2,c='tab:green')
    plt.plot(rs+rA,mon+c0+c1*rs+0.5*c2*rs**2,c='tab:blue')
    plt.xlim(re-0.15,re+0.2)
    plt.ylim(np.amin(np.concatenate(([0],dimer[3]))-1/627.5),np.amax(dimer[3]+1/627.5))

    #plt.ylim(-20,25)
    if c0>0:
        print("Unbound?, {}".format(Ehb))
        if plot: plt.show()
    elif re<rA:
        print("Shorter, {}".format(rA-re))
        if plot: plt.show()

    elif c2>0:
        print("Blue shift?!")
        if plot: plt.show()
    plt.close()
string=json.dumps(data,skipkeys=True,indent=4)
print(string)
file=open("./DerivParams_All.json",'w')
file.write(string)
file.close()
