import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import matplotlib as mpl
import json
import pandas as pd
from Useful import CD1,CD2,CD3,wRMSD,MAE,wMAE,RMSD
from WorkUpFits.ExpansionAndDerivatives import *

plt.rcParams.update({'font.size':12,'font.family':'arial'})

def MAD(params, xs, ys, func):
    return np.sum(np.abs(ys - func(xs, *params)))
def MSD(params, xs, ys, func):
    return np.sum(np.abs(ys - func(xs, *params)) ** 2)
def MLD(params, xs, ys, func):
    return np.sum(np.log10(np.abs(ys - func(xs, *params)) + 1))

files=["../Fitted_Potentials/d1_p2/Fits_NH3_NH3_Polar.json","../Fitted_Potentials/d1_p2/Fits_H2O_H2O_Polar.json","../Fitted_Potentials/d1_p2/Fits_HF_HF_Polar.json",
       "../Fitted_Potentials/d1_p2/Fits_H2O_NH3_Polar.json","../Fitted_Potentials/d1_p2/Fits_HF_H2O_Polar.json","../Fitted_Potentials/d1_p2/Fits_HF_NH3_Polar.json"]
params=[]
for i in range(6):
    try:
        file=open(files[i],'r')
        params.append(json.loads(file.read()))
        file.close()
    except:
        params.append({
    "fitgn": {"fun": 1,"x": [0,0,0,0,0,0,0,0,0,0]},
    "fitn": {"fun": 1,"x": [0,0,0,0,0,0,0,0,0,0]},
    "fitg": {"fun": 1,"x": [0,0,0,0,0,0,0,0,0,0]},
    "fit": {"fun": 1,"x": [0,0,0,0,0,0,0,0,0,0]}})
        continue

file = open("../Finite_Difference_Values/DerivParams_All.json", 'r')
string = file.read()
file.close()
data = json.loads(string)
data=pd.DataFrame(data).T
colors=[]
for entry in data.T:
    if data.T[entry]['donor']=='hf':
        if data.T[entry]['acceptor']=='hf':colors.append("tab:green")
        if data.T[entry]['acceptor']=='h2o':colors.append("tab:brown")
        if data.T[entry]['acceptor']=='nh3':colors.append("tab:cyan")
    elif data.T[entry]['donor'] == 'h2o':
        if data.T[entry]['acceptor'] == 'h2o': colors.append("tab:red")
        if data.T[entry]['acceptor'] == 'nh3': colors.append("tab:purple")
    else: colors.append('tab:blue')
dR=[]
dre=[]
dres=[]
don=["nh3","h2o","hf","h2o","hf","hf"]
acc=['nh3','h2o','hf','nh3','h2o','nh3']
ass=[]
das=[]
Res=[]
wss=[]
Es=[]
ras=[]
c0=[]
for i in range(6):
    temp=data[data['donor']==don[i]]
    temp=temp[temp['acceptor']==acc[i]]
    temp=temp[temp['re']>temp['ra']]
    ras.append(temp['ra'][0])
    dR.append(temp['R'].to_numpy(dtype='float')-temp['Re'].to_numpy(dtype='float'))
    c0.append(temp['c0'].to_numpy(dtype='float'))
    dre.append(temp['re'].to_numpy(dtype='float') - temp['ra'].to_numpy(dtype='float'))
    dres.append(np.linspace(0,np.amax(temp['re'].to_numpy(dtype='float')-temp['ra'].to_numpy(dtype='float')),1000))
    Es.append(temp['Ee'].to_numpy(dtype='float'))
    ass.append(temp['alpha'][0])
    das.append(temp['da'][0])
    Res.append(temp['Re'][0])
colors=['tab:blue','tab:red','tab:green','tab:purple','tab:brown','tab:cyan']
o=0.5
s=1e-8
type='fitddaaa'
print(type)
titles=["NH$_3\cdot$NH$_3$","H$_2$O$\cdot$H$_2$O","HF$\cdot$HF",
        "H$_2$O$\cdot$NH$_3$","HF$\cdot$H$_2$O","HF$\cdot$NH$_3$"]

kamds=[]
xs=[]
ys=[]
zzs=[]

fig,axes=plt.subplots(2,3,figsize=(6,4))
for i in range(6):
    s=1e-12
    temp = data[data['donor'] == don[i]]
    temp = temp[temp['acceptor'] == acc[i]]
    print(files[i])
    print(Res[i])
    ps = params[i][type]['x']
    a=ass[i]
    d=das[i]
    ka=2*a*a*d
    j = i % 3
    k = i // 3
    Rs=np.logspace(2,np.log10(Res[i]-0.1),1000)
    rs=Rr(Rs,ka,a,ps)
    rss=rs[Rs>=Res[i]]
    Rss=Rs[Rs>=Res[i]]
    zs=np.exp(-a*rs)
    zss = np.exp(-a * rss)

    o = np.argsort(dre[i])
    RS=dR[i][o]
    dr = dre[i][o]
    c0s = c0[i][o]
    Ess=Es[i][o]/627.5
    dr=dr[RS>=0]
    c0s=c0s[RS>=0]
    Ess=Ess[RS>=0]
    psPolar=np.copy(ps)
    psPolar[0]=0
    axes[k,j].scatter(dr*100,Ess*627.5,marker='.',c='k')
    axes[k,j].plot(rs*100,Evsr(rs,Rs,a,d,ps)*627.5,c='tab:gray')


    temp = c0dr(Rs, ka, a, psPolar)
    z = np.argmax(np.isnan(temp))
    if z != 0:
        Rs = Rs[:z]
        rs = rs[:z]
        zs = zs[:z]
    temp=c0dr(Rss,ka,a,psPolar)
    z=np.argmax(np.isnan(temp))
    if z!=0:
        Rss=Rss[:z]
        rss=rss[:z]
        zss=zss[:z]

    def func1(x,a,b,c): return (a*x*x+b*x)/(x+c)
    def func1i(x, a, b, c):return c * (a * c - b) * np.log(np.abs(x + c)) + a * x ** 2 / 2 + (b - a * c) * x - c * (a * c - b) * np.log(np.abs(c))
    fit1=opt.basinhopping(MAD,x0=np.array([1,0,0]),
                          minimizer_kwargs={'args':(rss,c0dr(Rss,ka,a,psPolar)+ka*ps[6]*(2*zss**2-zss),func1),
                                            'bounds':((-np.inf,np.inf),(-10,10),(0,1e-2))},
                          niter=500).x
    print(fit1)

    ri=np.amax(dres[i])
    if i>=1: s=0.05
    def func2(x,A,ri): return -A/(x-ri)-A/ri
    def func2i(x, A, ri):return -A * np.log(np.abs(x - ri)) - A * x / ri + A * np.log(np.abs(-ri))
    fit2=opt.basinhopping(MAD,x0=np.array([s,ri]),
                          minimizer_kwargs={'args':(rss,c0dr(Rss,ka,a,ps)-c0dr(Rss,ka,a,psPolar),func2),
                                              'bounds':((-np.inf,np.inf),(ri-1e-12,ri+s))})
    print(fit2.fun,fit2.x,np.amax(dres[i]))
    fit2=fit2.x

    temp = c1dr(Rs, ka, a, psPolar)
    z = np.argmax(np.isnan(temp))
    if z != 0:
        Rs = Rs[:z]
        rs = rs[:z]
        zs = zs[:z]
    temp=c1dr(Rss,ka,a,psPolar)
    z=np.argmax(np.isnan(temp))
    if z!=0:
        Rss=Rss[:z]
        rss=rss[:z]
        zss=zss[:z]

    def func4(x, a, b):return a * x * x/2 + b * x
    def func4i(x,a,b):return a*x**3/6 + b*x**2/2
    fit4 = opt.minimize(MAD, x0=np.array([1, 0]),
                        args=(rss, c1dr(Rss, ka, a, psPolar) + ka * (2 * zss ** 2 - zss), func4),
                        bounds=((-np.inf, np.inf), (-np.inf, np.inf))).x
    print(fit4)

    def func5(x, a, b):return a * x * x / 2 + b * x
    def func5i(x, a, b):return a * x ** 3 / 6 + b * x ** 2 / 2
    fit5 = opt.minimize(MAD, x0=np.array([1, 0]),
                        args=(rss, c1dr(Rss, ka, a, ps) - c1dr(Rss, ka, a, psPolar), func5),
                        bounds=((-np.inf, np.inf), (-np.inf, np.inf)))
    print(fit5.x)
    fit5 = fit5.x

    b,c,ro=fit1
    B,ri=fit2
    xi1,zeta1=fit4
    xi2,zeta2=fit5
    zeta=zeta1+zeta2
    xi=xi1+xi2
    C=ro*(b*ro-c)

    x=3*a*ka/2*ps[6]+b/2-ka/2
    y=-7*a*a*ka/6*ps[6]+a*ka+zeta/2
    z=15*a*a*a*ka/24*ps[6]-7/8*a*a+xi/6
    axes[k, j].plot(rs * 100, (ka*ps[6]*rs-x*rs**2-y*rs**3-z*rs**4
                               -np.log((rs+ro)**C/(ri-rs)**B*ri**B/ro**C)
                               +C/ro*rs+B/ri*rs
                               ) * 627.5, c='tab:gray',zorder=-10)

    kamds.append(ka*ps[6])
    xs.append(x)
    ys.append(y)
    zzs.append(z)

    dz=np.exp(-a*dr)
    ans=(-d*(dz-1)**2
                           -2*a*d*ps[6]*(dz**2-dz)
                           -2*a*d*(dz**2-dz)*dr
                           -func1i(dr,*fit1)
                           -func2i(dr,*fit2)
                           -func4i(dr,*fit4)*dr
                           -func5i(dr,*fit5)*dr)
    axes[k,j].plot(rs*100,(-d*(zs-1)**2
                           -2*a*d*ps[6]*(zs**2-zs)
                           -2*a*d*(zs**2-zs)*rs
                           -func1i(rs,*fit1)
                           -func2i(rs,*fit2)
                           -func4i(rs,*fit4)*rs
                           -func5i(rs,*fit5)*rs
                           )*627.5,c='k')
    axes[k, j].text(np.amax(dr) * 110 * 0.95, np.amax(Ess)*627.5 * 1.1 * 0.05,
                    "MaxE: {:.2f}\nMAE: {:.2f}\nkcal/mol".format(np.amax(np.abs(Ess -ans))* 627.5,
                                                                    MAE(Ess,ans)*627.5),
                    va='bottom', ha='right', fontsize=10)
    axes[k,j].set_title(titles[i],fontsize=12)
    axes[k,j].set_ylim(0,np.amax(Ess)*627.5*1.1)
    axes[k,j].set_xlim(0,np.amax(dr)*110)
    Rticks=np.linspace(Res[i],Res[i]+2,3)
    if k==1: axes[k,j].set_xlabel("$\Delta r_\mathrm{HB}$ (pm)")

#fig.suptitle(type)
axes[0,0].set_xticks([0,0.2,0.4],["$0.0$","$0.2$","$0.4$"])
axes[0,1].set_xticks([0,0.3,0.6],["$0.0$","$0.3$","$0.6$"])
axes[0,2].set_xticks([0,0.3,0.6],["$0.0$","$0.3$","$0.6$"])
axes[1,0].set_xticks([0,0.6,1.2],["$0.0$","$0.6$","$1.2$"])
axes[1,1].set_xticks([0,0.8,1.6],["$0.0$","$0.8$","$1.6$"])
axes[1,2].set_xticks([0,1.5,3.0],["$0.0$","$1.5$","$3.0$"])

axes[0,0].set_yticks([0,1,2,3],["$0$","$1$","$2$","$3$",])
axes[0,1].set_yticks([0,2,4],["$0$","$2$","$4$"])
axes[0,2].set_yticks([0,2,4],["$0$","$2$","$4$"])
axes[1,0].set_yticks([0,2,4,6],["$0$","$2$","$4$","$6$"])
axes[1,1].set_yticks([0,3,6,9],["$0$","$3$","$6$","$9$"])
axes[1,2].set_yticks([0,4,8,12],["$0$","$4$","$8$","$12$"])

axes[1,0].text(-0.4,0," ",ha='right',va='center',rotation='vertical')
plt.tight_layout()
axes[1,0].text(-0.4,9,"$E_\mathrm{HB}$  (kcal/mol)",ha='right',va='center',rotation='vertical')

plt.savefig("../Figures/Ess_fits_{}.png".format(type),dpi=600)
plt.show()

caption=["$k_A \\frac{\mu^{(0)}}{\mu^{(1)}}$ & (kcal/mol pm$^{-1}$)",
         "$x$ & (kcal/mol pm$^{-2}$)",
         '$y$ & (kcal/mol pm$^{-3}$)',
         '$z$ & (kcal/mol pm$^{-4}$)']
str=caption[0]
for i in range(6):
    str+=" & {:.4f}".format(kamds[i]*627.5/100)
print(str+"\\\\")
str=caption[1]
for i in range(6):
    str+=" & {:.4f}".format(xs[i]*627.5/100/100)
print(str+"\\\\")
str=caption[2]
for i in range(6):
    str+=" & {:.4f}".format(ys[i]*627.5/100/100/100)
print(str+"\\\\")
str=caption[3]
for i in range(6):
    str+=" & {:.4f}".format(zzs[i]*627.5/100/100/100/100)
print(str+'\\\\\hline')
