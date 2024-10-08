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

files=["../d1_p2/Fits_NH3_NH3_Polar.json","../d1_p2/Fits_H2O_H2O_Polar.json","../d1_p2/Fits_HF_HF_Polar.json",
       "../d1_p2/Fits_H2O_NH3_Polar.json","../d1_p2/Fits_HF_H2O_Polar.json","../d1_p2/Fits_HF_NH3_Polar.json"]
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

file = open("../DerivParams_All.json", 'r')
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
c1=[]
for i in range(6):
    temp=data[data['donor']==don[i]]
    temp=temp[temp['acceptor']==acc[i]]
    temp=temp[temp['re']>temp['ra']]
    ras.append(temp['ra'][0])
    dR.append(temp['R'].to_numpy(dtype='float')-temp['Re'].to_numpy(dtype='float'))
    c1.append(temp['c1'].to_numpy(dtype='float'))
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

ms=[]
zetas=[]
xis=[]

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
    c1s = c1[i][o]
    dr=dr[RS>=0]
    c1s=c1s[RS>=0]
    psPolar=np.copy(ps)
    psPolar[0]=0
    axes[k,j].scatter(dr*100,c1s*627.5/100,marker='.',c='k')
    axes[k,j].plot(rs*100,c1R(Rs,*ps)*627.5/100,c='tab:gray')


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


    def func4(x, a, b):
        return a * x * x/2 + b * x
    def func4i(x,a,b):
        return a*x**3/6 + b*x**2/2
    fit4 = opt.minimize(MAD, x0=np.array([1, 0]),
                        args=(rss, c1dr(Rss, ka, a, psPolar) + ka * (2 * zss ** 2 - zss), func4),
                        bounds=((-np.inf, np.inf), (-np.inf, np.inf))).x
    print(fit4)

    ri = np.amax(dres[i])
    if i >= 1: s = 0.05
    def func5(x, a, b):
        return a * x * x / 2 + b * x
    def func5i(x, a, b):
        return a * x ** 3 / 6 + b * x ** 2 / 2
    fit5 = opt.minimize(MAD, x0=np.array([1, 0]),
                        args=(rss, c1dr(Rss, ka, a, ps) - c1dr(Rss, ka, a, psPolar), func5),
                        bounds=((-np.inf, np.inf), (-np.inf, np.inf)))
    print(fit5.x)
    fit5 = fit5.x

    xi1,zeta1=fit4
    xi2,zeta2=fit5
    ms.append(2*a*d)
    zetas.append(zeta1+zeta2)
    xis.append(xi1+xi2)

    dz=np.exp(-a*dr)
    ans=2*a*d*(dz**2-dz)+ func4i(dr,*fit4) + func5i(dr, *fit5)
    axes[k,j].plot(rs * 100, (2*a*d*(zs**2-zs)
                              + func4i(rs,*fit4) + func5i(rs, *fit5)) * 627.5/100, c='k')
    axes[k, j].text(np.amax(rs) * 100 * 0.975, np.amin(c1s)*627.5/100 * 1.1 * 0.025,
                    "MaxE: {:.4f}\nMAE: {:.4f}\n $\\frac{{\mathrm{{kcal/mol}}}}{{\mathrm{{pm}}}}$".format(np.amax(np.abs(c1s -ans))* 627.5/100,
                                                                    MAE(c1s,ans)*627.5/100),
                    va='top', ha='right', fontsize=10)
    axes[k,j].set_title(titles[i],fontsize=12)
    axes[k,j].set_ylim(np.amin(c1s)*627.5*1.1/100,0)
    axes[k,j].set_xlim(0,np.amax(rs)*100)
    Rticks=np.linspace(Res[i],Res[i]+2,3)
    if k==1: axes[k,j].set_xlabel("$\Delta r_\mathrm{HB}$ (pm)")

#fig.suptitle(type)
axes[0,0].set_xticks([0,0.2,0.4],["$0.0$","$0.2$","$0.4$"],fontsize=10)
axes[0,1].set_xticks([0,0.3,0.6],["$0.0$","$0.3$","$0.6$"],fontsize=10)
axes[0,2].set_xticks([0,0.3,0.6],["$0.0$","$0.3$","$0.6$"],fontsize=10)
axes[1,0].set_xticks([0,0.6,1.2],["$0.0$","$0.6$","$1.2$"],fontsize=10)
axes[1,1].set_xticks([0,0.8,1.6],["$0.0$","$0.8$","$1.6$"],fontsize=10)
axes[1,2].set_xticks([0,1.5,3.0],["$0.0$","$1.5$","$3.0$"],fontsize=10)

axes[0,0].set_yticks([0,-0.02,-0.04],["$0.00$","$-0.02$","$-0.04$"],fontsize=10)
axes[0,1].set_yticks([0,-0.04,-0.08],["$0.00$","$-0.04$","$-0.08$"],fontsize=10)
axes[0,2].set_yticks([0,-0.04,-0.08],["$0.00$","$-0.04$","$-0.08$"],fontsize=10)
axes[1,0].set_yticks([0,-0.05,-0.10,-0.15],["$0.00$","$-0.05$","$-0.10$","$-0.15$"],fontsize=10)
axes[1,1].set_yticks([0,-0.07,-0.14,-0.21],["$0.00$","$-0.07$","$-0.14$","$-0.21$"],fontsize=10)
axes[1,2].set_yticks([0,-0.2,-0.4],["$0.00$","$-0.20$","$-0.04$"],fontsize=10)

axes[1,0].text(-0.7,0," ",ha='right',va='center',rotation='vertical')
plt.tight_layout()
axes[1,0].text(-0.7,0.05,"$c^{(1)}$  (kcal/mol pm$^{-1}$)",ha='right',va='center',rotation='vertical')
plt.subplots_adjust(left=0.14)
plt.savefig("../Figures/c1Vsdr1_fits_{}.png".format(type),dpi=600)
plt.show()

caption=["$2\\alpha D_A $ & (kcal/mol pm$^{-1}$)",
         "$\\zeta$ & (kcal/mol pm$^{-3}$)",
         "$\\xi$ & (kcal/mol pm$^{-4}$)"]
str=caption[0]
for i in range(6):
    str+=" & {:.4f}".format(ms[i]*627.5)
print(str+"\\\\\hline")
str=caption[1]
for i in range(6):
    str+=" & {:.4f}".format(zetas[i]*627.5/100/100/100)
print(str+"\\\\")
str=caption[2]
for i in range(6):
    str+=" & {:.4f}".format(xis[i]*627.5/100/100/100/100)
print(str+"\\\\\hline\hline")

