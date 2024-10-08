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
ad1=[]

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
    dr = dre[i][o]
    c0s = c0[i][o]

    psPolar=np.copy(ps)
    psPolar[0]=0

    axes[k,j].plot(rs*100,c0dr(Rs,ka,a,ps)*6.275,c='k')

    axes[k,j].plot(rs*100,(-ka*ps[6]*(2*zs**2-zs))*6.275,c='tab:blue',lw=5,zorder=-10)
    axes[k,j].plot(rs*100,(-ka*ps[6]*(2*zs**2-zs))*6.275,c='tab:orange',lw=2,zorder=-9)

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

    axes[k,j].plot(rs*100,(c0dr(Rs,ka,a,psPolar)+ka*ps[6]*(2*zs**2-zs))*6.275,c='tab:red',lw=2,zorder=-7)

    def func1(x,a,b,c): return (a*x*x+b*x)/(x+c)
    fit1=opt.basinhopping(MAD,x0=np.array([1,0,0]),
                          minimizer_kwargs={'args':(rss,c0dr(Rss,ka,a,psPolar)+ka*ps[6]*(2*zss**2-zss),func1),
                                            'bounds':((-np.inf,np.inf),(-10,10),(0,1e-2))},
                          niter=500).x
    print(fit1)

    ri=np.amax(dres[i])
    if i>=1: s=0.05
    def func2(x,A,ri): return -A/(x-ri)-A/ri
    fit2=opt.basinhopping(MAD,x0=np.array([s,ri]),
                          minimizer_kwargs={'args':(rss,c0dr(Rss,ka,a,ps)-c0dr(Rss,ka,a,psPolar),func2),
                                              'bounds':((-np.inf,np.inf),(ri-1e-12,ri+s))})
    print(fit2.fun,fit2.x,np.amax(dres[i]))
    fit2=fit2.x

    axes[k,j].plot(rs*100,func1(rs,*fit1)*6.275,c='tab:green',zorder=-8,lw=5)
    axes[k,j].plot(rs*100,(c0dr(Rs,ka,a,ps)-c0dr(Rs,ka,a,psPolar))*6.275,c='tab:olive',lw=2,zorder=-5)
    axes[k,j].plot(rs * 100, func2(rs, *fit2) * 6.275, c='tab:purple', lw=5, zorder=-6)
    axes[k,j].plot(rs * 100, (-ps[6] * ka*(2*zs**2-zs)
                              + func1(rs,*fit1) + func2(rs, *fit2)) * 6.275, c='tab:gray', lw=5, zorder=-1)

    axes[k,j].set_title(titles[i],fontsize=12)
    axes[k,j].set_ylim(-10,5)
    axes[k,j].set_xlim(0,np.amin([np.amax(rss),np.amax(dre[i])+0.01])*100)
    Rticks=np.linspace(Res[i],Res[i]+2,3)
    if k==1: axes[k,j].set_xlabel("$\Delta r_\mathrm{HB}$ (pm)")
    axes[k,j].set_yticks([-10,-5,0,5],["$-10$","$-5$","$0$","$5$"])

#fig.suptitle(type)
axes[0,0].set_xticks([0,0.2,0.4],["$0.0$","$0.2$","$0.4$"])
axes[0,1].set_xticks([0,0.3,0.6],["$0.0$","$0.3$","$0.6$"])
axes[0,2].set_xticks([0,0.3,0.6],["$0.0$","$0.3$","$0.6$"])
axes[1,0].set_xticks([0,0.5,1.0],["$0.0$","$0.5$","$1.0$"])
axes[1,1].set_xticks([0,0.8,1.6],["$0.0$","$0.8$","$1.6$"])
axes[1,2].set_xticks([0,1.5,3.0],["$0.0$","$1.5$","$3.0$"])
axes[1,0].text(-0.6,0," ",ha='right',va='center',rotation='vertical')
plt.tight_layout()
axes[1,0].text(-0.6,10,"$\partial_{\Delta r_\mathrm{HB}}c^{(0)}$  (kcal/mol pm$^{-1}$)",ha='right',va='center',rotation='vertical')
plt.savefig("../Figures/DipoleDerivs1_{}.png".format(type),dpi=600)
plt.show()
plt.close()

from matplotlib.lines import Line2D
ls=[Line2D((0,0),(0,0),c='k',lw=2),
    Line2D((0,0),(0,0),c='tab:orange',lw=2),
    Line2D((0,0),(0,0),c='tab:gray',lw=5),
    Line2D((0,0),(0,0),c='tab:blue',lw=5),
    Line2D((0,0),(0,0),c='tab:red',lw=2),
    Line2D((0,0),(0,0),c='tab:olive',lw=2),
    Line2D((0,0),(0,0),c='tab:green',lw=5),
    Line2D((0,0),(0,0),c='tab:purple',lw=5),
    ]

plt.legend(ls,["$\partial_{\Delta r_\mathrm{HB}} c^{(0)}$",
               "$\partial_{\Delta r_\mathrm{HB}} c^{(0)}_E$",
               "approx.$\partial_{\Delta r_\mathrm{HB}} c^{(0)}$",
               "approx.$\partial_{\Delta r_\mathrm{HB}} c^{(0)}_E$",
               "$\partial_{\Delta r_\mathrm{HB}} c^{(0)}_P$",
               "$\partial_{\Delta r_\mathrm{HB}} c^{(0)}_R$",
               "approx.$\partial_{\Delta r_\mathrm{HB}} c^{(0)}_P$",
               "approx.$\partial_{\Delta r_\mathrm{HB}} c^{(0)}_R$"],
           ncols=4)

plt.show()

from matplotlib.lines import Line2D
ls=[Line2D((0,0),(0,0),c='k',lw=2),
    Line2D((0,0),(0,0),c='tab:orange',lw=2),
    Line2D((0,0),(0,0),c='tab:red',lw=2),
    Line2D((0,0),(0,0),c='tab:olive',lw=2),
    Line2D((0,0),(0,0),c='tab:gray',lw=5),
    Line2D((0,0),(0,0),c='tab:blue',lw=5),
    Line2D((0,0),(0,0),c='tab:green',lw=5),
    Line2D((0,0),(0,0),c='tab:purple',lw=5),
    ]

plt.legend(ls,["$\partial_{\Delta r_\mathrm{HB}} c^{(0)}$",
               "$\partial_{\Delta r_\mathrm{HB}} c^{(0)}_E$",
               "$\partial_{\Delta r_\mathrm{HB}} c^{(0)}_P$",
               "$\partial_{\Delta r_\mathrm{HB}} c^{(0)}_R$",
               "approx.$\partial_{\Delta r_\mathrm{HB}} c^{(0)}$",
               "approx.$\partial_{\Delta r_\mathrm{HB}} c^{(0)}_E$",
               "approx.$\partial_{\Delta r_\mathrm{HB}} c^{(0)}_P$",
               "approx.$\partial_{\Delta r_\mathrm{HB}} c^{(0)}_R$"],
           ncols=2)

plt.show()