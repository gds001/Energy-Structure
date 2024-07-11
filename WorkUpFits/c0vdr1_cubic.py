import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import matplotlib as mpl
import json
import pandas as pd
from Useful import CD1,CD2,CD3,wRMSD,MAE,wMAE,RMSD,PointSlope,TrapSum
from Fitted_Potentials.d2_p0_fixed_dipole_ratios.ExpansionAndDerivatives import *
plt.rcParams.update({'font.size':10,
                    'font.family':'arial'})

np.seterr(all='ignore', divide='ignore', over='ignore', under='ignore', invalid='ignore')

dir="../Fitted_Potentials/d2_p0_fixed_dipole_ratios"
files=[dir+"/Fits_NH3_NH3_Polar.json",dir+"/Fits_H2O_H2O_Polar.json",dir+"/Fits_HF_HF_Polar.json",
       dir+"/Fits_H2O_NH3_Polar.json",dir+"/Fits_HF_H2O_Polar.json",dir+"/Fits_HF_NH3_Polar.json"]
params=[]
for i in range(6):
    file=open(files[i],'r')
    params.append(json.loads(file.read()))
    file.close()

file = open("../Finite_Difference_Values/DerivParams_All.json", 'r')
string = file.read()
file.close()
data = json.loads(string)
data=pd.DataFrame(data).T
name="All Points"
name="All Bound Points"
name="Dissociation Only"
if name=="Dissociation Only":data=data[data['R']>=data['Re']]
if name=="All Bound Points":data=data[data['Ee']>=0]
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
kas=[]
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
titles=["NH$_3\cdot$NH$_3$","H$_2$O$\cdot$H$_2$O","HF$\cdot$HF",
        "H$_2$O$\cdot$NH$_3$","HF$\cdot$H$_2$O","HF$\cdot$NH$_3$"]
type='fitddda'
print(type)
fig,axes=plt.subplots(2,3,figsize=(6,4))
qs=[]
ts=[]
for i in range(6):
    temp = data[data['donor'] == don[i]]
    temp = temp[temp['acceptor'] == acc[i]]
    print(files[i])
    #print(Res[i])
    ps = params[i][type]['x']
    #print(ps)
    a=ass[i]
    d=das[i]
    ka=2*a*a*d
    kas.append(ka)
    j = i % 3
    k = i // 3

    o = np.argsort(dre[i])
    dr = dre[i][o]
    c0s = c0[i][o]
    Ess=Es[i][o]
    Rss = np.logspace(2,np.log10(Res[i]-.1), 100000)
    rss = Rr(Rss,ka,a,ps)
    zss = np.exp(-a * rss)
    md=ps[6]
    md2=ps[7]

    axes[k, j].scatter(dr * 100, c0s * 627.5, marker='.', c='k')
    if i < 3:
        m,b=PointSlope(Res[i]+1000,Res[i]+1,Rr,c0dr2,[ka,a,ps])
    elif i < 5:
        m, b = PointSlope(Res[i] + 0.25, Res[i] + 1, Rr, c0dr2, [ka, a, ps])
    else:
        m=0
        Rs = np.logspace(2,np.log10(Res[i]), 10000)
        b=TrapSum(Rs,Rr,c0dr2,[ka,a,ps])
        print(b*6.275/100/100)
        #"""

    print(b * 6.275 / 100, m * 6.275 / 100 / 100)
    qs.append(b-3*a*ka*md-2*ka*md*md2)
    ts.append(m)
    axes[k, j].plot(rss * 100, (-md * ka * rss + b * rss ** 2 / 2 + m * rss ** 3 / 6) * 627.5, c='k')
    axes[k,j].set_title(titles[i])

    axes[k,j].text(np.amax(dr)*107.5,np.amin(c0s)*627.5*1.1*0.05,
                   "MAE: {:.2f} kcal/mol\nMax Error:\n{:.2f} kcal/mol".format(MAE(c0s,-md*ka*dr+b*dr**2/2+m*dr**3/6)*627.5,
                                                                              np.amax(np.abs(c0s-(-md*ka*dr+b*dr**2/2+m*dr**3/6)))*627.5,
                   ),
                   va='top',ha='right',fontsize=9)
    axes[k,j].set_ylim(np.amin(c0s)*627.5*1.1,0)
    axes[k,j].set_xlim(0,np.amax(dr)*110)
    if k==1: axes[k,j].set_xlabel("$\Delta r_\mathrm{HB}$ (pm)")
    if j==0: axes[k,j].set_ylabel("$c^{(0)}$  (kcal/mol)")
    print(np.amax(np.abs(c0s-(-md*ka*dr+b*dr**2/2+m*dr**3/6))*627.5))

axes[0,0].set_xticks([0.00,0.15,0.3,0.45],["$0.0$","$0.15$","$0.30$","$0.45$"])
axes[0,1].set_xticks([0.00,0.2,0.4,0.6],["$0.0$","$0.2$","$0.4$","$0.6$"])
axes[0,2].set_xticks([0.00,0.2,0.4,0.6],["$0.0$","$0.2$","$0.4$","$0.6$"])
axes[1,0].set_xticks([0,0.4,0.8,1.2],["$0.0$","$0.4$","$0.8$","$1.2$"])
axes[1,1].set_xticks([0,0.5,1,1.5],["$0.0$","$0.5$","$1.0$","$1.5$"])
axes[1,2].set_xticks([0,1,2,3],["$0.0$","$1.0$","$2.0$","$3.0$"])

axes[0,0].set_yticks([0.00,-1,-2,-3],["$0$","$-1$","$-2$","$-3$"])
axes[0,1].set_yticks([0.00,-2,-4,],["$0$","$-2$","$-4$"])
axes[0,2].set_yticks([0.00,-2,-4,],["$0$","$-2$","$-4$"])
axes[1,0].set_yticks([0.00,-2,-4,-6],["$0$","$-2$","$-4$","$-6$"])
axes[1,1].set_yticks([0.00,-3,-6,-9],["$0$","$-3$","$-6$","$-9$"])
axes[1,2].set_yticks([0.00,-4,-8,-12],["$0$","$-4$","$-8$","$-12$"])

#fig.suptitle(type)
plt.tight_layout()
plt.savefig("c0Vsdr1_cubic.png",dpi=600)
plt.show()

preamble=["$\\frac{\mu^{(0)}}{\mu^{(1)}}k_A$ & (kcal/mol pm$^{-1}$)" ,
          "$3/2\\alpha k_A\\frac{\mu^{(0)}}{\mu^{(1)}}$ & \multirow{3}{*}{(kcal/mol pm$^{-2}$)}",
          "$k_A\\frac{\mu^{(0)}\mu^{(2)}}{\left(\mu^{(1)}\\right)^2}$ &",
          "$\\frac12\left\langle\partial_{\Delta r_\mathrm{HB}}^2c^{(2)}_\mathrm{P}\\right\\rangle$ &",
          "$\\frac16\left\langle\partial_{\Delta r_\mathrm{HB}}^3c^{(2)}\\right\\rangle$ & (kcal/mol pm$^{-3}$)",
          "$k_A$ & (kcal/mol pm$^{-2}$)",
          "$3/2\\alpha k_A$ & \multirow{2}{*}{(kcal/mol pm$^{-3}$)}",
          "$ k_A\\frac{\mu^{(2)}}{\mu^{(1)}}$ & "]
postscript=["\\\\\hline","\\\\","\\\\","\\\\\hline",
            "\\\\\hline\hline","\\\\\hline","\\\\",
            "\\\\\hline"]

str=preamble[0]
for j in range(6):
    ps = params[j][type]['x']
    md=ps[6]
    str+="& {:.4f} ".format(md*kas[j]*6.275)
str+=postscript[0]+'\n'
str+=preamble[1]
for j in range(6):
    ps = params[j][type]['x']
    md=ps[6]
    str+="& {:.4f} ".format(3/2*ass[j]*kas[j]*md*6.275/100)
str+=postscript[1]+'\n'
str+=preamble[2]
for j in range(6):
    ps = params[j][type]['x']
    md=ps[6]
    md2=ps[7]
    str+="& {:.4f} ".format(kas[j]*md*md2*6.275/100)
str+=postscript[2]+'\n'
str+=preamble[3]
for j in range(6):
    if np.abs(qs[j])<1e-4:
        str+="& 0.0000 ".format()
    else:
        str+="& {:.4f} ".format(qs[j]/2*6.275/100)
str+=postscript[3]+'\n'
str+=preamble[4]
for j in range(6):
    str+="& {:.4f} ".format(ts[j]/6*6.275/100/100)
str+=postscript[4]+'\n'
str+=preamble[5]
for j in range(6):
    str+="& {:.4f} ".format(kas[j]*6.275/100)
str+=postscript[5]+'\n'
str+=preamble[6]
for j in range(6):
    str+="& {:.4f} ".format(3/2*ass[j]*kas[j]*6.275/100/100)
str+=postscript[6]+'\n'
str+=preamble[7]
for j in range(6):
    ps = params[j][type]['x']
    md=ps[6]
    md2=ps[7]
    str+="& {:.4f} ".format(kas[j]*md2*6.275/100/100)
str+=postscript[7]+'\n'

print(str)