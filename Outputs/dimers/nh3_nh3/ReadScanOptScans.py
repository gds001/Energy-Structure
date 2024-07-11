import numpy as np
from Useful import mDynepA_to_EhpA2

donor="../nh3/nh3"
file=open(donor+".log",'r')
donor=file.readlines()
file.close()
acceptor="../nh3/nh3"
file=open(acceptor+".log",'r')
acceptor=file.readlines()
nameog='nh3_nh3_linear_mp2_avtz_rab_aahb_opts'
file.close()
titleline=" Title:"
energyline=" E2 =   "
reline=" !       Rah "
raline=' !       R1 '
khbline=" Frc consts  --"
freqline=" Frequencies --"
saveline=" Normal termination of Gaussian 16"
ra=0
ka=0
ea=0
wa=0
i=0
while i<len(donor):
    line=donor[i][:-1]
    if line[:len(energyline)]==energyline:
        ea=line.split()[5]
        ea=float(ea.replace('D','E'))
    if line[:len(raline)]==raline:
        ra=float(line.split()[2])
    if line[:len(khbline)]==khbline:
        ka=float(line.split()[5])*mDynepA_to_EhpA2
    if line[:len(freqline)] == freqline:
        wa = float(line.split()[4])
    i+=1
eb=0
i=0
while i<len(acceptor):
    line=acceptor[i][:-1]
    if line[:len(energyline)]==energyline:
        eb=line.split()[5]
        eb=float(eb.replace('D','E'))
    if line[:len(saveline)]==saveline:
        break
    i+=1


res=[]
retemp=0
khbs=[]
khbtemp=0
whbs=[]
whbtemp=0
Rs=[]
Atemp=0
As=[]
Rtemp=0
es=[]
etemp=0
nmts=[]
nmt=90.
k=1
while True:
    name = nameog
    if k!=1:
        name=nameog+"_{}".format(k)
    k+=1
    try:
        file=open(name+".log",'r')
        lines=file.readlines()
        file.close()
    except: break
    i=0
    while i<len(lines):
        line=lines[i][:-1]
        if line[:len(titleline)]==titleline:
            Rtemp=line.split()[3]
            Atemp=line.split()[6]
        if line[:len(energyline)]==energyline:
            etemp=line.split()[5]
            etemp=float(etemp.replace('D','E'))
        if line[:len(reline)]==reline:
            retemp=float(line.split()[2])
        if line[:len(khbline)]==khbline:
            khbtemp=float(line.split()[5])*mDynepA_to_EhpA2
        if line[:len(freqline)]==freqline:
            whbtemp=float(line.split()[4])
        if line[:len(saveline)]==saveline:
            i+=1
            try: line=lines[i][:-1]
            except: break
            if line[:len(" Link1:  Proceeding to internal")]==" Link1:  Proceeding to internal":
                continue
            Rs.append(Rtemp)
            As.append(Atemp)
            es.append(ea+eb-etemp)
            res.append(retemp-ra)
            khbs.append(1-khbtemp/ka)
            whbs.append(wa-whbtemp)
        i+=1

file=open(nameog+".csv",'w')
for i in range(len(Rs)):
    file.write('{},{},{},{},{},{}\n'.format(Rs[i],As[i],es[i],res[i],khbs[i],whbs[i]))
file.close()