import numpy as np

name='hf_nh3_mp2_avtz_rahs_rab'
file=open(name+".log",'r')
lines=file.readlines()
file.close()

titleline=" Title:"
energyline=" E2 =   "
angleline=" !       nmt     "
qCharge=" ESP charges:"
saveline=" Normal termination of Gaussian 16"

rs=[]
rtemp=0
Rs=[]
Rtemp=[]
es=[]
etemp=0
nmts=[]
nmt=90.
qHs=[]
qH=1
qAs=[]
qA=0
i=0
while i<len(lines):
    line=lines[i]
    if line[:len(titleline)]==titleline:
        rtemp=line.split()[3]
        Rtemp=line.split()[6]
    if line[:len(energyline)]==energyline:
        etemp=line.split()[5]
        etemp=etemp.replace('D','E')
    if line[:len(angleline)]==angleline:
        nmt=line.split()[2]
    if line[:len(qCharge)]==qCharge:
        i+=2
        qF=lines[i].split()[2]
        i+=1
        qH=lines[i].split()[2]
        qA=float(qF)+float(qH)
    if line[:len(saveline)]==saveline:
        rs.append(rtemp)
        Rs.append(Rtemp)
        es.append(etemp)
        nmts.append(nmt)
        qHs.append(qH)
        qAs.append(qA)
    i+=1

file=open(name+".csv",'w')
for i in range(len(rs)):
    file.write('{},{},{},{},{},{}\n'.format(rs[i], Rs[i], nmts[i], es[i], qHs[i],qAs[i]))
file.close()

