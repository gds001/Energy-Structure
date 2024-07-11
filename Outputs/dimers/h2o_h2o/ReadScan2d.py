import numpy as np

name='h2o_h2o_linear_mp2_avtz_rahs_rab'
file=open(name+".log",'r')
lines=file.readlines()
file.close()

titleline=" Title:"
energyline=" E2 =   "
angleline=" !       nmt     "
saveline=" Normal termination of Gaussian 16"

rs=[]
rtemp=0
Rs=[]
Rtemp=[]
es=[]
etemp=0
nmts=[]
nmt=90.
i=0
for line in lines:
    line=line[:-1]
    if line[:len(titleline)]==titleline:
        rtemp=line.split()[3]
        Rtemp=line.split()[6]
    if line[:len(energyline)]==energyline:
        etemp=line.split()[5]
        etemp=etemp.replace('D','E')
    if line[:len(angleline)]==angleline:
        nmt=line.split()[2]
    if line[:len(saveline)]==saveline:
        rs.append(rtemp)
        Rs.append(Rtemp)
        es.append(etemp)
        nmts.append(nmt)
    i+=1

file=open(name+".csv",'w')
for i in range(len(rs)):
    file.write('{},{},{},{}\n'.format(rs[i],Rs[i],nmts[i],es[i]))
file.close()

