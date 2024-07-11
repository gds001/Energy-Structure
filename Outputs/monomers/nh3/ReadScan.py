import numpy as np

name='nh3_mp2_avtz_scan'
file=open(name+".log",'r')
lines=file.readlines()
file.close()

titleline=" Title:"
energyline=" E2 =   "
saveline=" Normal termination of Gaussian 16"

rs=[]
rtemp=0
es=[]
etemp=0
i=0
for line in lines:
    line=line[:-1]
    if line[:len(titleline)]==titleline:
        rtemp=line.split()[3]
    if line[:len(energyline)]==energyline:
        etemp=line.split()[5]
        etemp=etemp.replace('D','E')
    if line[:len(saveline)]==saveline:
        rs.append(rtemp)
        es.append(etemp)
    i+=1

file=open(name+".csv",'w')
for i in range(len(rs)):
    file.write('{},{}\n'.format(rs[i],es[i]))
file.close()

