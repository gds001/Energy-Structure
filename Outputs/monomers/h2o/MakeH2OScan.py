import numpy as np

endstring=" Normal termination of Gaussian 16"
r1string=" !       R1        "
a1string=" !       A1      "

zmatrix="""
O
H 1 R1
H 1 R2 2 A1

R2={:.4f}
A1={:.4f}
R1={:.4f}
"""
header="""%nproc=40
%mem=50gb
# mp2/aug-cc-pVTZ Pop=ChelpG

Title: R = {:.4f}

0 1"""
linker="""
--Link1--
"""

file="h2o_mp2_avtz"
f=open(file+'_opt_freq.log','r')
lines=f.readlines()
f.close()


r1=0
r2=0
a1=0
i=0
while True:
    line=lines[i]
    if line[:len(endstring)]==endstring: break
    if line[:len(r1string)]==r1string:
        r1=float(line.split()[2])
    if line[:len(a1string)]==a1string:
        a1=float(line.split()[2])
    i+=1

rs=np.linspace(-0.1,0.1,21)
rs=np.sort(rs)
rs+=r1
print(rs)



string=""
for r in rs:
    string+=header.format(r)
    string+=zmatrix.format(r1,a1,r)
    string+=linker

f=open(file+"_scan.gjf",'w')
f.write(string[:-10])
f.close()
