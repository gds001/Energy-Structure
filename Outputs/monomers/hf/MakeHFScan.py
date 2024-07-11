import numpy as np

endstring=" Normal termination of Gaussian 16"
monomerstring=" !       R1        "
anglestring=" !       A1      "

zmatrix="""
F
H 1 R1

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

file="hf_mp2_avtz"
f=open(file+'_opt_freq.log','r')
lines=f.readlines()
f.close()


dist=0
ang=0
i=0
while True:
    line=lines[i]
    if line[:len(endstring)]==endstring: break
    if line[:len(monomerstring)]==monomerstring:
        dist=float(line.split()[2])
    if line[:len(anglestring)]==anglestring:
        ang=float(line.split()[2])
    i+=1

rs=np.linspace(-0.1,0.1,21)
rs=np.sort(rs)
print(rs)
rs+=dist
print(rs)

string=""
for r in rs:
    string+=header.format(r)
    string+=zmatrix.format(r)
    string+=linker

f=open(file+"_scan.gjf",'w')
f.write(string[:-10])
f.close()
