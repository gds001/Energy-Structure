import numpy as np
import scipy.interpolate as inter

endstring = " Normal termination of Gaussian 16"
r1string = " !       R1        "
a1string = " !       A1      "
r2string = " !       R2  "
r3string = " !       R3  "
nmzstring = " !       nmt   "
a2string = " !       A2   "
a3string = " !       A3   "
d1string = " !       D1   "
rahstring = " !       Rah  "
rabstring = " !       Rab  "

zmatrix = """
 F1
 X2                   1     1. 
 H3                   1    Rah      2     90. 
 F4                   1    Rab      2     90.      3     0.       0
 X5                   4     1.       1     90.      2     90.      0
 H5                   4    R1       5     90.      1    D1       0
       Variables:
  R1                   {:.4f}                 
  D1                  {:.4f}          
  Rah                   {:.4f}      
  Rab                   {:.4f}       
"""
header = """%nproc=40
%mem=50gb
# mp2/aug-cc-pVTZ Pop=ChelpG

Title: Rah = {:.4f} Rab = {:.4f} Scans

0 1"""
linker = """
--Link1--
"""

file = "hf_hf_linear_mp2_avtz"
f = open(file + '_rab_opts.log', 'r')
lines = f.readlines()
f.close()

r1 = 0
rah = 0
rab = 0
a1 = 0
a2 = 0
r2 = 0
r3 = 0
a3 = 0
d1 = 0
nmz = 0
i = 0
rabs=[]
r1s = []
a1s = []
a2s = []
r2s = []
r3s = []
a3s = []
d1s = []
nmzs = []
while i<len(lines):
    line = lines[i]
    if line[:len(endstring)] == endstring:
        rabs.append(rab)
        r1s.append(r1)
        a1s.append(a1)
        a2s.append(a2)
        r2s.append(r2)
        r3s.append(r3)
        a3s.append(a3)
        d1s.append(d1)
        nmzs.append(nmz)
    if line[:len(r1string)] == r1string:
        r1 = float(line.split()[2])
    if line[:len(r2string)] == r2string:
        r2 = float(line.split()[2])
    if line[:len(r3string)] == r3string:
        r3 = float(line.split()[2])
    if line[:len(rahstring)] == rahstring:
        rah = float(line.split()[2])
    if line[:len(rabstring)] == rabstring:
        rab = float(line.split()[2])
    if line[:len(a1string)] == a1string:
        a1 = float(line.split()[2])
    if line[:len(a2string)] == a2string:
        a2 = float(line.split()[2])
    if line[:len(a3string)] == a3string:
        a3 = float(line.split()[2])
    if line[:len(nmzstring)] == nmzstring:
        nmz = float(line.split()[2])
    if line[:len(d1string)] == d1string:
        d1 = float(line.split()[2])
    i += 1

file = "../hf/hf_mp2_avtz"
f = open(file + '_opt_freq.log', 'r')
lines = f.readlines()
f.close()
file = "hf_hf_linear_mp2_avtz"
i=0
rahmon=0
while True:
    line = lines[i]
    if line[:len(endstring)] == endstring: break
    if line[:len(r1string)] == r1string:
        rahmon = float(line.split()[2])
    i += 1

rahs=np.linspace(-0.1,0.1,21)
rahs=np.sort(rahs)
rahs += rahmon
print(rahs)
print(rabs)
string = ""
for k in range(len(rabs)):
    for rah in rahs:
        if rabs[k] - rah < min(rahs):
            print("Tight", rah, rabs[k])
            continue
        string += header.format(rah,rabs[k])
        string += zmatrix.format(r1s[k], d1s[k], rah, rabs[k])
        string += linker

f = open(file + "_rahs_rab.gjf", 'w')
f.write(string[:-10])
f.close()


