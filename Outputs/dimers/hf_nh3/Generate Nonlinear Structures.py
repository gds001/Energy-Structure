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
 F
 X                    1     1. 
 H                    1    Rah      2     Aahb
 N                    1    Rab      2     90.      3    0.0        0
 H(Iso=2)                    4    R1       3    A1       2     0.       0
 H(Iso=2)                    4    R2       3    A2       2     D1     0
 H(Iso=2)                    4    R2       3    A2       2    -D1     0
       Variables:
  R1                   {:.4f}     
  R2                 {:.4f}              
  A1                  {:.4f}
  A2                   {:.4f}       
  Rah                   {:.4f}      
  D1                   120.0
               
  Rab                   {:.4f}     
  Aahb                  {:.4f}
"""
header = """%nproc=40
%mem=50gb
# mp2/aug-cc-pVTZ pop=ChelpG opt=(z-matrix,MaxCycles=1000,MaxStep=1,CalcFC) freq

Title: Rab = {:.4f} Aabh = {:.4f}

0 1"""
linker = """
--Link1--
"""

file = "hf_nh3_mp2_avtz"
f = open(file + '_opt_freq.log', 'r')
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
while True:
    line = lines[i]
    if line[:len(endstring)] == endstring: break
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
file = "hf_nh3_mp2_avtz"
rahmon = 0
while True:
    line = lines[i]
    if line[:len(endstring)] == endstring: break
    if line[:len(r1string)] == r1string:
        rahmon = float(line.split()[2])
    i += 1

r2s = np.concatenate((np.linspace(-0.5, -0.1, 5),
                      np.linspace(0, 0.9, 10),
                      np.linspace(1, 3, 11)))
r2s = np.sort(r2s)
r2s += rab
print(r2s)
aos=[-25,-20,-17.5,-15,-12.5,-10,-8,-6,-5,-4,-3,-2,-1,
     0,1,2,3,4,5,6,8,10,12.5,15,17.5,20,25]
string = ""
for ao in aos:
    for rab in r2s:
        if rab - rah < rah:
            print("Tight", rah, rab)
            continue
        string += header.format(rab,ao)
        string += zmatrix.format(r1,r1, a1,a1, rah, rab,ao+90)
        string += linker

f = open(file + "_rab_aahb_opts.gjf", 'w')
f.write(string[:-10])
f.close()

