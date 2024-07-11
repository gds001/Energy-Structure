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
 H3                   1    Rah      2    90.0
 O4                   1    Rab      2     90.      3     0.       0
 H5                   4    R1       1    A1       2    D1       0
 H6                   4    R1       1    A1       2    -D1      0
       Variables:
  R1                 {:.4f}                      
  A1                 {:.4f}                      
  D1                 {:.4f}         
  Rah                {:.4f}    
                    
  Rab                {:.4f} 
"""
header = """%nproc=28
%mem=50gb
# mp2/aug-cc-pVTZ pop=ChelpG opt=(z-matrix,MaxCycles=100,CalcFC) SCF=(MaxCycles=10000)

Title: Rab = {:.4f} Optimizations

0 1"""
linker = """
--Link1--
"""

file = "hf_h2o_linear_mp2_avtz"
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
file = "hf_h2o_linear_mp2_avtz"
i=0
rahmon=0
while True:
    line = lines[i]
    if line[:len(endstring)] == endstring: break
    if line[:len(r1string)] == r1string:
        rahmon = float(line.split()[2])
    i += 1

r2s = np.concatenate((np.linspace(-0.5, -0.1, 5),
                      np.linspace(0, 0.9, 10),
                      np.linspace(1, 2.8, 10),
                      np.linspace(3, 8, 6),
                      np.linspace(10, 18, 5)))
r2s=np.sort(r2s)
r2s += rab
print(r2s)

string = ""
for rab in r2s:
    if rab-rah<min([rah]):
        print("Tight",rah,rab)
        continue
    string += header.format(rab)
    string += zmatrix.format(r1, a1,d1, rah, rab)
    string += linker

f = open(file + "_rab_opts.gjf", 'w')
f.write(string[:-10])
f.close()
