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
rahstring = " !       Roh  "
rabstring = " !       Roo  "

zmatrix = """
 O
 X                    1     1. 
 H                    1    Roh      2    A0
 O                    1    Roo      2     90.      3     0.       0
 H                    1    R1       3    A1       2     0.       0
 H                    4    R2       1    A2       5    D1       0
 H                    4    R2       1    A2       5    -D1      0
       Variables:                     
  R1                  {:.4f}                       
  A1                  {:.4f}                       
  R2                  {:.4f}                     
  A2                 {:.4f}        
  D1                  {:.4f}  
  Roh                {:.4f} 
  
  Roo               {:.4f}     
  A0                 {:.4f}
"""
header = """%nproc=40
%mem=50gb
# mp2/aug-cc-pVTZ Pop=ChelpG opt=(z-matrix,MaxCycles=1000)

Title: Rah = {:.4f} Rab = {:.4f}

0 1"""
linker = """
--Link1--
"""

file = "h2o_h2o_linear_mp2_avtz"
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

file = "../h2o/h2o_mp2_avtz"
f = open(file + '_opt_freq.log', 'r')
lines = f.readlines()
f.close()
file = "h2o_h2o_linear_mp2_avtz"
i=0
rahmon=0
while True:
    line = lines[i]
    if line[:len(endstring)] == endstring: break
    if line[:len(r1string)] == r1string:
        rahmon = float(line.split()[2])
    i += 1

aos=np.linspace(-10,10,10)
r2s = np.concatenate((np.linspace(-0.5, -0.1, 5),
                      np.linspace(0, 0.9, 10),
                      np.linspace(1, 2.8, 10),
                      np.linspace(3, 8, 6),
                      np.linspace(10, 18, 5)))
r2s=np.sort(r2s)
r2s += rab
print(r2s)
aos=[90.0]
string = ""
for ao in aos:
    for rab in r2s:
        if rab-rah<rah:
            print("Tight",rah,rab)
            continue
        string += header.format(rab,ao)
        string += zmatrix.format(r1,a1,r2,a2,d1,rah,rab, ao)
        string += linker

f = open(file + "_rab_opts.gjf", 'w')
f.write(string[:-10])
f.close()

