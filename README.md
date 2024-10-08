# Energy-Structure-Spectra

This repository is to both store output files associated with the work <BLANK>
as well as provide a base of code so someone can expend this work to new potentials
and new molecuels.

Do you want to do some analysis on your own set of molecules or with your own potential? 
Well lucky here we can do that easily for you.

## Make New Molecule
1) Create appropriate monomer and dimer geometries, and corresponding folders
2) Run 1-D scan of A-H bond in monomer.  This is used to make monomer potential
3) Run 1-D relaxed scan of A-B heavy atom distance in linear dimer. Hold symmetry and angle.
4) Run 2-D rigid scan of A-H at each A-B geometry.
5) Run 2-D relaxed scan of A-B and H-A-B angle, with frequencies.  Deuterate all non-hydrogen bonding hydrogens.
6) Export all 2-D scans to .csv

## Make New Potential Surfaces
1) If coding up a new potential function, make sure all is coded.  Code up c coefficiets for later use.
2) If fitting to a new molecule, begin making a fitting procedure
   1) Make sure to ensure stability in your potential and fitting procedure.
   2) Utilize repetative optimizations (slowly relaxing variables) to speed up high dimensional fittings.
   3) Use a global searching method, such as basin hopping.
   4) If confused, please look at codes used for these Buckinghamm-like potentials.
3) If you have a new linear molecule, make sure to work up c coefficiets via. finite difference too.

## Analize and Fit Energy vs. Structure
This can be done in a large number of ways; the methods done in the published work focueses on a
Taylor expansion around R to infinity. If the Taylor expansion does not diverge at this point, it's
a good solution.  If it does diverge, another technique should be explored.  Below shows
how best go do this Taylor expansion.

1) Explore derivatives of c0 coefficient.
   1) Simply plot first, second, third, and more if need be.
   2) Determine where derivative limits vanish to zero (typically second).
   3) Add "derivative corrections" to describe non-zero higher derivatves for
      stronger hydrogen bonds.  This is done for the third derivative by fitting the second 
      derivative with a line, but could be done in a number of ways, to a number of degrees.
2) Explore derivatives of c1 coefficient.
   1) Limit derivatives appeared sufficient in this work, but it woudl be good for you to check especially if
      if you changed the potential type.
3) Include and simplify the derivatives of c0 and c1 in eqn. 9.  Plot and check these answers
   to make sure things work out.

## Descriptions of the files in this repository follows...

### 1) Outputs
A repo with all output files, associated input files, csvs that contain the important
information for scans in a faster-to-read style, and python scripts that both write
scans and read scans.

#### a) monomers
Monomers are scanned in order to get good information regarding the convalent potential
of the hydrogen bond donor.  In each one, there are a set of files

###### \*opt_freq.\*
Simple optimization file to get the minimum geometry of each molecule

###### \*scan.\*
Simple rigid (other degrees of freedom relaxed) 1-D scan of the X-H bond.

###### Make\*Scan.py
Python file that takes the optimized geometry and creates the rigid 1-D scan for the
X-H bond.  

###### ReadScan.py
Python file taht reads a Link1 concatenated scan and reads the scanning parameters 
into a .csv.  Here, only energy and structure are saved to .csv.

#### b) dimers
Dimers are scanned twice to get a lot of information to get training information on 
the linear hydrogen bond and test information on the non-linear hydrogen bond.

###### \*opt_freq.\*
The optimization of the linear hydrogen bond

###### \*rab_opts.\*
Semi-relaxed scan of heavy atom distance (hydrogen bond angle and Cs symmetry held)

###### \*rahs_rab.\*
Rigid scan of rah at each of the geometries from the rab scan

###### \*rab_aahb_opts*/
Semi-relaxed scan of hydrogen bond angle and heavy atom distance (Cs symmetry held).
This may be separated across many files. Subsequent files are constructed via. copy and 
paste, then removing the jobs that finished.

###### Generate Linear Structues.py
Python script to write the linear semi-rigid scans of rab.

###### Generate Nonlinear Structures.py
Python script to write the non-linear semi-rigid scans of rab and aahb

###### Scanr.py
Python script to take the linear semi-rigid scan and add a rah scan to each geometry

###### ReadScan2d.py
Python script to read the 2-D scan of the linear hydrogen bond in rah and rab

###### ReadScanOptScans.py
Python script to read the 2-D scan of the non-linear hydrogen bonds in rah and aahb.

### 2) Fitted_Potentials

The repo should contain a folder for all the different types of potentials explored in this work
as well as files to fit the linear hydrogen bonds to the fitted potentials and files that contain 
fitted parameters for the potetial.

#### a) general_buckingham_potential
This potential is a buckingham potential with dipole and polarizability expanded up to cubic.
This can be truncated to shorter order by bounding that dipole derivative to zero.

###### ExpansionsAndDerivatives.py
This file contains the function to compute the potential.  It should also contain functions
for each c coefficients, derivatives of the c coefficinets with respect to R, derivatives of things
with respect to rHB.

###### 2DFittingPolar*.py
This file contains code to fit 