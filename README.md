# Energy-Structure-Spectra

This repository is to both store output files associated with the work <BLANK>
as well as provide a base of code so someone can expend this work to new potentials
and new molecuels.

Descriptions of the files in this repository follows...

## 1) Outputs
A repo with all output files, associated input files, csvs that contain the important
information for scans in a faster-to-read style, and python scripts that both write
scans and read scans.

### a) monomers
Monomers are scanned in order to get good information regarding the convalent potential
of the hydrogen bond donor.  In each one, there are a set of files

##### \*opt_freq.\*
Simple optimization file to get the minimum geometry of each molecule

##### \*scan.\*
Simple rigid (other degrees of freedom relaxed) 1-D scan of the X-H bond.

##### Make\*Scan.py
Python file that takes the optimized geometry and creates the rigid 1-D scan for the
X-H bond.  

##### ReadScan.py
Python file taht reads a Link1 concatenated scan and reads the scanning parameters 
into a .csv.  Here, only energy and structure are saved to .csv.

### b) dimers
Dimers are scanned twice to get a lot of information to get training information on 
the linear hydrogen bond and test information on the non-linear hydrogen bond.

##### \*opt_freq.\*
The optimization of the linear hydrogen bond

##### \*rab_opts.\*
Semi-relaxed scan of heavy atom distance (hydrogen bond angle and Cs symmetry held)

##### \*rahs_rab.\*
Rigid scan of rah at each of the geometries from the rab scan

##### \*rab_aahb_opts*/
Semi-relaxed scan of hydrogen bond angle and heavy atom distance (Cs symmetry held).
This may be separated across many files. Subsequent files are constructed via. copy and 
paste, then removing the jobs that finished.

##### Generate Linear Structues.py
Python script to write the linear semi-rigid scans of rab.

##### Generate Nonlinear Structures.py
Python script to write the non-linear semi-rigid scans of rab and aahb

##### Scanr.py
Python script to take the linear semi-rigid scan and add a rah scan to each geometry

##### ReadScan2d.py
Python script to read the 2-D scan of the linear hydrogen bond in rah and rab

##### ReadScanOptScans.py
Python script to read the 2-D scan of the non-linear hydrogen bonds in rah and aahb.