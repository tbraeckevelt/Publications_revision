Here you can the input scripts for calculations using the M06_L XC functional.

VASP version 5.4.4 was used.

If convergence problems are encountered, it is recommended to use a preconverged wavefunction for the PBE functional. 
Furthermore, ALGO = A (conjugate gradient algorithm for orbitals) is often more stable than charge density mixing, in particular, if the system contains vacuum regions.
From personal experience I would say that convergence issues occur more often for M06_L than for SCAN, so maybe use as default ALGO = A and a preconverged wavefunction.


Meta-GGA calculations require POTCAR files that include information on the kinetic energy density of the core-electrons. 
To check whether a particular POTCAR contains this information, type:
    grep kinetic POTCAR
This should yield at least the following lines (for each element on the file):
    kinetic energy-density
    mkinetic energy-density pseudized
and for PAW datasets with partial core corrections:
    kinetic energy density (partial)

