Here you can find the input file for the all-in-one approach to calculate the RPA+HF energy

Make sure the computing nodes have enough memory or use enough nodes, for CsPbI3 with these settings more than 1 TB was needed.
Note that LOPTICS=.TRUE. (without LPEAD= .TRUE.) does not effect the energy in the all-in-one approach. However, in the four steps procedure you can use it to convergence the energy with respect to the kpoint density faster.

The total energy is reported in the OUTCAR file, in the paper we have reported the value proceded by: "HF+E_corr(extrapolated)="
