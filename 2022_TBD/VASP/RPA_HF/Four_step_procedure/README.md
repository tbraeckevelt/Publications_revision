Here you can find the input files for the four steps procedure to calculate the RPA+HF energy

The four steps are:
1. Inital DFT (PBE) run
2. HF using DFT orbitals
3. Calculations and exact diagonalization of unoccupied states from DFT WAVECAR
4. Calculations of the RPA correlation energy

Make sure the computing nodes have enough memory for steps 3 and 4 or use enough nodes, for CsPbI3 with these settings more than 1 TB was needed.
Due to the large number of unoccupied steps the WAVECAR file after step 3 can be quite large (e.g., for these setting  around 50 GB), make sure you have enough storage space.

The HF exchange energy is reported in the OUTCAR file of step 2 and the correlation energy is reported in the OUTCAR file of step 4 (the value proceded by: "converged value"). The total RPA+HF energy is the sum of both contributions.
