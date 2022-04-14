Here you can the input scripts for SOC calculations using the PBE XC functional.
INCAR_pre was used to find a good initial guess for the charge density and wavefunction.
INCAR_SOC performs the actual PBE+SOC calculation.
No dispersion method was used, but that is also possible via IVDW tag similar to the regular PBE calculations.

For the actual SOC calculation (not the pre-calculation), non-collinear calculations with VASP version 5.4.4 was performed.

ICHARG = 11 is recommended by VASP, but could lead to quite different results compared to ICHARG = 1, 
if ICHARG = 1 does not convergence after a lot of steps you could consider to restart after a few ICHARG = 1 steps, 
using that WAVECAR and ICHARG = 11 could already lead to better results instead of directly setting ICHARG = 11.
