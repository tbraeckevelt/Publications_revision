For all calculations the PBE PAW potentials were used

For all calculations a gamma centered kpoint grid was used, with:
gamma: 443 \
delta_Cs: 642 \
delta_FA: 434 \

We used VASP version 6.1.2

There are two different procedures implemented in VASP to calculate the RPA+HF energies: a all-in-one approach, and a four step procedure, more details can be found in the corresponding subfolder.
Normally the all-in-one approach was preferred, but sometimes the calculations hanged on the exact exchange calculations and for those calculations we opted for the four step procedure.
