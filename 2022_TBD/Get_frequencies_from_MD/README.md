in utils_pDOS.py you can find the PhononDensity class and all methods needed to go from CP2K velocities in xyz format to a h5 file to a PhononDensity instance. From which we can compute the entropy and free energy. Moreover functions to plot the phonon density and calculate the similarity between different PhononDensity instances, are also present in this file. \

When calculating the phonon density, there are a few options that you can set: \
bsize = None or (delta_f = None):     To specify the frequency sampling, in this work we always set bsize = 1200 \
remove_vel_cm = True:                 To remove the center of mass velocity \
normed_to_ndof = True:                To norm the phonon density of states to the total number of degrees of freedom \
par_sym_lst = None:                   To create parital phonon densities, only include the phonon density contributions of the elements specified in this list


main_pDOS.py shows how you could use these methods to calculate from a CP2K velocity output file or a h5 file the phonon density and plot it together with the free energy.
