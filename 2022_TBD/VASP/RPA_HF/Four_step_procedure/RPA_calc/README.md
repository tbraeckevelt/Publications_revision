Here you can find the input settings of the actual RPA correlation calculation, setting NBANDS to the maximum number of plane-waves (as in step3) tells VASP that you do not want to use the all-in-one approach.
Makes sure that the WAVECAR (typically very heavy) of step 3 is present in the folder, if LOPTICS = .TRUE.  in step 3, then also the WAVEDER file from that calculation should be present in this folder.