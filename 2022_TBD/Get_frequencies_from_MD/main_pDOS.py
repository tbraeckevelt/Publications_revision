import molmod
import numpy as np

from utils_pDOS import PhononDensity, Create_slice_h5_file_from_CP2K_vel_xyz, initialize_plot, write_out_plot

unit_THz = molmod.units.hertz * 10**12
unit_kjmol = molmod.units.kjmol

if __name__ == '__main__':

    #fig, axis = initialize_plot()
    
    phase = "gamma"
    tem = 300
    
    """File names"""
    xyz_file_name = 'traj_vel_'+phase+'_T'+str(tem)+'.xyz'
    h5_file_name = 'snap_'+phase+'_T'+str(tem)+'.h5'
    p_file_name = 'ps_'+phase+'_T'+str(tem)+'.p'
    
    """create h5 file from CP2K xyz velocity file"""
    dt_file_au = 2 * molmod.units.femtosecond
    Create_slice_h5_file_from_CP2K_vel_xyz(xyz_file_name, h5_file_name, dt_file_au, step = 10)
    
    """create phonon density object from h5 file"""
    pdos = PhononDensity.from_file(h5_file_name, bsize = 1200)
    pdos.calc_free_energy_and_entropy(tem)
    
    
    """write out or import phonon density object to or from pickle file"""
    pdos.write_pickle_file(p_file_name)
    pdos = PhononDensity.from_file(p_file_name)
    
    """print and plot properties"""
    print(pdos.free_energy/ unit_kjmol/64)
    print(pdos.temperature)
    print(pdos.entropy/ unit_kjmol/64)
    pdos.plot(path_plot = "PhononDos.pdf", nat_pfu = 5, xlim = 5)
    
    #write_out_plot(fig, axis, path_plot = "PhononDOS.pdf", xlim = 5)
