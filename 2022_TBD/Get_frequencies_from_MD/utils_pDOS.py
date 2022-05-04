import molmod
import yaff
import numpy as np
import pickle
import h5py
from ase.io import read
from ase import Atom

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

matplotlib.pyplot.rcParams['axes.axisbelow'] = True
matplotlib.pyplot.rcParams['font.family'] = 'sans-serif'
matplotlib.pyplot.rcParams['font.size'] = 14
matplotlib.pyplot.rcParams['axes.linewidth'] = 0.65

FIGURE_KWARGS = {
        'figsize': (5.5, 4.4),
        'tight_layout': False,
        }
PLOT_KWARGS = {
        #'marker': '.',
        #'markersize': 4,
        'linewidth': 1.5,
        'markeredgecolor': 'k',
        'markeredgewidth': 0.5,
        }

#unit_invcm = molmod.constants.lightspeed / molmod.units.centimeter
unit_THz = molmod.units.hertz * 10**12
unit_kjmol = molmod.units.kjmol

class PhononDensity():
    """Class to represent a phonon density of states"""

    def __init__(self, ndof, frequencies, spectrum, delta_f=None, bsize = None, remove_vel_cm = True, normed_to_ndof = True, par_sym_lst = None):
        assert len(frequencies.shape) == 1
        assert len(spectrum.shape) == 1
        self.frequencies = frequencies
        self.spectrum = spectrum
        self.delta_f = delta_f
        self.bsize = bsize
        self.ndof = ndof
        self.remove_vel_cm = remove_vel_cm
        self.normed_to_ndof = normed_to_ndof
        self.par_sym_lst = par_sym_lst
        self.temperature = None
        self.entropy = None
        self.free_energy = None

    def calc_free_energy_and_entropy(self, temperature):
        """Computes the entropy as a function of frequency"""
        
        self.temperature = temperature
        delta_f = self.frequencies[1]- self.frequencies[0]
        
        # compute entropy for each frequency and integrate
        free_energy_vib, entropy_vib  = _free_energy_and_entropy_quantum(self.frequencies, self.temperature)
        free_energy_spec = delta_f * self.spectrum * free_energy_vib
        entropy_spec = delta_f * self.spectrum * entropy_vib
        self.free_energy = np.sum(free_energy_spec)
        self.entropy = np.sum(entropy_spec)
        return self.free_energy, self.entropy

    def plot(self, fig = None, axis = None, path_plot = None, nat_pfu = None, xlim = None, ylim = None, inthz = True, label = "", lab_FE = True):
        """Plot the normalized spectrum and the entropy function"""
        
        if fig is None or axis is None:
            fig, axis = initialize_plot()
        
        if inthz:
            spectrum = self.spectrum * unit_THz
            frequencies = self.frequencies / unit_THz
        else:
            spectrum = self.spectrum
            frequencies = self.frequencies
        
        if nat_pfu is not None:
            number_of_fu = np.ceil(self.ndof/(3*nat_pfu))
            if nat_pfu == 1:
                number_of_fu+=1
            spectrum /= number_of_fu
            free_energy  = self.free_energy / number_of_fu
        else:
            free_energy  = self.free_energy
        
        if lab_FE:
            label += "F = " + str(np.round(free_energy / unit_kjmol, 2)) + " kJ/mol"
            
        axis.plot(frequencies, spectrum, label = label, **PLOT_KWARGS)
        
        if path_plot is not None:
            legend = True
            if label == "":
                legend = False
            write_out_plot(fig, axis, path_plot, xlim = xlim, ylim = ylim, inthz = inthz, legend = legend)

    @classmethod
    def from_file(cls, filepath, delta_f=None, bsize = None, remove_vel_cm = True, normed_to_ndof = True, par_sym_lst = None):
        if filepath[-3:]== '.h5':
            return cls._from_h5(filepath, delta_f = delta_f, bsize = bsize, remove_vel_cm = remove_vel_cm, normed_to_ndof = normed_to_ndof, par_sym_lst = par_sym_lst)
        elif filepath[-2:] == '.p':
            return cls._from_p(filepath, delta_f = delta_f, bsize = bsize, remove_vel_cm = remove_vel_cm, normed_to_ndof = normed_to_ndof, par_sym_lst = par_sym_lst)
        else:
            raise NotImplementedError('I do not know the extension of {}'.format(
                filepath))

    @classmethod
    def _from_h5(cls, filepath, delta_f=None, bsize = None, remove_vel_cm = True, normed_to_ndof = True, par_sym_lst = None):
        """Computes the phonon density from a trajectory of velocities

        By default, the FFTs are computed such that the frequency resolution
        (as specified by the delta_f argument) is at least 0.1 invcm. If
        delta_f is None, then a single FFT is used over the entire trajectory.

        When calculating the phonon density, there are a few options that you can set:
        bsize = None or (delta_f = None):     To specify the frequency sampling, in this work we always set bsize = 1200
        remove_vel_cm = True:                 To remove the center of mass velocity
        normed_to_ndof = True:                To norm the phonon density of states to the total number of degrees of freedom
        par_sym_lst = None:                   To create parital phonon densities, only include the phonon density contributions of the elements specified in this list
        """
        with h5py.File(filepath, 'r') as f:
            time_signal = np.array(list(f['trajectory']['vel'])) #time_signal[time, atom, component]
            time = np.array(list(f['trajectory']['time']))
            masses = np.array(list(f['system']['masses']))
            if par_sym_lst is not None:
                time_signal, masses = _get_partial_velocities_and_masses(time_signal, masses, par_sym_lst)

        sampling_period = time[1] - time[0]
        assert np.shape(time)[0] == len(time_signal)
        assert len(masses) == np.shape(time_signal)[1]
            
        if remove_vel_cm:
            tot_mass =  np.sum(masses)
            for i in range(len(time)):
                for j in range(3): 
                    time_signal[i,:,j] -= np.dot(masses,time_signal[i,:,j])/tot_mass #remove center of mass velocity
            ndof = time_signal.shape[1] * 3 - 3
        else:
            ndof = time_signal.shape[1] * 3

        if delta_f is not None:
            # determine block size
            size = _determine_blocksize(delta_f, sampling_period)
            if bsize is not None:
                assert bsize == size
        elif bsize is not None:
            size = bsize
        else:
            # use entire signal as single block
            size = time_signal.shape[0]
        n_blocks = time_signal.shape[0] // size
        assert n_blocks >= 1
        delta_f = _compute_frequency_resolution(size, sampling_period)

        # iterate over all blocks and compute spectrum
        spectra = None
        for i in range(n_blocks):
            start = i * size
            end = start + size
            frequencies, block_spectrum = _compute_block(
                    time_signal[start:end, :],
                    sampling_period,
                    masses,
                    remove_vel_cm
                    )
            if spectra is None:
                spectra = block_spectrum.copy()
            else:
                spectra += block_spectrum
        spectrum = spectra / (n_blocks)
        
        if normed_to_ndof:
            area_spectrum = np.trapz(spectrum, frequencies)
            spectrum *= ndof / area_spectrum
            
        return cls(ndof, frequencies, spectrum, delta_f = delta_f, bsize = bsize, remove_vel_cm = remove_vel_cm, normed_to_ndof = normed_to_ndof, par_sym_lst = par_sym_lst)

    @classmethod
    def _from_p(cls, filepath, delta_f=None, bsize = None, remove_vel_cm = True, normed_to_ndof = True, par_sym_lst = None):
        """Loads a pickled file and returns a ``PhononDensity`` instance"""
        pdos = pickle.load(open(filepath, 'rb'))
        assert isinstance(pdos, PhononDensity)
        if delta_f is not None:
            assert pdos.delta_f == delta_f 
        if bsize is not None:
            assert pdos.bsize == bsize
        assert pdos.remove_vel_cm == remove_vel_cm
        assert pdos.normed_to_ndof == normed_to_ndof
        assert pdos.par_sym_lst == par_sym_lst

        return pdos
        
    def write_pickle_file(self, filepath_out):
        pickle.dump(self, open(filepath_out, 'wb'))


def _free_energy_and_entropy_quantum(f, T):
    if (f > 0).all():
        h = molmod.constants.planck
        k = molmod.constants.boltzmann
        beta = 1 / (molmod.constants.boltzmann * T)
        q_quantum = np.exp(- (beta * h * f) / 2) / (1 - np.exp(- beta * h * f))
        f_quantum = - np.log(q_quantum) / beta
        s_quantum = -k * np.log(1 - np.exp(- beta * h * f)) + h * f / T * (np.exp(beta * h * f) - 1) ** (-1)
        return f_quantum, s_quantum
    else:
        raise ValueError('Entropy at 0Hz is infinite')

def _get_partial_velocities_and_masses(time_signal, masses, par_sym_lst):
    par_mass_lst=[]
    for sym in par_sym_lst:
        par_mass_lst.append(Atom(sym).mass/molmod.units.unified)
    mask = np.isin(masses, par_mass_lst)
    return time_signal[:,mask,:],  masses[mask]

def _compute_block(time_signal, sampling_period, masses, remove_vel_cm):
    """Computes the spectrum of the autocorrelation of a time signal

    This function computes the power spectrum of a time signal by computing
    the norm-squared DFT spectrum of the time signal using the FFT algorithm.

    Arguments
    ---------

    time_signal (ndarray of shape (nsteps, natoms, 3)):
        trajectory of samples of positions, velocities etc, in atomic units

    sampling_period (float):
        elapsed time between each of the samples, in atomic units.
    """
    nsteps = time_signal.shape[0]
    _all = np.reshape(time_signal, (nsteps, -1))
    N = 2 * nsteps - 1 # avoid periodic copies influencing the result
    out = np.fft.fft(_all, n=N, axis=0)
    power_spectrum = np.abs(out) ** 2
    
    masses_c = np.array([masses,masses,masses]).T
    masses_c_1 = np.reshape(masses_c, -1)
    
    spectrum = np.zeros(N)
    for freq in range(N):
        spectrum[freq] = np.dot(masses_c_1,power_spectrum[freq,:])

    # frequency axis of N-point DFT
    frequencies = 1 / N * np.arange(N) * (1 / sampling_period)

    # spectrum is even
    n = N // 2
    if remove_vel_cm:
        #print("removed spectrum(freq = 0): " + str(spectrum[0]))
        #print("Which is {:f} of the variance of the spectrum".format(spectrum[0]/np.var(spectrum)))
        return frequencies[1:n], spectrum[1:n]
    else:
        return frequencies[:n], spectrum[:n]


def _compute_frequency_resolution(nsteps, sampling_period):
    """Computes the resolution of the FFT spectrum"""
    fs = 1 / sampling_period
    return fs / (2 * nsteps - 1)


def _determine_blocksize(delta_f, sampling_period):
    """Partitions the time_signal into blocks to fix the frequency resolution

    The return value is an integer representing the required block size, as
    obtained through the relation:

        1 / sampling_period = (2n - 1) * delta_f

    Arguments
    ---------

    delta_f (float):
        desired frequency resolution, in atomic units.

    sampling_period (float):
        total time inbetween two consecutive samples, in atomic units

    """
    return int(np.ceil((1 / (sampling_period * delta_f) + 1) / 2))
    
def Create_slice_h5_file_from_CP2K_vel_xyz(fn, h5_fn, dt_file_au, begin=0, end=None, step=1):
    with h5py.File(h5_fn, 'w') as fout:

        traj = read(fn, index = slice (begin,end,step))
        
        vel_au = np.zeros([len(traj), traj[0].get_global_number_of_atoms(), 3])
        for i,atoms in enumerate(traj):
            vel_au[i,:,:] = atoms.get_positions()
            
        time = np.arange(len(traj)) * step * dt_file_au
        mass_au = traj[0].get_masses()/molmod.units.unified

        fout.create_dataset('system/masses',data = mass_au, dtype='float64')
        fout.create_dataset('trajectory/vel',data = vel_au, dtype='float64')
        fout.create_dataset('trajectory/time',data = time, dtype='float64')

        return fout
        
def initialize_plot():
    fig = matplotlib.figure.Figure(**FIGURE_KWARGS)
    axis = fig.add_axes([0, 0, 1, 1])
    
    return fig, axis
    
def write_out_plot(fig, axis, path_plot, xlim = None, ylim = None, inthz = True, legend = True):
    axis.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    axis.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    #axis.yaxis.set_ticklabels([])
    #axis.xaxis.set_ticklabels([])
    axis.get_yaxis().set_tick_params(which='both', direction='in')
    axis.get_xaxis().set_tick_params(which='both', direction='in')
    if xlim is None:
        axis.set_xlim(0, None)
    else: 
        axis.set_xlim([0, xlim])
    if ylim is None:
        axis.set_ylim(0, None)
    else:
        axis.set_ylim([0, ylim])
    if inthz is True:
        axis.set_xlabel('Frequency [THz]')
        axis.set_ylabel('Normalized phonon density [1/THz]')
    else:
        axis.set_xlabel('Frequency [a.u.]')
        axis.set_ylabel('Normalized phonon density [a.u]')
    if legend:
        axis.legend(loc="upper right")
    fig.savefig(path_plot, bbox_inches='tight')
    fig.clf()

def print_similarity_metric(pdos1, pdos2):

    assert (pdos1.frequencies == pdos2.frequencies).all()
    ad_spec = np.abs(pdos1.spectrum/pdos1.ndof - pdos2.spectrum/pdos2.ndof)/2

    difference_met = np.trapz(ad_spec,pdos1.frequencies)
    similarity_met = 1- difference_met
    print(similarity_met)
