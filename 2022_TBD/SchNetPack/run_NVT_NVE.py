import os, sys
import schnetpack as spk
from schnetpack.md import System
import torch

import numpy as np
from schnetpack.md.utils import HDF5Loader
from schnetpack import Properties
from schnetpack.md.utils.md_units import MDUnits


T = int(sys.argv[1])
phase = sys.argv[2]
data_frac = 1
model = 'allphases'
supercel = sys.argv[3]


def get_NVT_data(log_file, T):
    data = HDF5Loader(log_file)
    temperatures = data.get_temperature()
    temperature = temperatures[-1]
    pos = data.get_positions()[-1]
    return temperature, pos

def run_NVT(steps):

    #
    # Initialization
    #

    model_path = '../../../training/final/%s_10fs_rcut_6_128_50_6_0.0001_%d/best_model'%(model, data_frac)
    structure = 'snapshots/%s_%dK_%s.xyz'%(phase, T, supercel)

    md_device = 'cuda'
    n_replicas = 1

    # Initialize the system
    md_system = System(n_replicas, device=md_device)
    md_system.load_molecules_from_xyz(structure)

    from schnetpack.md import MaxwellBoltzmannInit
    md_initializer = MaxwellBoltzmannInit(T, remove_translation=True, remove_rotation=True)
    md_initializer.initialize_system(md_system)

    #
    # MD integrator & calculator
    #

    from schnetpack.md.integrators import VelocityVerlet
    from schnetpack.md.calculators import SchnetPackCalculator
    from schnetpack import Properties
    from schnetpack.md.neighbor_lists import TorchNeighborList, ASENeighborList

    # Setup the integrator
    md_integrator = VelocityVerlet(2, device=md_device)

    # Load the stored model
    md_model = torch.load(model_path, map_location=md_device).to(md_device)

    # Generate the calculator
    md_calculator = SchnetPackCalculator(
        md_model,
        required_properties=[Properties.energy, Properties.forces],
        force_handle=Properties.forces,
        position_conversion='A',
        force_conversion='eV/A',
        neighbor_list=ASENeighborList,
        cutoff_shell=3.0
    )

    # Thermostat
    from schnetpack.md.simulation_hooks import thermostats

    # Set temperature and thermostat constant
    bath_temperature = T # K
    time_constant = 100 # fs

    # Initialize the thermostat
    langevin = thermostats.LangevinThermostat(bath_temperature, time_constant)


    #
    # Logging
    #
    from schnetpack.md.simulation_hooks import logging_hooks

    log_file = 'results/phase_%s_model_%s_data_%d_%dK_%s.hdf5'%(phase, model, data_frac, T, supercel)
    buffer_size = 1

    # Set up data streams to store positions, momenta and all properties
    data_streams = [
        logging_hooks.MoleculeStream(),
        logging_hooks.PropertyStream()
    ]

    # Create the file logger
    file_logger = logging_hooks.FileLogger(
        log_file,
        buffer_size,
        data_streams=data_streams,
        every_n_steps=10
    )

    simulation_hooks = [
        langevin,
        file_logger
    ]

    #
    # Run MD
    #

    from schnetpack.md import Simulator

    n_steps = steps
    md_simulator = Simulator(md_system, md_integrator, md_calculator, simulator_hooks=simulation_hooks)
    md_simulator.simulate(n_steps)

    return log_file


def setup_system(system, pos):
    pos = torch.from_numpy(pos).float().to('cuda')
    system.positions[0, 0, :, :] = pos


def run_NVE(log_file, steps):

    temperature, pos = get_NVT_data(log_file, T)
    os.remove(log_file)

    #
    # Initialization
    #

    model_path = '../../../training/final/%s_10fs_rcut_6_128_50_6_0.0001_%d/best_model'%(model, data_frac)
    structure = 'snapshots/%s_%dK_%s.xyz'%(phase, T, supercel)

    md_device = 'cuda'
    n_replicas = 1

    # Initialize the system
    md_system = System(n_replicas, device=md_device)
    md_system.load_molecules_from_xyz(structure)

    from schnetpack.md import MaxwellBoltzmannInit
    md_initializer = MaxwellBoltzmannInit(T, remove_translation=True, remove_rotation=True)
    md_initializer.initialize_system(md_system)

    setup_system(md_system, pos)

    #
    # MD integrator & calculator
    #

    from schnetpack.md.integrators import VelocityVerlet
    from schnetpack.md.calculators import SchnetPackCalculator
    from schnetpack import Properties
    from schnetpack.md.neighbor_lists import TorchNeighborList, ASENeighborList

    # Setup the integrator
    md_integrator = VelocityVerlet(2, device=md_device)

    # Load the stored model
    md_model = torch.load(model_path, map_location=md_device).to(md_device)

    # Generate the calculator
    md_calculator = SchnetPackCalculator(
        md_model,
        required_properties=[Properties.energy, Properties.forces],
        force_handle=Properties.forces,
        position_conversion='A',
        force_conversion='eV/A',
        neighbor_list=ASENeighborList,
        cutoff_shell=3.0
    )

    #
    # Logging
    #
    from schnetpack.md.simulation_hooks import logging_hooks

    log_file = 'results/phase_%s_model_%s_data_%d_%dK_%s.hdf5'%(phase, model, data_frac, T, supercel)
    buffer_size = 1

    # Set up data streams to store positions, momenta and all properties
    data_streams = [
        logging_hooks.MoleculeStream(),
        logging_hooks.PropertyStream()
    ]

    # Create the file logger
    file_logger = logging_hooks.FileLogger(
        log_file,
        buffer_size,
        data_streams=data_streams,
        every_n_steps=10
    )

    simulation_hooks = [
        file_logger
    ]

    #
    # Run MD
    #

    from schnetpack.md import Simulator

    n_steps = steps
    md_simulator = Simulator(md_system, md_integrator, md_calculator, simulator_hooks=simulation_hooks)
    md_simulator.simulate(n_steps)


# Run inital NVT simulation
log_file = run_NVT(50000)
# Extract last frame and run NVE simulation
run_NVE(log_file, 500000)

