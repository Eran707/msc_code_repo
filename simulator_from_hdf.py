"""
This class functions to:

 A) "Extend" a previous simulation:
        Copies all data from a previous simulation HDF5 file into a new HDF file.
        The simulation will be setup based on the values stored in the last saved interval.
        The set timing function will be designed to add additional time to the total time of the base simulation.
 B) "Resume" a previous simulation:
        Used when a simulation has crashed before completion - often in the simulation of load shedding.
        Copies all the data from a previous simulation HDF file into a new HDF file
        Then runs the simulation to the total time described in the base simulation starting from the point of the crash
 C) "LastValues" of a simulation:
        Used when one wants to start a new simulation off the back of just the last values (and not entire dataset) of
        another experiment. This function will just copy the setup from the previous simulation and make a new file that
        starts with those initial conditions. This is useful when wanting to run simulations of another simulation
        that has already reached a steady state.


The controller_from_hdf.py file is used to interface with the SimulatorFromHDF class

Creator: EFS
Date:  April 2022
"""

import compartment
import common
import h5py
import numpy
import electrodiffusion


class SimulatorFromHDF:

    def __init__(self, old_file_name, new_file_name, amend_type='Extend'):
        """
        @param old_file_name: hdf5 file which must be in the same base directory as this code
        @param new_file_name: name of new hdf5 which will be saved in C
        @param amend_type: "Extend" or "Resume" or "LastValues"
        """
        # HDF file names and amendment type
        self.file_name = "C:\ " + new_file_name
        self.old_file_name = old_file_name
        self.amend_type = amend_type
        # Compartment structure initialization
        self.comp_arr, self.num_comps = [], 0
        # Surface area and pump defaults initialization
        self.static_atpase, self.static_sa = True, True
        self.p = 0
        # Electrodiffusion initialization
        self.ED_on, self.ed_arr, self.diff_constants = True, [], {}
        # External ion properties initialization
        self.na_o, self.k_o, self.cl_o, self.x_o, self.z_o, self.osm_o = 0, 0, 0, 0, 0, 0
        # Timing and interval initialization
        self.start_t, self.end_t, self.run_t, self.total_t, self.dt = 0, 0, 0, 0, 0
        self.one_percent_t, self.intervals, self.interval_num, self.steps = 0, 0, 1, 0
        # X-flux initialization
        self.xflux_setup, self.xflux_dict, self.xflux_count = True, {}, 0
        # z-flux initialization
        self.zflux_setup, self.sim_zflux_params = True, {"On": False}
        # Current initialization
        self.sim_current_params = {"On": False}
        # Synapse initialization
        self.syn_dict = {"On": False}
        # Hodgkin-Huxley initialization :
        self.hh_on, self.hh_comp_num, self.hh_t_on = False, 0, 0

        self.amend_sim()

    def amend_sim(self):

        if self.amend_type == "Extend":
            # Open old HDF file
            with h5py.File(self.file_name, mode='w') as hdf_new:
                try:
                    with h5py.File(self.old_file_name, mode='r') as hdf_old:
                        groups = list(hdf_old.keys())
                        for _ in groups:
                            hdf_old.copy(_, hdf_new)
                            hdf_old.close()
                except:
                    raise "could not open old file and copy data- check if it is in the current directory"

            # Copy all the data into new file





        elif self.amend_type == "Resume":
            print('code not complete')
        elif self.amend_type == "LastValues":
            print('code not complete')
