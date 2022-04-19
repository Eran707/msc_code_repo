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
import numpy as np
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
        self.extend_t = 0
        self.oldSim_timing_dict = {}
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

        if self.amend_type == "Extend":
            self._extend_sim()
        elif self.amend_type == "Resume":
            self._resume_sim()
        elif self.amend_type == "LastValues":
            self._last_values_sim()

    def _extend_sim(self):

        # Open old HDF file and copy all data into new file
        with h5py.File(self.file_name, mode='w') as hdf_new:
            try:
                with h5py.File(self.old_file_name, mode='r') as hdf_old:
                    groups = list(hdf_old.keys())
                    for _ in groups:
                        group = hdf_old.get(_)
                        hdf_old.copy(group, hdf_new)
                hdf_old.close()
            except:
                raise "could not open old file and copy data- check if it is in the current directory"

        # Open new file - find timing settings
        with h5py.File(self.file_name, mode='r') as hdf_new:

            T = hdf_new.get('TIMING')
            oldSim_total_t = T.get('TOTAL_T')[()]
            oldSim_dt = T.get("DT")[()]
            oldSim_intervals = T.get('INTERVALS')[()]
            oldSim_total_steps = oldSim_total_t / oldSim_dt
            oldSim_interval_step = oldSim_total_steps / oldSim_intervals
            oldSim_interval_arr = [round(oldSim_interval_step * i) for i in range(oldSim_intervals)]

            comps = hdf_new.get(groups[0])
            comp_list = list(comps.keys())
            # Use timing settings to find last interval and the interval steps
            # only find the second last interval for the first compartment
            comp = comps.get(comp_list[0])
            comp_intervals_list = list(comp.keys())
            oldSim_last_interval = oldSim_interval_arr[len(comp_intervals_list) - 2]

            self.oldSim_timing_dict = {"total_t": oldSim_total_t, "dt": oldSim_dt, "intervals": oldSim_intervals,
                                       "interval_step": oldSim_interval_step, "last_interval": oldSim_last_interval}
            print(self.oldSim_timing_dict)

            # find the last value for each compartment and initialize all compartments:

            for _ in comp_list:
                hdfcomp = comps.get(_)
                dataset_last = list(hdfcomp.get(str(oldSim_last_interval)))
                rad, length = dataset_last[1:3]
                na, k, cl, x, z = dataset_last[4:9]
                dataset_start = list(hdfcomp.get("0"))
                rad_start = dataset_start[1]
                sa_start = 2 * np.pi * rad_start * length
                comp_name = _
                comp = compartment.Compartment(compartment_name=comp_name, radius=rad, length=length,
                                               sa_value=sa_start, static_sa=True)
                comp.set_ion_properties(na_i=na, k_i=k, cl_i=cl, x_i=x, z_i=z, osmol_neutral_start=False)

                self.comp_arr.append(comp)

    def _resume_sim(self):
        print('code not complete')

    def _last_values_sim(self):
        print('code not complete')

    def add_compartment(self, comp=compartment.Compartment):
        print('code not complete')

    def set_timing(self, extend_t=30 ):

        if self.amend_type == "Extend":
            self.extend_t = extend_t
            self.dt = self.oldSim_timing_dict["dt"]
            self.run_t = self.oldSim_timing_dict["last_interval"] * self.oldSim_timing_dict["dt"]
            self.total_t = self.run_t + self.extend_t
            self.total_steps = self.total_t/self.dt
            self.intervals = round(self.total_steps/self.oldSim_timing_dict["interval_step"])
            self.interval_arr = [i*self.oldSim_timing_dict["interval_step"] for i in range(self.intervals)]




    def set_static_sa(self):
        print('code not complete')

    def run_simulation(self):

        while self.run_t < self.total_t:
            self.run_t += self.dt

            print('code not complete')

    def save_simulation_settings(self):
        print('code not complete')
