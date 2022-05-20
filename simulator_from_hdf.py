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
from common import F
import h5py
import numpy as np
import electrodiffusion
import time


class SimulatorFromHDF:

    def __init__(self, old_file_name, new_file_name, amend_type='Extend', already_extended = False):
        """
        @param old_file_name: hdf5 file which must be in the same base directory as this code
        @param new_file_name: name of new hdf5 which will be saved in C
        @param amend_type: "Extend" or "Resume" or "LastValues"
        @param already_extended" "True" if you are extending an already extended file
        """
        # HDF file names and amendment type
        self.file_name = "C:\ " + new_file_name
        self.old_file_name = old_file_name
        self.amend_type = amend_type
        self.already_extended = already_extended
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
        self.syn_arr = []
        self.syn_count =0
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
            if self.already_extended == False:
                oldSim_total_t = T.get('TOTAL_T')[()]
                oldSim_dt = T.get("DT")[()]
                oldSim_intervals = T.get('INTERVALS')[()]
                oldSim_interval_step = oldSim_total_steps / oldSim_intervals
            else:
                oldSim_total_t = T.get('EXTENDED_TOTAL_T')[()]
                oldSim_dt = T.get("EXTENDED_DT")[()]
                oldSim_intervals = T.get('EXTENDED_INTERVALS')[()]
                oldSim_interval_step = T.get('EXTENDED_INTERVAL_STEP')[()]

            oldSim_total_steps = oldSim_total_t / oldSim_dt

            oldSim_interval_arr = [round(oldSim_interval_step * i) for i in range(oldSim_intervals)]

            comps = hdf_new.get(groups[0])
            comp_list = list(comps.keys())
            # Use timing settings to find last interval and the interval steps
            # only find the second last interval for the first compartment
            comp = comps.get(comp_list[0])
            comp_intervals_list = list(comp.keys())
            oldSim_last_interval = len(comp_intervals_list) - 1
            oldSim_last_step = oldSim_interval_arr[oldSim_last_interval]
            self.oldSim_timing_dict = {"total_t": oldSim_total_t, "dt": oldSim_dt, "intervals": oldSim_intervals,
                                       "last_step": oldSim_last_step, "interval_step": oldSim_interval_step,
                                       "last_interval": oldSim_last_interval}
            print(self.oldSim_timing_dict)




            # find the last value for each compartment and initialize all compartments:

            for _ in comp_list:
                hdfcomp = comps.get(_)
                dataset_last = list(hdfcomp.get(str(oldSim_last_step)))
                rad, length = dataset_last[1:3]
                na, k, cl, x, z = dataset_last[4:9]
                dataset_start = list(hdfcomp.get("0"))
                rad_start = dataset_start[1]
                sa_start = 2 * np.pi * rad_start * length
                # RETRIEVE ATPASE SETTINGS
                try:
                    A = hdf_new.get('ATPASE_SETTINGS')
                    self.na_start = A.get("NA-START")[()]
                except:
                    self.na_start = dataset_start[4]
                comp_name = _
                comp = compartment.Compartment(compartment_name=comp_name, radius=rad, length=length,
                                               sa_value=sa_start, static_sa=True)
                comp.set_ion_properties(na_i=na, k_i=k, cl_i=cl, x_i=x, z_i=z, osmol_neutral_start=False)
                comp.na_i_start = self.na_start
                self.comp_arr.append(comp)

    def _resume_sim(self):

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
            try:
                oldSim_total_t = T.get('EXTENDED_TOTAL_T')[()]
            except:
                oldSim_total_t = T.get('TOTAL_T')[()]
            try:
                oldSim_dt = T.get("EXTENDED_DT")[()]
            except:
                oldSim_dt = T.get("DT")[()]
            try:
                oldSim_intervals = T.get('EXTENDED_INTERVALS')[()]
            except:
                oldSim_intervals = T.get('INTERVALS')[()]

            oldSim_total_steps = round(oldSim_total_t / oldSim_dt)
            oldSim_interval_step = round(oldSim_total_steps / oldSim_intervals)
            oldSim_interval_arr = [round(oldSim_interval_step * i) for i in range(oldSim_intervals)]

            comps = hdf_new.get(groups[0])
            comp_list = list(comps.keys())
            # Use timing settings to find last interval and the interval steps
            # only find the second last interval for the first compartment
            comp = comps.get(comp_list[-1])
            comp_intervals_list = list(comp.keys())
            oldSim_last_interval = len(comp_intervals_list) - 1
            oldSim_last_step = round(oldSim_interval_arr[oldSim_last_interval])
            self.oldSim_timing_dict = {"total_t": oldSim_total_t, "dt": oldSim_dt, "intervals": oldSim_intervals,
                                       "last_step": oldSim_last_step, "interval_step": oldSim_interval_step,
                                       "last_interval": oldSim_last_interval}
            print(self.oldSim_timing_dict)

            # find the last value for each compartment and initialize all compartments:
            try:
                S = hdf_new.get('SA_SETTINGS')
                sa_start_arr = S.get("SA-VALUES")[()]
            except:
                pass
            j=0
            for _ in comp_list:

                hdfcomp = comps.get(_)
                dataset_last = list(hdfcomp.get(str(oldSim_last_step)))
                rad, length = dataset_last[1:3]
                na, k, cl, x, z = dataset_last[4:9]
                dataset_start = list(hdfcomp.get("0"))
                rad_start = dataset_start[1]

                try:
                    A = hdf_new.get('ATPASE_SETTINGS')
                    self.na_start = A.get("NA-START")[()]
                except:
                    self.na_start = dataset_start[4]
                try:
                    sa_start = sa_start_arr[j]
                except:
                    sa_start = 2 * np.pi * rad_start * length


                comp_name = _
                comp = compartment.Compartment(compartment_name=comp_name, radius=rad, length=length,
                                               sa_value=sa_start, static_sa=True)
                comp.set_ion_properties(na_i=na, k_i=k, cl_i=cl, x_i=x, z_i=z, osmol_neutral_start=False)
                comp.na_i_start = self.na_start
                self.comp_arr.append(comp)
                j+=1

    def _last_values_sim(self):

        # Open old HDF file and copy only last values into new file

        comp_arr_lastvalues = [] #saves all the last values of the base experiment
        comp_arr_startvalues = [] #saves all the first values of the base experiment

        with h5py.File(self.old_file_name, mode='r') as hdf_old:
            C = hdf_old.get('COMPARTMENTS')
            C_group_arr = []
            comp_names_arr = list(C.keys())
            if self.already_extended == True:
                T = hdf_old.get('TIMING')
                prev_total_t = T.get('TOTAL_T')[()]
                prev_intervals = T.get('INTERVALS')[()]
                dt = T.get("EXTENDED_DT")[()]
                prev_total_steps = round(prev_total_t / dt)
                interval_step = prev_total_steps / prev_intervals
                total_t = T.get('EXTENDED_TOTAL_T')[()]
                intervals = T.get('EXTENDED_INTERVALS')[()]
                interval_arr = [int(interval_step * i) for i in range(intervals)]
            else:
                T = hdf_old.get('TIMING')
                total_t = T.get('TOTAL_T')[()]
                intervals = T.get('INTERVALS')[()]
                dt = T.get("DT")[()]
                total_steps = total_t / dt
                interval_step = total_steps / intervals
                interval_arr = [int(interval_step * i) for i in range(intervals)]

            try:
                S = hdf_old.get('SA_SETTINGS')
                sa_start_arr = S.get("SA-VALUES")[()]
            except:
                pass

            # Looping through old compartments and saving last dataset values


            for _ in comp_names_arr:
                comp = C.get(_)
                steplist = list(comp.keys())
                rev_interval_arr = interval_arr[::-1]
                for j in rev_interval_arr:
                    interval = j
                    if str(interval) in str(steplist):
                        last_time_point = str(interval)
                        break
                last_dataset = comp.get(last_time_point)
                comp_arr_lastvalues.append(list(last_dataset))
                first_time_point = str(interval_arr[0])
                first_dataset = comp.get(first_time_point)
                comp_arr_startvalues.append(list(first_dataset))

            try:
                A = hdf_old.get('ATPASE_SETTINGS')
                self.na_start = A.get("NA-START")[()]
            except:
                self.na_start = first_dataset[4]


        hdf_old.close()

        # writing last values to the new file
        comp_arr = []  # array of compartment objects
        with h5py.File(self.file_name, mode='w') as hdf_new:
            hdf_new.create_group('COMPARTMENTS')
            C = hdf_new.get("COMPARTMENTS")
            for i in range(len(comp_names_arr)):
                C.create_group(comp_names_arr[i])
                Comp = C.get(comp_names_arr[i])
                last_dataset = comp_arr_lastvalues[i]
                last_dataset[0] = 0
                Comp.create_dataset(name='0',data=last_dataset)
                ## creating a new compartment object and adding to the simulator
                rad = last_dataset[1]
                length = last_dataset[2]
                try:
                    initial_sa =sa_start_arr[i]
                except:
                    initial_sa = 2 * np.pi * comp_arr_startvalues[i][1] * comp_arr_startvalues[i][2]
                comp = compartment.Compartment(compartment_name=comp_names_arr[i],radius=rad,length=length,static_sa=True, sa_value=initial_sa)
                comp.set_ion_properties(na_i=last_dataset[4],k_i=last_dataset[5],cl_i=last_dataset[6],
                                        x_i=last_dataset[7],z_i=last_dataset[8],osmol_neutral_start=False)
                comp.vm = last_dataset[-3]
                comp.na_i_start = self.na_start
                self.comp_arr.append(comp)
        hdf_new.close()


    def set_electrodiffusion_properties(self, ED_on: bool = True,
                                    diff_constant_dict=None):
        """
        Function to define the electrodiffusion setup at the beginning of a simulation.
        @param ED_on: boolean, sets either electrodiffusion between compartments to be on or off.
        @param diff_constant_dict: dictionary, diffusion constants in the format {"ion":value,"ion":value,...}

        """
        if diff_constant_dict is None:
            diff_constant_dict = {"na": (1.33 / 2) * 1e-7, "k": (1.96 / 2) * 1e-7,
                                  "cl": (2.03 / 2) * 1e-7, "x": 0}
        self.ED_on = ED_on
        self.diff_constants = diff_constant_dict

        if self.ED_on:
            for e in range(len(self.comp_arr) - 1):
                name = self.comp_arr[e].name + ' -> ' + self.comp_arr[e + 1].name

                ed = electrodiffusion.Electrodiffusion(comp_a_name=self.comp_arr[e].name,
                                                       comp_a_length=self.comp_arr[e].length,
                                                       comp_b_name=self.comp_arr[e + 1].name,
                                                       comp_b_length=self.comp_arr[e + 1].length)

                ed.set_diff_constants(self.diff_constants)
                self.ed_arr.append(ed)
        if self.amend_type =='LastValues':
            with h5py.File(self.file_name, mode='a') as self.hdf:
                ED = self.hdf.create_group("ELECTRODIFFUSION")
                for _ in range(len(self.ed_arr)):
                    ED.create_group(self.comp_arr[_].name + ' -> ' + self.comp_arr[_+1].name)

    def set_timing(self, extend_t=None,dt=None,total_t=None,intervals=None):
        if self.amend_type == "Extend":
            self.extend_t = extend_t
            self.dt = self.oldSim_timing_dict["dt"]
            self.run_t = self.oldSim_timing_dict["last_step"] * self.oldSim_timing_dict["dt"]
            self.total_t = self.run_t + self.extend_t
            self.total_steps = self.total_t / self.dt
            self.start_step = self.run_t / self.dt + 1
            self.remaining_steps = self.total_steps - self.start_step
            self.intervals = round(self.total_steps / self.oldSim_timing_dict["interval_step"])
            self.interval_arr = [int(i * self.oldSim_timing_dict["interval_step"]) for i in range(self.intervals)]
            self.output_intervals = (0.001, 0.005, 0.01, 0.1, 0.25, 0.5, 0.75, 1)
            self.output_arr = [self.start_step + (round(self.output_intervals[a] * self.remaining_steps, 0)) for a in
                               range(len(self.output_intervals))]
            self.output_arr = tuple(self.output_arr)

            with h5py.File(self.file_name, mode='a') as self.hdf:
                if self.already_extended == False:
                    timing = self.hdf.get("TIMING")
                elif self.already_extended ==True:
                    try:
                        self.hdf.create_group("TIMING_2")
                        timing = self.hdf.get("TIMING_2")
                    except:
                        self.hdf.create_group("TIMING_3")
                        timing = self.hdf.get("TIMING_3")
                timing.create_dataset("EXTENDED_DT", data=self.dt)
                timing.create_dataset("EXTENDED_TOTAL_T", data=self.total_t)
                timing.create_dataset("EXTENDED_INTERVALS", data=self.intervals)
                timing.create_dataset("EXTENDED_INTERVAL_STEP", data=self.oldSim_timing_dict["interval_step"])
            for i in range(len(self.comp_arr)):
                self.comp_arr[i].dt = self.dt

        if self.amend_type == "Resume":
            self.dt = self.oldSim_timing_dict["dt"]
            self.run_t = self.oldSim_timing_dict["last_step"] * self.oldSim_timing_dict["dt"]
            self.total_t = self.oldSim_timing_dict["total_t"]
            self.total_steps = self.total_t / self.dt
            self.start_step = self.run_t / self.dt + 1
            self.remaining_steps = self.total_steps - self.start_step
            self.intervals = round(self.total_steps / self.oldSim_timing_dict["interval_step"])
            self.interval_arr = [int(i * self.oldSim_timing_dict["interval_step"]) for i in range(self.intervals)]
            self.output_intervals = (0.001, 0.005, 0.01, 0.1, 0.25, 0.5, 0.75, 1)
            self.output_arr = [self.start_step + (round(self.output_intervals[a] * self.remaining_steps, 0)) for a in
                               range(len(self.output_intervals))]
            self.output_arr = tuple(self.output_arr)

            with h5py.File(self.file_name, mode='a') as self.hdf:
                timing = self.hdf.get("TIMING")
                try:
                    timing.create_dataset("RESUME_T", data=self.run_t)
                except:
                    timing.create_dataset("RESUME_T2", data=self.run_t)

            for i in range(len(self.comp_arr)):
                self.comp_arr[i].dt = self.dt

        if self.amend_type == "LastValues":
            self.dt = dt
            self.total_t = total_t
            self.intervals = intervals
            self.total_steps = self.total_t / self.dt
            self.interval_step = int(self.total_steps / self.intervals)
            self.run_t = 0
            self.start_step = 1
            self.interval_arr = [i * self.interval_step for i in range(self.intervals)]
            self.output_intervals = (0.001, 0.005, 0.01, 0.1, 0.25, 0.5, 0.75, 1)
            self.output_arr = [round(self.total_steps * _) for _ in self.output_intervals]
            self.output_arr = tuple(self.output_arr)

            with h5py.File(self.file_name, mode='a') as self.hdf:
                timing = self.hdf.create_group("TIMING")
                T = self.hdf.get("TIMING")
                T.create_dataset("DT", data=self.dt)
                T.create_dataset("TOTAL_T", data=total_t)
                T.create_dataset("INTERVALS", data=intervals)

            for i in range(len(self.comp_arr)):
                self.comp_arr[i].dt = self.dt

    def set_external_ion_properties(self, na_o=145e-3, k_o=3.5e-3, cl_o=119e-3, x_o=29.5e-3, z_o=-0.85):
        """
        Capacity to change the extracellular bath properties before the simulation
        @param na_o: float,extracellular sodium concentration in Molar
        @param k_o: float,extracellular potassium concentration in Molar
        @param cl_o: float,extracellular chloride concentration in Molar
        @param x_o: float,extracellular impermeant anion concentration
        @param z_o : float,extracellular impermeant anion average charge
        """
        self.na_o, self.k_o, self.cl_o, self.x_o, self.z_o = na_o, k_o, cl_o, x_o, z_o
        self.x_o = -1 * (self.cl_o - self.na_o - self.k_o)
        self.osm_o = self.x_o + self.na_o + self.cl_o + self.k_o

    def set_sa(self,rad=0.5e-5, length=20e-5):
        """
        Function to preset surface area
        @param rad: radius in dm
        @param len: length in dm
        @return:
        """

        for i in range(len(self.comp_arr)):
            if self.comp_arr[i].name =='Comp0(Soma)':
                self.comp_arr[i].sa = 2 * np.pi * (2 * rad) * (2 * length)
            else:
                self.comp_arr[i].sa  = 2 * np.pi * (rad) * (length)
        return

    def set_static_sa(self, static_sa=True):
        """
        Function to fix surface area
        @param static_sa: boolean which sets the surface area of the compartment to static
            A static surface area means the cell membrane wrinkles or becomes more turgid as the volume changes.
            The implication of this is that the channel conductances aren't scaled by surface area.
        """
        self.static_sa = static_sa
        sa_values_data =[]
        for i in range(len(self.comp_arr)):
            self.comp_arr[i].static_sa = self.static_sa
            sa_values_data.append(self.comp_arr[i].sa)

        with h5py.File(self.file_name, mode='a') as self.hdf:
            static_settings = self.hdf.create_group("SA_SETTINGS")
            static_settings.create_dataset(name="SA-VALUES", data=sa_values_data)

    def set_na_i_start(self,na_i=14e-3):
        """
        Preset starting sodium concentrations
        @param na_i: starting sodium concentrations
        @return:
        """
        for _ in self.comp_arr:
            _.na_i_start = na_i
        return

    def set_atpase_static(self, static_atpase=True):
        """
        Function to change ATPase pump rate from dynamic (dependant on sodium concentration) to static (independent)
        @param static_atpase: boolean --> if true ATPase is static
        """
        self.static_atpase = static_atpase
        for i in range(len(self.comp_arr)):
            self.comp_arr[i].static_sa = self.static_atpase

        atpase_settings_data= []
        if static_atpase:
            atpase_settings_data.append(1)
        else:
            atpase_settings_data.append(0)

        atpase_settings_data.append(self.comp_arr[1].p_atpase)
        atpase_settings_data.append(self.comp_arr[1].na_i_start)

        with h5py.File(self.file_name, mode='a') as self.hdf:
            self.hdf.create_group("ATPASE_SETTINGS")
            atpase_settings = self.hdf.get("ATPASE_SETTINGS")
            atpase_settings.create_dataset(name="STATIC", data=atpase_settings_data[0])
            atpase_settings.create_dataset(name="P-ATPASE", data=atpase_settings_data[1])
            atpase_settings.create_dataset(name="NA-START", data=atpase_settings_data[2])

    def set_gkcc2(self, gkcc2=2e-3):
        self.gkcc2 = gkcc2/F
        for i in self.comp_arr:
            i.p_kcc2 = self.gkcc2

    def get_avg_osmo(self, excl_comp_name=''):
        """
        Determines the average osmolarity in all the compartments.
        @param excl_comp_name: exclude a particular compartments from the calculation
        - useful for calculating the average osmolarity when z is being fluxed in a particular compartment.
        """
        total_osmo = 0
        counter = 0
        for i in self.comp_arr:
            if i.name != excl_comp_name and i.name != "Comp0(Soma)":
                total_osmo += i.osm_i
                counter += 1
        avg_osmo = total_osmo / counter
        return avg_osmo

    def add_current(self, comp_name='', current_type='Inhibitory', start_t=0, duration=2 * 1e-3,
                    current_A=1e-3, dt=1e-6, total_t=1, frequency=1):
        """
        Function to add a current to a particular compartment
        @param comp_name: string, compartment name that current is being pulsed into
        @param current_type: string, "Excitatory" or "Inhibitory" current
        @param start_t: float, start time of the current (in seconds)
        @param duration: float, duration of current (in seconds)
        @param current_A: float, amount of current injected (in Amperes)
        @param total_t: float,total simulation time
        @param frequency: float, Hz (pulses/second)
        @return:
        """

        self.sim_current_params["On"] = True

        comp_num = 0
        for i in range(len(self.comp_arr)):
            if comp_name == self.comp_arr[i].name:
                comp_num = i
                self.sim_current_params["compartment"] = comp_num

        if current_type == "Inhibitory":
            self.current_type = 0
            self.sim_current_params["current_type"] = "inhibitory"
        elif current_type == "Excitatory":
            self.current_type = 1
            self.sim_current_params["current_type"] = "excitatory"

        self.sim_current_params["start_t"] = start_t + self.run_t
        self.sim_current_params["duration"] = duration
        self.sim_current_params["current_A"] = current_A
        self.sim_current_params["frequency"] = frequency

        self.pulse_start_t_arr = []
        self.pulse_end_t_arr = []
        remaining_time = total_t - start_t
        time_intervals = remaining_time / frequency

        for i in range(frequency):
            pulse_start_t = time_intervals * (i) + start_t + self.run_t
            pulse_end_t = pulse_start_t + duration
            self.pulse_start_t_arr.append(pulse_start_t + self.run_t)
            self.pulse_end_t_arr.append(pulse_end_t)

        self.comp_arr[comp_num].dt = dt
        self.comp_arr[comp_num].set_current(current_type=self.current_type, current_amplitude=current_A)
        current_settings_list = list(self.sim_current_params.values())
        current_settings_list = current_settings_list[1:]
        if current_settings_list[1] == 'excitatory':
            current_settings_list[1] = 1
        else:
            current_settings_list[1] =0
        with h5py.File(self.file_name, mode='a') as self.hdf:
            current = self.hdf.create_group("CURRENT")
            C = self.hdf.get("CURRENT")
            C.create_dataset("CURR", data=current_settings_list)

    def add_synapse(self, comp_name='', synapse_type='Inhibitory', start_t=0, duration=2 * 1e-3,
                    max_neurotransmitter=1e-3, synapse_conductance=1e-9):
        """
        function to add a synapse to a particular compartment
        @param synapse_type: string,either 'Inhibitory' (GABAergic) or 'Excitatory' (Glutamatergic)
        @param comp_name: string,compartment name on which to synapse onto
        @param start_t: float,start time for synaptic input
        @param duration: float,duration of synaptic input
        @param max_neurotransmitter: float,max neurotransmitter concentration in moles/liter
        @param synapse_conductance: float, conductance of the synapse channel in Siemens, default is 1nS
        """
        syn_dict = {}
        self.syn_count +=1
        syn_dict["On"] = True
        syn_dict["Synapse_num" ] = self.syn_count
        comp_num = 0
        for i in range(len(self.comp_arr)):
            if comp_name == self.comp_arr[i].name:
                comp_num = i
                syn_dict["compartment"] = comp_num

        if synapse_type == "Excitatory":
            syn_dict["synapse_type"] = 0
            syn_dict["alpha"] = 2e6  # ms-1.mM-1 --> s-1.M-1= Forward rate constant
            syn_dict["beta"] = 1e3  # ms-1 --> s-1 == Backward rate constant
        elif synapse_type == "Inhibitory":
            syn_dict["synapse_type"] = 1
            syn_dict["alpha"] = 0.5e6  # ms-1.mM-1 --> s-1.M-1= Forward rate constant
            syn_dict["beta"] = 0.1e3  # ms-1 --> s-1 == Backward rate constant

        syn_dict["start_t"] = start_t
        syn_dict["duration"] = duration
        syn_dict["end_t"] = start_t + duration
        syn_dict["max_neurotransmitter_conc"] = max_neurotransmitter
        syn_dict["synapse_conductance"] = synapse_conductance
        self.comp_arr[comp_num].set_synapse(syn_dict, self.dt)
        self.syn_arr.append(syn_dict)

        with h5py.File(self.file_name, mode='a') as self.hdf:
            synapse_name = "SYNAPSE-" + self.comp_arr[comp_num].name
            syn_data_arr = list(syn_dict.values())
            try:
                self.hdf.create_group("SYNAPSE-SETTINGS")
            except:
                pass
            synapse_group = self.hdf.get("SYNAPSE-SETTINGS")
            synapse_group.create_dataset(name=synapse_name, data=syn_data_arr)

        return

    def set_hh_on(self, comp="Comp0(Soma)", t_on=5):
        """
        Function to add Hodgkin-Huxley Channels to a compartment
        @param comp: string, name of the compartment which requires HH channels
            default is the name of the soma
        @param t_on: float, time when to add HH channels
            needs to be earlier than total simulation run time
        """
        self.hh_on = True
        self.hh_t_on = t_on
        for i in range(len(self.comp_arr)):
            if self.comp_arr[i].name == comp:
                self.hh_comp_num = i

    def run_simulation(self):
        self.start_t = time.time()
        if self.amend_type == "Resume" or self.amend_type =="Extend":
            self.interval_num = self.oldSim_timing_dict["last_interval"] + 1
        else:
            self.interval_num = 1
        self.steps = int(self.start_step)

        while self.run_t < self.total_t:

            for a in self.comp_arr:

                a.step()  # step for each compartment

                if a.xflux_switch and (a.xflux_params["start_t"] <= self.run_t < a.xflux_params["end_t"]):
                    a.x_flux()

                if a.zflux_switch and \
                        (a.zflux_params["start_t"] <= self.run_t < a.zflux_params[
                            "end_t"]):
                    a.z_flux()

                if a.synapse_on:
                    a.synapse_step(run_t=self.run_t)

                if a.current_on:
                    for i in range(len(self.pulse_start_t_arr)):
                        if (self.run_t >= self.pulse_start_t_arr[i]) and (self.run_t <= self.pulse_end_t_arr[i]):
                            a.current_step(dt=self.dt)

                    # electrodiffusion dictionary for each compartment

            if self.ED_on:
                for b in range(len(self.ed_arr)):
                    ed_conc_changes = self.ed_arr[b].calc_ed(self.dt, self.comp_arr[b].w,
                                                             self.comp_arr[b].get_ed_dict(),
                                                             self.comp_arr[b + 1].get_ed_dict())
                    self.comp_arr[b].ed_update(ed_conc_changes, "positive")
                    self.comp_arr[b + 1].ed_update(ed_conc_changes, "negative")

                    # appending the electrodiffusion concentrations for each compartment

            for d in self.comp_arr:
                if d.adjust_x_bool:
                    avg_osm = self.get_avg_osmo(d.name)
                    d.update_volumes(avg_osm)
                else:
                    d.update_volumes()
                    # updates of the volumes, arrays for each compartment

            for f in range(len(self.output_arr) - 1):
                if self.steps == self.output_arr[f]:
                    if f == 2:
                        self.one_percent_t = time.time() - self.start_t
                        self.hundred_percent_t = self.one_percent_t * 100
                        print(str(self.output_intervals[f] * 100) + " % complete in " + str(
                            round(self.one_percent_t, 2)) + " s")
                        print("Estimated time to complete :" + str(
                            round(self.hundred_percent_t / 60, 2)) + " minutes")
                    else:
                        print(str(self.output_intervals[f] * 100) + " % complete in " + str(
                            round(time.time() - self.start_t, 2)) + " s")

            if self.hh_on and self.run_t >= self.hh_t_on:
                self.comp_arr[self.hh_comp_num].hh_on = True

            if self.interval_num < len(self.interval_arr):
                if self.steps == self.interval_arr[self.interval_num]:
                    self.interval_num += 1
                    self.save_to_file()

            self.steps += 1

            self.run_t += self.dt

        print("100.0 % complete in " + str(
            round(time.time() - self.start_t, 2)) + " s")


    def save_to_file(self):
        """
        Function to save simulation to HDF file.
        This is an internal function and should not be called by the user
        """
        with h5py.File(self.file_name, mode='a') as self.hdf:
            for i in range(len(self.comp_arr)):
                group = self.hdf.get('COMPARTMENTS')
                subgroup = group.get(self.comp_arr[i].name)
                data_array = self.comp_arr[i].get_array(self.run_t)
                subgroup.create_dataset(name=str(self.steps), data=data_array)

            if self.ED_on:
                for j in range(len(self.ed_arr)):
                    group = self.hdf.get('ELECTRODIFFUSION')
                    subgroup = group.get(self.ed_arr[j].name)
                    data_array = self.ed_arr[j].ed_change_arr
                    subgroup.create_dataset(name=str(self.steps), data=data_array)


    def write_settings_to_file(self, basefile_name=""):
        settings_file_name = self.file_name + "_settings"
        f = open(settings_file_name, "w")
        f.write("SETTINGS")
        f.write("\n")
        f.write("=============================")
        f.write("\n")
        if basefile_name != "":
            f.write("HDF5 Base File (Template) " + basefile_name)
            f.write("\n")
        else:
            f.write("No Base/Template)")
            f.write("\n")
        f.write("HDF5 File name:" + self.file_name)
        f.write("\n")
        f.write("\n")
        f.write("==============================")
        f.write("\n")
        f.write("Dendrite Morphology")
        f.write("\n")
        f.write("==============================")
        f.write("\n")
        for i in self.comp_arr:
            f.write(i.name + " Len(um): " + str(i.length * 1e5) + " Rad(um): " + str(i.radius * 1e5))
            f.write("\n")
        f.write("\n")
        f.write("==============================")
        f.write("\n")
        f.write("General settings")
        f.write("\n")
        f.write("==============================")
        f.write("\n")
        f.write("ED on: " + str(self.ED_on))
        f.write("\n")
        f.write("Surface Area Constant: " + str(self.static_sa))
        f.write("\n")
        f.write("Total sim time(s):" + str(self.total_t))
        f.write("\n")
        f.write("Time step(s): " + str(self.dt))
        f.write("\n")
        f.write("#Intervals recorded in file: " + str(self.intervals))
        f.write("\n")
        f.write("Starting [Na]: " + str(self.comp_arr[1].na_i * 1e3))
        f.write("\n")
        f.write("Starting [K]: " + str(self.comp_arr[1].k_i * 1e3))
        f.write("\n")
        f.write("Starting [Cl]: " + str(self.comp_arr[1].cl_i * 1e3))
        f.write("\n")
        f.write("Starting [X]: " + str(self.comp_arr[1].x_i * 1e3))
        f.write("\n")
        f.write("\n")
        f.write("\n")
        f.write("==============================")
        f.write("\n")
        f.write("Pump settings")
        f.write("\n")
        f.write("==============================")
        f.write("\n")
        f.write("ATPase rate Constant: " + str(self.static_atpase))
        f.write("\n")
        f.write("ATPase starting [Na]: " + str(self.comp_arr[1].na_i_start * 1e3))
        f.write("\n")
        f.write("ATPase pump constant: " + str(self.comp_arr[1].p_atpase))
        f.write("\n")
        f.write("g_KCC2: " + str(self.comp_arr[1].p_kcc2))
        f.write("\n")
        f.write("\n")
        f.write("==============================")
        f.write("\n")
        f.write("[X]-flux settings")
        f.write("\n")
        f.write("==============================")
        f.write("\n")
        if self.xflux_count == 0:
            f.write("No [X]-flux")
            f.write("\n")
            f.write("\n")
        else:
            for i in range(len(self.xflux_names_arr)):
                f.write(self.xflux_names_arr[i])
                f.write("\n")
                f.write("Flux rate: " + str(self.xflux_data_arr[0]))
                f.write("\n")
                f.write("z: " + str(self.xflux_data_arr[1]))
                f.write("\n")
                f.write("Start time: " + str(self.xflux_data_arr[2]))
                f.write("\n")
                f.write("End time: " + str(self.xflux_data_arr[3]))
                f.write("\n")
                f.write("\n")
        f.write("==============================")
        f.write("\n")
        f.write("z-flux settings")
        f.write("\n")
        f.write("==============================")
        f.write("\n")
        if not self.sim_zflux_params["On"]:
            f.write("No z-flux")
            f.write("\n")
            f.write("\n")
        else:
            f.write(self.sim_zflux_params["Comp"])
            f.write("\n")
            f.write("z end: " + str(self.sim_zflux_params["z_end"]))
            f.write("\n")
            f.write("Start time: " + str(self.sim_zflux_params["start_t"]))
            f.write("\n")
            f.write("End time: " + str(self.sim_zflux_params["end_t"]))
            f.write("\n")
            f.write("Adjust x: " + str(self.sim_zflux_params["adjust_X"]))
            f.write("\n")
            f.write("\n")
        f.write("==============================")
        f.write("\n")
        f.write("Current settings")
        f.write("\n")
        f.write("==============================")
        f.write("\n")
        if not self.sim_current_params["On"]:
            f.write("No external current added")
            f.write("\n")

        else:
            f.write("Current onto: Comp" + str(self.sim_current_params["compartment"]))
            f.write("\n")
            f.write("Current type: " + self.sim_current_params["current_type"])
            f.write("\n")
            f.write("Start time: " + str(self.sim_current_params["start_t"]))
            f.write("\n")
            f.write("Duration: " + str(self.sim_current_params["duration"]))
            f.write("\n")
            f.write("Current Amplitude (A): " + str(self.sim_current_params["current_A"]))
            f.write("\n")
            f.write("Current Frequency (Hz): " + str(self.sim_current_params["frequency"]))
        f.write("\n")
        f.write("\n")
        f.write("==============================")
        f.write("\n")
        f.write("Synapse settings")
        f.write("\n")
        f.write("==============================")
        f.write("\n")
        if len(self.syn_arr) == 0:
            f.write("\n")
            f.write("No synapses added")
            f.write("\n")
        else:
            for j in range(len(self.syn_arr)):
                f.write("\n")
                f.write("Synapse number:" + str(self.syn_arr[j]["Synapse_num"]))
                f.write("\n")
                f.write("Synapse on:" + self.comp_arr[self.syn_arr[j]["compartment"]].name)
                f.write("\n")
                if self.syn_arr[j]["synapse_type"] == 0:
                    f.write("\n")
                    f.write("Synapse type: Excitatory")
                else:
                    f.write("\n")
                    f.write("Synapse type: Inhibitory")
                f.write("\n")
                f.write("alpha rate constant: " + str(self.syn_arr[j]["alpha"]))
                f.write("\n")
                f.write("beta rate constant: " + str(self.syn_arr[j]["beta"]))
                f.write("\n")
                f.write("start_t: " + str(self.syn_arr[j]["start_t"]))
                f.write("\n")
                f.write("end_t: " + str(self.syn_arr[j]["end_t"]))
                f.write("\n")
                f.write("max [NT](M) :" + str(self.syn_arr[j]["max_neurotransmitter_conc"]))
                f.write("\n")
                f.write("synaptic conductance: " + str(self.syn_arr[j]["synapse_conductance"]))
                f.write("\n")
        f.write("\n")
        f.write("==============================")
        f.write("\n")
        f.write("HH settings")
        f.write("\n")
        f.write("==============================")
        if not self.hh_on:
            f.write("\n")
            f.write("No HH channels ")
            f.write("\n")
        else:
            f.write("\n")
            f.write("HH channels enabled")
            f.write("\n")
        f.write("\n")
        f.write("===============================")

        f.close()
        print(settings_file_name + " created")
