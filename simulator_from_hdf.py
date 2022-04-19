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
import time

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
            oldSim_last_interval = len(comp_intervals_list) - 2
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
                self.na_start = dataset_start[4]
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
            for e in range(self.num_comps - 1):
                name = self.comp_arr[e].name + ' -> ' + self.comp_arr[e + 1].name

                ed = electrodiffusion.Electrodiffusion(comp_a_name=self.comp_arr[e].name,
                                                       comp_a_length=self.comp_arr[e].length,
                                                       comp_b_name=self.comp_arr[e + 1].name,
                                                       comp_b_length=self.comp_arr[e + 1].length)

                ed.set_diff_constants(self.diff_constants)
                self.ed_arr.append(ed)

    def set_timing(self, extend_t=30):

        if self.amend_type == "Extend":
            self.extend_t = extend_t
            self.dt = self.oldSim_timing_dict["dt"]
            self.run_t = self.oldSim_timing_dict["last_step"] * self.oldSim_timing_dict["dt"]
            self.total_t = self.run_t + self.extend_t
            self.total_steps = self.total_t / self.dt
            self.start_step = self.run_t / self.dt
            self.remaining_steps = self.total_steps - self.start_step
            self.intervals = round(self.total_steps / self.oldSim_timing_dict["interval_step"])
            self.interval_arr = [i * self.oldSim_timing_dict["interval_step"] for i in range(self.intervals)]
            self.output_intervals = (0.001, 0.005, 0.01, 0.1, 0.25, 0.5, 0.75, 1)
            self.output_arr = [self.start_step+ (round(self.output_intervals[a] * self.remaining_steps, 0)) for a in
                               range(len(self.output_intervals))]
            self.output_arr = tuple(self.output_arr)

            with h5py.File(self.file_name, mode='a') as self.hdf:
                timing = self.hdf.get("TIMING")
                timing.create_dataset("EXTENDED_DT", data=self.dt)
                timing.create_dataset("EXTENDED_TOTAL_T", data=self.total_t)
                timing.create_dataset("EXTENDED_INTERVALS", data=self.intervals)

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

    def set_static_sa(self, static_sa=True):
        """
        Function to fix surface area
        @param static_sa: boolean which sets the surface area of the compartment to static
            A static surface area means the cell membrane wrinkles or becomes more turgid as the volume changes.
            The implication of this is that the channel conductances aren't scaled by surface area.
        """
        self.static_sa = static_sa
        for i in range(len(self.comp_arr)):
            self.comp_arr[i].static_sa = self.static_sa

    def set_atpase_static(self, static_atpase=True):
        """
        Function to change ATPase pump rate from dynamic (dependant on sodium concentration) to static (independent)
        @param static_atpase: boolean --> if true ATPase is static
        """
        self.static_atpase = static_atpase
        for i in range(len(self.comp_arr)):
            self.comp_arr[i].static_sa = self.static_atpase
            self.comp_arr[i].na_i_start = self.na_start

    def run_simulation(self):

        self.start_t = time.time()
        self.interval_num = self.oldSim_timing_dict["last_interval"]
        self.steps = self.start_step

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
        f.write("ATPase rate Constant: " + str(self.static_atpase))
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
            f.write("Current onto: " + self.comp_arr[self.sim_current_params["compartment"].name])
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
        if not self.syn_dict["On"]:
            f.write("\n")
            f.write("No synapses added")
            f.write("\n")
        else:
            f.write("\n")
            f.write("Synapse on:" + self.comp_arr[self.syn_dict["compartment"]].name)
            f.write("\n")
            if self.syn_dict["synapse_type"] == 0:
                f.write("\n")
                f.write("Synapse type: Excitatory")
            else:
                f.write("\n")
                f.write("Synapse type: Inhibitory")
            f.write("\n")
            f.write("alpha rate constant: " + str(self.syn_dict["alpha"]))
            f.write("\n")
            f.write("beta rate constant: " + str(self.syn_dict["beta"]))
            f.write("\n")
            f.write("start_t: " + str(self.syn_dict["start_t"]))
            f.write("\n")
            f.write("end_t: " + str(self.syn_dict["end_t"]))
            f.write("\n")
            f.write("max [NT](M) :" + str(self.syn_dict["max_neurotransmitter_conc"]))
            f.write("\n")
            f.write("synaptic conductance: " + str(self.syn_dict["synapse_conductance"]))
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
