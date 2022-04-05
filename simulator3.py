"""

Control the functioning of the Neural Physiological Emulator
Main simulation class
-parent classes: directed with controller class and controller_from_HDF classes
-child clases: compartment and electrodiffusion classes


Creator: Eran Shorer
"""
import time
import h5py
import numpy as np
import pandas as pd
import compartment
import electrodiffusion


class Simulator:

    def __init__(self, file_name=""):
        """
        Initializing a new instance of the simulator
        @param file_name: string of the file name
        """

        self.file_name = file_name
        self.file_name = "\ " + file_name

        try:
            with h5py.File(self.file_name, mode='w') as self.hdf:
                self.hdf.create_group('COMPARTMENTS')
                self.hdf.create_group('ELECTRODIFFUSION')
                self.hdf.create_group("TIMING")
                self.hdf.create_group("X-FLUX-SETTINGS")
                self.hdf.create_group("CURRENT-SETTINGS")
                self.hdf.create_group("SYNAPSE-SETTINGS")
                print("simulation file ('" + file_name + "') created in base directory")
        except:
            raise Exception("File not created")

        self.num_comps = 0
        self.one_percent_t = 0
        self.interval_num = 1
        self.intervals = 0
        self.steps = 0
        self.ED_on = True
        self.diff_constants = {}
        self.static_atpase, self.static_sa = True, True
        self.comp_arr, self.ed_arr = [], []
        self.external_xflux_setup, self.xflux_setup, self.zflux_setup = True, True, True
        self.na_o, self.k_o, self.cl_o, self.x_o, self.z_o, self.osm_o = 0, 0, 0, 0, 0, 0
        self.p = 0
        self.start_t, self.end_t, self.run_t, self.total_t, self.dt = 0, 0, 0, 0, 0
        self.xflux_dict = {}
        self.xflux_count = 0
        self.xoflux_switch = False
        self.xoflux_params = {"start_t": 0, "end_t": 0, "xo_conc": 0, "zo": 0}
        self.xoflux_setup = True
        self.xo_start, self.cl_o_start, self.d_xoflux, self.xo_final, self.xo_flux, self.t_xoflux = 0, 0, 0, 0, 0, 0
        self.xoflux_points, self.dt_xoflux, self.xo_alpha, self.xo_beta = 0, 0, 0, 0
        self.sim_zflux_params = {"On": False}
        self.sim_current_params = {"On": False}
        self.syn_dict = {"On": False}
        self.hh_on = False
        self.hh_t_on = 0
        self.hh_comp_num = 0

    def add_compartment(self, comp=compartment):
        """
        Adding compartment object to the simulator
        @param comp: compartment object
        """

        new_comp = comp.get_array(time=0)
        with h5py.File(self.file_name, mode='a') as self.hdf:
            group = self.hdf.get('COMPARTMENTS')
            subgroup = group.create_group(name=comp.name)
            subgroup.create_dataset(name='0', data=new_comp)

        self.num_comps += 1
        self.comp_arr.append(comp)

    def add_default_multicompartment(self, number_of_comps=9, rad=0.5e-5, len=20e-5, soma=False):
        """
        Sets the simulation to run with the default multicompartment model
         -- 9 compartments + 1 soma
        @param number_of_comps: integer,#of compartments arranged linearly in a multicompartment setup
        @param rad: float, radius of a single compartment (in decimeters); default values are micrometers converted to decimeters
        @param len: floot, length of a single compartment (in decimeters)
        @param soma: boolean, value to add a soma to the end of the multicompartment chain. Soma length and radius is
        double that of a standard compartment. Default is to have hodgkin huxley channels on the soma off.
        """

        if soma:
            soma = compartment.Compartment("Comp0(Soma)", radius=2 * rad, length=2 * len, static_sa=True,
                                           sa_value=2 * np.pi * (2 * rad) * (2 * len), hh_on=False)
            # soma.set_ion_properties(na_i=0.013995241563512785, k_i=0.12286753014443351, cl_i=0.005171468255812758, x_i=0.15496634531836323)
            soma.set_ion_properties()
            self.add_compartment(soma)

        for i in range(number_of_comps):
            comp = compartment.Compartment("Comp" + str(i + 1), radius=rad, length=len, static_sa=True,
                                           sa_value=2 * np.pi * rad * len, hh_on=False)
            comp.set_ion_properties()
            self.add_compartment(comp)


    def set_electrodiffusion_properties(self, ED_on: bool = True,
                                        diff_constant_dict=None):
        """
        Function to define the electrodiffusion setup at the beginning of a simulation.
        @param ED_on: boolean, sets either electrodiffusion between compartments to be on or off.
        @param diff_constant_dict: dictionary, diffusion constants in the format {"ion":value,"ion":value,...}

        """
        if diff_constant_dict is None:
            diff_constant_dict = {"na": 1.33e-7, "k": 1.96e-7, "cl": 2.03e-7,
                                  "x": 0}
        self.ED_on = ED_on
        self.diff_constants = diff_constant_dict

        with h5py.File(self.file_name, mode='a') as self.hdf:

            ed_group = self.hdf.get('ELECTRODIFFUSION')

            if self.ED_on:
                for e in range(self.num_comps - 1):
                    name = self.comp_arr[e].name + ' -> ' + self.comp_arr[e + 1].name
                    ed_group.create_group(name)

                    ed = electrodiffusion.Electrodiffusion(comp_a_name=self.comp_arr[e].name,
                                                           comp_a_length=self.comp_arr[e].length,
                                                           comp_b_name=self.comp_arr[e + 1].name,
                                                           comp_b_length=self.comp_arr[e + 1].length)

                    ed.set_diff_constants(self.diff_constants)
                    self.ed_arr.append(ed)

            else:
                ed_group.create_group("NO ED")

            self.hdf.close()

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

    def set_z(self, comp="Comp8", z=-0.85, adjust_x=False):
        """
        Function which sets the impermeant anion average charge in a compartment

        @param comp: string, compartment name. "All" if all compartments are desired.
        @param z :, average charge of impermeant anion desired for compartment
        @param adjust_x: boolean, allows for the compartment to adjust IA concentration based on z value
                        to ensure volume remains constant
        """
        if comp == "All":
            for i in range(len(self.comp_arr)):
                self.comp_arr[i].z_i = z
                self.comp_arr[i].set_osmo_neutral_start()
                if adjust_x:
                    self.comp_arr[i].bool_adjust_x = True
        else:
            for i in range(len(self.comp_arr)):
                if comp == self.comp_arr[i].name:
                    self.comp_arr[i].z_i = z
                    self.comp_arr[i].set_osmo_neutral_start()
                    if adjust_x:
                        self.comp_arr[i].bool_adjust_x = True

    def set_atpase_static(self, static_atpase=True):
        """
        Function to change ATPase pump rate from dynamic (dependant on sodium concentration) to static (independent)
        @param static_atpase: boolean --> if true ATPase is static
        """
        self.static_atpase = static_atpase
        for i in range(len(self.comp_arr)):
            self.comp_arr[i].static_sa = self.static_atpase

    def set_sa_static(self, static_sa=True):
        """
        Function to fix surface area
        @param static_sa: boolean which sets the surface area of the compartment to static
            A static surface area means the cell membrane wrinkles or becomes more turgid as the volume changes.
            The implication of this is that the channel conductances aren't scaled by surface area.
        """
        self.static_sa = static_sa
        for i in range(len(self.comp_arr)):
            self.comp_arr[i].static_sa = self.static_sa

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

    def set_timing(self, total_t, time_step=1e-6, intervals=1000):
        """
        @param total_t: float, total simulation time (in seconds)
        @param time_step: float, time step of the simulation(in seconds) this determines how often calculations are repeated
            default is set to 1e-6 seconds.
        @param intervals: integer, determines how often to write simulation results to the HDF file
        """
        self.total_t, self.dt = total_t, time_step
        with h5py.File(self.file_name, mode='a') as self.hdf:
            timing = self.hdf.get("TIMING")
            timing.create_dataset("DT", data=time_step)
            timing.create_dataset("TOTAL_T", data=total_t)
            timing.create_dataset("INTERVALS", data=intervals)

        self.total_steps = self.total_t / self.dt
        self.intervals = intervals
        self.output_intervals = (0.001, 0.005, 0.01, 0.1, 0.25, 0.5, 0.75, 1)
        self.output_arr = [round(self.output_intervals[a] * self.total_steps, 0) for a in
                           range(len(self.output_intervals))]
        self.output_arr = tuple(self.output_arr)

        self.interval_step = self.total_steps / intervals
        self.interval_arr = [round(self.interval_step * i) for i in range(intervals)]
        self.interval_arr = tuple(self.interval_arr)

        for i in range(len(self.comp_arr)):
            self.comp_arr[i].dt = self.dt

    # print(self.interval_arr)

    def set_xflux(self, comps=None, start_t=0, end_t=0,
                  flux_rate=1e-5, z=-0.85):
        """
        function which defines how impermeant anion flux should be constructed
        @param comps: string, receives the name of a compartment object where intracellular impermeant anions
            need to be fluxed
        @param start_t: float, time in seconds to start the flux
        @param end_t: float, time in seconds when flux should end
        @param flux_rate: float, MOLES (not concentration) per second flux should occur at
        @param z: float, average charge of added impermeant anions
        """
        self.xflux_data_arr, self.xflux_names_arr = [], []  # array which will be sent to the HDF5 file

        for i in range(len(self.comp_arr)):
            if comps[0] == self.comp_arr[i].name:
                self.xflux_data_arr.append(i)

        self.xflux_data_arr.extend((flux_rate, z, start_t, end_t))

        for i in range(len(comps)):
            for j in self.comp_arr:
                if comps[i] == j.name:
                    j.xflux_switch = True
                    j.xflux_params["flux_rate"] = flux_rate
                    j.xflux_params["z"] = z
                    j.xflux_params["start_t"] = start_t
                    j.xflux_params["end_t"] = end_t

        self.xflux_names_arr.append("X-FLUX-" + str(self.xflux_count))  # names of the xflux
        self.xflux_count += 1
        with h5py.File(self.file_name, mode='a') as self.hdf:
            xflux_group = self.hdf.get("X-FLUX-SETTINGS")
            xflux_group.create_dataset(name=self.xflux_names_arr[-1], data=self.xflux_data_arr)

    def set_zflux(self, all_comps=False, comps=None, start_t=0, end_t=0, z_end=-1,adjust_x=False):
        """
        function which changes the average charge of impermeant anions in a compartment
        @param all_comps: boolean, if true impermeant anion change occurs in all compartments
        @param comps: string, receives the name of the compartment to change the average impermeant charge
            if all compartments is true this can be left out.
        @param start_t: float, time to start charge flux in seconds
        @param end_t: float, time to end charge flux in seconds
        @param adjust_x: boolean, if true changes the impermeant anion concentration of the compartment
            such that there is no volume change.
       """


        for i in range(len(comps)):
            for j in range(len(self.comp_arr)):
                if comps[i] == self.comp_arr[j].name or all_comps:
                    self.comp_arr[j].zflux_switch = True
                    self.comp_arr[j].zflux_params = {"start_t": start_t, "end_t": end_t, "z": z_end,}
                    self.comp_arr[j].adjust_x_bool = adjust_x
                    self.sim_zflux_params = {"On": True, "Comp": self.comp_arr[j].name, "z_end": z_end,
                                             "start_t": start_t, "end_t": end_t,
                                             "adjust_X": adjust_x}

    def set_xoflux(self, start_t=0, end_t=50, xo_conc=1e-3, z=-0.85):
        """
        CHANGE THE FLUX OF IMPERMEANTS OUTSIDE THE COMPARTMENT
        @Param start_t: float,time to flux external impermeant anions
        @param end_t: float,time to end flux of impermeants
        @param xo_conc: float,amount of impermeant anions to flux (in Mols/L)
        @param z: float,charge of impermeants
        """
        self.xoflux_switch = True
        self.xoflux_params = {"start_t": start_t, "end_t": end_t, "xo_conc": xo_conc, "zo": z}

    def xoflux(self):
        """
        implementation of impermeant anion change.
        Should not be called by public user. Internal function for the simulator
        """
        if self.xoflux_setup:
            # starting values for external x flux
            self.xo_start = self.x_o
            self.cl_o_start = self.cl_o
            self.xo_final = self.x_o + self.xoflux_params["xo_conc"]
            self.xoflux_points = (self.xoflux_params["end_t"] - self.xoflux_params["start_t"]) * (1 / self.dt)
            self.dt_xoflux = 4 / self.xoflux_points
            self.xo_alpha = 1
            self.xo_beta = -1
            self.xoflux_setup = False

        elif self.x_o <= self.xo_final:

            self.d_xoflux = self.xo_alpha - np.e ** (self.xo_beta * self.t_xoflux)
            self.xo_flux = self.d_xoflux * self.xoflux_params["xo_conc"]
            self.x_o = self.xo_start + self.xo_flux
            self.cl_o = self.cl_o_start - self.xo_flux  # balancing the charges added externally
            self.t_xoflux += self.dt_xoflux

        else:
            return

    def get_avg_osmo(self, excl_comp_name=''):
        """
        Determines the average osmolarity in all the compartments.
        @param excl_comp_name: exclude a particular compartments from the calculation
        - useful for calculating the average osmolarity when z is being fluxed in a particular compartment.

        """
        total_osmo = 0
        counter = 0
        for i in self.comp_arr:
            if i.name != excl_comp_name:
                total_osmo += i.osm_i
                counter += 1
        avg_osmo = total_osmo / counter
        return avg_osmo

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
        self.syn_dict["On"] = True
        comp_num = 0
        for i in range(len(self.comp_arr)):
            if comp_name == self.comp_arr[i].name:
                comp_num = i
                self.syn_dict["compartment"] = comp_num

        if synapse_type == "Excitatory":
            self.syn_dict["synapse_type"] = 0
            self.syn_dict["alpha"] = 2e6  # ms-1.mM-1 --> s-1.M-1= Forward rate constant
            self.syn_dict["beta"] = 1e3  # ms-1 --> s-1 == Backward rate constant
        elif synapse_type == "Inhibitory":
            self.syn_dict["synapse_type"] = 1
            self.syn_dict["alpha"] = 0.5e6  # ms-1.mM-1 --> s-1.M-1= Forward rate constant
            self.syn_dict["beta"] = 0.1e3  # ms-1 --> s-1 == Backward rate constant

        self.syn_dict["start_t"] = start_t
        self.syn_dict["duration"] = duration
        self.syn_dict["end_t"] = start_t + duration
        self.syn_dict["max_neurotransmitter_conc"] = max_neurotransmitter
        self.syn_dict["synapse_conductance"] = synapse_conductance
        self.comp_arr[comp_num].set_synapse(self.syn_dict, self.dt)

        with h5py.File(self.file_name, mode='a') as self.hdf:
            synapse_name = "SYNAPSE-" + self.comp_arr[comp_num].name
            syn_data_arr = list(self.syn_dict.values())
            synapse_group = self.hdf.get("SYNAPSE-SETTINGS")
            synapse_group.create_dataset(name=synapse_name, data=syn_data_arr)

        return

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
            self.sim_current_params["current type"] = "inhibitory"

        self.sim_current_params["start_t"] = start_t
        self.sim_current_params["duration"] = duration
        self.sim_current_params["current_A"] = current_A
        self.sim_current_params["frequency"] = frequency

        self.pulse_start_t_arr = []
        self.pulse_end_t_arr = []
        remaining_time = total_t - start_t
        time_intervals = remaining_time / frequency

        for i in range(frequency):
            pulse_start_t = time_intervals * (i) + start_t
            pulse_end_t = pulse_start_t + duration
            self.pulse_start_t_arr.append(pulse_start_t)
            self.pulse_end_t_arr.append(pulse_end_t)

        self.comp_arr[comp_num].dt = dt
        self.comp_arr[comp_num].set_current(current_type=self.current_type, current_amplitude=current_A)

    def run_simulation(self):
        """
        Main function to run the simulation.
        """
        self.start_t = time.time()
        self.interval_num = 1
        self.steps = 0
        self.run_t = 0
        while self.run_t < self.total_t:

            for a in self.comp_arr:

                a.step()  # step for each compartment

                if a.xflux_switch and (a.xflux_params["start_t"] <= self.run_t < a.xflux_params["end_t"]):
                    a.x_flux()

                if a.zflux_switch and \
                        (a.zflux_params["start_t"] <= self.run_t < a.zflux_params[
                            "end_t"]):
                    a.z_flux()

                if self.xoflux_switch and \
                        self.xoflux_params["start_t"] <= self.run_t < self.xoflux_params["end_t"]:
                    self.xoflux()

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
        self.end_t = time.time()

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

    def print_settings(self):
        print("HDF5 File name:" + self.file_name)
        print("==============================")
        print("Dendrite Morphology")
        print("==============================")
        for i in self.comp_arr:
            print(i.name + " Len(um): " + str(i.length*1e5) + " Rad(um): " + str(i.radius*1e5))
        print("")
        print("==============================")
        print("General settings")
        print("==============================")
        print("ED on: " + str(self.ED_on))
        print("Surface Area Constant: " + str(self.static_sa))
        print("ATPase rate Constant: "+ str(self.static_atpase))
        print("Total sim time(s):" + str(self.total_t))
        print("Time step(s): "+ str(self.dt))
        print("#Intervals recorded in file: " + str(self.intervals))
        print("Starting [Na]: " + str(self.comp_arr[1].na_i*1e3))
        print("Starting [K]: " + str(self.comp_arr[1].k_i * 1e3))
        print("Starting [Cl]: " + str(self.comp_arr[1].cl_i * 1e3))
        print("Starting [X]: " + str(self.comp_arr[1].x_i * 1e3))
        print("")
        print("==============================")
        print("[X]-flux settings")
        print("==============================")
        if self.xflux_count == 0:
            print("No [X]-flux")
        else:
            for i in range(len(self.xflux_names_arr)):
                print(self.xflux_names_arr[i])
                print("Flux rate: " + self.xflux_data_arr[1])
                print("z: " + self.xflux_data_arr[2])
                print("Start time: " + self.xflux_data_arr[3])
                print("End time: " + self.xflux_data_arr[4])
        print("")
        print("==============================")
        print("z-flux settings")
        print("==============================")
        if not self.sim_zflux_params["On"]:
            print("No z-flux")
        else:
            print(self.sim_zflux_params["Comp"])
            print("z end: " + str(self.sim_zflux_params["z_end"]))
            print("Start time: " + str(self.sim_zflux_params["start_t"]))
            print("End time: " + str(self.sim_zflux_params["end_t"]))
            print("Adjust x: " + str(self.sim_zflux_params["adjust_X"]))

        print("")
        print("==============================")
        print("Current settings")
        print("==============================")
        if self.sim_current_params["On"] == False:
            print("No external current added")
        else:
            print("Current onto: " + self.comp_arr[self.sim_current_params["compartment"].name])
            print("Current type: " + self.sim_current_params["current_type"])
            print("Start time: " + str(self.sim_current_params["start_t"]))
            print("Duration: " + str(self.sim_current_params["duration"]))
            print("Current Amplitude (A): " + str(self.sim_current_params["current_A"]))
            print("Current Frequency (Hz): " + str(self.sim_current_params["frequency"]) )

        print("")
        print("==============================")
        print("Synapse settings")
        print("==============================")
        if not self.syn_dict["On"]:
            print("No synapses added")
        else:
            print("Synapse on:" + self.comp_arr[self.syn_dict["compartment"]].name)
            if self.syn_dict["synapse_type"] ==0:
                print("Synapse type: Excitatory")
            else:
                print("Synapse type: Inhibitory")
            print("alpha rate constant: " + str(self.syn_dict["alpha"]))
            print("beta rate constant: "+ str(self.syn_dict["beta"]))
            print("start_t: " + str(self.syn_dict["start_t"]))
            print("end_t: " + str(self.syn_dict["end_t"]))
            print("max [NT](M) :" + str(self.syn_dict["max_neurotransmitter_conc"]))
            print("synaptic conductance: " + str(self.syn_dict["synapse_conductance"]))
        print("")
        print("==============================")
        print("HH settings")
        print("==============================")
        if not self.hh_on:
            print("No HH channels ")
        else:
            print("HH channels enabled")
        print("")
        print("===============================")