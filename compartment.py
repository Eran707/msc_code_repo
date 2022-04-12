# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 17:45:32 2021

@author: eshor

Class which defines the compartments object and related methods.

Class: Compartment : New compartment Methods: __int__ : Initializes compartment object set_ion_properties: define
intracellular ionic properties of the compartment (extracellular properties are imported) step: actions to take at
every time point in the simulation update_volumes: change the volume of the compartment (as well as ion
concentrations based on new volume) update_arrays: update the arrays for each parameter of the compartment ed_update:
make changes to the compartment based on the results of electrodiffusion get_ed_dict: sends the current status of the
compartment back to the simulation to be evaluated by electrodiffusion equations get_fin_vals: sends the final values
of the compartment back to the simulation get_df_array: sends the dataframe arrays back to the simulation


"""

##################################################################################
# IMPORTS


import numpy as np

from common import \
    gk, gna, gcl, p_atpase, p_kcc2, \
    pw, vw, RTF, cm, F, \
    hh_v, gna_bar, gk_bar, \
    m, b


##################################################################################
# COMPARTMENT CLASS

class Compartment:

    def __init__(self, compartment_name, radius=1e-5, length=10e-5, static_sa=True, sa_value=0, hh_on=False):

        self.static_sa = static_sa
        self.name = compartment_name
        self.radius = radius  # in dm
        self.length = length  # in dm
        self.w = np.pi * (self.radius ** 2) * self.length
        self.w_initial = self.w
        self.dw, self.w2 = 0, 0
        if static_sa and sa_value != 0:
            self.sa = sa_value
        else:
            self.sa = 2 * np.pi * self.radius * self.length

        self.hh_on = hh_on
        self.hh_n = 0
        self.hh_m = 0
        self.hh_h = 0
        self.hh_g_k = 0
        self.hh_d_k = 0
        self.hh_g_na = 0
        self.hh_d_na = 0

        self.static_atpase = True
        self.FinvC = F / cm

        self.j_kcc2 = 0
        self.j_p = 0

        self.p_kcc2 = p_kcc2
        self.p_atpase = p_atpase

        self.v, self.E_cl, self.E_k, self.E_na, self.drivingf_cl = 0, 0, 0, 0, 0
        self.na_i, self.k_i, self.cl_i, self.x_i, self.z_i, self.osm_i = 0, 0, 0, 0, 0, 0
        self.na_o, self.k_o, self.cl_o, self.osm_o = 0, 0, 0, 1

        self.na_i_start, self.x_start, self.z_start = 0, 0, 0

        self.x_default = 154.962e-3
        self.z_default = -0.85
        self.adjust_x_bool = False
        self.xflux_setup, self.zflux_setup, self.external_xflux_setup = True, True, True
        self.xflux_switch, self.zflux_switch = False, False
        self.xflux, self.xoflux = 0, 0
        self.t_xflux = 0
        self.xflux_params = {"start_t": 0, "end_t": 0,
                             "flux_rate": 60, "z": 0}
        self.zflux_params = {"start_t": 0, "end_t": 0, "z": 0}
        self.dt_xflux, self.flux_points, self.alpha, self.beta = 0, 0, 0, 0
        self.d_xflux, self.d_zflux, self.total_x_flux, self.static_xflux, self.x_final = 0, 0, 0, 0, 0
        self.osmo_final = 0
        self.z_diff, self.z_final, self.z_diff, self.z_inc, self.zflux = 0, 0, 0, 0, 0
        self.synapse_on = False
        self.current_on = False
        # Zeroing Delta values
        self.d_na_i, self.d_na_atpase, self.d_na_leak, self.d_na_hh = 0, 0, 0, 0
        self.d_k_i, self.d_k_atpase, self.d_k_leak, self.d_k_kcc2, self.d_k_hh = 0, 0, 0, 0, 0
        self.d_cl_i, self.d_cl_leak, self.d_cl_kcc2 = 0, 0, 0
        self.d_x_i = 0

        self.synapse_dict = {}
        self.molar_per_s = 0
        self.current_on = False
        self.current_type = 0

        self.r, self.d_r = 0, 0

        self.dt, self.syn_t_on, self.syn_t_off = 0, 0, 0  # Timing

    def set_ion_properties(self,
                           na_i=0.014,
                           k_i=0.1229,
                           cl_i=0.0052,
                           x_i= 0.1549,
                           z_i=-0.85,
                           osmol_neutral_start=False):
        """
        - Adjustment of starting concentrations to ensure starting electroneutrality
       old defaults: na_i=14.001840415288e-3, k_i=122.870162657e-3, cl_i=5.1653366e-3,
                           x_i=154.972660318083e-3, z_i=-0.85

                            na_i=0.013614929701933197,
                           k_i=0.12323046997950128,
                           cl_i=0.005161213856695764,
                           x_i= 0.15499338647116973,
                           z_i=-0.85,
                           osmol_neutral_start=False)
        """
        self.na_i, self.k_i, self.cl_i, self.x_i, self.z_i = na_i, k_i, cl_i, x_i, z_i  # Intracellular ion conc.
        self.na_i_start, self.x_start, self.z_start = na_i, x_i, z_i
        self.na_o = 145e-3
        self.k_o = 3.5e-3
        self.cl_o = 119e-3
        self.x_o = 29.5e-3
        self.osm_o = self.na_o + self.k_o + self.cl_o + self.x_o
        self.osm_initial = self.na_i + self.k_i + self.cl_i + self.x_i

        if osmol_neutral_start:
            self.k_i = self.cl_i - self.z_i * self.x_i - self.na_i
            # self.cl_i = self.k_i +self.z_i*self.x_i +self.na_i

        self.v = self.FinvC * (self.w / self.sa) * (
                    self.na_i + self.k_i + (self.z_i * self.x_i) - self.cl_i)  # starting membrane potential is -72.6mV
        self.E_na = -1 * RTF * np.log(self.na_i / self.na_o)
        self.E_k = -1 * RTF * np.log(self.k_i / self.k_o)
        self.E_cl = RTF * np.log(self.cl_i / self.cl_o)
        self.drivingf_cl = self.v - self.E_cl

    def set_osmo_neutral_start(self):
        self.k_i = self.cl_i - self.z_i * self.x_i - self.na_i

    def set_voltages(self, vm=72.6 * 1e-3):
        self.v = vm

    def set_synapse(self, syn_dict, dt=1e-6):
        self.synapse_on = True
        self.synapse_dict = syn_dict
        self.r = 0  # proportion of bound synapses
        self.d_r = 0  # delta r
        self.dt = dt

    def synapse_step(self, run_t):

        if self.synapse_dict['start_t'] < run_t < self.synapse_dict['end_t']:
            nt_conc = self.synapse_dict["max_neurotransmitter_conc"]
        else:
            nt_conc = 0

        self.d_r = (self.synapse_dict["alpha"] * nt_conc * (1 - self.r) - self.synapse_dict["beta"] * self.r) * self.dt
        self.r = self.r + self.d_r

        if self.synapse_dict["synapse_type"] == 0:  # excitatory synapse
            i_syn = self.synapse_dict["synapse_conductance"] * self.r * (self.E_na - self.v)
            i_syn = i_syn / F  # converting coloumb to mol
            i_syn = i_syn * self.dt  # getting the mol input for the timestep
            na_entry = i_syn / self.w
            self.na_i += na_entry

        elif self.synapse_dict["synapse_type"] == 1:  # inhibitory synapse
            i_syn = self.synapse_dict["synapse_conductance"] * self.r * (self.v - self.E_cl)
            i_syn = i_syn * 4 / 5  # CL- only contributes about 80% of the GABA current, HCO3- contributes the rest.
            i_syn = i_syn / F  # converting coloumb to mol
            i_syn = i_syn * self.dt  # getting the mol input for the timestep
            cl_entry = i_syn / self.w
            self.cl_i += cl_entry

    def set_current(self, current_type=1, current_amplitude=1e-9):

        self.molar_per_s = current_amplitude / F
        self.current_on = True
        self.current_type = current_type

    def current_step(self, dt=1e-6):

        if self.current_type == 0:  # Inihibitory current
            cl_current = self.molar_per_s / self.w
            cl_current = cl_current * dt
            self.cl_i += cl_current

        elif self.current_type == 1:  # Excitatory current
            na_current = self.molar_per_s / self.w
            na_current = na_current * dt
            self.na_i += na_current

    def step(self):
        """
        Perform a time step for the specific compartment.
        """
        # 1) Zeroing deltas
        self.d_na_i, self.d_k_i, self.d_cl_i, self.d_x_i = 0, 0, 0, 0

        # 2) Updating voltages
        self.v = self.FinvC * (self.w / self.sa) * (self.na_i + self.k_i + (self.z_i * self.x_i) - self.cl_i)
        self.E_na = -1 * RTF * np.log(self.na_i / self.na_o)
        self.E_k = -1 * RTF * np.log(self.k_i / self.k_o)
        self.E_cl = RTF * np.log(self.cl_i / self.cl_o)
        self.drivingf_cl = self.v - self.E_cl

        # 3) Update ATPase, KCC2 pump rate, assess HH if soma
        if not self.static_atpase:
            self.j_p = self.p_atpase * (self.na_i / self.na_o) ** 3
        elif self.static_atpase:
            self.j_p = self.p_atpase * (self.na_i_start / self.na_o) ** 3

        self.j_kcc2 = self.p_kcc2 * (self.E_k - self.E_cl)

        if self.hh_on:
            self.hh_g_na = self.calc_hh_g_na()
            self.hh_g_k = self.calc_hh_g_k()

        # 4) Solve ion flux equations for t+dt from t
        self.d_na_leak = - self.dt * self.sa / self.w * gna * (self.v + RTF * np.log(self.na_i / self.na_o))
        self.d_na_atpase = - self.dt * self.sa / self.w * (+3 * self.j_p)
        if self.hh_on:
            self.d_na_hh = -self.dt * self.sa / self.w * self.hh_g_na * (self.v + RTF * np.log(self.na_i / self.na_o))
            self.d_na_i = self.d_na_leak + self.d_na_atpase + self.d_na_hh  # inward na HH current
        else:
            self.d_na_i = self.d_na_leak + self.d_na_atpase

        self.d_k_leak = - self.dt * self.sa / self.w * gk * (self.v + RTF * np.log(self.k_i / self.k_o))
        self.d_k_atpase = - self.dt * self.sa / self.w * (- 2 * self.j_p)
        self.d_k_kcc2 = - self.dt * self.sa / self.w * (- self.j_kcc2)
        if self.hh_on:
            self.d_k_hh = -self.dt * self.sa / self.w * self.hh_g_k * (
                        self.v + RTF * np.log(self.k_i / self.k_o))  ## outward K+ current
            self.d_k_i = self.d_k_leak + self.d_k_atpase + self.d_k_kcc2 + self.d_k_hh
        else:
            self.d_k_i = self.d_k_leak + self.d_k_atpase + self.d_k_kcc2

        self.d_cl_leak = + self.dt * self.sa / self.w * gcl * (self.v + RTF * np.log(self.cl_o / self.cl_i))
        self.d_cl_kcc2 = + self.dt * self.sa / self.w * self.j_kcc2
        self.d_cl_i = self.d_cl_leak + self.d_cl_kcc2

        if self.cl_i < 0:
            print("Cl_i = " + str(self.cl_i))
            print("d_Cl_i = " + str(self.d_cl_i))
            raise Exception("chloride log can't have a negative number")

        if self.k_i < 0:
            print("k_i = " + str(self.k_i))
            print("d_k_i = " + str(self.d_k_i))
            raise Exception("[K+] <0 --  log can't have a negative number")

        # 5) Update ion concentrations
        self.na_i = self.na_i + self.d_na_i
        self.k_i = self.k_i + self.d_k_i
        self.cl_i = self.cl_i + self.d_cl_i

        # Note that X_i values may change based on wether there is adjustment with z change

    def update_volumes(self,avg_osmo = 0):
        """ Calculates the new compartment volume (dm3)
        Elongation should occur radially
        """
        if self.adjust_x_bool:
            self.osm_i = avg_osmo
            self.x_i = self.osm_i - self.na_i - self.k_i - self.cl_i
        else:
            self.osm_i = self.na_i + self.k_i + self.cl_i + self.x_i

        self.dw = self.dt * (vw * pw * self.sa * (self.osm_i - self.osm_o))
        self.w2 = self.w + self.dw

        self.na_i = self.na_i * self.w / self.w2
        self.k_i = self.k_i * self.w / self.w2
        self.cl_i = self.cl_i * self.w / self.w2
        self.x_i = self.x_i * self.w/self.w2

        self.radius = np.sqrt(self.w2 / (self.length * np.pi))

        if self.static_sa:
            self.sa = self.sa
        else:
            self.sa = 2 * np.pi * self.radius * self.length

        self.w = self.w2

    def ed_update(self, ed_change: dict, sign="positive"):
        """
        Receives a dictionary and update
        """

        if sign == "positive":
            self.na_i += (ed_change["na"] / self.length) / self.w
            self.cl_i += (ed_change["cl"] / self.length) / self.w
            self.k_i += (ed_change["k"] / self.length) / self.w
            self.x_i += (ed_change["x"] / self.length) / self.w
        elif sign == "negative":
            self.na_i -= (ed_change["na"] / self.length) / self.w
            self.cl_i -= (ed_change["cl"] / self.length) / self.w
            self.k_i -= (ed_change["k"] / self.length) / self.w
            self.x_i -= (ed_change["x"] / self.length) / self.w

    def get_ed_dict(self):
        ed_dict = {"na": self.na_i, "k": self.k_i, "cl": self.cl_i, "x": self.x_i, "Vm": self.v}
        return ed_dict

    def get_df_dict(self, time=0):
        df_dict = {"time": time, "name": self.name, "radius": self.radius, "length": self.length, "volume": self.w,
                   "na": self.na_i, "k": self.k_i,
                   "cl": self.cl_i, "x": self.x_i, "z": self.z_i, "vm": self.v, "e_k": self.E_k, "e_cl": self.E_cl}
        return df_dict

    def get_array(self, time=0):
        array = [time, self.radius, self.length, self.w,
                 self.na_i, self.k_i, self.cl_i, self.x_i, self.z_i,
                 self.d_na_i, self.d_na_leak, self.d_na_atpase,
                 self.d_k_i, self.d_k_leak, self.d_k_atpase, self.d_k_kcc2,
                 self.d_cl_i, self.d_cl_leak, self.d_cl_kcc2,
                 self.v, self.E_k, self.E_cl]
        return array

    def x_flux(self):
        """
        FLUX IMPERMEANTS INTO THE COMPARTMENT

        """
        if self.xflux_setup:
            # starting values for flux
            self.d_xflux = (self.xflux_params["flux_rate"]) * self.dt  # in MOLES PER TIME STEP
            self.x_i = self.x_i + (self.d_xflux / self.w)
            self.xflux_setup = False
        else:
            self.x_i = self.x_i + (self.d_xflux / self.w)

    def z_flux(self):
        """
        Changing the charge of intra-compartmental impermeants during the simulation
        """
        if self.zflux_setup:
            self.z_diff = self.zflux_params["z"] - self.z_i
            # self.x_mol_start =  self.w * self.x_i
            t_diff = (self.zflux_params["end_t"] - self.zflux_params["start_t"]) / self.dt
            self.z_inc = self.z_diff / t_diff
            self.zflux_setup = False
        else:
            self.z_i += self.z_inc

    def adjust_x(self):
        """
        X = z*m + b
        z is impermeant anion charge
        m and b are constants derived from the best fit curve for the relationship between x and z
            constants found in common.py
        """
        self.x_i = self.z_i * m + b

    def calc_hh_g_k(self):

        # alpha_n
        numerator = -self.v + hh_v[2]
        denominator = hh_v[1] * (np.exp((-self.v + hh_v[2]) / hh_v[3]) - 1)
        alpha_n = numerator / denominator

        # beta_n
        beta_n = hh_v[4] * np.exp((-self.v + hh_v[5]) / hh_v[6])

        # calculating n
        dn = alpha_n * (1 - self.hh_n) - beta_n * self.hh_n
        self.hh_n = self.hh_n + dn * self.dt

        # calculating gk
        self.hh_g_k = gk_bar * self.hh_n ** 4

        return self.hh_g_k

    def calc_hh_g_na(self):

        # alpha_m
        numerator = -self.v + hh_v[8]
        denominator = hh_v[7] * (np.exp((-self.v + hh_v[8]) / hh_v[9]) - 1)
        alpha_m = numerator / denominator

        # beta_m
        beta_m = hh_v[10] * np.exp((-self.v + hh_v[11]) / hh_v[12])

        # calculating m
        dm = alpha_m * (1 - self.hh_m) - beta_m * self.hh_m
        self.hh_m = self.hh_m + dm * self.dt

        #  alpha_h

        alpha_h = hh_v[13] * (np.exp((-self.v + hh_v[14]) / hh_v[15]))

        # beta_h

        beta_h = 1 / (hh_v[16] * (np.exp((-self.v + hh_v[17]) / hh_v[18]) + 1))

        # calculating h
        dh = alpha_h * (1 - self.hh_h) - beta_h * self.hh_h
        self.hh_h = self.hh_h + dh * self.dt

        # calculating gna
        self.hh_g_na = (self.hh_m ** 3) * self.hh_h * gna_bar
        return self.hh_g_na
