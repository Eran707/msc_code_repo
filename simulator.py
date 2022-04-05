# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 18:10:40 2021

Controls the processes of the multi-compartmental model

@author: E Shorer

"""

import pandas as pd

import compartment
import electrodiffusion
import matplotlib
#matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

import seaborn as sns
import graphing as gr

global comp_arr
global df_sim
global df_start

df_sim = pd.DataFrame()
df_start =pd.DataFrame()

dt = 1e-3  # 1ms time steps
total_t = 80  # s
run_t = 0  # current simulation timing
t_arr = []

comp_arr = []
ed_arr = []  # array of all the electro-diffusion links between the compartments
ed_dict_arr = []  # array of all the electro-diffusion dictionaries (constantly changing)
ed_conc_changes_arr = []

# print(df_sim)

# example of how this simulation should be run
"""
comp_1 = compartment.Compartment("comp_1")
comp_arr.append(comp_1)
comp_1.set_ion_properties(x_i=250e-3)
df_sim[comp_1.name] = comp_1.get_df_array()
df_start[comp_1.name] = comp_1.get_df_array()

comp_2 = compartment.Compartment("comp_2")
comp_arr.append(comp_2)
comp_2.set_ion_properties(x_i=154e-3)
df_sim[comp_2.name] = comp_2.get_df_array()

comp_3 = compartment.Compartment("comp_3")
comp_arr.append(comp_3)
comp_3.set_ion_properties(x_i=154e-3)
df_sim[comp_3.name] = comp_3.get_df_array()
"""

df_sim.index = ['Radius', 'Length', 'Volume', 'Na_i', 'K_i', 'Cl_i', 'X_i', 'z_i', 'ATP pump rate', 'KCC2 pump rate',
                'Vm', 'Ek', 'ECl']
df_start.index = ['Radius', 'Length', 'Volume', 'Na_i', 'K_i', 'Cl_i', 'X_i', 'z_i', 'ATP pump rate', 'KCC2 pump rate',
                'Vm', 'Ek', 'ECl']
#print(df_sim)


##Linking all the compartments with electrodiffusion

def multi_comp_ed_link():
    for e in range(len(comp_arr) - 1):
        ed_arr.append(electrodiffusion.Electrodiffusion(comp_arr[e], comp_arr[e + 1]))


multi_comp_ed_link()
# #
run_t = 0
constant_j_atp = False


ED_ON = True

while run_t < total_t:

    if ED_ON:

        ed_dict_arr = []  # array of all the electro-diffusion dictionaries (constantly changing)
        ed_conc_changes_arr = []


        for a in range(len(comp_arr)):
            comp_arr[a].step(dt)  # step for each compartment
            comp_arr[a].x_flux()
            ed_dict_arr.append(comp_arr[a].get_ed_dict())  # electrodiffusion dictionary for each compartment

        for b in range(len(comp_arr) - 1):
            ed_conc_changes_arr.append(
                ed_arr[b].calc_ed(dt, ed_dict_arr[b], ed_dict_arr[b + 1]))  # makes an array of all the ED conc changes

        for c in range(len(comp_arr) - 1):
            comp_arr[c].ed_update(ed_conc_changes_arr[c],
                                  "positive")  # appending the electrodiffusion concentrations for each compartment
            comp_arr[c + 1].ed_update(ed_conc_changes_arr[c], "negative")

        for d in range(len(comp_arr)):
            comp_arr[d].update_volumes()  # updates of the volumes, arrays, and dataframe for each compartment
            if run_t != 0:
                comp_arr[d].update_arrays()
            df_sim[comp_arr[d].name] = comp_arr[d].get_df_array()

    else: # if you want to run with normal diffusion not ED
        for a in range(len(comp_arr)):
            comp_arr[a].step(dt)
            comp_arr[a].x_flux()
            comp_arr[a].update_volumes()  # updates of the volumes, arrays, and dataframe for each compartment
            comp_arr[a].update_arrays()
            df_sim[comp_arr[a].name] = comp_arr[d].get_df_array()

    if run_t != 0:
        t_arr.append(run_t)
    run_t += dt

def get_vm_dict():
    vm_dict ={}
    for x in range(len(comp_arr)):
       vm_dict[comp_arr[x].name] = comp_arr[x].v_arr
    return vm_dict


vm_dict = get_vm_dict()
g1 = gr.graph(start_time=20)
fig1,ax1 = g1.graph_time_vm_allcomps(vm_dict,t_arr)
sns.despine()
fig1.show()


"""
fig_vm, (a1, a2, a3) = plt.subplots(1, 3,sharey=True)
plt.xlabel("Time (s)")
plt.ylabel("Voltage (mV)")
a1.set_title("Comp 1")
a2.set_title("Comp 2")
a3.set_title("Comp 3")

a1.plot(t_arr[45000:-1], comp_1.v_arr[45000:-1])
a2.plot(t_arr[45000:-1], comp_2.v_arr[45000:-1])
a3.plot(t_arr[45000:-1], comp_3.v_arr[45000:-1])
sns.despine()


fig_cl, (b1, b2, b3) = plt.subplots(1, 3, sharey=True)
plt.xlabel("Time (s)")
plt.ylabel("Intracellular Chloride (mM)")
b1.set_title("Comp 1")
b2.set_title("Comp 2")
b3.set_title("Comp 3")

b1.plot(t_arr[45000:-1], comp_1.cl_arr[45000:-1])
b2.plot(t_arr[45000:-1], comp_2.cl_arr[45000:-1])
b3.plot(t_arr[45000:-1], comp_3.cl_arr[45000:-1])
sns.despine()

fig_vm.show()
fig_cl.show()

print(df_sim)
print(comp_1.v_arr[0:50])
"""