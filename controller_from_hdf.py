"""
Version of the controller class to be able to start a simulation from where an HDF file started off
"""
import numpy as np

import compartment
import simulator3
import h5py

# Base file is the original file where the simulation was run

base_file_name = "T1_v2"
new_file_name = "Exp3-2"

comp_arr = []
with h5py.File(base_file_name, mode='r') as hdf:
    C = hdf.get('COMPARTMENTS')
    C_group_arr = []
    comp_names_arr = list(C.keys())

    T = hdf.get('TIMING')
    total_t = T.get('TOTAL_T')[()]
    intervals = T.get('INTERVALS')[()]
    dt = T.get("DT")[()]
    total_steps = total_t / dt
    interval_step = total_steps / intervals
    interval_arr = [round(interval_step * i) for i in range(intervals)]

    for i in comp_names_arr:
        temp_arr = []
        comp_name = i
        current_comp = C.get(comp_name)
        last_time_point = str(interval_arr[-1])
        last_dataset = current_comp.get(last_time_point)
        last_dataset = list(last_dataset)

        radius = last_dataset[1]
        length = last_dataset[2]
        na_i = last_dataset[4]
        k_i = last_dataset[5]
        cl_i = last_dataset[6]
        x_i = last_dataset[7]
        z_i = last_dataset[8]
        vm = last_dataset[-3]

        if comp_name == "Comp0(Soma)":
            comp = compartment.Compartment(comp_name, radius, length, static_sa=True, sa_value=2 * np.pi * 1e-5 * 40e-5)
            comp.set_ion_properties(na_i, k_i, cl_i, x_i, z_i, osmol_neutral_start=False, fixed_osm=False)
        else:
            comp = compartment.Compartment(comp_name, radius, length, static_sa=True,
                                           sa_value=2 * np.pi * 0.5e-5 * 20e-5)
            comp.set_ion_properties(na_i, k_i, cl_i, x_i, z_i, osmol_neutral_start=False)
        # comp.set_voltages(vm)
        comp_arr.append(comp)
    hdf.close()

sim = simulator3.Simulator(new_file_name)
for i in comp_arr:
    sim.add_compartment(i)
sim.set_electrodiffusion_properties(ED_on=True, diff_constant_dict={"na": (1.33 / 2) * 1e-7, "k": (1.96 / 2) * 1e-7,
                                                                    "cl": (2.03 / 2) * 1e-7, "x": 0})
sim.set_external_ion_properties()
sim.set_atpase_static(static_atpase=True)
sim.set_sa_static(static_sa=True)
#sim.set_z(comp="Comp8", z=-0.65, fixed_osm=True)
#sim.set_hh_on("Comp0(Soma)", t_on=1e-3)
total_t = 500
time_step = 1e-6
sim.set_timing(total_t=total_t, time_step=time_step, intervals=5000)
sim.set_zflux(comps=["Comp8"], start_t=100, end_t=150, z_end=-0.65, fixed_osm=False)
#sim.add_synapse("Comp8", "Excitatory", start_t=50e-3, duration=1e-3, max_neurotransmitter=4e-3, synapse_conductance=4e-9)
#sim.add_synapse("Comp4", "Inhibitory", start_t=52e-3, duration=1e-3, max_neurotransmitter=4e-3, synapse_conductance=4e-9)
#sim.add_current("Comp9", current_type="Excitatory", start_t=1e-3, duration=1e-3, current_A=0.1e-9, dt=time_step)
#sim.set_zflux(comps=["Comp8"], start_t=10, end_t=40, z_end=-1.05, fixed_osm=True)

# sim.set_xflux(comps=["Comp8"], start_t=20, end_t=40, flux_rate=3e-16)
# sim.add_current("Comp9", current_type="Excitatory", start_t=120e-3, duration=1e-3, current_A=0.1e-9, dt=time_step)
sim.run_simulation()
print("fin")
