"""
Main script to run simulation
"""

import simulator3

file_name = "T7_v2"

# 1) DEFINE SIMULATOR CLASS AND ADD COMPARTMENTS
sim = simulator3.Simulator(file_name)

sim.add_default_multicompartment(number_of_comps=9, rad=0.5e-5, len=20e-5, soma=True)
#sim.set_z("Comp0(Soma)", z=-0.65, adjust_x=True)
# 2) SET SIMULATION SETTINGS

sim.set_electrodiffusion_properties(ED_on=True, diff_constant_dict={"na": (1.33 / 2) * 1e-7, "k": (1.96 / 2) * 1e-7,
                                                                     "cl": (2.03 / 2) * 1e-7, "x": 0})
sim.set_external_ion_properties()
sim.set_atpase_static(static_atpase=True)
sim.set_sa_static(static_sa=True)
sim.set_zflux(comps=["Comp8"], start_t=350, end_t=400, z_end=-0.65, adjust_x = True)
total_t = 750
time_step = 1e-6
sim.set_timing(total_t=total_t, time_step=time_step, intervals=10000)
#sim.set_hh_on(comp="Comp0(Soma)", t_on=0)
sim.print_settings()
sim.run_simulation()
print("fin")


#sim.add_synapse("Comp1", "Excitatory", start_t=0.5, duration=5 * 1e-3, max_neurotransmitter=1e-3, synapse_conductance=1e-9)
# sim.add_current("Comp9","Excitatory",start_t=10e-3,duration=1e-3, current_A=0.1e-9,
# total_t=total_t,frequency=1)

# sim.set_xflux(comps=["Comp8"], start_t=120, end_t=180, z=-0.85, flux_rate=300*1e-3/60)


# sim.set_zflux()
# sim.set_xoflux()

# 4) RUN SIMULATION
