"""
Controller class interfaces with simulator_from_hdf to amend HDF simulation files
A) "Extend" a previous simulation
B) "Resume" a previous simulation
C) "LastValues" starts a simulation based on the previous values

Creator: EFS
Revised: April 2022
"""

import simulator_from_hdf


old_file_name = "Steady_105_HH"
new_file_name = "Steady_105_HH_v2"
amend_type = 'LastValues'

sim = simulator_from_hdf.SimulatorFromHDF(old_file_name, new_file_name, amend_type, already_extended= False)
sim.set_timing(total_t=150, dt=1e-6, intervals=1000)
#sim.set_timing(extend_t=100)
sim.set_electrodiffusion_properties(ED_on=True)
sim.set_external_ion_properties()
sim.set_static_sa()
sim.set_atpase_static(static_atpase=True)
#sim.add_synapse("Comp8", "Excitatory", start_t=150e-3, duration=5 * 1e-3, max_neurotransmitter=3e-3, synapse_conductance=6e-9)
#sim.add_synapse("Comp4", "Inhibitory", start_t=(150.1)*1e-3, duration=5 * 1e-3, max_neurotransmitter=3e-3, synapse_conductance=10e-9)
sim.set_hh_on("Comp0(Soma)", t_on=0)

#sim.set_gkcc2(gkcc2=2e-3)

#sim.comp_arr[7].adjust_x_bool = True
#sim.add_current(comp_name='Comp9', current_type="Excitatory", start_t=1e-3, duration=1e-3, current_A=0.1e-9)
sim.write_settings_to_file(old_file_name)
sim.run_simulation()

