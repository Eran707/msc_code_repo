"""
Controller class interfaces with simulator_from_hdf to amend HDF simulation files
A) "Extend" a previous simulation
B) "Resume" a previous simulation
C) "LastValues" starts a simulation based on the previous values

Creator: EFS
Revised: April 2022
"""

import simulator_from_hdf

old_file_name = "SS_gk80_HH"
new_file_name = "Exp6_gk80_65"
amend_type = 'LastValues'

sim = simulator_from_hdf.SimulatorFromHDF(old_file_name, new_file_name, amend_type, already_extended=False)
sim.set_timing(total_t=200e-3, dt=1e-6, intervals=10000)
# sim.set_timing(extend_t=100)
sim.set_electrodiffusion_properties(ED_on=True)
sim.set_external_ion_properties()
sim.set_sa()
sim.set_static_sa()
sim.set_na_i_start()
sim.set_atpase_static(static_atpase=True)
sim.set_gkcc2(comp_name="Comp4", gkcc2=8e-3)
sim.set_hh_on("Comp0(Soma)", t_on=0)
sim.add_synapse("Comp8", "Excitatory", start_t=50e-3, duration=5 * 1e-3,
                max_neurotransmitter=3e-3, synapse_conductance=6.5e-9)
sim.add_synapse("Comp4", "Inhibitory", start_t=(50.1)*1e-3, duration=5 * 1e-3,
                max_neurotransmitter=3e-3, synapse_conductance=10e-9)

sim.write_settings_to_file(old_file_name)
sim.run_simulation()
