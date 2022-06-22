"""
Main script to control the neural simulation.

# STEP 1: INITIALIZE SIMULATION AND DEFINE MORPHOLOGY:

    - Define the file name (simulation name) where the results of the simulation will be saved to
    - The file will be saved to the root directory
    - Either add the default multicompartment which adds a set number of defined compartments
        This function allows an option to add a somatic compartment which creates a compartment that is double the size
    - Or, add individual compartments

# STEP 2: SIMULATION SETTINGS:

    - set_electrodiffusion_properties: turn electrodiffusion on/off and define diffusion constants
    - set_external_ion_properties: option to change external ion concentrations away from default
    - set_atpase_static: use the initialization Na concentration to fix the Na-K-ATPase pump rate
    - set_sa_static: ensure that the surface area of each compartment does not change
    - set_hh: add hodgkin huxley channels to a compartment
    - set_timing: define the time step, total time, and intervals (number of times to save to file)

# STEP 3: ADDITIONAL ION FLUXES SETTINGS:

    - set_zflux: change impermeant anion valence
    - set_xflux: change impermeant anion quantity
    - add_current: add external current source
    - add_synapse: add a synapse

# STEP 4: RUN SIMULATION
    - write_settings_to_file: option to make a text file in the root directory which contains the simulation settings
    - run_simulation: runs the simulation

Creator: EF Shorer
"""

###################################################################
# 1) DEFINE SIMULATOR CLASS AND DEFINE MORPHOLOGY
###################################################################
import simulator
file_name = "SS_gk80_HH"

sim = simulator.Simulator(file_name)
sim.add_default_multicompartment(number_of_comps=9, rad=0.5e-5, len=20e-5, soma=True)
#sim.add_compartment()

##################################################################
# 2) SET SIMULATION SETTINGS
##################################################################
sim.set_electrodiffusion_properties(ED_on=True, diff_constant_dict={"na": (1.33 / 2) * 1e-7, "k": (1.96 / 2) * 1e-7,
                                                                    "cl": (2.03 / 2) * 1e-7, "x": 0})
sim.set_external_ion_properties()
sim.set_atpase_static(static_atpase=True)
sim.set_sa_static(static_sa=True)
# sim.set_hh_on(comp="Comp0(Soma)", t_on=0)

total_t = 30
time_step = 1e-6
sim.set_timing(total_t=total_t, time_step=time_step, intervals=1000)

##################################################################
# 3) ADDITIONAL ION FLUX SETTINGS
##################################################################

#sim.set_zflux(comps=["Comp4"], start_t=80e-4, end_t=120e-4, z_end=-0.65, adjust_x=True)
#sim.set_xflux(comps=['Comp8'], start_t=0.5, end_t=0.7, flux_rate=20e-16 / 60, z=-0.85)
#sim.add_synapse("Comp1", "Excitatory", start_t=0.5, duration=5 * 1e-3, max_neurotransmitter=1e-3, synapse_conductance=1e-9)
#sim.add_current("Comp9","Excitatory",start_t=10e-3,duration=1e-3, current_A=0.1e-9, total_t=total_t,frequency=1)



##################################################################
# 4) RUN SIMULATION
##################################################################

sim.write_settings_to_file()
sim.run_simulation()
print("fin")


