"""
Class to calculate the synapse concentration changes
Each synapse is an object
"""
import numpy as np
from constants import F


class synapse:

    def __init__(self, comp_name='', synapse_type='Inhibitory', start_t=0, duration=2e-3, max_neurotransmitter=1e-3):
        self.comp_name = comp_name
        self.synapse_type = synapse_type
        self.start_t = start_t,
        self.duration = duration
        self.nt_max = max_neurotransmitter
        self.alpha = 0.5  # ms-1.mM-1 == Forward rate constant
        self.beta = 0.1  # ms-1 == Backward rate constant
        self.r_initial = 0  # ratio of NT bound initially
        self.r_t = 0  # current ratio of NT bound
        self.r_infinity = (self.alpha * self.nt_max) / (self.alpha * self.nt_max + self.beta)
        self.tau = 1 / (self.alpha * self.nt_max + self.beta)
        self.g_synapse = 1 * 10 - 9

    def synapse_step(self, dt, run_t, V_membrane=0, V_ion=0):
        self.r_t += self.r_infinity + (self.r_initial - self.r_infinity) * np.e ** (-(run_t - self.start_t) / self.tau)

        I_syn = self.g_synapse * self.r_t * (V_membrane - V_ion)
        if self.synapse_type == 'Inhibitory':
            I_syn = I_syn * 4 / 5  # CL- only contributes about 80% of the GABA current, HCO3- contributes the rest.

        I_syn = I_syn / F  # converting coloumb to mol
        I_syn = I_syn * dt  # getting the mol input for the timestep

        return I_syn
