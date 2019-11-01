# SINGLE PAIR SIMULATION SCRIPT WITH THE SAME NOISE AT MUE AND D2D

import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, 'D:/Google Drive/trabalho/mestrado/dev/link-simulator')
sys.path.insert(1, 'D:/Google Drive/trabalho/mestrado/dev/sys-simulator')
import general as gen
import matplotlib.pyplot as plt
# import pygraphviz as pgv
import numpy as np
import scipy.spatial as spatial
from system_simulation import system_simulation_type_enum, system_simulation
from devices import base_station, mobile_user, d2d_user
from channel import pathloss_channel


# simulation parameters
lower_lim_uc_db = 10 # sinr lower bound for MUEs in dB
lower_lim_d2d_db = 7  # sinr lower bound for DUEs in dB
sinr_margin = 1     # safety margin 
n_orthogonal_resources = 1

# devices distribution 
n_mues = n_orthogonal_resources
n_d2d = 2*n_mues  # must be an even number
bs = base_station((0,0), radius=2000)

# channel specifications
prop_coef = 3.5       # propagation coefficient
channel = pathloss_channel(prop_coef)

# noise specifications
sigma = 1e-6        # awgn variance

# simulation
sim_type = system_simulation_type_enum.SINGLE_D2D_PAIR_DIFFERENT_NOISE
sim = system_simulation(sim_type, lower_lim_uc_db, lower_lim_d2d_db, sinr_margin, n_orthogonal_resources, n_mues, n_d2d, bs)

# simulation execution
n_loops = 3e3
radius = [r for r in range(10, 1500, 100)]
allocation_rates = sim.run_loops(int(n_loops), radius)

# plot stuff
plt.plot(radius,allocation_rates['d2d_dcrit'], '-o', label='Critical distance method')
plt.plot(radius,allocation_rates['d2d_sinr'], '-d', label='SINR method')
plt.xlabel('Radius [m]')
plt.ylabel('Allocation ratio')
plt.title('Single pair')
plt.legend()
plt.show()
