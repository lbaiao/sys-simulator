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
from system_simulation import system_simulation_type_enum, system_simulation, allocation_algorithm_enum
from devices import base_station, mobile_user, d2d_user
from channel import pathloss_channel


# simulation parameters
lower_lim_uc_db = 10 # sinr lower bound for MUEs in dB
lower_lim_d2d_db = 7  # sinr lower bound for DUEs in dB
sinr_margin = 1     # safety margin 
n_orthogonal_resources = 25

# devices distribution 
n_mues = n_orthogonal_resources
n_d2d = 2*n_mues  # must be an even number
bs = base_station((0,0), radius=500)

# channel specifications
prop_coef = 3.5       # propagation coefficient
channel = pathloss_channel(prop_coef)

# noise specifications
sigma = 1e-6        # awgn variance

# simulation
sim_type = system_simulation_type_enum.MULTIPLE_D2D_PAIRS

sim1 = system_simulation(sim_type, lower_lim_uc_db, lower_lim_d2d_db, sinr_margin, n_orthogonal_resources, n_mues, n_d2d, bs, allocation_algorithm=allocation_algorithm_enum.RANDOM)

sim2 = system_simulation(sim_type, lower_lim_uc_db, lower_lim_d2d_db, sinr_margin, n_orthogonal_resources, n_mues, n_d2d, bs, allocation_algorithm=allocation_algorithm_enum.LOWEST_AVAILABILITY)

# simulation execution
n_loops = 2.5e3
d2d_distances_list = np.linspace(0, bs.radius/2, 20)
allocation_avg_rnd, distances_appearances_rnd = sim1.run_loops(int(n_loops), d2d_distances_list)
allocation_avg_pbd, distances_appearances_pbd = sim2.run_loops(int(n_loops), d2d_distances_list)

plt.figure(1)
plt.plot(d2d_distances_list, allocation_avg_rnd, label='SR')
plt.plot(d2d_distances_list, allocation_avg_pbd, label='PBD')
plt.xlabel('Distância do enlace D2D [m]')
plt.ylabel('Razão de alocação')
plt.legend()

plt.figure(2)
plt.bar([2*(i+1) for i in range(len(d2d_distances_list))], distances_appearances_rnd, align='center', tick_label= [ f'{i:.3g}' for i in d2d_distances_list] )
plt.xlabel('Distâncias do enlace D2D [m]')
plt.ylabel('Ocorrências')
plt.title('SR')

plt.figure(3)
plt.bar([2*(i+1) for i in range(len(d2d_distances_list))], distances_appearances_pbd, align='center', tick_label= [ f'{i:.3g}' for i in d2d_distances_list] )
plt.xlabel('Distâncias do enlace D2D [m]')
plt.ylabel('Ocorrências')
plt.title('PBD')

plt.show()


