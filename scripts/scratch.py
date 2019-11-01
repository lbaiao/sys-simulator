import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, 'D:/Google Drive/trabalho/mestrado/dev/link-simulator')
sys.path.insert(1, 'D:/Google Drive/trabalho/mestrado/dev/sys-simulator')
import general as gen
import matplotlib.pyplot as plt
import sim_methods as sim
import pygraphviz as pgv
import numpy as np
import scipy.spatial as spatial
from devices import base_station, mobile_user, d2d_user
from channel import pathloss_channel


# simulation parameters
lower_lim_uc_db = 10 # sinr lower bound for MUEs in dB
lower_lim_uc = 10**(lower_lim_uc_db/10) # sinr lower bound for MUEs in dB
lower_lim_d2d_db = 7  # sinr lower bound for DUEs in dB
lower_lim_d2d = 10**(lower_lim_d2d_db/10)  # sinr lower bound for DUEs
sinr_margin = 1     # safety margin 
n_orthogonal_resources = 200

# devices distribution 
n_mues = n_orthogonal_resources
n_d2d = 2*n_mues  # must be an even number
bs = base_station((0,0), radius=250)

# channel specifications
prop_coef = 3.5       # propagation coefficient
channel = pathloss_channel(prop_coef)

# noise specifications
sigma = 1e-6        # awgn variance

mues = [mobile_user(x) for x in range(n_mues)]    # mobile users
dues = [d2d_user(x) for x in range(n_d2d)]   # d2d users

nodes_bs_distances_table = gen.distribute_nodes(mues+dues,bs)     # distribute nodes on the map

# calculate pathloss from BS to each other node
for node in mues+dues:    
    node.set_pathloss_to_bs(
        channel.calculate_pathloss(
            nodes_bs_distances_table[node.id]
        )
    )

# calculate pathloss between d2d pairs
d2d_distances_table = gen.get_distances_table(dues)
d2d_pairs_table, d2d_pairs_pathloss_table = gen.get_d2d_links(d2d_distances_table, dues, channel)

# set minimal allowed tx power for MUEs
noises_mues = [n**2/2 for n in sigma*np.random.randn(len(mues))]
sim.set_tx_power_mues_nodes(mues, lower_lim_uc, noises_mues, sinr_margin)

# allocate MUEs
mues_allocation_table = sim.allocate_resources_sequential(mues, n_orthogonal_resources)

# allocate DUEs
dues_ids = [id[0] for id in np.array(list(d2d_pairs_table.values()))[:,0]]
dues_txs = [due for due in dues if due.id in dues_ids]
d2d_request_access_table = sim.allocate_resources_sequential(dues_txs, n_orthogonal_resources)
noises_dues = [n**2/2 for n in sigma*np.random.randn(len(dues))]
# set max allowed tx power for DUEs
sim.set_tx_power_dues_nodes(dues_txs, lower_lim_d2d, noises_dues, sinr_margin)

# mues_sharing = [m for m in mues if m.id in [mues_allocation_table[i] for i in list(d2d_request_access_table.keys())]]
sharing_table = dict()
d2d_mue_distances_table = dict()
for key in list(d2d_request_access_table.keys()):
    sharing_table[d2d_request_access_table[key]] = mues_allocation_table[key]    
for key in list(sharing_table.keys()):    
    d2d_mue_distances_table[key] = spatial.distance.euclidean(
        next(m for m in mues if m.id == sharing_table[key]).position,
        next(d for d in dues_txs if d.id == key).position
    )

# calculate sinr_dd
d2d_mue_pathloss_table = dict()
for key in list(d2d_mue_distances_table.keys()):
    d2d_mue_pathloss_table[key] = channel.calculate_pathloss(d2d_mue_distances_table[key])
sinr_dd_table = dict()
for key in list(d2d_mue_distances_table.keys()):
    due = next(d for d in dues_txs if d.id == key)
    mue = next(m for m in mues if m.id == sharing_table[key])
    i = int(due.id[due.id.index(':')+1:])
    sinr_dd_table[key] = due.tx_power*d2d_mue_distances_table[key]/(noises_dues[i]+mue.tx_power*d2d_mue_pathloss_table[key])

# d2d allocation table
dues_allocation_table = dict()
for key in list(sinr_dd_table.keys()):
    dues_allocation_table[key] = sinr_dd_table[key] > lower_lim_d2d    

# TODO: implementar verificação de enlace por meio da distância crítica e comparar com o método da SINR          

# critical distance method
d2d_crit_distance_table = dict()
for d in dues_txs:
    d2d_crit_distance_table[d.id] = sim.d2d_critical_distance(sinr_margin,channel,lower_lim_d2d,lower_lim_uc,nodes_bs_distances_table[sharing_table[d.id]], d2d_mue_distances_table[d.id],nodes_bs_distances_table[d.id])
dues_allocation_table_dcrit = dict()
for d in dues_txs:
    dues_allocation_table_dcrit[d.id] = d2d_pairs_table[d.link_id][1] <=  nodes_bs_distances_table[d.id]

d2d_allocation_rate_sinr = np.sum(list(dues_allocation_table.values()))/len(dues_allocation_table)
d2d_allocation_rate_dcrit = np.sum(list(dues_allocation_table_dcrit.values()))/len(dues_allocation_table_dcrit)

# plot stuff
print(f'D2D request access table:\n{d2d_request_access_table}')
print(f'Sharing table:\n{sharing_table}')
print(f'DUES allocation table:\n{dues_allocation_table}')
print(f'DUES allocation table (critical distance method):\n{dues_allocation_table_dcrit}')
print(f'D2D allocation rate: {d2d_allocation_rate_sinr}')
print(f'D2D allocation rate (d_crit): {d2d_allocation_rate_dcrit}')

mues_x = [i.position[0] for i in mues]
mues_y = [i.position[1] for i in mues]
d2d_x = [i.position[0] for i in dues]
d2d_y = [i.position[1] for i in dues]

p1 = plt.plot(mues_x, mues_y, '*', label='Usuários móveis')
p2 = plt.plot(d2d_x, d2d_y, 'd', label='Usuários D2D')
p3 = plt.plot(bs.position[0], bs.position[1], 'o', label='BS')

for i in d2d_pairs_table.values():
    due1 = d2d_user.get_due_by_id(dues, i[0][0])
    due2 = d2d_user.get_due_by_id(dues, i[0][1])
    pos_x = [due1.position[0], due2.position[0]]
    pos_y = [due1.position[1], due2.position[1]]
    plt.plot(pos_x, pos_y, 'k-')

plt.legend()
plt.show()
