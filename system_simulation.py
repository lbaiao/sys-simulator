import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, 'D:/Google Drive/trabalho/mestrado/dev/link-simulator')
sys.path.insert(1, 'D:/Google Drive/trabalho/mestrado/dev/sys-simulator')
import random
import general as gen
import matplotlib.pyplot as plt
import numpy as np
import scipy.spatial as spatial
from devices import base_station, mobile_user, d2d_user, d2d_node_type
from channel import pathloss_channel
from enum import Enum
from timeit import default_timer as timer

class system_simulation_type_enum(Enum):
    """
    Enum that defines wich type of simulation will be run by the system_simulation object.
    VARYING_BS_RADIUS: multiple users, the bs radius changes on each iteration
    SINGLE_D2D_PAIR_SAME_NOISE: single mue, single d2d pair, the same noise is assumed to be on the MUE and DUE equipments
    SINGLE_D2D_PAIR_DIFFERENT_NOISE: single mue, single d2d pair, different noise signals, with the same probabilistic characteristics, are assumed to be on the MUE and DUE equipments
    """    
    VARYING_BS_RADIUS = 0,
    SINGLE_D2D_PAIR_SAME_NOISE = 1,
    SINGLE_D2D_PAIR_DIFFERENT_NOISE = 2,
    MULTIPLE_D2D_PAIRS = 3

class allocation_algorithm_enum(Enum):
    RANDOM = 0,
    LOWEST_AVAILABILITY = 1

class system_simulation:
    def __init__(self, sim_type, lower_lim_uc_db,lower_lim_d2d_db,sinr_margin,n_orthogonal_resources,n_mues,n_d2d,bs, **kwargs):
        if isinstance(sim_type, system_simulation_type_enum):
            self.simulation_type = sim_type
        else: raise Exception('Invalid system simulation type.')
        if sim_type == system_simulation_type_enum.VARYING_BS_RADIUS:
            self.run_loops = self.__run_loops_varying_bs_radius
        elif sim_type == system_simulation_type_enum.SINGLE_D2D_PAIR_SAME_NOISE or sim_type == system_simulation_type_enum.SINGLE_D2D_PAIR_DIFFERENT_NOISE:
            self.run_loops = self.__run_loops_single_pair
        elif sim_type == system_simulation_type_enum.MULTIPLE_D2D_PAIRS:
            self.run_loops = self.__run_loops_multiple_pairs
        self.lower_lim_uc_db = lower_lim_uc_db
        self.lower_lim_uc = 10**(lower_lim_uc_db/10)
        self.lower_lim_d2d_db = lower_lim_d2d_db
        self.lower_lim_d2d = 10**(lower_lim_d2d_db/10)
        self.sinr_margin = sinr_margin
        self.n_orthogonal_resources = n_orthogonal_resources
        self.n_mues = n_mues
        self.n_d2d = n_d2d
        self.bs = bs        
        if kwargs.get('allocation_algorithm') is not None:            
            if kwargs['allocation_algorithm'] == allocation_algorithm_enum.RANDOM:
                self.allocation_algorithm = self.__random_allocation
            if kwargs['allocation_algorithm'] == allocation_algorithm_enum.LOWEST_AVAILABILITY:            
                self.allocation_algorithm = self.__low_disponibility_priority
        else: self.allocation_algorithm = self.__random_allocation

    def get_snr_uc(self, p_uc, pathloss, noise):
        return p_uc*pathloss/noise

    def get_sinr_uc(self, mue, due, noise):
        return mue.tx_power*mue.pathloss_to_bs/(noise + due.tx_power*due.pathloss_to_bs)

    def is_snr_uc(self):
        pass

    def set_tx_power_mues_nodes(self,nodes_list, lower_snr_bound, noises, margin):
        for i in range(len(nodes_list)):
            tx_power = (1+margin)*lower_snr_bound*noises[i]/nodes_list[i].pathloss_to_bs
            nodes_list[i].set_tx_power(tx_power)    

    def set_tx_power_dues_nodes(self,nodes_list, lower_snr_bound, noises, margin):
        for i in range(len(nodes_list)):
            tx_power = margin*noises[i]/nodes_list[i].pathloss_to_bs
            nodes_list[i].set_tx_power(tx_power)    

    def allocate_resources_sequential(self,nodes, n_resources, **kwargs):
        #TODO: garantir que exista apenas 1 tabela de alocação no simulador inteiro??
        allocation_table = dict()
        n_nodes = len(nodes)
        for i in range(n_resources):
            if i>=n_nodes:
                break
            allocation_table[f'RB:{i}'] = nodes[i].id        
        return allocation_table

    # def allocate_resources_d2d_sequential(self,dues_txs,d2d_access_table,mues,mues_allocation_table):
    #     dues_allocation_table = dict()
    #     for node in dues_txs:
    #         pass

    def d2d_critical_distance(self,margin, channel, lower_d2d_sinr, lower_mue_sinr, uc_bs_distance, uc_d2d_distance, d2d_bs_distance):
        d_crit = d2d_bs_distance*(margin/(lower_d2d_sinr*(1+lower_mue_sinr*(uc_bs_distance/uc_d2d_distance)**channel.prop_coef*(1+margin))))**(1/channel.prop_coef)
        return d_crit

    def __run(self):
        # simulation parameters
        lower_lim_uc = 10**(self.lower_lim_uc_db/10) # sinr lower bound for MUEs in dB
        lower_lim_d2d = 10**(self.lower_lim_d2d_db/10)  # sinr lower bound for DUEs

        # channel specifications
        prop_coef = 3.5       # propagation coefficient
        channel = pathloss_channel(prop_coef)

        # noise specifications
        sigma = 1e-6        # awgn variance

        mues = [mobile_user(x) for x in range(self.n_mues)]    # mobile users
        dues = [d2d_user(x) for x in range(self.n_d2d)]   # d2d users

        nodes_bs_distances_table = gen.distribute_nodes(mues+dues,self.bs)     # distribute nodes on the map

        # calculate pathloss from BS to each other node
        for node in mues+dues:    
            node.set_pathloss_to_bs(
                channel.calculate_pathloss(
                    nodes_bs_distances_table[node.id]
                )
            )

        # calculate pathloss between d2d pairs
        d2d_distances_table = gen.get_distances_table(dues)
        d2d_pairs_table, _ = gen.get_d2d_links(d2d_distances_table, dues, channel)

        # set minimal allowed tx power for MUEs
        noises_mues = [n**2/2 for n in sigma*np.random.randn(len(mues))]
        self.set_tx_power_mues_nodes(mues, lower_lim_uc, noises_mues, self.sinr_margin)

        # allocate MUEs
        mues_allocation_table = self.allocate_resources_sequential(mues, self.n_orthogonal_resources)

        # allocate DUEs
        dues_ids = [id[0] for id in np.array(list(d2d_pairs_table.values()))[:,0]]
        dues_txs = [due for due in dues if due.id in dues_ids]
        d2d_request_access_table = self.allocate_resources_sequential(dues_txs, self.n_orthogonal_resources)
        noises_dues = [n**2/2 for n in sigma*np.random.randn(len(dues_txs))]
        # set max allowed tx power for DUEs
        self.set_tx_power_dues_nodes(dues_txs, lower_lim_d2d, noises_dues, self.sinr_margin)

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
            d2d_crit_distance_table[d.id] = self.d2d_critical_distance(self.sinr_margin,channel,lower_lim_d2d,lower_lim_uc,nodes_bs_distances_table[sharing_table[d.id]], d2d_mue_distances_table[d.id],nodes_bs_distances_table[d.id])
        dues_allocation_table_dcrit = dict()
        for d in dues_txs:
            dues_allocation_table_dcrit[d.id] = d2d_pairs_table[d.link_id][1] <=  nodes_bs_distances_table[d.id]

        d2d_allocation_rate_sinr = np.sum(list(dues_allocation_table.values()))/len(dues_allocation_table)
        d2d_allocation_rate_dcrit = np.sum(list(dues_allocation_table_dcrit.values()))/len(dues_allocation_table_dcrit)

        return d2d_allocation_rate_dcrit, d2d_allocation_rate_sinr

    def __run_single_pair(self, pair_distance):
        lower_lim_uc = 10**(self.lower_lim_uc_db/10) # sinr lower bound for MUEs in dB
        lower_lim_d2d = 10**(self.lower_lim_d2d_db/10)  # sinr lower bound for DUEs

        # channel specifications
        prop_coef = 3.5       # propagation coefficient
        channel = pathloss_channel(prop_coef)

        # noise specifications
        sigma = 1e-6        # awgn variance

        mues = [mobile_user(x) for x in range(self.n_mues)]    # mobile users
        dues = [d2d_user(x) for x in range(self.n_d2d)]   # d2d users

        _ = gen.distribute_nodes(mues,self.bs)     # distribute nodes on the map
        gen.distribute_pair_fixed_distance(dues, self.bs,pair_distance)        

        # calculate pathloss between d2d pairs
        d2d_pair_pathloss = channel.calculate_pathloss(pair_distance)        
        
        # set minimal allowed tx power for MUEs
        mue_bs_distance = spatial.distance.euclidean(mues[0].position, self.bs.position)
        noises_mues = [n**2/2 for n in sigma*np.random.randn(len(mues))]
        mues[0].set_pathloss_to_bs(channel.calculate_pathloss(mue_bs_distance))
        self.set_tx_power_mues_nodes(mues, lower_lim_uc, noises_mues, self.sinr_margin)

        if self.simulation_type == system_simulation_type_enum.SINGLE_D2D_PAIR_DIFFERENT_NOISE:
            noises_dues = [n**2/2 for n in [sigma*np.random.randn()]]
        else:
            noises_dues = noises_mues
            
        # set max allowed tx power for DUEs
        d2d_bs_distance = spatial.distance.euclidean(dues[0].position, self.bs.position)
        dues[0].set_pathloss_to_bs(channel.calculate_pathloss(d2d_bs_distance))        
        self.set_tx_power_dues_nodes([dues[0]], lower_lim_d2d, noises_dues, self.sinr_margin)

        # calculate sinr_dd
        d2d_mue_distance = spatial.distance.euclidean(mues[0].position, dues[0].position)
        d2d_mue_pathloss = channel.calculate_pathloss(d2d_mue_distance)
        sinr_dd = dues[0].tx_power*d2d_pair_pathloss/(noises_dues[0]+mues[0].tx_power*d2d_mue_pathloss)
        sinr_mue = self.get_sinr_uc(mues[0], dues[0], noises_mues[0])

        is_d2d_allocated_sinr = sinr_dd > lower_lim_d2d and sinr_mue > lower_lim_uc

        # critical distance method
        d2d_crit_distance = self.d2d_critical_distance(self.sinr_margin,channel,lower_lim_d2d,lower_lim_uc,mue_bs_distance, d2d_mue_distance, d2d_bs_distance)

        is_d2d_allocated_dcrit = pair_distance < d2d_crit_distance

        return is_d2d_allocated_dcrit, is_d2d_allocated_sinr

    def __critical_distance_allocation_matrix(self, mues, dues, channel):
        allocation_matrix = [ [ int( d.distance_d2d <= self.d2d_critical_distance(self.sinr_margin, channel, self.lower_lim_d2d, self.lower_lim_uc, m.distance_to_bs, spatial.distance.euclidean(m.position, d.position), d.distance_to_bs)) for d in dues ] for m in mues ]        
        return allocation_matrix

    def __random_allocation(self, allocation_matrix):
        allocation_matrix = np.array(allocation_matrix)
        n_mues = list(range(allocation_matrix.shape[0]))
        n_dues = list(range(allocation_matrix.shape[1]))
        allocation_vector = np.zeros(allocation_matrix.shape[1], dtype='int')
        for i in n_mues:
            for j in n_dues:
                if allocation_matrix[i][j] > 0:
                    allocation_vector[j] = 1
                    n_dues.pop(n_dues.index(j))
                    n_mues.pop(n_mues.index(i))
                    break
        return allocation_vector

    def __low_disponibility_priority(self, allocation_matrix):
        allocation_matrix = np.array(allocation_matrix)
        n_dues = allocation_matrix.shape[1]
        allocation_vector = np.zeros(allocation_matrix.shape[1], dtype='int')
        while(allocation_matrix.sum()>0):
            aux = allocation_matrix.sum(axis=0)
            min_index = list(aux).index(np.min(aux[aux>0]))
            aux2 = list(allocation_matrix[:,min_index])
            mue_index = next(aux2.index(j) for j in aux2 if j==1)
            allocation_vector[min_index] = 1
            allocation_matrix[mue_index,:] = 0
            allocation_matrix[:, min_index] = 0
        return allocation_vector

    def __run_multiple_pairs(self, pair_distances_list):
        lower_lim_uc = 10**(self.lower_lim_uc_db/10) # sinr lower bound for MUEs in dB
        lower_lim_d2d = 10**(self.lower_lim_d2d_db/10)  # sinr lower bound for DUEs

        # channel specifications
        prop_coef = 3.5       # propagation coefficient
        channel = pathloss_channel(prop_coef)

        # noise specifications
        sigma = 1e-6        # awgn variance        

        mues = [mobile_user(x) for x in range(self.n_mues)]    # mobile users

        dues_tx = list()
        dues_rx = list()

        for i in range(self.n_d2d//2):
            due_tx = d2d_user(i,type=d2d_node_type.TX) 
            due_rx = d2d_user(i,type=d2d_node_type.RX)   
            distance = random.choice(pair_distances_list)            
            due_tx.set_distance_d2d(distance)   
            due_rx.set_distance_d2d(distance)               
            dues_tx.append(due_tx) # d2d users tx   
            dues_rx.append(due_rx) # d2d users rx 

        gen.distribute_nodes(mues,self.bs)     # distribute nodes on the map

        gen.distribute_pair_fixed_distance_multiple(dues_tx, dues_rx, self.bs) 

        # set nodes distances to bs
        for n in dues_tx+dues_rx:
            n.set_distance_to_bs(spatial.distance.euclidean(n.position, self.bs.position))

        #TODO: calcular matriz de alocacao
        # critical distance method
        allocation_matrix = self.__critical_distance_allocation_matrix(mues, dues_tx, channel)

        allocation_vector = self.allocation_algorithm(allocation_matrix)
        
        dues_indexes = np.where(allocation_vector==1)[0]        
        distances = [d.distance_d2d for d in np.array(dues_tx)[dues_indexes]]
        distances_indexes = np.array(list(), dtype='int32')
        for i in distances:
            distances_indexes = np.concatenate((np.array(distances_indexes),np.where(pair_distances_list == i)[0]))

        distances_total_allocation = np.zeros(len(pair_distances_list))
        distances_appeareances = np.zeros(len(pair_distances_list))
        unique, counts = np.unique([list(pair_distances_list).index(d.distance_d2d) for d in dues_tx], return_counts=True) 
        aux = dict(zip(unique, counts))

        for i in aux.keys():
            distances_appeareances[i] = aux[i]

        for i in distances_indexes:
            distances_total_allocation[i] += 1

        return distances_total_allocation, distances_appeareances

    def __run_loops_multiple_pairs(self, n_loops, d2d_distances_list):
        start = timer()
        distances_allocation_ratio = np.zeros(len(d2d_distances_list))
        distances_appeareances = np.zeros(len(d2d_distances_list))
        for i in range(n_loops):
            print(f'\r  Loop {i+1}/{n_loops}',end='', flush=True)                
            aux1, aux2 = self.__run_multiple_pairs(d2d_distances_list)
            distances_allocation_ratio += aux1
            distances_appeareances += aux2
        print(f'\nElapsed time: {timer() - start}')        

        distances_allocation_ratio = [ distances_allocation_ratio[i]/distances_appeareances[i] for i in range(len(d2d_distances_list)) ]

        return distances_allocation_ratio, distances_appeareances

    def __run_loops_varying_bs_radius(self, n_loops, bs_radius_list):
        start = timer()
        allocation_rates = {'d2d_dcrit':list(), 'd2d_sinr':list()}
        for r in bs_radius_list:
            aux=list()
            self.bs.set_radius = r
            print(f'\nRadius: {r}m out of {bs_radius_list[-1]}m')            
            for j in range(n_loops):
                print(f'\r  Loop {j+1}/{n_loops}')                
                rate_d2d_dcrit, rate_d2d_sinr = self.__run()
                aux.append((rate_d2d_dcrit, rate_d2d_sinr))
            aux = np.array(aux)    
            allocation_rates['d2d_dcrit'].append(np.sum(aux[:, 0])/n_loops)
            allocation_rates['d2d_sinr'].append(np.sum(aux[:, 1])/n_loops)
        print(f'\nElapsed time: {timer() - start}')
        return allocation_rates
    
    def __run_loops_single_pair(self, n_loops, pair_distances_list):
        start = timer()
        allocation_rates = {'d2d_dcrit':list(), 'd2d_sinr':list()}
        for r in pair_distances_list:
            aux=list()            
            print(f'\nPair distance: {r}m out of {pair_distances_list[-1]}m')            
            for j in range(int(n_loops)):
                print(f'\r  Loop {j+1}/{n_loops}',end='', flush=True)                
                rate_d2d_dcrit, rate_d2d_sinr = self.__run_single_pair(r)
                aux.append((rate_d2d_dcrit, rate_d2d_sinr))
            aux = np.array(aux)    
            allocation_rates['d2d_dcrit'].append(np.sum(aux[:, 0])/n_loops)
            allocation_rates['d2d_sinr'].append(np.sum(aux[:, 1])/n_loops)
        print(f'\nElapsed time: {timer() - start}')
        return allocation_rates



        
        # test stuff
        # tx_x = [i.position[0] for i in dues_tx]
        # tx_y = [i.position[1] for i in dues_tx]
        # rx_x = [i.position[0] for i in dues_rx]
        # rx_y = [i.position[1] for i in dues_rx]
        # m_x = [i.position[0] for i in mues]
        # m_y = [i.position[1] for i in mues]

        # plt.plot(tx_x, tx_y, '*', label='dues tx')
        # plt.plot(rx_x, rx_y, 'd', label='dues rx')
        # plt.plot(m_x, m_y, 's', label='mues')
        # patch = plt.Circle(self.bs.position, self.bs.radius, edgecolor='red', facecolor='None', linewidth=5.0, zorder=10)
        # plt.legend()
        # ax = plt.gca()
        # ax.add_patch(patch)
        # plt.show()