import numpy as np

def get_snr_uc(p_uc, pathloss, noise):
    return p_uc*pathloss/noise

def is_snr_uc():
    pass

def set_tx_power_mues_nodes(nodes_list, lower_snr_bound, noises, margin):
    for i in range(len(nodes_list)):
        tx_power = (1+margin)*lower_snr_bound*noises[i]/nodes_list[i].pathloss_to_bs
        nodes_list[i].set_tx_power(tx_power)    

def set_tx_power_dues_nodes(nodes_list, lower_snr_bound, noises, margin):
    for i in range(len(nodes_list)):
        tx_power = margin*noises[i]/nodes_list[i].pathloss_to_bs
        nodes_list[i].set_tx_power(tx_power)    

def allocate_resources_sequential(nodes, n_resources, **kwargs):
    #TODO: garantir que exista apenas 1 tabela de alocação no simulador inteiro??
    allocation_table = dict()
    n_nodes = len(nodes)
    for i in range(n_resources):
        if i>=n_nodes:
            break
        allocation_table[f'RB:{i}'] = nodes[i].id        
    return allocation_table

# def allocate_resources_d2d_sequential(dues_txs,d2d_access_table,mues,mues_allocation_table):
#     dues_allocation_table = dict()
#     for node in dues_txs:
#         pass

def d2d_critical_distance(margin, channel, lower_d2d_sinr, lower_mue_sinr, uc_bs_distance, uc_d2d_distance, d2d_bs_distance):
    d_crit = d2d_bs_distance*(margin/(lower_d2d_sinr*(1+lower_mue_sinr*(uc_bs_distance/uc_d2d_distance)**channel.prop_coef*(1+margin))))**(1/channel.prop_coef)
    return d_crit
    