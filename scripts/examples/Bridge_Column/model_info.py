# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 15:38:08 2019

@author: Angshuman Deb
"""

import numpy as np
from pypbee.utility import Utility
from scipy.stats import beta, lognorm, norm

########################################################################################################################
# Possible Random Model Parameters for Bridge Column
########################################################################################################################

fpc_e = 1.3 * 5  # Caltrans SDC v1.7
Ec_e = (40000. * np.sqrt(fpc_e * 1000) + 1e6) / 1000  # Caltrans SDC v1.7
fy_e = 69.144  # Darwin
fu_e = 95.1969  # Darwin
Es_e = 29200  # Lee Mosalam (2006)
b_e = 0.01  # calculated
w_e = 143.96 * 0.001 / (12 ** 3)  # Caltrans SDC v1.7
xi_e = 0.01  # Self

fpc_COV = 0.16  # Restrepo
Ec_COV = 0.10  # Restrepo
fy_COV = 0.04  # Carreno
fu_COV = fy_COV / 1.4  # Restrepo
Es_COV = 0.033  # Lee Mosalam (2006)
b_COV = 0.30  # calculatedf
w_COV = 0.04  # Restrepo, Weldon, Conte
xi_COV = 0.50  # Astroza

fy_min = 60  # Darwin
fy_max = 85.4  # Darwin

fu_min = 80  # Darwin
fu_max = 116  # Darwin

b_shift = 0.005
xi_shift = 0.003

corr_fpc_Ec = 0.40  # Restrepo
corr_fy_fu = 0.85  # Lee Mosalam (2006)
corr_fy_b = 0.40  # calculated
corr_fu_b = 0.50  # calculated
corr_Es_b = -0.13  # calculated

a_fy, b_fy = Utility.get_standard_beta_dist_params(fy_e, fy_COV * abs(fy_e), fy_min, fy_max)
a_fu, b_fu = Utility.get_standard_beta_dist_params(fu_e, fu_COV * abs(fu_e), fu_min, fu_max)

num_bar_clusters = 4
num_cols_total = 1
num_secdef_per_col = 3  # two extremes and one in middle
col_edge_list = [1]  # for EDP definition

# Primary design parameters
primary_design_vars = {'1': [5.51, 0.02]}

col_dia = [5, 6, 7, 8]
rho_long = [0.01, 0.015, 0.02, 0.025, 0.03]

sample_size = 50
dist_params_sample_size = 10000
estimation_sample_size = 5

sampling_method = 'lhs'

########################################################################################################################
# Random model parameters
########################################################################################################################

prob_dist_list = list()

# prob_dists
# 0: fpc
# 1: Ec
# 2: fy
# 3: fu
# 4: Es
# 5: b
# 6: w
# 7: xi

prob_dist_list.append(norm(fpc_e, fpc_COV * abs(fpc_e)))
prob_dist_list.append(norm(Ec_e, Ec_COV * abs(Ec_e)))
prob_dist_list.append(beta(a_fy, b_fy, fy_min, fy_max - fy_min))
prob_dist_list.append(beta(a_fu, b_fu, fu_min, fu_max - fu_min))
prob_dist_list.append(norm(Es_e, Es_COV * abs(Es_e)))
prob_dist_list.append(lognorm(*Utility.get_lognormal_dist_params(b_e, b_COV * abs(b_e), loc=b_shift)))
prob_dist_list.append(norm(w_e, w_COV * abs(w_e)))
prob_dist_list.append(lognorm(*Utility.get_lognormal_dist_params(xi_e, xi_COV * abs(xi_e), loc=xi_shift)))

corr_matrix = np.eye(len(prob_dist_list))
corr_matrix[0, 1] = corr_fpc_Ec
corr_matrix[1, 0] = corr_fpc_Ec

corr_matrix[2, 3] = corr_fy_fu
corr_matrix[3, 2] = corr_fy_fu
corr_matrix[2, 5] = corr_fy_b
corr_matrix[5, 2] = corr_fy_b
corr_matrix[3, 5] = corr_fu_b
corr_matrix[5, 3] = corr_fu_b
corr_matrix[4, 5] = corr_Es_b
corr_matrix[5, 4] = corr_Es_b

rv_list = list()
diag_blocks = [[0, 1], [2, 3, 4, 5], [6], [7]]

# ----------------------------------------------------------------------------------------------------------------------
# Concrete parameters
# ----------------------------------------------------------------------------------------------------------------------

expected_value_list = [fpc_e, Ec_e]
name_list = list()
for i_col in list(range(1, num_cols_total + 1)):
    for i_sec in list(range(1, num_secdef_per_col + 1)):
        name_list.append(f'fpc_col_{i_col}_secdef_{i_sec}')
        name_list.append(f'Ec_col_{i_col}_secdef_{i_sec}')
rv_list_entry = {'prob_dist_ind_list': [0, 1],
                 'count': num_cols_total * num_secdef_per_col,
                 'name_list': name_list,
                 'expected_value_list': expected_value_list
                 }
rv_list.append(rv_list_entry)

# ----------------------------------------------------------------------------------------------------------------------
# Steel parameters
# ----------------------------------------------------------------------------------------------------------------------

expected_value_list = [fy_e, fu_e, Es_e, b_e]
name_list = list()
for i_col in list(range(1, num_cols_total + 1)):
    for i_sec in list(range(1, num_secdef_per_col + 1)):
        for i_bc in list(range(1, num_bar_clusters + 1)):
            name_list.append(f'fy_col_{i_col}_secdef_{i_sec}_barcluster_{i_bc}')
            name_list.append(f'fu_col_{i_col}_secdef_{i_sec}_barcluster_{i_bc}')
            name_list.append(f'Es_col_{i_col}_secdef_{i_sec}_barcluster_{i_bc}')
            name_list.append(f'b_col_{i_col}_secdef_{i_sec}_barcluster_{i_bc}')

rv_list_entry = {'prob_dist_ind_list': [2, 3, 4, 5],
                 'count': num_cols_total * num_secdef_per_col * num_bar_clusters,
                 'name_list': name_list,
                 'expected_value_list': expected_value_list
                 }
rv_list.append(rv_list_entry)

# ----------------------------------------------------------------------------------------------------------------------
# Mass and damping
# ----------------------------------------------------------------------------------------------------------------------

expected_value_list = [w_e, xi_e]
name_list = ['wconc_all', 'damp_ratio']
rv_list_entry = {'prob_dist_ind_list': [6, 7],
                 'count': 1,
                 'name_list': name_list,
                 'expected_value_list': expected_value_list
                 }
rv_list.append(rv_list_entry)

# ----------------------------------------------------------------------------------------------------------------------

random_model_params = dict()
random_model_params['rv_list'] = rv_list
random_model_params['prob_dist_list'] = prob_dist_list
random_model_params['corr_matrix'] = corr_matrix
random_model_params['diag_blocks'] = diag_blocks
random_model_params['sample_size'] = sample_size
random_model_params['sampling_method'] = sampling_method
random_model_params['estimation_sample_size'] = estimation_sample_size
random_model_params['dist_params_sample_size'] = dist_params_sample_size

########################################################################################################################
# Primary design parameters
########################################################################################################################

design_num = 2
for i_col_dia in range(len(col_dia)):
    for i_rho_long in range(len(rho_long)):
        primary_design_vars.update({f'{design_num}': [col_dia[i_col_dia], rho_long[i_rho_long]]})
        design_num += 1

n_dp = len(primary_design_vars.keys())
primary_design_params = dict()
primary_design_params['value_list_dict'] = primary_design_vars
primary_design_params['name_list'] = ['all_col_dia_in_ft', 'all_rho_long']
# to define grid point vs non-grid-point
primary_design_params['design_point_qualifier_dict'] = dict(zip(list(primary_design_vars.keys()),
                                                                ['as-designed non-gp'] + ['gp'] * (n_dp - 1)))

########################################################################################################################
# Other model parameters
########################################################################################################################

other_model_params = dict()
other_model_params['value_list'] = [num_bar_clusters, num_secdef_per_col]
other_model_params['name_list'] = ['num_bar_clusters', 'num_secdef_per_col']

########################################################################################################################
# Misc parameters
########################################################################################################################

damping_models = list()

damping_model = dict()
damping_model['name'] = 'rayleigh_damping_first_trans_mode'
damping_model['xi_i'] = ['damp_ratio']
damping_models.append(damping_model)

damping_model = dict()
damping_model['name'] = 'rayleigh_damping'
damping_model['xi_i'] = ['damp_ratio'] * 2
damping_model['i'] = [1, 2]
damping_model['w_i'] = [None, None]
damping_models.append(damping_model)

damping_model = dict()
damping_model['name'] = 'modal_damping'
damping_model['xi_i'] = ['damp_ratio'] * 4 + list(np.linspace(0.01, 0.10, 10))
damping_models.append(damping_model)

model_attributes = dict()
model_attributes['vertical_gm'] = [0, 1]
model_attributes['ssi_springs'] = [0]
model_attributes['base_hinge'] = [0]
model_attributes['steel_material'] = ['SteelMPF']
model_attributes['col_elem_type'] = ['1', '2', '3', '4', '5']
model_attributes['damping'] = ['rayleigh_damping_first_trans_mode', 'rayleigh_damping', 'modal_damping']
model_attributes['backfill_resp_skew'] = [1, 0]

########################################################################################################################
# model_params
########################################################################################################################

model_params = dict()
model_params['primary_design_params'] = primary_design_params
model_params['random_model_params'] = random_model_params
model_params['other_model_params'] = other_model_params
model_params['damping_models'] = damping_models
model_params['model_attributes'] = model_attributes
model_params['num_cols_total'] = num_cols_total
model_params['col_edge_list'] = col_edge_list

########################################################################################################################
# location
########################################################################################################################

location_info = dict()
location_info['mechanism'] = 'U'
location_info['latitude'] = 37.7531
location_info['longitude'] = -121.1427
location_info['v_s30'] = 216.8470
location_info['region'] = 'california'
