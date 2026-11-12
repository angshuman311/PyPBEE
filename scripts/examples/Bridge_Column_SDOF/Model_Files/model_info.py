# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 2026

@author: Angshuman Deb
"""

import numpy as np
from pypbee.utility import Utility
from scipy.stats import lognorm, norm

########################################################################################################################
# Equivalent SDOF of Bridge_Column (as-designed column : dia = 5.51 ft, rho_long = 0.02)
########################################################################################################################

W_e = 3143.  # kip, deck + half column
k_e = 363.  # kip/in, gives T_1 of Bridge_Column (0.94 s)
Fy_e = 465.  # kip, pushover of Bridge_Column
b_e = 0.01  # same as Bridge_Column
xi_e = 0.01  # same as Bridge_Column

W_COV = 0.04  # same as w in Bridge_Column
b_COV = 0.30  # same as Bridge_Column
xi_COV = 0.50  # same as Bridge_Column

b_shift = 0.005
xi_shift = 0.003

sample_size = 50
dist_params_sample_size = 10000
estimation_sample_size = 5

sampling_method = 'lhs'

########################################################################################################################
# Random model parameters
########################################################################################################################

prob_dist_list = list()

# prob_dists
# 0: W
# 1: b
# 2: xi

prob_dist_list.append(norm(W_e, W_COV * abs(W_e)))
prob_dist_list.append(lognorm(*Utility.get_lognormal_dist_params(b_e, b_COV * abs(b_e), loc=b_shift)))
prob_dist_list.append(lognorm(*Utility.get_lognormal_dist_params(xi_e, xi_COV * abs(xi_e), loc=xi_shift)))

corr_matrix = np.eye(len(prob_dist_list))
diag_blocks = [[0], [1], [2]]

rv_list = list()
rv_list_entry = {
    'prob_dist_ind_list': [0, 1, 2],
    'count': 1,
    'name_list': ['W', 'b', 'damp_ratio'],
    'expected_value_list': [W_e, b_e, xi_e]
}
rv_list.append(rv_list_entry)

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

primary_design_vars = {'1': [k_e, Fy_e]}

primary_design_params = dict()
primary_design_params['value_list_dict'] = primary_design_vars
primary_design_params['name_list'] = ['k_in_kip_per_in', 'Fy_in_kip']
primary_design_params['design_point_qualifier_dict'] = {'1': 'as-designed non-gp'}

########################################################################################################################
# Other model parameters
########################################################################################################################

other_model_params = dict()
other_model_params['value_list'] = []
other_model_params['name_list'] = []

########################################################################################################################
# Misc parameters
########################################################################################################################

# damping is in ops_model.py (1 dof : modal damping of PyPBEE needs more than 1 dof)
damping_models = list()

model_attributes = dict()
model_attributes['spring_material'] = ['Steel01']

########################################################################################################################
# model_params
########################################################################################################################

model_params = dict()
model_params['primary_design_params'] = primary_design_params
model_params['random_model_params'] = random_model_params
model_params['other_model_params'] = other_model_params
model_params['damping_models'] = damping_models
model_params['model_attributes'] = model_attributes

########################################################################################################################
# location
########################################################################################################################

location_info = dict()
location_info['mechanism'] = 'U'
location_info['latitude'] = 37.7531
location_info['longitude'] = -121.1427
location_info['v_s30'] = 216.8470
location_info['region'] = 'california'
