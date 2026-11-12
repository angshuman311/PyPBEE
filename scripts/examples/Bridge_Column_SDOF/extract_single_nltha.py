# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 2026

@author: Angshuman Deb
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from create_objects import name, model_files_path, gm_database_dir_path, sdof, nltha, prelim_analysis

# ---------------------
# For local extract run
# ---------------------

for_which_str = '100 1 1 1 3 7'
for_which = for_which_str.split()

element_response_list = list()
element_response_list.append([1, 'deformation'])
element_response_list.append([1, 'force'])

recorder_info = {'element_response_list': element_response_list}

nltha.extract_single_run(for_which, prelim_analysis.num_modes, gm_database_dir_path, recorder_info=recorder_info)

extract_dir_path = os.path.join(os.path.dirname(model_files_path), f'{name}_extract_OpenSeesPy')
sdof.structural_analysis_platform.run_file(extract_dir_path, 'run_nltha.py', 'local')

# ----------------------------------------------------------------------------
# Spring hysteresis
# ----------------------------------------------------------------------------
rec_dir_path = os.path.join(extract_dir_path, 'NLTHA_Recorder_Files')
u = np.loadtxt(os.path.join(rec_dir_path, 'element_response_elem_1_resp_deformation.txt'))[:, 1]
f = np.loadtxt(os.path.join(rec_dir_path, 'element_response_elem_1_resp_force.txt'))[:, -1]

fig, ax = plt.subplots()
ax.plot(u, f, 'k-', lw=0.75)
ax.grid(True, alpha=0.5)
ax.set_xlabel('Spring deformation [in]')
ax.set_ylabel('Spring force [kip]')
ax.set_title(f'for_which = {for_which}')
