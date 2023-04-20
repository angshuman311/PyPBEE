# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 15:40:20 2020

@author: Angshuman Deb
"""

from pypbee.utility import Utility
import matplotlib.pyplot as plt
from create_objects_local import im, gm_database_dir_path, nltha, prelim_analysis

# ---------------------
# For local extract run
# ---------------------

# for_which = ['100', '1', '1', '1', '6', '1']
# for_which = ['400', '1', '1', '13', '6', '39']
# for_which = ['400', '1', '1', '11', '6', '1']
for_which_str = '100 1 1 1 6 1'
for_which = for_which_str.split()

element_response_list = list()

element_response_list.append([2001, 'localForce'])

element_response_list.append([1501, 'force'])
element_response_list.append([1601, 'force'])
element_response_list.append([1501, 'deformation'])
element_response_list.append([1601, 'deformation'])

element_section_response_list = list()
element_section_response_list.append([2001, 1, 'force'])
element_section_response_list.append([2001, 8, 'force'])
element_section_response_list.append([2001, 1, 'deformation'])
element_section_response_list.append([2001, 8, 'deformation'])


element_section_fiber_response_list = list()
element_section_fiber_response_list.append([2001, 1, 2, 0, 30, 1001])
element_section_fiber_response_list.append([2001, 1, 3, 0, 30, 2001])

node_response_list = list()
node_response_list.append([200, 1, 'disp'])
node_response_list.append([200, 2, 'disp'])
node_response_list.append([2002, 1, 'accel'])
node_response_list.append([2002, 2, 'accel'])

recorder_info = {'element_response_list': element_response_list,
                 'element_section_response_list': element_section_response_list,
                 'element_section_fiber_response_list': element_section_fiber_response_list,
                 'node_response_list': node_response_list,
                 'addnl_file_list': ["generate_curv_diag"]
                 }

nltha.extract_single_run(for_which, prelim_analysis.num_modes, gm_database_dir_path, recorder_info=recorder_info, ai_end=0.6)

plt.close('all')

# ----------------------------------------------------------------------------
# Input ground motion pseudo spectral acceleration spectrum
# ----------------------------------------------------------------------------
_, ax, _ = im.plot_gm_psa_spectra(for_which, figkwargs={'num': 1, 'figsize': (3.56, 2.68)}, save_mat=False,
                                  minor_grid_alpha=0.5)
ax.set_xlim(0.05, 1.5)
ax.set_ylim(ax.get_ylim()[0], 3)
# plt.savefig('gms.png')

# ----------------------------------------------------------------------------
# Input ground motion acceleration history
# ----------------------------------------------------------------------------
gm_file_names, scale_fac = im.get_gm_rec_info(for_which)
Utility.plot_raw_ground_acc_time_history(gm_database_dir_path, gm_file_names=gm_file_names, scale_fac=scale_fac)
im.plot_gm_psa_spectra(for_which)
