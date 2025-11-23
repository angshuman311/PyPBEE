# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 18:34:41 2020

@author: Angshuman Deb
"""

import os
from pypbee.utility import Utility
from pypbee.structure import OSB
from pypbee.structural_analysis_platform import OpenSeesPy, OpenSeesTcl
from pypbee.prelim_analysis import PrelimAnalysis
from pypbee.avg_sa import AvgSa
from pypbee.sa_t import SaT
from pypbee.psha import PSHA
from pypbee.gms import GMS
from pypbee.edp import MaxColRebarStrain, MaxSpringDeformation
from pypbee.nltha import NLTHA
from pypbee.psdemha import PSDemHA
from pypbee.ds import DS
from pypbee.psdamha import PSDamHA
from pypbee.interp_exterp_model import InterpExterpModel
from pypbee.pygmm_extension.boore_atkinson_2008 import BooreAtkinson2008
import pygmm.baker_jayaram_2008
from scipy.stats import lognorm
import numpy as np

haz_lev_list = ['1', '2', '3', '4', '5', '6']
mrp_list = [72, 224, 475, 975, 2475, 4975]

name = 'Bridge_Column'
analysis_case = '100'
n_gm_list = [50] * len(mrp_list)
rng_seed = 'unique 6'  # for other Bridges

base_dir_path = r"E:\PyPBEE"
base_work_dir_path = rf"E:\PEER_UQ_Project_Work_Dir"
model_files_path = os.path.join(base_dir_path, 'scripts', 'examples', name)
model_work_dir_path = os.path.join(base_work_dir_path, f"{name}_Work_Dir")
local_opensees_path = os.path.join(base_dir_path, "vendor", "OpenSees_Windows", "OpenSees")
local_python_path = os.path.join(base_dir_path, "venv", "pypbee", "Scripts", "python.exe")
gm_database_dir_path = r"E:\NGA_PEER_EQ_DB"
local_bash_path = r"C:\Program Files\Git\git-bash"
model_params = Utility.import_attr_from_module(model_files_path, f"model_info", 'model_params')
location_info = Utility.import_attr_from_module(model_files_path, f"model_info", 'location_info')
structural_analysis_platform = OpenSeesPy(model_files_path, local_python_path)
osb = OSB(name, location_info, model_files_path, model_work_dir_path, model_params, structural_analysis_platform)
design_num_list = osb.get_design_num_list('all')

if analysis_case == '400':
    if name == 'Bridge_A':
        dp = [6.5, 0.015]  # Bridge_A
    if name == 'Bridge_B':
        dp = [6, 0.015]  # Bridge_B
    if name == 'Bridge_C':
        dp = [5.5, 0.01]  # Bridge_C
    if name == 'Bridge_MAOC':
        dp = [5, 0.02]  # Bridge_MAOC
    osb.add_design_point(dp)
    design_num_list = osb.get_design_num_list(osb.get_design_pts(['1']) + [dp])
    # design_num_list = osb.get_design_num_list(osb.get_design_pts(['1']))
    # design_num_list = osb.get_design_num_list([dp])
elif analysis_case == '100':
    design_num_list = osb.get_design_num_list('all')
else:
    design_num_list = ['1']

design_num_list = ['1']
prelim_analysis = PrelimAnalysis(osb, 8)
gmm = BooreAtkinson2008
correl_func = pygmm.baker_jayaram_2008.calc_correls
# im = SaT(osb, gmm, correl_func, 1.)
# im = SaT(osb, gmm, correl_func, 'T_1_trans')
im = AvgSa(osb, gmm, correl_func, ['T_1_trans', ], range_multiplier=[1, 2.5])
im_input = np.logspace(np.log10(1e-4), np.log10(5), 1000, endpoint=True)
# im_input = np.linspace(1e-4, 5, 1000, endpoint=True)
spectral_periods = np.logspace(np.log10(0.05), np.log10(5), 50, endpoint=True)
# spectral_periods = np.linspace(0.05, 5, 50, endpoint=True)
psha = PSHA(im)
gms = GMS(im)
gms_n_loop = 0
edp_list = [
    MaxColRebarStrain('compression', osb, '1', 'shared'),
    MaxColRebarStrain('tension', osb, '2', 'shared'),
    MaxColRebarStrain('range', osb, '3', 'shared'),
    MaxSpringDeformation('shear_key', 'compression', osb, '4', 'separate', normalize_with='$D3', small_value=1.0e-10, large_value=10.0)
]
delta_input_list = [np.logspace(np.log10(1e-7), np.log10(1.0), 1000)] * 3 + [np.logspace(np.log10(1e-7), np.log10(10.0), 1000)]
haz_req = dict()
haz_req['fit_dist'] = lognorm
haz_req['transf_on_dist_params_list'] = [lambda x, y, z: x, lambda x, y, z: y, lambda x, y, z: z]
haz_req['transf_to_dist_params_list'] = [lambda x, y, z: x, lambda x, y, z: y, lambda x, y, z: z]
haz_req['interp_exterp_function_list'] = [InterpExterpModel.piecewise_linear_interp_constant_extrap,
                                          lambda x, y, xq: 0 * xq,
                                          InterpExterpModel.power_law_regression_lin
                                          ]
[edp.set_haz_req(haz_req) for edp in edp_list]
nltha = NLTHA(edp_list, im)
nltha_setup_pool_size = 24
psdemha = PSDemHA(edp_list, im)
ds_list = list()
ds_list.append(
    DS(
        edp_list[0], lambda x: 0.004,
        haz_req={'normalized_fragility_dist': lognorm(0.326, 0, 1.02),
                 'estimation_sample_size': 5},
        ds_type='col_rebar_strain_damage'
    )
)
ds_list.append(
    DS(
        edp_list[1], lambda x: 0.03 + 700 * x[1] * x[2] / x[3] - 0.1 * x[8] / (x[4] * x[5]),
        haz_req={'normalized_fragility_dist': lognorm(0.201, 0, 1.05),
                 'estimation_sample_size': 5},
        ds_type='col_rebar_strain_damage'
    )
)
ds_list.append(
    DS(
        edp_list[2], lambda x: 0.11 + min(0.054, 3.2 * x[1]) - 0.0175 * abs(x[6] ** (1 / 3) - 2.93) - 0.054 * x[7],
        haz_req={'normalized_fragility_dist': lognorm(0.109, 0, 0.99),
                 'estimation_sample_size': 5},
        ds_type='col_rebar_strain_damage'
    )
)
ds_list.append(
    DS(
        edp_list[3], lambda x: 1.0,
        haz_req={'normalized_fragility_dist': lognorm(0.11, 0, 1.14),
                 'estimation_sample_size': 5},
        ds_type='spring_deformation_damage',
    )
)
psdamha = PSDamHA(ds_list, im=im, sol_type='numerical')

# f, a, _ = im.plot_gm_psa_spectra(
#     ['100', '1', '1', '1', '3'],
#     figkwargs={'num': 1, 'figsize': (3.56, 2.68)},
#     minor_grid_alpha=0.5,
#     plot_uhs=True,
# )
# f.savefig(r"E:\Dropbox\Personal\Journal Papers\Paper 5\gms.png", dpi=600)
#
# f, a, _ = im.plot_shc(
#     ['100', '1', '1', '1', '3'],
#     figkwargs={'num': 1, 'figsize': (3.56, 2.68)},
#     minor_grid_alpha=0.5,
#     plot_all=True,
# )
# f.savefig(r"E:\Dropbox\Personal\Journal Papers\Paper 5\shc.png", dpi=600)
#
# f, a, _ = edp_list[0].plot_conditional_demand_model_3d(
#     for_which=['100', '1', '1', '1'],
#     edp_str="col_1_edge_2",
#     haz_lev_list=haz_lev_list,
#     im=im,
#     im_lim=[0, 1],
#     edp_lim=[0.0001, 0.05],
#     figkwargs={'num': 1, 'figsize': (4.5, 3), 'constrained_layout': True},
# )
# f.savefig(r"E:\Dropbox\Personal\Journal Papers\Paper 5\cdm.png", dpi=600, transparent=True)
#
# f, a, _ = ds_list[0].plot_haz_deagg(
#     for_which=['100', '1', '1', '1'],
#     wrt="im",
# )
# f.savefig(r"E:\Dropbox\Personal\Journal Papers\Paper 5\haz_deagg.png", dpi=600)
#
# f, a, _ = osb.plot_scenario_rates(figkwargs={'num': 1, 'figsize': (4.5, 3), 'constrained_layout': True})
# f.savefig(r"E:\Dropbox\Personal\Journal Papers\Paper 5\rates.png", dpi=600, transparent=True)
#
# f, a, _ = ds_list[0].plot_mrp_histogram(
#     for_which=['200', '1', '1', '1'],
#     mrp_det=ds_list[0].get_damage_mrp_system_mean_estimate(['100', '1', '1', '1'])[0],
#     figkwargs={'num': 1, 'figsize': (4.5, 3), 'constrained_layout': True},
# )
# f.savefig(r"E:\Dropbox\Personal\Journal Papers\Paper 5\hist.png", dpi=600, transparent=True)