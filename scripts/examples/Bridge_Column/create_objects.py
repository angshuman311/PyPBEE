# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 18:34:41 2020

@author: Angshuman Deb
"""

import os
from pypbee.utility import Utility
from pypbee.structure import OSB
from pypbee.structural_analysis_platform import OpenSeesPy
from pypbee.prelim_analysis import PrelimAnalysis
from pypbee.sa_t import SaT
from pypbee.psha import PSHA
from pypbee.gms import GMS
from pypbee.edp import MaxColRebarStrain
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
n_gm_list = [10] * len(mrp_list)
rng_seed = 'unique 6'  # for other Bridges

base_dir_path = r"C:\Users\adeb\Work\PyPBEE"
base_work_dir_path = r"C:\Users\adeb\Work\PyPBEE_Work_Dir"
model_files_path = os.path.join(base_dir_path, 'scripts', 'examples', name, 'Model_Files')
model_work_dir_path = os.path.join(base_work_dir_path, f"{name}_Work_Dir")
local_opensees_path = os.path.join(base_dir_path, "vendor", "OpenSees_Windows", "OpenSees")
local_python_path = os.path.join(base_dir_path, "venv", "pypbee", "Scripts", "python.exe")
gm_database_dir_path = r"C:\Users\adeb\Work\NGA_PEER_EQ_DB"
local_bash_path = r"C:\Program Files\Git\git-bash"
model_params = Utility.import_attr_from_module(
    model_files_path, f"model_info", 'model_params'
)
location_info = Utility.import_attr_from_module(
    model_files_path, f"model_info", 'location_info'
)
structural_analysis_platform = OpenSeesPy(model_files_path, local_python_path)
osb = OSB(
    name,
    location_info,
    model_files_path,
    model_work_dir_path,
    model_params,
    structural_analysis_platform
)
gmm = BooreAtkinson2008
correl_func = pygmm.baker_jayaram_2008.calc_correls
im = SaT(osb, gmm, correl_func, 'T_1')
edp_list = [
    MaxColRebarStrain(
        max_what='compression',
        frame_structure=osb,
        tag='1', recorder_file_storage='shared'
    ),
    MaxColRebarStrain(
        max_what='tension',
        frame_structure=osb,
        tag='2',
        recorder_file_storage='shared'
    ),
]
haz_req = dict()
haz_req['fit_dist'] = lognorm
haz_req['transf_on_dist_params_list'] = [
    lambda x, y, z: x,
    lambda x, y, z: y,
    lambda x, y, z: z
]
haz_req['transf_to_dist_params_list'] = [
    lambda x, y, z: x,
    lambda x, y, z: y,
    lambda x, y, z: z
]
haz_req['interp_exterp_function_list'] = [
    InterpExterpModel.piecewise_linear_interp_constant_extrap,
    lambda x, y, xq: 0 * xq,
    InterpExterpModel.power_law_regression_lin
]
[edp.set_haz_req(haz_req) for edp in edp_list]

ds_list = list()
ds_list.append(
    DS(
        edp_list[0],
        lambda x: 0.004,
        haz_req={
            'normalized_fragility_dist': lognorm(0.326, 0, 1.02),
            'estimation_sample_size': 5
        },
        ds_type='col_rebar_strain_damage'
    )
)
ds_list.append(
    DS(
        edp_list[1],
        lambda x: 0.03 + 700 * x[1] * x[2] / x[3] - 0.1 * x[8] / (x[4] * x[5]),
        haz_req={
            'normalized_fragility_dist': lognorm(0.201, 0, 1.05),
            'estimation_sample_size': 5
        },
        ds_type='col_rebar_strain_damage'
    )
)

design_num_list = osb.get_design_num_list([[5.51, 0.02], [5, 0.01]])
num_modes = 4
im_input = np.logspace(np.log10(1e-4), np.log10(5), 1000, endpoint=True)
spectral_periods = np.logspace(np.log10(0.05), np.log10(5), 50, endpoint=True)
gms_n_loop = 0
nltha_setup_pool_size = 24
delta_input_list = [np.logspace(np.log10(1e-7), np.log10(1.0), 1000)] * 3 + [np.logspace(np.log10(1e-7), np.log10(10.0), 1000)]

prelim_analysis = PrelimAnalysis(osb, num_modes)
psha = PSHA(im)
gms = GMS(im)
nltha = NLTHA(edp_list, im)
psdemha = PSDemHA(edp_list, im)
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