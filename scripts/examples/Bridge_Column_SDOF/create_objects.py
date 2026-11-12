# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 2026

@author: Angshuman Deb
"""

import os
from pypbee.utility import Utility
from pypbee.structure import Structure
from pypbee.structural_analysis_platform import OpenSeesPy
from pypbee.prelim_analysis import PrelimAnalysis
from pypbee.sa_t import SaT
from pypbee.sa import sigma_c_watson_lamprey_boore_2007
from pypbee.psha import PSHA
from pypbee.gms import GMS
from pypbee.nltha import NLTHA
from pypbee.psdemha import PSDemHA
from pypbee.ds import DS
from pypbee.psdamha import PSDamHA
from pypbee.interp_exterp_model import InterpExterpModel
from pypbee.pygmm_extension.boore_atkinson_2008 import BooreAtkinson2008
import pygmm.baker_jayaram_2008
from scipy.stats import lognorm
from max_disp import MaxDisp

haz_lev_list = ['1', '2', '3', '4', '5', '6']
mrp_list = [72, 224, 475, 975, 2475, 4975]

name = 'Bridge_Column_SDOF'
analysis_case = '100'
n_gm_list = [10] * len(mrp_list)
rng_seed = 'unique 6'

base_dir_path = r"C:\Users\adeb\Work\PyPBEE"
base_work_dir_path = os.path.join(base_dir_path, "PyPBEE_Work_Dir")
model_files_path = os.path.join(base_dir_path, 'scripts', 'examples', name, 'Model_Files')
model_work_dir_path = os.path.join(base_work_dir_path, f"{name}_Work_Dir")
local_python_path = os.path.join(base_dir_path, "venv", "pypbee", "Scripts", "python.exe")
gm_database_dir_path = os.path.join(base_dir_path, "data", "NGA_PEER_EQ_DB")
model_params = Utility.import_attr_from_module(
    model_files_path, f"model_info", 'model_params'
)
location_info = Utility.import_attr_from_module(
    model_files_path, f"model_info", 'location_info'
)
structural_analysis_platform = OpenSeesPy(model_files_path, local_python_path)
sdof = Structure(
    name,
    location_info,
    model_files_path,
    model_work_dir_path,
    model_params,
    structural_analysis_platform
)
gmm = BooreAtkinson2008
correl_func = pygmm.baker_jayaram_2008.calc_correls
im = SaT(
    sdof, gmm, correl_func, 'T_1',
    gm_component_mode='single_horizontal',
    sigma_component_func=sigma_c_watson_lamprey_boore_2007
)
edp_list = [
    MaxDisp(
        node_tag=2, dof=1,
        structure=sdof,
        tag='1', recorder_file_storage='shared'
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
        lambda x: 5.5,
        haz_req={
            'normalized_fragility_dist': lognorm(0.082, 0, 1.03),
            'estimation_sample_size': 5
        },
        ds_type='disp_damage'
    )
)

design_num_list = sdof.get_design_num_list([[363., 465.]])
num_modes = 1

prelim_analysis = PrelimAnalysis(sdof, num_modes)
psha = PSHA(im)
gms = GMS(im)
nltha = NLTHA(edp_list, im)
psdemha = PSDemHA(edp_list, im)
psdamha = PSDamHA(ds_list, im=im, sol_type='numerical')
