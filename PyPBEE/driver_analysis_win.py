# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 14:53:16 2019

@author: Angshuman Deb
"""

from create_objects_local import osb, analysis_case, rng_seed, \
    gm_database_dir_path, local_bash_path, \
    design_num_list, haz_lev_list, mrp_list, n_gm_list, im_input, spectral_periods, delta_input_list, \
    prelim_analysis, psha, gms, nltha, psdemha, psdamha, \
    Utility

# script code to be put within if __name__ == "__main__" block so that nothing from main gets rerun when subprocess
# from within imported modules are running !!!! to avoid iPython problems with multiprocessing on Windows machines,
# set run -> config -> execute in external system terminal

if __name__ == "__main__":

    pool_size = 10
    gms_n_loop = 2

    # prelim_analysis.setup(analysis_case, design_num_list, rng_seed=rng_seed)
    # prelim_analysis.run(analysis_case, pool_size)
    # prelim_analysis.wrap_up(analysis_case)

    psha.setup(analysis_case, design_num_list)
    psha.run(analysis_case, pool_size, im_input=im_input)
    psha.wrap_up(analysis_case)

    #
    gms.setup(analysis_case, design_num_list, haz_lev_list, mrp_list, n_gm_list)
    gms.run(analysis_case, pool_size, spectral_periods=spectral_periods, n_loop=gms_n_loop, rng_seed=rng_seed, uhs=True)
    gms.wrap_up(analysis_case)
    #
    # nltha.extract_single_run(['100', '1', '1', '1', '1', '1'], prelim_analysis.num_modes, gm_database_dir_path)
    # nltha.stampede_upload_pre(analysis_case, design_num_list, haz_lev_list, 'adeb', 'DesignSafe-Conte', 8, 50,
    #                           f'/scratch/04236/adeb/stampede2/PEER_UQ_Project/Software/{osb.name}',
    #                           local_bash_path=local_bash_path,
    #                           nltha_include=False, extract=True, stampede_comp_env='stampede_knl',
    #                           new_file=True, keep_file=False,
    #                           )

    # nltha.setup(analysis_case, pool_size, design_num_list, haz_lev_list, gm_database_dir_path, stage=True)
    # nltha.restage(analysis_case, pool_size, gm_database_dir_path, local_opensees_path=local_opensees_path)
    # nltha.run(analysis_case, pool_size, n_batch=1, batch_index=1, restage=False)
    # nltha.wrap_up(analysis_case)
    # a, _, _ = nltha.get_post_analysis_lists(analysis_case)
    # nltha.collect_edps(analysis_case, pool_size, a)

    # psdemha.setup(analysis_case, design_num_list, haz_lev_list)
    # psdemha.run(analysis_case, pool_size, delta_input_list=delta_input_list)
    # psdemha.wrap_up(analysis_case)

    # psdamha.setup(analysis_case, design_num_list, haz_lev_list=haz_lev_list, n_gm_list=n_gm_list, rng_seed=rng_seed)
    # psdamha.run(analysis_case, pool_size, delta_input_list=delta_input_list)
    # psdamha.wrap_up(analysis_case)
    #