# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 14:53:16 2019

@author: Angshuman Deb
"""

from create_objects import analysis_case, rng_seed, \
    gm_database_dir_path, local_bash_path,\
    design_num_list, haz_lev_list, mrp_list, n_gm_list, \
    prelim_analysis, psha, gms, nltha, psdemha, psdamha

if __name__ == "__main__":

    pool_size = 12

    prelim_analysis.setup(analysis_case, design_num_list, rng_seed=rng_seed)
    prelim_analysis.run(analysis_case, pool_size)
    prelim_analysis.wrap_up(analysis_case)

    psha.setup(analysis_case, design_num_list)
    psha.run(analysis_case, pool_size)
    psha.wrap_up(analysis_case)

    gms.setup(analysis_case, design_num_list, haz_lev_list, mrp_list, n_gm_list)
    gms.run(analysis_case, pool_size, rng_seed=rng_seed, n_loop=2)
    gms.wrap_up(analysis_case)

    nltha.stampede_upload_pre(
        analysis_case, design_num_list, haz_lev_list, 'adeb', 'BCS25086', 8, 50,
        f'/scratch/04236/adeb/stampede3/PEER_UQ_Project/Software/{nltha.structure.name}',
        local_bash_path=local_bash_path,
        nltha_include=False, extract=True, stampede_comp_env='stampede_knl',
        new_file=True, keep_file=False,
    )

    nltha.setup(analysis_case, pool_size, design_num_list, haz_lev_list, gm_database_dir_path)
    nltha.run(analysis_case, pool_size)
    nltha.wrap_up(analysis_case)
    # # 2 analysis didn't converge, changed algorithm and tolerances
    # nltha.run(analysis_case, pool_size)
    # nltha.wrap_up(analysis_case)
    # # After all analyses converges:
    nltha.collect_edps(analysis_case, pool_size)

    psdemha.setup(analysis_case, design_num_list, haz_lev_list)
    psdemha.run(analysis_case, pool_size)
    psdemha.wrap_up(analysis_case)

    psdamha.setup(analysis_case, design_num_list, haz_lev_list, n_gm_list, rng_seed=rng_seed)
    psdamha.run(analysis_case, pool_size)
    psdamha.wrap_up(analysis_case)

    # Plotting results
    im = psha.im
    edp1 = nltha.edp_list[0]
    ds1 = psdamha.ds_list[0]
    im.plot_shc(for_which=['100', '1', '1', '1', '3'])
    im.plot_seismic_hazard_deagg(mrp=475, for_which=['100', '1', '1', '1', '3'])
    im.plot_gm_psa_spectra(for_which=['100', '1', '1', '1', '3'], plot_uhs=True)

    edp1.plot_conditional_demand_model_3d(
        for_which=['100', '1', '1', '1'],
        edp_str='col_1_edge_1',
        haz_lev_list=haz_lev_list,
        im=im,
        im_lim=(1e-8, 1.),
        edp_lim=(1e-8, 0.125)
    )

    nltha.extract_single_run(
        for_which=['100', '1', '1', '1', '3', '7'],
        num_modes=prelim_analysis.num_modes,
        gm_database_dir_path=gm_database_dir_path
    )

    # ds1.plot_mrp_histogram(for_which=['200', '1', '1'])