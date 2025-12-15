# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 14:53:16 2019

@author: Angshuman Deb
"""
import numpy as np

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

    nltha.setup(analysis_case, pool_size, design_num_list, haz_lev_list, gm_database_dir_path)
    nltha.run(analysis_case, pool_size)
    nltha.wrap_up(analysis_case)
    # 2 analysis didn't converge, changed algorithm and tolerances
    # nltha.run(analysis_case, pool_size)
    # nltha.wrap_up(analysis_case)
    # After all analyses converges:
    nltha.collect_edps(analysis_case, pool_size)

    psdemha.setup(analysis_case, design_num_list, haz_lev_list)
    psdemha.run(analysis_case, pool_size)
    psdemha.wrap_up(analysis_case)

    psdamha.setup(analysis_case, design_num_list, haz_lev_list, n_gm_list, rng_seed=rng_seed)
    psdamha.run(analysis_case, pool_size)
    psdamha.wrap_up(analysis_case)

    # Plotting results
    im = psha.im
    fig1, ax1, shc_data = im.plot_shc(for_which=['100', '1', '1', '1'])
    ax1.set_xlabel('Sa [g]')
    ax1.set_ylabel('MAF [-]')
    mrp = 475
    mrp_ind = np.argmin(np.abs(shc_data['shc'][:, 1] - 1/mrp))
    sa_at_mrp = shc_data['shc'][mrp_ind][0]
    ax1_xlim = ax1.get_xlim()
    ax1_ylim = ax1.get_ylim()
    ax1.plot(
        [ax1_xlim[0], sa_at_mrp, sa_at_mrp],
        [1/mrp, 1/mrp, ax1_ylim[0]],
        'r-', lw=0.75
    )
    ax1.text(ax1_xlim[0], 1/mrp, f"MRP = {mrp} yrs")
    ax1.text(sa_at_mrp, ax1_ylim[0], f"Sa = {sa_at_mrp:0.3f} g")

    fig2, ax2, _ = im.plot_seismic_hazard_deagg(
        mrp=475, for_which=['100', '1', '1', '1']
    )
    ax2.set_xlabel('Magnitude')
    ax2.set_ylabel('Source-to-site\ndistance [km]')
    ax2.set_zlabel('Contribution\nto hazard [-]')

    fig3, ax3, _ = im.plot_gm_psa_spectra(
        for_which=['100', '1', '1', '1', '3'], plot_uhs=True
    )
    ax3.set_xlabel('T [s]')
    ax3.set_ylabel('Sa [g]')

    edp2 = nltha.edp_list[1]
    fig4, ax4, _ = edp2.plot_conditional_demand_model_3d(
        for_which=['100', '1', '1', '1'],
        edp_str='col_1_edge_1',
        haz_lev_list=haz_lev_list,
        im=im,
        im_lim=(1e-8, 1.),
        edp_lim=(1e-8, 0.2)
    )
    ax4.set_xlabel('EDP_2 [-]')
    ax4.set_ylabel('Sa(T1) [g]')
    ax4.set_zlabel('P[EDP | IM]')

    fig5, axs5, _ = edp2.plot_conditional_demand_regression_model(
        for_which=['100', '1', '1', '1'],
        edp_str='col_1_edge_1',
        haz_lev_list=haz_lev_list,
        im=im,
        im_lim=[0, 1],
    )
    axs5[0].set_xlabel('Sa [g]')
    axs5[1].set_xlabel('Sa [g]')
    axs5[2].set_xlabel('Sa [g]')
    axs5[0].set_ylabel('lognorm shape [-]')
    axs5[1].set_ylabel('lognorm loc [-]')
    axs5[2].set_ylabel('lognorm scale [-]')
    fig5.set_constrained_layout(True)

    fig6, ax6, _ = edp2.plot_dhc(
        for_which=['100', '1', '1', '1'],
        edp_str='col_1_edge_1'
    )
    ax6.set_xlabel('EDP_2 [-]')
    ax6.set_ylabel('MAF [-]')

    ds2 = psdamha.ds_list[1]
    mrp_ds2 = 1 / ds2.get_psdamha_results(
        for_which=['100', '1', '1', '1'],
        ds_str='col_1_edge_1'
    )[0]

    fig7, axs7, _ = ds2.plot_haz_deagg(
        for_which=['100', '1', '1', '1'], wrt='im'
    )
    axs7[0].set_xlabel('IM: SaT [g]')

    axs7[0].set_ylabel('MAF IM [-]', color='blue')
    axs7[0].tick_params(axis='y', colors='blue')

    axs7[1].set_ylabel('P[IM | LSE] [-]', color='red')
    axs7[1].tick_params(axis='y', colors='red')

    nltha.extract_single_run(
        for_which=['100', '1', '1', '1', '3', '7'],
        num_modes=prelim_analysis.num_modes,
        gm_database_dir_path=gm_database_dir_path
    )