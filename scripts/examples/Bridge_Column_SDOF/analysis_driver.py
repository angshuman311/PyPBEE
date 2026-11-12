# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 2026

@author: Angshuman Deb
"""
import numpy as np

from create_objects import analysis_case, rng_seed, \
    gm_database_dir_path, \
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

    edp1 = nltha.edp_list[0]
    fig4, ax4, _ = edp1.plot_conditional_demand_model_3d(
        for_which=['100', '1', '1', '1'],
        edp_str='node_2',
        haz_lev_list=haz_lev_list,
        im=im,
        im_lim=(1e-8, 1.),
        edp_lim=(1e-8, 15.)
    )
    ax4.set_xlabel('EDP_1 [in]')
    ax4.set_ylabel('Sa(T1) [g]')
    ax4.set_zlabel('P[EDP | IM]')

    fig5, axs5, _ = edp1.plot_conditional_demand_regression_model(
        for_which=['100', '1', '1', '1'],
        edp_str='node_2',
        haz_lev_list=haz_lev_list,
        im=im,
        im_lim=[0, 1],
    )
    axs5[0].set_xlabel('Sa [g]')
    axs5[1].set_xlabel('Sa [g]')
    axs5[2].set_xlabel('Sa [g]')
    axs5[0].set_ylabel('lognorm shape [-]')
    axs5[1].set_ylabel('lognorm loc [-]')
    axs5[2].set_ylabel('lognorm scale [in]')
    fig5.set_constrained_layout(True)

    fig6, ax6, _ = edp1.plot_dhc(
        for_which=['100', '1', '1', '1'],
        edp_str='node_2'
    )
    ax6.set_xlabel('EDP_1 [in]')
    ax6.set_ylabel('MAF [-]')

    ds1 = psdamha.ds_list[0]
    fig7, ax7, _ = ds1.plot_denormalized_fragility(
        for_which=['100', '1', '1', '1'],
        ds_str='node_2',
    )
    ax7.set_xlabel('EDP_1 [in]')
    ax7.set_ylabel('P[DS | EDP]')
    mrp_ds1 = 1 / ds1.get_psdamha_results(
        for_which=['100', '1', '1', '1'],
        ds_str='node_2'
    )[0]

    fig8, axs8, _ = ds1.plot_haz_deagg(
        for_which=['100', '1', '1', '1'], wrt='im'
    )
    axs8[0].set_xlabel('IM: SaT [g]')

    axs8[0].set_ylabel('MAF IM [-]', color='blue')
    axs8[0].tick_params(axis='y', colors='blue')

    axs8[1].set_ylabel('P[IM | LSE] [-]', color='red')
    axs8[1].tick_params(axis='y', colors='red')

    # one horizontal component of each record : 1 or 2
    gm_records = im.get_gms_results(for_which=['100', '1', '1', '1', '3'])[0]
    for rsn, comp, file_name in zip(gm_records['NGA_RSN'], gm_records['component'], gm_records['file_name_1']):
        print(f"RSN {rsn}, component {comp}, {file_name}")
