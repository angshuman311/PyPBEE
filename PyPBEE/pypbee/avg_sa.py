# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 16:16:25 2019

@author: Angshuman Deb
"""

import numpy as np
from .utility import Utility
from .sa import Sa
import os
import sys
import scipy.io as sio
import matplotlib.pyplot as plt
from matplotlib import colors
import mpl_toolkits.mplot3d.art3d as art3d


class AvgSa(Sa):
    large_err = 1.0e99

    ####################################################################################################################
    # Constructor
    ####################################################################################################################

    def __init__(self, structure, gmm, correl_func, define_range, **kwargs):
        self.define_range = define_range
        self.range_multiplier = kwargs.get('range_multiplier', [1, 1])
        self.num_periods = kwargs.get('num_periods', 10)
        self.type_spacing = kwargs.get('type_spacing', 'log')
        super().__init__(structure, gmm, correl_func)

    ####################################################################################################################
    # Static methods
    ####################################################################################################################

    @staticmethod
    def evaluate_gmm_helper(med_sa_period_range, sig_ln_sa_period_range, moc):
        n_p = med_sa_period_range.shape[1]

        med_sa_avg = np.exp((1 / n_p) * np.sum(np.log(med_sa_period_range), axis=1))
        sig_ij = np.einsum('ij,il->jli', sig_ln_sa_period_range, sig_ln_sa_period_range)
        rho_ij_sig_ij = np.einsum('ijk,ij->ijk', sig_ij, moc)
        sig_ln_sa_avg = np.sqrt(((1 / n_p) ** 2) * np.sum(rho_ij_sig_ij, axis=(0, 1)))

        # return median and std
        return med_sa_avg, sig_ln_sa_avg

    ####################################################################################################################
    # Public Functionalities (instance methods)
    ####################################################################################################################

    def evaluate_period(self, for_which):
        structure = self.structure
        define_range = self.define_range
        range_multiplier = self.range_multiplier
        num_periods = self.num_periods
        type_spacing = self.type_spacing
        period_range = structure.get_period_range(for_which, define_range, range_multiplier, num_periods, type_spacing)
        return period_range

    # ------------------------------------------------------------------------------------------------------------------

    def setup_psa_spectra_gms(self, mrp, for_which, mean_req_exact, sig_req_exact, spectral_periods,
                              gm, sa_known, per_known, is_scaled):
        # Here zeros are ok as intialized values of indices because they are going to be filled up before they are used
        rec_per = np.zeros(len(spectral_periods)).astype(int)
        for i in range(len(spectral_periods)):
            rec_per[i] = np.argmin(np.abs(per_known - spectral_periods[i]))

        sa_known[sa_known == 0] = sys.float_info.epsilon
        sample_big = np.log(sa_known[:, rec_per])

        period_range = self.evaluate_period(for_which)
        target = self.get_inv_seismic_hazard(mrp, for_which)
        ind = np.in1d(spectral_periods, period_range)
        rec_avg = np.exp(np.sum(sample_big[:, ind], axis=1) / len(period_range))
        if is_scaled:
            scale_fac = target / rec_avg
        else:
            scale_fac = rec_avg / rec_avg

        return sample_big, scale_fac, mean_req_exact, sig_req_exact, gm

    # ------------------------------------------------------------------------------------------------------------------

    def evaluate_gmm(self, period_range):
        med_sa_period_range, sig_ln_sa_period_range = self._evaluate_gmm(period_range)
        moc = self.get_moc(period_range)
        med_sa_avg, sig_ln_sa_avg = AvgSa.evaluate_gmm_helper(med_sa_period_range, sig_ln_sa_period_range, moc)
        # return 1-d arrays
        return med_sa_avg, sig_ln_sa_avg

    # ------------------------------------------------------------------------------------------------------------------

    def _get_target_spectrum(self, mrp, for_which, med_sa_spectrum, sig_ln_sa_spectrum,
                             spectral_periods, uncondition):
        if uncondition:
            scenario_weights = self.get_seismic_hazard_deagg(mrp, for_which)
            mean_ln_sa_spectrum, sig_ln_sa_spectrum = Utility.get_exact_mixture_params(
                np.log(med_sa_spectrum),
                sig_ln_sa_spectrum,
                scenario_weights)  # (1, n_sp), (1, n_sp)
            med_sa_spectrum = np.exp(mean_ln_sa_spectrum)  # (1, n_sp)

        period_range = self.evaluate_period(for_which)

        n_s = med_sa_spectrum.shape[0]
        n_sp = len(spectral_periods)
        n_p = len(period_range)

        ind_period_range = np.in1d(spectral_periods, period_range)
        med_sa_period_range = med_sa_spectrum[:, ind_period_range]  # (n_s, n_p)
        sig_ln_sa_period_range = sig_ln_sa_spectrum[:, ind_period_range]  # (n_s, n_p)

        moc = self.get_moc(period_range)
        med_sa_avg, sig_ln_sa_avg = AvgSa.evaluate_gmm_helper(med_sa_period_range,
                                                              sig_ln_sa_period_range,
                                                              moc)  # (n_s,), (n_s,)

        rho_ij = np.zeros((n_sp, n_p))  # (n_sp, n_p)
        for j in range(n_p):
            rho_ij[:, j] = self.correl_func(spectral_periods, period_range[j])

        rho = (rho_ij @ sig_ln_sa_period_range.T).T  # ((n_sp, n_p) @ (n_p, n_s)).T = (n_s, n_sp)
        rho = rho / (n_p * sig_ln_sa_avg.reshape(n_s, 1))  # (n_s, n_sp) / (n_s, 1) = (n_s, n_sp)

        mean_req, cov_req = self.get_conditional_spectrum(mrp, for_which,
                                                          med_sa_spectrum,
                                                          sig_ln_sa_spectrum,
                                                          rho, self.get_moc(spectral_periods),
                                                          med_sa_avg, sig_ln_sa_avg)
        return mean_req, cov_req, spectral_periods

    # ------------------------------------------------------------------------------------------------------------------

    def plot_shc(self, for_which, **kwargs):
        save_dir_path = kwargs.get('save_dir_path', Utility.get_path(os.getcwd()))
        save_mat = kwargs.get('save_mat', False)
        lc = kwargs.get('lc', 'black')
        lw = kwargs.get('lw', 1.5)
        plot_all = kwargs.get('plot_all', False)

        fig, ax, export_mat_dict = super().plot_shc(for_which, **kwargs)
        structure = self.structure
        if plot_all:
            period_range = self.evaluate_period(for_which)
            shc = export_mat_dict['shc']
            shc_list = self.shc_at_periods(period_range, im_input=shc[:, 0])
            for itr in range(len(period_range)):
                shc_individual = shc_list[itr]
                ax.loglog(shc_individual[:, 0], shc_individual[:, 1], '-',
                          color=0.5 * (np.ones(3) - np.array(colors.to_rgb(lc))), alpha=0.5, linewidth=lw / 2)
                export_mat_dict[f'shc_individual_{itr + 1}'] = shc_individual
            if save_mat:
                name = structure.name
                for_which_str = '_'.join(for_which)
                file_name = f'plot_shc_{name}_{for_which_str}.mat'
                sio.savemat(Utility.get_path(save_dir_path, file_name), export_mat_dict)

        return fig, ax, export_mat_dict

    # ------------------------------------------------------------------------------------------------------------------

    def mark_period_on_axis(self, for_which, ax, **kwargs):
        patch_alpha = kwargs.get('period_patch_alpha', 0.25)
        patch_color = kwargs.get('period_patch_color', (0, 1, 0))
        period_range = self.evaluate_period(for_which)
        ylim = Utility.get_ylim(ax, (0.01, None))
        rect = plt.Rectangle(
            (period_range[0], ylim[0]), period_range[-1] - period_range[0], ylim[1] - ylim[0],
            color=patch_color, alpha=patch_alpha
        )
        ax.add_patch(rect)
        if '3D' in ax.__class__.__name__:
            art3d.pathpatch_2d_to_3d(rect, z=0)
        ax.set_ylim(ylim)
