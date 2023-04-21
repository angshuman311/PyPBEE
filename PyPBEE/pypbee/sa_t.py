# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 16:16:25 2019

@author: Angshuman Deb
"""

import numpy as np
from .structure import Structure
from .sa import Sa
from .utility import Utility
import os
import sys
import scipy.io as sio


class SaT(Sa):
    large_value = 1e99
    small_value = 1e-99

    ####################################################################################################################
    # Constructor
    ####################################################################################################################

    def __init__(self, structure, gmm, correl_func, period):
        self.period = period
        super().__init__(structure, gmm, correl_func)

    ####################################################################################################################
    # Public Functionalities (instance methods)
    ####################################################################################################################

    def evaluate_period(self, for_which):
        period = self.period
        if isinstance(period, str) and period.startswith('T_'):
            structure = self.structure
            eval_period = structure.get_period_range_extremes(for_which, [period, ])[0]
        elif isinstance(period, float):
            eval_period = period
        else:
            return

        rec_data_file_path = Utility.get_path(os.path.dirname(__file__), 'Rec_Data', 'rec_data.pickle')
        rec_data = Utility.pickle_load_dict(rec_data_file_path)
        per_known = rec_data['per_known']
        return per_known[np.argmin(np.abs(per_known - eval_period))]

    # ------------------------------------------------------------------------------------------------------------------

    def _cond_spectrum(self, mrp, for_which, med_sa_spectrum, sig_ln_sa_spectrum, spectral_periods):
        target = self.get_inv_seismic_hazard(mrp, for_which)
        period = self.evaluate_period(for_which)
        period_mask_lt = spectral_periods < period
        period_mask_gt = spectral_periods > period

        med_sa_spectrum = np.hstack(
            (med_sa_spectrum[period_mask_lt],
             target,
             med_sa_spectrum[period_mask_gt])
        )

        sig_ln_sa_spectrum = np.hstack(
            (sig_ln_sa_spectrum[period_mask_lt],
             0.0,
             sig_ln_sa_spectrum[period_mask_gt])
        )
        spectral_periods = Utility.merge_sorted_1d_arrays(spectral_periods, [period])
        return med_sa_spectrum, sig_ln_sa_spectrum, spectral_periods

    # ------------------------------------------------------------------------------------------------------------------

    def setup_psa_spectra_gms(self, mrp, for_which, mean_req_exact, sig_req_exact, spectral_periods,
                              gm, sa_known, per_known, is_scaled):
        target = self.get_inv_seismic_hazard(mrp, for_which)
        period = self.evaluate_period(for_which)
        period_mask_lt = spectral_periods < period
        period_mask_gt = spectral_periods > period

        mean_req_exact = np.hstack(
            (mean_req_exact[period_mask_lt],
             np.log(target),
             mean_req_exact[period_mask_gt])
        )

        sig_req_exact = np.hstack(
            (sig_req_exact[period_mask_lt],
             0.0,
             sig_req_exact[period_mask_gt])
        )

        gm = np.column_stack(
            (gm[:, period_mask_lt],
             target * np.ones(gm.shape[0]),
             gm[:, period_mask_gt])
        )

        spectral_periods = Utility.merge_sorted_1d_arrays(spectral_periods, [period])
        # Here zeros are ok as intialized values of indices because they are going to be filled up before they are used
        rec_per = np.zeros(len(spectral_periods)).astype(int)
        for i in range(len(spectral_periods)):
            rec_per[i] = np.argmin(np.abs(per_known - spectral_periods[i]))

        sa_known[sa_known == 0] = sys.float_info.epsilon
        sample_big = np.log(sa_known[:, rec_per])

        ind = np.in1d(spectral_periods, period)
        rec_val = np.exp(sample_big[:, ind]).flatten()
        if is_scaled:
            scale_fac = target / rec_val
        else:
            scale_fac = rec_val / rec_val

        return sample_big, scale_fac, mean_req_exact, sig_req_exact, gm

    # ------------------------------------------------------------------------------------------------------------------

    def evaluate_gmm(self, period):
        med_sa_period, sig_ln_sa_period = self._evaluate_gmm(period)
        # return 1-d arrays
        return med_sa_period.flatten(), sig_ln_sa_period.flatten()

    # ------------------------------------------------------------------------------------------------------------------

    def compute_seismic_hazard_integral(self, for_which, **kwargs):
        im_input = kwargs.get('im_input', np.array([]))

        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)

        period = self.evaluate_period(for_which)

        med_sa_period, sig_ln_sa_period = self.evaluate_gmm(period)  # get back 1-d arrays for scenarios
        to_return = self._compute_seismic_hazard_integral(med_sa_period, sig_ln_sa_period, im_input)

        return to_return

    # ------------------------------------------------------------------------------------------------------------------

    def get_target_spectrum(self, mrp, for_which, **kwargs):
        uncondition = kwargs.get('uncondition', True)
        spectral_periods = kwargs.get('spectral_periods', np.array([])).flatten()  # 1-d array
        if spectral_periods.size == 0:
            spectral_periods = np.logspace(np.log10(0.05), np.log10(5), 50, endpoint=True)

        period = self.evaluate_period(for_which)
        spectral_periods = Utility.merge_sorted_1d_arrays(spectral_periods, [period])
        med_sa_spectrum, sig_ln_sa_spectrum = self._evaluate_gmm(spectral_periods)  # (n_s, n_sp), (n_s, n_sp)

        mean_req, cov_req, spectral_periods = self._get_target_spectrum(mrp, for_which,
                                                                        med_sa_spectrum, sig_ln_sa_spectrum,
                                                                        spectral_periods,
                                                                        uncondition)

        return mean_req, cov_req, spectral_periods

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

        period = self.evaluate_period(for_which)

        ind_period = np.in1d(spectral_periods, period)
        med_sa_period = med_sa_spectrum[:, ind_period].flatten()  # (n_s,)
        sig_ln_sa_period = sig_ln_sa_spectrum[:, ind_period].flatten()  # (n_s,)

        period_mask = spectral_periods != period
        spectral_periods = spectral_periods[period_mask]

        med_sa_spectrum = med_sa_spectrum[:, period_mask]
        sig_ln_sa_spectrum = sig_ln_sa_spectrum[:, period_mask]

        n_s = med_sa_spectrum.shape[0]
        n_sp = len(spectral_periods)

        rho = np.zeros((n_s, n_sp))
        rho[0:, :] = self.correl_func(spectral_periods, period)  # (n_s, n_sp)

        mean_req, cov_req = self.get_conditional_spectrum(mrp, for_which,
                                                          med_sa_spectrum,
                                                          sig_ln_sa_spectrum,
                                                          rho, self.get_moc(spectral_periods),
                                                          med_sa_period, sig_ln_sa_period)
        return mean_req, cov_req, spectral_periods

    # ------------------------------------------------------------------------------------------------------------------

    def select_ground_motion_records(self, n_gm, mean_req, cov_req, spectral_periods, mrp, for_which, **kwargs):
        rng_seed = kwargs.get('rng_seed', None)
        max_scale = kwargs.get('max_scale', 4)
        min_scale = kwargs.get('min_scale', 1 / 3)
        dev_weights = kwargs.get('dev_weights', [1, 2])
        n_loop = kwargs.get('n_loop', 2)
        penalty = kwargs.get('penalty', 0)
        is_scaled = kwargs.get('is_scaled', True)
        not_allowed = kwargs.get('not_allowed', [])
        classify_pulse = kwargs.get('classify_pulse', True)
        sampling_method = kwargs.get('sampling_method', 'mcs')

        ground_motion_records = self._select_ground_motion_records(n_gm, mean_req, cov_req, mrp, for_which,
                                                                   spectral_periods,
                                                                   rng_seed, max_scale, min_scale, dev_weights, n_loop,
                                                                   penalty, is_scaled, not_allowed, classify_pulse,
                                                                   sampling_method)

        return ground_motion_records

    # ------------------------------------------------------------------------------------------------------------------

    def plot_gm_psa_spectra(self, for_which, **kwargs):
        lc = kwargs.get('lc', 'black')
        lw = kwargs.get('lw', 1)
        save_dir_path = kwargs.get('save_dir_path', Utility.get_path(os.getcwd()))
        save_mat = kwargs.get('save_mat', False)

        fig, ax, export_mat_dict = self._plot_gm_psa_spectra(for_which, **kwargs)

        period = self.evaluate_period(for_which)

        ylim = ax.get_ylim()
        ax.plot([period] * 2, ylim, color=lc, linestyle='-', linewidth=lw)
        ax.set_ylim(ylim)

        export_mat_dict['period'] = period

        if save_mat:
            structure = self.structure
            name = structure.name
            for_which_str = '_'.join(for_which)
            file_name = f'plot_gm_psa_spectra_{name}_{for_which_str}.mat'
            sio.savemat(Utility.get_path(save_dir_path, file_name), export_mat_dict)

        return fig, ax, export_mat_dict

    # ------------------------------------------------------------------------------------------------------------------

    def plot_conditional_spectra_pdfs(self, mrp, for_which, **kwargs):
        lc = kwargs.get('lc', 'black')
        lw = kwargs.get('lw', 1)
        save_dir_path = kwargs.get('save_dir_path', Utility.get_path(os.getcwd()))
        save_mat = kwargs.get('save_mat', False)

        fig, ax, export_mat_dict = self._plot_conditional_spectra_pdfs(mrp, for_which, **kwargs)

        structure = self.structure
        period = self.evaluate_period(for_which)

        ylim = ax.get_ylim()
        ax.plot([period] * 2, ylim, 0, color=lc, linestyle='-', linewidth=lw)
        ax.set_ylim(ylim)

        export_mat_dict['period'] = period

        if save_mat:
            name = structure.name
            for_which_str = '_'.join(for_which)
            file_name = f'plot_conditional_spectra_pdfs_{name}_{for_which_str}_mrp_{mrp}.mat'
            sio.savemat(Utility.get_path(save_dir_path, file_name), export_mat_dict)

        return fig, ax, export_mat_dict

    # ------------------------------------------------------------------------------------------------------------------
