# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 16:16:25 2019

@author: Angshuman Deb
"""

import numpy as np
from .sa import Sa
from .utility import Utility
import os
import sys


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

        n_s = med_sa_spectrum.shape[0]
        n_sp = len(spectral_periods)

        rho = np.zeros((n_s, n_sp))
        rho[0:, :] = self.correl_func(spectral_periods, period)  # (n_s, n_sp)

        mean_req, cov_req = self.get_conditional_spectrum(mrp, for_which,
                                                          med_sa_spectrum,
                                                          sig_ln_sa_spectrum,
                                                          rho, self.get_moc(spectral_periods),
                                                          med_sa_period, sig_ln_sa_period)

        cov_req[ind_period, :, :] = 0.
        cov_req[:, ind_period, :] = 0.

        return mean_req, cov_req, spectral_periods

    # ------------------------------------------------------------------------------------------------------------------

    def mark_period_on_axis(self, for_which, ax, **kwargs):
        lc = kwargs.get('period_lc', 'black')
        lw = kwargs.get('period_lw', 1)
        period = self.evaluate_period(for_which)
        ylim = Utility.get_ylim(ax, (0.01, None))
        if '3D' in ax.__class__.__name__:
            ax.plot([period] * 2, ylim, 0, color=lc, linestyle='-', linewidth=lw)
        else:
            ax.plot([period] * 2, ylim, color=lc, linestyle='-', linewidth=lw)
        ax.set_ylim(ylim)
