# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 16:16:25 2019

@author: Angshuman Deb
"""

from .im import IM
import numpy as np
from .structure import Structure
from .utility import Utility
from .multivariate_nataf import multivariate_nataf
from .mixture import mixture
from abc import abstractmethod
import pygmm
from scipy.stats import lognorm
from scipy.stats import skew
import os
import scipy.io as sio
from scipy.stats import norm
import matplotlib.pyplot as plt
from matplotlib import colors
from numba import jit


class Sa(IM):
    large_err = 1.0e99

    ####################################################################################################################
    # Constructor
    ####################################################################################################################

    def __init__(self, structure, gmm, correl_func):
        self.correl_func = correl_func  # function defining correlation b/w Sa at different structural periods
        super().__init__(structure, gmm)

    ####################################################################################################################
    # Abstract methods
    ####################################################################################################################

    @abstractmethod
    def setup_psa_spectra_gms(self, *args, **kwargs):
        pass

    # ------------------------------------------------------------------------------------------------------------------

    @abstractmethod
    def _get_target_spectrum(self, *args, **kwargs):
        pass

    @abstractmethod
    def evaluate_period(self, for_which):
        pass

    @abstractmethod
    def mark_period_on_axis(self, for_which, ax, **kwargs):
        pass

    ####################################################################################################################
    # Public Functionalities (instance methods)
    ####################################################################################################################

    def shc_at_periods(self, periods, **kwargs):
        im_input = kwargs.get('im_input', np.array([]))
        im_input = np.array(im_input).flatten()
        if im_input.size == 0:
            im_input = np.logspace(np.log10(1e-4), np.log10(5), 1000, endpoint=True)
        med_sa_periods, sig_ln_sa_periods = self._evaluate_gmm(periods)
        shc_list = []
        for itr in range(len(periods)):
            prob_dist_params_im = np.row_stack(
                (sig_ln_sa_periods[:, itr], np.abs(0 * med_sa_periods[:, itr]), med_sa_periods[:, itr]))
            shc_list.append(self._compute_seismic_hazard_integral(lognorm, prob_dist_params_im, im_input)['seismic_hazard_curve'])
        return shc_list

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
        deterministic = kwargs.get('deterministic', False)

        ground_motion_records = self._select_ground_motion_records(deterministic, n_gm, mean_req, cov_req, mrp, for_which,
                                                                   spectral_periods,
                                                                   rng_seed, max_scale, min_scale, dev_weights, n_loop,
                                                                   penalty, is_scaled, not_allowed, classify_pulse,
                                                                   sampling_method)

        return ground_motion_records

    # ------------------------------------------------------------------------------------------------------------------

    def get_target_spectrum(self, mrp, for_which, **kwargs):
        uncondition = kwargs.get('uncondition', True)
        spectral_periods = kwargs.get('spectral_periods', np.array([]))
        spectral_periods = np.array(spectral_periods).flatten()
        if spectral_periods.size == 0:
            spectral_periods = np.logspace(np.log10(0.05), np.log10(5), 50, endpoint=True)
        uhs = kwargs.get('uhs', False)

        spectral_periods = Utility.merge_sorted_1d_arrays(spectral_periods, self.evaluate_period(for_which))
        if uhs:
            temp = self.uniform_hazard_spectrum(mrp, spectral_periods=spectral_periods)
            return np.log(temp[:, 1][np.newaxis, :]), np.zeros(shape=(temp.shape[0], temp.shape[0], 1)), temp[:, 0]

        med_sa_spectrum, sig_ln_sa_spectrum = self._evaluate_gmm(spectral_periods)  # (n_s, n_sp), (n_s, n_sp)
        mean_req, cov_req, spectral_periods = self._get_target_spectrum(mrp, for_which,
                                                                        med_sa_spectrum, sig_ln_sa_spectrum,
                                                                        spectral_periods,
                                                                        uncondition)
        return mean_req, cov_req, spectral_periods

    # ------------------------------------------------------------------------------------------------------------------

    def compute_seismic_hazard_integral(self, for_which, **kwargs):
        im_input = kwargs.get('im_input', np.array([]))
        im_input = np.array(im_input).flatten()
        if im_input.size == 0:
            im_input = np.logspace(np.log10(1e-4), np.log10(5), 1000, endpoint=True)

        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)

        med_sa, sig_ln_sa = self.evaluate_gmm(self.evaluate_period(for_which))
        prob_dist_params_im = np.row_stack((sig_ln_sa, np.abs(0 * med_sa), med_sa))
        to_return = self._compute_seismic_hazard_integral(lognorm, prob_dist_params_im, im_input)

        return to_return

    # ------------------------------------------------------------------------------------------------------------------

    def uniform_hazard_spectrum(self, mrp, **kwargs):
        spectral_periods = kwargs.get('spectral_periods', np.array([]))
        spectral_periods = np.array(spectral_periods).flatten()
        if spectral_periods.size == 0:
            spectral_periods = np.logspace(np.log10(0.05), np.log10(5), 50, endpoint=True)

        im_input = kwargs.get('im_input', np.array([]))
        im_input = np.array(im_input).flatten()
        if im_input.size == 0:
            im_input = np.logspace(np.log10(1e-4), np.log10(5), 1000, endpoint=True)

        uhs = []
        shc_list = self.shc_at_periods(spectral_periods, im_input=im_input)
        for itr in range(len(spectral_periods)):
            uhs.append(IM.invert_seismic_hazard(shc_list[itr], mrp))
        return np.column_stack((spectral_periods, uhs))

    def get_moc(self, periods):
        n_p = len(periods)
        moc = np.zeros((n_p, n_p))
        correl_func = self.correl_func

        for j in range(n_p):
            moc[:, j] = correl_func(periods, periods[j])

        return moc

    # ------------------------------------------------------------------------------------------------------------------

    def _gmm(self, scenario_data, periods):
        gmm = self.gmm
        scenario = pygmm.Scenario(**scenario_data)
        return gmm(scenario).interp_spec_accels(periods), gmm(scenario).interp_ln_stds(periods)

    # ------------------------------------------------------------------------------------------------------------------

    def _evaluate_gmm(self, periods):

        if np.isscalar(periods):
            periods = np.array([periods])
        else:
            periods = np.array(periods)

        scenario_data = self.structure.get_site_hazard_info()['scenario_data']

        n_p = len(periods)
        n_s = len(list(scenario_data.values())[0])

        med_sa_periods = np.zeros((n_s, n_p))
        sig_ln_sa_periods = np.zeros((n_s, n_p))

        scenarios = [{key: value for key, value in zip(scenario_data.keys(), values)} for values in zip(*scenario_data.values())]

        for i in range(n_s):
            med_sa_periods[i, :], sig_ln_sa_periods[i, :] = self._gmm(scenarios[i], periods)

        to_return = list()
        to_return.append(med_sa_periods)
        to_return.append(sig_ln_sa_periods)

        return to_return

    # ------------------------------------------------------------------------------------------------------------------

    def get_conditional_spectrum(self, mrp, for_which, med_sa_spectrum, sig_ln_sa_spectrum,
                                 rho, rho_all,
                                 med_sa_star, sig_ln_sa_star):

        # med_sa_spectrum : (n_s, n_sp)
        # sig_ln_sa_spectrum : (n_s, n_sp)
        # rho : (n_s, n_sp)
        # rho_all : (n_sp, n_sp)
        # med_sa_star : (n_s,)
        # sig_ln_sa_star : (n_s,)

        sa_cond_value = self.get_inv_seismic_hazard(mrp, for_which)
        n_s = len(med_sa_star)
        n_sp = rho_all.shape[0]

        # Use the CMS method to estimate target means and covariances
        eps_bar = (np.log(sa_cond_value) - np.log(med_sa_star)) / sig_ln_sa_star  # (n_s,)
        mean_req = np.log(med_sa_spectrum) + sig_ln_sa_spectrum * eps_bar.reshape(n_s, 1) * rho  # (n_s, n_sp)

        sig_1 = np.einsum('ij,il->jli', sig_ln_sa_spectrum, sig_ln_sa_spectrum) * rho_all.reshape(n_sp, n_sp, 1)
        sig_2 = rho * sig_ln_sa_spectrum * sig_ln_sa_star.reshape(n_s, 1)
        sig_2_sig_2_t = np.einsum('ij,il->jli', sig_2, sig_2)
        term_2 = sig_2_sig_2_t * (1 / (sig_ln_sa_star ** 2)).reshape(1, 1, n_s)

        cov_req = sig_1 - term_2  # (n_sp, n_sp, n_s)
        np.einsum('iij->ji', cov_req)[...] = np.abs(np.einsum('iij->ji', cov_req))

        return mean_req, cov_req

    # ------------------------------------------------------------------------------------------------------------------

    def _select_ground_motion_records(self, deterministic, n_gm, mean_req, cov_req, mrp, for_which, spectral_periods,
                                      rng_seed, max_scale, min_scale, dev_weights, n_loop,
                                      penalty, is_scaled, not_allowed, classify_pulse, sampling_method):

        # Extract rec_data
        rec_data_file_path = Utility.get_path(os.path.dirname(__file__), 'Rec_Data', 'rec_data.pickle')
        rec_data = Utility.pickle_load_dict(rec_data_file_path)
        file_name_1 = rec_data['file_name_1']
        file_name_2 = rec_data['file_name_2']
        file_name_vert = rec_data['file_name_vert']
        per_known = rec_data['per_known']
        sa_known = rec_data['sa_known']
        soil_vs_30 = rec_data['soil_vs_30']
        not_allowed_rec_num = rec_data['not_allowed_rec_num']
        not_allowed_ind = np.unique([i - 1 for i in [*not_allowed_rec_num, *not_allowed]])
        pulse_like_records = rec_data['pulse_like_records']

        if deterministic:
            cov_req[...] = 0

        if np.all(cov_req == 0):
            deterministic = True
        else:
            deterministic = False

        if mean_req.shape[0] > 1 and cov_req.shape[2] > 1:
            is_mixture = True
        else:
            is_mixture = False

        if is_mixture:
            scenario_weights = self.get_seismic_hazard_deagg(mrp, for_which)
            threshold = 1.0e-5
            scenario_mask = scenario_weights >= threshold
            scenario_weights = scenario_weights[scenario_mask]
            scenario_weights = scenario_weights / np.sum(scenario_weights)
        else:
            scenario_weights = np.array([1.0])
            scenario_mask = np.array([True])

        if not deterministic:
            # 2 options:
            # 1. mixture CMS
            # 2. non-mixture CMS
            ind = np.where(np.all(cov_req[:, :, 0] == 0, axis=1))[0][0]
            cov_req = np.delete(np.delete(cov_req, ind, axis=0), ind, axis=1)
            mean_req = np.delete(mean_req, ind, axis=1)
            spectral_periods = np.delete(spectral_periods, ind)
            n_sp = len(spectral_periods)
            # Simulate realizations of ground motion spectra
            n_trials = 20
            if is_mixture:
                cov_req = cov_req[:, :, scenario_mask]
                mean_req = mean_req[scenario_mask, :]
                sig_req = np.sqrt(np.einsum('iik->ki', cov_req))
                mean_req_exact, sig_req_exact = Utility.get_exact_mixture_params(mean_req, sig_req, scenario_weights)
                mean_req_exact = mean_req_exact.flatten()
                sig_req_exact = sig_req_exact.flatten()
                prob_dist_params = np.hstack((np.einsum('ijd->idj', mean_req[:, :, np.newaxis]),
                                              np.einsum('ijd->idj', sig_req[:, :, np.newaxis])))
                mixture_list = list()
                for i_sp in range(n_sp):
                    mixture_list.append(mixture(prob_dist_params[:, :, i_sp], weights=scenario_weights, each=norm))
                m = multivariate_nataf(mixture_list, cov_req / np.einsum('ij,ik->jki', sig_req, sig_req))
                gms = np.exp(m.rvs(size=(n_gm, n_trials), random_state=rng_seed, method=sampling_method))
            else:
                mean_req_exact = mean_req[0, :]
                sig_req = np.sqrt(np.diag(cov_req[:, :, 0]))[np.newaxis, :]
                sig_req_exact = sig_req.flatten()
                prob_dist_list = list()
                for i_sp in range(n_sp):
                    prob_dist_list.append(norm(mean_req[0, i_sp], sig_req[0, i_sp]))
                corr_matrix = cov_req[:, :, 0] / np.einsum('i,j->ij', sig_req_exact, sig_req_exact)
                np.fill_diagonal(corr_matrix, 1.0)
                m = multivariate_nataf(prob_dist_list, corr_matrix, corr_matrix_z=corr_matrix)
                gms = np.exp(m.rvs(size=(n_gm, n_trials), random_state=rng_seed, method=sampling_method))
            # gms : (n_gm : 0, n_trials: 1, n_sp: 2)
            # to
            # gms : (n_gm, n_sp, n_trials)
            gms = gms.transpose((0, 2, 1))
            dev_mean_sim = dev_weights[0] * np.sum(np.square(np.mean(np.log(gms), axis=0).T - mean_req_exact), axis=1)
            dev_skew_sim = 0.1 * (dev_weights[0] + dev_weights[1]) * np.sum(np.square(skew(np.log(gms), axis=0).T), axis=1)
            dev_sig_sim = dev_weights[1] * np.sum(np.square(np.std(np.log(gms), axis=0).T - sig_req_exact), axis=1)
            dev_total_sim = dev_mean_sim + dev_sig_sim + dev_skew_sim
            gm = gms[:, :, np.argmin(np.abs(dev_total_sim))]
        else:
            # 3 options:
            # 1. mixture CMS
            # 2. non-mixture CMS
            # 3. UHS (non-mixture by default), still calling the variable mean_req (abuse of notation)
            if is_mixture:
                mean_req = mean_req[scenario_mask, :]
                mean_req_exact, sig_req_exact = Utility.get_exact_mixture_params(mean_req, np.zeros(mean_req.shape), scenario_weights)
                mean_req_exact = mean_req_exact.flatten()
                sig_req_exact = sig_req_exact.flatten()
            else:
                mean_req_exact = mean_req[0, :]
                sig_req_exact = np.zeros(mean_req_exact.shape)
            gm = np.repeat(np.exp(mean_req_exact)[np.newaxis, :], n_gm, axis=0)

        sample_big, scale_fac, mean_req_exact, sig_req_exact, gm = self.setup_psa_spectra_gms(mrp, for_which,
                                                                                              mean_req_exact,
                                                                                              sig_req_exact,
                                                                                              spectral_periods,
                                                                                              gm, sa_known, per_known,
                                                                                              is_scaled)

        # Find best matches to simulated spectra from ground motion database
        rec_id = np.zeros(n_gm).astype(int) - 1  # 1-d array (-1 because 0 is a valid index in python)
        sample_small = np.zeros((n_gm, sample_big.shape[1]))
        final_scale_fac = np.ones(n_gm)  # 1-d array

        if is_scaled:
            discard_mask = np.any([scale_fac > max_scale, scale_fac < min_scale, soil_vs_30 == -1], axis=0)
        else:
            discard_mask = (soil_vs_30 == -1)

        for i in range(n_gm):
            err = np.sum(
                np.square(
                    np.log(
                        np.exp(sample_big) * scale_fac.reshape(len(scale_fac), 1)
                    ) - np.log(gm[i, :]).reshape(1, gm.shape[1])
                ), axis=1
            )

            err[discard_mask] = self.large_err
            err[not_allowed_ind] = self.large_err
            err[rec_id[rec_id >= 0]] = self.large_err

            min_err = np.min(err)
            rec_id[i] = np.argmin(err)
            if min_err >= self.large_err:
                for_which_str = '_'.join(for_which)
                warning_str = f'Warning: Problem with GMS. No good mathces found for for_which = {for_which_str}!'
                print(warning_str)
            if is_scaled:
                final_scale_fac[i] = scale_fac[rec_id[i]]
            sample_small[i, :] = np.log(np.exp(sample_big[rec_id[i], :]) * scale_fac[rec_id[i]])

        # Greedy subset modification procedure
        rec_id, final_scale_fac = greedy(
            n_loop, rec_id, sample_small, sample_big,
            mean_req_exact, sig_req_exact, np.array(dev_weights), scale_fac, final_scale_fac,
            deterministic, penalty, is_scaled, max_scale, min_scale, soil_vs_30, not_allowed_ind, self.large_err
        )

        final_rsn = rec_id + 1
        final_scale_factors = final_scale_fac

        # Output data
        ground_motion_records = dict()
        ground_motion_records['gm_nbr'] = list()
        ground_motion_records['NGA_RSN'] = list()
        ground_motion_records['scale_factor'] = list()
        ground_motion_records['file_name_1'] = list()
        ground_motion_records['file_name_2'] = list()
        ground_motion_records['file_name_vert'] = list()

        if classify_pulse:
            ground_motion_records['is_pulse'] = list()

        for i in range(n_gm):
            ground_motion_records['gm_nbr'].append(i + 1)
            ground_motion_records['NGA_RSN'].append(final_rsn[i])
            ground_motion_records['scale_factor'].append(final_scale_factors[i])
            ground_motion_records['file_name_1'].append(file_name_1[final_rsn[i] - 1].replace('.at2', '.AT2'))
            ground_motion_records['file_name_2'].append(file_name_2[final_rsn[i] - 1].replace('.at2', '.AT2'))
            ground_motion_records['file_name_vert'].append(file_name_vert[final_rsn[i] - 1].replace('.at2', '.AT2'))
            if classify_pulse:
                if final_rsn[i] in pulse_like_records:
                    ground_motion_records['is_pulse'].append(1)
                else:
                    ground_motion_records['is_pulse'].append(0)

        return ground_motion_records

    # ------------------------------------------------------------------------------------------------------------------

    def plot_gm_psa_spectra(self, for_which, **kwargs):
        figkwargs = kwargs.get('figkwargs', dict())
        grid_alpha = kwargs.get('grid_alpha', 0.5)
        minor_grid_alpha = kwargs.get('minor_grid_alpha', 0.0)
        lc = kwargs.get('lc', 'black')
        lc1 = kwargs.get('lc1', 'red')
        lw = kwargs.get('lw', 1)
        txc = kwargs.get('txc', 'red')
        k = kwargs.get('k', 1.96)
        fig_ax = kwargs.get('fig_ax', None)
        plot_uhs = kwargs.get('plot_uhs', False)
        uhs_lc = kwargs.get('uhs_lc', 'deepskyblue')
        save_dir_path = kwargs.get('save_dir_path', Utility.get_path(os.getcwd()))
        save_mat = kwargs.get('save_mat', False)

        if fig_ax is None or np.array(fig_ax, dtype='object').any() is None:
            fig = plt.figure(**figkwargs)
            ax = fig.gca()
        else:
            fig = fig_ax[0]
            ax = fig_ax[1]
        ax.minorticks_on()
        ax.grid(True, which="major", alpha=grid_alpha)
        ax.grid(True, which="minor", alpha=minor_grid_alpha)

        orig_for_which = for_which
        dir_level_to_seek = 'Hazard_Level_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)

        if len(orig_for_which) > len(for_which):
            gm_num = int(orig_for_which[len(for_which)])
        else:
            gm_num = None

        gm_records, target_spectra, mrp, _ = self.get_gms_results(for_which)
        mean_req, cov_req, spectral_periods = target_spectra

        if mean_req.shape[0] > 1 and len(cov_req.shape) > 2 and cov_req.shape[2] > 1:
            sig_req = np.sqrt(np.einsum('iik->ki', cov_req))
            scenario_weights = self.get_seismic_hazard_deagg(mrp, for_which)
            mean_req_exact, sig_req_exact = Utility.get_exact_mixture_params(mean_req, sig_req, scenario_weights)
            mean_req_exact = mean_req_exact.flatten()
            sig_req_exact = sig_req_exact.flatten()
        else:
            try:
                mean_req_exact = mean_req[0, :]
                sig_req_exact = np.sqrt(np.diag(cov_req[:, :, 0]))
            except IndexError:
                mean_req_exact = mean_req
                sig_req_exact = np.sqrt(np.diag(cov_req[:, :]))

        rec_data_file_path = Utility.get_path(os.path.dirname(__file__), 'Rec_Data', 'rec_data.pickle')
        rec_data = Utility.pickle_load_dict(rec_data_file_path)
        per_known = rec_data['per_known']
        sa_known = rec_data['sa_known']

        # Scale already saved response spectra of eqkes in ensemble. Much faster
        psa = sa_known[np.array(gm_records['NGA_RSN']) - 1, :].T
        scale_fac = np.array(gm_records['scale_factor']).reshape(1, len(gm_records['scale_factor']))
        psa = psa * scale_fac

        for i_nbr in range(len(gm_records['gm_nbr'])):
            ax.plot(per_known, psa[:, i_nbr], color=0.5 * (np.ones(3) - np.array(colors.to_rgb(lc))),
                    alpha=0.5, linewidth=lw / 2)

        if gm_num is not None:
            ax.plot(per_known, psa[:, gm_num - 1], color=lc1, linewidth=lw)

        if plot_uhs:
            uhs = self.uniform_hazard_spectrum(mrp)
            ax.plot(uhs[:, 0], uhs[:, 1], color=uhs_lc, linewidth=lw * 1.5)

        sample_mean = np.mean(np.log(psa), axis=1)
        sample_median = np.exp(sample_mean)
        sample_std = np.std(np.log(psa), axis=1, ddof=1)
        sample_k_lower = np.exp(sample_mean - k * sample_std)
        sample_k_upper = np.exp(sample_mean + k * sample_std)
        ax.plot(per_known, sample_median, color=lc, linestyle='--', linewidth=lw)
        ax.plot(per_known, sample_k_lower, color=lc, linestyle='--', linewidth=lw)
        ax.plot(per_known, sample_k_upper, color=lc, linestyle='--', linewidth=lw)

        median_req_exact = np.exp(mean_req_exact)
        target_median = median_req_exact
        mean_req_exact = np.log(median_req_exact)
        ax.plot(spectral_periods, target_median, color=lc, linestyle='-', linewidth=lw)

        if not np.all(sig_req_exact == 0):
            target_std = sig_req_exact
            target_k_lower = np.exp(mean_req_exact - k * target_std)
            target_k_upper = np.exp(mean_req_exact + k * target_std)
            ax.plot(spectral_periods, target_k_lower, color=lc, linestyle='-', linewidth=lw)
            ax.plot(spectral_periods, target_k_upper, color=lc, linestyle='-', linewidth=lw)

        im_value = self.get_inv_seismic_hazard(mrp, for_which)
        ax.plot(spectral_periods, im_value * np.ones(len(spectral_periods)),
                color=txc, linestyle='-.', linewidth=lw / 2)
        ax.text(spectral_periods[-1], im_value, f'{im_value:.3}', color=txc, horizontalalignment='right',
                verticalalignment='bottom')

        ax.set_xlim([spectral_periods[0], spectral_periods[-1]])

        export_mat_dict = dict()
        export_mat_dict['per_known'] = per_known
        export_mat_dict['psa'] = psa
        if gm_num is not None:
            export_mat_dict['gm_num'] = gm_num
        export_mat_dict['spectral_periods'] = spectral_periods
        export_mat_dict['mean_req'] = mean_req_exact
        export_mat_dict['sig_req'] = sig_req_exact
        export_mat_dict['mrp'] = mrp
        export_mat_dict['im_value'] = im_value

        self.mark_period_on_axis(for_which, ax, **kwargs)
        export_mat_dict['period'] = self.evaluate_period(for_which)

        if save_mat:
            structure = self.structure
            name = structure.name
            for_which_str = '_'.join(for_which)
            file_name = f'plot_gm_psa_spectra_{name}_{for_which_str}.mat'
            sio.savemat(Utility.get_path(save_dir_path, file_name), export_mat_dict)

        return fig, ax, export_mat_dict

    # ------------------------------------------------------------------------------------------------------------------

    def plot_unconditional_spectra_pdfs(self, mrp, for_which, **kwargs):
        periods_lim = kwargs.get('periods_lim', (0.05, 5))
        scenario = kwargs.get('scenario', 'uncondition')
        save_dir_path = kwargs.get('save_dir_path', Utility.get_path(os.getcwd()))
        save_mat = kwargs.get('save_mat', False)
        spectral_periods = kwargs.get('spectral_periods', np.logspace(np.log10(periods_lim[0]),
                                                                      np.log10(periods_lim[-1]),
                                                                      50))

        spectral_periods = np.array(spectral_periods)

        structure = self.structure
        site_hazard_info = structure.get_site_hazard_info()
        scenario_data = site_hazard_info['scenario_data']
        scenario_weights = self.get_seismic_hazard_deagg(mrp, for_which)

        med_sa_spectrum, sig_ln_sa_spectrum = self._evaluate_gmm(spectral_periods)  # (n_s, n_sp), (n_s, n_sp)

        if scenario == 'uncondition':
            mean_ln_sa_spectrum, sig_ln_sa_spectrum = Utility.get_exact_mixture_params(
                np.log(med_sa_spectrum),
                sig_ln_sa_spectrum,
                scenario_weights)  # (1, n_sp), (1, n_sp)
            med_sa_spectrum = np.exp(mean_ln_sa_spectrum).ravel()  # (1, n_sp)
            sig_ln_sa_spectrum = sig_ln_sa_spectrum.ravel()

        if scenario == 'mean':
            numeric_params = [param.name for param in self.gmm.PARAMS if param.__class__.__name__ == 'NumericParameter']
            mean_scenario = {numeric_param: np.sum(scenario_data[numeric_param] * scenario_weights) for numeric_param in numeric_params}
            stacked = np.column_stack([scenario_data[numeric_param] for numeric_param in numeric_params])
            mean_array = np.array(list(mean_scenario.values()))[np.newaxis, :]
            mean_ind = np.argmin(np.sum((stacked - mean_array) ** 2, axis=1))
            mean_scenario = {**mean_scenario, **{key: scenario_data[key][mean_ind] for key in scenario_data.keys() if key not in numeric_params}}
            med_sa_spectrum, sig_ln_sa_spectrum = self._gmm(mean_scenario, spectral_periods)

        if scenario == 'median':
            ind = np.argmin(np.abs(np.cumsum(scenario_weights) - 0.50))
            med_sa_spectrum = med_sa_spectrum[ind, :]
            sig_ln_sa_spectrum = sig_ln_sa_spectrum[ind, :]

        if scenario == 'mode':
            ind = np.argmax(scenario_weights)
            med_sa_spectrum = med_sa_spectrum[ind, :]
            sig_ln_sa_spectrum = sig_ln_sa_spectrum[ind, :]

        if isinstance(scenario, int):
            med_sa_spectrum = med_sa_spectrum[scenario, :]
            sig_ln_sa_spectrum = sig_ln_sa_spectrum[scenario, :]

        if 'spectral_periods' in kwargs.keys():
            kwargs.pop('spectral_periods')
        fig, ax, export_mat_dict = Sa.plot_spectra_pdfs(med_sa_spectrum, sig_ln_sa_spectrum,
                                                        spectral_periods, **kwargs)

        if save_mat:
            name = structure.name
            for_which_str = '_'.join(for_which)
            file_name = f'plot_unconditional_spectra_pdfs_{name}_{for_which_str}_mrp_{mrp}.mat'
            sio.savemat(Utility.get_path(save_dir_path, file_name), export_mat_dict)

        return fig, ax, export_mat_dict

    # ------------------------------------------------------------------------------------------------------------------

    def plot_conditional_spectra_pdfs(self, mrp, for_which, **kwargs):
        periods_lim = kwargs.get('periods_lim', (0.05, 5))
        scenario = kwargs.get('scenario', 'uncondition')
        spectral_periods = kwargs.get('spectral_periods', np.logspace(np.log10(periods_lim[0]),
                                                                      np.log10(periods_lim[-1]),
                                                                      50))
        save_mat = kwargs.get('save_mat', False)
        save_dir_path = kwargs.get('save_dir_path', Utility.get_path(os.getcwd()))

        spectral_periods = np.array(spectral_periods)

        structure = self.structure
        site_hazard_info = structure.get_site_hazard_info()
        scenario_data = site_hazard_info['scenario_data']
        scenario_weights = self.get_seismic_hazard_deagg(mrp, for_which)

        ind = 0
        if scenario == 'uncondition':
            mean_req, cov_req, spectral_periods = self.get_target_spectrum(mrp, for_which,
                                                                           spectral_periods=spectral_periods,
                                                                           uncondition=True)
        else:
            mean_req, cov_req, spectral_periods = self.get_target_spectrum(mrp, for_which,
                                                                           spectral_periods=spectral_periods,
                                                                           uncondition=False)

            if scenario == 'mean':
                numeric_params = [param.name for param in self.gmm.PARAMS if param.__class__.__name__ == 'NumericParameter']
                mean_scenario = {numeric_param: np.sum(scenario_data[numeric_param] * scenario_weights) for numeric_param in numeric_params}
                stacked = np.column_stack([scenario_data[numeric_param] for numeric_param in numeric_params])
                mean_array = np.array(list(mean_scenario.values()))[np.newaxis, :]
                mean_ind = np.argmin(np.sum((stacked - mean_array) ** 2, axis=1))
                mean_scenario = {**mean_scenario, **{key: scenario_data[key][mean_ind] for key in scenario_data.keys() if key not in numeric_params}}
                med_sa_spectrum, sig_ln_sa_spectrum = self._gmm(mean_scenario, spectral_periods)
                mean_req, cov_req, spectral_periods = self._get_target_spectrum(mrp, for_which,
                                                                                med_sa_spectrum[np.newaxis, :],
                                                                                sig_ln_sa_spectrum[np.newaxis, :],
                                                                                spectral_periods, False)

            if scenario == 'median':
                ind = np.argmin(np.abs(np.cumsum(scenario_weights) - 0.50))

            if scenario == 'mode':
                ind = np.argmax(scenario_weights)

            if isinstance(scenario, int):
                ind = scenario

        med_sa_spectrum = np.exp(mean_req[ind, :])
        sig_ln_sa_spectrum = np.sqrt(np.diag(cov_req[:, :, ind]))

        if 'spectral_periods' in kwargs.keys():
            kwargs.pop('spectral_periods')
        fig, ax, export_mat_dict = Sa.plot_spectra_pdfs(med_sa_spectrum, sig_ln_sa_spectrum,
                                                        spectral_periods, **kwargs)

        self.mark_period_on_axis(for_which, ax, **kwargs)
        export_mat_dict['period'] = self.evaluate_period(for_which)

        if save_mat:
            name = structure.name
            for_which_str = '_'.join(for_which)
            file_name = f'plot_conditional_spectra_pdfs_{name}_{for_which_str}_mrp_{mrp}.mat'
            sio.savemat(Utility.get_path(save_dir_path, file_name), export_mat_dict)

        return fig, ax, export_mat_dict

    ####################################################################################################################
    # Static methods
    ####################################################################################################################

    @staticmethod
    def plot_spectra_pdfs(med_sa_spectrum, sig_ln_sa_spectrum, spectral_periods, **kwargs):
        figkwargs = kwargs.get('figkwargs', dict())
        minor_grid = kwargs.get('minor_grid', False)
        lc = kwargs.get('lc', 'black')
        lw = kwargs.get('lw', 1)
        k = kwargs.get('k', 3)
        n_cont = kwargs.get('n_cont', 100)
        periods_lim = kwargs.get('periods_lim', (0.05, 5))
        plot_periods = kwargs.get('plot_periods', np.logspace(np.log10(periods_lim[0]),
                                                              np.log10(periods_lim[-1]),
                                                              5))
        fig_ax = kwargs.get('fig_ax', None)

        if fig_ax is None or np.array(fig_ax, dtype='object').any() is None:
            fig = plt.figure(**figkwargs)  # also accepts tight_layout=bool
            ax = fig.add_subplot(projection='3d')
        else:
            fig = fig_ax[0]
            ax = fig_ax[1]
        if minor_grid:
            ax.minorticks_on()
        else:
            ax.minorticks_off()

        export_mat_dict = dict()

        ax.plot(spectral_periods, med_sa_spectrum, 0, color=lc, linestyle='-.', linewidth=lw)
        ax.plot(spectral_periods, med_sa_spectrum / np.exp(1.96 * sig_ln_sa_spectrum), 0,
                color=lc, linestyle='-.', linewidth=lw / 2)
        ax.plot(spectral_periods, med_sa_spectrum * np.exp(1.96 * sig_ln_sa_spectrum), 0,
                color=lc, linestyle='-.', linewidth=lw / 2)

        export_mat_dict['spectral_periods'] = spectral_periods
        export_mat_dict['med_sa_spectrum'] = med_sa_spectrum
        export_mat_dict['sig_ln_sa_spectrum'] = sig_ln_sa_spectrum
        export_mat_dict['plot_periods'] = plot_periods

        for itr in range(len(plot_periods)):
            ind = np.argmin(np.abs(spectral_periods - plot_periods[itr]))
            med_sa = med_sa_spectrum[ind]
            sig_ln_sa = sig_ln_sa_spectrum[ind]
            # ln_sa = np.logspace(np.log10(1e-4), np.log10(med_sa * np.exp(sig_ln_sa * k)), n_cont)
            ln_sa = np.linspace(1.0e-4, med_sa * np.exp(sig_ln_sa * k), n_cont)
            x = plot_periods[itr] * np.ones(n_cont)
            y = ln_sa
            z = lognorm(sig_ln_sa, 0, med_sa).pdf(y)
            ax.plot(x, y, z, color=lc, linestyle='-', linewidth=lw)
            export_mat_dict[f'ln_sa_period_{itr}'] = y
            export_mat_dict[f'pdf_ln_sa_period_{itr}'] = z

        ax.set_zlim((0, ax.get_zlim()[-1]))
        ax.set_xlim(periods_lim)

        return fig, ax, export_mat_dict

    # ------------------------------------------------------------------------------------------------------------------


@jit(nopython=True)
def greedy(n_loop, rec_id, sample_small, sample_big,
           mean_req_exact, sig_req_exact, dev_weights, scale_fac, final_scale_fac,
           deterministic, penalty, is_scaled, max_scale, min_scale, soil_vs_30, not_allowed_ind, large_err):
    n_gm = sample_small.shape[0]
    n_big = sample_big.shape[0]
    for k in range(n_loop):
        for i in range(n_gm):
            min_id = i
            min_dev = large_err
            # sample_small = np.delete(sample_small, i, axis=0)
            # rec_id = np.delete(rec_id, i, axis=0)
            sample_small = np.concatenate((sample_small[:i, :], sample_small[i + 1:, :]), axis=0)
            rec_id = np.concatenate((rec_id[:i], rec_id[i + 1:]), axis=0)

            # Try to add a new spectra to the subset list
            for j in range(n_big):
                temp_sample_small = np.zeros((sample_small.shape[0] + 1, sample_small.shape[1]))
                temp_sample_small[:-1, :] = sample_small
                if is_scaled:
                    # sample_small = np.vstack((sample_small, sample_big[j, :] + np.log(scale_fac[j])))
                    temp_sample_small[-1, :] = sample_big[j, :] + np.log(scale_fac[j])
                else:
                    # sample_small = np.vstack((sample_small, sample_big[j, :]))
                    temp_sample_small[-1, :] = sample_big[j, :]
                sample_small = temp_sample_small

                # Compute deviations from target
                sample_mean = np.sum(sample_small, axis=0) / sample_small.shape[0]
                dev_mean = sample_mean - mean_req_exact
                if not deterministic:
                    sample_std = np.sqrt(np.sum(np.square(sample_small - sample_mean), axis=0) / (sample_small.shape[0] - 1))
                    dev_sig = sample_std - sig_req_exact
                    dev_total = dev_weights[0] * np.sum(np.square(dev_mean)) + dev_weights[1] * np.sum(np.square(dev_sig))
                    # Penalize bad spectra (set penalty to zero if this is not required)
                    if penalty != 0:
                        for m in range(sample_small.shape[0]):
                            dev_total += np.sum(
                                np.exp(sample_small[m, :]) > np.exp(mean_req_exact + 3. * sig_req_exact)
                            ) * penalty
                else:
                    dev_total = dev_weights[0] * np.sum(np.square(dev_mean))

                if is_scaled:
                    if scale_fac[j] > max_scale or scale_fac[j] < min_scale or soil_vs_30[j] == -1 or np.any(not_allowed_ind == j):
                        dev_total += large_err
                else:
                    if soil_vs_30[j] == -1 or np.any(not_allowed_ind == j):
                        dev_total += large_err

                # Should cause improvement and record should not be repeated
                if dev_total < min_dev and not np.any(rec_id == j):
                    min_id = j
                    min_dev = dev_total
                # sample_small = np.delete(sample_small, sample_small.shape[0] - 1, axis=0)
                sample_small = sample_small[:-1, :]

            # Add new element in the right slot
            if is_scaled:
                final_scale_fac[i] = scale_fac[min_id]

            # sample_small = np.vstack((sample_small[:i, :], sample_big[min_id, :] + np.log(scale_fac[min_id]), sample_small[i:, :]))
            # rec_id = np.hstack((rec_id[:i], min_id, rec_id[i:]))
            to_add = sample_big[min_id, :] + np.log(scale_fac[min_id])
            sample_small = np.concatenate((sample_small[:i, :], np.expand_dims(to_add, axis=0), sample_small[i:, :]), axis=0)
            rec_id = np.concatenate((rec_id[:i], np.array([min_id]).astype(rec_id.dtype), rec_id[i:]), axis=0)
    return rec_id, final_scale_fac
