# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 12:36:55 2019

@author: Angshuman Deb
"""

from abc import ABC, abstractmethod
from .structure import Structure
from .utility import Utility
import numpy as np
import scipy.io as sio
import os
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d
import random
from benedict import benedict
from functools import partial
from collections.abc import Sequence
from pathos.multiprocessing import ProcessPool as Pool


class EDP(ABC):
    
    ####################################################################################################################
    # Constructor
    ####################################################################################################################
    
    def __init__(self, structure, tag, recorder_file_storage, haz_req, small_value, large_value):
        self.structure = structure
        self.tag = tag
        self.recorder_file_storage = recorder_file_storage
        self.haz_req = haz_req
        self.small_value = small_value
        self.large_value = large_value

    ####################################################################################################################
    # Setters (instance methods)
    ####################################################################################################################
    
    def set_haz_req(self, haz_req):
        self.haz_req = haz_req

    ####################################################################################################################
    # Static methods
    ####################################################################################################################

    @staticmethod
    def collect_parallel(tag, results_dir_path, edp_strs, run_case):

        for_which = list(run_case.astype(str))
        run_case_str = ' '.join(for_which)

        target_work_dir_path = Utility.get_path(results_dir_path, *for_which, f'EDP_{tag}_Results')
        nltha_status_file_path = Utility.get_path(results_dir_path, *for_which, 'NLTHA_STATUS.txt')
        status = 'FAIL'
        if os.path.isfile(nltha_status_file_path):
            with open(nltha_status_file_path, 'r') as fid:
                status = fid.readline().strip()

        to_return = dict()
        for edp_str in edp_strs:
            edp_file_name = Utility.get_path(target_work_dir_path, f'edp_{tag}_{edp_str}.txt')
            if os.path.isfile(edp_file_name) and status == 'SUCCESS':
                with open(Utility.get_path(
                        target_work_dir_path, f'edp_{tag}_{edp_str}.txt'), 'r') as fid:
                    edp_val = float(fid.readline().strip('\n'))
            else:
                edp_val = np.inf

            to_return[f'{edp_str}'] = edp_val

        print(f'Collecting edp {tag} for case {run_case_str} complete!')

        return f'{run_case_str}', to_return

    ####################################################################################################################
    # Abstract methods
    ####################################################################################################################

    @abstractmethod
    def write_evaluate_edp_and_helpers(self, *args, **kwargs):
        pass

    @abstractmethod
    def generate_recorder(self, *args, **kwargs):
        pass

    @abstractmethod
    def get_edp_strings(self, *args, **kwargs):
        pass

    ####################################################################################################################
    # Public functionalities
    ####################################################################################################################

    def collect(self, analysis_case, pool_size, **kwargs):
        analysis_list = kwargs.get('analysis_list', 'analysis_list_original.txt')
        structure = self.structure
        setup_dir_path = Utility.get_path(structure.model_work_dir_path, 'Work_Dir', 'NLTHA_Setup')

        if isinstance(analysis_list, str):
            list_file_path = Utility.get_path(setup_dir_path, analysis_case, analysis_list)
            analysis_list = Utility.read_numpy_array_from_txt_file(list_file_path, skiprows='find').astype(int)

        if analysis_list.size == 0:
            return

        results_dir_path = Utility.get_path(structure.model_work_dir_path, 'Work_Dir', 'NLTHA_Results')
        analysis_range = np.arange(1, analysis_list.shape[0] + 1)
        to_pass = [[analysis_list[analysis_index - 1, :] for analysis_index in analysis_range], ]

        design_num_list = list(np.unique(analysis_list[:, 1]).astype(str))

        edp_strs = self.get_edp_strings(list(analysis_list[0, :].astype(str)))
        tag = self.tag
        to_run = partial(EDP.collect_parallel, tag, results_dir_path, edp_strs)

        if pool_size > 1:
            my_pool = Pool(pool_size)
            edp_results_this = dict(my_pool.map(to_run, *to_pass))
            my_pool.close()
            my_pool.join()
            my_pool.clear()
        else:
            edp_results_this = dict(list(map(to_run, *to_pass)))

        for key in edp_results_this.keys():
            for edp_str in edp_strs:
                if edp_results_this[key][f'{edp_str}'] == np.inf:
                    warning_str = f"Warning: Invalid edp_{tag}_{edp_str} " \
                                  f"found for for_which = {key}!"
                    print(warning_str)
                    for_which = key.split()
                    nltha_status_file_path = Utility.get_path(
                        results_dir_path, *for_which, 'NLTHA_STATUS.txt')
                    with open(nltha_status_file_path, 'w') as fid:
                        fid.write('FAIL\n')
                if edp_results_this[key][f'{edp_str}'] == 0:
                    print(f"Warning: edp_{tag}_{edp_str} = 0.0 found for for_which = {key}!")
                    edp_results_this[key][f'{edp_str}'] = self.small_value

        for design_num in design_num_list:
            target_work_dir_path = Utility.get_path(results_dir_path.replace('_Results', '_EDP'), analysis_case,
                                                    design_num)

            if not os.path.isdir(target_work_dir_path):
                os.makedirs(target_work_dir_path)

            edp_results_file_path = Utility.get_path(target_work_dir_path, 'edp_results.pickle')

            if os.path.isfile(edp_results_file_path):
                edp_results = benedict(Utility.pickle_load_dict(edp_results_file_path))
            else:
                edp_results = benedict()

            edp_results_this_design_num = {key: edp_results_this[key] for key in edp_results_this.keys() if
                                           key.split()[Structure.get_dir_level_index('Design_Num')] == design_num}

            edp_results.merge({tag: edp_results_this_design_num})
            Utility.pickle_dump_dict(edp_results_file_path, edp_results)
        return

    def get_rec_save_dir_path(self):
        if self.recorder_file_storage == 'shared':
            rec_save_dir_path = Utility.get_path(self.structure.structural_analysis_platform.nltha_rec_dir_name, 'EDP')
        else:
            rec_save_dir_path = Utility.get_path(self.structure.structural_analysis_platform.nltha_rec_dir_name, f'EDP_{self.tag}')
        return rec_save_dir_path

    def compute_demand_hazard_integral(self, for_which, haz_lev_list, im, n_gm_list=None, **kwargs):
        delta_input = kwargs.get('delta_input', np.array([]))
        delta_input = np.array(delta_input).flatten()  # 1-d array
        min_max_scale_fac = kwargs.get('min_max_scale_fac', [1, 1])  # list

        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)

        mrp_list = list()
        n_gm_list_send = list()

        for haz_lev in haz_lev_list:
            _, _, mrp, n_gm = im.get_gms_results(for_which + [haz_lev])
            mrp_list.append(mrp)
            n_gm_list_send.append(n_gm)

        if isinstance(n_gm_list, Sequence):
            n_gm_list_send = n_gm_list

        seismic_hazard_curve, _ = im.get_psha_results(for_which)

        inv_seismic_hazard_at_haz_lev_list = list()
        for itr in range(len(haz_lev_list)):
            inv_seismic_hazard_at_haz_lev_list.append(
                seismic_hazard_curve[np.argmin(np.abs(seismic_hazard_curve[:, 1] - (1 / mrp_list[itr]))), 0])

        im_vals = seismic_hazard_curve[:, 0]
        im_vals_integration = np.exp(np.mean(np.log(np.column_stack([im_vals[:-1], im_vals[1:]])), axis=1))
        nu_im = seismic_hazard_curve[:, 1]
        dnu_im = np.abs(np.diff(nu_im))

        to_return = benedict()
        edp_strs = self.get_edp_strings(for_which)
        for edp_str in edp_strs:
            dhc_save, nu_edp_im_save = self._compute_demand_hazard_integral(for_which, haz_lev_list,
                                                                            n_gm_list_send,
                                                                            inv_seismic_hazard_at_haz_lev_list,
                                                                            im_vals_integration, dnu_im,
                                                                            edp_str, delta_input,
                                                                            min_max_scale_fac)

            to_return.merge(
                {f'{edp_str}': {'demand_hazard_curve': dhc_save,
                                'demand_hazard_curve_deagg_im': nu_edp_im_save
                                }
                 }
            )

        to_return.merge(
            {'seismic_hazard_curve': np.column_stack(
                [im_vals_integration,
                 np.exp(np.mean(np.log(np.column_stack([nu_im[:-1], nu_im[1:]])), axis=1)),
                 dnu_im]
            )
            }
        )

        to_return.merge({'mrp_list': mrp_list})

        return to_return

    # ------------------------------------------------------------------------------------------------------------------

    def get_inv_demand_hazard(self, mrp, for_which, edp_str):
        dhc, _, _, _ = self.get_psdemha_results(for_which, edp_str)
        return dhc[np.argmin(np.abs(dhc[:, 1] - 1 / mrp)), 0]

    # ------------------------------------------------------------------------------------------------------------------

    def get_demand_hazard_deagg_im(self, mrp, for_which, edp_str):
        dhc, dhc_deagg_im, _, _ = self.get_psdemha_results(for_which, edp_str)
        ind = np.argmin(np.abs(dhc[:, 1] - 1 / mrp))
        # return 1-d array
        return dhc_deagg_im[ind, :] / dhc[ind, 1]

    # ------------------------------------------------------------------------------------------------------------------

    def get_psdemha_results(self, for_which, edp_str=None):
        model_work_dir_path = self.structure.model_work_dir_path
        dir_category = 'PSDemHA_Results'
        file_path = Utility.get_path(model_work_dir_path, 'Work_Dir', dir_category,
                                     *Structure.get_modified_for_which(for_which, 'Design_Num'), 'psdemha_results.pickle')
        psdemha_results = Utility.pickle_load_dict(file_path)
        if len(for_which) == (Structure.get_dir_level_index('Design_Num') + 1) or edp_str is None:
            return psdemha_results[self.tag]

        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)
        run_case_str = ' '.join(for_which)
        psdemha_results_for_which = psdemha_results[self.tag][run_case_str]

        to_return = list()
        to_return.append(psdemha_results_for_which[f'{edp_str}']['demand_hazard_curve'])
        to_return.append(psdemha_results_for_which[f'{edp_str}']['demand_hazard_curve_deagg_im'])
        to_return.append(psdemha_results_for_which['seismic_hazard_curve'])
        to_return.append(psdemha_results_for_which['mrp_list'])

        return to_return

    # ------------------------------------------------------------------------------------------------------------------

    def get_edp_at_hazard_level(self, for_which, haz_lev_list, n_gm_list, edp_str):
        tag = self.tag
        edp_vals = list()

        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)

        edp_results = self.get_edp_results(for_which)

        random.seed(999)
        for haz_lev in haz_lev_list:
            temp = [edp_results[key][f'{edp_str}'] for key in edp_results.keys() if
                    key.split()[Structure.get_dir_level_index('Hazard_Level_Num')] == haz_lev and
                    Structure.get_modified_for_which(key.split(), dir_level_to_seek) == for_which]
            temp = [temp_content for temp_content in temp if temp_content > 0.0]
            temp = [temp_content for temp_content in temp if temp_content < np.inf]
            if len(temp) != n_gm_list[haz_lev_list.index(haz_lev)]:
                warning_str = f"Warning: Number of edp_{tag}_{edp_str} values found for haz_lev = {haz_lev} " \
                              f"doesn't match " \
                              f"n_gm = {n_gm_list[haz_lev_list.index(haz_lev)]} for for_which = {for_which}!"
                print(warning_str)
                try:
                    temp = random.sample(temp, n_gm_list[haz_lev_list.index(haz_lev)])
                except ValueError:
                    pass
            edp_vals.append(np.array(temp))

        return edp_vals

    # ------------------------------------------------------------------------------------------------------------------

    def get_conditional_demand_params(self, for_which, haz_lev_list, n_gm_list,
                                      inv_seismic_hazard_at_haz_lev_list, im_vals_integration,
                                      edp_str):
        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)
        edp_vals = self.get_edp_at_hazard_level(for_which, haz_lev_list, n_gm_list, edp_str)

        haz_req = self.haz_req
        fit_dist = haz_req['fit_dist']
        transf_on_dist_params_list = haz_req['transf_on_dist_params_list']
        transf_to_dist_params_list = haz_req['transf_to_dist_params_list']
        interp_exterp_function_list = haz_req['interp_exterp_function_list']

        dist_params = np.zeros((len(haz_lev_list), len(transf_on_dist_params_list)))
        transf_params = np.zeros((len(haz_lev_list), len(transf_on_dist_params_list)))

        for itr in range(len(haz_lev_list)):
            dist_params[itr, :] = Utility.fit_dist_to_data(fit_dist, edp_vals[itr])

        for itr in range(transf_params.shape[1]):
            transf_params[:, itr] = transf_on_dist_params_list[itr](*dist_params.T)

        transf_params_interp_exterp = np.zeros((len(im_vals_integration), transf_params.shape[1]))
        dist_params_interp_exterp = np.zeros((len(im_vals_integration), transf_params.shape[1]))

        for itr in range(transf_params.shape[1]):
            transf_params_interp_exterp[:, itr] = interp_exterp_function_list[itr](
                inv_seismic_hazard_at_haz_lev_list, transf_params[:, itr], im_vals_integration)

        for itr in range(transf_params.shape[1]):
            dist_params_interp_exterp[:, itr] = transf_to_dist_params_list[itr](*transf_params_interp_exterp.T)

        to_return = dict()
        to_return['dist_params'] = dist_params
        to_return['transf_params'] = transf_params
        to_return['dist_params_interp_exterp'] = dist_params_interp_exterp
        to_return['transf_params_interp_exterp'] = transf_params_interp_exterp
        to_return['edp_vals'] = edp_vals

        return to_return

    # ------------------------------------------------------------------------------------------------------------------

    def _compute_demand_hazard_integral(self, for_which, haz_lev_list, n_gm_list, inv_seismic_hazard_at_haz_lev_list,
                                        im_vals_integration, dnu_im, edp_str, delta_input, min_max_scale_fac):

        fit_dist = self.haz_req['fit_dist']
        conditional_demand_params = self.get_conditional_demand_params(for_which,
                                                                       haz_lev_list,
                                                                       n_gm_list,
                                                                       inv_seismic_hazard_at_haz_lev_list,
                                                                       im_vals_integration, edp_str)

        dist_params_interp_exterp = conditional_demand_params['dist_params_interp_exterp']
        edp_vals = conditional_demand_params['edp_vals']

        if delta_input.size == 0:
            min_edp = min_max_scale_fac[0] * min([min(temp) for temp in edp_vals])
            max_edp = min_max_scale_fac[1] * max([max(temp) for temp in edp_vals])
            delta_input = np.logspace(np.log10(min_edp), np.log10(max_edp), 1000, endpoint=True)
        delta_input = delta_input.reshape(len(delta_input), 1)  # 2-d vector

        nu_edp_im = (1 - fit_dist.cdf(delta_input, *dist_params_interp_exterp.T)) * dnu_im
        nu_edp = np.sum(nu_edp_im, axis=1)
        dhc = np.column_stack([delta_input, nu_edp])

        dhc_save = dhc[dhc[:, 1] != 0, :]
        nu_edp_im_save = nu_edp_im[dhc[:, 1] != 0, :]

        return dhc_save, nu_edp_im_save

    # ------------------------------------------------------------------------------------------------------------------

    def get_edp_results(self, for_which):
        model_work_dir_path = self.structure.model_work_dir_path
        tag = self.tag
        dir_category = 'NLTHA_EDP'
        dir_level_to_seek = 'Design_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)
        edp_results = Utility.pickle_load_dict(
            Utility.get_path(
                model_work_dir_path, 'Work_Dir', dir_category, *for_which, 'edp_results.pickle'
            )
        )[tag]

        return edp_results

    # ------------------------------------------------------------------------------------------------------------------

    def find_anomalous_run_cases(self, for_which, threshold):

        dir_level_to_seek = 'Design_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)

        edp_results = self.get_edp_results(for_which)
        run_cases = [np.array(key.split()).astype(int) for key in edp_results.keys()
                     if np.any(np.array(list(edp_results[key].values())) > threshold)]
        run_cases = np.array(run_cases)

        if run_cases.shape[0] == 0:
            return run_cases.reshape(0, len(Structure.get_all_dir_levels()))
        else:
            return run_cases

    # ------------------------------------------------------------------------------------------------------------------

    def plot_conditional_demand_regression_model(self, for_which, edp_str,
                                                 haz_lev_list, im, **kwargs):

        figkwargs = kwargs.get('figkwargs', dict())
        lc = kwargs.get('lc', 'blue')
        mc = kwargs.get('mc', 'red')
        lw = kwargs.get('lw', 1)
        grid_alpha = kwargs.get('grid_alpha', 0.5)
        minor_grid_alpha = kwargs.get('minor_grid_alpha', 0.0)
        im_lim = kwargs.get('im_lim', list())
        save_dir_path = kwargs.get('save_dir_path', Utility.get_path(os.getcwd()))
        save_mat = kwargs.get('save_mat', False)
        fig_ax = kwargs.get('fig_ax', None)

        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)

        mrp_list = list()
        n_gm_list = list()

        for haz_lev in haz_lev_list:
            _, _, mrp, n_gm = im.get_gms_results(for_which + [haz_lev])
            mrp_list.append(mrp)
            n_gm_list.append(n_gm)

        haz_req = self.haz_req
        num_params = len(haz_req['interp_exterp_function_list'])

        if fig_ax is None or np.array(fig_ax, dtype='object').any() is None:
            fig, axs = plt.subplots(1, num_params, **figkwargs)  # also accepts tight_layout=bool
        else:
            fig = fig_ax[0]
            axs = fig_ax[1]
        for ax in axs:
            ax.minorticks_on()
            ax.grid(True, which="major", alpha=grid_alpha)
            ax.grid(True, which="minor", alpha=minor_grid_alpha)

        inv_seismic_hazard_at_haz_lev_list = list()
        for itr in range(len(haz_lev_list)):
            inv_seismic_hazard_at_haz_lev_list.append(im.get_inv_seismic_hazard(mrp_list[itr], for_which))

        seismic_hazard_curve, _ = im.get_psha_results(for_which)
        im_vals = seismic_hazard_curve[:, 0]
        im_vals_integration = np.exp(np.mean(np.log(np.column_stack([im_vals[:-1], im_vals[1:]])), axis=1))

        if len(im_lim) == 0:
            im_lim = [im_vals_integration[0], im_vals_integration[-1]]

        conditional_demand_params = self.get_conditional_demand_params(for_which,
                                                                       haz_lev_list,
                                                                       n_gm_list,
                                                                       inv_seismic_hazard_at_haz_lev_list,
                                                                       im_vals_integration, edp_str)

        transf_params = conditional_demand_params['transf_params']
        transf_params_interp_exterp = conditional_demand_params['transf_params_interp_exterp']

        for itr in range(transf_params.shape[1]):
            axs[itr].plot(inv_seismic_hazard_at_haz_lev_list, transf_params[:, itr], 'o', color=mc)
            axs[itr].plot(*Utility.get_xy_within_lims(im_vals_integration,
                                                      transf_params_interp_exterp[:, itr],
                                                      im_lim, [-np.inf, np.inf]), '-', color=lc, linewidth=lw)
            axs[itr].set_xlim(im_lim)

        export_mat_dict = dict()
        export_mat_dict[f'inv_seismic_hazard_at_haz_lev_list'] = inv_seismic_hazard_at_haz_lev_list
        export_mat_dict[f'im_vals_integration'] = im_vals_integration
        export_mat_dict[f'transf_params'] = transf_params
        export_mat_dict[f'transf_params_interp_exterp'] = transf_params_interp_exterp

        if save_mat:
            structure = self.structure
            name = structure.name
            tag = self.tag
            for_which_str = '_'.join(for_which)
            file_name = f'plot_conditional_demand_regression_model_{name}_{for_which_str}_' \
                        f'edp_{tag}_{edp_str}.mat'
            sio.savemat(Utility.get_path(save_dir_path, file_name), export_mat_dict)

        return fig, axs, export_mat_dict

    # ------------------------------------------------------------------------------------------------------------------

    def plot_dhc(self, for_which, edp_str, **kwargs):

        figkwargs = kwargs.get('figkwargs', dict())
        save_dir_path = kwargs.get('save_dir_path', Utility.get_path(os.getcwd()))
        save_mat = kwargs.get('save_mat', False)
        grid_alpha = kwargs.get('grid_alpha', 0.5)
        minor_grid_alpha = kwargs.get('minor_grid_alpha', 0.0)
        lc = kwargs.get('lc', 'black')
        lw = kwargs.get('lw', 1.2)
        fig_ax = kwargs.get('fig_ax', None)

        if fig_ax is None or np.array(fig_ax, dtype='object').any() is None:
            fig = plt.figure(**figkwargs)
            ax = fig.gca()
        else:
            fig = fig_ax[0]
            ax = fig_ax[1]
        ax.minorticks_on()
        ax.grid(True, which="major", alpha=grid_alpha)
        ax.grid(True, which="minor", alpha=minor_grid_alpha)

        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)

        dhc = self.get_psdemha_results(for_which, edp_str)[0]
        ax.loglog(dhc[:, 0], dhc[:, 1], color=lc, linewidth=lw)

        export_mat_dict = dict()
        export_mat_dict['dhc'] = dhc

        if save_mat:
            structure = self.structure
            name = structure.name
            tag = self.tag
            for_which_str = '_'.join(for_which)
            file_name = f'plot_dhc_{name}_{for_which_str}_edp_{tag}_{edp_str}.mat'
            sio.savemat(Utility.get_path(save_dir_path, file_name), export_mat_dict)

        return fig, ax, export_mat_dict

    # ------------------------------------------------------------------------------------------------------------------

    def plot_conditional_demand_model_3d(self, for_which, edp_str, haz_lev_list, im, **kwargs):

        figkwargs = kwargs.get('figkwargs', dict())
        mc = kwargs.get('mc', 'red')
        ms = kwargs.get('ms', 2)
        lw = kwargs.get('lw', 1)
        lc = kwargs.get('lc', 'blue')
        save_dir_path = kwargs.get('save_dir_path', Utility.get_path(os.getcwd()))
        save_mat = kwargs.get('save_mat', False)
        edp_lim = kwargs.get('edp_lim', [self.small_value, self.large_value])
        im_lim = kwargs.get('im_lim', list())
        n_bins = kwargs.get('n_bins', 10)
        first_pctl = kwargs.get('first_pctl', 0.025)
        second_pctl = kwargs.get('second_pctl', 0.50)
        third_pctl = kwargs.get('third_pctl', 0.975)
        minor_grid = kwargs.get('minor_grid', False)
        patch_color = kwargs.get('patch_color', 'cyan')
        patch_alpha = kwargs.get('patch_alpha', 0.25)
        n_cont = kwargs.get('n_cont', 100)
        fig_ax = kwargs.get('fig_ax', None)

        if fig_ax is None or np.array(fig_ax, dtype='object').any() is None:
            fig = plt.figure(**figkwargs)
            ax = fig.add_subplot(projection='3d')
        else:
            fig = fig_ax[0]
            ax = fig_ax[1]
        if minor_grid:
            ax.minorticks_on()
        else:
            ax.minorticks_off()

        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)

        mrp_list = list()
        n_gm_list = list()

        for haz_lev in haz_lev_list:
            _, _, mrp, n_gm = im.get_gms_results(for_which + [haz_lev])
            mrp_list.append(mrp)
            n_gm_list.append(n_gm)

        haz_req = self.haz_req
        fit_dist = haz_req['fit_dist']

        inv_seismic_hazard_at_haz_lev_list = list()
        for itr in range(len(haz_lev_list)):
            inv_seismic_hazard_at_haz_lev_list.append(im.get_inv_seismic_hazard(mrp_list[itr], for_which))

        seismic_hazard_curve, _ = im.get_psha_results(for_which)
        im_vals = seismic_hazard_curve[:, 0]
        im_vals_integration = np.exp(np.mean(np.log(np.column_stack([im_vals[:-1], im_vals[1:]])), axis=1))

        if len(im_lim) == 0:
            im_lim = [im_vals_integration[0], im_vals_integration[-1]]

        conditional_demand_params = self.get_conditional_demand_params(for_which,
                                                                       haz_lev_list,
                                                                       n_gm_list,
                                                                       inv_seismic_hazard_at_haz_lev_list,
                                                                       im_vals_integration, edp_str)

        edp_vals = conditional_demand_params['edp_vals']
        dist_params_interp_exterp = conditional_demand_params['dist_params_interp_exterp']
        dist_params = conditional_demand_params['dist_params']

        first_pctl_values = fit_dist.ppf(first_pctl, *dist_params_interp_exterp.T)
        second_pctl_values = fit_dist.ppf(second_pctl, *dist_params_interp_exterp.T)
        second_pctl_values_haz_lev = fit_dist.ppf(second_pctl, *dist_params.T)
        third_pctl_values = fit_dist.ppf(third_pctl, *dist_params_interp_exterp.T)

        ax.plot(*Utility.get_xy_within_lims(first_pctl_values, im_vals_integration, edp_lim, im_lim), 0,
                '-.', color=lc, linewidth=lw)
        ax.plot(*Utility.get_xy_within_lims(second_pctl_values, im_vals_integration, edp_lim, im_lim), 0,
                '-.', color=lc, linewidth=lw)
        ax.plot(second_pctl_values_haz_lev, inv_seismic_hazard_at_haz_lev_list, 0, 'd', color=mc,
                markeredgecolor='black', markersize=ms * 2)
        ax.plot(*Utility.get_xy_within_lims(third_pctl_values, im_vals_integration, edp_lim, im_lim), 0,
                '-.', color=lc, linewidth=lw)

        export_mat_dict = dict()

        for itr in range(len(haz_lev_list)):
            x = edp_vals[itr]
            export_mat_dict[f'edp_vals_haz_lev_{haz_lev_list[itr]}'] = x
            y = np.ones(x.size) * inv_seismic_hazard_at_haz_lev_list[itr]
            ax.plot(*Utility.get_xy_within_lims(x, y, edp_lim, im_lim), 0, '.', color=mc, markersize=ms)
            hist, bin_edges = np.histogram(x, bins=n_bins, density=True)
            export_mat_dict[f'hist_haz_lev_{haz_lev_list[itr]}'] = hist
            export_mat_dict[f'bin_edges_haz_lev_{haz_lev_list[itr]}'] = bin_edges
            for i_bin in range(len(hist)):
                rect = plt.Rectangle((bin_edges[i_bin], 0),
                                     bin_edges[i_bin + 1] - bin_edges[i_bin],
                                     hist[i_bin],
                                     edgecolor='black', facecolor='None', linewidth=lw / 2)
                ax.add_patch(rect)
                art3d.pathpatch_2d_to_3d(rect, z=inv_seismic_hazard_at_haz_lev_list[itr], zdir="y")
            x_cont = np.logspace(np.log10(edp_lim[0]), np.log10(edp_lim[-1]), n_cont)
            y_cont = np.ones(x_cont.size) * inv_seismic_hazard_at_haz_lev_list[itr]
            z_cont = fit_dist.pdf(x_cont, *dist_params[itr, :])
            ax.plot(x_cont, y_cont, z_cont, '-', color=lc, linewidth=lw)
            export_mat_dict[f'cont_pdf_curves_haz_lev_{haz_lev_list[itr]}'] = np.column_stack([x_cont, z_cont])

        rect = plt.Rectangle((edp_lim[0], inv_seismic_hazard_at_haz_lev_list[0]),
                             edp_lim[-1] - edp_lim[0],
                             inv_seismic_hazard_at_haz_lev_list[-1] - inv_seismic_hazard_at_haz_lev_list[0],
                             color=patch_color, alpha=patch_alpha)
        ax.add_patch(rect)
        art3d.pathpatch_2d_to_3d(rect, z=0)

        zlim = ax.get_zlim()
        ax.set_zlim([0, zlim[-1]])
        ax.set_xlim(edp_lim)
        ax.set_ylim(im_lim)

        export_mat_dict[f'inv_seismic_hazard_at_haz_lev_list'] = inv_seismic_hazard_at_haz_lev_list
        export_mat_dict[f'im_vals_integration'] = im_vals_integration
        export_mat_dict[f'first_pctl_values'] = first_pctl_values
        export_mat_dict[f'second_pctl_values'] = second_pctl_values
        export_mat_dict[f'second_pctl_values_haz_lev'] = second_pctl_values_haz_lev
        export_mat_dict[f'third_pctl_values'] = third_pctl_values
        export_mat_dict[f'fit_dist_name'] = fit_dist.name
        export_mat_dict[f'seismic_hazard_curve'] = seismic_hazard_curve
        export_mat_dict[f'dist_params_interp_exterp'] = dist_params_interp_exterp
        export_mat_dict[f'dist_params'] = dist_params

        if save_mat:
            structure = self.structure
            name = structure.name
            tag = self.tag
            for_which_str = '_'.join(for_which)
            file_name = f'plot_conditional_demand_model_3d_{name}_{for_which_str}_' \
                        f'edp_{tag}_{edp_str}.mat'
            sio.savemat(Utility.get_path(save_dir_path, file_name), export_mat_dict)

        return fig, ax, export_mat_dict

    # ------------------------------------------------------------------------------------------------------------------


class FrameMaxDeformation(EDP, ABC):
    def __init__(self, max_what, frame_structure, tag, recorder_file_storage, haz_req, small_value, large_value):
        super().__init__(frame_structure, tag, recorder_file_storage, haz_req, small_value, large_value)
        self.max_what = max_what
        self.resp = 'deformation'
        self.normalize_with = 1.0
        self.resp_index = 1

    @abstractmethod
    def write_generate_edp_recorders(self, for_which, rec_gen_file_open_mode):
        pass

    def write_evaluate_edp_and_helpers(self):
        max_what = self.max_what
        tag = self.tag
        resp = self.resp
        resp_index = self.resp_index
        self.structure.structural_analysis_platform.write_evaluate_edp_and_helpers(
            tag, max_what, resp, resp_index)
        return

    def generate_recorder(self, for_which, **kwargs):
        rec_gen_file_open_mode = kwargs.get('rec_gen_file_open_mode', 'a+')
        gen_rec = kwargs.get('gen_rec', True)

        frame_structure = self.structure
        model_work_dir_path = frame_structure.model_work_dir_path
        dir_category = 'NLTHA_Results'
        dir_level_to_seek = 'Ground_Motion_Rec_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)
        work_dir_path = Utility.get_path(model_work_dir_path, 'Work_Dir', dir_category, Utility.get_path(*for_which))
        if gen_rec:
            self.write_generate_edp_recorders(for_which, rec_gen_file_open_mode)
        frame_structure.structural_analysis_platform.write_evaluate_edps(work_dir_path, self, for_which, gen_rec, rec_gen_file_open_mode)
        return


class MaxColRebarStrain(FrameMaxDeformation):

    ####################################################################################################################
    # Constructor
    ####################################################################################################################

    def __init__(self, max_what, frame_structure, tag, recorder_file_storage, **kwargs):
        haz_req = kwargs.get('haz_req', {})
        small_value = kwargs.get('small_value', 1.0e-10)
        large_value = kwargs.get('large_value', 1.0e0)
        super().__init__(max_what, frame_structure, tag, recorder_file_storage, haz_req, small_value, large_value)
        self.resp = 'strain'
        self.resp_index = 2

    ####################################################################################################################
    # Public functionalities (instance methods)
    ####################################################################################################################

    def write_generate_edp_recorders(self, for_which, rec_gen_file_open_mode):
        frame_structure = self.structure
        model_work_dir_path = frame_structure.model_work_dir_path
        dir_category = 'NLTHA_Results'
        rec_save_dir_path = self.get_rec_save_dir_path()
        dir_level_to_seek = 'Ground_Motion_Rec_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)
        work_dir_path = Utility.get_path(model_work_dir_path, 'Work_Dir', dir_category, Utility.get_path(*for_which))
        frame_structure.write_col_rebar_strain_recorder(for_which, work_dir_path,
                                                        rec_gen_file_open_mode, rec_save_dir_path)
        return

    def get_edp_strings(self, for_which):
        model_params = self.structure.model_params
        num_cols_total = model_params['num_cols_total']
        col_edge_list = model_params['col_edge_list']
        edp_strs = []
        for i_col in range(1, num_cols_total + 1):
            for i_edge in col_edge_list:
                col_id = i_col
                edge_id = i_edge
                edp_strs += [f'col_{col_id}_edge_{edge_id}']
        return edp_strs

    # ------------------------------------------------------------------------------------------------------------------


class MaxSpringDeformation(FrameMaxDeformation):

    ####################################################################################################################
    # Constructor
    ####################################################################################################################

    def __init__(self, spring_type, max_what, frame_structure, tag, recorder_file_storage, **kwargs):
        haz_req = kwargs.get('haz_req', {})
        normalize_with = kwargs.get('normalize_with', 1.0)
        small_value = kwargs.get('small_value', 0.)
        large_value = kwargs.get('large_value', 1.0e99)
        resp_index = kwargs.get('resp_index', 1)
        super().__init__(max_what, frame_structure, tag, recorder_file_storage, haz_req, small_value, large_value)
        self.spring_type = spring_type
        self.normalize_with = normalize_with
        self.resp_index = resp_index

    ####################################################################################################################
    # Public functionalities (instance methods)
    ####################################################################################################################

    def write_generate_edp_recorders(self, for_which, rec_gen_file_open_mode):
        frame_structure = self.structure
        model_work_dir_path = frame_structure.model_work_dir_path
        dir_category = 'NLTHA_Results'
        rec_save_dir_path = self.get_rec_save_dir_path()
        dir_level_to_seek = 'Ground_Motion_Rec_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)
        work_dir_path = Utility.get_path(model_work_dir_path, 'Work_Dir', dir_category, Utility.get_path(*for_which))
        frame_structure.write_spring_deformation_recorder(for_which, self.spring_type, work_dir_path,
                                                          rec_gen_file_open_mode, rec_save_dir_path)
        return

    def get_edp_strings(self, for_which):
        frame_structure = self.structure
        spring_type = self.spring_type
        edp_strs = []
        num_springs = len(frame_structure.get_spring_elem_tags(for_which, spring_type))
        for i_spring in range(1, num_springs + 1):
            spring_id = i_spring
            edp_strs += [f'{spring_type}_{spring_id}']
        return edp_strs

    # ------------------------------------------------------------------------------------------------------------------
