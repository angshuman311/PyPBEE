# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 15:37:53 2019

@author: Angshuman Deb
"""

from .utility import Utility
from .multivariate_nataf import multivariate_nataf
from .mixture import mixture
import numpy as np
import os
from scipy.linalg import block_diag
import scipy.io as sio
import matplotlib.pyplot as plt
from collections.abc import Sequence
import matplotlib.tri as tri
import shutil
from copy import deepcopy


class Structure:

    ####################################################################################################################
    # Constructor
    ####################################################################################################################

    def __init__(self, name, location_info, model_files_path, model_work_dir_path, model_params, structural_analysis_platform):
        self.name = name
        self.location_info = location_info
        self.model_files_path = Utility.get_path(model_files_path)
        self.model_work_dir_path = Utility.get_path(model_work_dir_path)
        self.model_params = model_params
        self.structural_analysis_platform = structural_analysis_platform

    ####################################################################################################################
    # Static methods
    ####################################################################################################################

    @staticmethod
    def get_all_dir_levels():
        return ['Analysis_Case', 'Design_Num', 'Model_Option_Num', 'Model_Realization_Num',
                'Hazard_Level_Num', 'Ground_Motion_Rec_Num']

    # ------------------------------------------------------------------------------------------------------------------

    @staticmethod
    def get_dir_level_index(dir_level_to_seek):
        all_dir_levels = Structure.get_all_dir_levels()
        return all_dir_levels.index(dir_level_to_seek)

    # ------------------------------------------------------------------------------------------------------------------

    @staticmethod
    def get_analysis_cases(tag):
        to_return = list()
        if tag == 'Deterministic_FE_Model':
            to_return.extend(['100'])
        if tag == 'Incl_FE_Model_Param_Uncertainty':
            to_return.extend(['200'])
        if tag == 'Deterministic_FE_Model_Incl_Model_Form_Error':
            to_return.extend(['300'])
        if tag == 'Incl_Prob_Dist_Param_Estimation_Uncertainty':
            to_return.extend(['400'])
        if tag == 'All':
            to_return.extend(['100', '200', '300', '400'])

        return to_return

    # ------------------------------------------------------------------------------------------------------------------

    @staticmethod
    def get_analysis_type(case):
        if case == '100':
            return 'Deterministic_FE_Model'
        if case == '200':
            return 'Incl_FE_Model_Param_Uncertainty'
        if case == '300':
            return 'Deterministic_FE_Model_Incl_Model_Form_Error'
        if case == '400':
            return 'Incl_Prob_Dist_Param_Estimation_Uncertainty'

    # ------------------------------------------------------------------------------------------------------------------

    @staticmethod
    def get_modified_for_which(for_which, dir_level_to_seek):
        return for_which[0: Structure.get_dir_level_index(dir_level_to_seek) + 1]

    ####################################################################################################################
    # Public functionalities (instance methods)
    ####################################################################################################################

    def get_design_num_list(self, values):
        if isinstance(values, Sequence) and not isinstance(values, str):
            d = self.model_params['primary_design_params']['value_list_dict']
            design_num_list = [key for key in d.keys() for value in values if d[key] == value]
            return design_num_list
        elif isinstance(values, str) and values == 'all':
            return list(self.model_params['primary_design_params']['value_list_dict'].keys())

    # ------------------------------------------------------------------------------------------------------------------

    def get_design_pts(self, design_num_list):
        d = self.model_params['primary_design_params']['value_list_dict']
        if isinstance(design_num_list, Sequence) and not isinstance(design_num_list, str):
            return [d[num] for num in design_num_list]
        elif isinstance(design_num_list, str) and design_num_list == 'all':
            return [d[key] for key in d.keys()]

    def add_design_point(self, value, gp=False):
        if value is None:
            return
        d = self.model_params['primary_design_params']['value_list_dict']
        if value in d.values():
            return
        max_design_num = np.max(np.array(list(d.keys())).astype(int))
        d[f'{max_design_num + 1}'] = value
        d0 = self.model_params['primary_design_params']['design_point_qualifier_dict']
        if gp:
            d0[f'{max_design_num + 1}'] = ['gp']
        else:
            d0[f'{max_design_num + 1}'] = ['non-gp']

    def set_site_hazard_info(self, im):
        if im.gmm.NAME == 'Boore and Atkinson (2008)':
            site_hazard_info_dir_path = Utility.get_path(self.model_files_path, 'Site_Hazard')

            # Get scenario data
            scenario_data = Utility.read_numpy_array_from_txt_file(
                Utility.get_path(site_hazard_info_dir_path, 'deagg_data_small_sa.txt'), skiprows='find')

            # Eliminate 0 contrib scenarios
            scenario_data = scenario_data[scenario_data[:, -1] != 0, :]

            # Eliminate dist > 300 km scenarios
            scenario_data = scenario_data[scenario_data[:, 0] <= 300, :]

            # scenario_data = [dist mag contrib]
            scenario_data = np.array([scenario_data[:, 0], scenario_data[:, 1], scenario_data[:, -1]]).T

            with open(Utility.get_path(site_hazard_info_dir_path, 'deagg_summary_small_sa.txt'), 'r') as fid:
                curr_line = fid.readline()

            hazard_prob_exceedance = float(curr_line.split(')')[0].split()[-1])
            # TODO: include T in formula
            hazard_mar = -np.log(1 - hazard_prob_exceedance)

            scenario_data[:, -1] = (scenario_data[:, -1] / 100) * hazard_mar

            param_names = im.get_gmm_param_names()

            site_hazard_info = dict()
            site_hazard_info['scenario_rate'] = scenario_data[:, -1]
            site_hazard_info['scenario_data'] = dict()
            site_hazard_info['scenario_data'][param_names[0]] = scenario_data[:, 1]
            site_hazard_info['scenario_data'][param_names[1]] = scenario_data[:, 0]
            site_hazard_info['scenario_data'][param_names[2]] = np.array([self.location_info['v_s30']] * scenario_data.shape[0])
            site_hazard_info['scenario_data'][param_names[3]] = np.array([self.location_info['mechanism']] * scenario_data.shape[0])
            site_hazard_info['scenario_data']['region'] = np.array([self.location_info['region']] * scenario_data.shape[0])

            Utility.pickle_dump_dict(Utility.get_path(site_hazard_info_dir_path, 'site_hazard_info.pickle'),
                                     site_hazard_info)

        return

    # ------------------------------------------------------------------------------------------------------------------

    def get_site_hazard_info(self):
        site_hazard_info_file_path = Utility.get_path(self.model_files_path, 'Site_Hazard',
                                                      'site_hazard_info.pickle')
        return Utility.pickle_load_dict(site_hazard_info_file_path)

    # ------------------------------------------------------------------------------------------------------------------

    def generate_random_model_param_vals(self, analysis_case, design_num, **kwargs):
        rng_seed_recvd = kwargs.get('rng_seed', None)
        save = kwargs.get('save', True)
        sample_size = kwargs.get('sample_size', None)
        sampling_method = kwargs.get('sampling_method', None)

        model_params = self.model_params
        model_work_dir_path = self.model_work_dir_path
        dir_category = 'Random_Model_Param_Vals'

        if sample_size is None:
            sample_size = model_params['random_model_params']['sample_size']
        if sampling_method is None:
            sampling_method = model_params['random_model_params']['sampling_method']

        prob_dist_list = self.model_params['random_model_params']['prob_dist_list']
        corr_matrix = self.model_params['random_model_params']['corr_matrix']
        dist_params_sample_size = model_params['random_model_params']['dist_params_sample_size']
        estimation_sample_size = model_params['random_model_params']['estimation_sample_size']
        diag_blocks = model_params['random_model_params']['diag_blocks']

        if analysis_case in Structure.get_analysis_cases('Deterministic_FE_Model'):
            expected_value_list = self.get_flattened_rv_expected_value_list()
            random_model_param_vals = np.array(expected_value_list)[np.newaxis, :, np.newaxis]

        elif analysis_case in Structure.get_analysis_cases('Incl_FE_Model_Param_Uncertainty'):
            m0 = multivariate_nataf(prob_dist_list, corr_matrix)
            m = multivariate_nataf(self.get_flattened_rv_prob_dist_list(),
                                   self.get_corr_matrix_full(),
                                   corr_matrix_z=self.get_corr_matrix_full(corr_matrix=m0.corr_matrix_z))
            random_model_param_vals = m.rvs(size=sample_size,
                                            random_state=rng_seed_recvd, method=sampling_method)[:, :, np.newaxis]

        elif analysis_case in Structure.get_analysis_cases('Incl_Prob_Dist_Param_Estimation_Uncertainty'):
            dist_param_vals_list = list()
            corr_matrix_realizations = np.array([]).reshape((0, 0, dist_params_sample_size))
            prob_dists = np.array(prob_dist_list)
            itr_group = 0
            for diag_block in diag_blocks:
                prob_dist_block = prob_dists[np.array(diag_block)]
                corr_matrix_block = Utility.get_indexed_matrix(corr_matrix, diag_block)
                if rng_seed_recvd is None:
                    rng_seed = None
                else:
                    rng_seed = Utility.get_num([rng_seed_recvd, itr_group])
                dist_param_vals_list_temp, corr_matrix_realizations_temp = Utility.generate_random_dist_param_vals(
                    prob_dist_block,
                    corr_matrix_block,
                    estimation_sample_size,
                    dist_params_sample_size,
                    rng_seed=rng_seed
                )
                itr_group += 1
                dist_param_vals_list.extend(dist_param_vals_list_temp)
                corr_matrix_realizations = Utility.block_diag_3d(corr_matrix_realizations,
                                                                 corr_matrix_realizations_temp)
            itr_group = 0
            iter_model_option = 0
            prob_dist_list_temp = list()
            for iter_dist in range(len(prob_dist_list)):
                prob_dist_list_temp.append(
                    mixture(dist_param_vals_list[iter_dist], each=prob_dist_list[iter_dist].dist))
            m0 = multivariate_nataf(prob_dist_list_temp, corr_matrix_realizations)
            m = multivariate_nataf(self.get_flattened_rv_prob_dist_list(prob_dist_list=prob_dist_list_temp),
                                   self.get_corr_matrix_full(corr_matrix=m0.corr_matrix),
                                   corr_matrix_z=self.get_corr_matrix_full(corr_matrix=m0.corr_matrix_z))
            if rng_seed_recvd is None:
                rng_seed = None
            else:
                rng_seed = Utility.get_num([rng_seed_recvd, itr_group, iter_model_option])
            random_model_param_vals = m.rvs(size=sample_size, random_state=rng_seed,
                                            method=sampling_method)[:, :, np.newaxis]

        elif analysis_case in Structure.get_analysis_cases('Deterministic_FE_Model_Incl_Model_Form_Error'):
            num_model_options = self.get_num_model_options()
            expected_value_list = self.get_flattened_rv_expected_value_list()
            random_model_param_vals = np.array(expected_value_list)[np.newaxis, :, np.newaxis].repeat(
                num_model_options, axis=2)
        else:
            return

        if save:
            # Save random_model_param_vals
            file_dir_path = Utility.get_path(model_work_dir_path, 'Work_Dir', dir_category,
                                             analysis_case, design_num)
            if not os.path.isdir(file_dir_path):
                os.makedirs(file_dir_path)
            file_path = Utility.get_path(file_dir_path, 'random_model_param_vals.pickle')
            Utility.pickle_dump_dict(file_path, {'random_model_param_vals': random_model_param_vals})
        return random_model_param_vals

    # ------------------------------------------------------------------------------------------------------------------

    def get_random_model_param_vals(self, for_which):
        model_work_dir_path = self.model_work_dir_path
        dir_category = 'Random_Model_Param_Vals'

        dir_level_to_seek = 'Design_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)

        work_dir_path = Utility.get_path(model_work_dir_path, 'Work_Dir', dir_category, Utility.get_path(*for_which))

        file_path = Utility.get_path(work_dir_path, 'random_model_param_vals.pickle')

        return Utility.pickle_load_dict(file_path)['random_model_param_vals']

    # ------------------------------------------------------------------------------------------------------------------

    def get_num_model_options(self):
        model_attributes = self.model_params['model_attributes']
        prod = 1
        for item in list(model_attributes.values()):
            prod *= len(item)
        return prod

    # ------------------------------------------------------------------------------------------------------------------

    def remove_from_work_dir(self, what, **kwargs):

        # what = 'all', 'startswith <tag>', 'contains <tag>', [list of thing(s)]

        local_bash_path = kwargs.get('local_bash_path', '')
        comp_env = kwargs.get('comp_env', 'local')

        if type(what) == str and what == 'all':
            stuff_to_remove = [f'*']

        elif type(what) == str and what.startswith('startswith '):
            tag = what.split()[1]
            stuff_to_remove = [f'{tag}*']

        elif type(what) == str and what.startswith('contains '):
            tag = what.split()[1]
            stuff_to_remove = [f'*{tag}*']

        elif type(what) == list:
            stuff_to_remove = what

        else:
            return

        find_string = ' -o '.join([f'-name "{stuff}"' for stuff in stuff_to_remove])

        model_work_dir_path = self.model_work_dir_path
        work_dir_path = Utility.get_path(model_work_dir_path, 'Work_Dir')

        exec_commands = [f'find "{work_dir_path}" \\( {find_string} \\) -exec echo -n \'"{{}}" \' \\;'
                         f' | xargs -n 50 -P 1 rm -rf']
        Utility.shell_exec('REMOVE', '', exec_commands, comp_env=comp_env, local_bash_path=local_bash_path)

        return

    # ------------------------------------------------------------------------------------------------------------------

    def tar_work_dir(self, **kwargs):

        local_bash_path = kwargs.get('local_bash_path', '')
        exclude_list = kwargs.get('exclude_list', [])
        include_list = kwargs.get('include_list', [])
        if len(exclude_list) != 0 and len(include_list) != 0:
            raise ValueError("Cannot specify both exclude and include lists!")
        comp_env = kwargs.get('comp_env', 'local')
        new_file = kwargs.get('new_file', False)
        flags = kwargs.get('flags', ['c', 'z', 'v', 'f'])

        if new_file:
            tar_file_name = self.name + '_Work_Dir_' + Utility.get_curr_time_stamp() + '.tar'
        else:
            tar_file_name = self.name + '_Work_Dir.tar'

        tar_string = f"tar -{''.join(flags)} {tar_file_name}"

        exclude_string = ''
        for item in exclude_list:
            exclude_string += f' --exclude=\'{item}\''

        tar_string += exclude_string
        if len(include_list) == 0:
            tar_string += ' Work_Dir'
        else:
            for item in include_list:
                tar_string += f' Work_Dir/{item}'

        exec_commands = [tar_string]
        Utility.shell_exec('TAR', self.model_work_dir_path, exec_commands,
                           comp_env=comp_env, local_bash_path=local_bash_path)

        return tar_file_name

    # ------------------------------------------------------------------------------------------------------------------

    def get_material_hysteresis(self, for_which, mat_tag, input_data, num_incr):
        input_data = ' '.join([str(float(data)) for data in input_data])
        model_work_dir_path = self.model_work_dir_path
        dir_category = 'Prelim_Analysis_Results'
        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)
        mat_data_dir_path = Utility.get_path(model_work_dir_path, 'Work_Dir', dir_category, *for_which, self.structural_analysis_platform.model_info_dir_name)
        return self.structural_analysis_platform.perform_mat_test(mat_data_dir_path, mat_tag, input_data, num_incr)

    # ------------------------------------------------------------------------------------------------------------------

    def get_periods(self, for_which, modes):
        model_work_dir_path = self.model_work_dir_path
        dir_category = 'Prelim_Analysis_Results'

        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)

        work_dir_path = Utility.get_path(model_work_dir_path, 'Work_Dir', dir_category, Utility.get_path(*for_which))
        periods = Utility.read_numpy_array_from_txt_file(
            Utility.get_path(work_dir_path, self.structural_analysis_platform.model_info_dir_name, 'periods.txt'),
            skiprows='find').flatten()

        if isinstance(modes, Sequence):
            mode_indices = (np.array(modes) - 1).astype(int)
        else:
            mode_indices = int(modes - 1)
        return periods[mode_indices]

    def get_period_range_extremes(self, for_which, define_range):
        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)
        period_extremes = list()
        for itr in define_range:
            itr_split = itr.split('_')
            if len(itr_split) == 2 and itr_split[0] == 'T' and itr_split[1].isdigit():
                period = self.get_periods(for_which, int(itr_split[1]))
            else:
                raise NotImplementedError(f"period getter not defined for {itr}")
            period_extremes.append(period)
        if len(period_extremes) == 1:
            period_extremes = [period_extremes[0], period_extremes[0]]
        return period_extremes

    def get_period_range(self, for_which, define_range, range_multiplier, num_periods, type_spacing):
        period_extremes = self.get_period_range_extremes(for_which, define_range)
        # return 1-d array
        if type_spacing == 'log':
            return np.logspace(np.log10(range_multiplier[0] * period_extremes[0]),
                               np.log10(range_multiplier[1] * period_extremes[1]), num_periods, endpoint=True)
        else:
            return np.linspace(range_multiplier[0] * period_extremes[0],
                               range_multiplier[1] * period_extremes[1], num_periods, endpoint=True)

    # ------------------------------------------------------------------------------------------------------------------

    def get_flattened_rv_name_list(self):
        rv_list = self.model_params['random_model_params']['rv_list']
        flattened_list = list()
        for itr in range(len(rv_list)):
            flattened_list.extend(rv_list[itr]['name_list'])
        return flattened_list

    # ------------------------------------------------------------------------------------------------------------------

    def get_flattened_rv_expected_value_list(self):
        rv_list = self.model_params['random_model_params']['rv_list']
        flattened_list = list()
        for itr in range(len(rv_list)):
            flattened_list.extend(rv_list[itr]['expected_value_list'] * rv_list[itr]['count'])
        return flattened_list

    # ------------------------------------------------------------------------------------------------------------------

    def get_flattened_rv_prob_dist_list(self, prob_dist_list=None):
        rv_list = self.model_params['random_model_params']['rv_list']
        flattened_list = list()
        if prob_dist_list is None:
            prob_dists = np.array(self.model_params['random_model_params']['prob_dist_list'])
        else:
            prob_dists = np.array(prob_dist_list)
        for itr in range(len(rv_list)):
            prob_dist_ind = np.array(rv_list[itr]['prob_dist_ind_list'])
            temp = list(prob_dists[prob_dist_ind])
            flattened_list.extend(temp * rv_list[itr]['count'])
        return flattened_list

    # ------------------------------------------------------------------------------------------------------------------

    def get_corr_matrix_full(self, corr_matrix=None):
        if corr_matrix is None:
            corr_matrix = self.model_params['random_model_params']['corr_matrix']
        corr_matrix_full = np.array([]).reshape(0, 0)
        rv_list = self.model_params['random_model_params']['rv_list']
        for itr in range(len(rv_list)):
            temp = Utility.get_indexed_matrix(corr_matrix, rv_list[itr]['prob_dist_ind_list'])
            corr_matrix_full = block_diag(corr_matrix_full, *([temp] * rv_list[itr]['count']))
        return corr_matrix_full

    # ------------------------------------------------------------------------------------------------------------------

    def get_rv_index(self, key, whole=False):
        name_list = self.get_flattened_rv_name_list()
        if not whole:
            index_list = [itr for itr in range(len(name_list)) if key in name_list[itr]]
        else:
            index_list = [itr for itr in range(len(name_list)) if key == name_list[itr]]
        return index_list

    # ------------------------------------------------------------------------------------------------------------------

    def get_prob_dist_list_from_flattened(self, rv_indices):
        prob_dist_list_flattened = self.get_flattened_rv_prob_dist_list()
        prob_dist_list = list()
        for ind in rv_indices:
            prob_dist_list.append(prob_dist_list_flattened[ind])
        corr_matrix_full = self.get_corr_matrix_full()
        return prob_dist_list, Utility.get_indexed_matrix(corr_matrix_full, rv_indices)

    # ------------------------------------------------------------------------------------------------------------------

    def plot_scenario_rates(self, **kwargs):
        figkwargs = kwargs.get('figkwargs', dict())
        cmap = kwargs.get('cmap', 'viridis')
        grid_alpha = kwargs.get('grid_alpha', 0.5)
        minor_grid = kwargs.get('minor_grid', False)
        save_dir_path = kwargs.get('save_dir_path', Utility.get_path(os.getcwd()))
        save_mat = kwargs.get('save_mat', False)
        fig_ax = kwargs.get('fig_ax', None)

        site_hazard_info = self.get_site_hazard_info()
        rate = site_hazard_info['scenario_rate']
        mag = site_hazard_info['scenario_data']['mag']
        dist = site_hazard_info['scenario_data'][[item for item in site_hazard_info['scenario_data'].keys() if 'dist' in item][0]]

        fig, ax, export_mat_dict = Utility.plot_3d_bar(mag, dist, rate,
                                                       figkwargs=figkwargs, cmap=cmap,
                                                       grid_alpha=grid_alpha, minor_grid=minor_grid,
                                                       fig_ax=fig_ax)

        export_mat_dict['mag'] = export_mat_dict.pop('x')
        export_mat_dict['dist'] = export_mat_dict.pop('y')
        export_mat_dict['rate'] = export_mat_dict.pop('dz')

        if save_mat:
            name = self.name
            file_name = f'plot_scenario_rates_{name}.mat'
            sio.savemat(Utility.get_path(save_dir_path, file_name), export_mat_dict)

        return fig, ax, export_mat_dict

    # ------------------------------------------------------------------------------------------------------------------

    def plot_rv_distribution_functions(self, rv_index, **kwargs):
        figkwargs = kwargs.get('figkwargs', dict())
        grid_alpha = kwargs.get('grid_alpha', 0.5)
        minor_grid_alpha = kwargs.get('minor_grid_alpha', 0.0)
        save_dir_path = kwargs.get('save_dir_path', Utility.get_path(os.getcwd()))
        save_mat = kwargs.get('save_mat', False)
        k = kwargs.get('k', 5)
        n_cont = kwargs.get('n_cont', 100)
        plot_cdf = kwargs.get('plot_cdf', True)
        color_pdf = kwargs.get('color_pdf', 'blue')
        color_cdf = kwargs.get('color_cdf', 'red')
        lw = kwargs.get('lw', 1)
        scilimits = kwargs.get('scilimits', (-4, 4))
        fig_ax = kwargs.get('fig_ax', None)

        if fig_ax is None or np.array(fig_ax, dtype='object').any() is None:
            fig, ax1 = plt.subplots(**figkwargs)
            ax2 = ax1.twinx()
        else:
            fig = fig_ax[0]
            axs = fig_ax[1]
            ax1 = axs[0]
            if len(axs) == 2:
                ax2 = axs[1]
            else:
                ax2 = ax1.twinx()

        prob_dist = self.get_flattened_rv_prob_dist_list()[rv_index]
        dist_name = prob_dist.dist.name
        mean = prob_dist.mean()
        std = prob_dist.std()
        x1 = mean - k * std
        x2 = mean + k * std
        x = np.linspace(x1, x2, n_cont)
        pdf_x = prob_dist.pdf(x)
        cdf_x = prob_dist.cdf(x)

        ax1.ticklabel_format(axis='both', scilimits=scilimits)
        if plot_cdf:
            ax1.minorticks_on()
            ax1.grid(True, which="major", alpha=grid_alpha, color=color_pdf)
            ax1.grid(True, which="minor", alpha=minor_grid_alpha, color=color_pdf)
            ax1.tick_params(axis='y', labelcolor=color_pdf)
            ax2.minorticks_on()
            ax2.grid(True, which="major", alpha=grid_alpha, color=color_cdf)
            ax2.grid(True, which="minor", alpha=minor_grid_alpha, color=color_cdf)
            ax2.tick_params(axis='y', labelcolor=color_cdf)
            ax2.plot(x, cdf_x, color=color_cdf, linewidth=lw)
            axs = [ax1, ax2]
        else:
            ax1.minorticks_on()
            ax1.grid(True, which="major", alpha=grid_alpha, color='black')
            ax1.grid(True, which="minor", alpha=minor_grid_alpha, color='black')
            ax1.tick_params(axis='y', labelcolor='black')
            axs = [ax1]
            ax2.axis('off')

        ax1.plot(x, pdf_x, color=color_pdf, linewidth=lw)

        export_mat_dict = dict()
        export_mat_dict['x'] = x
        export_mat_dict['pdf_x'] = pdf_x
        export_mat_dict['cdf_x'] = cdf_x
        export_mat_dict['dist_name'] = dist_name

        if save_mat:
            name = self.name
            file_name = f'plot_rv_distribution_functions_{name}_rv_{rv_index}.mat'
            sio.savemat(Utility.get_path(save_dir_path, file_name), export_mat_dict)

        return fig, axs, export_mat_dict

    # ------------------------------------------------------------------------------------------------------------------

    def plot_random_model_param_realizations(self, for_which, rv_indices, **kwargs):
        rng_seed_recvd = kwargs.get('rng_seed', None)
        analysis_case = for_which[Structure.get_dir_level_index('Analysis_Case')]
        sample_size = kwargs.get('sample_size', 1000)
        sampling_method = kwargs.get('sampling_method', 'mcs')
        dist_func = kwargs.get('dist_func', 'pdf')

        prob_dist_list, corr_matrix = self.get_prob_dist_list_from_flattened(rv_indices)
        m = multivariate_nataf(prob_dist_list, corr_matrix)
        data = None

        if 'Deterministic' not in Structure.get_analysis_type(analysis_case):

            if len(for_which) > Structure.get_dir_level_index('Design_Num') + 1:
                dir_level_to_seek = 'Model_Option_Num'
                for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)
                realizations = self.get_random_model_param_vals(for_which)[:, :, int(for_which[-1]) - 1]
                data = realizations[:, rv_indices]

            if analysis_case in Structure.get_analysis_cases('Incl_FE_Model_Param_Uncertainty'):
                if data is None:
                    data = m.rvs(sample_size, rng_seed_recvd, sampling_method)

            elif analysis_case in Structure.get_analysis_cases('Incl_Prob_Dist_Param_Estimation_Uncertainty'):
                itr_group = 0
                if rng_seed_recvd is None:
                    rng_seed = None
                else:
                    rng_seed = Utility.get_num([rng_seed_recvd, itr_group])

                prob_dist_list_temp = list()
                num_rv = len(prob_dist_list)
                dist_param_vals_list, corr_matrix_realztns = Utility.generate_random_dist_param_vals(
                    prob_dist_list,
                    corr_matrix,
                    self.model_params['random_model_params']['estimation_sample_size'],
                    self.model_params['random_model_params']['dist_params_sample_size'],
                    rng_seed=rng_seed
                )
                for iter_num_rv in range(num_rv):
                    prob_dist_list_temp.append(
                        mixture(dist_param_vals_list[iter_num_rv], each=prob_dist_list[iter_num_rv].dist))
                m = multivariate_nataf(prob_dist_list_temp, corr_matrix_realztns)
                prob_dist_list = prob_dist_list_temp

                if data is None:
                    data = m.rvs(sample_size, rng_seed_recvd, sampling_method)
            else:
                return

            if data is None:
                return

            figkwargs = kwargs.get('figkwargs', dict())
            grid_alpha = kwargs.get('grid_alpha', 0.5)
            minor_grid_alpha = kwargs.get('minor_grid_alpha', 0.0)
            save_dir_path = kwargs.get('save_dir_path', Utility.get_path(os.getcwd()))
            save_mat = kwargs.get('save_mat', False)
            k = kwargs.get('k', 5)
            n_cont = kwargs.get('n_cont', 100)
            n_bins = kwargs.get('n_bins', 10)
            ms = kwargs.get('ms', 10)
            lc = kwargs.get('lc', 'black')
            lw = kwargs.get('lw', 1)
            hc = kwargs.get('hc', 'gray')
            mc = kwargs.get('mc', 'red')
            same_lim = kwargs.get('same_lim', False)
            fig_ax = kwargs.get('fig_ax', None)

            if len(rv_indices) == 1:
                if fig_ax is None or np.array(fig_ax, dtype='object').any() is None:
                    fig = plt.figure(**figkwargs)
                    axs = [fig.gca()]
                else:
                    fig = fig_ax[0]
                    axs = fig_ax[1]
                fig, axs[0], export_mat_dict = Utility.scatter_hist_1d(data, prob_dist_list[0],
                                                                       figkwargs=figkwargs,
                                                                       grid_alpha=grid_alpha,
                                                                       minor_grid_alpha=minor_grid_alpha,
                                                                       k=k, n_cont=n_cont,
                                                                       n_bins=n_bins, ms=ms, mc=mc,
                                                                       lw=lw, lc=lc, hc=hc,
                                                                       dist_func=dist_func,
                                                                       fig_ax=[fig, axs[0]])

                export_mat_dict[f'{dist_func}_plot_1'] = export_mat_dict.pop(f'{dist_func}_plot')
                export_mat_dict['dist_name_1'] = export_mat_dict.pop(f'dist_name')

            elif len(rv_indices) == 2:
                fig, axs, export_mat_dict = Utility.scatter_hist_2d(data, m,
                                                                    figkwargs=figkwargs,
                                                                    grid_alpha=grid_alpha,
                                                                    minor_grid_alpha=minor_grid_alpha,
                                                                    k=k, n_cont=n_cont,
                                                                    n_bins=n_bins, ms=ms, mc=mc,
                                                                    lw=lw, lc=lc, hc=hc,
                                                                    same_lim=same_lim, dist_func=dist_func,
                                                                    fig_ax=fig_ax)
                if sampling_method == 'lhs':
                    axs[0].grid(False, which="both")
                    grid, ind = m.get_2d_lhs_grid(sample_size, n_cont)
                    grid = grid[-1]
                    axs[0] = Utility.plot_grid(axs[0], grid[:, :, 0], grid[:, :, 1], ind, lw, lc, grid_alpha)

            elif len(rv_indices) == 3:
                if fig_ax is None or np.array(fig_ax, dtype='object').any() is None:
                    fig = plt.figure(**figkwargs)
                    axs = list()
                    axs.append(fig.add_subplot(2, 2, 1, projection='3d'))
                    axs.append(fig.add_subplot(2, 2, 3))
                    axs.append(fig.add_subplot(2, 2, 4))
                    axs.append(fig.add_subplot(2, 2, 2))
                else:
                    fig = fig_ax[0]
                    axs = fig_ax[1]

                if minor_grid_alpha > 0:
                    axs[0].minorticks_on()
                else:
                    axs[0].minorticks_off()
                axs[0].plot(data[:, 0], data[:, 1], data[:, 2], '.', markersize=ms, color=mc)

                export_mat_dict = dict()
                axes_lim = list()
                for itr in range(len(prob_dist_list)):
                    fig, axs[itr + 1], export_mat_dict = Utility.scatter_hist_1d(data[:, itr], prob_dist_list[itr],
                                                                                 grid_alpha=grid_alpha,
                                                                                 minor_grid_alpha=minor_grid_alpha,
                                                                                 k=k, n_cont=n_cont,
                                                                                 n_bins=n_bins, ms=ms, mc=mc,
                                                                                 lw=lw, lc=lc, hc=hc,
                                                                                 dist_func=dist_func,
                                                                                 fig_ax=[fig, axs[itr + 1]])
                    axes_lim.append(axs[itr + 1].get_xlim())
                    export_mat_dict[f'{dist_func}_plot_{itr + 1}'] = export_mat_dict.pop(f'{dist_func}_plot')
                    export_mat_dict[f'dist_name_{itr + 1}'] = export_mat_dict.pop(f'dist_name')

                axes_lim = np.array(axes_lim)
                if same_lim:
                    axes_lim[0:, :] = np.array([np.min(axes_lim[:, 0]), np.max(axes_lim[:, 1])])
                axs[0].set_xlim(axes_lim[0, :])
                axs[0].set_ylim(axes_lim[1, :])
                axs[0].set_zlim(axes_lim[2, :])
                axs[1].set_xlim(axes_lim[0, :])
                axs[2].set_xlim(axes_lim[1, :])
                axs[3].set_xlim(axes_lim[2, :])

                export_mat_dict[f'data'] = data

            else:
                return

            if save_mat:
                name = self.name
                for_which_str = '_'.join(for_which)
                file_name = f"plot_random_model_param_realizations_{name}_{for_which_str}_" \
                            f"rv_{'_'.join(map(str, rv_indices))}.mat"
                sio.savemat(Utility.get_path(save_dir_path, file_name), export_mat_dict)

            return fig, axs, export_mat_dict

        else:
            return

    # ------------------------------------------------------------------------------------------------------------------

    def get_mixture_prob_dist_list(self, rv_indices, rng_seed=None):
        prob_dist_list, corr_matrix = self.get_prob_dist_list_from_flattened(rv_indices)
        dist_param_vals_list, corr_matrix_realztns = Utility.generate_random_dist_param_vals(
            prob_dist_list,
            corr_matrix,
            self.model_params['random_model_params']['estimation_sample_size'],
            self.model_params['random_model_params']['dist_params_sample_size'],
            rng_seed=rng_seed
        )
        num_rv = len(prob_dist_list)
        prob_dist_list_temp = list()
        for iter_num_rv in range(num_rv):
            prob_dist_list_temp.append(
                mixture(dist_param_vals_list[iter_num_rv], each=prob_dist_list[iter_num_rv].dist))
        return prob_dist_list_temp, multivariate_nataf(prob_dist_list_temp, corr_matrix_realztns)

    # ------------------------------------------------------------------------------------------------------------------

    def create_model_options(self, analysis_case):
        model_attributes = self.model_params['model_attributes']
        if 'Incl_Model_Form_Error' not in Structure.get_analysis_type(analysis_case):
            return [[item[0] for item in list(model_attributes.values())]]
        opts = 1
        opt_inds = list()
        model_attr_vals = list(model_attributes.values())
        for item in model_attr_vals:
            opt_inds.append(list(range(len(item))))
            opts *= len(item)

        cases = Utility.get_multi_level_iterator_cases(*opt_inds).astype(int)
        model_opts = list()
        for opt in range(opts):
            model_opts.append([model_attr_vals[itr][cases[opt, itr]] for itr in range(len(model_attr_vals))])
        return model_opts

    # ------------------------------------------------------------------------------------------------------------------

    def get_model_option_num_list(self, analysis_case, model_attributes):
        model_opts = self.create_model_options(analysis_case)
        inds = np.arange(len(model_opts), dtype=int)
        attr_keys = list(self.model_params['model_attributes'].keys())
        to_ret = list()
        for key in model_attributes.keys():
            for ind in inds:
                if model_opts[ind][attr_keys.index(key)] in model_attributes[key] and (ind + 1) not in to_ret:
                    to_ret.append(ind + 1)
        return to_ret

    # ------------------------------------------------------------------------------------------------------------------

    def fetch_random_model_param_vals(self, for_which, **kwargs):
        random_model_param_vals = kwargs.get('random_model_param_vals', self.get_random_model_param_vals(for_which))
        for_which = Structure.get_modified_for_which(for_which, 'Model_Realization_Num')
        iter_model_realztn = int(for_which[Structure.get_dir_level_index('Model_Realization_Num')])
        iter_model_option = int(for_which[Structure.get_dir_level_index('Model_Option_Num')])
        return random_model_param_vals[iter_model_realztn - 1, :, iter_model_option - 1]

    def write_model_param_vals(self, dir_category, for_which, **kwargs):
        random_model_param_vals = kwargs.get('random_model_param_vals', self.get_random_model_param_vals(for_which))
        model_param_vals = self.fetch_random_model_param_vals(for_which,
                                                              random_model_param_vals=random_model_param_vals)
        for_which = Structure.get_modified_for_which(for_which, 'Model_Realization_Num')
        design_num = for_which[Structure.get_dir_level_index('Design_Num')]
        name_value_pairs = [
            (self.get_flattened_rv_name_list(), model_param_vals),
            (self.model_params['primary_design_params']['name_list'], self.model_params['primary_design_params']['value_list_dict'][design_num]),
            (self.model_params['other_model_params']['name_list'], self.model_params['other_model_params']['value_list']),
        ]
        work_dir_path = Utility.get_path(self.model_work_dir_path, 'Work_Dir', dir_category, *for_which)
        self.structural_analysis_platform.write_model_param_vals(work_dir_path, name_value_pairs)

    # ------------------------------------------------------------------------------------------------------------------

    def fetch_model_attributes(self, for_which, **kwargs):
        model_attributes_values = kwargs.get('model_attributes_values', self.create_model_options(for_which[0]))
        for_which = Structure.get_modified_for_which(for_which, 'Model_Option_Num')
        iter_model_option = int(for_which[Structure.get_dir_level_index('Model_Option_Num')])
        return dict(zip(self.model_params['model_attributes'].keys(), model_attributes_values[iter_model_option - 1]))

    def write_model_attributes(self, dir_category, for_which, **kwargs):
        model_attributes_values = kwargs.get('model_attributes_values', self.create_model_options(for_which[0]))
        name_value_pairs = self.fetch_model_attributes(for_which, model_attributes_values=model_attributes_values)
        name_value_pairs = [
            (list(name_value_pairs.keys()), list(name_value_pairs.values()))
        ]
        for_which = Structure.get_modified_for_which(for_which, 'Model_Realization_Num')
        work_dir_path = Utility.get_path(self.model_work_dir_path, 'Work_Dir', dir_category, *for_which)
        self.structural_analysis_platform.write_model_attributes(work_dir_path, name_value_pairs)

    # ------------------------------------------------------------------------------------------------------------------

    def get_damping_model(self, for_which):
        damping_model_name = self.fetch_model_attributes(for_which)["damping"]
        damping_models = self.model_params['damping_models']
        damping_model = [item for item in damping_models if item['name'] == damping_model_name][0]
        return damping_model

    def get_rayleigh_damping_params(self, for_which, damping_model=None):
        if damping_model is None:
            damping_model = self.get_damping_model(for_which)
        if damping_model['name'] == 'rayleigh_damping':
            xi_i = damping_model['xi_i']
            i = damping_model['i']
            w_i = damping_model['w_i']
            rmpvs = self.fetch_random_model_param_vals(for_which)
            damp_ratio_vals = list()
            for xi in xi_i:
                if type(xi) is str:
                    rv_index = self.get_rv_index(xi, whole=True)[0]
                    damp_ratio_vals.append(rmpvs[rv_index])
                else:
                    damp_ratio_vals.append(xi)
            if i[0] is not None:
                w_i[0] = 2 * np.pi / self.get_periods(for_which, i[0])
            if i[1] is not None:
                w_i[1] = 2 * np.pi / self.get_periods(for_which, i[1])
            alpha, beta = Utility.get_rayleigh_damping_params(damp_ratio_vals[0], w_i[0], damp_ratio_vals[1], w_i[1])
            return alpha, beta
        else:
            raise ValueError(f"Damping model used in for_which = {for_which} is not of type Rayleigh!")

    def get_modal_damping_params(self, for_which, damping_model=None):
        if damping_model is None:
            damping_model = self.get_damping_model(for_which)
        if damping_model['name'] == 'modal_damping':
            xi_i = damping_model['xi_i']
            rmpvs = self.fetch_random_model_param_vals(for_which)
            damp_ratio_vals = list()
            for xi in xi_i:
                if type(xi) is str:
                    rv_index = self.get_rv_index(xi, whole=True)[0]
                    damp_ratio_vals.append(rmpvs[rv_index])
                else:
                    damp_ratio_vals.append(xi)
            return len(xi_i), damp_ratio_vals
        else:
            raise ValueError(f"Damping model used in for_which = {for_which} is not of type Modal!")

    def plot_rayleigh_damping_model(self, for_which, freq, freq_type, **kwargs):
        figkwargs = kwargs.get('figkwargs', dict())
        lc = kwargs.get('lc', 'black')
        lw = kwargs.get('lw', 1.0)
        mc = kwargs.get('mc', 'red')
        grid_alpha = kwargs.get('grid_alpha', 0.5)
        minor_grid_alpha = kwargs.get('minor_grid_alpha', 0.0)
        save_dir_path = kwargs.get('save_dir_path', Utility.get_path(os.getcwd()))
        save_mat = kwargs.get('save_mat', False)
        fig_ax = kwargs.get('fig_ax', None)

        if fig_ax is None or np.array(fig_ax, dtype='object').any() is None:
            fig = plt.figure(**figkwargs)  # also accepts tight_layout=bool
            ax = fig.gca()
        else:
            fig = fig_ax[0]
            ax = fig_ax[1]
        ax.minorticks_on()
        ax.grid(True, which="major", alpha=grid_alpha)
        ax.grid(True, which="minor", alpha=minor_grid_alpha)

        periods = self.get_periods(for_which, [1, 2, 3, 4, 5, 6, 7, 8])
        freq_periods = 2 * np.pi / periods
        if freq_type == 'circular':
            freq_input = 2 * np.pi * np.array(freq)
            freq_periods_plot = 1 / periods
        elif freq_type == 'angular':
            freq_input = np.array(freq)
            freq_periods_plot = freq_periods
        else:
            raise ValueError(f"freq_type should be in ['circular', 'angular']")

        alpha, beta = self.get_rayleigh_damping_params(for_which)
        ax.plot(freq, alpha / 2 / freq_input + beta * freq_input / 2, color=lc, linewidth=lw)
        ax.plot(freq_periods_plot, alpha / 2 / freq_periods + beta * freq_periods / 2, '*', color=mc)

        export_mat_dict = dict()
        export_mat_dict['freq'] = freq
        export_mat_dict['freq_periods_plot'] = freq_periods_plot
        export_mat_dict['freq_type'] = freq_type
        export_mat_dict['alpha'] = alpha
        export_mat_dict['beta'] = beta

        if save_mat:
            name = self.name
            for_which_str = '_'.join(for_which)
            file_name = f'plot_rayleigh_damping_model_{name}_{for_which_str}.mat'
            sio.savemat(Utility.get_path(save_dir_path, file_name), export_mat_dict)

        return fig, ax, export_mat_dict

    def write_damping_file(self, for_which):
        model_work_dir_path = self.model_work_dir_path
        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)
        dir_category = 'Prelim_Analysis_Results'
        work_dir_path = Utility.get_path(model_work_dir_path, 'Work_Dir', dir_category, Utility.get_path(*for_which))
        try:
            damping_model = self.get_damping_model(for_which)
            if 'rayleigh_damping' in damping_model['name']:
                alpha, beta = self.get_rayleigh_damping_params(for_which, damping_model)
                self.structural_analysis_platform.write_damping_file('rayleigh_damping', work_dir_path, alpha=alpha,
                                                                     beta=beta)
            elif 'modal_damping' in damping_model['name']:
                n_modes, damp_ratios = self.get_modal_damping_params(for_which, damping_model)
                self.structural_analysis_platform.write_damping_file('modal_damping', work_dir_path, n_modes=n_modes,
                                                                     damp_ratios=damp_ratios)
            else:
                raise NotImplementedError(f"{damping_model['name']} is NOT an implemented damping model.")
        except KeyError:
            print(f"Damping information not found for_which = {for_which}! Proceeding with empty damping file.")
            self.structural_analysis_platform.write_damping_file(None, work_dir_path)
            pass

    # ------------------------------------------------------------------------------------------------------------------

    def get_m_alpha(self, d_line_list):
        d_vals = np.array(self.get_design_pts(d_line_list))
        line_coeffs = np.polyfit(d_vals[:, 0], d_vals[:, 1], 1)
        m, alpha = line_coeffs[:]
        return m, alpha

    # ------------------------------------------------------------------------------------------------------------------

    def get_d_from_x(self, x, d_line_list):
        scalar = True if np.isscalar(x) else False
        if scalar:
            x = np.array([x])
        else:
            x = np.array(x)
        m, alpha = self.get_m_alpha(d_line_list)
        d = list()
        for itr in range(len(x)):
            d.append((np.linalg.inv(np.array([[1 / m, 1], [-m, 1]])) @ np.array([[x[itr]], [alpha]])).flatten())
        return np.array(d).flatten() if scalar else np.array(d)

    # ------------------------------------------------------------------------------------------------------------------

    def get_x_from_d(self, d, d_line_list):
        d = np.array(d)
        scalar = True if d.ndim == 1 else False
        if scalar:
            d = d[np.newaxis, :]
        m, alpha = self.get_m_alpha(d_line_list)
        x = d[:, 1] + 1 / m * d[:, 0]
        return x.item() if scalar else x

    # ------------------------------------------------------------------------------------------------------------------

    def triangulate_regular_grid_design_space(self):
        model_params = self.model_params
        points = model_params['primary_design_params']['value_list_dict']
        dp_qual = model_params['primary_design_params']['design_point_qualifier_dict']
        points = np.array([points[key] for key in points.keys() if 'non-gp' not in dp_qual[key]])

        x_points = np.sort(np.unique(points[:, 0]))
        y_points = np.sort(np.unique(points[:, 1]))
        x_len = len(x_points) - 1
        y_len = len(y_points) - 1
        grid_points = np.array([[x, y] for x in x_points for y in y_points])
        a = [[i + j * (y_len + 1), (i + 1) + j * (y_len + 1), i + (j + 1) * (y_len + 1)]
             for i in range(y_len) for j in range(x_len)]
        b = [[i + (j + 1) * (y_len + 1), (i + 1) + j * (y_len + 1), (i + 1) + (j + 1) * (y_len + 1)]
             for i in range(y_len) for j in range(x_len)]
        # noinspection PyTypeChecker
        triangulation = tri.Triangulation(grid_points[:, 0],
                                          grid_points[:, 1],
                                          triangles=np.array(a + b))
        return triangulation, grid_points

    # ------------------------------------------------------------------------------------------------------------------

    def plot_design_space(self, **kwargs):
        figkwargs = kwargs.get('figkwargs', dict())
        grid_alpha = kwargs.get('grid_alpha', 0.5)
        minor_grid_alpha = kwargs.get('minor_grid_alpha', 0.0)
        save_dir_path = kwargs.get('save_dir_path', Utility.get_path(os.getcwd()))
        save_mat = kwargs.get('save_mat', False)
        lw = kwargs.get('lw', 1)
        lc = kwargs.get('lc', 'black')
        mc_gp = kwargs.get('mc_gp', 'red')
        ms_gp = kwargs.get('ms_gp', 6)
        mt_gp = kwargs.get('mt_gp', 'o')
        mc_ngp = kwargs.get('mc_ngp', 'red')
        ms_ngp = kwargs.get('ms_ngp', 6)
        mt_ngp = kwargs.get('mt_ngp', 'd')
        mc_ad = kwargs.get('mc_ad', 'red')
        ms_ad = kwargs.get('ms_ad', 18)
        mt_ad = kwargs.get('mt_ad', '*')
        mc_mps = kwargs.get('mc_mps', list())
        ms_mps = kwargs.get('ms_mps', list())
        mt_mps = kwargs.get('mt_mps', list())
        mark_points = kwargs.get('mark_points', list())
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

        tr, grid_data = self.triangulate_regular_grid_design_space()
        ax.set_xticks(np.unique(grid_data[:, 0]))
        ax.set_yticks(np.unique(grid_data[:, 1]))
        dp_qual = self.model_params['primary_design_params']['design_point_qualifier_dict']
        as_designed_key = [key for key in dp_qual.keys() if 'as-designed' in dp_qual[key]]
        non_gp_keys = [key for key in dp_qual.keys() if 'non-gp' in dp_qual[key] and 'as-designed' not in dp_qual[key]]
        as_designed = self.model_params['primary_design_params']['value_list_dict'][as_designed_key[0]]
        non_gp_dps = np.array(
            [self.model_params['primary_design_params']['value_list_dict'][key] for key in non_gp_keys])
        triangles = tr.triangles

        # for itr in range(triangles.shape[0]):
        #     ax.plot(grid_data[np.append(triangles[itr], triangles[itr][0]), 0],
        #             grid_data[np.append(triangles[itr], triangles[itr][0]), 1], '-', lw=lw, color=lc)

        ax.triplot(tr, lw=lw, color=lc)

        ax.plot(grid_data[:, 0], grid_data[:, 1], mt_gp,
                ms=ms_gp, markerfacecolor=mc_gp, markeredgecolor='black', markeredgewidth=0.5)
        ax.plot(*as_designed, mt_ad, ms=ms_ad, markerfacecolor=mc_ad, markeredgecolor='black', markeredgewidth=0.5)
        if non_gp_dps.size != 0:
            ax.plot(non_gp_dps[:, 0], non_gp_dps[:, 1], mt_ngp, ms=ms_ngp, markerfacecolor=mc_ngp,
                    markeredgecolor='black', markeredgewidth=0.5)

        for itr in range(len(mark_points)):
            ax.plot(*mark_points[itr], mt_mps[itr], ms=ms_mps[itr], markerfacecolor=mc_mps[itr],
                    markeredgecolor='black', markeredgewidth=0.5)

        dx = np.mean(np.diff(np.unique(grid_data[:, 0])))
        dy = np.mean(np.diff(np.unique(grid_data[:, 1])))

        ax.set_xlim([np.min(grid_data[:, 0]) - dx, np.max(grid_data[:, 0]) + dx])
        ax.set_ylim([np.min(grid_data[:, 1]) - dy, np.max(grid_data[:, 1]) + dy])

        export_mat_dict = dict()
        export_mat_dict['triangles'] = triangles
        export_mat_dict['grid_data'] = grid_data
        export_mat_dict['as_designed'] = as_designed

        if save_mat:
            file_name = f'plot_design_space_{self.name}.mat'
            sio.savemat(Utility.get_path(save_dir_path, file_name), export_mat_dict)

        return fig, ax, export_mat_dict

    # ------------------------------------------------------------------------------------------------------------------

    def extract_prelim_analysis(self, for_which, extract_dest_path, num_modes):
        model_work_dir_path = self.model_work_dir_path
        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)
        dir_category = 'Prelim_Analysis_Results'
        work_dir_path = Utility.get_path(model_work_dir_path, 'Work_Dir', dir_category, Utility.get_path(*for_which))
        self.write_damping_file(for_which)
        prelim_analysis_files = self.structural_analysis_platform.get_prelim_analysis_files()
        for file in prelim_analysis_files:
            shutil.copy(Utility.get_path(work_dir_path, file), extract_dest_path)
        temp_platform = deepcopy(self.structural_analysis_platform)
        temp_platform.model_files_path = extract_dest_path
        temp_platform.write_run_prelim_analysis(extract_dest_path, num_modes)

    def extract_nltha(self, for_which, extract_dest_path):
        model_work_dir_path = self.model_work_dir_path
        dir_category = 'NLTHA_Results'
        dir_level_to_seek = 'Ground_Motion_Rec_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)
        work_dir_path = Utility.get_path(model_work_dir_path, 'Work_Dir', dir_category, Utility.get_path(*for_which))
        nltha_files = self.structural_analysis_platform.get_nltha_files()
        for file in nltha_files:
            shutil.copy(Utility.get_path(work_dir_path, file), extract_dest_path)
        temp_platform = deepcopy(self.structural_analysis_platform)
        temp_platform.model_files_path = extract_dest_path
        temp_platform.write_run_nltha(
            extract_dest_path,
            extract_dest_path,
            write_model_files=1,
            generate_additional_recorders=True,
            delete_rec=False,
        )


class FrameStructure(Structure):
    ####################################################################################################################
    # Constructor
    ####################################################################################################################

    def __init__(self, name, location_info, model_files_path, model_work_dir_path, model_params, structural_analysis_platform):
        super().__init__(name, location_info, model_files_path, model_work_dir_path, model_params, structural_analysis_platform)

    @staticmethod
    def read_secdef_file(file_path, delim):
        with open(file_path, 'r') as section_file_fid:
            section_file_lines = [[item.strip(', "\'') for item in line.strip().strip(';').split(delim)]
                                  for line in section_file_fid.readlines() if not line.isspace()]
        return section_file_lines

    @staticmethod
    def get_mat_tags_from_secdef(secdef):
        mat_tags = []
        for item in secdef:
            if item[0] in ['patch', 'layer']:
                mat_tags.append(int(item[2]))
            if item[0] == 'fiber':
                mat_tags.append(int(item[-1]))
        return mat_tags

    ####################################################################################################################
    # Public functionalities (instance methods)
    ####################################################################################################################

    def generate_rebar_stress_strain_recorder(self, col_id, edge_id, col_edge_sec_info, rec_save_dir_path, rec_gen_file_fid):
        section_file_lines = col_edge_sec_info['secdef']
        ele_tag = col_edge_sec_info['ele_tag']
        sec_num = col_edge_sec_info['sec_num']
        fib_ctr = 0
        for curr_line in section_file_lines:
            first_word = curr_line[0]

            if first_word == 'fiber':
                # code block
                pass
            if first_word == 'patch':
                # code block
                pass
            if first_word == 'layer':
                second_word = curr_line[1]

                if second_word == 'straight':
                    # code block
                    pass
                if second_word == 'circ':
                    mat_tag = int(curr_line[2])
                    if mat_tag in col_edge_sec_info['rebar_mat_tags']:
                        fib_ctr = self.write_circ_layer_stress_strain_recorder(
                            curr_line, ele_tag, sec_num, col_id, edge_id, fib_ctr, rec_save_dir_path, rec_gen_file_fid)
        return

    # ------------------------------------------------------------------------------------------------------------------

    def write_circ_layer_stress_strain_recorder(self, curr_line, ele_tag, sec_num, col_id, edge_id, fib_ctr,
                                                rec_save_dir_path, rec_gen_file_fid):
        mat_tag = curr_line[2]
        num_fiber = int(curr_line[3])
        y_center = float(curr_line[5])
        z_center = float(curr_line[6])
        radius = float(curr_line[7])
        if len(curr_line) > 8:
            start_ang = float(curr_line[8])
            end_ang = float(curr_line[9])
        else:
            start_ang = 0.0
            end_ang = 360.0 - 360.0 / num_fiber

        ang = start_ang
        ang_incr = (end_ang - start_ang) / (num_fiber - 1)

        # Write recorders for rebars of current layer
        for i_fib in range(num_fiber):
            fib_ctr += 1
            y = y_center + radius * np.cos(np.deg2rad(ang))
            z = z_center + radius * np.sin(np.deg2rad(ang))
            ang += ang_incr
            recorder_file_name = f'fib_stress_strain_col_{col_id}_edge_{edge_id}_elem_{ele_tag}' \
                                 f'_mat_{mat_tag}_fib_{fib_ctr}.txt'
            recorder_file_path = Utility.get_path(rec_save_dir_path, recorder_file_name)
            new_recorder_line = self.structural_analysis_platform.get_stress_strain_recorder_command(
                recorder_file_path, ele_tag, sec_num, y, z, mat_tag)
            rec_gen_file_fid.write(new_recorder_line)
        return fib_ctr

    def write_col_rebar_strain_recorder(self, for_which, work_dir_path,
                                        rec_gen_file_open_mode, rec_save_dir_path):
        rec_gen_file_fid = self.structural_analysis_platform.open_edp_recorder_file(
            work_dir_path,
            rec_gen_file_open_mode,
            rec_save_dir_path
        )
        col_edge_sec_info = self.get_col_edge_sec_info(for_which)
        num_cols_total = self.model_params['num_cols_total']
        col_edge_list = self.model_params['col_edge_list']

        for i_col in range(1, num_cols_total + 1):
            col_id = i_col
            for i_edge in col_edge_list:
                edge_id = i_edge

                self.generate_rebar_stress_strain_recorder(
                    col_id, edge_id,
                    col_edge_sec_info[col_id, edge_id],
                    rec_save_dir_path,
                    rec_gen_file_fid,
                )
        rec_gen_file_fid.close()

    def get_col_edge_sec_info(self, for_which):
        dir_category = 'Prelim_Analysis_Results'
        for_which = Structure.get_modified_for_which(for_which, 'Model_Realization_Num')
        prelim_analysis_work_dir_path = Utility.get_path(self.model_work_dir_path, 'Work_Dir', dir_category, *for_which)
        model_info_files_dir_path = Utility.get_path(prelim_analysis_work_dir_path, self.structural_analysis_platform.model_info_dir_name)

        col_rebar_mat_tags = self.get_col_rebar_mat_tags(for_which)

        col_edge_info = {}
        num_cols_total = self.model_params['num_cols_total']
        col_edge_list = self.model_params['col_edge_list']

        for col_id in range(1, num_cols_total + 1):
            col_elem_sec_info = self.get_col_elem_sec_info(for_which, col_id)
            for edge_id in col_edge_list:
                dict_key = (col_id, edge_id)

                ele_tag, sec_tag, sec_num = None, None, None
                if len(col_elem_sec_info) == 1 and edge_id == 1:
                    ele_tag = col_elem_sec_info[0][0]
                    sec_tag = col_elem_sec_info[0][1]
                    sec_num = 1

                if len(col_elem_sec_info) == 1 and edge_id == 2:
                    ele_tag = col_elem_sec_info[0][0]
                    sec_tag = col_elem_sec_info[0][-1]
                    sec_num = len(col_elem_sec_info[0]) - 1

                if len(col_elem_sec_info) > 1 and edge_id == 1:
                    ele_tag = col_elem_sec_info[0][0]
                    sec_tag = col_elem_sec_info[0][1]
                    sec_num = 1

                if len(col_elem_sec_info) > 1 and edge_id == 2:
                    ele_tag = col_elem_sec_info[-1][0]
                    sec_tag = col_elem_sec_info[-1][-1]
                    sec_num = len(col_elem_sec_info[-1]) - 1

                secdef_file_path = Utility.get_path(model_info_files_dir_path, f'fib_secdef_sec_{sec_tag}.txt')
                section_file_lines_space = FrameStructure.read_secdef_file(secdef_file_path, delim=' ')
                section_file_lines_comma = FrameStructure.read_secdef_file(secdef_file_path, delim=',')
                section_file_lines = None
                for itr in range(len(section_file_lines_comma)):
                    if len(section_file_lines_comma[itr]) > len(section_file_lines_space[itr]):
                        section_file_lines = section_file_lines_comma
                        break
                    elif len(section_file_lines_comma[itr]) < len(section_file_lines_space[itr]):
                        section_file_lines = section_file_lines_space
                        break
                    else:
                        pass
                if section_file_lines is None:
                    # both are same, choose any one
                    section_file_lines = section_file_lines_space

                section_mat_tags = FrameStructure.get_mat_tags_from_secdef(section_file_lines)

                col_edge_info[dict_key] = {
                    'ele_tag': ele_tag,
                    'sec_num': sec_num,
                    'secdef': section_file_lines,
                    'rebar_mat_tags': [mat_tag for mat_tag in col_rebar_mat_tags if mat_tag in section_mat_tags]
                }
        return col_edge_info

    def generate_spring_elem_deformation_recorder(self, ele_tag, spring_type, spring_id,
                                                  rec_save_dir_path, rec_gen_file_fid):
        recorder_file_name = f'spring_deformation_{spring_type}_{spring_id}_elem_{ele_tag}.txt'
        recorder_file_path = Utility.get_path(rec_save_dir_path, recorder_file_name)
        new_recorder_line = self.structural_analysis_platform.get_elem_deformation_recorder_command(recorder_file_path, ele_tag)
        rec_gen_file_fid.write(new_recorder_line)
        return

    def write_spring_deformation_recorder(self, for_which, spring_type, work_dir_path,
                                          rec_gen_file_open_mode, rec_save_dir_path):
        rec_gen_file_fid = self.structural_analysis_platform.open_edp_recorder_file(
            work_dir_path,
            rec_gen_file_open_mode,
            rec_save_dir_path
        )
        spring_elem_tags = self.get_spring_elem_tags(for_which, spring_type)
        for i_spring in range(1, len(spring_elem_tags) + 1):
            self.generate_spring_elem_deformation_recorder(
                spring_elem_tags[i_spring - 1], spring_type, i_spring, rec_save_dir_path, rec_gen_file_fid)
        rec_gen_file_fid.close()

    def get_col_rebar_mat_tags(self, for_which):

        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)

        model_work_dir_path = self.model_work_dir_path
        dir_category = 'Prelim_Analysis_Results'
        file_path = Utility.get_path(model_work_dir_path, 'Work_Dir', dir_category, Utility.get_path(*for_which),
                                     self.structural_analysis_platform.model_info_dir_name, 'col_rebar_mat_info.txt')

        col_rebar_mat_tags = Utility.read_numpy_array_from_txt_file(file_path, skiprows='find')
        return col_rebar_mat_tags.flatten().astype(int)  # 1-d array

    # ------------------------------------------------------------------------------------------------------------------

    def get_spring_elem_tags(self, for_which, spring_type):

        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)

        model_work_dir_path = self.model_work_dir_path
        dir_category = 'Prelim_Analysis_Results'
        file_path = Utility.get_path(model_work_dir_path, 'Work_Dir', dir_category, Utility.get_path(*for_which),
                                     self.structural_analysis_platform.model_info_dir_name, f'{spring_type}_spring_elem_tags.txt')

        spring_elem_tags = Utility.read_numpy_array_from_txt_file(file_path, skiprows='find')
        return spring_elem_tags.flatten().astype(int)  # 1-d array

    # ------------------------------------------------------------------------------------------------------------------

    def get_col_elem_sec_info(self, for_which, col_id):

        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)

        model_work_dir_path = self.model_work_dir_path
        dir_category = 'Prelim_Analysis_Results'
        file_path = Utility.get_path(model_work_dir_path, 'Work_Dir', dir_category, Utility.get_path(*for_which),
                                     self.structural_analysis_platform.model_info_dir_name, f'col_elem_sec_info_col_{col_id}.txt')
        with open(file_path, 'r') as fid:
            lines = [line for line in fid.readlines() if not line.isspace()]
        return [list(np.array(line.strip(';\n ').split()).astype(int)) for line in lines]


class OSB(FrameStructure):

    ####################################################################################################################
    # Constructor
    ####################################################################################################################

    def __init__(self, name, location_info, model_files_path, model_work_dir_path, model_params, structural_analysis_platform):
        super().__init__(name, location_info, model_files_path, model_work_dir_path, model_params, structural_analysis_platform)

    ####################################################################################################################
    # Public functionalities (instance methods)
    ####################################################################################################################

    def get_first_trans_mode_info(self, for_which, **kwargs):
        # NOTE:
        # BridgeCoordsInOpenSees = R*BridgeCoordsInLTV
        # BridgeCoordsInLTV = R'*BridgeCoordsInOpenSees
        # Assume L,T,V as basis
        # 1st row of R = OSB's global X in L,T,V
        # 2nd row of R = OSB's global Y in L,T,V
        # 3rd row of R = OSB's global Z in L,T,V

        crd_transf_matrix = kwargs.get('crd_transf_matrix', np.eye(3))

        model_work_dir_path = self.model_work_dir_path
        dir_category = 'Prelim_Analysis_Results'

        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)

        work_dir_path = Utility.get_path(model_work_dir_path, 'Work_Dir', dir_category, Utility.get_path(*for_which))

        primary_nodes_deck = Utility.read_numpy_array_from_txt_file(
            Utility.get_path(work_dir_path, self.structural_analysis_platform.model_info_dir_name, 'primary_nodes_deck.txt'), skiprows='find').flatten()
        periods = Utility.read_numpy_array_from_txt_file(
            Utility.get_path(work_dir_path, self.structural_analysis_platform.model_info_dir_name, 'periods.txt'),
            skiprows='find').flatten()

        mode_shapes = list()
        for i_mode in range(len(periods)):
            temp = Utility.read_numpy_array_from_txt_file(
                Utility.get_path(work_dir_path, self.structural_analysis_platform.model_info_dir_name, f'mode_shape_{i_mode + 1}.txt'),
                skiprows='find')
            temp[:, 1:4] = (crd_transf_matrix.T @ temp[:, 1:4].T).T
            mode_shapes.append(temp)

        norm_mode = np.zeros((len(periods), 2))

        # Find norm of mode shapes in L and T directions
        for i_mode in range(len(periods)):
            norm_mode[i_mode, 0] = np.linalg.norm(mode_shapes[i_mode][:, 1])
            norm_mode[i_mode, 1] = np.linalg.norm(mode_shapes[i_mode][:, 2])

        # Isolate modes whose T-norm is greater than L-norm
        temp = np.divide(norm_mode[:, 1], norm_mode[:, 0]) > 1.0
        trans_mode_index = np.flatnonzero(temp)  # starting zero

        # Find location of primary nodes of deck in EigenVector matrix
        primary_nodes_deck_loc = np.flatnonzero(np.in1d(mode_shapes[0][:, 0], primary_nodes_deck) == 1)

        # Find T-disp of primary nodes of deck
        disp_trans_primary_nodes_deck = np.zeros((len(primary_nodes_deck_loc), len(trans_mode_index)))
        for i_mode in range(len(trans_mode_index)):
            for i_primary_nodes_deck in range(len(primary_nodes_deck)):
                disp_trans_primary_nodes_deck[i_primary_nodes_deck, i_mode] = mode_shapes[trans_mode_index[i_mode]][
                    primary_nodes_deck_loc[i_primary_nodes_deck]][2]

        # Sum the sign of T-disp of primary nodes of deck for each isolated T-mode
        temp = np.abs(np.sum(np.sign(disp_trans_primary_nodes_deck), axis=0))

        # First T-mode is the mode for which most primary nodes of deck are transversely displaced in the same direction
        # return transModeNum (= trans_mode_index + 1), trans_mode_period
        trans_mode_index = trans_mode_index[np.argmin(np.abs(temp - len(primary_nodes_deck)))]
        trans_mode_period = periods[trans_mode_index]
        return trans_mode_index + 1, trans_mode_period

    # ------------------------------------------------------------------------------------------------------------------

    def get_period_range_extremes(self, for_which, define_range):
        try:
            period_extremes = super().get_period_range_extremes(for_which, define_range)
            return period_extremes
        except NotImplementedError:
            dir_level_to_seek = 'Model_Realization_Num'
            for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)
            period_extremes = list()
            for itr in define_range:
                if itr == 'T_1_trans':
                    _, period = self.get_first_trans_mode_info(for_which)
                else:
                    raise NotImplementedError(f"period getter not defined for {itr}")
                period_extremes.append(period)
            if len(period_extremes) == 1:
                period_extremes = [period_extremes[0], period_extremes[0]]
            return period_extremes

    # ------------------------------------------------------------------------------------------------------------------

    def get_rayleigh_damping_params(self, for_which, damping_model=None):
        try:
            alpha, beta = super().get_rayleigh_damping_params(for_which, damping_model)
            return alpha, beta
        except ValueError:
            if damping_model is None:
                damping_model = self.get_damping_model(for_which)
            if damping_model['name'] == 'rayleigh_damping_first_trans_mode':
                xi_i = damping_model['xi_i'][0]
                rmpvs = self.fetch_random_model_param_vals(for_which)
                if type(xi_i) is str:
                    rv_index = self.get_rv_index(xi_i, whole=True)[0]
                    damp_ratio_val = rmpvs[rv_index]
                else:
                    damp_ratio_val = xi_i
                alpha, beta = self.get_rayleigh_damping_first_trans_mode_params(for_which, damp_ratio_val)
                return alpha, beta
            else:
                raise ValueError(f"Damping model used in for_which = {for_which} is not of type Rayleigh!")

    # ------------------------------------------------------------------------------------------------------------------

    def get_rayleigh_damping_first_trans_mode_params(self, for_which, xi_1):
        _, damp_period = self.get_first_trans_mode_info(for_which)
        xi_2 = 1.10 * xi_1 if xi_1 > 0.05 else 0.05
        w_1 = 2 * np.pi / damp_period
        w_2 = 1.10 * xi_2 * w_1 / xi_1
        alpha, beta = Utility.get_rayleigh_damping_params(xi_1, w_1, xi_2, w_2)
        return alpha, beta

    # ------------------------------------------------------------------------------------------------------------------
