# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 11:14:02 2019

@author: Angshuman Deb
"""

from .structure import Structure
from .utility import Utility
from abc import ABC, abstractmethod
import numpy as np
import os
import scipy.io as sio
import matplotlib.pyplot as plt


class IM(ABC):

    ####################################################################################################################
    # Constructor
    ####################################################################################################################

    def __init__(self, structure, gmm):
        self.structure = structure
        self.gmm = gmm

    ####################################################################################################################
    # Abstract methods
    ####################################################################################################################

    @abstractmethod
    def evaluate_gmm(self, *args, **kwargs):
        pass

    # ------------------------------------------------------------------------------------------------------------------

    @abstractmethod
    def select_ground_motion_records(self, *args, **kwargs):
        pass

    # ------------------------------------------------------------------------------------------------------------------

    @abstractmethod
    def compute_seismic_hazard_integral(self, *args, **kwargs):
        pass

    ####################################################################################################################
    # Public functionalities (instance methods)
    ####################################################################################################################

    def _compute_seismic_hazard_integral(self, prob_dist_im, prob_dist_params_im, im_input):
        im_input = im_input.reshape(len(im_input), 1)  # 2-d vector
        structure = self.structure
        site_hazard_info = structure.get_site_hazard_info()
        scenario_rate = site_hazard_info['scenario_rate']  # 1-d array
        prob_dist_im = prob_dist_im(*prob_dist_params_im)
        nu_im_scenario = (1 - prob_dist_im.cdf(im_input)) * scenario_rate
        nu_im = np.sum(nu_im_scenario, axis=1)  # sum along (over) scenarios

        to_return = {'seismic_hazard_curve': np.column_stack((im_input, nu_im)),
                     'seismic_hazard_curve_deagg': nu_im_scenario}

        return to_return

    # ------------------------------------------------------------------------------------------------------------------

    def get_inv_seismic_hazard(self, mrp, for_which):
        shc, _ = self.get_psha_results(for_which)
        return shc[np.argmin(np.abs(shc[:, 1] - 1 / mrp)), 0]

    # ------------------------------------------------------------------------------------------------------------------

    def get_seismic_hazard_deagg(self, mrp, for_which):
        shc, shc_deagg = self.get_psha_results(for_which)
        ind = np.argmin(np.abs(shc[:, 1] - 1 / mrp))
        # return 1-d array
        # return shc_deagg[ind, :] / shc[ind, 1]
        return (shc_deagg[ind, :] - shc_deagg[ind + 1, :]) / (shc[ind, 1] - shc[ind + 1, 1])

    # ------------------------------------------------------------------------------------------------------------------

    def get_psha_results(self, for_which):
        model_work_dir_path = self.structure.model_work_dir_path
        dir_category = 'PSHA_Results'
        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)
        file_path = Utility.get_path(model_work_dir_path, 'Work_Dir', dir_category,
                                     *Structure.get_modified_for_which(for_which, 'Design_Num'), 'psha_results.pickle')
        run_case_str = ' '.join(for_which)
        psha_results = Utility.pickle_load_dict(file_path)[run_case_str]
        return psha_results['seismic_hazard_curve'], psha_results['seismic_hazard_curve_deagg']

    # ------------------------------------------------------------------------------------------------------------------

    def get_gms_results(self, for_which):
        model_work_dir_path = self.structure.model_work_dir_path
        dir_category = 'GMS_Results'
        dir_level_to_seek = 'Hazard_Level_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)
        file_path = Utility.get_path(model_work_dir_path, 'Work_Dir', dir_category,
                                     *Structure.get_modified_for_which(for_which, 'Design_Num'), 'gms_results.pickle')
        run_case_str = ' '.join(for_which)
        gms_results = Utility.pickle_load_dict(file_path)[run_case_str]

        to_return = list()
        to_return.append(gms_results['ground_motion_records'])
        to_return.append(gms_results['target_spectra'])
        to_return.append(gms_results['mrp'])
        to_return.append(gms_results['n_gm'])

        return to_return

    # ------------------------------------------------------------------------------------------------------------------

    def get_gm_rec_info(self, for_which):
        dir_level_to_seek = 'Ground_Motion_Rec_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)

        gm_records, _, _, _ = self.get_gms_results(for_which)
        gm_index = int(for_which[Structure.get_dir_level_index('Ground_Motion_Rec_Num')]) - 1
        gm_file_names = [gm_records[key][gm_index] for key in gm_records.keys() if 'file_name_' in key]
        scale_fac = gm_records['scale_factor'][gm_index]

        return gm_file_names, scale_fac

    # ------------------------------------------------------------------------------------------------------------------

    def plot_shc(self, for_which, **kwargs):
        figkwargs = kwargs.get('figkwargs', dict())
        grid_alpha = kwargs.get('grid_alpha', 0.5)
        minor_grid_alpha = kwargs.get('minor_grid_alpha', 0.0)
        save_dir_path = kwargs.get('save_dir_path', Utility.get_path(os.getcwd()))
        save_mat = kwargs.get('save_mat', False)
        lc = kwargs.get('lc', 'black')
        lw = kwargs.get('lw', 1.5)
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

        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)

        shc = self.get_psha_results(for_which)[0]
        ax.loglog(shc[:, 0], shc[:, 1], '-', color=lc, linewidth=lw)

        export_mat_dict = dict()
        export_mat_dict['shc'] = shc

        if save_mat:
            structure = self.structure
            name = structure.name
            for_which_str = '_'.join(for_which)
            file_name = f'plot_shc_{name}_{for_which_str}.mat'
            sio.savemat(Utility.get_path(save_dir_path, file_name), export_mat_dict)

        return fig, ax, export_mat_dict

    # ------------------------------------------------------------------------------------------------------------------

    def plot_seismic_hazard_deagg(self, mrp, for_which, **kwargs):
        figkwargs = kwargs.get('figkwargs', dict())
        cmap = kwargs.get('cmap', 'viridis')
        grid_alpha = kwargs.get('grid_alpha', 0.5)
        minor_grid = kwargs.get('minor_grid', False)
        save_dir_path = kwargs.get('save_dir_path', Utility.get_path(os.getcwd()))
        save_mat = kwargs.get('save_mat', False)
        fig_ax = kwargs.get('fig_ax', None)

        structure = self.structure
        site_hazard_info = structure.get_site_hazard_info()
        mag = site_hazard_info['scenario_data']['mag']
        dist = site_hazard_info['scenario_data']['dist']
        rate = self.get_seismic_hazard_deagg(mrp, for_which)

        fig, ax, export_mat_dict = Utility.plot_3d_bar(mag, dist, rate,
                                                       figkwargs=figkwargs, cmap=cmap,
                                                       grid_alpha=grid_alpha, minor_grid=minor_grid,
                                                       fig_ax=fig_ax)

        export_mat_dict['mag'] = export_mat_dict.pop('x')
        export_mat_dict['dist'] = export_mat_dict.pop('y')
        export_mat_dict['rate'] = export_mat_dict.pop('dz')

        if save_mat:
            name = structure.name
            for_which_str = '_'.join(for_which)
            file_name = f'plot_seismic_hazard_deagg_{name}_{for_which_str}_mrp_{mrp}.mat'
            sio.savemat(Utility.get_path(save_dir_path, file_name), export_mat_dict)

        return fig, ax, export_mat_dict

    # ------------------------------------------------------------------------------------------------------------------

    def get_gmm_param_names(self):
        gmm = self.gmm
        params = gmm.PARAMS
        reqd_param_names = [param.name for param in params]
        return reqd_param_names

    # ------------------------------------------------------------------------------------------------------------------

    @staticmethod
    def invert_seismic_hazard(shc, mrp):
        return shc[np.argmin(np.abs(shc[:, 1] - 1 / mrp)), 0]
