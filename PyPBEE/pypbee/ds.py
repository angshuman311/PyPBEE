# -*- coding: utf-8 -*-
"""
Created on Sat Dec  7 18:32:34 2019

@author: Angshuman Deb
"""

from .structure import Structure
from .utility import Utility
from .mixture import mixture
from .interp_exterp_model import InterpExterpModel
from .im import IM
import numpy as np
import os
import scipy.io as sio
from scipy.stats import lognorm
import sys
from functools import partial
from scipy.interpolate import interp1d
from benedict import benedict
from collections.abc import Sequence
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


class DS:

    ####################################################################################################################
    # Constructor
    ####################################################################################################################

    def __init__(self, edp, predictor, haz_req, ds_type):
        self.structure = edp.structure
        self.tag = edp.tag
        self.edp = edp
        self.predictor = predictor
        self.haz_req = haz_req
        self.ds_type = ds_type

    ####################################################################################################################
    # Setters (instance methods)
    ####################################################################################################################

    def set_haz_req(self, haz_req):
        self.haz_req = haz_req

    ####################################################################################################################
    # Public functionalities
    ####################################################################################################################

    def compute_damage_hazard_integral(self, for_which, im=None, haz_lev_list=None, n_gm_list=None, **kwargs):

        sol_type = kwargs.get('sol_type', 'numerical')

        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)
        run_case_str = ' '.join(for_which)

        edp = self.edp

        to_return = benedict()
        edp_strs = self.edp.get_edp_strings(for_which)

        if sol_type == 'numerical':
            if im is None:
                psdemha_results = edp.get_psdemha_results(Structure.get_modified_for_which(for_which, 'Design_Num'))
                psdemha_results_for_which = psdemha_results[run_case_str]
            elif isinstance(im, IM) and isinstance(haz_lev_list, Sequence):
                delta_input = kwargs.get('delta_input', np.array([]))
                delta_input = np.array(delta_input).flatten()  # 1-d array
                min_max_scale_fac = kwargs.get('min_max_scale_fac', [1, 1])  # list
                psdemha_results_for_which = edp.compute_demand_hazard_integral(for_which, haz_lev_list, im, n_gm_list,
                                                                               delta_input=delta_input,
                                                                               min_max_scale_fac=min_max_scale_fac)
            else:
                raise ValueError("Cannot proceed to damage hazard assessment!")

            for edp_str in edp_strs:
                ds_str = edp_str
                demand_hazard_curve = psdemha_results_for_which[f'{ds_str}']['demand_hazard_curve']
                demand_hazard_curve_deagg_im = psdemha_results_for_which[f'{ds_str}']['demand_hazard_curve_deagg_im']
                temp = self._compute_damage_hazard_integral(for_which, demand_hazard_curve,
                                                            demand_hazard_curve_deagg_im, ds_str)
                nu_ds, nu_ds_edp, nu_ds_im, demand_hazard_curve = temp
                to_return.merge(
                    {f'{ds_str}': {'damage_hazard': nu_ds,
                                   'damage_hazard_deagg_edp': nu_ds_edp,
                                   'damage_hazard_deagg_im': nu_ds_im,
                                   'demand_hazard_curve': demand_hazard_curve
                                   }
                     }
                )

            to_return.merge({'seismic_hazard_curve': psdemha_results_for_which['seismic_hazard_curve']})
            to_return.merge({'mrp_list': psdemha_results_for_which['mrp_list']})

        elif sol_type in ['cf-cornell', 'cf-vamvatsikos']:
            if isinstance(im, IM) and isinstance(haz_lev_list, Sequence):
                mrp_list = list()
                n_gm_list_send = list()
                shc, _ = im.get_psha_results(for_which)
                for haz_lev in haz_lev_list:
                    _, _, mrp, n_gm = im.get_gms_results(for_which + [haz_lev])
                    mrp_list.append(mrp)
                    n_gm_list_send.append(n_gm)
                if isinstance(n_gm_list, Sequence):
                    n_gm_list_send = n_gm_list
                for edp_str in edp_strs:
                    ds_str = edp_str
                    nu_ds = self.get_cf_sol(for_which, sol_type, ds_str, haz_lev_list, mrp_list,
                                            n_gm_list_send, shc)
                    to_return.merge(
                        {f'{ds_str}': {'damage_hazard': np.array([nu_ds])}}
                    )
            else:
                raise ValueError("Cannot proceed to damage hazard assessment!")

        else:
            raise ValueError(f"Solution type = {sol_type} is not implemented! Incorrect solution type!")

        return to_return

    def get_damage_hazard_system(self, for_which, psdamha_results=None):

        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)
        run_case_str = ' '.join(for_which)

        if psdamha_results is None:
            psdamha_results = self.get_psdamha_results(Structure.get_modified_for_which(for_which, 'Design_Num'))
            assert isinstance(psdamha_results, dict)
            psdamha_results_for_which = psdamha_results[run_case_str]
        else:
            psdamha_results_for_which = psdamha_results[run_case_str]

        dam_haz_sys = list()
        edp_strs = self.edp.get_edp_strings(for_which)

        for edp_str in edp_strs:
            dam_haz = psdamha_results_for_which[edp_str]['damage_hazard']
            dam_haz_sys.append(dam_haz)

        dam_haz_sys = np.array(dam_haz_sys)  # 2d array (regions, iter_frag_dist_realztn)
        return np.max(dam_haz_sys, axis=0), list(np.array(edp_strs)[np.argmax(dam_haz_sys, axis=0)])

    def eval_predictor(self, for_which, ds_str):

        predictor = self.predictor
        model_work_dir_path = self.structure.model_work_dir_path

        dir_category = 'Prelim_Analysis_Results'

        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)

        pred_info_file_path = Utility.get_path(
            model_work_dir_path, 'Work_Dir', dir_category, *for_which,
            self.structure.structural_analysis_platform.model_info_dir_name, f'predictor_info_{self.ds_type}_{ds_str}.txt')
        if os.path.isfile(pred_info_file_path):
            x = np.loadtxt(pred_info_file_path).flatten()
        else:
            x = np.array([])
        x = np.hstack((0, x))
        return predictor(x)

    # ------------------------------------------------------------------------------------------------------------------

    def get_denormalized_fragility_dist_params(self, for_which, ds_str):

        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)

        predictor_val = self.eval_predictor(for_which, ds_str)

        normalized_fragility_dist = self.haz_req['normalized_fragility_dist']
        normalized_fragility_dist_param_vals = self.get_normalized_fragility_dist_param_vals(for_which)  # 2-d array

        dist_name = normalized_fragility_dist.dist.name

        if dist_name == 'lognorm':
            normalized_fragility_dist_param_vals[:, 2] *= predictor_val
            return normalized_fragility_dist_param_vals  # 2-d array

    # ------------------------------------------------------------------------------------------------------------------

    def generate_normalized_fragility_dist_param_vals(self, analysis_case, **kwargs):
        rng_seed = kwargs.get('rng_seed', None)

        if analysis_case in Structure.get_analysis_cases('Incl_Prob_Dist_Param_Estimation_Uncertainty'):
            haz_req = self.haz_req
            normalized_fragility_dist = haz_req['normalized_fragility_dist']
            estimation_sample_size = haz_req['estimation_sample_size']
            structure = self.structure
            dist_params_sample_size = structure.model_params['random_model_params']['sample_size']
            dist_param_vals_list, _ = Utility.generate_random_dist_param_vals(
                [normalized_fragility_dist], np.eye(1), estimation_sample_size, dist_params_sample_size,
                rng_seed=rng_seed)
            return dist_param_vals_list[0]  # 2-d array

        else:
            return np.array(Utility.get_all_dist_params(self.haz_req['normalized_fragility_dist']))[np.newaxis, :]

    # ------------------------------------------------------------------------------------------------------------------

    def get_normalized_fragility_dist_param_vals(self, for_which):
        model_work_dir_path = self.structure.model_work_dir_path
        tag = self.tag
        dir_category = 'PSDamHA_Results'

        dir_level_to_seek = 'Design_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)

        file_path = Utility.get_path(model_work_dir_path, 'Work_Dir', dir_category,
                                     *for_which,
                                     'normalized_fragility_dist_param_vals.pickle')

        return Utility.pickle_load_dict(file_path)[tag]  # 2-d array

    # ------------------------------------------------------------------------------------------------------------------

    def _compute_damage_hazard_integral(self, for_which, demand_hazard_curve, demand_hazard_curve_deagg_im, ds_str):
        dem_vals = demand_hazard_curve[:, 0]
        dem_vals_integration = np.exp(np.mean(np.log(np.column_stack([dem_vals[:-1], dem_vals[1:]])), axis=1))
        nu_edp = demand_hazard_curve[:, 1]
        dnu_edp = np.abs(np.diff(nu_edp))

        denormalized_fragility_dist_params = self.get_denormalized_fragility_dist_params(for_which,
                                                                                         ds_str)  # 2-d array
        denormalized_fragility_dist = Utility.get_many_dists(self.haz_req['normalized_fragility_dist'].dist,
                                                             denormalized_fragility_dist_params)

        # nu_ds_edp will be (n_dem_vals, n_dist_params)
        nu_ds_edp = denormalized_fragility_dist.cdf(dem_vals_integration[:, np.newaxis]) * dnu_edp[:, np.newaxis]

        # nu_ds will be (n_dist_params,)
        nu_ds = np.sum(nu_ds_edp, axis=0)

        # nu_ds_im will be (n_dist_params, n_im_vals).T => (n_im_vals, n_dist_params)
        nu_ds_im = np.einsum('ij,ik->jk',
                             denormalized_fragility_dist.cdf(dem_vals_integration[:, np.newaxis]),
                             np.abs(np.diff(demand_hazard_curve_deagg_im, axis=0))).T

        to_return = list()
        to_return.append(nu_ds)
        to_return.append(nu_ds_edp)
        to_return.append(nu_ds_im)
        to_return.append(np.column_stack(
            [dem_vals_integration, np.exp(np.mean(np.log(np.column_stack([nu_edp[:-1], nu_edp[1:]])), axis=1)),
             dnu_edp]))
        return to_return

    # ------------------------------------------------------------------------------------------------------------------

    def get_psdamha_results(self, for_which, ds_str=None):
        model_work_dir_path = self.structure.model_work_dir_path
        dir_category = 'PSDamHA_Results'
        file_path = Utility.get_path(
            model_work_dir_path, 'Work_Dir', dir_category,
            *Structure.get_modified_for_which(for_which, 'Design_Num'), 'psdamha_results.pickle')
        try:
            psdamha_results = Utility.pickle_load_dict(file_path)
        except FileNotFoundError:
            raise FileNotFoundError

        if len(for_which) == (Structure.get_dir_level_index('Design_Num') + 1) or ds_str is None:
            return psdamha_results[self.tag]

        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)
        run_case_str = ' '.join(for_which)
        psdamha_results_for_which = psdamha_results[self.tag][run_case_str]

        to_return = list()
        to_return.append(psdamha_results_for_which[f'{ds_str}']['damage_hazard'])
        to_return.append(psdamha_results_for_which[f'{ds_str}']['damage_hazard_deagg_edp'])
        to_return.append(psdamha_results_for_which[f'{ds_str}']['damage_hazard_deagg_im'])
        to_return.append(psdamha_results_for_which[f'{ds_str}']['demand_hazard_curve'])
        to_return.append(psdamha_results_for_which['seismic_hazard_curve'])
        to_return.append(psdamha_results_for_which['mrp_list'])

        return to_return

    # ------------------------------------------------------------------------------------------------------------------

    def interp_damage_mrp_system_mean_estimate(self, for_which, **kwargs):
        surf_type = kwargs.get('surf_type', 'piecewise_linear')
        d_list = kwargs.get('d_list', [])
        da_list = kwargs.get('da_list', [])

        dir_level_to_seek = 'Design_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)
        analysis_case = for_which[Structure.get_dir_level_index('Analysis_Case')]
        design_num = for_which[Structure.get_dir_level_index('Design_Num')]
        structure = self.structure

        if surf_type not in ['piecewise_linear', 'bilinear_contours']:
            raise ValueError("Invalid surf_type!")

        if surf_type == 'bilinear_contours':
            if len(d_list) < 2:
                raise ValueError("Invalid d_list!")
            if len(da_list) != 2:
                raise ValueError("Invalid da_list!")

        d = np.array(structure.get_design_pts([design_num]))
        if surf_type == 'piecewise_linear':
            mrp_tri_interpolator, _ = self.get_mrp_tri_interpolator_with_data(analysis_case)
            mrp = mrp_tri_interpolator(d[0, 0], d[0, 1]).item()
            return mrp
        elif surf_type == 'bilinear_contours':
            mrp = self.mrp_bilinear_contours_interpolator(d[0, 0], d[0, 1], analysis_case, d_list, da_list)
            return mrp

    # ------------------------------------------------------------------------------------------------------------------

    def get_damage_mrp_system_mean_estimate(self, for_which):
        dir_level_to_seek = 'Design_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)
        analysis_case = for_which[Structure.get_dir_level_index('Analysis_Case')]
        design_num = for_which[Structure.get_dir_level_index('Design_Num')]
        structure = self.structure

        psdamha_results = self.get_psdamha_results(for_which)
        random_model_param_vals = structure.get_random_model_param_vals([analysis_case, design_num])
        # number of realizations for prob dist params of FE model params
        iter_model_option_max = random_model_param_vals.shape[2]
        range_iter_model_option = range(1, iter_model_option_max + 1)
        # number of realizations for FE model params
        iter_model_realztn_max = random_model_param_vals.shape[0]
        range_iter_model_realztn = range(1, iter_model_realztn_max + 1)

        haz_vals_design_num = list()  # will be a 3d array (model_opt, model_realztn, frag_dist_realztn)

        for iter_model_option in range_iter_model_option:
            temp = list()
            for iter_model_realztn in range_iter_model_realztn:
                dam_haz_sys, _ = self.get_damage_hazard_system(
                    [analysis_case, design_num, str(iter_model_option), str(iter_model_realztn)],
                    psdamha_results=psdamha_results
                )
                temp.append(dam_haz_sys)
            haz_vals_design_num.append(temp)

        haz_vals_design_num = np.array(haz_vals_design_num)
        mrp_design_num = 1 / haz_vals_design_num

        return np.mean(mrp_design_num), mrp_design_num

    # ------------------------------------------------------------------------------------------------------------------

    def get_d_star(self, analysis_case, d_list, tgt):
        d_list_mrp_data = self.get_d_list_mrp_data(analysis_case, d_list)
        structure = self.structure
        x = structure.get_x_from_d(d_list_mrp_data[:, :-1], d_list)

        x_star = InterpExterpModel.piecewise_power_interp_continuous_extrap(d_list_mrp_data[:, -1], x, tgt)
        d_star = structure.get_d_from_x(x_star, d_list)

        return d_star, x_star

    # ------------------------------------------------------------------------------------------------------------------

    def get_damage_mrp_interp_surface(self, analysis_case, mesh_size, surf_type, **kwargs):
        d_list = kwargs.get('d_list', [])
        da_list = kwargs.get('da_list', [])

        if surf_type == 'piecewise_linear':
            mrp_tri_interpolator, data = self.get_mrp_tri_interpolator_with_data(analysis_case)
            points = data[:, :-1]
            x = np.linspace(min(points[1:, 0]), max(points[1:, 0]), mesh_size)
            y = np.linspace(min(points[1:, 1]), max(points[1:, 1]), mesh_size)
            xx, yy = np.meshgrid(x, y)
            surface = np.dstack((xx, yy, mrp_tri_interpolator(xx, yy)))
            return data, surface

        elif surf_type == 'bilinear_contours':
            if len(d_list) < 2:
                raise ValueError("Invalid d_list!")
            if len(da_list) != 2:
                raise ValueError("Invalid da_list!")

            structure = self.structure
            d_all = np.array(structure.get_design_pts('all'))
            x = np.linspace(min(d_all[1:, 0]), max(d_all[1:, 0]), mesh_size)
            y = np.linspace(min(d_all[1:, 1]), max(d_all[1:, 1]), mesh_size)
            xx, yy = np.meshgrid(x, y)
            surface = np.dstack((xx, yy, self.mrp_bilinear_contours_interpolator(xx, yy,
                                                                                 analysis_case,
                                                                                 d_list, da_list)))
            d_list_mrp_data = self.get_d_list_mrp_data(analysis_case, d_list)
            da_list_mrp_data = self.get_d_list_mrp_data(analysis_case, da_list)
            data = np.row_stack((d_list_mrp_data, da_list_mrp_data))
            return data, surface

    # ------------------------------------------------------------------------------------------------------------------

    def mrp_bilinear_contours_interpolator(self, xx, yy, analysis_case, d_list, da_list):
        if len(d_list) < 2:
            raise ValueError("Invalid d_list!")
        if len(da_list) != 2:
            raise ValueError("Invalid da_list!")

        if np.isscalar(xx) and np.isscalar(yy):
            scalar = True
            xx = np.array([xx])
            yy = np.array([yy])
        else:
            scalar = False
            xx = np.array(xx)
            yy = np.array(yy)

        zz = np.zeros(xx.shape)

        ppl_interpolator_d, ppl_interpolator_x = \
            self.get_damage_mrp_system_mean_estimate_ppl_interpolator_along_line(analysis_case, d_list)

        structure = self.structure

        m_contours = self.get_m_contours(analysis_case, d_list, da_list)
        da_list_mrp_data = self.get_d_list_mrp_data(analysis_case, da_list)
        m, alpha = structure.get_m_alpha(d_list)
        pos_wrt_line = np.sign(yy - (m * xx) - alpha)
        pos_da_wrt_line = np.sign(da_list_mrp_data[:, 1] - m * da_list_mrp_data[:, 0] - alpha)

        ind_0 = pos_wrt_line == pos_da_wrt_line[0]
        ind_1 = pos_wrt_line == pos_da_wrt_line[1]
        ind_2 = pos_wrt_line == 0

        mat_inv = np.linalg.inv(np.array([[m, -1], [m_contours[0], -1]]))
        temp_d = mat_inv @ np.row_stack(([-alpha] * np.sum(ind_0), m_contours[0] * xx[ind_0] - yy[ind_0]))
        temp_x = structure.get_x_from_d(temp_d.T, d_list)
        zz[ind_0] = ppl_interpolator_x(temp_x)

        mat_inv = np.linalg.inv(np.array([[m, -1], [m_contours[1], -1]]))
        temp_d = mat_inv @ np.row_stack(([-alpha] * np.sum(ind_1), m_contours[1] * xx[ind_1] - yy[ind_1]))
        temp_x = structure.get_x_from_d(temp_d.T, d_list)
        zz[ind_1] = ppl_interpolator_x(temp_x)

        zz[ind_2] = ppl_interpolator_x(structure.get_x_from_d(np.column_stack((xx[ind_2], yy[ind_2])), d_list))

        return zz.item() if scalar else zz

    # ------------------------------------------------------------------------------------------------------------------

    def get_m_contours(self, analysis_case, d_list, da_list, n_cont=1000):
        ppl_interpolator_d, ppl_interpolator_x = \
            self.get_damage_mrp_system_mean_estimate_ppl_interpolator_along_line(analysis_case, d_list)

        d_list_mrp_data = self.get_d_list_mrp_data(analysis_case, d_list)
        da_list_mrp_data = self.get_d_list_mrp_data(analysis_case, da_list)

        structure = self.structure

        d_cont = np.column_stack((np.linspace(np.min(d_list_mrp_data[:, 0]),
                                              np.max(d_list_mrp_data[:, 0]),
                                              n_cont),
                                  np.linspace(np.min(d_list_mrp_data[:, 1]),
                                              np.max(d_list_mrp_data[:, 1]),
                                              n_cont)
                                  ))

        x_cont = structure.get_x_from_d(d_cont, d_list)
        mrp_ppl_values = ppl_interpolator_d(d_cont[:, 0], d_cont[:, 1])
        xa_prime = np.exp(
            interp1d(np.log(mrp_ppl_values), np.log(x_cont), kind='linear', fill_value='extrapolate')(
                np.log(da_list_mrp_data[:, -1])))

        da_prime_list_mrp_data = np.column_stack((structure.get_d_from_x(xa_prime, d_list), da_list_mrp_data[:, -1]))

        m_contours = (da_list_mrp_data[:, 1] - da_prime_list_mrp_data[:, 1]) / \
                     (da_list_mrp_data[:, 0] - da_prime_list_mrp_data[:, 0])

        return m_contours

    # ------------------------------------------------------------------------------------------------------------------

    def get_d_list_mrp_data(self, analysis_case, d_list):
        values = list()
        for design_num in d_list:
            values.append(self.get_damage_mrp_system_mean_estimate([analysis_case, design_num])[0])
        values = np.array(values)
        return np.column_stack((np.array(self.structure.get_design_pts(d_list)), values))

    # ------------------------------------------------------------------------------------------------------------------

    def get_damage_mrp_system_mean_estimate_ppl_interpolator_along_line(self, analysis_case, d_list):
        d_list_mrp_data = self.get_d_list_mrp_data(analysis_case, d_list)
        structure = self.structure
        x = structure.get_x_from_d(d_list_mrp_data[:, :-1], d_list)

        ppl_interpolator = partial(InterpExterpModel.piecewise_power_interp_continuous_extrap,
                                   x, d_list_mrp_data[:, -1])

        def return_func(dp_x, dp_y):
            xq = structure.get_x_from_d(np.squeeze(np.column_stack((dp_x, dp_y))), d_list)
            return ppl_interpolator(xq)

        return return_func, ppl_interpolator

    # ------------------------------------------------------------------------------------------------------------------

    def get_mrp_tri_interpolator_with_data(self, analysis_case):
        structure = self.structure
        model_params = structure.model_params
        points = model_params['primary_design_params']['value_list_dict']
        dp_qual = model_params['primary_design_params']['design_point_qualifier_dict']
        design_num_list = [key for key in points.keys() if 'non-gp' not in dp_qual[key]]
        points = np.array([points[key] for key in points.keys() if 'non-gp' not in dp_qual[key]])
        values = list()
        for design_num in design_num_list:
            values.append(self.get_damage_mrp_system_mean_estimate([analysis_case, design_num])[0])
        triangulation, grid_points = structure.triangulate_regular_grid_design_space()
        z = np.array(values)[[np.where((points == grid_point).all(axis=1))[0].item()
                              for grid_point in grid_points]]
        return tri.LinearTriInterpolator(triangulation, z), np.column_stack((points, values))

    # ------------------------------------------------------------------------------------------------------------------

    def get_cf_sol(self, for_which, sol_type, ds_str, haz_lev_list, mrp_list, n_gm_list, shc):

        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)

        w = np.zeros(shc.shape[0])

        if sol_type == 'cf-cornell':
            ind_low = np.argmin(np.abs(shc[:, 1] - 1 / mrp_list[0]))
            ind_high = np.argmin(np.abs(shc[:, 1] - 1 / mrp_list[-1]))
            w[ind_low:ind_high] = 1.0
            p = np.polyfit(np.log(shc[:, 0]), np.log(shc[:, 1]), 1, w=w)
            k2 = 0
            k1 = -p[0]
            k0 = np.exp(p[1])
        elif sol_type == 'cf-vamvatsikos':
            ind_low = np.argmin(np.abs(shc[:, 1] - 1 / mrp_list[0]))
            ind_high = np.argmin(np.abs(shc[:, 1] - 1 / mrp_list[-1]))
            im_low = 0.7 * shc[ind_low, 0]
            im_high = 1.3 * shc[ind_high, 0]
            ind_low = np.argmin(np.abs(shc[:, 0] - im_low))
            ind_high = np.argmin(np.abs(shc[:, 0] - im_high))

            w[ind_low:ind_high] = 1.0

            p = np.polyfit(np.log(np.array([shc[0, 0], shc[ind_low, 0]])),
                           np.log(np.array([sys.float_info.epsilon, 1])), 1)
            w[:ind_low] = np.exp(p[0] * np.log(shc[:ind_low, 0]) + p[1])

            p = np.polyfit(np.log(np.array([shc[ind_high, 0], shc[-1, 0]])),
                           np.log(np.array([1, sys.float_info.epsilon])), 1)
            w[ind_high:] = np.exp(p[0] * np.log(shc[ind_high:, 0]) + p[1])

            p = np.polyfit(np.log(shc[:, 0]), np.log(shc[:, 1]), 2, w=w)
            k2 = -p[0]
            k1 = -p[1]
            k0 = np.exp(p[2])
        else:
            raise ValueError(f"Solution type = {sol_type} is not implemented! Incorrect solution type!")

        dist_params = np.zeros((len(haz_lev_list), 3))
        edp_vals = self.edp.get_edp_at_hazard_level(for_which, haz_lev_list, n_gm_list, ds_str)
        for itr in range(len(haz_lev_list)):
            dist_params[itr, :] = lognorm.fit(edp_vals[itr], floc=0)

        x_haz_lev = [shc[np.argmin(np.abs(shc[:, 1] - 1 / mrp)), 0] for mrp in mrp_list]

        p = np.polyfit(np.log(x_haz_lev), np.log(dist_params[:, 2]), 1)
        b = p[0]
        a = np.exp(p[1])

        p = np.polyfit(x_haz_lev, dist_params[:, 0], 0)
        zeta_edp = p[0]

        denormalized_fragility_dist_params = self.get_denormalized_fragility_dist_params(for_which,
                                                                                         ds_str)  # 2-d array
        normalized_fragility_dist_name = self.haz_req['normalized_fragility_dist'].dist.name
        if normalized_fragility_dist_name == 'lognorm' and \
                denormalized_fragility_dist_params.shape[0] == 1:
            zeta_c = denormalized_fragility_dist_params[0, 0]
            eta_c = denormalized_fragility_dist_params[0, 2]
        elif normalized_fragility_dist_name == 'lognorm':
            temp = mixture(denormalized_fragility_dist_params, each=self.haz_req['normalized_fragility_dist'].dist)
            zeta_c, _, eta_c = lognorm.fit(temp.rvs(1000000, random_state=999), floc=0)
        else:
            raise ValueError(f"No implementation for normalized fragility dist of type: "
                             f"{normalized_fragility_dist_name}!")

        x_eta_c = (eta_c / a) ** (1 / b)

        q = 1 / (1 + 2 * k2 * (zeta_edp ** 2 / b ** 2))
        phi = 1 / (1 + 2 * k2 * (zeta_edp ** 2 + zeta_c ** 2) / (b ** 2))
        nu_im_x_eta_c = k0 * np.exp(-k2 * ((np.log(x_eta_c)) ** 2) - k1 * np.log(x_eta_c))
        exponent = (1 / (2 * b ** 2)) * q * k1 ** 2 * (zeta_edp ** 2 + phi * (zeta_c ** 2))
        nu_ds = np.sqrt(phi) * (k0 ** (1 - phi)) * (nu_im_x_eta_c ** phi) * np.exp(exponent)

        return nu_ds

    # ------------------------------------------------------------------------------------------------------------------

    def plot_normalized_fragility(self, for_which, **kwargs):
        figkwargs = kwargs.get('figkwargs', dict())
        grid_alpha = kwargs.get('grid_alpha', 0.5)
        minor_grid_alpha = kwargs.get('minor_grid_alpha', 0.0)
        save_dir_path = kwargs.get('save_dir_path', Utility.get_path(os.getcwd()))
        save_mat = kwargs.get('save_mat', False)
        n_cont = kwargs.get('n_cont', 100)
        lc = kwargs.get('lc', 'red')
        lw = kwargs.get('lw', 1)
        k = kwargs.get('k', 5)
        rng_seed = kwargs.get('rng_seed', None)
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

        analysis_case = for_which[Structure.get_dir_level_index('Analysis_Case')]
        if len(for_which) > Structure.get_dir_level_index('Analysis_Case') + 1:
            dir_level_to_seek = 'Design_Num'
            for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)
            normalized_frag_dist_param_vals = self.get_normalized_fragility_dist_param_vals(for_which)
        else:
            normalized_frag_dist_param_vals = self.generate_normalized_fragility_dist_param_vals(analysis_case,
                                                                                                 rng_seed=rng_seed)
        normalized_fragility_dist = self.haz_req['normalized_fragility_dist']

        export_mat_dict = dict()
        export_mat_dict['dist_name'] = normalized_fragility_dist.dist.name
        export_x = list()
        export_cdf_x = list()

        if analysis_case in Structure.get_analysis_cases('Incl_Prob_Dist_Param_Estimation_Uncertainty'):
            normalized_fragility_dist_temp = mixture(normalized_frag_dist_param_vals,
                                                     each=normalized_fragility_dist.dist)
        else:
            normalized_fragility_dist_temp = normalized_fragility_dist.dist(
                *normalized_frag_dist_param_vals.squeeze())
        mean = normalized_fragility_dist_temp.mean()
        std = normalized_fragility_dist_temp.std()
        x1 = mean - k * std
        x2 = mean + k * std
        x = np.linspace(x1, x2, n_cont)
        cdf_x = normalized_fragility_dist_temp.cdf(x)
        ax.plot(x, cdf_x, linewidth=lw, color=lc)
        export_x.append(x)
        export_cdf_x.append(cdf_x)

        export_mat_dict[f'x'] = np.column_stack(export_x)
        export_mat_dict[f'cdf_x'] = np.column_stack(export_cdf_x)

        if save_mat:
            structure = self.structure
            name = structure.name
            tag = self.tag
            file_name = f'plot_normalized_fragility_{name}_{analysis_case}_ds_{tag}.mat'
            sio.savemat(Utility.get_path(save_dir_path, file_name), export_mat_dict)

        return fig, ax, export_mat_dict

    # ------------------------------------------------------------------------------------------------------------------

    def plot_damage_mrp_surface(self, analysis_case, **kwargs):
        figkwargs = kwargs.get('figkwargs', dict())
        cmap = kwargs.get('cmap', 'viridis')
        minor_grid = kwargs.get('minor_grid', False)
        save_dir_path = kwargs.get('save_dir_path', Utility.get_path(os.getcwd()))
        save_mat = kwargs.get('save_mat', False)
        mesh_size = kwargs.get('mesh_size', 100)
        scilimits = kwargs.get('scilimits', (-4, 4))
        fig_ax = kwargs.get('fig_ax', None)
        surf_type = kwargs.get('surf_type', 'piecewise_linear')
        d_list = kwargs.get('d_list', [])
        da_list = kwargs.get('da_list', [])

        if surf_type not in ['piecewise_linear', 'bilinear_contours']:
            raise ValueError("Invalid surf_type!")

        if surf_type == 'bilinear_contours':
            if len(d_list) < 2:
                raise ValueError("Invalid d_list!")
            if len(da_list) != 2:
                raise ValueError("Invalid da_list!")

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
        ax.ticklabel_format(axis='both', scilimits=scilimits)

        data, surface = self.get_damage_mrp_interp_surface(analysis_case, mesh_size, surf_type,
                                                           d_list=d_list, da_list=da_list)

        ax.plot_surface(surface[:, :, 0], surface[:, :, 1], surface[:, :, -1],
                        linewidth=0, cmap=cmap, edgecolor='none', antialiased=False)

        ax_in = inset_axes(ax, width="50%", height="5%", loc='upper left')

        sm = cm.ScalarMappable(cmap=cmap)
        sm.set_array(surface[:, :, -1].flatten())
        cbar = fig.colorbar(sm, cax=ax_in, orientation='horizontal')
        cbar.ax.ticklabel_format(axis='both', scilimits=scilimits)

        export_mat_dict = dict()
        export_mat_dict['data'] = data
        export_mat_dict['surface'] = surface

        if save_mat:
            structure = self.structure
            name = structure.name
            tag = self.tag
            file_name = f'plot_damage_mrp_surface_{name}_{analysis_case}_ds_{tag}.mat'
            sio.savemat(Utility.get_path(save_dir_path, file_name), export_mat_dict)

        return fig, ax, export_mat_dict

    # ------------------------------------------------------------------------------------------------------------------

    def plot_denormalized_fragility(self, for_which, ds_str, **kwargs):
        figkwargs = kwargs.get('figkwargs', dict())
        grid_alpha = kwargs.get('grid_alpha', 0.5)
        minor_grid_alpha = kwargs.get('minor_grid_alpha', 0.0)
        save_dir_path = kwargs.get('save_dir_path', Utility.get_path(os.getcwd()))
        save_mat = kwargs.get('save_mat', False)
        n_cont = kwargs.get('n_cont', 100)
        lc = kwargs.get('lc', 'red')
        lw = kwargs.get('lw', 1)
        k = kwargs.get('k', 5)
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
        denormalized_fragility_dist_params = self.get_denormalized_fragility_dist_params(for_which,
                                                                                         ds_str)  # 2-d array

        normalized_fragility_dist_name = self.haz_req['normalized_fragility_dist'].dist.name
        if normalized_fragility_dist_name == 'lognorm':
            denormalized_fragility_dist = mixture(denormalized_fragility_dist_params,
                                                  each=self.haz_req['normalized_fragility_dist'].dist)
        else:
            raise ValueError(f"No implementation for normalized fragility dist of type: "
                             f"{normalized_fragility_dist_name}!")

        mean = denormalized_fragility_dist.mean()
        std = denormalized_fragility_dist.std()
        x1 = mean - k * std
        x2 = mean + k * std
        x = np.linspace(x1, x2, n_cont)
        cdf_x = denormalized_fragility_dist.cdf(x)
        ax.plot(x, cdf_x, linewidth=lw, color=lc)

        export_mat_dict = dict()
        export_mat_dict[f'x'] = x
        export_mat_dict[f'cdf_x'] = cdf_x
        export_mat_dict['dist_name'] = denormalized_fragility_dist.dist.name

        if save_mat:
            structure = self.structure
            name = structure.name
            tag = self.tag
            for_which_str = '_'.join(for_which)
            file_name = f'plot_denormalized_fragility_{name}_{for_which_str}_' \
                        f'ds_{tag}_{ds_str}.mat'
            sio.savemat(Utility.get_path(save_dir_path, file_name), export_mat_dict)

        return fig, ax, export_mat_dict

    # ------------------------------------------------------------------------------------------------------------------

    def plot_haz_deagg(self, for_which, wrt, **kwargs):
        figkwargs = kwargs.get('figkwargs', dict())
        grid_alpha = kwargs.get('grid_alpha', 0.5)
        minor_grid_alpha = kwargs.get('minor_grid_alpha', 0.0)
        save_dir_path = kwargs.get('save_dir_path', Utility.get_path(os.getcwd()))
        save_mat = kwargs.get('save_mat', False)
        lc1 = kwargs.get('lc1', 'blue')
        lc2 = kwargs.get('lc2', 'red')
        lw1 = kwargs.get('lw1', 1)
        lw2 = kwargs.get('lw2', 1)
        txc = kwargs.get('txc', 'black')
        patch_color = kwargs.get('patch_color', 'cyan')
        patch_alpha = kwargs.get('patch_alpha', 0.25)
        im_lim = kwargs.get('im_lim', [])
        edp_lim = kwargs.get('edp_lim', [])
        fig_ax = kwargs.get('fig_ax', None)

        if fig_ax is None or np.array(fig_ax, dtype='object').any() is None:
            fig, ax1 = plt.subplots(**figkwargs)
            ax2 = ax1.twinx()
        else:
            fig = fig_ax[0]
            axs = fig_ax[1]
            ax1 = axs[0]
            ax2 = axs[1]
        color1 = lc1
        color2 = lc2
        ax1.minorticks_on()
        ax1.grid(True, which="major", alpha=grid_alpha, color=color1)
        ax1.grid(True, which="minor", alpha=minor_grid_alpha, color=color1)
        ax2.minorticks_on()
        ax2.grid(True, which="major", alpha=grid_alpha, color=color2)
        ax2.grid(True, which="minor", alpha=minor_grid_alpha, color=color2)

        dir_level_to_seek = 'Model_Realization_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)

        dam_haz, ds_str = self.get_damage_hazard_system(for_which)

        ind_max = np.argmax(dam_haz)
        ds_str = ds_str[ind_max]
        dam_haz = dam_haz[ind_max]
        _, dam_haz_deagg_edp, dam_haz_deagg_im, dhc, shc, mrp_list = self.get_psdamha_results(for_which, ds_str)
        dam_haz_deagg_edp = dam_haz_deagg_edp[:, ind_max]
        dam_haz_deagg_im = dam_haz_deagg_im[:, ind_max]

        export_mat_dict = dict()

        if wrt == 'edp':
            ax1.semilogy(dhc[:, 0], dhc[:, 1], color=color1, linewidth=lw1)
            ax2.plot(dhc[1:, 0], dam_haz_deagg_edp[1:] / dam_haz / np.diff(dhc[:, 0]),
                     color=color2, linewidth=lw2)
            ax1.tick_params(axis='y', labelcolor=color1)
            ax2.tick_params(axis='y', labelcolor=color2)
            if len(edp_lim) == 0:
                edp_lim = ax1.get_xlim()
            ax1.set_xlim(edp_lim)
            export_mat_dict['dhc'] = dhc
            export_mat_dict['edp_deagg'] = np.column_stack([dhc[1:, 0],
                                                            dam_haz_deagg_edp[1:] / dam_haz / np.diff(dhc[:, 0])])

        elif wrt == 'im':
            ax1.semilogy(shc[:, 0], shc[:, 1], color=color1, linewidth=lw1)
            ax2.plot(shc[1:, 0], dam_haz_deagg_im[1:] / dam_haz / np.diff(shc[:, 0]),
                     color=color2, linewidth=lw2)
            ax1.tick_params(axis='y', labelcolor=color1)
            ax2.tick_params(axis='y', labelcolor=color2)
            if len(im_lim) == 0:
                im_lim = ax1.get_xlim()
            ax1.set_xlim(im_lim)
            ylim = ax1.get_ylim()
            mrp_min = min(mrp_list)
            mrp_max = max(mrp_list)
            im_min = shc[np.argmin(np.abs(shc[:, 1] - 1 / mrp_min)), 0]
            im_max = shc[np.argmin(np.abs(shc[:, 1] - 1 / mrp_max)), 0]
            rect = plt.Rectangle(
                (im_min, ylim[0]), im_max - im_min, ylim[1] - ylim[0], color=patch_color,
                alpha=patch_alpha)
            ax1.add_patch(rect)
            ax1.set_ylim(ylim)

            im_same_haz = shc[np.argmin(np.abs(shc[:, 1] - dam_haz)), 0]
            ax1.semilogy([im_lim[0], im_same_haz, im_same_haz], [dam_haz, dam_haz, ylim[0]],
                         color=txc, linestyle='-.', linewidth=lw1 / 2)
            ax1.text(im_lim[0], dam_haz, f'1/{int(np.floor(1 / dam_haz))}', color=txc, horizontalalignment='left',
                     verticalalignment='bottom')

            export_mat_dict['shc'] = shc
            export_mat_dict['im_deagg'] = np.column_stack([shc[1:, 0],
                                                           dam_haz_deagg_im[1:] / dam_haz / np.diff(shc[:, 0])])
            export_mat_dict['im_patch_lim'] = [im_min, im_max]
            export_mat_dict['dam_haz'] = dam_haz
            export_mat_dict['im_same_haz'] = im_same_haz

        if save_mat:
            structure = self.structure
            name = structure.name
            tag = self.tag
            for_which_str = '_'.join(for_which)
            file_name = f'plot_haz_deagg_wrt_{wrt}_{name}_{for_which_str}_' \
                        f'ds_{tag}.mat'
            sio.savemat(Utility.get_path(save_dir_path, file_name), export_mat_dict)

        axs = [ax1, ax2]
        return fig, axs, export_mat_dict

    # ------------------------------------------------------------------------------------------------------------------

    def plot_mrp_histogram(self, for_which, mrp_det, **kwargs):
        figkwargs = kwargs.get('figkwargs', dict())
        grid_alpha = kwargs.get('grid_alpha', 0.5)
        minor_grid_alpha = kwargs.get('minor_grid_alpha', 0.0)
        save_dir_path = kwargs.get('save_dir_path', Utility.get_path(os.getcwd()))
        save_mat = kwargs.get('save_mat', False)
        k = kwargs.get('k', 5)
        n_bins = kwargs.get('n_bins', 10)
        ms = kwargs.get('ms', 10)
        data_ms = kwargs.get('data_ms', 5)
        lc = kwargs.get('lc', 'black')
        lw = kwargs.get('lw', 1)
        hc = kwargs.get('hc', 'gray')
        mc = kwargs.get('mc', 'red')
        fs = kwargs.get('fs', 16)
        scilimits = kwargs.get('scilimits', (-4, 4))
        # leg_loc = kwargs.get('leg_loc', 'best')
        hatch_d_tgt = kwargs.get('hatch_d_tgt', 5)
        hatch_c_tgt = kwargs.get('hatch_c_tgt', lc)
        # lc_tgt = kwargs.get('lc_tgt', 'black')
        hatch_d_det = kwargs.get('hatch_d_det', 3)
        hatch_c_det = kwargs.get('hatch_c_det', lc)
        # lc_det = kwargs.get('lc_det', 'black')
        # ls_det = kwargs.get('ls_det', '--')
        # ls_tgt = kwargs.get('ls_tgt', '-.')
        tgt_mrp = kwargs.get('tgt_mrp', None)
        fig_ax = kwargs.get('fig_ax', None)
        model_opts_exclude = kwargs.get('model_opts_exclude', None)

        if fig_ax is None or np.array(fig_ax, dtype='object').any() is None:
            fig, ax1 = plt.subplots(**figkwargs)
            ax2 = ax1.twinx()
        else:
            fig = fig_ax[0]
            axs = fig_ax[1]
            ax1 = axs[0]
            ax2 = axs[1]

        ax1.ticklabel_format(axis='both', scilimits=scilimits)
        ax2.ticklabel_format(axis='both', scilimits=scilimits)

        color1 = 'black'
        color2 = lc
        ax1.minorticks_on()
        ax1.grid(True, which="major", alpha=grid_alpha, color=color1)
        ax1.grid(True, which="minor", alpha=minor_grid_alpha, color=color1)
        ax2.minorticks_on()
        ax2.grid(True, which="major", alpha=grid_alpha, color=color2)
        ax2.grid(True, which="minor", alpha=minor_grid_alpha, color=color2)

        analysis_case = for_which[Structure.get_dir_level_index('Analysis_Case')]

        allowed_cases = Structure.get_analysis_cases('All')
        allowed_cases.remove(Structure.get_analysis_cases('Deterministic_FE_Model')[0])

        _, mrp_vals = self.get_damage_mrp_system_mean_estimate(Structure.get_modified_for_which(for_which, 'Design_Num'))
        iter_model_option_max, iter_model_realztn_max, iter_frag_dist_param_realztn_max = mrp_vals.shape

        try:
            model_option_ind = int(for_which[Structure.get_dir_level_index('Model_Option_Num')]) - 1
            include_mask = np.zeros((iter_model_option_max,), dtype=bool)
            include_mask[model_option_ind] = True
        except IndexError:
            include_mask = np.ones((iter_model_option_max,), dtype=bool)
            if model_opts_exclude is not None and len(model_opts_exclude) != 0:
                model_opts_exclude_ind = np.array(model_opts_exclude) - 1
                include_mask[model_opts_exclude_ind] = False
        data = mrp_vals[include_mask, :, :]

        if analysis_case in allowed_cases:
            fig, ax1, export_mat_dict = Utility.scatter_hist_1d(data.ravel(), grid_alpha=grid_alpha,
                                                                minor_grid_alpha=minor_grid_alpha,
                                                                n_bins=n_bins, ms=data_ms, mc=mc,
                                                                lw=lw, lc=lc, hc=hc,
                                                                fig_ax=[fig, ax1])
            data_mean = np.mean(data.ravel())

            term_1 = np.mean(np.var(data, axis=2).ravel())
            term_2 = np.mean(np.var(np.mean(data, axis=2), axis=1))
            term_3 = np.var(np.mean(data, axis=(1, 2)))
            tot_var = term_1 + term_2 + term_3
            data_std = np.sqrt(tot_var)
            # data_std = np.std(data)

            ax1.plot(data_mean, 0.0, 'd', color=mc, markersize=ms, zorder=10, clip_on=False,
                     markeredgecolor='black', markeredgewidth=0.5,
                     # label=f'E[MRP]: {int(np.floor(data_mean))} yrs',
                     )

            ax1.plot(mrp_det, 0.0, 'o', color=mc, markersize=ms, zorder=10, clip_on=False,
                     markeredgecolor='black', markeredgewidth=0.5,
                     # label=r'$\mathrm{MRP^{Det}: }$' + f'{int(np.floor(mrp_det))}' + ' yrs',
                     )

            export_mat_dict['mrp_w/'] = data_mean
            export_mat_dict['mrp_w/o'] = mrp_det

            ax1.set_xlim(
                (np.max([0, data_mean - k * data_std]),
                 np.max([1.1 * np.max(data.ravel()), data_mean + k * data_std])))

            x1 = ax1.get_xlim()[0]
            if tgt_mrp is not None:
                x2 = tgt_mrp
                ax1.axvspan(x1, x2, clip_path=export_mat_dict['patches'][0], fill=False, hatch='/' * hatch_d_tgt,
                            color=hatch_c_tgt, lw=lw / 2)
                ax1.set_xlim(
                    (np.max([0, data_mean - k * data_std]),
                     np.max([1.1 * np.max(data.ravel()), data_mean + k * data_std, 1.1 * tgt_mrp])))
                prob_le_tgt = np.sum(data.ravel() < tgt_mrp) / data.size
                export_mat_dict['prob_mrp_le_tgt'] = prob_le_tgt
                # plotting probability lines
                # ax2.plot(([tgt_mrp] * 2) + [ax1.get_xlim()[-1]], [0, prob_le_tgt, prob_le_tgt],
                #          linestyle=ls_tgt, color=color2, lw=lw, zorder=10)
                # ax2.plot(tgt_mrp, prob_le_tgt, 'x',
                #          markeredgecolor=color2, markeredgewidth=lw, markerfacecolor="None",
                #          ms=ms, zorder=10)

            x2 = mrp_det
            ax1.axvspan(x1, x2, clip_path=export_mat_dict['patches'][0], fill=False, hatch="-" * hatch_d_det,
                        color=hatch_c_det, lw=lw / 2)
            prob_le_det = np.sum(data.ravel() < mrp_det) / data.size
            export_mat_dict['prob_mrp_le_det'] = prob_le_det
            # plotting probability lines
            # ax2.plot(([mrp_det] * 2) + [ax1.get_xlim()[-1]], [0, prob_le_det, prob_le_det],
            #          linestyle=ls_det, color=color2, lw=lw, zorder=10)
            # ax2.plot(mrp_det, prob_le_det, 'x',
            #          markeredgecolor=color2, markeredgewidth=lw, markerfacecolor="None",
            #          ms=ms, zorder=10)

            export_mat_dict['mrp_cov'] = data_std / np.abs(data_mean)
            export_mat_dict['mrp_ratio'] = data_mean / mrp_det

            ecdf_x, ecdf_y = Utility.ecdf(data.ravel())
            ecdf_x = np.hstack((ax2.get_xlim()[0], ecdf_x, ax2.get_xlim()[-1]))
            ecdf_y = np.hstack((0, ecdf_y, 1))

            ax2.step(ecdf_x, ecdf_y, where='post', color=color2, lw=lw)
            ax2.tick_params(axis='y', labelcolor=color2)
            ylim1 = ax1.get_ylim()
            ylim2 = ax2.get_ylim()

            ax1.set_ylim((-(ylim1[-1] - np.max(export_mat_dict['bin_hgts'])), ylim1[-1]))
            ax2.set_ylim(ylim2)

            if tgt_mrp is not None:
                ax2.vlines(tgt_mrp, *ax2.get_ylim(), colors='black', linestyles='-', linewidths=3, zorder=10)
                # ax2.annotate('Target ',
                #              xy=[tgt_mrp, ylim2[-1]],
                #              xytext=[tgt_mrp * 0.7, ylim2[-1] * 1.1],
                #              fontsize=fs,
                #              arrowprops=dict(fc='black',
                #                              ec=None,
                #                              arrowstyle='-|>'),
                #              va='bottom', ha='right')
            ax2.vlines(mrp_det, *ax2.get_ylim(), colors='black', linestyles='-', linewidths=3, zorder=10)
            # ax2.annotate(' Deterministic',
            #              xy=[mrp_det, ylim2[-1]],
            #              xytext=[mrp_det * 1.3, ylim2[-1] * 1.1],
            #              fontsize=fs,
            #              arrowprops=dict(fc='black',
            #                              ec=None,
            #                              arrowstyle='-|>'),
            #              va='bottom', ha='left')

            textstr = r'$\mathrm{\frac{E[MRP]}{MRP^{Det}}}$' + f'≈ {data_mean / mrp_det:0.2}  '
            xlim = ax2.get_xlim()
            ax2.text(xlim[-1],
                     0.5 * (ylim2[0] + ylim2[-1]),
                     textstr,
                     fontsize=fs,
                     ha='right',
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=1.0),
                     zorder=100,
                     )
            # handles, labels = ax1.get_legend_handles_labels()
            # # legend = plt.legend(handles, labels, loc='upper right')
            # legend = plt.legend(handles, labels, loc=leg_loc, fontsize=fs)
            # ax2.add_artist(legend)
            # legend.set_zorder(100)

            if tgt_mrp is not None:
                nu_data = 1 / data.ravel()

                pe_tgt = np.mean(1 - np.exp(-nu_data * tgt_mrp))
                pe_tgt_det = 1 - np.exp(-(1 / mrp_det) * tgt_mrp)

                pe_50 = np.mean(1 - np.exp(-nu_data * 50))
                pe_50_det = np.mean(1 - np.exp(-(1 / mrp_det) * 50))

                pe_annual = np.mean(1 - np.exp(-nu_data * 1))
                pe_annual_det = np.mean(1 - np.exp(-(1 / mrp_det) * 1))

                export_mat_dict['prob_lse_in_tgt_exp'] = pe_tgt
                export_mat_dict['prob_lse_in_tgt_exp_det'] = pe_tgt_det

                export_mat_dict['prob_lse_in_50_years'] = pe_50
                export_mat_dict['prob_lse_in_50_years_det'] = pe_50_det

                export_mat_dict['prob_lse_annual'] = pe_annual
                export_mat_dict['prob_lse_annual_det'] = pe_annual_det

                export_mat_dict['pc_change_tgt_exp'] = (pe_tgt - pe_tgt_det) / pe_tgt_det * 100
                export_mat_dict['pc_change_50'] = (pe_50 - pe_50_det) / pe_50_det * 100
                export_mat_dict['pc_change_annual'] = (pe_annual - pe_annual_det) / pe_annual_det * 100

            export_mat_dict['mrp_mean_estimate_w_uncertainty'] = data_mean
            export_mat_dict['mrp_mean_estimate_wo_uncertainty'] = mrp_det
            export_mat_dict['mrp_ecdf'] = np.column_stack((ecdf_x, ecdf_y))

            # be = export_mat_dict['bin_edges']
            # bh = export_mat_dict['bin_hgts']
            # if tgt_mrp is not None:
            # be_lte = be[be <= tgt_mrp]
            # x = np.diff(be_lte)
            # export_mat_dict['prob_le_tgt_area'] = np.sum(x * bh[:len(x)]) + (tgt_mrp - be_lte[-1]) * bh[len(x)]
            # be_lte = be[be <= mrp_det]
            # x = np.diff(be_lte)
            # export_mat_dict['prob_le_det_area'] = np.sum(x * bh[:len(x)]) + (mrp_det - be_lte[-1]) * bh[len(x)]

        else:
            return

        axs = [ax1, ax2]
        for ax in axs:
            [item.set_fontsize(fs) for item in
             ([ax.title, ax.xaxis.label, ax.yaxis.label, ax.xaxis.get_offset_text(), ax.yaxis.get_offset_text()]
              + ax.get_xticklabels() + ax.get_yticklabels())]

        if save_mat:
            for_which_str = '_'.join(for_which)
            file_name = f'plot_mrp_histogram_ds_{self.tag}_{self.structure.name}_{for_which_str}.mat'
            sio.savemat(Utility.get_path(save_dir_path, file_name), export_mat_dict)

        return fig, axs, export_mat_dict

    # ------------------------------------------------------------------------------------------------------------------
