# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 16:31:34 2019

@author: Angshuman Deb
"""

from .analysis import Analysis
from .utility import Utility
from .structure import Structure
import os
import numpy as np
from time import time
from pathos.multiprocessing import ProcessPool as Pool
from functools import partial
import scipy.io as sio
from benedict import benedict
import matplotlib.pyplot as plt
import matplotlib.colors as colors


class PSDamHA(Analysis):

    ####################################################################################################################
    # Constructor
    ####################################################################################################################

    def __init__(self, ds_list, im=None, **kwargs):
        comp_env = kwargs.get('comp_env', 'local')
        self.sol_type = kwargs.get('sol_type', 'numerical')
        self.ds_list = ds_list
        self.im = im
        super().__init__(self.__class__.__name__, ds_list[0].structure, comp_env)

    ####################################################################################################################
    # Static methods
    ####################################################################################################################

    @staticmethod
    def run_parallel(setup_dir_path, compute_damage_hazard_integral_list, sol_type,
                     im, haz_lev_list, n_gm_list, delta_input_list, min_max_scale_fac_list,
                     node_index, batch_index, analysis_index, run_case):

        start_time = time()

        for_which = list(run_case.astype(str))
        analysis_case = for_which[Structure.get_dir_level_index('Analysis_Case')]
        run_case_str = ' '.join(for_which)

        psdamha_results = list()
        for itr in range(len(compute_damage_hazard_integral_list)):
            temp = compute_damage_hazard_integral_list[itr](for_which, im, haz_lev_list, n_gm_list,
                                                            min_max_scale_fac=min_max_scale_fac_list[itr],
                                                            delta_input=delta_input_list[itr],
                                                            sol_type=sol_type)
            psdamha_results.append(temp)

        status_file_name = f'status_analysis_index_{int(analysis_index)}' \
                           f'_node_{int(node_index)}_batch_{int(batch_index)}.txt'
        status_file_path = Utility.get_path(setup_dir_path, analysis_case, status_file_name)

        with open(status_file_path, 'w') as fid:
            fid.write(f"{run_case_str} : COMPLETED in {time() - start_time} seconds\n")

        return f'{run_case_str}', psdamha_results

    ####################################################################################################################
    # Public functionalities (instance methods)
    ####################################################################################################################

    def stage(self, analysis_case, design_num_list, rng_seed):

        ds_list = self.ds_list

        results_dir_path = self.get_results_dir_path()

        set_rng_seed, seed_multiplier = Utility.setup_rng(rng_seed)

        for design_num in design_num_list:

            target_work_dir_path = Utility.get_path(results_dir_path, analysis_case, design_num)

            if not os.path.isdir(target_work_dir_path):
                os.makedirs(target_work_dir_path)

            normalized_fragility_dist_param_vals = dict()

            for ds in ds_list:
                if set_rng_seed:
                    rng_seed_send = seed_multiplier * Utility.get_num([int(design_num), ds.tag])
                else:
                    rng_seed_send = rng_seed

                normalized_fragility_dist_param_vals[ds.tag] = \
                    ds.generate_normalized_fragility_dist_param_vals(analysis_case, rng_seed=rng_seed_send)

            normalized_fragility_dist_param_vals_file_path = Utility.get_path(
                target_work_dir_path, 'normalized_fragility_dist_param_vals.pickle')
            Utility.pickle_dump_dict(normalized_fragility_dist_param_vals_file_path,
                                     normalized_fragility_dist_param_vals)

        return

    # ------------------------------------------------------------------------------------------------------------------

    def setup(self, analysis_case, design_num_list, haz_lev_list=None, n_gm_list=None, **kwargs):
        rng_seed = kwargs.get('rng_seed', None)
        run_time_limit = kwargs.get('run_time_limit', '01:00:00')
        allocation_name = kwargs.get('allocation_name', '')
        n_batch = kwargs.get('n_batch', 1)
        n_job = kwargs.get('n_job', 1)
        pool_size_factor = kwargs.get('pool_size_factor', 1)

        structure = self.structure

        self.stage(analysis_case, design_num_list, rng_seed)

        setup_dir_path = self.get_setup_dir_path()
        if not os.path.isdir(Utility.get_path(setup_dir_path, analysis_case)):
            os.makedirs(Utility.get_path(setup_dir_path, analysis_case))

        temp = list()
        for design_num in design_num_list:
            temp.append(
                Utility.get_multi_level_iterator_cases(
                    int(analysis_case),
                    int(design_num),
                    list(range(1, structure.get_random_model_param_vals([analysis_case, design_num]).shape[2] + 1)),
                    list(range(1, structure.get_random_model_param_vals([analysis_case, design_num]).shape[0] + 1)),
                )
            )
        analysis_list = np.vstack(temp)

        list_file_path = Utility.get_path(setup_dir_path, analysis_case, 'analysis_list.txt')
        Utility.save_array_as_text_file(list_file_path, analysis_list, fmt='%d',
                                        header=f'Total analyses to run = {analysis_list.shape[0]:d}')

        self.set_analysis_setup_info(
            analysis_case,
            {'design_num_list': design_num_list,
             'haz_lev_list': haz_lev_list,
             'n_gm_list': n_gm_list,
             }
        )

        if not self.comp_env == 'local':
            self.setup_non_local_run(analysis_case, n_batch, n_job, pool_size_factor, run_time_limit, allocation_name)

        return

    # ------------------------------------------------------------------------------------------------------------------

    def run(self, analysis_case, pool_size, **kwargs):
        ds_list = self.ds_list
        n_batch = kwargs.get('n_batch', 1)
        n_job = kwargs.get('n_job', 1)
        node_index = kwargs.get('node_index', 1)
        batch_index = kwargs.get('batch_index', 1)
        delta_input_list = kwargs.get('delta_input_list', [np.array([])] * len(ds_list))
        min_max_scale_fac_list = kwargs.get('min_max_scale_fac_list', [[1, 1]] * len(ds_list))

        compute_damage_hazard_integral_list = [ds.compute_damage_hazard_integral for ds in ds_list]
        im = self.im

        setup_dir_path = self.get_setup_dir_path()
        results_dir_path = self.get_results_dir_path()
        analysis_list, analysis_range, n_sim = self.get_analysis_info(analysis_case, 'analysis_list.txt',
                                                                      n_batch, n_job, node_index, batch_index)

        analysis_setup_info = self.get_analysis_setup_info(analysis_case)
        design_num_list = list(np.unique(analysis_list[:, 1]).astype(str))
        haz_lev_list = analysis_setup_info['haz_lev_list']
        n_gm_list = analysis_setup_info['n_gm_list']

        to_run = partial(PSDamHA.run_parallel, setup_dir_path, compute_damage_hazard_integral_list, self.sol_type,
                         im, haz_lev_list, n_gm_list, delta_input_list, min_max_scale_fac_list,
                         node_index, batch_index)

        to_pass = [[analysis_index for analysis_index in analysis_range],
                   [analysis_list[analysis_index - 1, :] for analysis_index in analysis_range], ]

        if pool_size > 1:
            my_pool = Pool(pool_size)
            psdamha_results = dict(my_pool.map(to_run, *to_pass))
            my_pool.close()
            my_pool.join()
            my_pool.clear()
        else:
            psdamha_results = dict(list(map(to_run, *to_pass)))

        ds_tags = [ds.tag for ds in ds_list]

        for design_num in design_num_list:
            target_work_dir_path = Utility.get_path(results_dir_path, analysis_case, design_num)
            psdamha_results_design_num = {key: psdamha_results[key] for key in psdamha_results.keys() if
                                          key.split()[Structure.get_dir_level_index('Design_Num')] == design_num}
            psdamha_results_ds = {
                ds_tag: {analysis_key: psdamha_results_design_num[analysis_key][ds_tags.index(ds_tag)]
                         for analysis_key in psdamha_results_design_num.keys()}
                for ds_tag in ds_tags
            }

            if n_batch == 1 and n_job == 1:
                Utility.pickle_dump_dict(Utility.get_path(target_work_dir_path, f'psdamha_results.pickle'),
                                         psdamha_results_ds)
            else:
                Utility.pickle_dump_dict(
                    Utility.get_path(
                        target_work_dir_path, f'psdamha_results_node_{node_index}_batch_{batch_index}.pickle'
                    ), psdamha_results_ds
                )

        # Write final_status_node_{}_batch_{}.txt
        contains_str = f'_node_{int(node_index)}_batch_{int(batch_index)}'
        Utility.combine_contents_of_files(Utility.get_path(setup_dir_path, analysis_case),
                                          [f'contains status_analysis_index_ {contains_str}'],
                                          f'final_status{contains_str}.txt',
                                          file_ext='.txt', delete_original=True)

        return

    # ------------------------------------------------------------------------------------------------------------------

    def wrap_up(self, analysis_case):
        setup_dir_path = self.get_setup_dir_path()
        results_dir_path = self.get_results_dir_path()

        # Write final_status.txt
        Utility.combine_contents_of_files(Utility.get_path(setup_dir_path, analysis_case),
                                          ['startswith final_status_', 'startswith status_analysis_index_'],
                                          'final_status.txt',
                                          file_ext='.txt', delete_original=True)

        analysis_list, _, _ = self.get_analysis_info(analysis_case, 'analysis_list.txt', 1, 1, 1, 1)

        # Combine results from all nodes and batches
        design_num_list = list(np.unique(analysis_list[:, 1]).astype(str))

        for design_num in design_num_list:
            target_work_dir_path = Utility.get_path(results_dir_path, analysis_case, design_num)
            files = os.listdir(target_work_dir_path)
            files = [file for file in files if file.startswith('psdamha_results') and file.endswith(
                '.pickle') and 'node' in file and 'batch' in file]
            if len(files) > 0:
                list_of_dicts = [Utility.pickle_load_dict(Utility.get_path(target_work_dir_path, file)) for file in
                                 files]

                psdamha_results_file_path = Utility.get_path(target_work_dir_path, 'psdamha_results.pickle')
                if os.path.isfile(psdamha_results_file_path):
                    psdamha_results = benedict(Utility.pickle_load_dict(psdamha_results_file_path))
                else:
                    psdamha_results = benedict()

                for d in list_of_dicts:
                    psdamha_results.merge(d)

                Utility.pickle_dump_dict(psdamha_results_file_path, psdamha_results)

                for file in files:
                    os.remove(Utility.get_path(target_work_dir_path, file))

        return

    # ------------------------------------------------------------------------------------------------------------------

    def get_d_star(self, analysis_case, target_mrp_list, d_list, **kwargs):
        exclude_ds_tags = kwargs.get('exclude_ds_tags', list())
        ds_list = self.ds_list
        ds_list = [ds for ds in ds_list if ds.tag not in exclude_ds_tags]

        x_star_list = []
        d_star_list = []
        for (ds, tgt) in zip(ds_list, target_mrp_list):
            d_star_temp, x_star_temp = ds.get_d_star(analysis_case, d_list, tgt)
            x_star_list.append(x_star_temp)
            d_star_list.append(d_star_temp)
        x_star_list = np.array(x_star_list)
        d_star_list = np.array(d_star_list)

        x_star = np.max(x_star_list)
        ind = np.argmax(x_star_list)
        d_star = d_star_list[ind, :]

        return d_star, x_star

    # ------------------------------------------------------------------------------------------------------------------

    def plot_feasible_domain(self, analysis_case, target_mrp_list, color_list, **kwargs):
        mesh_size = kwargs.get('mesh_size', 100)
        save_dir_path = kwargs.get('save_dir_path', Utility.get_path(os.getcwd()))
        save_mat = kwargs.get('save_mat', False)
        figkwargs = kwargs.get('figkwargs', dict())
        grid_alpha = kwargs.get('grid_alpha', 0.5)
        minor_grid_alpha = kwargs.get('minor_grid_alpha', 0.0)
        patch_alpha = kwargs.get('patch_alpha', 0.25)
        lw = kwargs.get('lw', 2)
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
        exclude_ds_tags = kwargs.get('exclude_ds_tags', list())
        fig_ax = kwargs.get('fig_ax', None)
        domain_type = kwargs.get('domain_type', 'full_fledged')
        d_list = kwargs.get('d_list', [])
        da_list = kwargs.get('da_list', [])

        if domain_type not in ['full_fledged', 'simplified_bilinear']:
            raise ValueError("Invalid domain_type!")

        surf_type = 'piecewise_linear'
        if domain_type == 'simplified_bilinear':
            if len(d_list) < 2:
                raise ValueError("Invalid d_list!")
            surf_type = 'bilinear_contours'
            if len(da_list) != 2:
                raise ValueError("Invalid da_list!")

        if fig_ax is None or np.array(fig_ax, dtype='object').any() is None:
            fig = plt.figure(**figkwargs)
            ax = fig.gca()
        else:
            fig = fig_ax[0]
            ax = fig_ax[1]
        ax.minorticks_on()
        ax.grid(True, which="major", alpha=grid_alpha)
        ax.grid(True, which="minor", alpha=minor_grid_alpha)

        ds_list = self.ds_list
        ds_list = [ds for ds in ds_list if ds.tag not in exclude_ds_tags]
        structure = self.structure
        model_params = structure.model_params
        design_vector = np.array(list(model_params['primary_design_params']['value_list_dict'].values()))
        design_num_list = list(model_params['primary_design_params']['value_list_dict'].keys())

        export_mat_dict = dict()
        export_mat_dict['design_vector'] = np.array(design_vector)
        export_mat_dict['design_num_list'] = np.array(design_num_list)
        export_mat_dict['target_mrp_list'] = np.array(target_mrp_list).astype(float)
        export_mat_dict['mesh_size'] = mesh_size

        safe_zone = np.zeros((mesh_size, mesh_size))

        for itr in range(len(ds_list)):
            ds = ds_list[itr]
            color = color_list[itr]
            target_mrp = target_mrp_list[itr]
            data, surface = ds.get_damage_mrp_interp_surface(analysis_case, mesh_size, surf_type,
                                                             d_list=d_list, da_list=da_list)
            ax.contour(surface[:, :, 0], surface[:, :, 1], surface[:, :, -1], [target_mrp],
                       colors=color, linewidths=lw)
            safe_zone[surface[:, :, -1] < target_mrp] = 1
            tag = ds.tag
            export_mat_dict[f'mrp_ds_{tag}'] = data[:, -1]
            export_mat_dict[f'surface_ds_{tag}'] = surface

        cmap = colors.ListedColormap([[0, 1, 0, 1], [1, 1, 1, 1]])

        ax.imshow(safe_zone,
                  extent=[min(design_vector[:, 0]),
                          max(design_vector[:, 0]),
                          min(design_vector[:, 1]),
                          max(design_vector[:, 1])],
                  aspect='auto',
                  origin='lower',
                  alpha=patch_alpha,
                  cmap=cmap)

        fig, ax, _ = structure.plot_design_space(
            grid_alpha=grid_alpha, minor_grid_alpha=minor_grid_alpha,
            save_mat=False, lw=0,
            mc_gp=mc_gp, ms_gp=ms_gp, mt_gp=mt_gp,
            mc_ngp=mc_ngp, ms_ngp=ms_ngp, mt_ngp=mt_ngp,
            mc_ad=mc_ad, ms_ad=ms_ad, mt_ad=mt_ad,
            mc_mps=mc_mps, ms_mps=ms_mps, mt_mps=mt_mps, mark_points=mark_points,
            fig_ax=[fig, ax]
        )

        export_mat_dict['safe_zone'] = safe_zone

        if save_mat:
            name = structure.name
            file_name = f'plot_feasible_domain_{name}_{analysis_case}.mat'
            sio.savemat(Utility.get_path(save_dir_path, file_name), export_mat_dict)

        return fig, ax, export_mat_dict

    # ------------------------------------------------------------------------------------------------------------------
