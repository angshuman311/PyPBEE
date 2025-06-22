# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 20:36:32 2019

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
from benedict import benedict


class PSDemHA(Analysis):

    ####################################################################################################################
    # Constructor
    ####################################################################################################################

    def __init__(self, edp_list, im, **kwargs):
        comp_env = kwargs.get('comp_env', 'local')
        self.edp_list = edp_list
        self.im = im
        super().__init__(self.__class__.__name__, im.structure, comp_env)

    ####################################################################################################################
    # Static methods
    ####################################################################################################################

    @staticmethod
    def run_parallel(setup_dir_path, compute_demand_hazard_integral_list,
                     haz_lev_list, n_gm_list, im, delta_input_list, min_max_scale_fac_list,
                     node_index, batch_index, analysis_index, run_case):

        start_time = time()

        for_which = list(run_case.astype(str))
        analysis_case = for_which[Structure.get_dir_level_index('Analysis_Case')]
        run_case_str = ' '.join(for_which)

        psdemha_results = list()
        for itr in range(len(compute_demand_hazard_integral_list)):
            temp = compute_demand_hazard_integral_list[itr](for_which, haz_lev_list, im, n_gm_list,
                                                            min_max_scale_fac=min_max_scale_fac_list[itr],
                                                            delta_input=delta_input_list[itr])
            psdemha_results.append(temp)

        status_file_name = f'status_analysis_index_{int(analysis_index)}' \
                           f'_node_{int(node_index)}_batch_{int(batch_index)}.txt'
        status_file_path = Utility.get_path(setup_dir_path, analysis_case, status_file_name)

        with open(status_file_path, 'w') as fid:
            fid.write(f"{run_case_str} : COMPLETED in {time() - start_time} seconds\n")

        return f'{run_case_str}', psdemha_results

    ####################################################################################################################
    # Public functionalities (instance methods)
    ####################################################################################################################

    def stage(self, analysis_case, design_num_list):

        results_dir_path = self.get_results_dir_path()

        for design_num in design_num_list:
            target_work_dir_path = Utility.get_path(results_dir_path, analysis_case, design_num)

            if not os.path.isdir(target_work_dir_path):
                os.makedirs(target_work_dir_path)

        return

    # ------------------------------------------------------------------------------------------------------------------

    def setup(self, analysis_case, design_num_list, haz_lev_list, n_gm_list=None, **kwargs):
        run_time_limit = kwargs.get('run_time_limit', '01:00:00')
        allocation_name = kwargs.get('allocation_name', '')
        n_batch = kwargs.get('n_batch', 1)
        n_job = kwargs.get('n_job', 1)
        pool_size_factor = kwargs.get('pool_size_factor', 1)

        structure = self.structure

        self.stage(analysis_case, design_num_list)

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
        edp_list = self.edp_list
        n_batch = kwargs.get('n_batch', 1)
        n_job = kwargs.get('n_job', 1)
        node_index = kwargs.get('node_index', 1)
        batch_index = kwargs.get('batch_index', 1)
        delta_input_list = kwargs.get('delta_input_list', [np.array([])] * len(edp_list))
        min_max_scale_fac_list = kwargs.get('min_max_scale_fac_list', [[1, 1]] * len(edp_list))

        compute_demand_hazard_integral_list = [edp.compute_demand_hazard_integral for edp in edp_list]
        im = self.im

        setup_dir_path = self.get_setup_dir_path()
        results_dir_path = self.get_results_dir_path()
        analysis_list, analysis_range, n_sim = self.get_analysis_info(analysis_case, 'analysis_list.txt',
                                                                      n_batch, n_job, node_index, batch_index)

        analysis_setup_info = self.get_analysis_setup_info(analysis_case)
        design_num_list = list(np.unique(analysis_list[:, 1]).astype(str))
        haz_lev_list = analysis_setup_info['haz_lev_list']
        n_gm_list = analysis_setup_info['n_gm_list']

        to_run = partial(PSDemHA.run_parallel, setup_dir_path, compute_demand_hazard_integral_list,
                         haz_lev_list, n_gm_list, im, delta_input_list, min_max_scale_fac_list,
                         node_index, batch_index)

        to_pass = [[analysis_index for analysis_index in analysis_range],
                   [analysis_list[analysis_index - 1, :] for analysis_index in analysis_range], ]

        if pool_size > 1:
            my_pool = Pool(pool_size)
            psdemha_results = dict(my_pool.map(to_run, *to_pass))
            my_pool.close()
            my_pool.join()
            my_pool.clear()
        else:
            psdemha_results = dict(list(map(to_run, *to_pass)))

        edp_tags = [edp.tag for edp in edp_list]

        for design_num in design_num_list:
            target_work_dir_path = Utility.get_path(results_dir_path, analysis_case, design_num)
            psdemha_results_design_num = {key: psdemha_results[key] for key in psdemha_results.keys() if
                                          key.split()[Structure.get_dir_level_index('Design_Num')] == design_num}
            psdemha_results_edp = {
                edp_tag: {analysis_key: psdemha_results_design_num[analysis_key][edp_tags.index(edp_tag)]
                          for analysis_key in psdemha_results_design_num.keys()}
                for edp_tag in edp_tags
            }

            if n_batch == 1 and n_job == 1:
                Utility.pickle_dump_dict(Utility.get_path(target_work_dir_path, f'psdemha_results.pickle'),
                                         psdemha_results_edp)
            else:
                Utility.pickle_dump_dict(
                    Utility.get_path(
                        target_work_dir_path, f'psdemha_results_node_{node_index}_batch_{batch_index}.pickle'
                    ), psdemha_results_edp
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
            files = [file for file in files if file.startswith('psdemha_results') and file.endswith(
                '.pickle') and 'node' in file and 'batch' in file]
            if len(files) > 0:
                list_of_dicts = [Utility.pickle_load_dict(Utility.get_path(target_work_dir_path, file)) for file in
                                 files]

                psdemha_results_file_path = Utility.get_path(target_work_dir_path, 'psdemha_results.pickle')
                if os.path.isfile(psdemha_results_file_path):
                    psdemha_results = benedict(Utility.pickle_load_dict(psdemha_results_file_path))
                else:
                    psdemha_results = benedict()

                for d in list_of_dicts:
                    psdemha_results.merge(d)

                Utility.pickle_dump_dict(psdemha_results_file_path, psdemha_results)

                for file in files:
                    os.remove(Utility.get_path(target_work_dir_path, file))

        return

    # ------------------------------------------------------------------------------------------------------------------
