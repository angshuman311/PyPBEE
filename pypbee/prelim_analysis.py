# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 14:01:56 2019

@author: Angshuman Deb
"""

from .analysis import Analysis
from .utility import Utility
import os
import numpy as np
from time import time
from functools import partial
from .structure import Structure
from pathos.multiprocessing import ProcessPool as Pool


class PrelimAnalysis(Analysis):

    ####################################################################################################################
    # Constructor
    ####################################################################################################################

    def __init__(self, structure, num_modes, **kwargs):
        comp_env = kwargs.get('comp_env', 'local')
        self.num_modes = num_modes
        super().__init__('Prelim_Analysis', structure, comp_env)

    ####################################################################################################################
    # Static methods
    ####################################################################################################################

    @staticmethod
    def run_parallel(results_dir_path, setup_dir_path, run_prelim_analysis, comp_env,
                     node_index, batch_index, analysis_index, run_case):
        start_time = time()

        for_which = list(run_case.astype(str))
        analysis_case = for_which[Structure.get_dir_level_index('Analysis_Case')]
        run_case_str = ' '.join(for_which)

        target_work_dir_path = Utility.get_path(results_dir_path, *for_which)

        run_prelim_analysis(target_work_dir_path, comp_env)

        status_file_name = f'status_analysis_index_{int(analysis_index)}' \
                           f'_node_{int(node_index)}_batch_{int(batch_index)}.txt'
        status_file_path = Utility.get_path(setup_dir_path, analysis_case, status_file_name)

        with open(status_file_path, 'w') as fid:
            fid.write(f"{run_case_str} : COMPLETED in {time() - start_time} seconds\n")

        return

    ####################################################################################################################
    # Public functionalities (instance methods)
    ####################################################################################################################

    def stage(self, analysis_case, design_num_list, rng_seed):

        structure = self.structure

        results_dir_path = self.get_results_dir_path()

        set_rng_seed, seed_multiplier = Utility.setup_rng(rng_seed)

        iter_model_option_max = [0] * len(design_num_list)
        iter_model_realztn_max = [0] * len(design_num_list)

        model_attributes_values = structure.create_model_options(analysis_case)

        for design_num in design_num_list:

            if set_rng_seed:
                rng_seed_send = seed_multiplier * Utility.get_num([int(design_num)])
            else:
                rng_seed_send = rng_seed

            random_model_param_vals = structure.generate_random_model_param_vals(
                analysis_case, design_num, rng_seed=rng_seed_send)

            # random_model_param_vals = structure.get_random_model_param_vals([analysis_case, design_num])

            # number of realizations for prob dist params of FE model params
            iter_model_option_max[design_num_list.index(design_num)] = random_model_param_vals.shape[2]
            range_iter_model_option = range(
                1, iter_model_option_max[design_num_list.index(design_num)] + 1)

            # number of realizations for FE model params
            iter_model_realztn_max[design_num_list.index(design_num)] = random_model_param_vals.shape[0]
            range_iter_model_realztn = range(1, iter_model_realztn_max[design_num_list.index(design_num)] + 1)

            for iter_model_option in range_iter_model_option:
                for iter_model_realztn in range_iter_model_realztn:

                    work_dir_path = Utility.get_path(results_dir_path, analysis_case,
                                                     design_num, str(iter_model_option),
                                                     str(iter_model_realztn))

                    if not os.path.isdir(work_dir_path):
                        os.makedirs(work_dir_path)

                    for_which = [analysis_case, design_num, str(iter_model_option), str(iter_model_realztn)]
                    structure.write_model_param_vals('Prelim_Analysis_Results', for_which, random_model_param_vals=random_model_param_vals)
                    structure.write_model_attributes('Prelim_Analysis_Results', for_which, model_attributes_values=model_attributes_values)
                    structure.structural_analysis_platform.write_run_prelim_analysis(work_dir_path, self.num_modes)
        return

    # ------------------------------------------------------------------------------------------------------------------

    def setup(self, analysis_case, design_num_list, **kwargs):
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
                    list(range(1, structure.get_random_model_param_vals([analysis_case, design_num]).shape[0] + 1))
                )
            )
        analysis_list = np.vstack(temp)
        list_file_path = Utility.get_path(setup_dir_path, analysis_case, 'analysis_list.txt')
        Utility.save_array_as_text_file(list_file_path, analysis_list, fmt='%d',
                                        header=f'Total analyses to run = {analysis_list.shape[0]:d}')

        self.set_analysis_setup_info(analysis_case, {'design_num_list': design_num_list})

        if not self.comp_env == 'local':
            self.setup_non_local_run(analysis_case, n_batch, n_job, pool_size_factor, run_time_limit, allocation_name)

        return

    # ------------------------------------------------------------------------------------------------------------------

    def run(self, analysis_case, pool_size, **kwargs):
        n_batch = kwargs.get('n_batch', 1)
        n_job = kwargs.get('n_job', 1)
        node_index = kwargs.get('node_index', 1)
        batch_index = kwargs.get('batch_index', 1)

        structure = self.structure
        run_prelim_analysis = structure.structural_analysis_platform.run_prelim_analysis
        results_dir_path = self.get_results_dir_path()
        setup_dir_path = self.get_setup_dir_path()
        analysis_list, analysis_range, n_sim = self.get_analysis_info(analysis_case, 'analysis_list.txt',
                                                                      n_batch, n_job, node_index, batch_index)

        to_run = partial(PrelimAnalysis.run_parallel, results_dir_path, setup_dir_path, run_prelim_analysis,
                         self.comp_env, node_index, batch_index)

        to_pass = [[analysis_index for analysis_index in analysis_range],
                   [analysis_list[analysis_index - 1, :] for analysis_index in analysis_range], ]

        if pool_size > 1:
            my_pool = Pool(pool_size)
            my_pool.map(to_run, *to_pass)
            my_pool.close()
            my_pool.join()
            my_pool.clear()
        else:
            list(map(to_run, *to_pass))

        for analysis_index in analysis_range:
            run_case = analysis_list[analysis_index - 1, :]
            for_which = list(run_case.astype(str))
            structure.write_damping_file(for_which)

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

        # Write final_status.txt
        Utility.combine_contents_of_files(Utility.get_path(setup_dir_path, analysis_case),
                                          ['startswith final_status_', 'startswith status_analysis_index_'],
                                          'final_status.txt',
                                          file_ext='.txt', delete_original=True)

        return

    # ------------------------------------------------------------------------------------------------------------------
