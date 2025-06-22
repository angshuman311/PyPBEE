# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 14:01:56 2019

@author: Angshuman Deb
"""

from .analysis import Analysis
from .utility import Utility
from .structure import Structure
import os
import numpy as np
from time import time
from functools import partial
from pathos.multiprocessing import ProcessPool as Pool
from benedict import benedict


class PSHA(Analysis):

    ####################################################################################################################
    # Constructor
    ####################################################################################################################

    def __init__(self, im, **kwargs):
        comp_env = kwargs.get('comp_env', 'local')
        self.im = im
        super().__init__(self.__class__.__name__, im.structure, comp_env)

    ####################################################################################################################
    # Static methods
    ####################################################################################################################

    @staticmethod
    def run_parallel(setup_dir_path, compute_seismic_hazard_integral, im_input, node_index, batch_index,
                     analysis_index, run_case):
        start_time = time()

        for_which = list(run_case.astype(str))
        analysis_case = for_which[Structure.get_dir_level_index('Analysis_Case')]
        run_case_str = ' '.join(for_which)

        psha_results = compute_seismic_hazard_integral(for_which, im_input=im_input)

        status_file_name = f'status_analysis_index_{int(analysis_index)}' \
                           f'_node_{int(node_index)}_batch_{int(batch_index)}.txt'
        status_file_path = Utility.get_path(setup_dir_path, analysis_case, status_file_name)

        with open(status_file_path, 'w') as fid:
            fid.write(f"{run_case_str} : COMPLETED in {time() - start_time} seconds\n")

        return f'{run_case_str}', psha_results

    ####################################################################################################################
    # Public functionalities (instance methods)
    ####################################################################################################################

    def stage(self, analysis_case, design_num_list):

        self.structure.set_site_hazard_info(self.im)

        results_dir_path = self.get_results_dir_path()

        for design_num in design_num_list:
            target_work_dir_path = Utility.get_path(results_dir_path, analysis_case, design_num)

            if not os.path.isdir(target_work_dir_path):
                os.makedirs(target_work_dir_path)

        return

    # ------------------------------------------------------------------------------------------------------------------

    def setup(self, analysis_case, design_num_list, **kwargs):
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
        im_input = kwargs.get('im_input', np.array([]))

        compute_seismic_hazard_integral = self.im.compute_seismic_hazard_integral

        setup_dir_path = self.get_setup_dir_path()
        results_dir_path = self.get_results_dir_path()
        analysis_list, analysis_range, n_sim = self.get_analysis_info(analysis_case, 'analysis_list.txt',
                                                                      n_batch, n_job, node_index, batch_index)

        design_num_list = list(np.unique(analysis_list[:, 1]).astype(str))

        to_run = partial(PSHA.run_parallel, setup_dir_path, compute_seismic_hazard_integral, im_input,
                         node_index, batch_index)

        to_pass = [[analysis_index for analysis_index in analysis_range],
                   [analysis_list[analysis_index - 1, :] for analysis_index in analysis_range], ]

        if pool_size > 1:
            my_pool = Pool(pool_size)
            psha_results = dict(my_pool.map(to_run, *to_pass))
            my_pool.close()
            my_pool.join()
            my_pool.clear()
        else:
            psha_results = dict(list(map(to_run, *to_pass)))

        for design_num in design_num_list:
            target_work_dir_path = Utility.get_path(results_dir_path, analysis_case, design_num)
            psha_results_design_num = {key: psha_results[key] for key in psha_results.keys() if
                                       key.split()[Structure.get_dir_level_index('Design_Num')] == design_num}
            if n_batch == 1 and n_job == 1:
                Utility.pickle_dump_dict(Utility.get_path(target_work_dir_path, f'psha_results.pickle'),
                                         psha_results_design_num)
            else:
                Utility.pickle_dump_dict(
                    Utility.get_path(
                        target_work_dir_path, f'psha_results_node_{node_index}_batch_{batch_index}.pickle'
                    ), psha_results_design_num
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
            files = [file for file in files if file.startswith('psha_results') and file.endswith(
                '.pickle') and 'node' in file and 'batch' in file]
            if len(files) > 0:
                list_of_dicts = [Utility.pickle_load_dict(Utility.get_path(target_work_dir_path, file)) for file in
                                 files]

                psha_results_file_path = Utility.get_path(target_work_dir_path, 'psha_results.pickle')
                if os.path.isfile(psha_results_file_path):
                    psha_results = benedict(Utility.pickle_load_dict(psha_results_file_path))
                else:
                    psha_results = benedict()

                for d in list_of_dicts:
                    psha_results.merge(d)

                Utility.pickle_dump_dict(psha_results_file_path, psha_results)

                for file in files:
                    os.remove(Utility.get_path(target_work_dir_path, file))

        return

    # ------------------------------------------------------------------------------------------------------------------
