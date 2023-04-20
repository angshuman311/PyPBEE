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
import shutil
from copy import deepcopy
# import platform


class NLTHA(Analysis):

    ####################################################################################################################
    # Constructor
    ####################################################################################################################

    def __init__(self, edp_list, im, **kwargs):
        comp_env = kwargs.get('comp_env', 'local')
        self.edp_list = edp_list
        self.im = im
        super().__init__(self.__class__.__name__, edp_list[0].structure, comp_env)

    ####################################################################################################################
    # Static methods
    ####################################################################################################################

    @staticmethod
    def stage_parallel(results_dir_path, model_fetch_dir_path, write_run_nltha,
                       generate_recorder_list, rec_gen_file_open_mode_list,
                       gen_rec_list, get_gm_rec_info, gm_database_dir_path, ai_end,
                       run_case, model_fetch_case, gm_fetch_case):

        for_which = list(run_case.astype(str))
        model_fetch_for_which = list(model_fetch_case.astype(str))
        gm_fetch_for_which = list(gm_fetch_case.astype(str))

        target_work_dir_path = Utility.get_path(results_dir_path, *for_which)
        work_dir_path_prelim_analysis = Utility.get_path(model_fetch_dir_path, *model_fetch_for_which)

        # make target work dir if it doesn't exist
        if not os.path.isdir(target_work_dir_path):
            os.makedirs(target_work_dir_path)

        status_file_path = Utility.get_path(target_work_dir_path, 'NLTHA_STATUS.txt')
        if os.path.isfile(status_file_path):
            os.remove(status_file_path)

        # generate edp recorders
        for i_edp in range(len(generate_recorder_list)):
            generate_recorder_list[i_edp](for_which, rec_gen_file_open_mode=rec_gen_file_open_mode_list[i_edp],
                                          gen_rec=gen_rec_list[i_edp])

        # write run_nltha
        write_run_nltha(target_work_dir_path, work_dir_path_prelim_analysis)

        # fetch ground motions
        if gm_database_dir_path is not None:
            gm_file_names, scale_fac = get_gm_rec_info([*gm_fetch_for_which,
                                                        for_which[
                                                            Structure.get_dir_level_index('Ground_Motion_Rec_Num')]])
            Utility.fetch_selected_ground_motion_records(
                gm_database_dir_path, gm_file_names, scale_fac, target_work_dir_path, ai_end)

        run_case_str = ' '.join(for_which)
        print(f'Staging analysis case {run_case_str} complete!')

        return

    # ------------------------------------------------------------------------------------------------------------------

    @staticmethod
    def run_parallel(results_dir_path, setup_dir_path, run_nltha, comp_env, node_index, batch_index,
                     analysis_index, run_case):
        start_time = time()

        for_which = list(run_case.astype(str))
        analysis_case = for_which[Structure.get_dir_level_index('Analysis_Case')]
        run_case_str = ' '.join(for_which)

        target_work_dir_path = Utility.get_path(results_dir_path, *for_which)

        run_nltha(target_work_dir_path, comp_env)

        status = 'FAIL'
        nltha_status_file_path = Utility.get_path(target_work_dir_path, 'NLTHA_STATUS.txt')
        if os.path.isfile(nltha_status_file_path):
            with open(nltha_status_file_path, 'r') as fid:
                status = fid.readline().strip()

        status_file_name = f'status_analysis_index_{int(analysis_index)}' \
                           f'_node_{int(node_index)}_batch_{int(batch_index)}.txt'
        status_file_path = Utility.get_path(setup_dir_path, analysis_case, status_file_name)

        with open(status_file_path, 'w') as fid:
            if status == 'SUCCESS':
                fid.write(f"{run_case_str} : CONVERGED and COMPLETED in {time() - start_time} seconds\n")
            else:
                fid.write(f"{run_case_str} : FAILED in {time() - start_time} seconds\n")

        return

    # ------------------------------------------------------------------------------------------------------------------

    @staticmethod
    def run_parallel_prelim_analysis(results_dir_path, comp_env, run_prelim_analysis, run_case):
        for_which = list(run_case.astype(str))
        target_work_dir_path = Utility.get_path(results_dir_path, *for_which)
        run_prelim_analysis(target_work_dir_path, comp_env)
        return

    ####################################################################################################################
    # Public functionalities (instance methods)
    ####################################################################################################################

    def stage(self, pool_size, analysis_list, model_fetch_list, gm_fetch_list, gm_database_dir_path, **kwargs):

        ai_end = kwargs.get('ai_end', 0.7)

        edp_list = self.edp_list
        results_dir_path = self.get_results_dir_path()
        model_fetch_dir_path = results_dir_path.replace(self.analysis_type, 'Prelim_Analysis')

        for edp in edp_list:
            edp.write_evaluate_edp_and_helpers()

        generate_recorder_list = list()
        rec_gen_file_open_mode_list = ['w'] + ['a+'] * (len(edp_list) - 1)

        gen_rec_list = np.array([True] * len(edp_list))
        ind_shared = np.array([edp.recorder_file_storage == 'shared' for edp in edp_list])
        ind_shared[0] = False
        gen_rec_list[ind_shared] = False
        gen_rec_list = list(gen_rec_list)

        for edp in edp_list:
            generate_recorder_list.append(edp.generate_recorder)

        get_gm_rec_info = self.im.get_gm_rec_info

        n_sim = analysis_list.shape[0]

        to_run = partial(NLTHA.stage_parallel, results_dir_path, model_fetch_dir_path,
                         self.structure.structural_analysis_platform.write_run_nltha,
                         generate_recorder_list, rec_gen_file_open_mode_list,
                         gen_rec_list, get_gm_rec_info, gm_database_dir_path, ai_end)

        to_pass = [[analysis_list[itr, :] for itr in range(n_sim)],
                   [model_fetch_list[itr, :] for itr in range(n_sim)],
                   [gm_fetch_list[itr, :] for itr in range(n_sim)], ]

        if pool_size > 1:
            my_pool = Pool(pool_size)
            my_pool.map(to_run, *to_pass)
            my_pool.close()
            my_pool.join()
            my_pool.clear()
        else:
            list(map(to_run, *to_pass))

        return

    # ------------------------------------------------------------------------------------------------------------------

    def setup(self, analysis_case, pool_size, design_num_list, haz_lev_list, gm_database_dir_path,
              **kwargs):
        run_time_limit = kwargs.get('run_time_limit', '01:00:00')
        allocation_name = kwargs.get('allocation_name', '')
        n_batch = kwargs.get('n_batch', 1)
        n_job = kwargs.get('n_job', 1)
        pool_size_factor = kwargs.get('pool_size_factor', 1)
        stage = kwargs.get('stage', True)
        ai_end = kwargs.get('ai_end', 0.7)

        structure = self.structure
        im = self.im

        setup_dir_path = self.get_setup_dir_path()
        if not os.path.isdir(Utility.get_path(setup_dir_path, analysis_case)):
            os.makedirs(Utility.get_path(setup_dir_path, analysis_case))

        temp = list()
        for design_num in design_num_list:
            for haz_lev in haz_lev_list:
                _, _, _, n_gm = im.get_gms_results([analysis_case, design_num, '1', '1', haz_lev])
                temp.append(
                    Utility.get_multi_level_iterator_cases(
                        int(analysis_case),
                        int(design_num),
                        list(range(1, structure.get_random_model_param_vals([analysis_case, design_num]).shape[2] + 1)),
                        list(range(1, structure.get_random_model_param_vals([analysis_case, design_num]).shape[0] + 1)),
                        int(haz_lev),
                        list(range(1, n_gm + 1))
                    )
                )
        analysis_list = np.vstack(temp)
        model_fetch_list = analysis_list[:, :Structure.get_dir_level_index('Model_Realization_Num') + 1]
        gm_fetch_list = analysis_list[:, :Structure.get_dir_level_index('Hazard_Level_Num') + 1]

        list_file_path = Utility.get_path(setup_dir_path, analysis_case, 'analysis_list.txt')
        list_file_path_orig = Utility.get_path(setup_dir_path, analysis_case, 'analysis_list_original.txt')
        Utility.save_array_as_text_file(list_file_path, analysis_list, fmt='%d',
                                        header=f'Total analyses to run = {analysis_list.shape[0]:d}')
        Utility.save_array_as_text_file(list_file_path_orig, analysis_list, fmt='%d',
                                        header=f'Total analyses to run = {analysis_list.shape[0]:d}')

        if stage:
            self.stage(pool_size, analysis_list.astype(int), model_fetch_list.astype(int),
                       gm_fetch_list.astype(int), gm_database_dir_path, ai_end=ai_end)

        self.set_analysis_setup_info(
            analysis_case,
            {'design_num_list': design_num_list,
             'haz_lev_list': haz_lev_list,
             'allocation_name': allocation_name
             }
        )

        if not self.comp_env == 'local':
            self.setup_non_local_run(analysis_case, n_batch, n_job, pool_size_factor, run_time_limit, allocation_name)

        return

    # ------------------------------------------------------------------------------------------------------------------

    def run(self, analysis_case, pool_size, **kwargs):
        n_batch = kwargs.get('n_batch', 1)
        n_job = kwargs.get('n_job', 1)
        node_index = kwargs.get('node_index', 1)
        batch_index = kwargs.get('batch_index', 1)
        design_nbrs = kwargs.get('design_nbrs', list())
        model_opts = kwargs.get('model_opts', list())
        model_realizations = kwargs.get('model_realizations', list())
        haz_levs = kwargs.get('haz_levs', list())
        gm_nbrs = kwargs.get('gm_nbrs', list())
        restage = kwargs.get('restage', False)

        include_list = [design_nbrs, model_opts, model_realizations, haz_levs, gm_nbrs]

        setup_dir_path = self.get_setup_dir_path()
        results_dir_path = self.get_results_dir_path()
        analysis_list, analysis_range, n_sim = self.get_analysis_info(analysis_case, 'analysis_list.txt', n_batch,
                                                                      n_job, node_index, batch_index)

        structure = self.structure
        run_nltha = structure.structural_analysis_platform.run_nltha
        analysis_ind = np.ones((analysis_list.shape[0],), dtype=bool)

        dir_levels = Structure.get_all_dir_levels()[Structure.get_dir_level_index('Analysis_Case') + 1:]

        for itr in range(len(dir_levels)):
            if len(include_list[itr]) == 0:
                include_list[itr] = np.unique(analysis_list[:, Structure.get_dir_level_index(dir_levels[itr])])
            else:
                include_list[itr] = np.array(include_list[itr]).astype(int)

        for itr in range(len(dir_levels)):
            if not np.array_equal(np.sort(include_list[itr]),
                                  np.unique(analysis_list[:, Structure.get_dir_level_index(dir_levels[itr])])):
                temp = ~np.ones((analysis_list.shape[0],), dtype=bool)
                for incl_item in include_list[itr]:
                    temp |= analysis_list[:, itr + 1] == incl_item
                analysis_ind &= temp

        if not np.all(analysis_ind):
            analysis_list = analysis_list[analysis_ind, :]
            analysis_range = np.arange(1, len(analysis_list) + 1)
            restage = True

        if restage:
            self.restage(analysis_case, pool_size, None, analysis_list=analysis_list)

        to_run = partial(NLTHA.run_parallel, results_dir_path, setup_dir_path, run_nltha,
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

        analysis_ran, analysis_failed, analysis_skipped = self.get_post_analysis_lists(analysis_case)

        shutil.copy(Utility.get_path(setup_dir_path, analysis_case, 'analysis_list.txt'),
                    Utility.get_path(setup_dir_path, analysis_case, 'analysis_list_last.txt'))

        analysis_list_rerun = np.vstack([analysis_skipped, analysis_failed])
        if analysis_list_rerun.size != 0:
            analysis_list_rerun = np.unique(analysis_list_rerun, axis=0)
        self.setup_rerun(analysis_case, analysis_list_rerun)
        return

    # ------------------------------------------------------------------------------------------------------------------

    def setup_rerun(self, analysis_case, analysis_list_rerun, **kwargs):
        stage = kwargs.get('stage', False)
        pool_size = kwargs.get('pool_size', 1)
        gm_database_dir_path = kwargs.get('gm_database_dir_path', None)

        setup_dir_path = self.get_setup_dir_path()
        new_list_file_path = Utility.get_path(setup_dir_path, analysis_case, 'analysis_list.txt')
        Utility.save_array_as_text_file(new_list_file_path, analysis_list_rerun, fmt='%d',
                                        header=f'Total analyses to run = {analysis_list_rerun.shape[0]:d}')

        num_rerun = analysis_list_rerun.shape[0]
        if num_rerun > 0:
            rerun_stampede = input(f'{num_rerun} NLTHA did not converge! Re-run on Stampede2 [y/n] ? : ')
            if rerun_stampede == 'y' or rerun_stampede == 'Y':
                comp_env = self.comp_env
                analysis_setup_info = self.get_analysis_setup_info(analysis_case)
                allocation_name = analysis_setup_info['allocation_name']
                comp_env_rerun = input(f'Enter comp_env to re-run (default = \'{comp_env}\'): ').strip(
                    '"\'') or comp_env
                self.comp_env = comp_env_rerun
                n_batch = int(input('Enter number of batches to re-run (default = 1): ').strip('"\'') or '1')
                n_job = int(input('Enter number of nodes per batch to re-run (default = 1): ').strip('"\'') or '1')
                pool_size_factor = input('Enter pool_size_factor to re-run (default = 1): ').strip('"\'') or '1'
                if '/' in pool_size_factor:
                    frac = pool_size_factor.split('/')
                    numerator = float(frac[0])
                    denominator = float(frac[1])
                    pool_size_factor = numerator / denominator
                else:
                    pool_size_factor = float(pool_size_factor)
                run_time_limit = input('Enter run_time_limit to re-run (default = \'01:00:00\'): ').strip(
                    '"\'') or '01:00:00'
                allocation_name = input(f'Enter allocation_name to re-run (default = \'{allocation_name}\'): ').strip(
                    '"\'') or allocation_name
                self.setup_non_local_run(analysis_case, n_batch, n_job, pool_size_factor, run_time_limit,
                                         allocation_name)
                if stage:
                    self.restage(analysis_case, pool_size, gm_database_dir_path)

        return

    # ------------------------------------------------------------------------------------------------------------------

    def recheck_status(self, analysis_case, **kwargs):
        print_status_file = kwargs.get('print_status_file', True)

        setup_dir_path = self.get_setup_dir_path()
        results_dir_path = self.get_results_dir_path()

        analysis_list, _, n_sim = self.get_analysis_info(analysis_case, 'analysis_list_original.txt', 1, 1, 1, 1)

        conv_count = 0

        write_lines = list()

        for i_sim in range(n_sim):
            run_case = analysis_list[i_sim, :]
            for_which = list(run_case.astype(str))
            run_case_str = ' '.join(for_which)
            target_work_dir_path = Utility.get_path(results_dir_path, *for_which)
            status = 'FAIL'
            nltha_status_file_path = Utility.get_path(target_work_dir_path, 'NLTHA_STATUS.txt')
            if os.path.isfile(nltha_status_file_path):
                with open(nltha_status_file_path, 'r') as nltha_status_fid:
                    status = nltha_status_fid.readline().strip()

            if status == 'SUCCESS':
                conv_count += 1
                if print_status_file:
                    write_lines.append(f"{run_case_str} : CONVERGED and COMPLETED in unknown seconds\n")
            else:
                if print_status_file:
                    write_lines.append(f"{run_case_str} : FAILED in unknown seconds\n")

        if print_status_file:
            status_file_path = Utility.get_path(setup_dir_path, analysis_case, 'final_status.txt')
            fid = open(status_file_path, 'w')
            fid.writelines(write_lines)
            fid.close()

        return conv_count, n_sim

    # ------------------------------------------------------------------------------------------------------------------

    def restage(self, analysis_case, pool_size, gm_database_dir_path, **kwargs):
        analysis_list = kwargs.get('analysis_list', None)
        ai_end = kwargs.get('ai_end', 0.7)

        structure = self.structure
        results_dir_path = self.get_results_dir_path().replace(self.analysis_type, 'Prelim_Analysis')

        if analysis_list is None:
            analysis_list, _, _ = self.get_analysis_info(analysis_case, 'analysis_list.txt', 1, 1, 1, 1)

        if analysis_list.size == 0:
            return

        model_fetch_list = analysis_list[:, :Structure.get_dir_level_index('Model_Realization_Num') + 1]
        gm_fetch_list = analysis_list[:, :Structure.get_dir_level_index('Hazard_Level_Num') + 1]

        prelim_analysis_list = np.unique(model_fetch_list, axis=0)

        to_run = partial(NLTHA.run_parallel_prelim_analysis, results_dir_path, self.comp_env,
                         structure.structural_analysis_platform.run_prelim_analysis)
        to_pass = [[prelim_analysis_list[analysis_index - 1, :]
                    for analysis_index in range(1, prelim_analysis_list.shape[0] + 1)], ]

        if pool_size > 1:
            my_pool = Pool(pool_size)
            my_pool.map(to_run, *to_pass)
            my_pool.close()
            my_pool.join()
            my_pool.clear()
        else:
            list(map(to_run, *to_pass))

        for analysis_index in range(1, prelim_analysis_list.shape[0] + 1):
            run_case = prelim_analysis_list[analysis_index - 1, :]
            for_which = list(run_case.astype(str))
            structure.write_damping_file(for_which)

        self.stage(pool_size, analysis_list.astype(int), model_fetch_list.astype(int),
                   gm_fetch_list.astype(int), gm_database_dir_path, ai_end=ai_end)

        return

    # ------------------------------------------------------------------------------------------------------------------

    def find_anomalous_run_cases(self, analysis_case, threshold_list, **kwargs):
        save_analysis_list = kwargs.get('save_analysis_list', False)
        edp_list = self.edp_list
        analysis_list, analysis_range, _ = self.get_analysis_info(analysis_case,
                                                                  'analysis_list_original.txt', 1, 1, 1, 1)
        design_num_list = list(np.unique(analysis_list[:, 1]).astype(str))
        run_cases = np.array([]).reshape(0, len(Structure.get_all_dir_levels()))
        for itr in range(len(edp_list)):
            edp = edp_list[itr]
            for design_num in design_num_list:
                run_cases = np.vstack([run_cases,
                                       edp.find_anomalous_run_cases([analysis_case, design_num],
                                                                    threshold_list[itr])])

        if run_cases.shape[0] == 0:
            run_cases = run_cases.reshape(0, len(Structure.get_all_dir_levels()))
        else:
            run_cases = np.unique(run_cases, axis=0)

        if save_analysis_list:
            save_path = Utility.get_path(self.get_setup_dir_path(), analysis_case, 'analysis_list_anomalous.txt')
            Utility.save_array_as_text_file(save_path, run_cases, fmt='%d',
                                            header=f'Total analyses to run = {run_cases.shape[0]:d}')

        return run_cases

    # ------------------------------------------------------------------------------------------------------------------

    def extract_single_run(self, for_which, num_modes, gm_database_dir_path, **kwargs):

        recorder_info = kwargs.get('recorder_info', {})
        addnl_rec_file_open_mode = kwargs.get('addnl_rec_file_open_mode', 'w')
        ai_end = kwargs.get('ai_end', 0.7)

        structure = self.structure
        name = structure.name
        model_files_path = structure.model_files_path

        for_which1 = Structure.get_modified_for_which(for_which, 'Model_Realization_Num')
        for_which2 = Structure.get_modified_for_which(for_which, 'Hazard_Level_Num')

        base_dir_path = os.path.dirname(model_files_path)
        analysis_platform = structure.structural_analysis_platform.__class__.__name__
        extract_dest_path = Utility.get_path(base_dir_path, f'{name}_extract_{analysis_platform}')
        if not os.path.isdir(extract_dest_path):
            os.makedirs(extract_dest_path)

        # PrelimAnalysis extract
        NLTHA.run_parallel_prelim_analysis(
            self.get_results_dir_path().replace(self.analysis_type, 'Prelim_Analysis'),
            'local',
            structure.structural_analysis_platform.run_prelim_analysis,
            np.array(for_which1).astype(int))
        structure.extract_prelim_analysis(for_which, extract_dest_path, num_modes)

        # NLTHA extract
        self.stage(1,
                   np.array([for_which]).astype(int),
                   np.array([for_which1]).astype(int),
                   np.array([for_which2]).astype(int),
                   None,
                   ai_end=ai_end)

        # GM files fetch
        gm_file_names, scale_fac = self.im.get_gm_rec_info(for_which)
        Utility.fetch_selected_ground_motion_records(gm_database_dir_path, gm_file_names, scale_fac, extract_dest_path, ai_end)
        structure.extract_nltha(for_which, extract_dest_path)

        # ModelFiles extract
        stuff_in_model_files_path = os.listdir(model_files_path)
        [shutil.copy(Utility.get_path(model_files_path, file), extract_dest_path)
         for file in stuff_in_model_files_path if
         os.path.isfile(Utility.get_path(model_files_path, file)) and file.endswith(structure.structural_analysis_platform.file_ext)]

        # AdditionalRecorders
        structure.structural_analysis_platform.write_additional_recorders(recorder_info, extract_dest_path, addnl_rec_file_open_mode)

        # Leave a bat file for easy execution
        structure.structural_analysis_platform.bat_file(extract_dest_path)
        return

    # ------------------------------------------------------------------------------------------------------------------

    def insert_extracted_run(self, for_which):
        structure = self.structure
        name = structure.name
        model_files_path = structure.model_files_path

        base_dir_path = os.path.dirname(model_files_path)
        analysis_platform = structure.structural_analysis_platform.__class__.__name__
        src_path = Utility.get_path(base_dir_path, f'{name}_extract_{analysis_platform}')

        results_dir_path = self.get_results_dir_path()
        work_dir_path = Utility.get_path(results_dir_path, *for_which)

        edp_tags = [edp.tag for edp in self.edp_list]
        for tag in edp_tags:
            stuff_in_edp_results = os.listdir(Utility.get_path(src_path, f'EDP_{tag}_Results'))
            dest_path = Utility.get_path(work_dir_path, f'EDP_{tag}_Results')
            if not os.path.isdir(dest_path):
                os.makedirs(dest_path)
            [shutil.copy(Utility.get_path(src_path, f'EDP_{tag}_Results', file),
                         dest_path)
             for file in stuff_in_edp_results if
             os.path.isfile(Utility.get_path(src_path, f'EDP_{tag}_Results', file))]

        with open(Utility.get_path(work_dir_path, 'NLTHA_STATUS.txt'), 'w') as fid:
            fid.write('SUCCESS\n')

        analysis_case = for_which[Structure.get_dir_level_index('Analysis_Case')]
        final_status_file_path = Utility.get_path(self.get_setup_dir_path(), analysis_case, 'final_status.txt')
        if os.path.isfile(final_status_file_path):
            with open(final_status_file_path, 'r') as fid:
                lines = fid.readlines()
            for_which_str = ' '.join(for_which)
            inds = [i for i in range(len(lines)) if
                    lines[i].startswith(f"{for_which_str} : FAILED in ")
                    ]
            for ind in inds:
                lines[ind] = f"{for_which_str} : CONVERGED and COMPLETED in unknown seconds\n"
            if len(inds) == 0:
                lines.append(f"{for_which_str} : CONVERGED and COMPLETED in unknown seconds\n")
            with open(final_status_file_path, 'w') as fid:
                fid.writelines(lines)

        [edp.collect(analysis_case, 1, analysis_list=np.array([[int(item) for item in for_which]]))
         for edp in self.edp_list]

        return

    # ------------------------------------------------------------------------------------------------------------------

    def stampede_upload_pre(self, analysis_case, design_num_list, haz_lev_list, username, allocation_name, n_batch, n_job, upload_to, **kwargs):
        local_bash_path = kwargs.get('local_bash_path', '')
        comp_env = kwargs.get('comp_env', 'local')
        new_file = kwargs.get('new_file', False)
        keep_file = kwargs.get('keep_file', False)
        nltha_include = kwargs.get('nltha_include', False)
        extract = kwargs.get('extract', True)
        stampede_comp_env = kwargs.get('stampede_comp_env', 'stampede_knl')
        nltha_setup_pool_size = kwargs.get('nltha_setup_pool_size', 50)
        gm_database_dir_path = kwargs.get('gm_database_dir_path', Utility.get_path(os.path.dirname(upload_to), 'NGA_PEER_EQ_DB'))
        nltha_run_pool_size_factor = kwargs.get('nltha_run_pool_size_factor', 0.85)

        structure = self.structure
        include_list = [f"GMS_Results/{analysis_case}/{{{','.join(design_num_list)}}}",
                        f"GMS_Setup/{analysis_case}",
                        f"Prelim_Analysis_Results/{analysis_case}/{{{','.join(design_num_list)}}}",
                        f"Prelim_Analysis_Setup/{analysis_case}",
                        f"PSHA_Results/{analysis_case}/{{{','.join(design_num_list)}}}",
                        f"PSHA_Setup/{analysis_case}",
                        f"Random_Model_Param_Vals/{analysis_case}/{{{','.join(design_num_list)}}}"]
        if nltha_include:
            include_list.append(f"NLTHA_EDP/{analysis_case}/{{{','.join(design_num_list)}}}")
            include_list.append(f"NLTHA_Setup/{analysis_case}")

        if len(design_num_list) == 1:
            include_list = [item.replace('{', '').replace('}', '') for item in include_list]
        temp_file_name = structure.tar_work_dir(local_bash_path=local_bash_path, comp_env=comp_env, new_file=True,
                                                flags=['cvf'],
                                                include_list=include_list)

        shutil.move(Utility.get_path(structure.model_work_dir_path, temp_file_name),
                    structure.model_files_path)

        exec_commands = [f'tar -uvf {temp_file_name} * --exclude={structure.name}_Work_Dir*.tar']

        Utility.shell_exec('APPENDTAR', structure.model_files_path, exec_commands,
                           comp_env=comp_env, local_bash_path=local_bash_path)

        stampede_files_dir = Utility.get_path(os.getcwd(), 'Stampede_Files')
        os.makedirs(stampede_files_dir, exist_ok=True)
        objects_dict = {
            'analysis_case': analysis_case,
            'nltha_setup_pool_size': nltha_setup_pool_size,
            'gm_database_dir_path': gm_database_dir_path,
            'design_num_list': design_num_list,
            'haz_lev_list': haz_lev_list,
            'allocation_name': allocation_name,
            'n_batch': n_batch,
            'n_job': n_job,
            'nltha_run_pool_size_factor': nltha_run_pool_size_factor,
            'nltha': self.get_nltha_for_stampede_upload(upload_to, stampede_comp_env),
        }
        Utility.pickle_dump_dict(Utility.get_path(stampede_files_dir, 'objects.pickle'), objects_dict)

        with open(Utility.get_path(stampede_files_dir, 'setup_nltha.py'), 'w') as fid:
            fid.write('from pypbee.utility import Utility\n')
            fid.write('\n')
            fid.write(f"objects = Utility.pickle_load_dict(\'{Utility.get_path(upload_to, 'objects.pickle')}\')\n")
            fid.write('analysis_case = objects[\'analysis_case\']\n')
            fid.write('nltha_setup_pool_size = objects[\'nltha_setup_pool_size\']\n')
            fid.write('gm_database_dir_path = objects[\'gm_database_dir_path\']\n')
            fid.write('design_num_list = objects[\'design_num_list\']\n')
            fid.write('haz_lev_list = objects[\'haz_lev_list\']\n')
            fid.write('allocation_name = objects[\'allocation_name\']\n')
            fid.write('n_batch = objects[\'n_batch\']\n')
            fid.write('n_job = objects[\'n_job\']\n')
            fid.write('nltha_run_pool_size_factor = objects[\'nltha_run_pool_size_factor\']\n')
            fid.write('nltha = objects[\'nltha\']\n')
            fid.write('kwargs = {\n')
            fid.write('    \'run_time_limit\': \'01:00:00\',\n')
            fid.write('    \'allocation_name\': allocation_name,\n')
            fid.write('    \'n_batch\': n_batch,\n')
            fid.write('    \'n_job\': n_job,\n')
            fid.write('    \'pool_size_factor\': nltha_run_pool_size_factor,\n')
            fid.write('}\n')
            fid.write('nltha.setup(\n')
            fid.write('    analysis_case,\n')
            fid.write('    nltha_setup_pool_size,\n')
            fid.write('    design_num_list,\n')
            fid.write('    haz_lev_list,\n')
            fid.write('    gm_database_dir_path,\n')
            fid.write('    **kwargs\n')
            fid.write(')\n')

        with open(Utility.get_path(stampede_files_dir, 'run_on_stampede.py'), 'w') as fid:
            fid.write('from pypbee.utility import Utility\n')
            fid.write('\n')
            fid.write(f"objects = Utility.pickle_load_dict(\'{Utility.get_path(upload_to, 'objects.pickle')}\')\n")
            fid.write('nltha = objects[\'nltha\']\n')
            fid.write('\n')
            fid.write('\n')
            fid.write('def run_on_stampede_nltha(analysis_case, n_batch, n_job, node_index, batch_index, pool_size):\n')
            fid.write('    nltha.run(\n')
            fid.write('        analysis_case,\n')
            fid.write('        pool_size,\n')
            fid.write('        n_batch=n_batch,\n')
            fid.write('        n_job=n_job,\n')
            fid.write('        node_index=node_index,\n')
            fid.write('        batch_index=batch_index\n')
            fid.write('    )\n')

        shutil.move(Utility.get_path(structure.model_files_path, temp_file_name), stampede_files_dir)

        Utility.shell_exec('APPENDTAR', stampede_files_dir, exec_commands, comp_env=comp_env, local_bash_path=local_bash_path)

        shutil.move(Utility.get_path(stampede_files_dir, temp_file_name), Utility.get_path(os.getcwd()))

        shutil.rmtree(stampede_files_dir)

        Utility.shell_exec('COMPRESS', '',
                           ['echo "compressing ..."', f'gzip {temp_file_name}'],
                           comp_env=comp_env, local_bash_path=local_bash_path)

        if new_file:
            new_file_name = f'{structure.name}_nltha_{analysis_case}_stampede_upload_pre_{Utility.get_curr_time_stamp()}.tar'
        else:
            new_file_name = f'{structure.name}_nltha_{analysis_case}_stampede_upload_pre.tar'

        os.replace(temp_file_name + '.gz', new_file_name)

        src = new_file_name
        server = f"{username}@stampede2.tacc.utexas.edu"

        Utility.shell_exec('UPLOAD', '',
                           [f'scp {src} {server}:{upload_to}'],
                           comp_env=comp_env, local_bash_path=local_bash_path)

        if extract:
            Utility.shell_exec('EXTRACT', '',
                               [f"ssh -T {server} << EOF",
                                f'cd {upload_to}',
                                f'tar -xzvf {new_file_name}',
                                f'rm -rf {new_file_name}\n' if not keep_file else '',
                                'EOF'],
                               comp_env=comp_env, local_bash_path=local_bash_path)

        if not keep_file:
            os.remove(new_file_name)

        return

    # ------------------------------------------------------------------------------------------------------------------

    def get_nltha_for_stampede_upload(self, upload_to, comp_env):
        new_structure = deepcopy(self.structure)
        new_structure.model_files_path = upload_to
        new_structure.model_work_dir_path = upload_to
        new_im = deepcopy(self.im)
        new_im.structure = new_structure
        new_edps = deepcopy(self.edp_list)
        for edp in new_edps:
            edp.structure = new_structure
            edp.haz_req = {}
        return NLTHA(new_edps, new_im, comp_env=comp_env)

    def get_post_analysis_lists(self, analysis_case):
        setup_dir_path = self.get_setup_dir_path()
        final_status_file_path = Utility.get_path(setup_dir_path, analysis_case, 'final_status.txt')

        with open(final_status_file_path, 'r') as fid:
            analysis_failed = np.array(
                [np.array(line.split(':')[0].strip().split()).astype(int)
                 for line in fid.readlines() if not line.isspace() and 'FAILED' in line])

        if analysis_failed.size == 0:
            analysis_failed = analysis_failed.reshape(0, len(Structure.get_all_dir_levels()))

        analysis_ran = Utility.read_numpy_array_from_txt_file(final_status_file_path, comments=':').astype(int)

        analysis_list, _, n_sim = self.get_analysis_info(analysis_case, 'analysis_list.txt', 1, 1, 1, 1)

        if n_sim > 0:
            analysis_skipped = analysis_list[Utility.in_nd(analysis_list, analysis_ran) == 0, :]
        else:
            analysis_skipped = np.array([]).reshape(0, len(Structure.get_all_dir_levels()))

        return analysis_ran, analysis_failed, analysis_skipped

    # ------------------------------------------------------------------------------------------------------------------

    def collect_edps(self, analysis_case, pool_size, analysis_list='analysis_list_original.txt'):
        edp_list = self.edp_list
        [edp.collect(analysis_case, pool_size, analysis_list=analysis_list) for edp in edp_list]

    # ------------------------------------------------------------------------------------------------------------------

    def stampede_download_post(self, analysis_case, pool_size, **kwargs):
        local_bash_path = kwargs.get('local_bash_path', '')
        comp_env = kwargs.get('comp_env', 'local')
        new_file = kwargs.get('new_file', False)
        analysis_list = kwargs.get('analysis_list', 'analysis_list_original.txt')

        self.collect_edps(analysis_case, pool_size, analysis_list=analysis_list)
        if isinstance(analysis_list, str):
            analysis_list, _, _ = self.get_analysis_info(analysis_case, analysis_list, 1, 1, 1, 1)
        design_num_list = list(np.unique(analysis_list[:, 1]).astype(str))

        structure = self.structure
        include_list = [f"GMS_Results/{analysis_case}/{{{','.join(design_num_list)}}}",
                        f"GMS_Setup/{analysis_case}",
                        f"Prelim_Analysis_Results/{analysis_case}/{{{','.join(design_num_list)}}}",
                        f"Prelim_Analysis_Setup/{analysis_case}",
                        f"PSHA_Results/{analysis_case}/{{{','.join(design_num_list)}}}",
                        f"PSHA_Setup/{analysis_case}",
                        f"Random_Model_Param_Vals/{analysis_case}/{{{','.join(design_num_list)}}}",
                        f"NLTHA_EDP/{analysis_case}/{{{','.join(design_num_list)}}}",
                        f"NLTHA_Setup/{analysis_case}"]
        if len(design_num_list) == 1:
            include_list = [item.replace('{', '').replace('}', '') for item in include_list]
        temp_file_name = structure.tar_work_dir(local_bash_path=local_bash_path, comp_env=comp_env, new_file=True,
                                                flags=['czvf'],
                                                include_list=include_list)

        if new_file:
            new_file_name = f'{structure.name}_nltha_case_{analysis_case}_' \
                            f'stampede_download_post_{Utility.get_curr_time_stamp()}.tar'
        else:
            new_file_name = f'{structure.name}_nltha_case_{analysis_case}_stampede_download_post.tar'

        os.replace(temp_file_name, new_file_name)
        return

    # ------------------------------------------------------------------------------------------------------------------

    def get_converged_analysis_count(self, analysis_case):
        conv_count, n_sim = self.recheck_status(analysis_case, print_status_file=False)
        print(f'{conv_count}/{n_sim} done!')
        return conv_count

    # ------------------------------------------------------------------------------------------------------------------
