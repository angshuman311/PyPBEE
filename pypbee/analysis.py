# -*- coding: utf-8 -*-
"""
Abstract class : Analysis

@author: Angshuman Deb
"""

from abc import ABC, abstractmethod
import numpy as np
from .utility import Utility
import os


class Analysis(ABC):

    stampede_knl_n = 68
    stampede_skx_n = 48

    ####################################################################################################################
    # Constructor
    ####################################################################################################################

    def __init__(self, analysis_type, structure, comp_env):
        self.analysis_type = analysis_type
        self.structure = structure
        self.comp_env = comp_env

    ####################################################################################################################
    # Abstract methods
    ####################################################################################################################

    @abstractmethod
    def setup(self, *args, **kwargs):
        pass

    # ------------------------------------------------------------------------------------------------------------------

    @abstractmethod
    def run(self, *args, **kwargs):
        pass

    # ------------------------------------------------------------------------------------------------------------------

    @abstractmethod
    def wrap_up(self, *args, **kwargs):
        pass

    ####################################################################################################################
    # Public functionalities (instance methods)
    ####################################################################################################################

    def get_setup_dir_path(self):
        setup_dir_path = Utility.get_path(self.structure.model_work_dir_path,
                                          'Work_Dir', self.analysis_type + '_Setup')
        return setup_dir_path

    # ------------------------------------------------------------------------------------------------------------------

    def get_results_dir_path(self):
        results_dir_path = Utility.get_path(self.structure.model_work_dir_path,
                                            'Work_Dir', self.analysis_type + '_Results')
        return results_dir_path

    # ------------------------------------------------------------------------------------------------------------------

    def get_analysis_setup_info(self, analysis_case):
        analysis_setup_info = Utility.pickle_load_dict(
            Utility.get_path(
                self.get_setup_dir_path(), analysis_case, 'analysis_setup_info.pickle'
            )
        )
        return analysis_setup_info

    # ------------------------------------------------------------------------------------------------------------------

    def set_analysis_setup_info(self, analysis_case, analysis_setup_info):
        Utility.pickle_dump_dict(
            Utility.get_path(
                self.get_setup_dir_path(), analysis_case, 'analysis_setup_info.pickle'
            ), analysis_setup_info
        )
        return

    # ------------------------------------------------------------------------------------------------------------------

    def generate_batch_file(self, analysis_case, n_batch, n_job, pool_size_factor):

        structure = self.structure
        name = structure.name
        model_work_dir_path = structure.model_work_dir_path
        analysis_type = self.analysis_type
        comp_env = self.comp_env
        slurm_file_name = f'slurm_file_{name}_{analysis_type.lower()}_case_{analysis_case}'
        batch_file_name_prefix = f'batch_file_{name}_{analysis_type.lower()}_case_{analysis_case}'
        err_file_name_prefix = f'err_file_{name}_{analysis_type.lower()}_case_{analysis_case}'
        out_file_name_prefix = f'out_file_{name}_{analysis_type.lower()}_case_{analysis_case}'

        pool_size = int(getattr(self, comp_env + '_n') * pool_size_factor)
        pool_size = max(1, pool_size)

        for iter_batch in range(1, n_batch + 1):
            with open(f'{model_work_dir_path}/{batch_file_name_prefix}_batch_{iter_batch}', 'w') as batch_fid:
                batch_fid.write(f'#!/bin/bash\n')

                write_line = f'if [ ! -d "{model_work_dir_path}/Work_Dir/EO_Files" ]; '
                write_line += f'then mkdir "{model_work_dir_path}/Work_Dir/EO_Files"; fi;\n'
                batch_fid.write(write_line)

                batch_fid.write(f'for i in {{{n_job * (iter_batch - 1) + 1}..{n_job * iter_batch}}}; do\n')

                write_line = f'sbatch -J run$i '
                write_line += f'-e "{model_work_dir_path}/Work_Dir/EO_Files/{err_file_name_prefix}_job_$i" '
                write_line += f'-o "{model_work_dir_path}/Work_Dir/EO_Files/{out_file_name_prefix}_job_$i" '
                write_line += f'{slurm_file_name} $i {iter_batch} {pool_size}\n'
                batch_fid.write(write_line)

                batch_fid.write(f'done\n')
            os.system(f'chmod +x "{model_work_dir_path}/{batch_file_name_prefix}_batch_{iter_batch}"')
        return slurm_file_name

    # ------------------------------------------------------------------------------------------------------------------

    def generate_slurm_file(self, analysis_case, slurm_file_name, run_time_limit, allocation_name, n_batch, n_job):

        structure = self.structure
        model_work_dir_path = structure.model_work_dir_path
        analysis_type = self.analysis_type
        comp_env = self.comp_env

        with open(f'{model_work_dir_path}/{slurm_file_name}', 'w') as slurm_file_fid:
            slurm_file_fid.write(f'#!/bin/bash\n')
            slurm_file_fid.write(f'#SBATCH -N 1\n')
            n = getattr(self, comp_env + '_n')
            slurm_file_fid.write(f'#SBATCH -n {n}\n')
            if comp_env == 'stampede_knl':
                slurm_file_fid.write(f'#SBATCH -p normal\n')
            elif comp_env == 'stampede_skx':
                slurm_file_fid.write(f'#SBATCH -p skx-normal\n')
            slurm_file_fid.write(f'#SBATCH -t {run_time_limit}\n')
            slurm_file_fid.write(f'#SBATCH -A {allocation_name}\n')
            slurm_file_fid.write(f'\n')
            slurm_file_fid.write(f'module load python3\n')

            python_base_path = Utility.get_path(os.getcwd())
            write_line = f'python3 -c "import sys; sys.path.insert(1, \'{python_base_path}\'); '
            write_line += f'from run_on_stampede import run_on_stampede_{analysis_type.lower()} as run_on_stampede; '
            write_line += f'run_on_stampede(\'{analysis_case}\', {n_batch}, {n_job}, $1, $2, $3)"\n'
            slurm_file_fid.write(write_line)
        return

    # ------------------------------------------------------------------------------------------------------------------

    def get_analysis_info(self, analysis_case, list_file_name, n_batch, n_job, node_index, batch_index):

        list_file_path = Utility.get_path(self.get_setup_dir_path(), analysis_case, list_file_name)

        if node_index > n_job:
            if node_index % n_job == 0:
                node_index = n_job
            else:
                node_index = node_index % n_job

        analysis_list = Utility.read_numpy_array_from_txt_file(list_file_path, skiprows='find')

        if analysis_list.size != 0:
            tot_num_analysis = analysis_list.shape[0]
        else:
            tot_num_analysis = 0

        n_sim = tot_num_analysis

        analysis_per_batch = np.zeros(n_batch)
        analysis_per_batch[0: analysis_per_batch.size - 1] = np.floor(n_sim / n_batch)
        analysis_per_batch[analysis_per_batch.size - 1] = n_sim - (n_batch - 1) * np.floor(n_sim / n_batch)

        analysis_this_batch = analysis_per_batch[batch_index - 1]

        analysis_per_job = np.zeros(n_job)
        analysis_per_job[0: analysis_per_job.size - 1] = np.floor(analysis_this_batch / n_job)
        analysis_per_job[analysis_per_job.size - 1] = analysis_this_batch - (n_job - 1) * np.floor(
            analysis_this_batch / n_job)

        analysis_range_i = sum(analysis_per_batch[0: batch_index - 1]) + sum(analysis_per_job[0: node_index - 1]) + 1
        analysis_range_f = sum(analysis_per_batch[0: batch_index - 1]) + sum(analysis_per_job[0: node_index])

        analysis_range = np.arange(analysis_range_i, analysis_range_f + 1)

        return analysis_list.astype(int), analysis_range.astype(int), len(analysis_range)

    # ------------------------------------------------------------------------------------------------------------------

    def setup_non_local_run(self, analysis_case, n_batch, n_job,
                            pool_size_factor, run_time_limit, allocation_name):
        # Write slurm job file and batch files:
        slurm_file_name = self.generate_batch_file(analysis_case, n_batch, n_job, pool_size_factor)
        self.generate_slurm_file(analysis_case, slurm_file_name, run_time_limit, allocation_name, n_batch, n_job)
        return

    # ------------------------------------------------------------------------------------------------------------------
