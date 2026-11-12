# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 2026

@author: Angshuman Deb
"""

from pypbee.edp import EDP
from pypbee.structure import Structure
from pypbee.utility import Utility


class MaxDisp(EDP):

    ####################################################################################################################
    # Constructor
    ####################################################################################################################

    def __init__(self, node_tag, dof, structure, tag, recorder_file_storage, **kwargs):
        haz_req = kwargs.get('haz_req', {})
        small_value = kwargs.get('small_value', 1.0e-10)
        large_value = kwargs.get('large_value', 1.0e2)
        super().__init__(structure, tag, recorder_file_storage, haz_req, small_value, large_value)
        self.node_tag = node_tag
        self.dof = dof
        self.resp = 'disp'
        self.normalize_with = 1.0

    ####################################################################################################################
    # Public functionalities (instance methods)
    ####################################################################################################################

    def write_evaluate_edp_and_helpers(self):
        model_files_path = self.structure.model_files_path
        tag = self.tag
        resp = self.resp

        lines = [
            'import numpy as np\n\n\n',
            f'def get_{resp}_vals(file_path):\n',
            '    return np.loadtxt(file_path)[:, 1]\n',
        ]
        with open(Utility.get_path(model_files_path, f'get_{resp}_vals.py'), 'w') as fid:
            fid.writelines(lines)

        lines = [
            'import numpy as np\n\n\n',
            f'def evaluate_edp_{tag}(x):\n',
            '    return np.max(np.abs(x))\n',
        ]
        with open(Utility.get_path(model_files_path, f'evaluate_edp_{tag}.py'), 'w') as fid:
            fid.writelines(lines)
        return

    def generate_recorder(self, for_which, **kwargs):
        rec_gen_file_open_mode = kwargs.get('rec_gen_file_open_mode', 'a+')
        gen_rec = kwargs.get('gen_rec', True)

        structure = self.structure
        platform = structure.structural_analysis_platform
        dir_category = 'NLTHA_Results'
        dir_level_to_seek = 'Ground_Motion_Rec_Num'
        for_which = Structure.get_modified_for_which(for_which, dir_level_to_seek)
        work_dir_path = Utility.get_path(structure.model_work_dir_path, 'Work_Dir', dir_category, *for_which)

        if gen_rec:
            rec_save_dir_path = self.get_rec_save_dir_path()
            rec_file_path = Utility.get_path(rec_save_dir_path, f'{self.resp}_node_{self.node_tag}.txt')
            rec_gen_file_fid = platform.open_edp_recorder_file(work_dir_path, rec_gen_file_open_mode, rec_save_dir_path)
            rec_gen_file_fid.write(
                f"    ops.recorder('Node', '-file', '{rec_file_path}', '-time', "
                f"'-node', {self.node_tag}, '-dof', {self.dof}, '{self.resp}')\n"
            )
            rec_gen_file_fid.close()
        platform.write_evaluate_edps(work_dir_path, self, for_which, gen_rec, rec_gen_file_open_mode)
        return

    def get_edp_strings(self, for_which):
        return [f'node_{self.node_tag}']
