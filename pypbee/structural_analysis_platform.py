# -*- coding: utf-8 -*-
"""
Abstract class : StructuralAnalysisPlatform
Derived class : OpenSeesTcl
Derived class : OpenSeesPy

@author: Angshuman Deb
"""

import json
import os.path
import platform
import subprocess
import numpy as np
from .utility import Utility
from abc import ABC, abstractmethod


class StructuralAnalysisPlatform(ABC):
    def __init__(self, model_files_path):
        self.model_files_path = Utility.get_path(model_files_path)
        self.file_ext = None
        self.nltha_rec_dir_name = 'NLTHA_Recorder_Files'
        self.model_info_dir_name = 'Model_Info_Files'

    @abstractmethod
    def bat_file(self, dest_path):
        pass

    @abstractmethod
    def write_additional_recorders(self, recorder_info, dest_path, file_open_mode):
        pass

    @abstractmethod
    def open_edp_recorder_file(self, work_dir_path, rec_gen_file_open_mode, rec_save_dir_path):
        pass

    @abstractmethod
    def write_evaluate_edp_and_helpers(self, edp_tag, max_what, resp, resp_num):
        pass

    @abstractmethod
    def write_evaluate_edps(self, work_dir_path, edp, for_which, gen_rec, mode):
        pass

    @abstractmethod
    def write_run_prelim_analysis(self, work_dir_path, num_modes):
        pass

    @abstractmethod
    def write_run_nltha(self, target_work_dir_path, work_dir_path_prelim_analysis, **kwargs):
        pass

    @abstractmethod
    def get_prelim_analysis_files(self):
        pass

    @abstractmethod
    def get_nltha_files(self):
        pass

    @abstractmethod
    def run_prelim_analysis(self, target_work_dir_path, comp_env):
        pass

    @abstractmethod
    def run_nltha(self, target_work_dir_path, comp_env):
        pass

    @abstractmethod
    def run_file(self, target_work_dir_path, run_file_name, comp_env):
        pass

    @abstractmethod
    def perform_mat_test(self, mat_data_dir_path, mat_tag, input_data, num_incr):
        pass

    @staticmethod
    @abstractmethod
    def write_model_param_vals(work_dir_path, name_value_pairs):
        pass

    @staticmethod
    @abstractmethod
    def write_model_attributes(work_dir_path, name_value_pairs):
        pass

    @staticmethod
    @abstractmethod
    def write_damping_file(damping_model, work_dir_path, **kwargs):
        pass

    @staticmethod
    @abstractmethod
    def get_stress_strain_recorder_command(recorder_file_path, ele_tag, sec_num, y, z, mat_tag):
        pass

    @staticmethod
    @abstractmethod
    def get_elem_deformation_recorder_command(recorder_file_path, ele_tag):
        pass

    @staticmethod
    @abstractmethod
    def get_recorder_files_to_read(rec_gen_file_path):
        pass


class OpenSeesTcl(StructuralAnalysisPlatform):
    def __init__(self, model_files_path, local_opensees_path):
        super().__init__(model_files_path)
        self.local_opensees_path = Utility.get_path(local_opensees_path)
        self.file_ext = '.tcl'

    def bat_file(self, dest_path):
        # OpenSees.bat
        if platform.system() == 'Windows':
            with open(Utility.get_path(dest_path, 'OpenSees.bat'), 'w') as fid:
                fid.write(f'"{self.local_opensees_path}"')

    def write_additional_recorders(self, recorder_info, dest_path, file_open_mode):
        element_response_list = recorder_info.get('element_response_list', [])
        element_section_response_list = recorder_info.get('element_section_response_list', [])
        element_section_fiber_response_list = recorder_info.get('element_section_fiber_response_list', [])
        node_response_list = recorder_info.get('node_response_list', [])
        addnl_file_list = recorder_info.get('addnl_file_list', [])

        if 'a' in file_open_mode:
            with open(Utility.get_path(dest_path, 'generate_additional_recorders.tcl'), 'r') as fid:
                lines = fid.readlines()
                if not lines[-1].endswith('\n'):
                    lines[-1] = lines[-1] + '\n'
        else:
            lines = []

        nltha_dir = self.nltha_rec_dir_name
        for item in element_response_list:
            ele_tag, what = item
            file_name = Utility.get_path(nltha_dir, f'element_response_elem_{ele_tag}_resp_{what}.txt')
            recorder_line = f'recorder Element -file "{file_name}" -time -ele {ele_tag} {what}\n'
            lines.append(recorder_line)

        for item in element_section_response_list:
            ele_tag, sec_num, what = item
            file_name = Utility.get_path(nltha_dir, f'element_section_response_elem_{ele_tag}_secnum_{sec_num}_resp_{what}.txt')
            recorder_line = f'recorder Element -file "{file_name}" -time -ele {ele_tag} section {sec_num} {what}\n'
            lines.append(recorder_line)

        for item in element_section_fiber_response_list:
            ele_tag, sec_num, fib_num, y, z, mat_tag = item
            file_name = Utility.get_path(nltha_dir, f'element_section_fiber_response_elem_{ele_tag}_secnum_{sec_num}_fib_{fib_num}_mat_{mat_tag}_resp_stressStrain.txt')
            recorder_line = f'recorder Element -file "{file_name}" -time -ele {ele_tag} section {sec_num} fiber {y} {z} {mat_tag} stressStrain\n'
            lines.append(recorder_line)

        for item in node_response_list:
            node_tag, dof, what = item
            file_name = Utility.get_path(nltha_dir, f'node_response_node_{node_tag}_dof_{dof}_resp_{what}.txt')
            recorder_line = f'recorder Node -file "{file_name}" -time -node {node_tag} -dof {dof} {what}\n'
            lines.append(recorder_line)

        for item in addnl_file_list:
            item = item.split('.')[0] + self.file_ext
            lines.append(f'source "{item}"\n')

        with open(Utility.get_path(dest_path, 'generate_additional_recorders.tcl'), 'w') as fid:
            fid.writelines(lines)

    def open_edp_recorder_file(self, work_dir_path, rec_gen_file_open_mode, rec_save_dir_path):
        rec_gen_file_fid = open(Utility.get_path(work_dir_path, 'generate_edp_recorders' + self.file_ext), rec_gen_file_open_mode)
        new_recorder_line = f'file mkdir "{rec_save_dir_path}";\n'
        rec_gen_file_fid.write(new_recorder_line)
        return rec_gen_file_fid

    def write_evaluate_edp_and_helpers(self, edp_tag, max_what, resp, resp_num):
        model_files_path = self.model_files_path
        lines = list()
        lines.append('proc max_list list {\n')
        lines.append('    set index 0\n')
        lines.append('    set max_index $index\n')
        lines.append('    set max_val [lindex $list 0]\n')
        lines.append('    foreach val $list {\n')
        lines.append('        if {$val > $max_val} {\n')
        lines.append('            set max_index $index\n')
        lines.append('            set max_val $val\n')
        lines.append('        }\n')
        lines.append('        incr index\n')
        lines.append('    }\n')
        lines.append('    return [list $max_val $max_index]\n')
        lines.append('}\n')
        max_list_file_path = Utility.get_path(model_files_path, 'max_list.tcl')
        with open(max_list_file_path, 'w') as fid:
            fid.writelines(lines)

        if max_what in ['compression', 'range']:
            lines = list()
            lines.append('proc min_list list {\n')
            lines.append('    set index 0\n')
            lines.append('    set min_index $index\n')
            lines.append('    set min_val [lindex $list 0]\n')
            lines.append('    foreach val $list {\n')
            lines.append('        if {$val < $min_val} {\n')
            lines.append('            set min_index $index\n')
            lines.append('            set min_val $val\n')
            lines.append('        }\n')
            lines.append('        incr index\n')
            lines.append('    }\n')
            lines.append('    return [list $min_val $min_index]\n')
            lines.append('}\n')
            min_list_file_path = Utility.get_path(model_files_path, 'min_list.tcl')
            with open(min_list_file_path, 'w') as fid:
                fid.writelines(lines)

        lines = list()
        lines.append(f'proc get_{resp}_vals file_path {{\n')
        lines.append(f'    set data_fid [open $file_path "r"]\n')
        lines.append(f'    set data [read $data_fid]\n')
        lines.append(f'    close $data_fid\n')
        lines.append(f'    set data_new [split $data "\\n"]\n')
        lines.append(f'    set {resp}_vals {{}}\n')
        lines.append(f'    for {{set k 0}} {{$k <= [expr [llength $data_new] - 2]}} {{incr k 1}} {{\n')
        lines.append(f'        set data_t [lindex $data_new $k]\n')
        lines.append(f'        lappend {resp}_vals [lindex $data_t [expr {resp_num}]]\n')
        lines.append(f'    }}\n')
        lines.append(f'    return ${resp}_vals\n')
        lines.append(f'}}\n')
        get_resp_vals_file_path = Utility.get_path(model_files_path, f'get_{resp}_vals.tcl')
        with open(get_resp_vals_file_path, 'w') as fid:
            fid.writelines(lines)

        lines = list()
        if max_what == 'tension':
            lines.append(f'proc evaluate_edp_{edp_tag} {resp}_vals {{\n')
            lines.append(f'    set temp [max_list ${resp}_vals]\n')
            lines.append(f'    return [lindex $temp 0]\n')
            lines.append(f'}}\n')
        if max_what == 'compression':
            lines.append(f'proc evaluate_edp_{edp_tag} {resp}_vals {{\n')
            lines.append(f'    set temp [min_list ${resp}_vals]\n')
            lines.append(f'    return [expr -([lindex $temp 0])]\n')
            lines.append(f'}}\n')
        if max_what == 'range':
            lines.append(f'proc evaluate_edp_{edp_tag} {resp}_vals {{\n')
            lines.append(f'    set temp [max_list ${resp}_vals]\n')
            lines.append(f'    set max_{resp} [lindex $temp 0]\n')
            lines.append(f'    set {resp}_vals [lrange ${resp}_vals [lindex $temp 1] end]\n')
            lines.append(f'    set temp [min_list ${resp}_vals]\n')
            lines.append(f'    set min_{resp} [lindex $temp 0]\n')
            lines.append(f'    return [expr $max_{resp} - $min_{resp}]\n')
            lines.append(f'}}\n')

        eval_file_path = Utility.get_path(model_files_path, f'evaluate_edp_{edp_tag}.tcl')
        with open(eval_file_path, 'w') as fid:
            fid.writelines(lines)

    def write_evaluate_edps(self, work_dir_path, edp, for_which, gen_rec, mode):
        rec_gen_file_path = Utility.get_path(work_dir_path, 'generate_edp_recorders' + self.file_ext)
        eval_edp_file_path = Utility.get_path(work_dir_path, 'evaluate_edps' + self.file_ext)

        recorder_files_to_read = OpenSeesTcl.get_recorder_files_to_read(rec_gen_file_path)
        if edp.recorder_file_storage == 'shared':
            recorder_files_to_read = [file for file in recorder_files_to_read
                                      if f'/EDP/' in file]
        else:
            recorder_files_to_read = [file for file in recorder_files_to_read
                                      if f'/EDP_{edp.tag}/' in file]

        lines = list()
        lines.append(f'source "{self.model_files_path}/evaluate_edp_{edp.tag}.tcl"\n')
        lines.append(f'file mkdir "EDP_{edp.tag}_Results"\n')
        if gen_rec:
            lines.append(f'source "{self.model_files_path}/get_{edp.resp}_vals.tcl"\n')
            lines.append(f'source "{self.model_files_path}/max_list.tcl"\n')
        if edp.max_what in ['compression', 'range']:
            lines.append(f'source "{self.model_files_path}/min_list.tcl"\n')

        edp_strings = edp.get_edp_strings(for_which)
        for edp_string in edp_strings:
            edp_var_name = f'edp_{edp.tag}_' + edp_string
            files_to_eval = [file for file in recorder_files_to_read if
                             edp_string in file]

            var_names_to_eval = [os.path.basename(file).split('.')[0] for file in files_to_eval]
            if gen_rec:
                for (var_name, file) in zip(var_names_to_eval, files_to_eval):
                    lines.append(f'set {var_name} [get_{edp.resp}_vals "{file}"]\n')

            lines.append(f'set {edp_var_name} {{0}}\n')
            for var_name in var_names_to_eval:
                lines.append(f'lappend {edp_var_name} [evaluate_edp_{edp.tag} ${var_name}]\n')

            lines.append(f'set temp [max_list ${edp_var_name}]\n')
            lines.append(f'set {edp_var_name} [lindex $temp 0]\n')
            lines.append(f'set {edp_var_name} [expr ${edp_var_name} / {edp.normalize_with}]\n')
            lines.append(f'set edp_file [open "EDP_{edp.tag}_Results/{edp_var_name}.txt" "w"]\n')
            lines.append(f'puts $edp_file ${edp_var_name}\n')
            lines.append(f'close $edp_file\n')

        with open(eval_edp_file_path, mode) as eval_file_fid:
            eval_file_fid.writelines(lines)
        return

    def write_run_prelim_analysis(self, work_dir_path, num_modes):
        model_files_path = self.model_files_path
        lines = ['wipe;\n',
                 f'set model_files_path "{model_files_path}";\n',
                 f'set model_info_directory "{self.model_info_dir_name}";\n',
                 'file mkdir $model_info_directory;\n',
                 'set write_model_files 1;\n',
                 f'set modes [expr int({num_modes})];\n',
                 'source "model_param_vals.tcl";\n',
                 'source "model_attributes.tcl";\n',
                 'source "$model_files_path/main.tcl";\n',
                 'source "$model_files_path/analysis_gravity.tcl";\n',
                 'source "$model_files_path/analysis_prelim.tcl";\n',
                 'wipe;\n']
        with open(Utility.get_path(work_dir_path, 'run_prelim_analysis.tcl'), 'w') as fid:
            fid.writelines(lines)

    def write_run_nltha(self, target_work_dir_path, work_dir_path_prelim_analysis, **kwargs):
        write_model_files = kwargs.get('write_model_files', 0)
        generate_additional_recorders = kwargs.get('generate_additional_recorders', False)
        delete_rec = kwargs.get('delete_rec', True)
        lines = [
            'wipe;\n',
            'set analysis_status_file [open "NLTHA_STATUS.txt" "w"];\n',
            'puts $analysis_status_file "FAIL";\n',
            'close $analysis_status_file;\n',
            f'set model_files_path "{self.model_files_path}";\n',
            f'set work_dir_path_prelim_analysis "{work_dir_path_prelim_analysis}";\n'
        ]
        if write_model_files == 0:
            lines += [
                'set write_model_files 0;\n',
            ]
        else:
            lines += [
                'set write_model_files 1;\n',
                f'set model_info_directory "{self.model_info_dir_name}";\n',
                'file mkdir $model_info_directory;\n',
            ]
        lines += [
            f'source "$work_dir_path_prelim_analysis/model_param_vals.tcl";\n',
            f'source "$work_dir_path_prelim_analysis/model_attributes.tcl";\n',
            'source "$model_files_path/main.tcl";\n',
            'source "$model_files_path/analysis_gravity.tcl";\n',
            'source "$work_dir_path_prelim_analysis/generate_damping.tcl";\n',
            'source "generate_edp_recorders.tcl";\n',
        ]
        if generate_additional_recorders:
            lines += [
                'source "generate_additional_recorders.tcl";\n'
            ]
        lines += [
            f'set gm_file_loc ".";\n',
            'source "$model_files_path/analysis_time_hist.tcl";\n',
            'remove recorders;\n',
            'if {$NLTHA_OK == 1} {\n',
            '    source "evaluate_edps.tcl";\n',
            '    set analysis_status_file [open "NLTHA_STATUS.txt" "w"];\n',
            '    puts $analysis_status_file "SUCCESS"\n',
            '    close $analysis_status_file;\n'
            '}\n',
        ]
        if delete_rec:
            lines += [
                f'file delete -force -- "{self.nltha_rec_dir_name}";\n',
            ]
        lines += [
            'wipe;\n'
        ]
        with open(Utility.get_path(target_work_dir_path, 'run_nltha.tcl'), 'w') as fid:
            fid.writelines(lines)

    def get_prelim_analysis_files(self):
        return [item + self.file_ext for item in ['generate_damping', 'model_param_vals', 'model_attributes']]

    def get_nltha_files(self):
        return ['evaluate_edps.tcl', 'generate_edp_recorders.tcl']

    def run_prelim_analysis(self, target_work_dir_path, comp_env):
        self.run_file(target_work_dir_path, 'run_prelim_analysis.tcl', comp_env)

    def run_nltha(self, target_work_dir_path, comp_env):
        self.run_file(target_work_dir_path, 'run_nltha.tcl', comp_env)

    def run_file(self, target_work_dir_path, run_file_name, comp_env):
        curr_dir = Utility.get_path(os.getcwd())
        os.chdir(target_work_dir_path)
        if comp_env == 'local' and platform.system() == 'Windows':
            subprocess.call(
                # f"\"{self.local_opensees_path}\" \"{run_file_name}\"", shell=True,
                [self.local_opensees_path, run_file_name], shell=True,
                # stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
            )
        else:
            os.system(f"\"{self.local_opensees_path}\" \"{run_file_name}\"")
        if os.path.isfile('gmon.out'):
            os.remove('gmon.out')
        os.chdir(curr_dir)

    def perform_mat_test(self, mat_data_dir_path, mat_tag, input_data, num_incr):
        material_data_file_path = Utility.get_path(mat_data_dir_path, 'material_data.txt')
        with open(material_data_file_path, 'r') as fid:
            lines = fid.readlines()
            lines = [line for line in lines if not line.strip().isspace() and line.strip()]
            material_data = ' '.join([line for line in lines if line.split()[2] == str(mat_tag)][0].split()) + '\n'

        lines = list()
        lines.append(f'wipe;\n')
        lines.append(f'model testUniaxial;\n')
        lines.append(f'set mat_tag {mat_tag};\n')
        lines.append(f'set strain_history {{{input_data}}};\n')
        lines.append(f'set file_out "hysteresis_mat_$mat_tag.txt";\n')
        lines.append(f'set out [open $file_out w];\n')
        lines.append(material_data)

        lines += ['uniaxialTest $mat_tag;\n',
                  'set strain 0.0;\n',
                  'set i_time 0;\n',
                  'set strain [expr $strain];\n',
                  'strainUniaxialTest $strain;\n',
                  'set stress [stressUniaxialTest];\n',
                  'set tangent [tangUniaxialTest];\n',
                  'set i_time [expr $i_time+1];\n',
                  'puts $out "$strain $stress";\n',
                  f'set num_incr {int(num_incr)};\n',
                  'foreach {strain_extreme_val} $strain_history {\n',
                  '    set strain_incr [expr ($strain_extreme_val - $strain)/$num_incr];\n',
                  '    for {set i 0} {$i < $num_incr} {incr i 1} {\n',
                  '	       set strain [expr $strain+$strain_incr];\n',
                  '        strainUniaxialTest $strain;\n',
                  '        set stress [stressUniaxialTest];\n',
                  '        set tangent [tangUniaxialTest];\n',
                  '        set i_time [expr $i_time+1];\n',
                  '        puts $out "$strain $stress";\n',
                  '    }\n',
                  '}\n',
                  'close $out;\n',
                  'puts "MATERIAL TESTER RAN SUCCESSFULLY!";\n',
                  'wipe;\n'
                  ]
        with open(Utility.get_path(mat_data_dir_path, 'mat_test.tcl'), 'w+') as fid:
            fid.writelines(lines)
        self.run_file(mat_data_dir_path, 'mat_test.tcl', 'local')
        out = Utility.read_numpy_array_from_txt_file(
            Utility.get_path(mat_data_dir_path, f'hysteresis_mat_{mat_tag}.txt'))
        os.remove(Utility.get_path(mat_data_dir_path, f'hysteresis_mat_{mat_tag}.txt'))
        os.remove(Utility.get_path(mat_data_dir_path, 'mat_test.tcl'))
        return out

    @staticmethod
    def write_model_param_vals(work_dir_path, name_value_pairs):
        OpenSeesTcl.write_tcl_input_params_file(work_dir_path, name_value_pairs, "model_param_vals.tcl")

    @staticmethod
    def write_model_attributes(work_dir_path, name_value_pairs):
        OpenSeesTcl.write_tcl_input_params_file(work_dir_path, name_value_pairs, "model_attributes.tcl")

    @staticmethod
    def write_damping_file(damping_model, work_dir_path, **kwargs):
        lines = []
        if damping_model == 'rayleigh_damping':
            alpha = kwargs['alpha']
            beta = kwargs['beta']
            lines.extend([f'rayleigh {alpha} 0. 0. {beta}\n'])
        if damping_model == 'modal_damping':
            n_modes = kwargs['n_modes']
            damp_ratios = kwargs['damp_ratios']
            lines.extend(
                [
                    f'eigen {int(n_modes)}\n',
                    f"modalDamping {' '.join(np.array(damp_ratios).astype(str))}\n"
                    f"wipeAnalysis\n"
                ]
            )
        with open(Utility.get_path(work_dir_path, 'generate_damping.tcl'), 'w') as fid:
            fid.writelines(lines)
        return

    @staticmethod
    def get_stress_strain_recorder_command(recorder_file_path, ele_tag, sec_num, y, z, mat_tag):
        new_recorder_line = f'recorder Element -file "{recorder_file_path}" -time -ele {ele_tag} ' \
                            f'section {sec_num} fiber {y} {z} {mat_tag} stressStrain\n'
        return new_recorder_line

    @staticmethod
    def get_elem_deformation_recorder_command(recorder_file_path, ele_tag):
        new_recorder_line = f'recorder Element -file "{recorder_file_path}" -time -ele {ele_tag} deformation;\n'
        return new_recorder_line

    @staticmethod
    def write_tcl_input_params_file(work_dir_path, name_value_pairs, filename):
        with open(Utility.get_path(work_dir_path, filename), 'w') as fid:
            for name, value in name_value_pairs:
                Utility.write_tcl_input_file(
                    fid, value, name
                )

    @staticmethod
    def get_recorder_files_to_read(rec_gen_file_path):
        with open(rec_gen_file_path, 'r') as fid:
            lines = fid.readlines()
        lines = [line.split() for line in lines if 'recorder' in line and '-file' in line]
        return [line[line.index('-file') + 1] for line in lines]


class OpenSeesPy(StructuralAnalysisPlatform):
    def __init__(self, model_files_path, local_python_path):
        super().__init__(model_files_path)
        self.local_python_path = Utility.get_path(local_python_path)
        self.file_ext = '.py'

    def bat_file(self, dest_path):
        pass

    def write_additional_recorders(self, recorder_info, dest_path, file_open_mode):
        element_response_list = recorder_info.get('element_response_list', [])
        element_section_response_list = recorder_info.get('element_section_response_list', [])
        element_section_fiber_response_list = recorder_info.get('element_section_fiber_response_list', [])
        node_response_list = recorder_info.get('node_response_list', [])
        addnl_file_list = recorder_info.get('addnl_file_list', [])

        if 'a' in file_open_mode:
            with open(Utility.get_path(dest_path, 'generate_additional_recorders.py'), 'r') as fid:
                lines = fid.readlines()
                if not lines[-1].endswith('\n'):
                    lines[-1] = lines[-1] + '\n'
        else:
            lines = [
                'import openseespy.opensees as ops\n',
                '\n',
                '\n',
                'def generate_additional_recorders():\n',
            ]

        importlib_line = 'import importlib.util\n'
        if len(addnl_file_list) > 0 and importlib_line not in lines:
            lines = [importlib_line, *lines]

        nltha_dir = self.nltha_rec_dir_name
        for item in element_response_list:
            ele_tag, what = item
            file_name = Utility.get_path(nltha_dir, f'element_response_elem_{ele_tag}_resp_{what}.txt')
            recorder_line = f'    ops.recorder(\'Element\', \'-file\', \'{file_name}\', \'-time\', \'-ele\', {ele_tag}, \'{what}\')\n'
            lines.append(recorder_line)

        for item in element_section_response_list:
            ele_tag, sec_num, what = item
            file_name = Utility.get_path(nltha_dir, f'element_section_response_elem_{ele_tag}_secnum_{sec_num}_resp_{what}.txt')
            recorder_line = f'    ops.recorder(\'Element\', \'-file\', \'{file_name}\', \'-time\', \'-ele\', {ele_tag}, \'section\', {sec_num}, \'{what}\')\n'
            lines.append(recorder_line)

        for item in element_section_fiber_response_list:
            ele_tag, sec_num, fib_num, y, z, mat_tag = item
            file_name = Utility.get_path(nltha_dir, f'element_section_fiber_response_elem_{ele_tag}_secnum_{sec_num}_fib_{fib_num}_mat_{mat_tag}_resp_stressStrain.txt')
            recorder_line = f'    ops.recorder(\'Element\', \'-file\', \'{file_name}\', \'-time\', \'-ele\', {ele_tag}, \'section\', {sec_num}, \'fiber\', {y}, {z}, {mat_tag}, \'stressStrain\')\n'
            lines.append(recorder_line)

        for item in node_response_list:
            node_tag, dof, what = item
            file_name = Utility.get_path(nltha_dir, f'node_response_node_{node_tag}_dof_{dof}_resp_{what}.txt')
            recorder_line = f'    ops.recorder(\'Node\', \'-file\', \'{file_name}\', \'-time\', \'-node\', {node_tag}, \'-dof\', {dof}, \'{what}\')\n'
            lines.append(recorder_line)

        for item in addnl_file_list:
            module_name = item.split('.')[0]
            lines.extend(OpenSeesPy.import_module(dest_path, module_name, module_name, tab='    '))
            lines.append(f'    {module_name}.{module_name}()\n')

        with open(Utility.get_path(dest_path, 'generate_additional_recorders.py'), 'w') as fid:
            fid.writelines(lines)

    def open_edp_recorder_file(self, work_dir_path, rec_gen_file_open_mode, rec_save_dir_path):
        rec_gen_file_fid = open(Utility.get_path(work_dir_path, 'generate_edp_recorders' + self.file_ext), rec_gen_file_open_mode)
        if 'w' in rec_gen_file_open_mode:
            lines = [
                f'import os\n',
                f'import openseespy.opensees as ops\n',
                f'\n\n',
                f'def generate_edp_recorders():\n',
                f'    os.makedirs(\'{rec_save_dir_path}\', exist_ok=True)\n',
            ]
        else:
            lines = [
                f'    os.makedirs(\'{rec_save_dir_path}\', exist_ok=True)\n',
            ]
        rec_gen_file_fid.writelines(lines)
        return rec_gen_file_fid

    def write_evaluate_edp_and_helpers(self, edp_tag, max_what, resp, resp_num):
        get_resp_vals = [
            f"def get_{resp}_vals(file_path):\n",
            f"    data = np.loadtxt(file_path, delimiter=' ')\n",
            f"    return data[:, {resp_num}]\n",
        ]

        edp_eval_func_name = f'evaluate_edp_{edp_tag}'

        max_list = [
            f"def {edp_eval_func_name}(x):\n",
            f"    return np.max(x)\n",
        ]

        abs_min_list = [
            f"def {edp_eval_func_name}(x):\n",
            f"    return np.abs(np.min(x))\n",
        ]

        max_range = [
            f"def {edp_eval_func_name}(x):\n",
            f"    max_x = np.max(x)\n",
            f"    max_x_i = np.argmax(x)\n",
            f"    x = x[max_x_i:]\n",
            f"    min_x_following = np.min(x)\n",
            f"    return max_x - min_x_following\n",
        ]

        model_files_path = self.model_files_path

        lines = ['import numpy as np\n\n\n'] + get_resp_vals
        get_resp_vals_file_path = Utility.get_path(model_files_path, f'get_{resp}_vals.py')
        with open(get_resp_vals_file_path, 'w') as fid:
            fid.writelines(lines)

        lines = ['import numpy as np\n\n\n']
        if max_what == 'tension':
            lines += max_list
        if max_what == 'compression':
            lines += abs_min_list
        if max_what == 'range':
            lines += max_range
        eval_file_path = Utility.get_path(model_files_path, f'{edp_eval_func_name}.py')
        with open(eval_file_path, 'w') as fid:
            fid.writelines(lines)

    def write_evaluate_edps(self, work_dir_path, edp, for_which, gen_rec, mode):
        rec_gen_file_path = Utility.get_path(work_dir_path, 'generate_edp_recorders' + self.file_ext)
        eval_edp_file_path = Utility.get_path(work_dir_path, 'evaluate_edps' + self.file_ext)

        recorder_files_to_read = OpenSeesPy.get_recorder_files_to_read(rec_gen_file_path)
        if edp.recorder_file_storage == 'shared':
            recorder_files_to_read = [file for file in recorder_files_to_read
                                      if f'/EDP/' in file]
        else:
            recorder_files_to_read = [file for file in recorder_files_to_read
                                      if f'/EDP_{edp.tag}/' in file]
        lines = list()
        if 'w' in mode:
            lines.extend(['import importlib.util\n', 'import os\n\n\n', 'def evaluate_edps():\n'])
            lines.extend(OpenSeesPy.import_module(self.model_files_path, 'ops_model', 'model', tab='    '))
        lines.extend(OpenSeesPy.import_module(self.model_files_path, f'evaluate_edp_{edp.tag}', f'edp_{edp.tag}', tab='    '))
        if gen_rec:
            lines.extend(OpenSeesPy.import_module(self.model_files_path, f'get_{edp.resp}_vals', f'resp_{edp.resp}', tab='    '))
        lines.append(f'    os.makedirs(\'EDP_{edp.tag}_Results\', exist_ok=True)\n')
        edp_strings = edp.get_edp_strings(for_which)
        if isinstance(edp.normalize_with, str):
            normalize_with = 'model.' + edp.normalize_with.strip('$')
        else:
            normalize_with = f'{edp.normalize_with}'
        for i, edp_string in enumerate(edp_strings):
            edp_var_name = f'edp_{edp.tag}_' + edp_string
            files_to_eval = [file for file in recorder_files_to_read if
                             edp_string in file]
            var_names_to_eval = [os.path.basename(file).split('.')[0] for file in files_to_eval]
            if gen_rec:
                for (var_name, file) in zip(var_names_to_eval, files_to_eval):
                    lines.append(f'    {var_name} = resp_{edp.resp}.get_{edp.resp}_vals("{file}")\n')
            lines.append(f'    {edp_var_name} = []\n')
            for var_name in var_names_to_eval:
                lines.append(f'    {edp_var_name}.append(edp_{edp.tag}.evaluate_edp_{edp.tag}({var_name}))\n')
            lines.append(f'    {edp_var_name} = max({edp_var_name}) / {normalize_with}\n')
            lines.append(f'    with open(\'EDP_{edp.tag}_Results/{edp_var_name}.txt\', \'w\') as edp_file:\n')
            lines.append(f'        edp_file.write(str({edp_var_name}) + \'\\n\')\n')

        with open(eval_edp_file_path, mode) as eval_file_fid:
            eval_file_fid.writelines(lines)
        return

    def write_run_prelim_analysis(self, work_dir_path, num_modes):
        model_files_path = self.model_files_path
        model_info_directory = Utility.get_path(work_dir_path, self.model_info_dir_name)
        module_name = 'ops_model'
        lines = [
            f'import os\n',
            f'import json\n',
            f'import importlib.util\n\n',
        ]
        lines.extend(OpenSeesPy.import_module(self.model_files_path, module_name, 'model'))
        lines.extend(
            [
                f'model_files_path = \'{model_files_path}\'\n',
                f'model_info_directory = \'{model_info_directory}\'\n',
                'os.makedirs(model_info_directory, exist_ok=True)\n',
            ]
        )
        lines.extend(OpenSeesPy.read_json('model_param_vals.json', 'model_param_vals'))
        lines.extend(OpenSeesPy.read_json('model_attributes.json', 'model_attributes'))
        names = [
            'model_files_path',
            'model_info_directory',
            'write_model_files',
            'modes',
            'model_param_vals',
            'model_attributes',
        ]
        values = [
            'model_files_path',
            'model_info_directory',
            1,
            num_modes,
            'model_param_vals',
            'model_attributes',
        ]
        lines.extend(OpenSeesPy.create_dict(names, values, 'model_input_dict'))
        lines.extend(
            [
                "model.main(model_input_dict)\n",
                'model.analysis_gravity(model_input_dict)\n',
                'model.analysis_prelim(model_input_dict)\n',
            ]
        )
        with open(Utility.get_path(work_dir_path, 'run_prelim_analysis.py'), 'w') as fid:
            fid.writelines(lines)

    def write_run_nltha(self, target_work_dir_path, work_dir_path_prelim_analysis, **kwargs):
        write_model_files = kwargs.get('write_model_files', 0)
        generate_additional_recorders = kwargs.get('generate_additional_recorders', False)
        delete_rec = kwargs.get('delete_rec', True)

        lines = [
            'import json\n',
            'import importlib.util\n',
            'import shutil\n',
            'import openseespy.opensees as ops\n\n',
            'with open("NLTHA_STATUS.txt", "w") as fid:\n',
            '    fid.write("FAIL\\n")\n',
        ]
        lines.extend(OpenSeesPy.import_module(self.model_files_path, 'ops_model', 'model'))
        lines.extend(OpenSeesPy.import_module(work_dir_path_prelim_analysis, 'generate_damping', 'damp'))
        lines.extend(OpenSeesPy.import_module('', 'generate_edp_recorders', 'edp_rec'))
        lines.extend(OpenSeesPy.import_module('', 'evaluate_edps', 'edp_eval'))
        lines.extend(OpenSeesPy.read_json(Utility.get_path(work_dir_path_prelim_analysis, 'model_param_vals.json'), 'model_param_vals'))
        lines.extend(OpenSeesPy.read_json(Utility.get_path(work_dir_path_prelim_analysis, 'model_attributes.json'), 'model_attributes'))
        names = [
            'write_model_files',
            'model_param_vals',
            'model_attributes',
            'gm_file_loc',
        ]
        values = [
            write_model_files,
            'model_param_vals',
            'model_attributes',
            '\'.\'',
        ]
        if write_model_files == 1:
            model_info_directory = Utility.get_path(target_work_dir_path, self.model_info_dir_name)
            lines = ['import os\n', *lines]
            lines += [
                f'model_info_directory = \'{model_info_directory}\'\n',
                'os.makedirs(model_info_directory, exist_ok=True)\n',
            ]
            names += ['model_info_directory']
            values += ['model_info_directory']
        lines.extend(OpenSeesPy.create_dict(names, values, 'model_input_dict'))
        lines.extend(
            [
                "model.main(model_input_dict)\n",
                'model.analysis_gravity(model_input_dict)\n',
                'damp.generate_damping()\n',
                'edp_rec.generate_edp_recorders()\n',
            ]
        )
        if generate_additional_recorders:
            lines.extend(OpenSeesPy.import_module(self.model_files_path, 'generate_additional_recorders', 'addnl_rec'))
            lines += [
                'addnl_rec.generate_additional_recorders()\n',
            ]

        lines += [
            'nltha_ok = model.analysis_time_hist(model_input_dict)\n',
            'ops.remove(\'recorders\')\n',
            'if nltha_ok:\n',
            '    edp_eval.evaluate_edps()\n',
            '    with open("NLTHA_STATUS.txt", "w") as fid:\n',
            '        fid.write("SUCCESS\\n")\n',
        ]
        if delete_rec:
            lines += [
                f'shutil.rmtree("{self.nltha_rec_dir_name}")\n',
            ]
        with open(Utility.get_path(target_work_dir_path, 'run_nltha.py'), 'w') as fid:
            fid.writelines(lines)

    def get_prelim_analysis_files(self):
        return ['generate_damping.py', 'model_param_vals.json', 'model_attributes.json']

    def get_nltha_files(self):
        return ['evaluate_edps.py', 'generate_edp_recorders.py']

    def run_prelim_analysis(self, target_work_dir_path, comp_env):
        self.run_file(target_work_dir_path, 'run_prelim_analysis.py', comp_env)

    def run_nltha(self, target_work_dir_path, comp_env):
        self.run_file(target_work_dir_path, 'run_nltha.py', comp_env)

    def run_file(self, target_work_dir_path, run_file_name, comp_env):
        curr_dir = Utility.get_path(os.getcwd())
        os.chdir(target_work_dir_path)
        if comp_env == 'local' and platform.system() == 'Windows':
            subprocess.call(
                [self.local_python_path, run_file_name], shell=True
            )
        os.chdir(curr_dir)

    def perform_mat_test(self, mat_data_dir_path, mat_tag, input_data, num_incr):
        raise NotImplementedError("openseespy implementation pending")

    @staticmethod
    def write_model_param_vals(work_dir_path, name_value_pairs):
        OpenSeesPy.write_json_input_params_file(work_dir_path, name_value_pairs, "model_param_vals.json")

    @staticmethod
    def write_model_attributes(work_dir_path, name_value_pairs):
        OpenSeesPy.write_json_input_params_file(work_dir_path, name_value_pairs, "model_attributes.json")

    @staticmethod
    def write_damping_file(damping_model, work_dir_path, **kwargs):
        lines = [
                f'import openseespy.opensees as ops\n',
                f'\n',
                f'\n',
                f'def generate_damping():\n',
        ]
        if damping_model == 'rayleigh_damping':
            alpha = kwargs['alpha']
            beta = kwargs['beta']
            lines.extend([
                f'    ops.rayleigh({alpha}, 0., 0., {beta})\n',
            ])
        if damping_model == 'modal_damping':
            n_modes = kwargs['n_modes']
            damp_ratios = kwargs['damp_ratios']
            lines.extend([
                f'    ops.eigen({int(n_modes)})\n',
                f"    ops.modalDamping({', '.join(np.array(damp_ratios).astype(str))})\n",
                f'    ops.wipeAnalysis()\n',
            ])
        if damping_model is None:
            lines.append('    pass\n')
        with open(Utility.get_path(work_dir_path, 'generate_damping.py'), 'w') as fid:
            fid.writelines(lines)
        return

    @staticmethod
    def get_stress_strain_recorder_command(recorder_file_path, ele_tag, sec_num, y, z, mat_tag):
        new_recorder_line = f'    ops.recorder(\'Element\', \'-file\', \'{recorder_file_path}\', \'-time\', \'-ele\', {ele_tag}, ' \
                            f'\'section\', {sec_num}, \'fiber\', {y}, {z}, {mat_tag}, \'stressStrain\')\n'
        return new_recorder_line

    @staticmethod
    def get_elem_deformation_recorder_command(recorder_file_path, ele_tag):
        new_recorder_line = f'    ops.recorder(\'Element\', \'-file\', \'{recorder_file_path}\', \'-time\', \'-ele\', {ele_tag}, \'deformation\')\n'
        return new_recorder_line

    @staticmethod
    def get_recorder_files_to_read(rec_gen_file_path):
        with open(rec_gen_file_path, 'r') as fid:
            lines = fid.readlines()
        lines = [[item.strip(' ,') for item in line.split(',')] for line in lines if 'ops.recorder(' in line and '-file' in line]
        return [line[line.index('\'-file\'') + 1].strip('\'') for line in lines]

    @staticmethod
    def import_module(module_path, module_name, alias, tab=''):
        # Assumes importlib.util is imported
        module_path = Utility.get_path(module_path, module_name + '.py')
        lines = [
            f"{tab}spec = importlib.util.spec_from_file_location('{module_name}', '{module_path}')\n",
            f"{tab}{alias} = importlib.util.module_from_spec(spec)\n",
            f"{tab}spec.loader.exec_module({alias})\n",
        ]
        return lines

    @staticmethod
    def read_json(file_path, var_name):
        # Assumes json is imported
        lines = [
            f'with open(\'{file_path}\', \'r\') as file:\n',
            f'    {var_name} = json.load(file)\n',
        ]
        return lines

    @staticmethod
    def create_dict(names, values, dict_name):
        lines = [f'{dict_name} = {{\n']
        for name, value in zip(names, values):
            lines.append(f'    \'{name}\': {value},\n')
        lines.append('}\n')
        return lines

    @staticmethod
    def write_json_input_params_file(work_dir_path, name_value_pairs, filename):
        to_dump_names = []
        to_dump_values = []
        for name, value in name_value_pairs:
            to_dump_names.extend(name)
            to_dump_values.extend(value)
        to_dump = dict(zip(to_dump_names, to_dump_values))
        with open(Utility.get_path(work_dir_path, filename), 'w') as fid:
            json.dump(to_dump, fid, indent=4)

    @staticmethod
    def unique_append(new_commands, lines, ex_commands):
        for command in new_commands:
            if command not in ex_commands:
                lines.append(command)
        return lines
