import numpy as np


def get_disp_vals(file_path):
    return np.loadtxt(file_path)[:, 1]
