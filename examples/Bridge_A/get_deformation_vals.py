import numpy as np


def get_deformation_vals(file_path):
    data = np.loadtxt(file_path, delimiter=' ')
    return data[:, 1]
