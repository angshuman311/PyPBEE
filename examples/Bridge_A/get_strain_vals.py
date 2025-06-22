import numpy as np


def get_strain_vals(file_path):
    data = np.loadtxt(file_path, delimiter=' ')
    return data[:, 2]
