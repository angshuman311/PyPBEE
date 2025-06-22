import numpy as np


def evaluate_edp_3(x):
    max_x = np.max(x)
    max_x_i = np.argmax(x)
    x = x[max_x_i:]
    min_x_following = np.min(x)
    return max_x - min_x_following
