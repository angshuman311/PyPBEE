# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 21:27:13 2019

@author: Angshuman Deb
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from numpy.linalg import lstsq


class InterpExterpModel:

    ####################################################################################################################
    # Static methods
    ####################################################################################################################

    @staticmethod
    def piecewise_power_interp_continuous_extrap(x, y, xq):
        return np.exp(interp1d(np.log(x), np.log(y), kind='linear', fill_value='extrapolate')(np.log(xq)))

    # ------------------------------------------------------------------------------------------------------------------

    @staticmethod
    def piecewise_linear_interp_constant_extrap(x, y, xq):
        return interp1d(np.log(x), y, kind='linear', fill_value=(y[0], y[-1]), bounds_error=False)(np.log(xq))
        # return interp1d(x, y, kind='linear', fill_value=(y[0], y[-1]), bounds_error=False)(xq)

    # ------------------------------------------------------------------------------------------------------------------

    @staticmethod
    def power_law_regression_lin_params(x, y):
        matrix_a = np.vstack([np.log(x), np.ones(len(x))]).T
        b, log_a = lstsq(matrix_a, np.log(y), rcond=None)[0]
        return b, np.exp(log_a)

    # ------------------------------------------------------------------------------------------------------------------

    @staticmethod
    def power_law_regression_lin(x, y, xq):
        b, a = InterpExterpModel.power_law_regression_lin_params(x, y)
        return a * (xq ** b)

    # ------------------------------------------------------------------------------------------------------------------

    @staticmethod
    def inc_dec_regression(x, y, xq):
        a = np.amax(y)
        b = x[np.argmax(y)]

        def func(_x, c):
            return a * (_x / b) * c / (c - 1 + (_x / b) ** c)

        if list(y).index(a) == len(y) - 1:
            # y_max is the last element of y 
            return InterpExterpModel.piecewise_linear_interp_linconst_extrap(x, y, xq)

        p0 = [2]
        bounds = (1, np.inf)
        popt = curve_fit(func, x, y, p0=p0, bounds=bounds)[0]
        return func(xq, *popt)

    # ------------------------------------------------------------------------------------------------------------------

    @staticmethod
    def piecewise_linear_interp_linconst_extrap(x, y, xq):
        xq_1 = xq[xq < x[0]]
        xq_2 = xq[xq >= x[0]]
        yq_1 = (y[0] / x[0]) * xq_1
        yq_2 = InterpExterpModel.piecewise_linear_interp_constant_extrap(x, y, xq_2)
        return np.hstack([yq_1, yq_2])

    # ------------------------------------------------------------------------------------------------------------------

    @staticmethod
    def power_law_regression_nonlin(x, y, xq):
        def func(_x, a, b):
            return a * (_x ** b)

        b_guess, a_guess = InterpExterpModel.power_law_regression_lin_params(x, y)

        p0 = [a_guess, b_guess]
        bounds = ([0, 1], np.inf)
        try:
            popt = curve_fit(func, x, y, p0=p0, bounds=bounds)[0]
        except ValueError:
            return InterpExterpModel.power_law_regression_lin(x, y, xq)
        return func(xq, *popt)

    # ------------------------------------------------------------------------------------------------------------------

    @staticmethod
    def constant(x, y, xq):
        matrix_a = np.vstack([np.ones(len(x))]).T
        c = lstsq(matrix_a, y, rcond=None)[0]
        return c * np.ones(len(xq))

    # ------------------------------------------------------------------------------------------------------------------

    @staticmethod
    def linear_origin(x, y, xq):
        matrix_a = np.vstack([x]).T
        m = lstsq(matrix_a, y, rcond=None)[0]
        return m * xq

    # ------------------------------------------------------------------------------------------------------------------
