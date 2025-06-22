# -*- coding: utf-8 -*-
"""
Created on Sat Dec  7 18:32:34 2019

@author: Angshuman Deb
"""

import numpy as np
from pyDOE2 import lhs
import sys
from scipy import integrate
from scipy import interpolate


class Mixture:

    def __init__(self, name):
        self.name = name

    def __call__(self, prob_dists, **kwargs):
        return FrozenMixture(self, prob_dists, **kwargs)

    def mean(self, prob_dists, **kwargs):
        return FrozenMixture(self, prob_dists, **kwargs).mean()

    def std(self, prob_dists, **kwargs):
        return FrozenMixture(self, prob_dists, **kwargs).std()

    def support(self, prob_dists, **kwargs):
        return FrozenMixture(self, prob_dists, **kwargs).support()

    def pdf(self, x, prob_dists, **kwargs):
        return FrozenMixture(self, prob_dists, **kwargs).pdf(x)

    def cdf(self, x, prob_dists, **kwargs):
        return FrozenMixture(self, prob_dists, **kwargs).cdf(x)

    def ppf(self, x, prob_dists, **kwargs):
        return FrozenMixture(self, prob_dists, **kwargs).ppf(x)

    def entropy(self, prob_dists, **kwargs):
        return FrozenMixture(self, prob_dists, **kwargs).entropy()

    def rvs(self, prob_dists, **kwargs):
        size = kwargs.get('size', None)
        random_state = kwargs.get('random_state', None)
        method = kwargs.get('method', 'mcs')
        each = kwargs.get('each', 'arbitrary')
        weights = kwargs.get('weights', np.array([1.0 / len(prob_dists)] * len(prob_dists)))
        return FrozenMixture(self, prob_dists, each=each,
                             weights=weights).rvs(size=size, random_state=random_state, method=method)


mixture = Mixture('mixture')


# Frozen Mixture
class FrozenMixture:

    ####################################################################################################################
    # Constructor
    ####################################################################################################################

    def __init__(self, dist, prob_dists, **kwargs):
        # prob_dists : either a 2d array or a list

        self.dist = dist
        self.each = kwargs.get('each', 'arbitrary')
        self.mixture_size = len(prob_dists)
        self.weights = kwargs.get('weights', np.array([1.0 / self.mixture_size] * self.mixture_size))

        if self.each == 'arbitrary':
            # prob_dists is a list
            self.prob_dists = prob_dists
            self.mean_each = np.array([prob_dist.mean() for prob_dist in self.prob_dists])
            self.std_each = np.array([prob_dist.std() for prob_dist in self.prob_dists])
            self.support_each = np.array([prob_dist.support() for prob_dist in self.prob_dists])
            self.args = [FrozenMixture.get_all_dist_params(prob_dist) for prob_dist in self.prob_dists]
            self.kwds = dict()
        else:
            # prob_dists is a 2d array
            self.prob_dists = self.each(*prob_dists.T)
            self.mean_each = self.prob_dists.mean()
            self.std_each = self.prob_dists.std()
            self.support_each = np.array(self.prob_dists.support()).T
            self.args = [list(item) for item in prob_dists]
            self.kwds = dict()

        self._mean = np.sum(self.mean_each * self.weights)
        self._std = np.sum(
            self.weights * (np.square(self.std_each) + np.square(self.mean_each - self._mean))) ** 0.5
        self._support = (np.min(self.support_each[:, 0]), np.max(self.support_each[:, 1]))
        self._k = 8
        self._approx_support = FrozenMixture.get_approximate_support(self._support, self._mean, self._std, self._k)
        self._ppf_data = None

    ####################################################################################################################
    # Public functionalities
    ####################################################################################################################

    def mean(self):
        return self._mean

    # ------------------------------------------------------------------------------------------------------------------

    def std(self):
        return self._std

    # ------------------------------------------------------------------------------------------------------------------

    def median(self):
        return self.ppf(0.5)

    # ------------------------------------------------------------------------------------------------------------------

    def support(self):
        return self._support

    # ------------------------------------------------------------------------------------------------------------------

    def pdf(self, x):
        scalar = True if np.isscalar(x) else False
        x = np.array(x)
        x_temp = x.ravel()[:, np.newaxis]
        temp = self.weights
        for itr in range(len(x_temp.shape)):
            temp = temp[:, np.newaxis]

        if self.each == 'arbitrary':
            pdf_x = list()
            for itr in range(self.mixture_size):
                pdf_x.append(self.prob_dists[itr].pdf(x_temp))
            pdf_x = np.array(pdf_x)
        else:
            pdf_x_temp = self.prob_dists.pdf(x_temp)
            pdf_x = np.einsum('ijd->jid', pdf_x_temp[:, :, np.newaxis])

        pdf_x = pdf_x * temp
        pdf_x = np.sum(pdf_x, axis=0).reshape(x.shape)
        return pdf_x.item() if scalar else pdf_x

    # ------------------------------------------------------------------------------------------------------------------

    def cdf(self, x):
        scalar = True if np.isscalar(x) else False
        x = np.array(x)
        x_temp = x.ravel()[:, np.newaxis]
        temp = self.weights
        for itr in range(len(x_temp.shape)):
            temp = temp[:, np.newaxis]

        if self.each == 'arbitrary':
            cdf_x = list()
            for itr in range(self.mixture_size):
                cdf_x.append(self.prob_dists[itr].cdf(x_temp))
            cdf_x = np.array(cdf_x)
        else:
            cdf_x_temp = self.prob_dists.cdf(x_temp)
            cdf_x = np.einsum('ijd->jid', cdf_x_temp[:, :, np.newaxis])

        cdf_x = cdf_x * temp
        cdf_x = np.sum(cdf_x, axis=0).reshape(x.shape)
        return cdf_x.item() if scalar else cdf_x

    # ------------------------------------------------------------------------------------------------------------------

    def entropy(self):
        eps = sys.float_info.epsilon
        entr = integrate.quad(FrozenMixture.entropy_integrand,
                              self._approx_support[0] + eps,
                              self._approx_support[-1] - eps,
                              args=(self,))
        return entr[0]

    # ------------------------------------------------------------------------------------------------------------------

    # def rvs(self, size=None, random_state=None):
    #     np.random.seed(random_state)
    #     mixture_inds = np.random.choice(self.mixture_size, size=size, p=self.weights)
    #     if np.isscalar(mixture_inds):
    #         mixture_inds = np.array([mixture_inds])
    #     mixture_inds = mixture_inds.flatten()
    #     if self.each == 'arbitrary':
    #         rvs_temp = np.array([self.prob_dists[ind].rvs() for ind in mixture_inds])
    #     else:
    #         rvs_temp = np.array([self.prob_dists.rvs()[ind] for ind in mixture_inds])
    #     return rvs_temp.reshape(size) if size is not None else rvs_temp.item()

    # ------------------------------------------------------------------------------------------------------------------

    def rvs(self, size=None, random_state=None, method='mcs'):
        if size is None:
            size = 1
            none = True
        else:
            none = False

        x = np.zeros(size)
        if method == 'mcs':
            np.random.seed(random_state)
            x = np.random.uniform(0.0, 1.0, size)
        elif method == 'lhs':
            if not np.isscalar(size):
                n_samples = size[0]
                if len(size) == 2:
                    n_set = size[1]
                elif len(size) > 2:
                    raise ValueError("Non-scalar size for this type of sampling method "
                                     "should not be of length more than 2 and should represent "
                                     "(number of realizations of random variable, number of sets of realizations)!")
                else:
                    n_set = 1
            else:
                n_samples = size
                n_set = 1

            for itr_set in range(n_set):
                if random_state is None:
                    rng_seed = None
                else:
                    if itr_set > 0:
                        rng_seed = int(''.join([str(i) for i in [random_state, itr_set]]))
                    else:
                        rng_seed = random_state
                temp = lhs(1, n_samples, random_state=rng_seed)
                if len(x.shape) == 2:
                    x[:, itr_set] = temp
                else:
                    x = temp

        rvs_temp = self.ppf(x)
        return rvs_temp.reshape(size) if not none else rvs_temp.item()

    # ------------------------------------------------------------------------------------------------------------------

    def ppf(self, p):
        if self._ppf_data is None:
            self.set_ppf_data()

        ppf_interp = interpolate.interp1d(self._ppf_data[:, 1], self._ppf_data[:, 0], fill_value="extrapolate")
        x = ppf_interp(p)

        if not np.isscalar(p):
            scalar = False
            x[p < 0.0] = np.nan
            x[p > 1.0] = np.nan
            x[p == 0.0] = self._support[0]
            x[p == 1.0] = self._support[1]
        else:
            scalar = True
            if p < 0.0:
                x = np.nan
            elif p > 1.0:
                x = np.nan
            elif p == 0.0:
                x = self._support[0]
            elif p == 1.0:
                x = self._support[1]
        return x.item() if scalar else x

    # ------------------------------------------------------------------------------------------------------------------

    def set_ppf_data(self):
        # p = np.linspace(self.cdf(self._approx_support[0]), self.cdf(self._approx_support[-1]), 500)
        # x = self.ppf_helper_bisection(p,
        #                               np.full(len(p), self._approx_support[0]),
        #                               np.full(len(p), self._approx_support[-1]))

        x = np.linspace(*self._approx_support, 1000)
        p = self.cdf(x)
        self._ppf_data = np.column_stack((x, p))

    # ------------------------------------------------------------------------------------------------------------------

    def ppf_helper_bisection(self, p, a, b, max_iterations=100, tolerance=1e-4):

        if p.size == 0:
            return np.empty(0)

        f_a = self.cdf(a) - p

        for i in range(max_iterations):

            x = a + (b - a) / 2
            f_x = self.cdf(x) - p

            if np.sum(np.abs(f_x) < tolerance) == len(f_x):
                return x

            ind1 = np.sign(f_a) * np.sign(f_x) > 0
            ind2 = ~ind1

            a[ind1] = x[ind1]
            f_a[ind1] = f_x[ind1]
            b[ind2] = x[ind2]

        return a + (b - a) / 2

    ####################################################################################################################
    # Static methods
    ####################################################################################################################

    @staticmethod
    def entropy_integrand(x, m):
        px = m.pdf(x)
        return -px * np.log(px)

    # ------------------------------------------------------------------------------------------------------------------

    @staticmethod
    def get_approximate_support(supp, m, s, k):
        return max(m - k * s, supp[0]), min(m + k * s, supp[-1])

    # ------------------------------------------------------------------------------------------------------------------

    @staticmethod
    def get_all_dist_params(prob_dist):
        temp = prob_dist.dist.__dict__['_parse_args'](*prob_dist.args, **prob_dist.kwds)
        return [*temp[0], *temp[1:len(temp)]]

    # ------------------------------------------------------------------------------------------------------------------
