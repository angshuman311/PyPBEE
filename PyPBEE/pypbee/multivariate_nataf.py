# -*- coding: utf-8 -*-
"""
Created on Sat Dec  7 18:32:34 2019

@author: Angshuman Deb
"""

import numpy as np
from scipy.stats import norm
from scipy.stats import multivariate_normal
from scipy.optimize import root_scalar
from pyDOE2 import lhs
import sys


class MultivariateNataf:

    def __init__(self):
        pass

    def __call__(self, prob_dist_list, corr_matrix, corr_matrix_z=None):
        return FrozenMultivariateNataf(self, prob_dist_list, corr_matrix, corr_matrix_z)

    def pdf(self, x, prob_dist_list, corr_matrix, corr_matrix_z=None):
        return FrozenMultivariateNataf(self, prob_dist_list, corr_matrix, corr_matrix_z).pdf(x)

    def cdf(self, x, prob_dist_list, corr_matrix, corr_matrix_z=None):
        return FrozenMultivariateNataf(self, prob_dist_list, corr_matrix, corr_matrix_z).cdf(x)

    def rvs(self, prob_dist_list, corr_matrix, corr_matrix_z=None, size=None, random_state=None, method='mcs'):
        return FrozenMultivariateNataf(self, prob_dist_list, corr_matrix, corr_matrix_z).rvs(size=size,
                                                                                             random_state=random_state,
                                                                                             method=method)


multivariate_nataf = MultivariateNataf()


class FrozenMultivariateNataf:

    ####################################################################################################################
    # Constructor
    ####################################################################################################################

    def __init__(self, dist, prob_dist_list, corr_matrix, corr_matrix_z=None):
        self._dist = dist
        self.prob_dist_list = prob_dist_list
        self.ndim = len(self.prob_dist_list)
        self._k = 8
        self._mean = np.array([prob_dist.mean() for prob_dist in self.prob_dist_list])
        self._std = np.array([prob_dist.std() for prob_dist in self.prob_dist_list])

        if corr_matrix.shape[0] != corr_matrix.shape[1] or corr_matrix.shape[0] != self.ndim:
            raise ValueError("Invalid corr_matrix!")
        if len(corr_matrix.shape) > 2:
            # implies each face for components of mixture
            temp = [prob_dist.dist.name for prob_dist in self.prob_dist_list]
            if temp != ['mixture'] * self.ndim:
                raise ValueError("Expecting all prob dists for this type of corr_matrix input to be mixtures!")
            temp = [prob_dist.mixture_size for prob_dist in self.prob_dist_list]
            if temp[1:] != temp[:-1]:
                raise ValueError("Expecting all mixtures for this type of corr_matrix input to have same size!")
            if temp != [corr_matrix.shape[2]] * self.ndim:
                raise ValueError(f"Expecting all mixtures for this type of corr_matrix input to have size = "
                                 f"{corr_matrix.shape[2]}!")
            temp = np.array([prob_dist.weights for prob_dist in self.prob_dist_list])
            if (temp != temp[0]).any():
                raise ValueError("Expecting all mixtures for this type of corr_matrix input to have "
                                 "equal weights!")

            # corr_matrix : (n_rv, n_rv, n_s)
            # total_mean : (n_rv,)
            total_mean = np.array([prob_dist.mean() for prob_dist in self.prob_dist_list])
            # total_sig : (n_rv,)
            total_sig = np.array([prob_dist.std() for prob_dist in self.prob_dist_list])
            # mean_all : (n_s, n_rv)
            mean_all = np.array([prob_dist.mean_each for prob_dist in self.prob_dist_list]).T
            # sig_all : (n_s, n_rv)
            sig_all = np.array([prob_dist.std_each for prob_dist in self.prob_dist_list]).T
            # weights : (n_s,)
            weights = self.prob_dist_list[0].weights
            self.corr_matrix = FrozenMultivariateNataf.corr_clipped(
                FrozenMultivariateNataf.get_total_correlation(corr_matrix,
                                                              total_mean,
                                                              total_sig,
                                                              mean_all,
                                                              sig_all,
                                                              weights)
            )
        else:
            self.corr_matrix = FrozenMultivariateNataf.corr_clipped(corr_matrix)

        temp = np.diag(self._std)
        self._covariance = temp @ self.corr_matrix @ temp

        if corr_matrix_z is None:
            self.corr_matrix_z = FrozenMultivariateNataf.corr_clipped(
                FrozenMultivariateNataf.get_corr_matrix_z(self.prob_dist_list, self.corr_matrix, self._k))
        else:
            self.corr_matrix_z = FrozenMultivariateNataf.corr_clipped(corr_matrix_z)

        eig_val, eig_vec = np.linalg.eig(self.corr_matrix_z)
        self.eig_val = np.real(eig_val)
        self.eig_vec = np.real(eig_vec)

    ####################################################################################################################
    # Public functionalities
    ####################################################################################################################

    def mean(self):
        return self._mean

    # ------------------------------------------------------------------------------------------------------------------

    def std(self):
        return self._std

    # ------------------------------------------------------------------------------------------------------------------

    def covariance(self):
        return self._covariance

    # ------------------------------------------------------------------------------------------------------------------

    def pdf(self, x):
        return self._nataf_pdf(x)

    # ------------------------------------------------------------------------------------------------------------------

    def cdf(self, x):
        return self._nataf_cdf(x)

    # ------------------------------------------------------------------------------------------------------------------

    def rvs(self, size=None, random_state=None, method='mcs'):
        if size is None:
            size = 1
        try:
            size = np.array(size).item()
        except ValueError:
            pass
        rvs_temp = self.sampler(size, random_state, method)
        if size == 1 and self.ndim == 1:
            return rvs_temp.item()
        if size == 1 and self.ndim > 1:
            return rvs_temp.flatten()
        if np.isscalar(size) and self.ndim == 1:
            return rvs_temp.flatten()
        return rvs_temp.squeeze()

    # ------------------------------------------------------------------------------------------------------------------

    def _nataf_pdf(self, x):
        # x : (i, j, k, ..., ndim)
        z = np.ones(x.shape)
        term_1 = np.ones(x.shape[:-1])
        term_3 = np.ones(term_1.shape)
        for itr in range(self.ndim):
            term_1 *= self.prob_dist_list[itr].pdf(x[..., itr])
            temp = self.prob_dist_list[itr].cdf(x[..., itr])
            temp[temp > 1.0] = 1.0
            z_itr = norm.ppf(temp)
            z_itr[z_itr >= self._k] = self._k
            z_itr[z_itr <= -self._k] = -self._k
            term_3 *= norm.pdf(z_itr)
            z[..., itr] = z_itr
        mvn = multivariate_normal(mean=np.zeros(self.ndim),
                                  cov=self.corr_matrix_z)
        term_2 = mvn.pdf(z)
        return term_1 * term_2 / term_3

    # ------------------------------------------------------------------------------------------------------------------

    def _nataf_cdf(self, x):
        # x : (i, j, k, ..., ndim)
        z = np.ones(x.shape)
        for itr in range(self.ndim):
            temp = self.prob_dist_list[itr].cdf(x[..., itr])
            temp[temp > 1.0] = 1.0
            z_itr = norm.ppf(temp)
            z_itr[z_itr >= self._k] = self._k
            z_itr[z_itr <= -self._k] = -self._k
            z[..., itr] = z_itr
        mvn = multivariate_normal(mean=np.zeros(self.ndim),
                                  cov=self.corr_matrix_z)
        return mvn.cdf(z)

    # ------------------------------------------------------------------------------------------------------------------

    def sampler(self, size, random_state, method, return_all=False):
        if method == 'lhs':
            if not np.isscalar(size):
                n_samples = size[0]
                if len(size) == 2:
                    n_set = size[1]
                elif len(size) > 2:
                    raise ValueError("Non-scalar size for this type of sampling method "
                                     "should not be of length more than 2 and should represent "
                                     "(number of realizations of random vector, number of sets of realizations)!")
                else:
                    n_set = 1
            else:
                n_samples = size
                n_set = 1
        else:
            n_samples = None
            n_set = None

        if np.isscalar(size):
            x = np.zeros((size, self.ndim))
        else:
            x = np.zeros((*size, self.ndim))

        if method == 'lhs':
            for itr_set in range(n_set):
                if random_state is None:
                    rng_seed = None
                else:
                    if itr_set > 0:
                        rng_seed = int(''.join([str(i) for i in [random_state, itr_set]]))
                    else:
                        rng_seed = random_state
                criterion = "m"
                temp = lhs(self.ndim, n_samples, random_state=rng_seed, criterion=criterion)
                if len(x.shape) == 3:
                    x[:, itr_set, :] = temp
                else:
                    x = temp

        elif method == 'mcs':
            np.random.seed(random_state)
            x = np.random.uniform(0.0, 1.0, size=x.shape)
        else:
            return

        if return_all:
            return self.transform_samples(x)
        else:
            return self.transform_samples(x)[-1]

    def transform_samples(self, x):
        transformed_samples = np.zeros_like(x)
        independent_samples = norm.ppf(x)
        dependent_samples = independent_samples @ np.diag(np.sqrt(self.eig_val)) @ self.eig_vec.T
        uniform_dependent_samples = norm.cdf(dependent_samples)

        for itr_dim in range(self.ndim):
            transformed_samples[..., itr_dim] = self.prob_dist_list[itr_dim].ppf(
                uniform_dependent_samples[..., itr_dim])
        return x, independent_samples, dependent_samples, uniform_dependent_samples, transformed_samples

    def get_2d_lhs_grid(self, n, n_cont):
        eps = sys.float_info.epsilon
        x_vec = np.linspace(eps, 1 - eps, n + 1)
        x_vec_grid = np.linspace(eps, 1 - eps, n_cont)
        for itr in range(len(x_vec)):
            if not np.any(x_vec_grid == x_vec[itr]):
                x_vec_grid = np.hstack(
                    (x_vec_grid[x_vec_grid < x_vec[itr]],
                     x_vec[itr],
                     x_vec_grid[x_vec_grid > x_vec[itr]]))
        ind = [np.where(x_vec_grid == item)[0][0] for item in x_vec]
        xx, yy = np.meshgrid(x_vec_grid, x_vec_grid)
        x = np.zeros((*xx.shape, self.ndim))
        x[..., 0] = xx
        x[..., 1] = yy
        return self.transform_samples(x), ind

    ####################################################################################################################
    # Static methods
    ####################################################################################################################

    @staticmethod
    def _nataf_func_solve_1(x, *args):
        rho_ij_prime = x
        rho_ij, term_1, z, z_i, z_j = args
        n = len(z)
        mean = np.array([0, 0])
        cov = np.array([[1, rho_ij_prime], [rho_ij_prime, 1]])
        mvn = multivariate_normal(mean=mean, cov=cov)
        term_2 = mvn.pdf(np.column_stack([z_i.flatten(), z_j.flatten()])).reshape(n, n)
        integrand = term_1 * term_2
        integral = np.trapz(np.trapz(integrand, z, axis=1), z)
        return integral - rho_ij

    # ------------------------------------------------------------------------------------------------------------------

    @staticmethod
    def _nataf_func_solve_2(x, *args):
        rho_ij_prime = x
        rho_ij, term_1, eta, xi = args
        mvn = multivariate_normal(mean=np.array([0, 0]), cov=np.array([[1, rho_ij_prime], [rho_ij_prime, 1]]))
        integral = np.sum(term_1 * mvn.pdf(np.column_stack([eta, xi])))
        return integral - rho_ij

    # ------------------------------------------------------------------------------------------------------------------

    @staticmethod
    def corr_clipped(corr, threshold=1e-15):
        x_new, clipped = FrozenMultivariateNataf.clip_evals(corr, value=threshold)
        if not clipped:
            return corr

        # cov2corr
        x_std = np.sqrt(np.diag(x_new))
        x_new = x_new / np.outer(x_std, x_std)
        return x_new

    # ------------------------------------------------------------------------------------------------------------------

    @staticmethod
    def clip_evals(x, value=0.0):
        evals, evecs = np.linalg.eigh(x)
        clipped = np.any(evals < value)
        x_new = np.dot(evecs * np.maximum(evals, value), evecs.T)
        return x_new, clipped

    # ------------------------------------------------------------------------------------------------------------------

    @staticmethod
    def get_total_correlation(corr_matrix, total_mean, total_sig, mean_all, sig_all, weights):
        # corr_matrix : (n_rv, n_rv, n_s)
        # total_mean : (n_rv,)
        # total_sig : (n_rv,)
        # mean_all : (n_s, n_rv)
        # sig_all : (n_s, n_rv)
        # weights : (n_s,)

        weights = weights[np.newaxis, np.newaxis, :]
        term_1 = corr_matrix * np.einsum('ki,kj->ijk', sig_all, sig_all)
        temp = mean_all - total_mean[np.newaxis, :]
        term_2 = np.einsum('ki,kj->ijk', temp, temp)
        total_covariance = np.sum((term_1 + term_2) * weights, axis=2)
        total_corr = total_covariance / np.einsum('i,j->ij', total_sig, total_sig)
        # zero_corr_mask = corr_matrix[:, :, 0] == 0
        # total_corr[zero_corr_mask] = 0.0
        return total_corr

    # ------------------------------------------------------------------------------------------------------------------

    @staticmethod
    def get_corr_matrix_z(prob_dist_list, corr_matrix, k):
        n_dim = len(prob_dist_list)
        corr_matrix_z = np.eye(n_dim)

        # method 1
        # z = np.linspace(-6.5, 6.5, 100)
        # # z = np.linspace(-8, 8, 1024)
        # z_i, z_j = np.meshgrid(z, z)
        #

        # method 2
        n = 100
        z_max = k
        z_min = -z_max
        points, weights = np.polynomial.legendre.leggauss(n)
        points = - (0.5 * (points + 1) * (z_max - z_min) + z_min)
        weights = weights * (0.5 * (z_max - z_min))
        xi = np.tile(points, [n, 1])
        xi = xi.flatten(order='F')
        eta = np.tile(points, n)
        first = np.tile(weights, n)
        first = np.reshape(first, [n, n])
        second = np.transpose(first)
        weights2d = first * second
        w2d = weights2d.flatten()
        #

        for i in range(n_dim):
            for j in range(n_dim):
                if i > j and corr_matrix[i, j] != 0.0:
                    mu_i = prob_dist_list[i].mean()
                    mu_j = prob_dist_list[j].mean()
                    sig_i = prob_dist_list[i].std()
                    sig_j = prob_dist_list[j].std()

                    # method 1
                    # x_i = prob_dist_list[i].ppf(norm.cdf(z_i))
                    # x_j = prob_dist_list[j].ppf(norm.cdf(z_j))
                    # term_1 = ((x_i - mu_i) / sig_i) * ((x_j - mu_j) / sig_j)
                    # args_1 = (corr_matrix[i, j], term_1, z, z_i, z_j)
                    #

                    # method 2
                    x_i = prob_dist_list[i].ppf(norm.cdf(points))  # x_i = ...(eta), do to x_i what you did to eta
                    x_j = prob_dist_list[j].ppf(norm.cdf(points))  # x_j = ...(xi), do to x_j what you did to xi
                    x_i = np.tile(x_i, n)
                    x_j = np.tile(x_j, [n, 1]).flatten(order='F')
                    term_1 = ((x_i - mu_i) / sig_i) * ((x_j - mu_j) / sig_j) * w2d
                    args_2 = (corr_matrix[i, j], term_1, eta, xi)
                    #

                    solve = True
                    if prob_dist_list[i].dist.name in ['norm', 'truncnorm']:
                        if prob_dist_list[j].dist.name in ['norm', 'truncnorm']:
                            corr_matrix_z[i, j] = corr_matrix[i, j]
                            solve = False
                        if prob_dist_list[j].dist.name == 'lognorm':
                            cov_j = sig_j / np.abs(mu_j)
                            corr_matrix_z[i, j] = cov_j / np.sqrt(np.log(1 + (cov_j ** 2))) * corr_matrix[i, j]
                            solve = False
                    if prob_dist_list[i].dist.name in ['lognorm']:
                        if prob_dist_list[j].dist.name in ['norm', 'truncnorm']:
                            cov_i = sig_i / np.abs(mu_i)
                            corr_matrix_z[i, j] = cov_i / np.sqrt(np.log(1 + (cov_i ** 2))) * corr_matrix[i, j]
                            solve = False
                        if prob_dist_list[j].dist.name == 'lognorm':
                            cov_i = sig_i / np.abs(mu_i)
                            cov_j = sig_j / np.abs(mu_j)
                            fac_num = np.log(1 + corr_matrix[i, j] * cov_i * cov_j)
                            fac_den = corr_matrix[i, j] * np.sqrt(np.log(1 + (cov_i ** 2)) * np.log(1 + (cov_j ** 2)))
                            corr_matrix_z[i, j] = (fac_num / fac_den) * corr_matrix[i, j]
                            solve = False
                    if solve:
                        # sol_1 = root_scalar(FrozenMultivariateNataf._nataf_func_solve_1, args=args_1,
                        #                     method='secant', x0=args_1[0], x1=args_1[0] * 0.9)
                        sol_2 = root_scalar(FrozenMultivariateNataf._nataf_func_solve_2, args=args_2,
                                            method='secant', x0=args_2[0], x1=args_2[0] * 0.9)
                        corr_matrix_z[i, j] = sol_2.root
                    corr_matrix_z[j, i] = corr_matrix_z[i, j]
        return corr_matrix_z

    # ------------------------------------------------------------------------------------------------------------------
