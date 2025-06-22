# -*- coding: utf-8 -*-
"""Boore and Atkinson (2008) model."""

import numpy as np
import os
import pygmm
import pygmm.model as model
from pygmm.types import ArrayLike

__author__ = 'Angshuman Deb'


class BooreAtkinson2008(model.GroundMotionModel):
    """Boore and Atkinson (2008) model.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    """

    NAME = 'Boore and Atkinson (2008)'
    ABBREV = 'BA08'

    # Load the coefficients for the model
    fname = os.path.join(os.path.dirname(__file__), 'data', 'boore_atkinson_2008.csv')
    COEFF = np.recfromcsv(fname, skip_header=2, case_sensitive=True).view(np.recarray)
    PERIODS = COEFF['period']
    
    INDEX_PGV = 0
    INDEX_PGA = 1
    INDICES_PSA = np.arange(INDEX_PGA + 1, len(PERIODS))

    PARAMS = [model.NumericParameter('mag', True),
              model.NumericParameter('dist_jb', True),
              model.NumericParameter('v_s30', True),
              model.CategoricalParameter('mechanism', False, ['U', 'SS', 'NS', 'RS'], 'U')
              ]

    def __init__(self, scenario: model.Scenario):
        """Initialize the model."""
        super(BooreAtkinson2008, self).__init__(scenario)
        pga_ref = np.exp(self._calc_ln_resp(np.nan)[self.INDEX_PGA])
        self._ln_resp = self._calc_ln_resp(pga_ref)
        self._ln_std = self._calc_ln_std()

    def _calc_ln_resp(self, pga_ref: ArrayLike) -> np.ndarray:
        """Calculate the natural logarithm of the response.

        Returns
        -------
        ln_resp : class:`np.array`:
            natural log of the response
        """
        s = self._scenario
        c = self.COEFF
        
        # Compute the event term
        ########################
        if s.mechanism == 'SS':
            event = np.array(c.e_02)
        elif s.mechanism == 'NS':
            event = np.array(c.e_03)
        elif s.mechanism == 'RS':
            event = np.array(c.e_04)
        else:
            # Unspecified
            event = np.array(c.e_01)

        mask = s.mag <= c.m_h
        event[mask] += (c.e_05 * (s.mag - c.m_h) + c.e_06 * (s.mag - c.m_h) ** 2)[mask]
        event[~mask] += (c.e_07 * (s.mag - c.m_h))[~mask]
        
        # Compute the distance terms
        ############################
        dist = np.sqrt(s.dist_jb ** 2 + c.h ** 2)
        path = (c.c_01 + c.c_02 * (s.mag - c.m_ref)) * np.log(dist / c.r_ref) + c.c_03 * (dist - c.r_ref)
        
        if np.isnan(pga_ref):
            # Reference condition. No site effect
            site = 0
        else:
            # Compute the site term
            f_lin = c.b_lin * np.log(s.v_s30 / c.v_ref)
            
            # Add the nonlinearity to the site term
            mask_1 = s.v_s30 <= c.v_1
            mask_2 = np.all([s.v_s30 > c.v_1, s.v_s30 <= c.v_2], axis=0) 
            mask_3 = np.all([s.v_s30 > c.v_2, s.v_s30 <= c.v_ref], axis=0) 
            
            b_nl = np.zeros(len(c.b_lin))
            b_nl[mask_1] = c.b_1[mask_1]
            b_nl[mask_2] = (c.b_2 + (c.b_1 - c.b_2) * np.log(s.v_s30 / c.v_2) / np.log(c.v_1 / c.v_2))[mask_2]
            b_nl[mask_3] = (c.b_2 * np.log(s.v_s30 / c.v_ref) / np.log(c.v_2 / c.v_ref))[mask_3]
            
            delta_x = np.log(c.a_2 / c.a_1)
            delta_y = b_nl * np.log(c.a_2 / c.pga_low)
            c_term = (3.0 * delta_y - b_nl * delta_x) / (delta_x ** 2)
            d_term = - (2.0 * delta_y - b_nl * delta_x) / (delta_x ** 3)
            
            mask_1 = pga_ref <= c.a_1
            mask_2 = np.all([pga_ref > c.a_1, pga_ref <= c.a_2], axis=0)
            mask_3 = pga_ref > c.a_2
            
            f_nl = np.zeros(len(f_lin))
            f_nl[mask_1] = (b_nl * np.log(c.pga_low / 0.1))[mask_1]
            f_nl[mask_2] = (b_nl * np.log(c.pga_low / 0.1) +
                            c_term * (np.log(pga_ref / c.a_1)) ** 2 +
                            d_term * (np.log(pga_ref / c.a_1)) ** 3)[mask_2]
            f_nl[mask_3] = (b_nl * np.log(pga_ref / 0.1))[mask_3]
            
            site = f_lin + f_nl

        ln_resp = event + path + site
        return ln_resp

    def _calc_ln_std(self) -> np.ndarray:
        """Calculate the logarithmic standard deviation.

        Returns
        -------
        ln_std : class:`np.array`:
            natural log standard deviation

        """
        c = self.COEFF
        s = self._scenario
        
        return (s.mechanism == 'U') * c.sig_tu + (s.mechanism != 'U') * c.sig_tm


pygmm.__all__.append('BooreAtkinson2008')
pygmm.models.append(BooreAtkinson2008)
