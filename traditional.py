#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2017 Srinivasa Rao Zinka
# License: New BSD License.

from __future__ import division
import numpy as np

def Butterworth_g(LAs, Ws):
    n = np.ceil(np.log(10**(0.1 * LAs) - 1) / (2 * np.log(Ws)))
    i = np.arange(n + 2)
    i1 = (2 * i - 1) * np.pi / 2 / n
    g = 2 * np.sin(i1)
    g[0] = 1
    g[-1] = 1
    return n, g

def Chebyshev_g(LAr, LAs, Ws):
    n = int(np.ceil(np.arccosh(
        np.sqrt((10**(0.1 * LAs) - 1) / (10**(0.1 * LAr) - 1))) / np.arccosh(Ws)))
    beta = np.log(1 / np.tanh(LAr / 17.37))
    gamma = np.sinh(beta / 2 / n)
    g = np.ones(n + 2)
    g[1] = (2 / gamma) * np.sin(np.pi / 2 / n)
    g[-1] = 1 / np.tanh(beta / 4)**2 if n % 2 == 0 else 1

    for i in range(2, n + 1):
        num = 4 * np.sin((2 * i - 1) * np.pi / 2 / n) * \
            np.sin((2 * i - 3) * np.pi / 2 / n)
        den = gamma**2 + np.sin((i - 1) * np.pi / n)**2
        g[i] = num / den / g[i - 1]
    return n, g

def g2LC(g, Z0=50, fc=1e9, ladder_type='LC'): # only for LPFs
    if(ladder_type == 'LC'):
        g[1:-1:2] = g[1:-1:2] * Z0 / (2 * np.pi * fc)
        g[2:-1:2] = g[2:-1:2] / Z0 / (2 * np.pi * fc)
        g[0] = Z0
        g[-1] = g[-1] * Z0
        g2LC = g
    if(ladder_type == 'CL'):
        g[1:-1:2] = g[1:-1:2] / Z0 / (2 * np.pi * fc)
        g[2:-1:2] = g[2:-1:2] * Z0 / (2 * np.pi * fc)
        g[0] = Z0
        g[-1] = g[-1] * Z0
        g2LC = g
    return g2LC

def stepped_impedance(g2LC, Z0_high, Z0_low, lg_high, lg_low, fc, ladder_type='LC', optimize=False):
    LC_ang = np.zeros_like(g2LC)
    if(ladder_type == 'LC'):
        g2LC[1:-1:2] = (lg_high / 2 / np.pi) * np.arcsin(2 *
                                                         np.pi * fc * g2LC[1:-1:2] / Z0_high)
        LC_ang[1:-1:2] = g2LC[1:-1:2] * (2 * np.pi / lg_high)
        g2LC[2:-1:2] = (lg_low / 2 / np.pi) * np.arcsin(2 *
                                                        np.pi * fc * g2LC[2:-1:2] * Z0_low)
        LC_ang[2:-1:2] = g2LC[2:-1:2] * (2 * np.pi / lg_low)
        LC_len = g2LC
    if(ladder_type == 'CL'):
        pass
    if(optimize):
        pass # implement (5.4) from J. S. Hong

    return LC_len, LC_ang

def stub_LPF():
    # Implement Fig. 5.3 from J. S. Hong
    pass

def g2K(g, LC_values=np.ones_like(g), Z0=1, ZL=1, ladder_type='L'):
    pass

if __name__ == '__main__':

    # print(Butterworth_g(LAs=40, Ws=2))
    # Chebyshev_g(LAr=0.1, LAs=40, Ws=2)

    # ==========================================================================
    # J. S. Hong, Ch.5, LPF Excercise
    # ==========================================================================

    # g = Chebyshev_g(LAr=0.1, LAs=40, Ws=6)[1]
    # g2LC = g2LC(g, Z0=50, fc=1e9, ladder_type='LC')
    # LC_len, LC_ang = stepped_impedance(
    #     g2LC, Z0_high=93, Z0_low=24, lg_high=118e-3, lg_low=105e-3, fc=1e9)
    # print(LC_ang * 180 / np.pi)
    # print(LC_len)

    # ==========================================================================
    # J. S. Hong, Ch.4, g to K conversion, Fig. 3.18
    # ==========================================================================

    g = Chebyshev_g(LAr=0.1, LAs=40, Ws=6)[1]



    # ==========================================================================
    # Notes Section
    # ==========================================================================
    # TODO implement elliptic from QUCS
    # TODO
