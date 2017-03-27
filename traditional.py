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
    return g

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
    return g

if __name__ == '__main__':

    print(Butterworth_g(LAs=40, Ws=2))
    print(Chebyshev_g(LAr=0.1, LAs=40, Ws=2))

# ==============================================================================
# Notes Section
# ==============================================================================
# TODO implement elliptic from QUCS  
