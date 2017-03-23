#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2017 Srinivasa Rao Zinka
# License: New BSD License.

from __future__ import division
import numpy as np
# import sympy as sp
# import scipy.signal as signal
# import matplotlib.pyplot as plt
# import warnings
# import Tkinter as ti
# import tkFileDialog as tkdlg

# adjusting "matplotlib" label fonts
from matplotlib import rc
rc('text', usetex=True)

import filtertool as ft

N = 7
poles = np.array([1.2576, -0.1546 - 0.9218j, -0.1546 + 0.9218j])
#     poles = np.array([-3,3])
eps = 6.0251j; eps_R = 1
#    poles = np.array([])
F, P = ft.Chebyshev_gen(N, poles)
ft.plot_rational(F, P, x_min= -1.1, x_max=1.1, x_num=1000)

F = ft.w_to_s(F, coef_norm=True)
P = ft.w_to_s(P, coef_norm=True)

print 'F:', '\n', F; print 'P:', '\n', P
E, roots_E = ft.poly_E(eps, eps_R, F, P)
print 'E:', '\n', E

#    F1, P1, E1 = dual_band(F, P, E, 0.4)

ft.plot_mag(eps, eps_R, F, P, E, w_min= -2, w_max=2, w_num=500, dB=True,
         dB_limit= -40, plot=True)
ft.plot_delay(roots_E)

# From now onwards, unlike the Cameron's example, this filter is doubly terminated
M, Rs, Rl = ft.coupling_N(F, P, E, eps, eps_R)
print 'M:', '\n', M.real
print 'Rs:', Rs
print 'Rl:', Rl

ft.MN_to_Sparam(M, Rs, Rl, w_min= -3, w_max=3, w_num=500, dB=True, dB_limit= -40)
