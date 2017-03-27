#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2017 Srinivasa Rao Zinka
# License: New BSD License.

from __future__ import division
import numpy as np
import filtertool as ft

M = np.zeros((8, 8))
Rs = Rl = 1.02828
M[0,1] = M[1,0] = 0.80799
M[0,3] = M[3,0] = -0.10066
M[1,2] = M[2,1] = 0.6514
M[2,3] = M[3,2] = 0.52837
M[3,4] = M[4,3] = 0.53064
M[4,5] = M[5,4] = 0.49184
M[4,7] = M[7,4] = -0.2346
M[5,6] = M[6,5] = 0.73748
M[6,7] = M[7,6] = 0.77967

print(M)
ft.MN_to_Sparam(M, Rs, Rl, w_min= -3, w_max=3, w_num=500, dB=True, dB_limit= -100)

# Actual coupling coef and quality factors
FBW = 0.07063
M = M*FBW
Rs = Rs/FBW
Rl = Rl/FBW
print(M)
print(Rs)
