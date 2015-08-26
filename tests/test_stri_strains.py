import sys; sys.path.append('../tde')
import numpy as np
import matplotlib.pyplot as plt
import tde

# tri dislocation parameters
pr = 0.25
ss = -1.
ts = 0.
ds = 0.
N = 30

sx, sy, sz = np.meshgrid( np.linspace(0, 100, N), np.linspace(0, 100, N), 0)

X = np.array([40., 60., 40.])
Y = np.array([50., 50., 30.])
Z = np.array([0., 0., 20.])

S = tde.calc_tri_strains(sx.ravel(order='F'), sy.ravel(order='F'), 
                         sz.ravel(order='F'), X, Y, Z, pr, ss, ts, ds)


