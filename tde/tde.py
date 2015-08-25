from __future__ import division
import numpy as np
from numpy import sin, cos, pi

import strain_functions as sf

try:
    import numexpr as ne
    ne_flag = True
except ImportError:
    ne_flag = False
    pass

try:
    from numba import jit, vectorize
    nu_flag = True
except:
    nu_flag = False
    pass

#Vec = namedtuple('Vec', ['xx','yy','zz','xy','xz','yz'])

def swap(a, b):
    b, a  =  a, b

    return b, a


def calc_tri_strains(sx=None, sy=None, sz=None, x=None, y=None, 
                     z=None, pr=0.25, ss=0., ts=0., ds=0.):
    '''
    docs
    '''

    slip_vec = calc_slip_vector(x, y, z)

    S = {
        'xx': np.zeros(len(sx)),
        'yy': np.zeros(len(sx)),
        'zz': np.zeros(len(sx)),
        'xy': np.zeros(len(sx)),
        'xz': np.zeros(len(sx)),
        'yz': np.zeros(len(sx)),
    }

    x = np.append(x, x[0])
    y = np.append(y, y[0])
    z = np.append(z, z[0])

    for i_tri in [0, 1, 2]:
        strike, dip, lss, ts, lds = get_edge_params(i_tri, x, y, z, slip_vec)

    return S


def calc_slip_vector(x, y, z, ss=0., ts=0., ds=0.):
    # Calculate the slip vector in XYZ coordinates
    v0 = np.array([x[0], y[0], z[0]])
    v1 = np.array([x[1], y[1], z[1]])
    v2 = np.array([x[2], y[2], z[2]])

    norm_vec = np.cross( (v1 - v0), (v2-v0))
    norm_vec /= np.linalg.norm(norm_vec)

    if norm_vec[2] < 0: # Enforce clockwise circulation
        norm_vec *= -1
        x[1], x[2] = x[2], x[1]
        y[1], y[2] = y[2], y[1]
        z[1], z[2] = z[2], z[1]

    strike_vec = np.array([-np.sin( np.arctan2( norm_vec[1], norm_vec[0])),
                            np.cos( np.arctan2( norm_vec[1], norm_vec[0])), 0])
    dip_vec = np.cross( norm_vec, strike_vec)
    slip_comp = np.array([ss, ds, ts])
    slip_vec = np.array([strike_vec[:], dip_vec[:], norm_vec[:]]).dot(slip_comp)

    return slip_vec


def get_edge_params(i_tri, x, y, z, slip_vec):

    # calculate strike and dip (trend and plunge?) of current leg
    x_dist = x[i_tri + 1] - x[i_tri]
    y_dist = y[i_tri + 1] - y[i_tri]
    z_dist = z[i_tri + 1] - z[i_tri]

    strike_rad = np.arctan2( y_dist, x_dist)
    strike = 180. / np.pi * strike_rad

    seg_map_length = np.sqrt(x_dist**2 + y_dist**2)
    rx, ry = rotate_xy_vec(x_dist, y_dist, -strike)
    dip = 180. / np.pi * np.arctan2( z_dist, rx)

    if dip >= 0.:
        beta = np.pi / 180. * (90. - dip)
        if beta > (np.pi/2.):
            beta = np.pi/2. - beta
    else:
        beta = -np.pi/180. * (90. + dip)
        if beta < (-np.pi/2.):
            beta = np.pi/2. - np.abs(beta)

    ss_vec = np.array([ np.cos(strike_rad), np.sin(strike_rad), 0.])
    ts_vec = np.array([-np.sin(strike_rad), np.cos(strike_rad), 0.])
    ds_vec = np.cross(ss_vec, ts_vec)
    lss = np.dot(slip_vec, ss_vec)
    lts = np.dot(slip_vec, ts_vec)
    lds = np.dot(slip_vec, ds_vec)

    return strike, dip, lss, lts, lds

def get_edge_strains():

    sx1, sy1 = rotate_xy_vec(sx-x[i_tri], sy-y[i_tri], -strike)
    



def rotate_xy_vec(x, y, alpha):
    '''Rotates a vector by an ange alpha'''
    alpha_rad = np.pi / 180. * alpha

    xp = np.cos(alpha_rad) * x - np.sin(alpha_rad) * y
    yp = np.sin(alpha_rad) * x + np.cos(alpha_rad) * y

    return xp, yp

def advs(y1, y2, y3, a, b, nu, B1, B2, B3):
    

    E = {
        'e11': 0.,
        'e22': 0.,
        'e33': 0.,
        'e12': 0.,
        'e13': 0.,
        'e23': 0.
    }
    return E


'''
.* -> _*_ (_ is space)
./ -> _/_ (_ is space)
.^ -> ** 
sin -> np.sin
cos -> np.cos
pi -> np.pi
'''

def cot(x):
    return 1 / np.tan(x)



