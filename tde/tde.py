from __future__ import division
import numpy as np
from numpy import sin, cos, pi

import disp_functions as df

try:
    import strain_functions_cy as sf
except ImportError:
    import strain_functions as sf


'''
Functions modified from Brendan Meade's "tde" MATLAB routines
(https://github.com/brendanjmeade/tde) detailed in
Computers and Geosciences (2006).

I would love to have more thorough documentation but I didn't spend the time
to follow the derivation well enough to know what is really going on. In
most cases, names of variables, functions and so forth either follow
directly from Meade's code, or is a guess of mine what is going on.

Ported by Richard Styron, 2015
'''


eps = 7./3 - 4./3 - 1 # get machine epsilon

def calc_tri_strains(sx=None, sy=None, sz=None, x=None, y=None, 
                     z=None, pr=0.25, ss=0., ts=0., ds=0.):
    '''
    Calculates strains due to slip on a triangular dislocation in an
    elastic half space utilizing the symbolically differentiated
    displacement gradient tensor derived from the expressions for
    the displacements due to an angular dislocation in an elastic half
    space (Comninou and Dunders, 1975).
    
    Arguments
     sx : x-coordinates of observation points
     sy : y-coordinates of observation points
     sz : z-coordinates of observation points
     x  : x-coordinates of triangle vertices.
     y  : y-coordinates of triangle vertices.
     z  : z-coordinates of triangle vertices.
     pr : Poisson's ratio
     ss : strike slip displacement
     ts : tensile slip displacement
     ds : dip slip displacement
    
    Returns
     E  : structure containing the strains (E.xx, E.yy, E.zz, E.xy, E.xz, E.yz)
    
    This paper should and related code should be cited as:
    Brendan J. Meade, Algorithms for the calculation of exact 
    displacements, strains, and stresses for Triangular Dislocation 
    Elements in a uniform elastic half space, Computers & 
    Geosciences (2007), doi:10.1016/j.cageo.2006.12.003.
    
    Use at your own risk and please let me know of any bugs/errors.
    
    Copyright (c) 2006 Brendan Meade

    Ported to Python by Richard Styron 2015.

    '''

    xx = np.copy(x)
    yy = np.copy(y)
    zz = np.copy(z)

    slip_vec = calc_slip_vector(xx, yy, zz, ss, ts, ds)

    E = {
        'xx': np.zeros(len(sx)),
        'yy': np.zeros(len(sx)),
        'zz': np.zeros(len(sx)),
        'xy': np.zeros(len(sx)),
        'xz': np.zeros(len(sx)),
        'yz': np.zeros(len(sx))
    }

    xx = np.append(xx, xx[0]) # for indexing during loops
    yy = np.append(yy, yy[0])
    zz = np.append(zz, zz[0])

    for i_tri in [0, 1, 2]:
        strike, dip, beta, lss, lts, lds = get_edge_params(i_tri, xx, yy, zz, 
                                                          slip_vec)
        
        e = get_edge_strains(sx, sy, sz, xx, yy, zz, i_tri, beta, pr, lss, lts,
                             lds, strike)

        E['xx'] += e['11']
        E['yy'] += e['22']
        E['zz'] += e['33']
        E['xy'] += e['12']
        E['xz'] += e['13']
        E['yz'] += e['23']

    return E


def calc_tri_displacements(sx=None, sy=None, sz=None, x=None, y=None, z=None, 
                           pr=0.25, ss=0., ts=0., ds=0.):
    '''
    Calculates displacements due to slip on a triangular dislocation in an
    elastic half space utilizing the Comninou and Dunders (1975) expressions
    for the displacements due to an angular dislocation in an elastic half
    space.
    
    Arguments
     sx : x-coordinates of observation points
     sy : y-coordinates of observation points
     sz : z-coordinates of observation points
     x  : x-coordinates of triangle vertices.
     y  : y-coordinates of triangle vertices.
     z  : z-coordinates of triangle vertices.
     pr : Poisson's ratio
     ss : strike slip displacement
     ts : tensile slip displacement
     ds : dip slip displacement
    
    Returns
     U  : Dictionary containing the displacements (U['x'], U['y'], U['z'])
    
    This paper should and related code should be cited as:
    Brendan J. Meade, Algorithms for the calculation of exact 
    displacements, strains, and stresses for Triangular Dislocation 
    Elements in a uniform elastic half space, Computers & 
    Geosciences (2007), doi:10.1016/j.cageo.2006.12.003.

    That code Copyright (c) 2006 Brendan Meade

    Ported to Python by Richard Styron, Aug 2015

    '''

    xx = np.copy(x)
    yy = np.copy(y)
    zz = np.copy(z)


    slip_vec = calc_slip_vector(xx, yy, zz, ss, ts, ds)

    U = {
        'x': np.zeros(len(sx)),
        'y': np.zeros(len(sx)),
        'z': np.zeros(len(sx))
    }

    xx = np.append(xx, xx[0]) # for indexing during loops
    yy = np.append(yy, yy[0])
    zz = np.append(zz, zz[0])

    for i_tri in [0, 1, 2]:
        strike, dip, beta, lss, lts, lds = get_edge_params(i_tri, xx, yy, zz, 
                                                           slip_vec)
        
        uxn, uyn, uzn = get_edge_displacements(sx, sy, sz, xx, yy, zz, i_tri, 
                                              beta, pr, lss, lts, lds, strike)

        U['x'] += uxn
        U['y'] += uyn
        U['z'] += uzn

    U = offset_underlying_points(xx, yy, zz, sx, sy, sz, slip_vec, U)

    return U



def calc_slip_vector(x, y, z, ss=0., ts=0., ds=0.):
    '''
    Calculates a 3D slip vector for displacement on a triangle, given
    strike-slip, tensile-slip and dip-slip displacement components.

    Arguments
     x  : x-coordinates of triangle vertices.
     y  : y-coordinates of triangle vertices.
     z  : z-coordinates of triangle vertices.
     ss : strike slip displacement
     ts : tensile slip displacement
     ds : dip slip displacement


    Note: Will modify the ordering of triangle coordinates variables passed,
    even outside of the function scope, if 'clockwise circulation' does not
    hold. If this function is used outside of `calc_tri_displacement` or
    `calc_tri_strains`, make sure to send it a copy of the triangle coordinate
    variables instead of the original, or strange bugs could result.
    '''
    # Calculate the slip vector in XYZ coordinates
    v0 = np.array([x[0], y[0], z[0]])
    v1 = np.array([x[1], y[1], z[1]])
    v2 = np.array([x[2], y[2], z[2]])

    norm_vec = np.cross( (v1 - v0), (v2 - v0))

    norm_vec = norm_vec / np.linalg.norm(norm_vec)

    if norm_vec[2] < 0: # Enforce clockwise circulation
        norm_vec *= -1
        x[1], x[2] = x[2], x[1]
        y[1], y[2] = y[2], y[1]
        z[1], z[2] = z[2], z[1]

    strike_vec = np.array([-np.sin( np.arctan2( norm_vec[1], norm_vec[0])),
                            np.cos( np.arctan2( norm_vec[1], norm_vec[0])), 0])
    
    dip_vec = np.cross( norm_vec, strike_vec)
    
    slip_comp = np.array([ss, ds, ts])
    
    fault_mat = np.array([[strike_vec[0], dip_vec[0], norm_vec[0]],
                          [strike_vec[1], dip_vec[1], norm_vec[1]],
                          [strike_vec[2], dip_vec[2], norm_vec[2]]])

    slip_vec = fault_mat.dot(slip_comp)

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

    return strike, dip, beta, lss, lts, lds


def get_edge_strains(sx, sy, sz, x, y, z, i_tri, beta, pr, lss, lts,lds, 
                     strike):
    
    sx1, sy1 = rotate_xy_vec(sx-x[i_tri], sy-y[i_tri], -strike)
    a = sf.advs(sx1, sy1, sz-z[i_tri], z[i_tri], beta, pr, lss, lts, lds)

    sx2, sy2 = rotate_xy_vec(sx-x[i_tri+1], sy-y[i_tri+1], -strike)
    b = sf.advs(sx2, sy2, sz-z[i_tri+1], z[i_tri+1], beta, pr, lss, lts, lds)

    bxx = a['11']-b['11'] 
    byy = a['22']-b['22'] 
    bzz = a['33']-b['33'] 
    bxy = a['12']-b['12'] 
    bxz = a['13']-b['13'] 
    byz = a['23']-b['23'] 

    g = np.pi / 180. * strike

    e = {}

    e['11'] = ((np.cos(g) * bxx - np.sin(g) * bxy) * np.cos(g) -
               (np.cos(g) * bxy - np.sin(g) * byy) * np.sin(g))
    e['12'] = ((np.cos(g) * bxx - np.sin(g) * bxy) * np.sin(g) +
               (np.cos(g) * bxy - np.sin(g) * byy) * np.cos(g))
    e['13'] = np.cos(g) * bxz - np.sin(g) * byz
    e['22'] = ((np.sin(g) * bxx + np.cos(g) * bxy) * np.sin(g) +
               (np.sin(g) * bxy + np.cos(g) * byy) * np.cos(g))
    e['23'] = np.sin(g) * bxz + np.cos(g) * byz
    e['33'] = bzz

    return e
               

def get_edge_displacements(sx, sy, sz, x, y, z, i_tri, beta, pr, lss, lts,lds, 
                           strike):

    sx1, sy1 = rotate_xy_vec(sx-x[i_tri], sy-y[i_tri], -strike)
    ux1, uy1, uz1 = df.adv(sx1, sy1, sz-z[i_tri], z[i_tri], beta, pr, lss, lts, 
                           lds)

    sx2, sy2 = rotate_xy_vec(sx-x[i_tri+1], sy-y[i_tri+1], -strike)
    ux2, uy2, uz2 = df.adv(sx2, sy2, sz-z[i_tri+1], z[i_tri+1], beta, pr, lss, 
                           lts, lds)

    uxn, uyn = rotate_xy_vec((ux1 - ux2), (uy1-uy2), strike)
    uzn = uz1 - uz2

    return uxn, uyn, uzn
    

def rotate_xy_vec(x, y, alpha):
    '''Rotates a vector by an ange alpha'''
    alpha_rad = np.pi / 180. * alpha

    xp = np.cos(alpha_rad) * x - np.sin(alpha_rad) * y
    yp = np.sin(alpha_rad) * x + np.cos(alpha_rad) * y

    return xp, yp


def offset_underlying_points(x, y, z, sx, sy, sz, slip_vec, U):

    poly = list(zip(x,y))

    for i in range(len(sx)):
        xi, yi, zi = sx[i], sy[i], sz[i]

        if line_plane_intersect(x, y, z, xi, yi, zi)[2] < 0:
            if point_in_poly(xi, yi, poly):
                U['x'][i] -= slip_vec[0]
                U['y'][i] -= slip_vec[1]
                U['z'][i] -= slip_vec[2]
    return U


def point_in_poly(x,y,poly):
    # from http://www.ariel.com.au/a/python-point-int-poly.html

    n = len(poly)
    inside = False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside


def line_plane_intersect(x, y, z, sx, sy, sz):
    '''
    Calculate the intersection of a line and a plane using a parametric
    representation of the plane. This is hardcoded for a vertical line.
    '''

    numerator = np.array([[1., 1., 1., 1.],
                          [x[0], x[1], x[2], sx],
                          [y[0], y[1], y[2], sy],
                          [z[0], z[1], z[2], sz]])

    numerator = np.linalg.det(numerator)

    denominator = np.array([[1., 1., 1., 1.],
                            [x[0], x[1], x[2], 0],
                            [y[0], y[1], y[2], 0],
                            [z[0], z[1], z[2], -sz]])

    denominator = np.linalg.det(denominator)

    if denominator == 0:
        denominator = np.spacing(1)

    t = numerator / denominator
    d = np.array([sx, sy, sz]) - t * (np.array([sx, sy, 0])- 
                                      np.array([sx, sy, sz]))
    return d

 
                          
def strain_to_stress(E, lamda, mu):
    '''Calculate stresses and invariants given a strain tensor and elastic
       moduli lambda and mu.

       Returns a dict S with the stress components and invariants.
    '''

    S = {}

    
    S['xx'] = 2 * mu * E['xx'] + lamda * (E['xx']+E['yy']+E['zz'])
    S['yy'] = 2 * mu * E['yy'] + lamda * (E['xx']+E['yy']+E['zz'])
    S['zz'] = 2 * mu * E['zz'] + lamda * (E['xx']+E['yy']+E['zz'])
    S['xy'] = 2 * mu * E['xy']
    S['xz'] = 2 * mu * E['xz']
    S['yz'] = 2 * mu * E['yz']
    
    S['I1'] = S['xx'] + S['yy'] + S['zz']
    S['I2'] = (-(S['xx'] * S['yy'] + S['yy'] * S['zz'] + S['xx'] * S['zz']) + 
               S['xy'] * S['xy'] + S['xz'] * S['xz'] + S['yz'] * S['yz'])
    S['I3'] = (S['xx'] * S['yy'] * S['zz'] + 2 * S['xy'] * S['xz'] * S['yz'] - 
               (S['xx'] * S['yz'] * S['yz'] + S['yy'] * S['xz'] * S['xz'] + 
                S['zz'] * S['xy'] * S['xy']))

    return S
