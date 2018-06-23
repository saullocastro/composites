"""
Composite Lamina Module (:mod:`composites.lamina`)
==========================================================

.. currentmodule:: composites.lamina

"""
from __future__ import division, absolute_import

import numpy as np
from numpy import cos, sin

from .matlamina import MatLamina


class Lamina(object):
    """
    =========  ===========================================================
    attribute  description
    =========  ===========================================================
    plyid      id of the composite lamina
    matobj     a pointer to a MatLamina object
    h          ply thickness
    theta      ply angle in degrees
    L          transformation matrix for displacements to laminate csys
    R          transformation matrix for stresses to laminate csys
    T          transformation matrix for stresses to lamina csys
    QL         constitutive matrix for plane-stress in laminate csys
    =========  ===========================================================

    References:
    -----------
    .. [1] Reddy, J. N., Mechanics of Laminated Composite Plates and
       Shells - Theory and Analysys. Second Edition. CRC PRESS, 2004.

    """
    def __init__(self):
        self.plyid = None
        self.matobj = None
        self.h = None
        self.theta = None
        self.L = None
        self.R = None
        self.T = None
        self.QL = None

    def rebuild(self):
        thetarad = np.deg2rad(self.theta)
        cost = cos(thetarad)
        sint = sin(thetarad)
        sin2t = sin(2*thetarad)

        cos2   = cost**2
        cos3   = cost**3
        cos4   = cost**4
        sin2   = sint**2
        sin3   = sint**3
        sin4   = sint**4
        sincos = sint*cost
        self.L = np.array([[ cost, sint, 0],
                           [-sint, cost, 0],
                           [   0,     0, 1]], dtype=np.float64)
        #STRESS
        #to lamina
        self.R = np.array(
            [[   cos2,   sin2, 0,   0,    0,     sin2t],
             [   sin2,   cos2, 0,   0,    0,    -sin2t],
             [      0,      0, 1,   0,    0,         0],
             [      0,      0, 0, cost, -sint,         0],
             [      0,      0, 0, sint,  cost,         0],
             [-sincos, sincos, 0,   0,    0, cos2-sin2]],dtype=np.float64)
        #to laminate
        self.T = np.array(
            [[  cos2,    sin2, 0,    0,   0,    -sin2t],
             [  sin2,    cos2, 0,    0,   0,     sin2t],
             [     0,       0, 1,    0,   0,         0],
             [     0,       0, 0,  cost, sint,         0],
             [     0,       0, 0, -sint, cost,         0],
             [sincos, -sincos, 0,    0,   0, cos2-sin2]],dtype=np.float64)
        # STRAINS
        # different from stress due to:
        #     2*e12 = e6    2*e13 = e5    2*e23 = e4
        # to laminate
        # self.Rstrain = np.transpose(self.Tstress)
        # to lamina
        # self.Tstrain = np.transpose(self.Rstress)
        if isinstance(self.matobj, MatLamina):
            e1   = self.matobj.e1
            e2   = self.matobj.e2
            nu12 = self.matobj.nu12
            nu21 = self.matobj.nu21
            g12  = self.matobj.g12
            g13  = self.matobj.g13
            g23  = self.matobj.g23
        else:
            e1   = self.matobj.e
            e2   = self.matobj.e
            nu12 = self.matobj.nu
            nu21 = self.matobj.nu
            g12  = self.matobj.g

        # plane stress
        q11  = e1/(1-nu12*nu21)
        q12  = nu12*e2/(1-nu12*nu21)
        q22  = e2/(1-nu12*nu21)
        q44  = g23
        q55  = g13
        q16 = 0
        q26 = 0
        q66  = g12

        q11L = q11*cos4 + 2*(q12 + 2*q66)*sin2*cos2 + q22*sin4
        q12L = (q11 + q22 - 4*q66)*sin2*cos2 + q12*(sin4 + cos4)
        q22L = q11*sin4 + 2*(q12 + 2*q66)*sin2*cos2 + q22*cos4
        q16L = (q11 - q12 - 2*q66)*sint*cos3 + (q12 - q22 + 2*q66)*sin3*cost
        q26L = (q11 - q12 - 2*q66)*sin3*cost + (q12 - q22 + 2*q66)*sint*cos3
        q66L = (q11 + q22 - 2*q12 - 2*q66)*sin2*cos2 + q66*(sin4 + cos4)
        q44L = q44*cos2 + q55*sin2
        q45L = (q55 - q44)*sincos
        q55L = q55*cos2 + q44*sin2

        self.QL = np.array([[q11L, q12L, q16L,    0,    0],
                            [q12L, q22L, q26L,    0,    0],
                            [q16L, q26L, q66L,    0,    0],
                            [   0,    0,    0, q44L, q45L],
                            [   0,    0,    0, q45L, q55L]], dtype=np.float64)

        #TODO add the thermal coeficient terms when calculating the
        #     stresses... to take into account eventual thermal expansions /
        #     contractions

