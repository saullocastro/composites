"""
Composite Laminate Module (:mod:`composites.laminate`)
==============================================================

.. currentmodule:: composites.laminate

"""
from __future__ import division, absolute_import

import numpy as np

from .lamina import Lamina
from .matlamina import read_laminaprop


def read_stack(stack, plyt=None, laminaprop=None, rho=None, plyts=None, laminaprops=None,
        rhos=None, offset=0., calc_scf=True):
    """Read a laminate stacking sequence data.

    An ``Laminate`` object is returned based on the inputs given.

    Parameters
    ----------
    stack : list
        Angles of the stacking sequence in degrees.
    plyt : float, optional
        When all plies have the same thickness, ``plyt`` can be supplied.
    laminaprop : tuple, optional
        When all plies have the same material properties, ``laminaprop``
        can be supplied.
    rho : float, optional
        Uniform material density to be used for all plies.
    plyts : list, optional
        A list of floats with the thickness of each ply.
    laminaprops : list, optional
        A list of tuples with a laminaprop for each ply.
    rhos : list, optional
        A list of floats with the material density of each ply.
    offset : float, optional
        Offset along the normal axis about the mid-surface, which influences
        the laminate properties.
    calc_scf : bool, optional
        If True, use :method:`.Laminate.calc_scf` to compute shear correction
        factors, otherwise the default value of 5/6 is used

    Notes
    -----
    ``plyt`` or ``plyts`` must be supplied
    ``laminaprop`` or ``laminaprops`` must be supplied

    For orthotropic plies, the ``laminaprop`` should be::

        laminaprop = (E11, E22, nu12, G12, G13, G23)

    For isotropic pliey, the ``laminaprop`` should be::

        laminaprop = (E, E, nu)

    """
    lam = Laminate()
    lam.offset = offset
    lam.stack = stack

    if plyts is None:
        if plyt is None:
            raise ValueError('plyt or plyts must be supplied')
        else:
            plyts = [plyt for i in stack]

    if laminaprops is None:
        if laminaprop is None:
            raise ValueError('laminaprop or laminaprops must be supplied')
        else:
            laminaprops = [laminaprop for i in stack]

    if rhos is None:
        rhos = [rho for i in stack]

    lam.plies = []
    for plyt, laminaprop, theta, rho in zip(plyts, laminaprops, stack, rhos):
        laminaprop = laminaprop
        ply = Lamina()
        ply.theta = float(theta)
        ply.h = plyt
        ply.matobj = read_laminaprop(laminaprop, rho)
        lam.plies.append(ply)

    lam.rebuild()
    lam.calc_constitutive_matrix()
    if calc_scf:
        lam.calc_scf()

    return lam


def read_lamination_parameters(thickness, laminaprop, rho,
                               xiA1, xiA2, xiA3, xiA4,
                               xiB1, xiB2, xiB3, xiB4,
                               xiD1, xiD2, xiD3, xiD4,
                               xiE1, xiE2, xiE3, xiE4):
    r"""Calculates a laminate based on the lamination parameters.

    The lamination parameters:
    `\xi_{A1} \cdots \xi_{A4}`,  `\xi_{B1} \cdots \xi_{B4}`,
    `\xi_{C1} \cdots \xi_{C4}`,  `\xi_{D1} \cdots \xi_{D4}`,
    `\xi_{E1} \cdots \xi_{E4}`

    are used to calculate the laminate constitutive matrix.

    Parameters
    ----------
    thickness : float
        The total thickness of the laminate
    laminaprop : tuple
        The laminaprop tuple used to define the laminate material.
    rho : float
        Material density.
    xiA1 to xiD4 : float
        The 16 lamination parameters used to define the laminate.

    Returns
    -------
    lam : Laminate
        laminate with the ABD and ABDE matrices already calculated

    """
    lam = Laminate()
    lam.h = thickness
    lam.matobj = read_laminaprop(laminaprop, rho)
    lam.xiA = np.array([1, xiA1, xiA2, xiA3, xiA4], dtype=np.float64)
    lam.xiB = np.array([0, xiB1, xiB2, xiB3, xiB4], dtype=np.float64)
    lam.xiD = np.array([1, xiD1, xiD2, xiD3, xiD4], dtype=np.float64)
    lam.xiE = np.array([1, xiE1, xiE2, xiE3, xiE4], dtype=np.float64)

    lam.calc_ABDE_from_lamination_parameters()
    return lam


class Laminate(object):
    r"""
    =========  ===========================================================
    attribute  description
    =========  ===========================================================
    plies      list of plies
    h          total thickness of the laminate
    offset     offset at the normal direction
    e1         equivalent laminate modulus in 1 direction
    e2         equivalent laminate modulus in 2 direction
    g12        equivalent laminate shear modulus in 12 direction
    nu12       equivalent laminate Poisson ratio in 12 direction
    nu21       equivalent laminate Poisson ratio in 21 direction
    xiA        laminate parameters for extensional matrix A
    xiB        laminate parameters for extension-bending matrix B
    xiD        laminate parameters for bending matrix D
    A          laminate extension matrix
    B          laminate extension-bending matrix
    D          laminate bending matrix
    E          laminate transferse shear matrix
    ABD        laminate ABD matrix
    ABDE       laminate ABD matrix with transverse shear terms
    scf_k13    shear correction factor 13
    scf_k23    shear correction factor 23
    =========  ===========================================================

    """
    def __init__(self):
        self.plies = []
        self.matobj = None
        self.h = None
        self.offset = 0.
        self.e1 = None
        self.e2 = None
        self.e3 = None
        self.nu12 = None
        self.g12 = None
        self.g13 = None
        self.g23 = None
        self.xiA = None
        self.xiB = None
        self.xiD = None
        self.xiE = None
        self.A = None
        self.B = None
        self.D = None
        self.E = None
        self.ABD = None
        self.ABDE = None
        self.scf_k13 = 5/6.
        self.scf_k23 = 5/6.


    def rebuild(self):
        lam_thick = 0.
        rho_avg = 0.
        for ply in self.plies:
            ply.rebuild()
            lam_thick += ply.h
            rho = ply.matobj.rho if ply.matobj.rho is not None else 0.
            rho_avg += rho * ply.h
        self.h = lam_thick
        self.rho = rho_avg / lam_thick


    def calc_scf(self):
        """Calculate improved shear correction factors

        Reference:

            Vlachoutsis, S. "Shear correction factors for plates and shells",
            Int. Journal for Numerical Methods in Engineering, Vol. 33,
            1537-1552, 1992.

            http://onlinelibrary.wiley.com/doi/10.1002/nme.1620330712/full


        Using "one shear correction factor" (see reference), assuming:

        - constant G13, G23, E1, E2, nu12, nu21 within each ply
        - g1 calculated using z at the middle of each ply
        - zn1 = Laminate.offset

        Returns
        -------

        k13, k23 : tuple
            Shear correction factors. Also updates attributes: `scf_k13`
            and `scf_k23`.

        """
        D1 = 0
        R1 = 0
        den1 = 0

        D2 = 0
        R2 = 0
        den2 = 0

        offset = self.offset
        zbot = -self.h/2 + offset
        z1 = zbot

        for ply in self.plies:
            z2 = z1 + ply.h
            e1 = (ply.matobj.e1 * np.cos(np.deg2rad(ply.theta)) +
                  ply.matobj.e2 * np.sin(np.deg2rad(ply.theta)))
            e2 = (ply.matobj.e2 * np.cos(np.deg2rad(ply.theta)) +
                  ply.matobj.e1 * np.sin(np.deg2rad(ply.theta)))
            nu12 = (ply.matobj.nu12 * np.cos(np.deg2rad(ply.theta)) +
                  ply.matobj.nu21 * np.sin(np.deg2rad(ply.theta)))
            nu21 = (ply.matobj.nu21 * np.cos(np.deg2rad(ply.theta)) +
              ply.matobj.nu12 * np.sin(np.deg2rad(ply.theta)))

            D1 += e1 / (1 - nu12*nu21)
            R1 += D1*((z2 - offset)**3/3 - (z1 - offset)**3/3)
            g13 = ply.matobj.g13
            d1 = g13 * ply.h
            den1 += d1 * (self.h / ply.h) * D1**2*(15*offset*z1**4 + 30*offset*z1**2*zbot*(2*offset - zbot) - 15*offset*z2**4 + 30*offset*z2**2*zbot*(-2*offset + zbot) - 3*z1**5 + 10*z1**3*(-2*offset**2 - 2*offset*zbot + zbot**2) - 15*z1*zbot**2*(4*offset**2 - 4*offset*zbot + zbot**2) + 3*z2**5 + 10*z2**3*(2*offset**2 + 2*offset*zbot - zbot**2) + 15*z2*zbot**2*(4*offset**2 - 4*offset*zbot + zbot**2))/(60*g13)

            D2 += e2 / (1 - nu12*nu21)
            R2 += D2*((z2 - self.offset)**3/3 - (z1 - self.offset)**3/3)
            g23 = ply.matobj.g23
            d2 = g23 * ply.h
            den2 += d2 * (self.h / ply.h) * D2**2*(15*offset*z1**4 + 30*offset*z1**2*zbot*(2*offset - zbot) - 15*offset*z2**4 + 30*offset*z2**2*zbot*(-2*offset + zbot) - 3*z1**5 + 10*z1**3*(-2*offset**2 - 2*offset*zbot + zbot**2) - 15*z1*zbot**2*(4*offset**2 - 4*offset*zbot + zbot**2) + 3*z2**5 + 10*z2**3*(2*offset**2 + 2*offset*zbot - zbot**2) + 15*z2*zbot**2*(4*offset**2 - 4*offset*zbot + zbot**2))/(60*g23)

            z1 = z2

        self.scf_k13 = R1**2 / den1
        self.scf_k23 = R2**2 / den2

        return self.scf_k13, self.scf_k23


    def calc_equivalent_modulus(self):
        """Calculates the equivalent laminate properties.

        The following attributes are calculated:
            e1, e2, g12, nu12, nu21

        """
        AI = np.matrix(self.ABD, dtype=np.float64).I
        a11, a12, a22, a33 = AI[0,0], AI[0,1], AI[1,1], AI[2,2]
        self.e1 = 1./(self.h*a11)
        self.e2 = 1./(self.h*a22)
        self.g12 = 1./(self.h*a33)
        self.nu12 = - a12 / a11
        self.nu21 = - a12 / a22


    def calc_lamination_parameters(self):
        """Calculate the lamination parameters.

        The following attributes are calculated:
            xiA, xiB, xiD, xiE

        """
        if len(self.plies) == 0:
            if self.xiA is None:
                raise ValueError('Laminate with 0 plies!')
            else:
                return
        xiA1, xiA2, xiA3, xiA4 = 0, 0, 0, 0
        xiB1, xiB2, xiB3, xiB4 = 0, 0, 0, 0
        xiD1, xiD2, xiD3, xiD4 = 0, 0, 0, 0
        xiE1, xiE2, xiE3, xiE4 = 0, 0, 0, 0

        lam_thick = sum([ply.h for ply in self.plies])
        self.h = lam_thick

        h0 = -lam_thick/2. + self.offset
        for ply in self.plies:
            if self.matobj is None:
                self.matobj = ply.matobj
            else:
                assert np.allclose(self.matobj.u, ply.matobj.u), "Plies with different materials"
            hk_1 = h0
            h0 += ply.h
            hk = h0

            Afac = ply.h / lam_thick
            Bfac = (2. / lam_thick**2) * (hk**2 - hk_1**2)
            Dfac = (4. / lam_thick**3) * (hk**3 - hk_1**3)
            Efac = (1. / lam_thick) * (hk - hk_1)

            thetarad = np.deg2rad(ply.theta)
            cos2t = np.cos(2*thetarad)
            sin2t = np.sin(2*thetarad)
            cos4t = np.cos(4*thetarad)
            sin4t = np.sin(4*thetarad)

            xiA1 += Afac * cos2t
            xiA2 += Afac * sin2t
            xiA3 += Afac * cos4t
            xiA4 += Afac * sin4t

            xiB1 += Bfac * cos2t
            xiB2 += Bfac * sin2t
            xiB3 += Bfac * cos4t
            xiB4 += Bfac * sin4t

            xiD1 += Dfac * cos2t
            xiD2 += Dfac * sin2t
            xiD3 += Dfac * cos4t
            xiD4 += Dfac * sin4t

            xiE1 += Efac * cos2t
            xiE2 += Efac * sin2t
            xiE3 += Efac * cos4t
            xiE4 += Efac * sin4t

        self.xiA = np.array([1, xiA1, xiA2, xiA3, xiA4], dtype=np.float64)
        self.xiB = np.array([0, xiB1, xiB2, xiB3, xiB4], dtype=np.float64)
        self.xiD = np.array([1, xiD1, xiD2, xiD3, xiD4], dtype=np.float64)
        self.xiE = np.array([1, xiE1, xiE2, xiE3, xiE4], dtype=np.float64)


    def calc_ABDE_from_lamination_parameters(self):
        """Use the ABDE matrix based on lamination parameters.

        Given the lamination parameters ``xiA``, ``xiB``, ``xiC`` and ``xiD``,
        the ABD matrix is calculated.

        """
        # dummies used to unpack vector results
        du1, du2, du3, du4, du5, du6 = 0, 0, 0, 0, 0, 0
        # A matrix terms
        A11,A22,A12, du1,du2,du3, A66,A16,A26 =\
            (self.h       ) * np.dot(self.matobj.u, self.xiA)
        # B matrix terms
        B11,B22,B12, du1,du2,du3, B66,B16,B26 =\
            (self.h**2/4. ) * np.dot(self.matobj.u, self.xiB)
        # D matrix terms
        D11,D22,D12, du1,du2,du3, D66,D16,D26 =\
            (self.h**3/12.) * np.dot(self.matobj.u, self.xiD)
        # E matrix terms
        du1,du2,du3, E44,E55,E45, du4,du5,du6 =\
            (self.h       ) * np.dot(self.matobj.u, self.xiE)

        self.A = np.array([[A11, A12, A16],
                           [A12, A22, A26],
                           [A16, A26, A66]], dtype=np.float64)

        self.B = np.array([[B11, B12, B16],
                           [B12, B22, B26],
                           [B16, B26, B66]], dtype=np.float64)

        self.D = np.array([[D11, D12, D16],
                           [D12, D22, D26],
                           [D16, D26, D66]], dtype=np.float64)

        # printing E acoordingly to Reddy definition for E44, E45 and E55
        self.E = np.array([[E55, E45],
                           [E45, E44]], dtype=np.float64)

        self.ABD = np.array([[A11, A12, A16, B11, B12, B16],
                             [A12, A22, A26, B12, B22, B26],
                             [A16, A26, A66, B16, B26, B66],
                             [B11, B12, B16, D11, D12, D16],
                             [B12, B22, B26, D12, D22, D26],
                             [B16, B26, B66, D16, D26, D66]], dtype=np.float64)

        # printing ABDE acoordingly to Reddy definition for E44, E45 and E55
        self.ABDE = np.array([[A11, A12, A16, B11, B12, B16, 0, 0],
                              [A12, A22, A26, B12, B22, B26, 0, 0],
                              [A16, A26, A66, B16, B26, B66, 0, 0],
                              [B11, B12, B16, D11, D12, D16, 0, 0],
                              [B12, B22, B26, D12, D22, D26, 0, 0],
                              [B16, B26, B66, D16, D26, D66, 0, 0],
                              [0, 0, 0, 0, 0, 0, E55, E45],
                              [0, 0, 0, 0, 0, 0, E45, E44]],
                               dtype=np.float64)


    def calc_constitutive_matrix(self):
        """Calculates the laminate constitutive matrix

        This is the commonly called ``ABD`` matrix with ``shape=(6, 6)`` when
        the classical laminated plate theory is used, or the ``ABDE`` matrix
        when the first-order shear deformation theory is used, containing the
        transverse shear terms.

        """
        self.A_general = np.zeros([5,5], dtype=np.float64)
        self.B_general = np.zeros([5,5], dtype=np.float64)
        self.D_general = np.zeros([5,5], dtype=np.float64)

        lam_thick = sum([ply.h for ply in self.plies])
        self.h = lam_thick

        h0 = -lam_thick/2 + self.offset
        for ply in self.plies:
            hk_1 = h0
            h0 += ply.h
            hk = h0
            self.A_general += ply.QL*(hk - hk_1)
            self.B_general += 1/2.*ply.QL*(hk**2 - hk_1**2)
            self.D_general += 1/3.*ply.QL*(hk**3 - hk_1**3)

        self.A = self.A_general[0:3, 0:3]
        self.B = self.B_general[0:3, 0:3]
        self.D = self.D_general[0:3, 0:3]
        self.E = self.A_general[3:5, 3:5]

        conc1 = np.concatenate([self.A, self.B], axis=1)
        conc2 = np.concatenate([self.B, self.D], axis=1)

        self.ABD = np.concatenate([conc1, conc2], axis=0)
        self.ABDE = np.zeros((8, 8), dtype=np.float64)
        self.ABDE[0:6, 0:6] = self.ABD
        self.ABDE[6:8, 6:8] = self.E


    def force_balanced_LP(self):
        r"""Force balanced lamination parameters

        The lamination parameters `\xi_{A2}` and `\xi_{A4}` are set to null
        to force a balanced laminate.

        """
        dummy, xiA1, xiA2, xiA3, xiA4 = self.xiA
        self.xiA = np.array([1, xiA1, 0, xiA3, 0], dtype=np.float64)
        self.calc_ABDE_from_lamination_parameters()


    def force_symmetric_LP(self):
        r"""Force symmetric lamination parameters

        The lamination parameters `\xi_{Bi}` are set to null
        to force a symmetric laminate.

        """
        self.xiB = np.zeros(5)
        self.calc_ABDE_from_lamination_parameters()


    def force_orthotropic(self):
        r"""Force an orthotropic laminate

        The terms
        `A_{13}`, `A_{23}`, `A_{31}`, `A_{32}`,
        `B_{13}`, `B_{23}`, `B_{31}`, `B_{32}`,
        `D_{13}`, `D_{23}`, `D_{31}`, `D_{32}` are set to zero to force an
        orthotropic laminate.

        """
        if self.offset != 0.:
            raise RuntimeError(
                    'Laminates with offset cannot be forced orthotropic!')
        self.A[0, 2] = 0.
        self.A[1, 2] = 0.
        self.A[2, 0] = 0.
        self.A[2, 1] = 0.

        self.B[0, 2] = 0.
        self.B[1, 2] = 0.
        self.B[2, 0] = 0.
        self.B[2, 1] = 0.

        self.D[0, 2] = 0.
        self.D[1, 2] = 0.
        self.D[2, 0] = 0.
        self.D[2, 1] = 0.

        self.ABD[0, 2] = 0. # A16
        self.ABD[1, 2] = 0. # A26
        self.ABD[2, 0] = 0. # A61
        self.ABD[2, 1] = 0. # A62

        self.ABD[0, 5] = 0. # B16
        self.ABD[5, 0] = 0. # B61
        self.ABD[1, 5] = 0. # B26
        self.ABD[5, 1] = 0. # B62

        self.ABD[3, 2] = 0. # B16
        self.ABD[2, 3] = 0. # B61
        self.ABD[4, 2] = 0. # B26
        self.ABD[2, 4] = 0. # B62

        self.ABD[3, 5] = 0. # D16
        self.ABD[4, 5] = 0. # D26
        self.ABD[5, 3] = 0. # D61
        self.ABD[5, 4] = 0. # D62

        self.ABDE[0, 2] = 0. # A16
        self.ABDE[1, 2] = 0. # A26
        self.ABDE[2, 0] = 0. # A61
        self.ABDE[2, 1] = 0. # A62

        self.ABDE[0, 5] = 0. # B16
        self.ABDE[5, 0] = 0. # B61
        self.ABDE[1, 5] = 0. # B26
        self.ABDE[5, 1] = 0. # B62

        self.ABDE[3, 2] = 0. # B16
        self.ABDE[2, 3] = 0. # B61
        self.ABDE[4, 2] = 0. # B26
        self.ABDE[2, 4] = 0. # B62

        self.ABDE[3, 5] = 0. # D16
        self.ABDE[4, 5] = 0. # D26
        self.ABDE[5, 3] = 0. # D61
        self.ABDE[5, 4] = 0. # D62


    def force_symmetric(self):
        """Force a symmetric laminate

        The `B` terms of the constitutive matrix are set to zero.

        """
        if self.offset != 0.:
            raise RuntimeError(
                    'Laminates with offset cannot be forced symmetric!')
        self.B = np.zeros((3,3))
        self.ABD[0:3, 3:6] = 0
        self.ABD[3:6, 0:3] = 0

        self.ABDE[0:3, 3:6] = 0
        self.ABDE[3:6, 0:3] = 0

