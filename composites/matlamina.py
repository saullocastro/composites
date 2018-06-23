"""
Composite Matlamina Module (:mod:`composites.matlamina`)
================================================================

.. currentmodule:: composites.matlamina

"""
import numpy as np

from .logger import error


class MatLamina(object):
    r"""
    Orthotropic material lamina

    ==========  ==========================================================
    attributes  description
    ==========  ==========================================================
    e1          Young Modulus in direction 1
    e2          Young Modulus in direction 2
    g12         in-plane shear modulus
    g13         transverse shear modulus for plane 1-Z
    g23         transverse shear modulus for plane 2-Z
    nu12        Poisson's ratio 12
    nu13        Poisson's ratio 13
    nu23        Poisson's ratio 23
    nu21        Poisson's ratio 21: use formula nu12/e1 = nu21/e2
    nu31        Poisson's ratio 31: use formula nu31/e3 = nu13/e1
    nu32        Poisson's ratio 32: use formula nu23/e2 = nu32/e3
    rho         especific mass (mass / volume)
    a1          thermal expansion coeffiecient in direction 1
    a2          thermal expansion coeffiecient in direction 2
    a3          thermal expansion coeffiecient in direction 3
    tref        reference temperature
    st1,st2     allowable tensile stresses for directions 1 and 2
    sc1,sc2     allowable compressive stresses for directions 1 and 2
    ss12        allowable in-plane stress for shear
    strn        allowable strain for direction 1
    q11         lamina constitutive constant 11
    q12         lamina constitutive constant 12
    q13         lamina constitutive constant 13
    q21         lamina constitutive constant 21
    q22         lamina constitutive constant 22
    q23         lamina constitutive constant 23
    q31         lamina constitutive constant 31
    q32         lamina constitutive constant 32
    q33         lamina constitutive constant 33
    q44         lamina constitutive constant 44
    q55         lamina constitutive constant 55
    q66         lamina constitutive constant 66
    u           matrix with lamina invariants
    c           matrix with lamina stiffness coefficients
    ==========  ==========================================================

    Notes
    -----
    For isotropic materials when the user defines `\nu` and `E`, `G` will be
    recaculated based on equation: `G = E/(2 \times (1+\nu))`; in a lower
    priority if the user defines `\nu` and `G`, `E` will be recaculated based
    on equation: `E = 2 \times (1+\nu) \times G`.

    """
    def __init__(self):
        super(MatLamina, self).__init__()
        self.e1 = None
        self.e2 = None
        self.e3 = 0
        self.g12 = None
        self.g13 = None
        self.g23 = None
        self.nu12 = None
        self.nu13 = 0
        self.nu21 = None
        self.nu23 = 0
        self.nu31 = 0
        self.nu32 = 0
        self.rho = None
        self.a1 = None
        self.a2 = None
        self.a3 = None
        self.st1 = None
        self.st2 = None
        self.sc1 = None
        self.sc2 = None
        self.ss12 = None
        self.strn = None
        self.q11 = None
        self.q12 = None
        self.q13 = None
        self.q21 = None
        self.q22 = None
        self.q23 = None
        self.q31 = None
        self.q32 = None
        self.q33 = None
        self.q44 = None
        self.q55 = None
        self.q66 = None
        self.u = None


    def rebuild(self):
        #
        # from references:
        #   Reddy, J. N., Mechanics of laminated composite plates and shells.
        #   Theory and analysis. Second Edition. CRC Press, 2004.
        e1 = self.e1
        e2 = self.e2
        e3 = self.e3
        nu12 = self.nu12
        nu21 = self.nu21
        nu13 = self.nu13
        nu31 = self.nu31
        nu23 = self.nu23
        nu32 = self.nu32
        delta = (1-nu12*nu21-nu23*nu32-nu31*nu13-2*nu21*nu32*nu13)/(e1*e2)
        c11 = (1    - nu23*nu23)/(delta*e2)
        c12 = (nu21 + nu31*nu23)/(delta*e2)
        c13 = (nu31 + nu21*nu32)/(delta*e2)
        c22 = (1    - nu13*nu31)/(delta*e1)
        c23 = (nu32 + nu12*nu31)/(delta*e1)
        c33 = e3*(1    - nu12*nu21)/(delta*e1*e2)
        c44 = self.g23
        c55 = self.g13
        c66 = self.g12
        self.c = np.array(
            [[c11, c12, c13,   0,   0,   0],
             [c12, c22, c23,   0,   0,   0],
             [c13, c23, c33,   0,   0,   0],
             [  0,   0,   0, c44,   0,   0],
             [  0,   0,   0,   0, c55,   0],
             [  0,   0,   0,   0,   0, c66]], dtype=np.float64)

        #
        # from references:
        #   hansen_hvejsen_2007 page 43
        #
        #   Guerdal Z., R. T. Haftka and P. Hajela (1999), Design and
        #   Optimization of Laminated Composite Materials, Wiley-Interscience.
        den = (1 - self.nu12 * self.nu21
                 - self.nu13 * self.nu31
                 - self.nu23 * self.nu32
                 - self.nu12 * self.nu23 * self.nu31
                 - self.nu13 * self.nu21 * self.nu32)
        den = np.array(den, dtype=np.float64)
        self.q11 = self.e1*(1         - self.nu23 * self.nu32) / den
        self.q12 = self.e1*(self.nu21 + self.nu23 * self.nu31) / den
        self.q13 = self.e1*(self.nu31 + self.nu21 * self.nu32) / den
        self.q21 = self.e2*(self.nu12 + self.nu13 * self.nu32) / den
        self.q22 = self.e2*(1         - self.nu13 * self.nu31) / den
        self.q23 = self.e2*(self.nu32 + self.nu12 * self.nu31) / den
        self.q31 = self.e3*(self.nu13 + self.nu12 * self.nu32) / den
        self.q32 = self.e3*(self.nu23 + self.nu13 * self.nu21) / den
        self.q33 = self.e3*(1         - self.nu12 * self.nu21) / den
        self.q44 = self.g12
        self.q55 = self.g23
        self.q66 = self.g13
        #
        # from reference:
        #   Jones R. M. (1999), Mechanics of Composite Materials, second edn,
        #   Taylor & Francis, Inc., 325 Chestnut Street, Philadelphia,
        #   PA 19106. ISBN 1-56032-712-X
        # slightly changed to include the transverse shear terms u6 ans
        # u7, taken from ABAQUS Example Problems Manual, vol1,
        # example 1.2.2 Laminated composite shell: buckling of a
        # cylindrical panel with a circular hole
        #
        u1 = (3*self.q11 + 3*self.q22 + 2*self.q12 + 4*self.q44) / 8.
        u2 = (self.q11 - self.q22) / 2.
        u3 = (self.q11 + self.q22 - 2*self.q12 - 4*self.q44) / 8.
        u4 = (self.q11 + self.q22 + 6*self.q12 - 4*self.q44) / 8.
        u5 = (u1-u4) / 2.
        u6 = (self.q55 + self.q66) / 2.
        u7 = (self.q55 - self.q66) / 2.
        self.u = np.array(
            [[u1,  u2,    0,  u3,   0],                     # q11
             [u1, -u2,    0,  u3,   0],                     # q22
             [u4,   0,    0, -u3,   0],                     # q12
             [u6,  u7,    0,   0,   0],                     # q55
             [u6, -u7,    0,   0,   0],                     # q66
             [ 0,   0,  -u7,   0,   0],                     # q56
             [u5,   0,    0, -u3,   0],                     # q44
             [ 0,   0, u2/2,   0,  u3],                     # q14
             [ 0,   0, u2/2,   0, -u3]], dtype=np.float64)  # q24


    def read_inputs(self, inputs={}):
        if len(inputs) > 0:
            self = user_setattr(self, inputs)
        if not self.nu21:
            nu21 = np.array(self.nu12*self.e2/self.e1, dtype=np.float64)
            self.nu21 = nu21
        if not self.nu12:
            nu12 = np.array(self.nu21*self.e1/self.e2, dtype=np.float64)
            self.nu12 = nu12


def read_laminaprop(laminaprop=None, rho=None):
    """Returns a ``MatLamina`` object based on an input ``laminaprop`` tuple

    Parameters
    ----------
    laminaprop : list or tuple
        Tuple containing the folliwing entries:

            (e1, e2, nu12, g12, g13, g23, e3, nu13, nu23)

        for othotropic materials the user can only supply:

            (e1, e2, nu12, g12, g13, g23)

        for isotropic materials the user can only supply:

            (e1, e2, nu12)

        ======  ==============================
        symbol  value
        ======  ==============================
        e1      Young Module in direction 1
        e2      Young Module in direction 2
        nu12    12 Poisson's ratio
        g12     12 Shear Modulus
        g13     13 Shear Modulus
        g23     13 Shear Modulus
        e3      Young Module in direction 3
        nu13    13 Poisson's ratio
        nu23    23 Poisson's ratio
        ======  ==============================


    rho : float, optional
        Material density


    Returns
    -------
    matlam : MatLamina
        A :class:`.MatLamina` object.

    """
    matlam = MatLamina()

    #laminaProp = (e1, e2, nu12, g12, g13, g23, e3, nu13, nu23)
    if laminaprop == None:
        error('laminaprop must be a tuple in the following format:\n\t'
              +'(e1, e2, nu12, g12, g13, g23, e3, nu13, nu23)')
    if len(laminaprop) == 3: #ISOTROPIC
        e = laminaprop[0]
        nu = laminaprop[2]
        g = e/(2*(1+nu))
        laminaprop = (e, e, nu, g, g, g, e, nu, nu)
    nu12 = laminaprop[2]

    if len(laminaprop) < 9:
        e2 = laminaprop[1]
        laminaprop = tuple(list(laminaprop)[:6] + [e2, nu12, nu12])
    matlam.e1 = laminaprop[0]
    matlam.e2 = laminaprop[1]
    matlam.e3 = laminaprop[6]
    matlam.nu12 = laminaprop[2]
    matlam.nu13 = laminaprop[7]
    matlam.nu23 = laminaprop[8]
    matlam.nu21 = matlam.nu12 * matlam.e2 / matlam.e1
    matlam.nu31 = matlam.nu13 * matlam.e3 / matlam.e1
    matlam.nu32 = matlam.nu23 * matlam.e3 / matlam.e2
    matlam.g12 = laminaprop[3]
    matlam.g13 = laminaprop[4]
    matlam.g23 = laminaprop[5]
    matlam.rho = rho

    matlam.rebuild()

    return matlam

