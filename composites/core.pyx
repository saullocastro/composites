#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
#cython: nonecheck=False
#cython: infer_types=False
#cython: overflowcheck=False
"""
Composites Core Module (:mod:`composites.core`)
==============================================================

.. currentmodule:: composites.core

"""
import copyreg
import warnings

import numpy as np

DOUBLE = np.float64


cdef class LaminationParameters:
    r"""Lamination parameters

    Attributes
    ----------
    xiA1, xiA2, xiA3, xiA4 : float
        Lamination parameters `\xi_{Ai}` (in-plane)
    xiB1, xiB2, xiB3, xiB4 : float
        Lamination parameters `\xi_{Bi}` (in-plane coupling with bending)
    xiD1, xiD2, xiD3, xiD4 : float
        Lamination parameters `\xi_{Di}` (bending)
    xiAtrans1, xiAtrans2 : float
        Lamination parameters `\xi_{{A_{trans}}i}` (transverse shear)

    """
    def __init__(LaminationParameters self):
        self.xiA1=0; self.xiA2=0; self.xiA3=0; self.xiA4=0
        self.xiB1=0; self.xiB2=0; self.xiB3=0; self.xiB4=0
        self.xiD1=0; self.xiD2=0; self.xiD3=0; self.xiD4=0
        self.xiAtrans1=0; self.xiAtrans2=0


cdef class MatLamina:
    r"""
    Orthotropic material lamina

    Attributes
    ----------

    e1 : float
        Young Modulus in direction 1
    e2 : float
        Young Modulus in direction 2
    g12 : float
        in-plane shear modulus
    g13 : float
        transverse shear modulus for plane 1-Z
    g23 : float
        transverse shear modulus for plane 2-Z
    nu12 :
        Poisson's ratio 12
    nu13 :
        Poisson's ratio 13
    nu23 :
        Poisson's ratio 23
    nu21 :
        Poisson's ratio 21: use formula nu12/e1 = nu21/e2
    nu31 :
        Poisson's ratio 31: use formula nu31/e3 = nu13/e1
    nu32 :
        Poisson's ratio 32: use formula nu23/e2 = nu32/e3
    rho :
        especific mass (mass / volume)
    a1 :
        thermal expansion coeffiecient in direction 1
    a2 :
        thermal expansion coeffiecient in direction 2
    a3 :
        thermal expansion coeffiecient in direction 3
    tref :
        reference temperature
    st1,st2 :
        allowable tensile stresses for directions 1 and 2
    sc1,sc2 :
        allowable compressive stresses for directions 1 and 2
    ss12 :
        allowable in-plane stress for shear
    q11 :
        lamina constitutive constant 11
    q12 :
        lamina constitutive constant 12
    q13 :
        lamina constitutive constant 13
    q21 :
        lamina constitutive constant 21
    q22 :
        lamina constitutive constant 22
    q23 :
        lamina constitutive constant 23
    q31 :
        lamina constitutive constant 31
    q32 :
        lamina constitutive constant 32
    q33 :
        lamina constitutive constant 33
    q44 :
        lamina constitutive constant 44
    q55 :
        lamina constitutive constant 55
    q66 :
        lamina constitutive constant 66
    ci :
        lamina stiffness constants
    ui :
        lamina material invariants

    Notes
    -----
    For isotropic materials when the user defines `\nu` and `E`, `G` will be
    recaculated based on equation: `G = E/(2 \times (1+\nu))`; in a lower
    priority if the user defines `\nu` and `G`, `E` will be recaculated based
    on equation: `E = 2 \times (1+\nu) \times G`.

    """
    def __init__(MatLamina self):
        pass

    cpdef void rebuild(MatLamina self):
        r"""Update constitutive and invariant terms

        Reference:

            Reddy, J. N., Mechanics of laminated composite plates and shells.
            Theory and analysis. Second Edition. CRC Press, 2004.

        """
        cdef double e1, e2, e3, nu12, nu21, nu13, nu31, nu23, nu32, delta, den
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
        self.c11 = (1    - nu23*nu23)/(delta*e2)
        self.c12 = (nu21 + nu31*nu23)/(delta*e2)
        self.c13 = (nu31 + nu21*nu32)/(delta*e2)
        self.c22 = (1    - nu13*nu31)/(delta*e1)
        self.c23 = (nu32 + nu12*nu31)/(delta*e1)
        self.c33 = e3*(1    - nu12*nu21)/(delta*e1*e2)
        self.c44 = self.g23
        self.c55 = self.g13
        self.c66 = self.g12

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
        self.q11 = self.e1*(1         - self.nu23 * self.nu32) / den
        self.q12 = self.e1*(self.nu21 + self.nu23 * self.nu31) / den
        self.q13 = self.e1*(self.nu31 + self.nu21 * self.nu32) / den
        self.q21 = self.e2*(self.nu12 + self.nu13 * self.nu32) / den
        self.q22 = self.e2*(1         - self.nu13 * self.nu31) / den
        self.q23 = self.e2*(self.nu32 + self.nu12 * self.nu31) / den
        self.q31 = self.e3*(self.nu13 + self.nu12 * self.nu32) / den
        self.q32 = self.e3*(self.nu23 + self.nu13 * self.nu21) / den
        self.q33 = self.e3*(1         - self.nu12 * self.nu21) / den
        self.q66 = self.g12
        self.q44 = self.g23
        self.q55 = self.g13
        #
        # from reference:
        #   Jones R. M. (1999), Mechanics of Composite Materials, second edn,
        #   Taylor & Francis, Inc., 325 Chestnut Street, Philadelphia,
        #   PA 19106. ISBN 1-56032-712-X
        # slightly changed to include the transverse shear terms u6 and u7,
        #   taken from ABAQUS Example Problems Manual, vol1, example 1.2.2
        #   Laminated composite shell: buckling of a
        #   cylindrical panel with a circular hole
        #
        self.u1 = (3*self.q11 + 3*self.q22 + 2*self.q12 + 4*self.q66) / 8.
        self.u2 = (self.q11 - self.q22) / 2.
        self.u3 = (self.q11 + self.q22 - 2*self.q12 - 4*self.q66) / 8.
        self.u4 = (self.q11 + self.q22 + 6*self.q12 - 4*self.q66) / 8.
        self.u5 = (self.u1 - self.u4) / 2.
        self.u6 = (self.q44 + self.q55) / 2.
        self.u7 = (self.q44 - self.q55) / 2.

    cpdef void trace_normalize_plane_stress(MatLamina self):
        r"""Trace-normalize the lamina properties for plane stress

        Modify the original :class:`.MatLamina` object with a
        trace-normalization performed after calculating the trace according to
        Eq. 1 of reference:

            Melo, J. D. D., Bi, J., and Tsai, S. W., 2017, “A Novel
            Invariant-Based Design Approach to Carbon Fiber Reinforced
            Laminates,” Compos. Struct., 159, pp. 44–52.

        The trace calculated as `tr = Q_{11} + Q_{22} + 2Q_{66}`.  The
        universal in-plane stress stiffness components
        `Q_{11},Q_{12},Q_{22},Q_{44},Q_{55},Q_{66}` are divided by `tr`, and
        the invariants `U_1,U_2,U_3,U_4,U_5,U_6,U_7` are calculated with the
        normalized stiffnesses, such they also become trace-normalized
        invariants. These can be accessed using the ``u1,u2,u3,u4,u5,u6,u7``
        attributes.

        """
        cdef double tr
        tr = self.q11 + self.q22 + 2*self.q66
        self.q11 /= tr
        self.q12 /= tr
        self.q22 /= tr
        self.q44 /= tr
        self.q55 /= tr
        self.q66 /= tr
        self.u1 /= tr
        self.u2 /= tr
        self.u3 /= tr
        self.u4 /= tr
        self.u5 /= tr
        self.u6 /= tr
        self.u7 /= tr

    cpdef double [:, ::1] get_constitutive_matrix(MatLamina self):
        r"""Return the constitutive matrix
        """
        return np.array(
            [[self.c11, self.c12, self.c13,   0,   0,   0],
             [self.c12, self.c22, self.c23,   0,   0,   0],
             [self.c13, self.c23, self.c33,   0,   0,   0],
             [  0,   0,   0, self.c44,   0,   0],
             [  0,   0,   0,   0, self.c55,   0],
             [  0,   0,   0,   0,   0, self.c66]], dtype=DOUBLE)

    cpdef double [:, ::1] get_invariant_matrix(MatLamina self):
        r"""Return the invariant matrix
        """
        return np.array(
            [[self.u1,  self.u2,    0,  self.u3,   0],            # q11
             [self.u1, -self.u2,    0,  self.u3,   0],            # q22
             [self.u4,   0,    0, -self.u3,   0],                 # q12
             [self.u5,   0,    0, -self.u3,   0],                 # q66
             [ 0,   0, self.u2/2.,   0,  self.u3],                # q16
             [ 0,   0, self.u2/2.,   0, -self.u3],                # q26
             [self.u6,  self.u7,    0,   0,   0],                 # q44
             [ 0,   0, -self.u7,    0,   0],                      # q45
             [self.u6, -self.u7,    0,   0,   0]], dtype=DOUBLE)  # q55


cdef class Lamina:
    r"""
    Attributes
    ----------

    plyid : int
        Identificaiton of the composite lamina
    matlamina : :class:`.MatLamina` object
        A :class:`.MatLamina` object
    h : float
        Ply thickness
    thetadeg : float
        Ply angle in degrees

    """
    def __init__(Lamina self):
        pass

    cpdef void rebuild(Lamina self):
        r"""Update constitutive matrices

        Reference:

            Reddy, J. N., Mechanics of Laminated Composite Plates and
            Shells - Theory and Analysys. Second Edition. CRC PRESS, 2004.
        """
        cdef double thetarad, e1, e2, nu12, nu21, g12, g13, g23
        cdef double q11, q12, q22, q44, q55, q16, q26, q66
        cdef double cos2, cos3, cos4, sin2, sin3, sin4, sincos
        thetarad = deg2rad(self.thetadeg)
        self.cost = cos(thetarad)
        self.cos2t = cos(2*thetarad)
        self.cos4t = cos(4*thetarad)
        self.sint = sin(thetarad)
        self.sin2t = sin(2*thetarad)
        self.sin4t = sin(4*thetarad)
        cos2 = self.cost**2
        cos3 = self.cost**3
        cos4 = self.cost**4
        sin2 = self.sint**2
        sin3 = self.sint**3
        sin4 = self.sint**4
        sincos = self.sint*self.cost
        # STRAINS
        # different from stress due to:
        #     2*e12 = e6    2*e13 = e5    2*e23 = e4
        # to laminate
        # self.Rstrain = np.transpose(self.Tstress)
        # to lamina
        # self.Tstrain = np.transpose(self.Rstress)
        e1   = self.matlamina.e1
        e2   = self.matlamina.e2
        nu12 = self.matlamina.nu12
        nu21 = self.matlamina.nu21
        g12  = self.matlamina.g12
        g13  = self.matlamina.g13
        g23  = self.matlamina.g23

        # plane stress
        #TODO plane strain
        q11  = e1/(1-nu12*nu21)
        q12  = nu12*e2/(1-nu12*nu21)
        q22  = e2/(1-nu12*nu21)
        q44  = g23
        q55  = g13
        q16 = 0
        q26 = 0
        q66  = g12

        self.q11L = q11*cos4 + 2*(q12 + 2*q66)*sin2*cos2 + q22*sin4
        self.q12L = (q11 + q22 - 4*q66)*sin2*cos2 + q12*(sin4 + cos4)
        self.q22L = q11*sin4 + 2*(q12 + 2*q66)*sin2*cos2 + q22*cos4
        self.q16L = (q11 - q12 - 2*q66)*self.sint*cos3 + (q12 - q22 + 2*q66)*sin3*self.cost
        self.q26L = (q11 - q12 - 2*q66)*sin3*self.cost + (q12 - q22 + 2*q66)*self.sint*cos3
        self.q66L = (q11 + q22 - 2*q12 - 2*q66)*sin2*cos2 + q66*(sin4 + cos4)
        self.q44L = q44*cos2 + q55*sin2
        self.q45L = (q55 - q44)*sincos
        self.q55L = q55*cos2 + q44*sin2

        #TODO add the thermal coeficient terms when calculating the
        #     stresses... to take into account eventual thermal expansions or
        #     contractions

    cpdef double [:, ::1] get_transf_matrix_displ_to_laminate(Lamina self):
        r"""Return displacement transformation matrix from lamina to laminate"""
        return np.array([[ self.cost, self.sint, 0],
                         [-self.sint, self.cost, 0],
                         [   0,     0, 1]], dtype=DOUBLE)

    cpdef double [:, ::1] get_constitutive_matrix(Lamina self):
        r"""Return the constitutive matrix"""
        return np.array([[self.q11L, self.q12L, self.q16L,    0,    0],
                         [self.q12L, self.q22L, self.q26L,    0,    0],
                         [self.q16L, self.q26L, self.q66L,    0,    0],
                         [   0,    0,    0, self.q44L, self.q45L],
                         [   0,    0,    0, self.q45L, self.q55L]], dtype=DOUBLE)

    cpdef double [:, ::1] get_transf_matrix_stress_to_lamina(Lamina self):
        r"""Return stress transformation matrix from laminate to lamina"""
        cdef double cos2, sin2, sincos
        cos2 = self.cost**2
        sin2 = self.sint**2
        sincos = self.sint*self.cost
        return np.array(
            [[ cos2, sin2, 0, 0, 0, self.sin2t],
             [ sin2, cos2, 0, 0, 0, -self.sin2t],
             [ 0, 0, 1, 0, 0, 0],
             [ 0, 0, 0, self.cost, -self.sint, 0],
             [ 0, 0, 0, self.sint,  self.cost, 0],
             [-sincos, sincos, 0, 0, 0, cos2-sin2]], dtype=DOUBLE)

    cpdef double [:, ::1] get_transf_matrix_stress_to_laminate(Lamina self):
        r"""Return stress transformation matrix from lamina to laminate"""
        cdef double cos2, sin2, sincos
        cos2 = self.cost**2
        sin2 = self.sint**2
        sincos = self.sint*self.cost
        return np.array(
            [[ cos2, sin2, 0, 0,   0, -self.sin2t],
             [ sin2, cos2, 0, 0,   0, self.sin2t],
             [ 0, 0, 1, 0, 0, 0],
             [ 0, 0, 0,  self.cost, self.sint, 0],
             [ 0, 0, 0, -self.sint, self.cost, 0],
             [sincos, -sincos, 0, 0, 0, cos2-sin2]], dtype=DOUBLE)


def _singular_Cs_error(int k, Lamina ply):
    return ValueError('Ply %d (plyid=%d, thetadeg=%g) has a singular '
            'transverse shear constitutive matrix Cs = [[q44L, q45L], '
            '[q45L, q55L]] = [[%g, %g], [%g, %g]]; g13 and g23 must be '
            'positive' % (k, ply.plyid, ply.thetadeg, ply.q44L, ply.q45L,
                ply.q45L, ply.q55L))


# NOTE names of the "cdef public" attributes of Laminate saved when pickling,
#      they must follow core.pxd (checked in tests/test_pickle.py)
_LAMINATE_STATE = tuple(
    ['%s%s' % (m, ij) for m in 'ABDEFH' for ij in ('11', '12', '16', '22', '26', '66')]
    + ['%s%s' % (m, ij) for m in ('A', 'Abar', 'Abarbar', 'D', 'F') for ij in ('44', '45', '55')]
    + ['e1', 'e2', 'g12', 'nu12', 'nu21', 'scf_k13', 'scf_k23', 'h', 'offset',
       'intrho', 'intrhoz', 'intrhoz2', 'plies', 'stack', 'shear_correction'])


cdef class Laminate:
    r"""
    Attributes
    ----------

    plies : list
        List of plies
    stack : list
        List of angles for each ply
    h : float
        Total thickness of the laminate
    offset : float
        Offset at the normal direction
    e1, e2 : float
        Equivalent laminate moduli in directions 1 and 2
    g12 : float
        Equivalent laminate shear modulus in the 12 direction
    nu12, nu21 : float
        Equivalent laminate Poisson ratios in the 12 and 21 directions
    A44, A45, A55 : float
        Transverse shear stiffnesses of the first-order shear deformation
        theory (FSDT), **with the shear correction already applied**
        according to ``shear_correction``. They are ready to be used in
        the element and no shear correction factor should be applied to
        them downstream. See :meth:`.calc_transverse_shear_stiffness`.
    Abar44, Abar45, Abar55 : float
        Constant-strain (uncorrected) transverse shear stiffnesses
        `\bar{A}_{ts} = \sum_k C_s^{(k)} h_k`. These are the terms to be used
        together with ``Dtrans`` and ``Ftrans`` in the third-order shear
        deformation theory (TSDT), which needs no shear correction.
    Abarbar44, Abarbar45, Abarbar55 : float
        Constant-stress transverse shear stiffnesses `\bar{\bar{A}}_{ts} = h^2
        [\sum_k (C_s^{(k)})^{-1} h_k]^{-1}`, for comparison only. Equal to
        ``nan`` when a ply has a singular `C_s^{(k)}`.
    shear_correction : str or None
        Method used to obtain ``A44``, ``A45``, ``A55`` from the ply data, see
        :meth:`.calc_transverse_shear_stiffness`. Default is ``'rohwer'``.
    scf_k13, scf_k23 : float
        Reported shear correction ratios ``A55/Abar55`` and ``A44/Abar44``.
        They are informative only, the correction is already inside ``A44``,
        ``A45``, ``A55``.
    intrho : float
        Integral `\int_{-h/2+offset}^{+h/2+offset} \rho(z) dz`, used in
        equivalent single layer finite element mass matrices
    intrhoz : float
        Integral `\int_{-h/2+offset}^{+h/2+offset} \rho(z)z dz`, used in
        equivalent single layer finite element mass matrices
    intrhoz2 : float
        Integral `\int_{-h/2+offset}^{+h/2+offset} \rho(z)z^2 dz`, used in
        equivalent single layer finite element mass matrices

    """
    def __init__(Laminate self):
        self.h = 0.
        self.e1 = 0.
        self.e2 = 0.
        self.g12 = 0.
        self.nu12 = 0.
        self.nu21 = 0.
        self.offset = 0.
        self.scf_k13 = 5/6.
        self.scf_k23 = 5/6.
        self.intrho = 0.
        self.intrhoz = 0.
        self.intrhoz2 = 0.
        self.plies = []
        self.stack = []
        self.shear_correction = 'rohwer'
        self._ts_ready = False

    def __reduce__(Laminate self):
        # NOTE the transverse shear distribution cache (_ts_z, _ts_fcoef) is
        #      left out and recomputed on demand after unpickling
        state = {name: getattr(self, name) for name in _LAMINATE_STATE}
        if hasattr(self, '__dict__'):
            state['__dict__'] = self.__dict__
        return copyreg.__newobj__, (type(self), ), state

    def __setstate__(Laminate self, dict state):
        state = dict(state)
        inst_dict = state.pop('__dict__', None)
        for name, value in state.items():
            setattr(self, name, value)
        if inst_dict:
            self.__dict__.update(inst_dict)
        self._ts_ready = False

    cdef double [:, ::1] get_A(Laminate self):
        return np.array([[self.A11, self.A12, self.A16],
                         [self.A12, self.A22, self.A26],
                         [self.A16, self.A26, self.A66]], dtype=DOUBLE)
    cdef double [:, ::1] get_B(Laminate self):
        return np.array([[self.B11, self.B12, self.B16],
                         [self.B12, self.B22, self.B26],
                         [self.B16, self.B26, self.B66]], dtype=DOUBLE)
    cdef double [:, ::1] get_D(Laminate self):
        return np.array([[self.D11, self.D12, self.D16],
                         [self.D12, self.D22, self.D26],
                         [self.D16, self.D26, self.D66]], dtype=DOUBLE)
    cdef double [:, ::1] get_E(Laminate self):
        return np.array([[self.E11, self.E12, self.E16],
                         [self.E12, self.E22, self.E26],
                         [self.E16, self.E26, self.E66]], dtype=DOUBLE)
    cdef double [:, ::1] get_F(Laminate self):
        return np.array([[self.F11, self.F12, self.F16],
                         [self.F12, self.F22, self.F26],
                         [self.F16, self.F26, self.F66]], dtype=DOUBLE)
    cdef double [:, ::1] get_H(Laminate self):
        return np.array([[self.H11, self.H12, self.H16],
                         [self.H12, self.H22, self.H26],
                         [self.H16, self.H26, self.H66]], dtype=DOUBLE)
    cdef double [:, ::1] get_Ats(Laminate self):
        return np.array([[self.A44, self.A45],
                         [self.A45, self.A55]], dtype=DOUBLE)
    cdef double [:, ::1] get_Atrans(Laminate self):
        # NOTE kept for cimporting packages, same as get_Ats
        return self.get_Ats()
    cdef double [:, ::1] get_Abar_ts(Laminate self):
        return np.array([[self.Abar44, self.Abar45],
                         [self.Abar45, self.Abar55]], dtype=DOUBLE)
    cdef double [:, ::1] get_Abarbar_ts(Laminate self):
        return np.array([[self.Abarbar44, self.Abarbar45],
                         [self.Abarbar45, self.Abarbar55]], dtype=DOUBLE)
    cdef double [:, ::1] get_Dtrans(Laminate self):
        return np.array([[self.D44, self.D45],
                         [self.D45, self.D55]], dtype=DOUBLE)
    cdef double [:, ::1] get_Ftrans(Laminate self):
        return np.array([[self.F44, self.F45],
                         [self.F45, self.F55]], dtype=DOUBLE)
    cdef double [:, ::1] get_ABD(Laminate self):
        return np.array([[self.A11, self.A12, self.A16, self.B11, self.B12, self.B16],
                         [self.A12, self.A22, self.A26, self.B12, self.B22, self.B26],
                         [self.A16, self.A26, self.A66, self.B16, self.B26, self.B66],
                         [self.B11, self.B12, self.B16, self.D11, self.D12, self.D16],
                         [self.B12, self.B22, self.B26, self.D12, self.D22, self.D26],
                         [self.B16, self.B26, self.B66, self.D16, self.D26, self.D66]], dtype=DOUBLE)
    @property
    def A(self):
        return np.asarray(self.get_A())
    @property
    def B(self):
        return np.asarray(self.get_B())
    @property
    def D(self):
        return np.asarray(self.get_D())
    @property
    def E(self):
        return np.asarray(self.get_E())
    @property
    def F(self):
        return np.asarray(self.get_F())
    @property
    def H(self):
        return np.asarray(self.get_H())
    @property
    def Ats(self):
        r"""Transverse shear stiffness matrix ``[[A44, A45], [A45, A55]]``

        Index 4 corresponds to `yz` and index 5 to `xz`, such that `\{Q_y,
        Q_x\}^T = A_{ts} \{\gamma_{yz}, \gamma_{xz}\}^T`. The shear correction
        is already applied, see :meth:`.calc_transverse_shear_stiffness`.

        """
        return np.asarray(self.get_Ats())
    @property
    def Atrans(self):
        r"""Deprecated, use :attr:`.Ats` instead

        Returns the same corrected matrix as :attr:`.Ats`.

        """
        warnings.warn("'Laminate.Atrans' is deprecated, use 'Laminate.Ats' "
                      "instead, which contains the shear correction",
                      DeprecationWarning, stacklevel=2)
        return np.asarray(self.get_Ats())
    @property
    def Abar_ts(self):
        r"""Constant-strain ``[[Abar44, Abar45], [Abar45, Abar55]]``

        Uncorrected transverse shear stiffness, to be used with ``Dtrans`` and
        ``Ftrans`` in the third-order shear deformation theory (TSDT).

        """
        return np.asarray(self.get_Abar_ts())
    @property
    def Abarbar_ts(self):
        r"""Constant-stress ``[[Abarbar44, Abarbar45], [Abarbar45, Abarbar55]]``"""
        return np.asarray(self.get_Abarbar_ts())
    @property
    def Dtrans(self):
        return np.asarray(self.get_Dtrans())
    @property
    def Ftrans(self):
        return np.asarray(self.get_Ftrans())
    @property
    def ABD(self):
        return np.asarray(self.get_ABD())


    cpdef tuple calc_scf(Laminate self):
        r"""Recompute the transverse shear stiffness and return the ratios

        .. deprecated:: 0.9.1
            Use :meth:`.calc_transverse_shear_stiffness`, called automatically
            by :meth:`.calc_constitutive_matrix`, and read ``scf_k13`` and
            ``scf_k23``.

        The ratios are informative only, since the correction is already
        applied to ``A44``, ``A45``, ``A55`` by
        :meth:`.calc_transverse_shear_stiffness`. If ``shear_correction`` is
        ``None``, it is set to ``'rohwer'`` before recomputing.

        Returns
        -------
        scf_k13, scf_k23 : tuple of float
            The ratios ``A55/Abar55`` and ``A44/Abar44``, also stored in the
            attributes ``scf_k13`` and ``scf_k23``.

        """
        warnings.warn("'Laminate.calc_scf' is deprecated, the transverse "
                      "shear stiffness is computed by "
                      "'calc_transverse_shear_stiffness', called by "
                      "'calc_constitutive_matrix'", DeprecationWarning,
                      stacklevel=2)
        if self.shear_correction is None:
            self.shear_correction = 'rohwer'
        self.calc_transverse_shear_stiffness()
        return self.scf_k13, self.scf_k23


    cdef int _calc_transverse_shear_distribution(Laminate self) except -1:
        r"""Coefficients of the distribution matrix `f^{(k)}(z)` (Rohwer, 1988)

        For each ply `k`, stores the 2x2 matrices `F_0`, `F_1`, `F_2` in
        ``_ts_fcoef[k, 0:3]``, such that `f^{(k)}(z) = F_0 + z F_1 + z^2 F_2`
        and `\{\tau_{yz}, \tau_{xz}\}^T = f^{(k)}(z) \{Q_y, Q_x\}^T`. The ply
        interfaces are stored in ``_ts_z``.

        The ABD matrix is recomputed from the plies, such that the
        distribution is consistent with them even if the laminate stiffness
        attributes were modified afterwards, e.g. by :meth:`.make_symmetric`.

        """
        cdef int i, j, k, N
        cdef double h, za, zb, dz1, dz2, dz3, zk
        cdef double [:, ::1] ABD, Hs
        cdef double [::1] z, acc_x, acc_y, p1, q1, p2, q2
        cdef double [::1] p1_prev, q1_prev, p2_prev, q2_prev
        cdef double [:, :, :, ::1] fcoef
        cdef double [3] c1, c2, c3
        cdef Lamina ply

        self._ts_ready = False
        N = <int>len(self.plies)
        if N == 0:
            raise ValueError('Laminate with 0 plies!')

        h = 0.
        for ply in self.plies:
            h += ply.h
        z = np.zeros(N + 1, dtype=DOUBLE)
        z[0] = -h/2. + self.offset

        # ABD = [[A, B], [B, D]] from the plies
        ABD = np.zeros((6, 6), dtype=DOUBLE)
        for k in range(N):
            ply = self.plies[k]
            za = z[k]
            zb = za + ply.h
            z[k+1] = zb
            dz1 = zb - za
            dz2 = (zb*zb - za*za)/2.
            dz3 = (zb*zb*zb - za*za*za)/3.
            c1[0] = ply.q11L; c1[1] = ply.q12L; c1[2] = ply.q16L
            c2[0] = ply.q12L; c2[1] = ply.q22L; c2[2] = ply.q26L
            c3[0] = ply.q16L; c3[1] = ply.q26L; c3[2] = ply.q66L
            for j in range(3):
                ABD[0, j] += c1[j]*dz1
                ABD[1, j] += c2[j]*dz1
                ABD[2, j] += c3[j]*dz1
                ABD[0, j+3] += c1[j]*dz2
                ABD[1, j+3] += c2[j]*dz2
                ABD[2, j+3] += c3[j]*dz2
                ABD[3, j] += c1[j]*dz2
                ABD[4, j] += c2[j]*dz2
                ABD[5, j] += c3[j]*dz2
                ABD[3, j+3] += c1[j]*dz3
                ABD[4, j+3] += c2[j]*dz3
                ABD[5, j+3] += c3[j]*dz3

        # Hstar = inv(ABD), with a symmetric diagonal scaling such that the
        # conditioning check does not depend on the units
        ABDnp = np.asarray(ABD)
        diag = np.diag(ABDnp)
        if not np.all(diag > 0):
            raise ValueError('The ABD matrix of the laminate is singular '
                             '(non-positive diagonal term), the transverse '
                             'shear distribution cannot be computed')
        scale = 1./np.sqrt(diag)
        M = scale[:, None]*ABDnp*scale[None, :]
        if not np.all(np.isfinite(M)) or np.linalg.cond(M) > 1e12:
            raise ValueError('The ABD matrix of the laminate is singular or '
                             'ill-conditioned, the transverse shear '
                             'distribution cannot be computed')
        Hs = np.ascontiguousarray(scale[:, None]*np.linalg.inv(M)*scale[None, :])

        fcoef = np.zeros((N, 3, 2, 2), dtype=DOUBLE)
        acc_x = np.zeros(6, dtype=DOUBLE)
        acc_y = np.zeros(6, dtype=DOUBLE)
        p1 = np.zeros(6, dtype=DOUBLE)
        q1 = np.zeros(6, dtype=DOUBLE)
        p2 = np.zeros(6, dtype=DOUBLE)
        q2 = np.zeros(6, dtype=DOUBLE)
        p1_prev = np.zeros(6, dtype=DOUBLE)
        q1_prev = np.zeros(6, dtype=DOUBLE)
        p2_prev = np.zeros(6, dtype=DOUBLE)
        q2_prev = np.zeros(6, dtype=DOUBLE)

        for k in range(N):
            ply = self.plies[k]
            c1[0] = ply.q11L; c1[1] = ply.q12L; c1[2] = ply.q16L
            c2[0] = ply.q12L; c2[1] = ply.q22L; c2[2] = ply.q26L
            zk = z[k]
            for j in range(6):
                # c Pi(z) = z c Hstar[:3, :] + z^2/2 c Hstar[3:, :]
                p1[j] = c1[0]*Hs[0, j] + c1[1]*Hs[1, j] + c1[2]*Hs[2, j]
                q1[j] = c1[0]*Hs[3, j] + c1[1]*Hs[4, j] + c1[2]*Hs[5, j]
                p2[j] = c2[0]*Hs[0, j] + c2[1]*Hs[1, j] + c2[2]*Hs[2, j]
                q2[j] = c2[0]*Hs[3, j] + c2[1]*Hs[4, j] + c2[2]*Hs[5, j]
                # a^(k) = sum_{i=1..k} (c^(i) - c^(i-1)) Pi(z_i)
                acc_x[j] += zk*(p1[j] - p1_prev[j]) + zk*zk/2.*(q1[j] - q1_prev[j])
                acc_y[j] += zk*(p2[j] - p2_prev[j]) + zk*zk/2.*(q2[j] - q2_prev[j])
                p1_prev[j] = p1[j]
                q1_prev[j] = q1[j]
                p2_prev[j] = p2[j]
                q2_prev[j] = q2[j]
            # tau_yz row: (a_y - c_2 Pi(z)) Ly, picking components 4 (Q_y) and 5 (Q_x)
            fcoef[k, 0, 0, 0] = acc_y[4]
            fcoef[k, 1, 0, 0] = -p2[4]
            fcoef[k, 2, 0, 0] = -q2[4]/2.
            fcoef[k, 0, 0, 1] = acc_y[5]
            fcoef[k, 1, 0, 1] = -p2[5]
            fcoef[k, 2, 0, 1] = -q2[5]/2.
            # tau_xz row: (a_x - c_1 Pi(z)) Lx, picking components 5 (Q_y) and 3 (Q_x)
            fcoef[k, 0, 1, 0] = acc_x[5]
            fcoef[k, 1, 1, 0] = -p1[5]
            fcoef[k, 2, 1, 0] = -q1[5]/2.
            fcoef[k, 0, 1, 1] = acc_x[3]
            fcoef[k, 1, 1, 1] = -p1[3]
            fcoef[k, 2, 1, 1] = -q1[3]/2.

        self._ts_z = z
        self._ts_fcoef = fcoef
        self._ts_ready = True
        return 0


    cpdef void calc_transverse_shear_stiffness(Laminate self) except *:
        r"""Update the transverse shear stiffnesses ``A44``, ``A45``, ``A55``

        Called at the end of :meth:`.calc_constitutive_matrix`, and computed
        once per laminate. The attributes ``A44``, ``A45``, ``A55`` are the
        transverse shear stiffnesses of the first-order shear deformation
        theory (FSDT) **with the shear correction already applied**, ready to
        be used in the element. No shear correction factor should be applied
        to them downstream.

        Conventions: `\{\tau_{yz}, \tau_{xz}\}^T`, `\{Q_y, Q_x\}^T`,
        `A_{ts} = [[A_{44}, A_{45}], [A_{45}, A_{55}]]`, index 4 corresponds
        to `yz` and index 5 to `xz`. The coordinate `z` is measured from the
        reference surface, with the plies running from `z_1 = -h/2 +
        offset` to `z_{N+1} = +h/2 + offset`, and `C_s^{(k)} = [[q_{44L},
        q_{45L}], [q_{45L}, q_{55L}]]`, which is in general a full matrix.

        The method is selected by the attribute ``shear_correction``:

        - ``'rohwer'`` (default): equilibrium approach of Rohwer (1988). The
          transverse shear stresses are obtained from the equilibrium of two
          cylindrical bending states, with zero tractions at the bottom and
          top faces and continuity at every interface, `\{\tau_{yz},
          \tau_{xz}\}^T = f^{(k)}(z) \{Q_y, Q_x\}^T`. The 2x2 stiffness is
          obtained from the complementary energy:

          .. math::

              A_{ts} = \left[ \sum_k \int_{z_k}^{z_{k+1}} f^{(k)T}
              \left(C_s^{(k)}\right)^{-1} f^{(k)} dz \right]^{-1}

          The integrand is a polynomial of degree 4 in `z` within each ply,
          such that the 3-point Gauss-Legendre rule used per ply is exact.
          The method is valid for arbitrary anisotropic and unsymmetric
          laminates, and the result does not depend on ``offset``. The result
          is not invariant to a rotation of the element axes, because the two
          cylindrical bending states are tied to the `x` and `y` axes. This is
          inherent to the method: for a `[0, 90]_s` CFRP laminate, deviations
          up to about 7 % are observed at 45 degrees.

        - ``'vlachoutsis'``: the scalar factors `k_{13}`, `k_{23}` of
          Vlachoutsis (1992) are applied to the constant-strain stiffness,
          ``A55 = k13*Abar55``, ``A44 = k23*Abar44``, and the ad-hoc ``A45 =
          (k13 + k23)/2*Abar45``, which cannot represent the coupling of
          angle-ply laminates with ``Abar45 = 0``. The factors use `C_{11}`
          and `C_{22}` of each ply in laminate axes and the direction-wise
          neutral surfaces, being exact only for specially orthotropic plies.

        - ``'constant'``: `k = 5/6`, i.e. ``A44 = 5/6*Abar44``, ``A45 =
          5/6*Abar45``, ``A55 = 5/6*Abar55``.

        - ``None``: no correction, ``A44 = Abar44``, ``A45 = Abar45``, ``A55 =
          Abar55``.

        The following attributes are also updated: ``Abar44``, ``Abar45``,
        ``Abar55`` (constant strain), ``Abarbar44``, ``Abarbar45``,
        ``Abarbar55`` (constant stress, ``nan`` if a ply has a singular
        `C_s^{(k)}`), and the ratios ``scf_k13 = A55/Abar55`` and ``scf_k23 =
        A44/Abar44``, which are informative only. For ``'rohwer'`` the
        through-thickness distribution used by
        :meth:`.calc_transverse_shear_stress` is also stored.

        References:

            Rohwer, K. "Improved transverse shear stiffness for layered
            finite elements", DFVLR-FB 88-32, 1988.

            Vlachoutsis, S. "Shear correction factors for plates and shells",
            Int. Journal for Numerical Methods in Engineering, Vol. 33,
            1537-1552, 1992.

        Raises
        ------
        ValueError
            If ``shear_correction`` is not recognized; if a ply has a
            singular transverse shear constitutive matrix `C_s^{(k)}`, e.g.
            ``g13 = 0`` or ``g23 = 0`` (``'rohwer'`` and ``'vlachoutsis'``);
            or if the ABD matrix of the laminate is singular (``'rohwer'``).

        """
        cdef int k, ig, N, alpha, singular_ply
        cdef double h, det, za, zb, zm, dz, zg, zg2, wdz
        cdef double i44, i45, i55
        cdef double f00, f01, f10, f11, g00, g01, g10, g11
        cdef double S00, S01, S11, detS
        cdef double Sbb44, Sbb45, Sbb55
        cdef double Dk, Gk, num, den, zn, R, d, I, gz, gza, kappa
        cdef double k13, k23
        cdef double [:, :, :, ::1] fc
        cdef double [::1] z
        cdef double [3] gp, gw
        cdef Lamina ply

        # 3-point Gauss-Legendre, exact for polynomials up to degree 5
        gp[0] = -0.7745966692414834; gp[1] = 0.; gp[2] = 0.7745966692414834
        gw[0] = 5./9.; gw[1] = 8./9.; gw[2] = 5./9.

        mode = self.shear_correction
        if not (mode is None or mode in ('rohwer', 'vlachoutsis', 'constant')):
            raise ValueError("shear_correction must be 'rohwer', "
                             "'vlachoutsis', 'constant' or None, got %r"
                             % (mode,))

        self._ts_ready = False
        self.A44 = 0; self.A45 = 0; self.A55 = 0
        self.Abar44 = 0; self.Abar45 = 0; self.Abar55 = 0
        self.Abarbar44 = 0; self.Abarbar45 = 0; self.Abarbar55 = 0
        N = <int>len(self.plies)
        if N == 0:
            return

        # constant-strain and constant-stress stiffnesses
        h = 0.
        singular_ply = -1
        Sbb44 = 0; Sbb45 = 0; Sbb55 = 0
        z = np.zeros(N + 1, dtype=DOUBLE)
        for k in range(N):
            ply = self.plies[k]
            h += ply.h
            z[k+1] = h
            self.Abar44 += ply.q44L*ply.h
            self.Abar45 += ply.q45L*ply.h
            self.Abar55 += ply.q55L*ply.h
            det = ply.q44L*ply.q55L - ply.q45L*ply.q45L
            if (ply.q44L <= 0 or ply.q55L <= 0
                    or det <= 1e-12*ply.q44L*ply.q55L):
                if singular_ply < 0:
                    singular_ply = k
            else:
                Sbb44 += ply.q55L/det*ply.h
                Sbb45 += -ply.q45L/det*ply.h
                Sbb55 += ply.q44L/det*ply.h
        for k in range(N + 1):
            z[k] += -h/2. + self.offset
        if singular_ply < 0:
            detS = Sbb44*Sbb55 - Sbb45*Sbb45
            self.Abarbar44 = h*h*Sbb55/detS
            self.Abarbar45 = -h*h*Sbb45/detS
            self.Abarbar55 = h*h*Sbb44/detS
        else:
            self.Abarbar44 = np.nan
            self.Abarbar45 = np.nan
            self.Abarbar55 = np.nan

        if mode is None:
            self.A44 = self.Abar44
            self.A45 = self.Abar45
            self.A55 = self.Abar55

        elif mode == 'constant':
            self.A44 = 5/6.*self.Abar44
            self.A45 = 5/6.*self.Abar45
            self.A55 = 5/6.*self.Abar55

        elif mode == 'rohwer':
            if singular_ply >= 0:
                raise _singular_Cs_error(singular_ply, self.plies[singular_ply])
            self._calc_transverse_shear_distribution()
            fc = self._ts_fcoef
            z = self._ts_z
            S00 = 0; S01 = 0; S11 = 0
            for k in range(N):
                ply = self.plies[k]
                det = ply.q44L*ply.q55L - ply.q45L*ply.q45L
                i44 = ply.q55L/det
                i45 = -ply.q45L/det
                i55 = ply.q44L/det
                za = z[k]
                zb = z[k+1]
                zm = (za + zb)/2.
                dz = (zb - za)/2.
                for ig in range(3):
                    zg = zm + dz*gp[ig]
                    zg2 = zg*zg
                    wdz = gw[ig]*dz
                    f00 = fc[k, 0, 0, 0] + zg*fc[k, 1, 0, 0] + zg2*fc[k, 2, 0, 0]
                    f01 = fc[k, 0, 0, 1] + zg*fc[k, 1, 0, 1] + zg2*fc[k, 2, 0, 1]
                    f10 = fc[k, 0, 1, 0] + zg*fc[k, 1, 1, 0] + zg2*fc[k, 2, 1, 0]
                    f11 = fc[k, 0, 1, 1] + zg*fc[k, 1, 1, 1] + zg2*fc[k, 2, 1, 1]
                    # g = inv(Cs) f
                    g00 = i44*f00 + i45*f10
                    g01 = i44*f01 + i45*f11
                    g10 = i45*f00 + i55*f10
                    g11 = i45*f01 + i55*f11
                    # S += f^T inv(Cs) f
                    S00 += wdz*(f00*g00 + f10*g10)
                    S01 += wdz*(f00*g01 + f10*g11)
                    S11 += wdz*(f01*g01 + f11*g11)
            detS = S00*S11 - S01*S01
            if not detS > 0:
                raise ValueError('Singular transverse shear compliance, the '
                                 'transverse shear stiffness cannot be '
                                 'computed')
            self.A44 = S11/detS
            self.A45 = -S01/detS
            self.A55 = S00/detS

        elif mode == 'vlachoutsis':
            k13 = 0; k23 = 0
            for alpha in range(2):
                # alpha = 0: direction 1 (x), uses C11 and G13 = q55L
                # alpha = 1: direction 2 (y), uses C22 and G23 = q44L
                num = 0; den = 0; d = 0
                for k in range(N):
                    ply = self.plies[k]
                    Gk = ply.q55L if alpha == 0 else ply.q44L
                    if not Gk > 0:
                        raise _singular_Cs_error(k, ply)
                    Dk = ply.q11L if alpha == 0 else ply.q22L
                    za = z[k]
                    zb = z[k+1]
                    num += Dk*(zb*zb - za*za)/2.
                    den += Dk*(zb - za)
                    d += Gk*ply.h
                if not den > 0:
                    raise ValueError('Vlachoutsis shear correction factors '
                                     'require positive in-plane stiffnesses')
                zn = num/den
                R = 0; I = 0; gza = 0
                for k in range(N):
                    ply = self.plies[k]
                    Gk = ply.q55L if alpha == 0 else ply.q44L
                    Dk = ply.q11L if alpha == 0 else ply.q22L
                    za = z[k]
                    zb = z[k+1]
                    R += Dk*((zb - zn)**3 - (za - zn)**3)/3.
                    zm = (za + zb)/2.
                    dz = (zb - za)/2.
                    for ig in range(3):
                        zg = zm + dz*gp[ig]
                        gz = gza - Dk*(0.5*(zg*zg - za*za) - zn*(zg - za))
                        I += gw[ig]*dz*gz*gz/Gk
                    # g at the top interface of this ply
                    gza = gza - Dk*(0.5*(zb*zb - za*za) - zn*(zb - za))
                kappa = R*R/(d*I)
                if alpha == 0:
                    k13 = kappa
                else:
                    k23 = kappa
            self.A44 = k23*self.Abar44
            self.A45 = (k13 + k23)/2.*self.Abar45
            self.A55 = k13*self.Abar55

        self.scf_k13 = self.A55/self.Abar55 if self.Abar55 != 0 else np.nan
        self.scf_k23 = self.A44/self.Abar44 if self.Abar44 != 0 else np.nan


    cpdef tuple calc_transverse_shear_stress(Laminate self, double z,
            double Qy, double Qx):
        r"""Transverse shear stresses at a given height

        Evaluates the equilibrium distribution of Rohwer (1988):

        .. math::

            \begin{Bmatrix} \tau_{yz} \\ \tau_{xz} \end{Bmatrix} =
            f^{(k)}(z) \begin{Bmatrix} Q_y \\ Q_x \end{Bmatrix}

        where `f^{(k)}(z)` is quadratic within each ply, vanishes at the bottom
        and top faces and is continuous at the ply interfaces. This is the
        consistent way to recover `\tau_{xz}` and `\tau_{yz}`, e.g. for
        failure criteria, since `C_s \gamma` is constant within each ply and
        non-zero at the free surfaces. For a homogeneous plate it gives the
        parabola `\tau_{xz} = 3 Q_x/(2h) (1 - 4 \bar{z}^2/h^2)`, with `\bar{z}`
        measured from the mid-surface.

        The distribution only depends on the in-plane stiffnesses of the plies
        and is the same for every ``shear_correction``. It is computed once by
        :meth:`.calc_transverse_shear_stiffness` when ``shear_correction``
        is ``'rohwer'``, or on the first call otherwise, and it is reset by
        :meth:`.calc_constitutive_matrix`, which must be called again if the
        plies are modified.

        Parameters
        ----------
        z : float
            Height measured from the reference surface, within `[-h/2 +
            offset, +h/2 + offset]`. At a ply interface both plies give the
            same result.
        Qy, Qx : float
            Transverse shear forces per unit length, `Q_y` and `Q_x`, e.g.
            ``{Qy, Qx} = Ats @ {gamma_yz, gamma_xz}``.

        Returns
        -------
        tau_yz, tau_xz : tuple of float
            Transverse shear stresses.

        Raises
        ------
        ValueError
            If ``z`` is outside the laminate, or if the ABD matrix of the
            laminate is singular.

        """
        cdef int k, lo, hi, mid, N
        cdef double tol
        cdef double [::1] zi
        cdef double [:, :, :, ::1] fc

        if not self._ts_ready:
            self._calc_transverse_shear_distribution()
        zi = self._ts_z
        fc = self._ts_fcoef
        N = <int>zi.shape[0] - 1
        tol = 1e-12*(zi[N] - zi[0])
        if not (zi[0] - tol <= z <= zi[N] + tol):
            raise ValueError('z=%g is outside the laminate, [%g, %g]'
                             % (z, zi[0], zi[N]))
        # last ply k with zi[k] <= z
        lo = 0
        hi = N - 1
        while lo < hi:
            mid = (lo + hi + 1)//2
            if zi[mid] <= z:
                lo = mid
            else:
                hi = mid - 1
        k = lo
        return ((fc[k, 0, 0, 0] + z*fc[k, 1, 0, 0] + z*z*fc[k, 2, 0, 0])*Qy
              + (fc[k, 0, 0, 1] + z*fc[k, 1, 0, 1] + z*z*fc[k, 2, 0, 1])*Qx,
                (fc[k, 0, 1, 0] + z*fc[k, 1, 1, 0] + z*z*fc[k, 2, 1, 0])*Qy
              + (fc[k, 0, 1, 1] + z*fc[k, 1, 1, 1] + z*z*fc[k, 2, 1, 1])*Qx)


    cpdef void calc_equivalent_properties(Laminate self):
        r"""Calculate the equivalent laminate properties

        The following attributes are updated:

            ``e1``, ``e2``, ``g12``, ```u12``, ``nu21``

        """
        AI = np.linalg.inv(self.get_ABD())
        a11, a12, a22, a33 = AI[0,0], AI[0,1], AI[1,1], AI[2,2]
        self.e1 = 1./(self.h*a11)
        self.e2 = 1./(self.h*a22)
        self.g12 = 1./(self.h*a33)
        self.nu12 = - a12 / a11
        self.nu21 = - a12 / a22


    cpdef void calc_constitutive_matrix(Laminate self):
        """Calculate the laminate constitutive terms

        This is the commonly called ``ABD`` matrix with ``shape=(6, 6)`` when
        the classical laminated plate theory is used, or the ``ABD`` matrix
        when the first-order shear deformation theory is used, containing the
        transverse shear terms.

        The transverse shear stiffnesses ``A44``, ``A45``, ``A55`` are
        calculated at the end by :meth:`.calc_transverse_shear_stiffness`,
        with the shear correction selected by ``shear_correction`` already
        applied.

        """
        cdef double h0, hk_1, hk, tmp_hk, tmp_hk_1
        self.h = 0.
        self.intrho = 0.
        self.intrhoz = 0.
        self.intrhoz2 = 0.
        for ply in self.plies:
            self.h += ply.h
        h0 = -self.h/2. + self.offset
        self.A11 = 0; self.A12 = 0; self.A16 = 0; self.A22 = 0; self.A26 = 0; self.A66 = 0
        self.B11 = 0; self.B12 = 0; self.B16 = 0; self.B22 = 0; self.B26 = 0; self.B66 = 0
        self.D11 = 0; self.D12 = 0; self.D16 = 0; self.D22 = 0; self.D26 = 0; self.D66 = 0
        self.E11 = 0; self.E12 = 0; self.E16 = 0; self.E22 = 0; self.E26 = 0; self.E66 = 0
        self.F11 = 0; self.F12 = 0; self.F16 = 0; self.F22 = 0; self.F26 = 0; self.F66 = 0
        self.H11 = 0; self.H12 = 0; self.H16 = 0; self.H22 = 0; self.H26 = 0; self.H66 = 0
        self.D44 = 0; self.D45 = 0; self.D55 = 0
        self.F44 = 0; self.F45 = 0; self.F55 = 0
        for ply in self.plies:
            hk_1 = h0
            h0 += ply.h
            hk = h0

            self.intrho += ply.matlamina.rho*(hk - hk_1)
            self.intrhoz += ply.matlamina.rho*(hk*hk/2. - hk_1*hk_1/2.)
            self.intrhoz2 += ply.matlamina.rho*(hk*hk*hk/3. - hk_1*hk_1*hk_1/3.)

            self.A11 += ply.q11L*(hk - hk_1)
            self.A12 += ply.q12L*(hk - hk_1)
            self.A16 += ply.q16L*(hk - hk_1)
            self.A22 += ply.q22L*(hk - hk_1)
            self.A26 += ply.q26L*(hk - hk_1)
            self.A66 += ply.q66L*(hk - hk_1)

            tmp_hk = hk*hk
            tmp_hk_1 = hk_1*hk_1
            self.B11 += 1/2.*ply.q11L*(tmp_hk - tmp_hk_1)
            self.B12 += 1/2.*ply.q12L*(tmp_hk - tmp_hk_1)
            self.B16 += 1/2.*ply.q16L*(tmp_hk - tmp_hk_1)
            self.B22 += 1/2.*ply.q22L*(tmp_hk - tmp_hk_1)
            self.B26 += 1/2.*ply.q26L*(tmp_hk - tmp_hk_1)
            self.B66 += 1/2.*ply.q66L*(tmp_hk - tmp_hk_1)

            tmp_hk = hk*hk*hk
            tmp_hk_1 = hk_1*hk_1*hk_1
            self.D11 += 1/3.*ply.q11L*(tmp_hk - tmp_hk_1)
            self.D12 += 1/3.*ply.q12L*(tmp_hk - tmp_hk_1)
            self.D16 += 1/3.*ply.q16L*(tmp_hk - tmp_hk_1)
            self.D22 += 1/3.*ply.q22L*(tmp_hk - tmp_hk_1)
            self.D26 += 1/3.*ply.q26L*(tmp_hk - tmp_hk_1)
            self.D66 += 1/3.*ply.q66L*(tmp_hk - tmp_hk_1)

            self.D44 += 1/3.*ply.q44L*(tmp_hk - tmp_hk_1)
            self.D45 += 1/3.*ply.q45L*(tmp_hk - tmp_hk_1)
            self.D55 += 1/3.*ply.q55L*(tmp_hk - tmp_hk_1)

            tmp_hk = hk*hk*hk*hk
            tmp_hk_1 = hk_1*hk_1*hk_1*hk_1
            self.E11 += 1/4.*ply.q11L*(tmp_hk - tmp_hk_1)
            self.E12 += 1/4.*ply.q12L*(tmp_hk - tmp_hk_1)
            self.E16 += 1/4.*ply.q16L*(tmp_hk - tmp_hk_1)
            self.E22 += 1/4.*ply.q22L*(tmp_hk - tmp_hk_1)
            self.E26 += 1/4.*ply.q26L*(tmp_hk - tmp_hk_1)
            self.E66 += 1/4.*ply.q66L*(tmp_hk - tmp_hk_1)

            tmp_hk = hk*hk*hk*hk*hk
            tmp_hk_1 = hk_1*hk_1*hk_1*hk_1*hk_1
            self.F11 += 1/5.*ply.q11L*(tmp_hk - tmp_hk_1)
            self.F12 += 1/5.*ply.q12L*(tmp_hk - tmp_hk_1)
            self.F16 += 1/5.*ply.q16L*(tmp_hk - tmp_hk_1)
            self.F22 += 1/5.*ply.q22L*(tmp_hk - tmp_hk_1)
            self.F26 += 1/5.*ply.q26L*(tmp_hk - tmp_hk_1)
            self.F66 += 1/5.*ply.q66L*(tmp_hk - tmp_hk_1)

            self.F44 += 1/5.*ply.q44L*(tmp_hk - tmp_hk_1)
            self.F45 += 1/5.*ply.q45L*(tmp_hk - tmp_hk_1)
            self.F55 += 1/5.*ply.q55L*(tmp_hk - tmp_hk_1)

            tmp_hk = hk*hk*hk*hk*hk*hk*hk
            tmp_hk_1 = hk_1*hk_1*hk_1*hk_1*hk_1*hk_1*hk_1
            self.H11 += 1/7.*ply.q11L*(tmp_hk - tmp_hk_1)
            self.H12 += 1/7.*ply.q12L*(tmp_hk - tmp_hk_1)
            self.H16 += 1/7.*ply.q16L*(tmp_hk - tmp_hk_1)
            self.H22 += 1/7.*ply.q22L*(tmp_hk - tmp_hk_1)
            self.H26 += 1/7.*ply.q26L*(tmp_hk - tmp_hk_1)
            self.H66 += 1/7.*ply.q66L*(tmp_hk - tmp_hk_1)

        self.calc_transverse_shear_stiffness()


    cpdef void make_balanced(Laminate self):
        r"""Make a balanced laminate

        The attributes `A_{16}`, `A_{26}`, `B_{16}`, `B_{26}` are set to zero
        to make a balanced laminate.

        """
        if self.offset != 0.:
            raise RuntimeError('Laminates with offset cannot be made balanced!')
        self.A16 = 0.
        self.A26 = 0.
        self.B16 = 0.
        self.B26 = 0.


    cpdef void make_orthotropic(Laminate self):
        r"""Make an orthotropic laminate

        The attributes `A_{16}`, `A_{26}`, `B_{16}`, `B_{26}`, `D_{16}`,
        `D_{26}` are set to zero to make an orthotropic laminate.

        """
        if self.offset != 0.:
            raise RuntimeError('Laminates with offset cannot be made orthotropic!')
        self.A16 = 0.
        self.A26 = 0.
        self.B16 = 0.
        self.B26 = 0.
        self.D16 = 0.
        self.D26 = 0.


    cpdef void make_symmetric(Laminate self):
        """Make a symmetric laminate

        The `B_{ij}` terms of the constitutive matrix are set to zero.

        """
        if self.offset != 0.:
            raise RuntimeError(
                    'Laminates with offset cannot be made symmetric!')
        self.B11 = 0
        self.B12 = 0
        self.B16 = 0
        self.B22 = 0
        self.B26 = 0
        self.B66 = 0

        # TSDT
        self.E11 = 0
        self.E12 = 0
        self.E16 = 0
        self.E22 = 0
        self.E26 = 0
        self.E66 = 0

        self.F16 = 0
        self.F26 = 0
        self.H16 = 0
        self.H26 = 0

        self.D45 = 0
        self.F45 = 0


    cpdef void make_smeared(Laminate self):
        r"""Make a laminated with smeared properties

        The `B_{ij}` terms of the constitutive matrix are set to zero.

        The `D_{ij}` terms are calculated from the membrane terms `A_{ij}`
        according to `D_{ij} = (h^2 A_{ij})/12`, where `h` is the
        laminate thickness.

        """
        if self.offset != 0.:
            raise NotImplementedError(
                    'Laminates with offset cannot be made smeared!')

        self.B11 = 0
        self.B12 = 0
        self.B16 = 0
        self.B22 = 0
        self.B26 = 0
        self.B66 = 0

        self.D11 = self.h**2/12 * self.A11
        self.D12 = self.h**2/12 * self.A12
        self.D16 = self.h**2/12 * self.A16
        self.D22 = self.h**2/12 * self.A22
        self.D26 = self.h**2/12 * self.A26
        self.D66 = self.h**2/12 * self.A66


    cpdef LaminationParameters calc_lamination_parameters(Laminate self):
        r"""Calculate the lamination parameters.

        The following attributes are calculated:

            ``xiA``, ``xiB``, ``xiD``, ``xiE``

        """
        cdef double h0, hk, hk_1, h, zbar1, zbar2, Afac, Bfac, Dfac, Efac
        cdef LaminationParameters lp = LaminationParameters()

        if len(self.plies) == 0:
            raise ValueError('Laminate with 0 plies!')

        h = 0.
        for ply in self.plies:
            h += ply.h

        h0 = -h/2. + self.offset
        for ply in self.plies:
            ply.rebuild()
            hk_1 = h0
            h0 += ply.h
            hk = h0
            zbar2 = hk/h
            zbar1 = hk_1/h

            Afac = zbar2 - zbar1
            Bfac = 2*(zbar2*zbar2 - zbar1*zbar1)
            Dfac = 4*(zbar2*zbar2*zbar2 - zbar1*zbar1*zbar1)
            Efac = zbar2 - zbar1

            lp.xiA1 += Afac * ply.cos2t
            lp.xiA2 += Afac * ply.sin2t
            lp.xiA3 += Afac * ply.cos4t
            lp.xiA4 += Afac * ply.sin4t

            lp.xiB1 += Bfac * ply.cos2t
            lp.xiB2 += Bfac * ply.sin2t
            lp.xiB3 += Bfac * ply.cos4t
            lp.xiB4 += Bfac * ply.sin4t

            lp.xiD1 += Dfac * ply.cos2t
            lp.xiD2 += Dfac * ply.sin2t
            lp.xiD3 += Dfac * ply.cos4t
            lp.xiD4 += Dfac * ply.sin4t

            lp.xiAtrans1 += Efac * ply.cos2t
            lp.xiAtrans2 += Efac * ply.sin2t

        return lp


cpdef LaminationParameters make_balanced_LP(LaminationParameters lp):
    r"""Make balanced lamination parameters

    The lamination parameters `\xi_{A2}` and `\xi_{A4}` are set to null to
    make a balanced laminate.

    """
    lp.xiA2 = 0
    lp.xiA4 = 0
    return lp


cpdef LaminationParameters make_symmetric_LP(LaminationParameters lp):
    r"""Make symmetric lamination parameters

    The lamination parameters `\xi_{Bi}` are set to null to make a symmetric
    laminate.

    """
    lp.xiB1 = 0
    lp.xiB2 = 0
    lp.xiB3 = 0
    lp.xiB4 = 0
    return lp


cpdef LaminationParameters make_orthotropic_LP(LaminationParameters lp):
    r"""Make orthotropic lamination parameters

    The lamination parameters `\xi_{A2}`, `\xi_{A4}`, `\xi_{B2}`, `\xi_{B4}`,
    `\xi_{D2}` and `\xi_{D4}` are set to null to make an orthotropic laminate.
    The `\xi_{D2}` and `\xi_{D4}` are related to the bend-twist coupling and
    become often very small for balanced laminates with a large amount of
    plies.

    """
    lp.xiA2 = 0
    lp.xiA4 = 0
    lp.xiB2 = 0
    lp.xiB4 = 0
    lp.xiD2 = 0
    lp.xiD4 = 0
    return lp


cpdef Laminate laminate_from_LaminationParameters(double thickness, MatLamina
        mat, LaminationParameters lp):
    r"""Return a :class:`.Laminate` object based in the thickness, material and
    lamination parameters

    Parameters
    ----------
    thickness : float
        The total thickness of the laminate
    mat : :class:`.MatLamina` object
        Material object
    lp : :class:`.LaminationParameters` object
        The container class with all lamination parameters already defined

    Returns
    -------
    lam : :class:`.Laminate`
        laminate with the constitutive matrices already calculated

    Notes
    -----
    Since the through-thickness distribution of the plies is not known from
    the lamination parameters, no shear correction can be computed. The
    transverse shear stiffnesses ``A44``, ``A45``, ``A55`` are therefore
    equal to the constant-strain ``Abar44``, ``Abar45``, ``Abar55``,
    ``shear_correction`` is ``None``, ``scf_k13 = scf_k23 = 1`` and the
    constant-stress ``Abarbar44``, ``Abarbar45``, ``Abarbar55`` are ``nan``.

    """
    lam = Laminate()
    lam.h = thickness

    lam.A11 = lam.h*(mat.u1 + mat.u2*lp.xiA1 + 0*lp.xiA2 + mat.u3*lp.xiA3 + 0*lp.xiA4)
    lam.A12 = lam.h*(mat.u4 + 0*lp.xiA1 + 0*lp.xiA2 + (-1)*mat.u3*lp.xiA3 + 0*lp.xiA4)
    lam.A22 = lam.h*(mat.u1 + (-1)*mat.u2*lp.xiA1 + 0*lp.xiA2 + mat.u3*lp.xiA3 + 0*lp.xiA4)
    lam.A16 = lam.h*(0 + 0*lp.xiA1 + mat.u2/2.*lp.xiA2 + 0*lp.xiA3 + mat.u3*lp.xiA4)
    lam.A26 = lam.h*(0 + 0*lp.xiA1 + mat.u2/2.*lp.xiA2 + 0*lp.xiA3 + (-1)*mat.u3*lp.xiA4)
    lam.A66 = lam.h*(mat.u5 + 0*lp.xiA1 + 0*lp.xiA2 + (-1)*mat.u3*lp.xiA3 + 0*lp.xiA4)

    lam.B11 = lam.h*lam.h/4.*(mat.u2*lp.xiB1 + 0*lp.xiB2 + mat.u3*lp.xiB3 + 0*lp.xiB4)
    lam.B12 = lam.h*lam.h/4.*(0*lp.xiB1 + 0*lp.xiB2 + (-1)*mat.u3*lp.xiB3 + 0*lp.xiB4)
    lam.B22 = lam.h*lam.h/4.*((-1)*mat.u2*lp.xiB1 + 0*lp.xiB2 + mat.u3*lp.xiB3 + 0*lp.xiB4)
    lam.B16 = lam.h*lam.h/4.*(0*lp.xiB1 + mat.u2/2.*lp.xiB2 + 0*lp.xiB3 + mat.u3*lp.xiB4)
    lam.B26 = lam.h*lam.h/4.*(0*lp.xiB1 + mat.u2/2.*lp.xiB2 + 0*lp.xiB3 + (-1)*mat.u3*lp.xiB4)
    lam.B66 = lam.h*lam.h/4.*(0*lp.xiB1 + 0*lp.xiB2 + (-1)*mat.u3*lp.xiB3 + 0*lp.xiB4)

    lam.D11 = lam.h*lam.h*lam.h/12.*(mat.u1 + mat.u2*lp.xiD1 + 0*lp.xiD2 + mat.u3*lp.xiD3 + 0*lp.xiD4)
    lam.D12 = lam.h*lam.h*lam.h/12.*(mat.u4 + 0*lp.xiD1 + 0*lp.xiD2 + (-1)*mat.u3*lp.xiD3 + 0*lp.xiD4)
    lam.D22 = lam.h*lam.h*lam.h/12.*(mat.u1 + (-1)*mat.u2*lp.xiD1 + 0*lp.xiD2 + mat.u3*lp.xiD3 + 0*lp.xiD4)
    lam.D16 = lam.h*lam.h*lam.h/12.*(0 + 0*lp.xiD1 + mat.u2/2.*lp.xiD2 + 0*lp.xiD3 + mat.u3*lp.xiD4)
    lam.D26 = lam.h*lam.h*lam.h/12.*(0 + 0*lp.xiD1 + mat.u2/2.*lp.xiD2 + 0*lp.xiD3 + (-1)*mat.u3*lp.xiD4)
    lam.D66 = lam.h*lam.h*lam.h/12.*(mat.u5 + 0*lp.xiD1 + 0*lp.xiD2 + (-1)*mat.u3*lp.xiD3 + 0*lp.xiD4)

    lam.Abar44 = lam.h*(mat.u6 + mat.u7*lp.xiAtrans1 + 0*lp.xiAtrans2)
    lam.Abar45 = lam.h*(0 + 0*lp.xiAtrans1 + (-1)*mat.u7*lp.xiAtrans2)
    lam.Abar55 = lam.h*(mat.u6 + (-1)*mat.u7*lp.xiAtrans1 + 0*lp.xiAtrans2)
    # NOTE the through-thickness ply distribution is not known, such that no
    #      shear correction can be computed
    lam.shear_correction = None
    lam.A44 = lam.Abar44
    lam.A45 = lam.Abar45
    lam.A55 = lam.Abar55
    lam.Abarbar44 = np.nan
    lam.Abarbar45 = np.nan
    lam.Abarbar55 = np.nan
    lam.scf_k13 = 1.
    lam.scf_k23 = 1.

    return lam


cpdef Laminate laminate_from_lamination_parameters(double thickness, MatLamina
        matlamina, double xiA1, double xiA2, double xiA3, double xiA4,
        double xiB1, double xiB2, double xiB3, double xiB4,
        double xiD1, double xiD2, double xiD3, double xiD4,
        double xiAtrans1=0, double xiAtrans2=0):
    r"""Return a :class:`.Laminate` object based in the thickness, material and
    lamination parameters

    Note that `\xi_{E1}` and `\xi_{E2}` are optional and usually equal to zero,
    becoming important only when the transverse shear modulus is different in
    the two directions, i.e.  when `G_{13} \ne G{23}`.

    Parameters
    ----------
    thickness : float
        The total thickness of the laminate
    matlamina : :class:`.MatLamina` object
        Material object
    xiAj, xiBj, xiDj, xiEj : float
        The 14 lamination parameters according to the first-order shear
        deformation theory: `\xi_{A1} \cdots \xi_{A4}`, `\xi_{B1} \cdots
        \xi_{B4}`, `\xi_{D1} \cdots \xi_{D4}`, `\xi_{E1}` and `\xi_{E2}`


    Returns
    -------
    lam : :class:`.Laminate`
        laminate with the constitutive matrices already calculated. See
        :func:`.laminate_from_LaminationParameters` for the transverse shear
        stiffnesses.

    """
    lp = LaminationParameters()
    lp.xiA1 = xiA1
    lp.xiA2 = xiA2
    lp.xiA3 = xiA3
    lp.xiA4 = xiA4
    lp.xiB1 = xiB1
    lp.xiB2 = xiB2
    lp.xiB3 = xiB3
    lp.xiB4 = xiB4
    lp.xiD1 = xiD1
    lp.xiD2 = xiD2
    lp.xiD3 = xiD3
    lp.xiD4 = xiD4
    lp.xiAtrans1 = xiAtrans1
    lp.xiAtrans2 = xiAtrans2
    return laminate_from_LaminationParameters(thickness, matlamina, lp)


_GRADABD_STATE = ('gradAij', 'gradBij', 'gradDij', 'gradAtransij')


cdef class GradABD:
    r"""Container to store the gradients of the ABD matrices with respect to
    the lamination parameters

    Attributes
    ==========

    gradAij, gradBij, gradDij, gradAtransij : tuple of 2D np.array objects
        The shapes of these gradient matrices are:

            gradAij: (6, 5)
            gradBij: (6, 5)
            gradDij: (6, 5)
            gradAtransij: (3, 3)

        They contain the gradients of each laminate stiffness with respect to
        the thickness and respective lamination parameters. The transverse
        shear terms are the constant-strain (uncorrected) ``Abar44``,
        ``Abar45``, ``Abar55``, as given by
        :func:`.laminate_from_LaminationParameters`. The rows and
        columns correspond to::

            gradAij
            -------

                h xiA1 xiA2 xiA3 xiA4
            A11
            A12
            A16
            A22
            A26
            A66

            gradBij
            -------

                h xiB1 xiB2 xiB3 xiB4
            B11
            B12
            B16
            B22
            B26
            B66

            gradDij
            -------

                h xiD1 xiD2 xiD3 xiD4
            D11
            D12
            D16
            D22
            D26
            D66

            gradAtransij
            -------

                h xiAtrans1 xiAtrans2
            A44
            A45
            A55

    """
    def __init__(GradABD self):
        self.gradAij = np.zeros((6, 5), dtype=DOUBLE)
        self.gradBij = np.zeros((6, 5), dtype=DOUBLE)
        self.gradDij = np.zeros((6, 5), dtype=DOUBLE)
        self.gradAtransij = np.zeros((3, 3), dtype=DOUBLE)

    def __reduce__(GradABD self):
        state = {}
        for name in _GRADABD_STATE:
            try:
                state[name] = np.asarray(getattr(self, name)).copy()
            except AttributeError: # memoryview not initialized
                pass
        return copyreg.__newobj__, (type(self), ), state

    def __setstate__(GradABD self, dict state):
        for name, value in state.items():
            setattr(self, name, np.ascontiguousarray(value, dtype=DOUBLE))

    cpdef void calc_LP_grad(GradABD self, double thickness, MatLamina mat, LaminationParameters lp):
        r"""Gradients of the shell stiffnesses with respect to the thickness and
        lamination parameters

        Parameters
        ----------
        thickness : float
            The total thickness of the laminate
        mat : :class:`.MatLamina` object
            Material object
        lp : :class:`.LaminationParameters` object
            The container class with all lamination parameters already defined

        Returns
        -------
        None
            The attributes of the object are updated.


        """
        cdef int i, j
        cdef double h
        cdef double [:, ::1] gradinv

        h = thickness
        gradinv = np.array([[mat.u2, 0, mat.u3, 0],
                            [0, 0, -mat.u3, 0],
                            [0, mat.u2/2., 0, mat.u3],
                            [-mat.u2, 0, mat.u3, 0],
                            [0, mat.u2/2., 0, -mat.u3],
                            [0, 0, -mat.u3, 0]], dtype=DOUBLE)

        # d(A11 A12 A16 A22 A26 A66) / dh
        self.gradAij[0, 0] = (mat.u1 + mat.u2*lp.xiA1 + 0*lp.xiA2 + mat.u3*lp.xiA3 + 0*lp.xiA4)
        self.gradAij[1, 0] = (mat.u4 + 0*lp.xiA1 + 0*lp.xiA2 + (-1)*mat.u3*lp.xiA3 + 0*lp.xiA4)
        self.gradAij[2, 0] = (0 + 0*lp.xiA1 + mat.u2/2.*lp.xiA2 + 0*lp.xiA3 + mat.u3*lp.xiA4)
        self.gradAij[3, 0] = (mat.u1 + (-1)*mat.u2*lp.xiA1 + 0*lp.xiA2 + mat.u3*lp.xiA3 + 0*lp.xiA4)
        self.gradAij[4, 0] = (0 + 0*lp.xiA1 + mat.u2/2.*lp.xiA2 + 0*lp.xiA3 + (-1)*mat.u3*lp.xiA4)
        self.gradAij[5, 0] = (mat.u5 + 0*lp.xiA1 + 0*lp.xiA2 + (-1)*mat.u3*lp.xiA3 + 0*lp.xiA4)

        # d(A11 A12 A16 A22 A26 A66) / d(xiA1, xiA2, xiA3, xiA4)
        for i in range(5):
            for j in range(4):
                self.gradAij[i, j+1] = h*gradinv[i, j]

        # d(B11 B12 B16 B22 B26 B66) / dh
        self.gradBij[0, 0] = h/2.*(mat.u2*lp.xiB1 + 0*lp.xiB2 + mat.u3*lp.xiB3 + 0*lp.xiB4)
        self.gradBij[1, 0] = h/2.*(0*lp.xiB1 + 0*lp.xiB2 + (-1)*mat.u3*lp.xiB3 + 0*lp.xiB4)
        self.gradBij[2, 0] = h/2.*(0*lp.xiB1 + mat.u2/2.*lp.xiB2 + 0*lp.xiB3 + mat.u3*lp.xiB4)
        self.gradBij[3, 0] = h/2.*((-1)*mat.u2*lp.xiB1 + 0*lp.xiB2 + mat.u3*lp.xiB3 + 0*lp.xiB4)
        self.gradBij[4, 0] = h/2.*(0*lp.xiB1 + mat.u2/2.*lp.xiB2 + 0*lp.xiB3 + (-1)*mat.u3*lp.xiB4)
        self.gradBij[5, 0] = h/2.*(0*lp.xiB1 + 0*lp.xiB2 + (-1)*mat.u3*lp.xiB3 + 0*lp.xiB4)

        # d(B11 B12 B16 B22 B26 B66) / d(xiB1, xiB2, xiB3, xiB4)
        for i in range(5):
            for j in range(4):
                self.gradBij[i, j+1] = h*h/4.*gradinv[i, j]

        # d(D11 D12 D16 D22 D26 D66) / dh
        self.gradDij[0, 0] = h*h/4.*(mat.u1 + mat.u2*lp.xiD1 + 0*lp.xiD2 + mat.u3*lp.xiD3 + 0*lp.xiD4)
        self.gradDij[1, 0] = h*h/4.*(mat.u4 + 0*lp.xiD1 + 0*lp.xiD2 + (-1)*mat.u3*lp.xiD3 + 0*lp.xiD4)
        self.gradDij[2, 0] = h*h/4.*(0 + 0*lp.xiD1 + mat.u2/2.*lp.xiD2 + 0*lp.xiD3 + mat.u3*lp.xiD4)
        self.gradDij[3, 0] = h*h/4.*(mat.u1 + (-1)*mat.u2*lp.xiD1 + 0*lp.xiD2 + mat.u3*lp.xiD3 + 0*lp.xiD4)
        self.gradDij[4, 0] = h*h/4.*(0 + 0*lp.xiD1 + mat.u2/2.*lp.xiD2 + 0*lp.xiD3 + (-1)*mat.u3*lp.xiD4)
        self.gradDij[5, 0] = h*h/4.*(mat.u5 + 0*lp.xiD1 + 0*lp.xiD2 + (-1)*mat.u3*lp.xiD3 + 0*lp.xiD4)

        # d(D11 D12 D16 D22 D26 D66) / d(xiD1, xiD2, xiD3, xiD4)
        for i in range(5):
            for j in range(4):
                self.gradDij[i, j+1] = h*h*h/12.*gradinv[i, j]

        # d(A44, A45, A55) / dh
        self.gradAtransij[0, 0] = (mat.u6 + mat.u7*lp.xiAtrans1 + 0*lp.xiAtrans2)
        self.gradAtransij[1, 0] = (0 + 0*lp.xiAtrans1 + (-1)*mat.u7*lp.xiAtrans2)
        self.gradAtransij[2, 0] = (mat.u6 + (-1)*mat.u7*lp.xiAtrans1 + 0*lp.xiAtrans2)

        # d(A44, A45, A55) / d(xiAtrans1, xiAtrans2)
        self.gradAtransij[0, 1] = h*mat.u7
        self.gradAtransij[1, 2] = h*(-mat.u7)
        self.gradAtransij[2, 1] = h*(-mat.u7)


cpdef Laminate n_double_laminate(double thickness, int n, double[::1] angles_deg, MatLamina matlamina):
    r"""Create a N-double laminated plate using a faster code

    This code is considerably faster than :func:`composites.utils.n_double_laminate`.

    An N-double laminate consists of `[\pm\phi_1,\pm\phi_2, \cdots,
    \pm\phi_n]`. With the principle of homogenization, at the limit where many
    plies are used we have that `B=0`.  Based on the double-double laminate as
    described by:

        Shrivastava, S., Sharma, N., Tsai, S. W., and Mohite, P. M., 2020,
        “D and DD-Drop Layup Optimization of Aircraft Wing Panels under
        Multi-Load Case Design Environment,” Compos. Struct., 248(January), p.
        112518.

    Parameters
    ----------
    thickness : float
        Total plate thickness.
    n : int
        Number of angle pairs,m defines the "N" in "N-double".
    angles_deg : array-like
        List of `\phi_n` of the N-double laminate.
    matlamina : MatLamina
        See :class:`.MatLamina` for details.

    """
    cdef int i
    cdef double tr

    lam = Laminate()
    lam.h = thickness

    tr = matlamina.q11 + matlamina.q22 + 2*matlamina.q66
    matlamina.trace_normalize_plane_stress()
    A11_star = matlamina.u1
    A12_star = matlamina.u4
    A22_star = matlamina.u1
    A16_star = 0
    A26_star = 0
    A66_star = matlamina.u5
    for i in range(n):
        angle_deg = angles_deg[i]
        angle = deg2rad(angle_deg)
        A11_star += matlamina.u2*(cos(2*angle) + cos(-2*angle))/(2*n)
        A11_star += matlamina.u3*(cos(4*angle) + cos(-4*angle))/(2*n)

        A12_star -= matlamina.u3*(cos(4*angle) + cos(-4*angle))/(2*n)

        A22_star -= matlamina.u2*(cos(2*angle) + cos(-2*angle))/(2*n)
        A22_star += matlamina.u3*(cos(4*angle) + cos(-4*angle))/(2*n)

        A66_star -= matlamina.u3*(cos(4*angle) + cos(-4*angle))/(2*n)

    lam.A11 = tr*A11_star*lam.h
    lam.A12 = tr*A12_star*lam.h
    lam.A16 = tr*A16_star*lam.h
    lam.A22 = tr*A22_star*lam.h
    lam.A26 = tr*A26_star*lam.h
    lam.A66 = tr*A66_star*lam.h
    lam.D11 = tr*A11_star*lam.h**3/12.
    lam.D12 = tr*A12_star*lam.h**3/12.
    lam.D16 = tr*A16_star*lam.h**3/12.
    lam.D22 = tr*A22_star*lam.h**3/12.
    lam.D26 = tr*A26_star*lam.h**3/12.
    lam.D66 = tr*A66_star*lam.h**3/12.
    lam.calc_equivalent_properties()

    return lam
