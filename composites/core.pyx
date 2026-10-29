#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
#cython: nonecheck=False
#cython: initializedcheck=False
#cython: infer_types=False
#cython: overflowcheck=False
"""
===============================================
Composites Core Module (:mod:`composites.core`)
===============================================

.. currentmodule:: composites.core

Classes of the material (:class:`.MatLamina`), of the ply (:class:`.Lamina`)
and of the laminate (:class:`.Laminate`), of the lamination parameters
(:class:`.LaminationParameters`) and of their gradients (:class:`.GradABD`),
and functions that create a :class:`.Laminate` from lamination parameters or
from the angles of an N-double laminate. The functions of
:mod:`composites.utils` create these objects from the ply data.

"""
import copyreg
import warnings

import numpy as np

DOUBLE = np.float64


cdef class LaminationParameters:
    r"""Lamination parameters

    Container of the 14 lamination parameters, all dimensionless and
    initialized to zero. It is returned by
    :meth:`.Laminate.calc_lamination_parameters` and used by
    :func:`.laminate_from_LaminationParameters` and
    :meth:`.GradABD.calc_LP_grad`. The deprecated properties ``xiAtrans1``
    and ``xiAtrans2`` are aliases of ``xiAts1`` and ``xiAts2``.

    Attributes
    ----------
    xiA1, xiA2, xiA3, xiA4 : float
        Lamination parameters `\xi_{Ai}` (in-plane), dimensionless
    xiB1, xiB2, xiB3, xiB4 : float
        Lamination parameters `\xi_{Bi}` (in-plane coupling with bending),
        dimensionless
    xiD1, xiD2, xiD3, xiD4 : float
        Lamination parameters `\xi_{Di}` (bending), dimensionless
    xiAts1, xiAts2 : float
        Lamination parameters `\xi_{{A_{ts}}i}` (transverse shear),
        dimensionless

    Notes
    -----
    The lamination parameters are defined in
    :meth:`.Laminate.calc_lamination_parameters`, and the stiffnesses that
    they give in :func:`.laminate_from_LaminationParameters`.

    """
    def __init__(LaminationParameters self):
        self.xiA1=0; self.xiA2=0; self.xiA3=0; self.xiA4=0
        self.xiB1=0; self.xiB2=0; self.xiB3=0; self.xiB4=0
        self.xiD1=0; self.xiD2=0; self.xiD3=0; self.xiD4=0
        self.xiAts1=0; self.xiAts2=0

    @property
    def xiAtrans1(self):
        r"""Deprecated, use :attr:`.xiAts1` instead

        Issues a :class:`DeprecationWarning` and returns ``xiAts1``.

        Returns
        -------
        xiAts1 : float
            Lamination parameter `\xi_{A_{ts}1}`, dimensionless.

        """
        _warn_deprecated('LaminationParameters.xiAtrans1',
                         'LaminationParameters.xiAts1')
        return self.xiAts1
    @xiAtrans1.setter
    def xiAtrans1(self, double value):
        r"""Deprecated, use :attr:`.xiAts1` instead

        Issues a :class:`DeprecationWarning` and sets ``xiAts1``.

        Parameters
        ----------
        value : float
            Lamination parameter `\xi_{A_{ts}1}`, dimensionless.

        """
        _warn_deprecated('LaminationParameters.xiAtrans1',
                         'LaminationParameters.xiAts1')
        self.xiAts1 = value
    @property
    def xiAtrans2(self):
        r"""Deprecated, use :attr:`.xiAts2` instead

        Issues a :class:`DeprecationWarning` and returns ``xiAts2``.

        Returns
        -------
        xiAts2 : float
            Lamination parameter `\xi_{A_{ts}2}`, dimensionless.

        """
        _warn_deprecated('LaminationParameters.xiAtrans2',
                         'LaminationParameters.xiAts2')
        return self.xiAts2
    @xiAtrans2.setter
    def xiAtrans2(self, double value):
        r"""Deprecated, use :attr:`.xiAts2` instead

        Issues a :class:`DeprecationWarning` and sets ``xiAts2``.

        Parameters
        ----------
        value : float
            Lamination parameter `\xi_{A_{ts}2}`, dimensionless.

        """
        _warn_deprecated('LaminationParameters.xiAtrans2',
                         'LaminationParameters.xiAts2')
        self.xiAts2 = value


def _warn_deprecated(old, new):
    # NOTE stacklevel=1 since the compiled functions add no Python frames,
    #      pointing to the line of the caller
    warnings.warn("'%s' is deprecated, use '%s' instead" % (old, new),
                  DeprecationWarning, stacklevel=1)


cdef class MatLamina:
    r"""
    Orthotropic material lamina

    Container of the elastic, thermal, density and strength properties of a
    material, created with all attributes equal to zero, except
    ``plane_stress = True``. After the elastic properties are set,
    :meth:`rebuild` computes the constitutive constants ``q11`` to ``q66``,
    the stiffness constants ``c11`` to ``c66`` and the invariants ``u1`` to
    ``u7``, all with the dimension of force per unit area;
    :meth:`trace_normalize_plane_stress` divides some of them by the trace,
    such that they become dimensionless. Usually created by
    :func:`.read_laminaprop`. By default, the constitutive constants and the
    invariants are those of plane stress, used by the plate and shell
    theories of :class:`.Laminate`, also when ``e3``, ``nu13`` and ``nu23``
    are given, see :meth:`rebuild`.

    Attributes
    ----------

    e1 : float
        Young Modulus in direction 1, force per unit area
    e2 : float
        Young Modulus in direction 2, force per unit area
    e3 : float
        Young Modulus in direction 3, force per unit area
    g12 : float
        in-plane shear modulus, force per unit area
    g13 : float
        transverse shear modulus for plane 1-Z, force per unit area
    g23 : float
        transverse shear modulus for plane 2-Z, force per unit area
    nu12 : float
        Poisson's ratio 12, dimensionless
    nu13 : float
        Poisson's ratio 13, dimensionless
    nu23 : float
        Poisson's ratio 23, dimensionless
    nu21 : float
        Poisson's ratio 21: use formula nu12/e1 = nu21/e2, dimensionless
    nu31 : float
        Poisson's ratio 31: use formula nu31/e3 = nu13/e1, dimensionless
    nu32 : float
        Poisson's ratio 32: use formula nu23/e2 = nu32/e3, dimensionless
    rho : float
        especific mass (mass / volume)
    a1 : float
        thermal expansion coeffiecient in direction 1, per unit temperature
    a2 : float
        thermal expansion coeffiecient in direction 2, per unit temperature
    a3 : float
        thermal expansion coeffiecient in direction 3, per unit temperature
    tref : float
        reference temperature
    st1,st2 : float
        allowable tensile stresses for directions 1 and 2, force per unit
        area
    sc1,sc2 : float
        allowable compressive stresses for directions 1 and 2, force per unit
        area
    ss12 : float
        allowable in-plane stress for shear, force per unit area
    plane_stress : bool
        ``True`` when the constitutive constants and the invariants are
        those of plane stress, ``False`` when they are three-dimensional, as
        selected in the last call of :meth:`rebuild`. Read-only; it is
        pickled and saved to JSON with the other attributes.
    q11 : float
        lamina constitutive constant 11
    q12 : float
        lamina constitutive constant 12
    q13 : float
        lamina constitutive constant 13
    q21 : float
        lamina constitutive constant 21
    q22 : float
        lamina constitutive constant 22
    q23 : float
        lamina constitutive constant 23
    q31 : float
        lamina constitutive constant 31
    q32 : float
        lamina constitutive constant 32
    q33 : float
        lamina constitutive constant 33
    q44 : float
        lamina constitutive constant 44
    q55 : float
        lamina constitutive constant 55
    q66 : float
        lamina constitutive constant 66
    c11, c12, c13, c22, c23, c33, c44, c55, c66 : float
        lamina stiffness constants
    u1, u2, u3, u4, u5, u6, u7 : float
        lamina material invariants

    Notes
    -----
    :class:`.MatLamina` recomputes neither `G` nor `E`: :meth:`rebuild`
    uses the moduli and the Poisson's ratios as they are set. For an
    isotropic material given as ``laminaprop = (e, nu)``,
    :func:`.read_laminaprop` sets ``g12 = g13 = g23 = e/(2*(1 + nu))``,
    i.e. `G = E/(2(1 + \nu))`; with 6 or 9 entries of ``laminaprop`` the
    shear moduli are taken as given. Otherwise, the user must set the shear
    moduli consistently with `E` and `\nu`.

    """
    def __init__(MatLamina self):
        self.plane_stress = True

    def _set_plane_stress(MatLamina self, bint plane_stress):
        # NOTE used by composites.json_io to restore the read-only attribute
        #      together with the stored constitutive constants
        self.plane_stress = plane_stress

    cpdef void rebuild(MatLamina self, bint plane_stress=True) except *:
        r"""Update constitutive and invariant terms

        Computes, from the moduli ``e1``, ``e2``, ``e3``, ``g12``, ``g13``,
        ``g23`` and the Poisson's ratios ``nu12``, ``nu21``, ``nu13``,
        ``nu31``, ``nu23``, ``nu32``, the stiffness constants ``c11`` to
        ``c66``, the constitutive constants ``q11`` to ``q66`` and the
        invariants ``u1`` to ``u7``, with the dimension of force per unit
        area. The moduli may be zero, e.g. ``e1 = e2 = g12 = 0`` for an
        antiplane core, see :func:`.read_laminaprop` and the Notes.

        Parameters
        ----------
        plane_stress : bool, optional
            ``True`` (default): ``q11`` to ``q33`` and the invariants ``u1``
            to ``u5`` are those of plane stress, `\sigma_{33} = 0`, used by
            :class:`.Laminate`, whether or not ``e3``, ``nu13`` and ``nu23``
            are given. ``False``: they are three-dimensional, e.g. for
            continuum shell elements, see the Notes. The choice is stored in
            the attribute ``plane_stress``. The stiffness constants ``c11``
            to ``c66`` are three-dimensional in both cases.

        Notes
        -----
        With ``plane_stress=True``, the default, the constitutive constants
        are the plane-stress reduced stiffnesses of Jones [3]_, Sec. 2.5,
        Eq. (2.66), p. 72, and of Reddy [1]_, Sec. 2.4, Eq. (2.4.4), p. 101:

        .. math::

            q_{11} = \frac{e_1}{1 - \nu_{12}\nu_{21}}, \quad
            q_{12} = q_{21} = \frac{\nu_{21} e_1}{1 - \nu_{12}\nu_{21}},
            \quad
            q_{22} = \frac{e_2}{1 - \nu_{12}\nu_{21}} \\
            q_{44} = g_{23}, \quad q_{55} = g_{13}, \quad q_{66} = g_{12}

        and `q_{13} = q_{23} = q_{31} = q_{32} = q_{33} = 0`. They follow
        from the three-dimensional stiffnesses with `\sigma_{33} = 0`,
        `q_{ij} = c_{ij} - c_{i3} c_{j3}/c_{33}`, `i, j = 1, 2` (Jones [3]_,
        Eq. (2.64), p. 72), and do not depend on ``e3``, ``nu13`` and
        ``nu23``: "the reduced stiffnesses involve six independent
        engineering constants: `E_1`, `E_2`, `\nu_{12}`, `G_{12}`, `G_{13}`,
        and `G_{23}`" (Reddy [1]_, p. 101). They are the stiffnesses of
        :meth:`.Lamina.rebuild`, such that
        :func:`.laminate_from_LaminationParameters` with the lamination
        parameters of :meth:`.Laminate.calc_lamination_parameters` gives the
        `A_{ij}`, `B_{ij}` and `D_{ij}` of :func:`.laminated_plate` for any
        ``laminaprop``.

        With ``plane_stress=False`` and `\Delta = 1 - \nu_{12}\nu_{21} -
        \nu_{13}\nu_{31} - \nu_{23}\nu_{32} - \nu_{12}\nu_{23}\nu_{31} -
        \nu_{13}\nu_{21}\nu_{32}`, the constitutive constants are the
        three-dimensional stiffnesses of the orthotropic material in its
        principal axes [2]_:

        .. math::

            q_{11} = e_1 (1 - \nu_{23}\nu_{32})/\Delta, \quad
            q_{12} = e_1 (\nu_{21} + \nu_{23}\nu_{31})/\Delta, \quad
            q_{13} = e_1 (\nu_{31} + \nu_{21}\nu_{32})/\Delta \\
            q_{21} = e_2 (\nu_{12} + \nu_{13}\nu_{32})/\Delta, \quad
            q_{22} = e_2 (1 - \nu_{13}\nu_{31})/\Delta, \quad
            q_{23} = e_2 (\nu_{32} + \nu_{12}\nu_{31})/\Delta \\
            q_{31} = e_3 (\nu_{13} + \nu_{12}\nu_{23})/\Delta, \quad
            q_{32} = e_3 (\nu_{23} + \nu_{13}\nu_{21})/\Delta, \quad
            q_{33} = e_3 (1 - \nu_{12}\nu_{21})/\Delta \\
            q_{44} = g_{23}, \quad q_{55} = g_{13}, \quad q_{66} = g_{12}

        which are symmetric, `q_{ij} = q_{ji}`, when the Poisson's ratios
        satisfy the reciprocity `\nu_{ij}/e_i = \nu_{ji}/e_j`, see
        :func:`.read_laminaprop`; the matrix of `q_{ij}`, `i, j = 1, 2, 3`, is
        then the inverse of the upper-left block of the compliance matrix
        given below, i.e. `q_{ij} = c_{ij}`. With `e_3 = 0`, as given by
        :func:`.read_laminaprop` for 2 or 6 entries, `\nu_{31} = \nu_{32} =
        0` and they are equal to the plane-stress ones.

        The stiffness constants are computed with `\delta = (1 -
        \nu_{12}\nu_{21} - \nu_{23}\nu_{32} - \nu_{31}\nu_{13} -
        2\nu_{21}\nu_{32}\nu_{13})/(e_1 e_2)` [1]_:

        .. math::

            c_{11} = (1 - \nu_{23}\nu_{32})/(\delta e_2), \quad
            c_{12} = (\nu_{21} + \nu_{31}\nu_{23})/(\delta e_2), \quad
            c_{13} = (\nu_{31} + \nu_{21}\nu_{32})/(\delta e_2) \\
            c_{22} = (1 - \nu_{13}\nu_{31})/(\delta e_1), \quad
            c_{23} = (\nu_{32} + \nu_{12}\nu_{31})/(\delta e_1), \quad
            c_{33} = e_3 (1 - \nu_{12}\nu_{21})/(\delta e_1 e_2) \\
            c_{44} = g_{23}, \quad c_{55} = g_{13}, \quad c_{66} = g_{12}

        With the reciprocity, `\delta e_1 e_2 = \Delta` and `c_{ij} =
        q_{ij}`. The matrix of the stiffness constants, see
        :meth:`get_constitutive_matrix`, is then the inverse of the
        orthotropic compliance matrix, with `S_{11} = 1/e_1`, `S_{22} =
        1/e_2`, `S_{33} = 1/e_3`, `S_{12} = -\nu_{12}/e_1`, `S_{13} =
        -\nu_{13}/e_1`, `S_{23} = -\nu_{23}/e_2`, `S_{44} = 1/g_{23}`,
        `S_{55} = 1/g_{13}` and `S_{66} = 1/g_{12}` [1]_.

        A zero modulus is the limit of a positive-definite material, where
        the stiffnesses stay finite while the compliances diverge. With
        `e_1 = 0` or `e_2 = 0`, the stiffness constants are computed as
        `c_{11} = e_1 (1 - \nu_{23}\nu_{32})/(\delta e_1 e_2)`, ..., `c_{33}
        = e_3 (1 - \nu_{12}\nu_{21})/(\delta e_1 e_2)`, with `\delta e_1 e_2
        = 1 - \nu_{12}\nu_{21} - \nu_{23}\nu_{32} - \nu_{31}\nu_{13} -
        2\nu_{21}\nu_{32}\nu_{13}`, the same values without dividing by the
        moduli, and ``ValueError`` is raised if this determinant is not
        positive; the compliance matrix is then undefined and the stiffness
        matrix singular. E.g. an antiplane core given as ``(0, 0, nu12, 0,
        g13, g23, e3, 0, 0)``, see :func:`.read_laminaprop`, has `c_{33} =
        e_3`, `c_{44} = g_{23}`, `c_{55} = g_{13}` and all other constants
        zero, and its plane-stress `q_{11} = q_{12} = q_{22} = q_{66} = 0`
        and `u_1 = \dots = u_5 = 0`, such that it adds only transverse shear
        stiffness to a laminate.

        The invariants `u_1` to `u_5` are those of Tsai and Pagano, as given
        by Jones [3]_, Sec. 2.7, Eq. (2.99), p. 86, and `u_6`, `u_7` the
        transverse shear invariants of the ABAQUS manual [4]_:

        .. math::

            u_1 = (3 q_{11} + 3 q_{22} + 2 q_{12} + 4 q_{66})/8, \quad
            u_2 = (q_{11} - q_{22})/2 \\
            u_3 = (q_{11} + q_{22} - 2 q_{12} - 4 q_{66})/8, \quad
            u_4 = (q_{11} + q_{22} + 6 q_{12} - 4 q_{66})/8 \\
            u_5 = (u_1 - u_4)/2, \quad
            u_6 = (q_{44} + q_{55})/2, \quad
            u_7 = (q_{44} - q_{55})/2

        They give the stiffnesses of a ply at any angle, see
        :meth:`get_invariant_matrix`.

        With ``plane_stress=True``, the invariants are those of the
        plane-stress reduced stiffnesses, as in Jones [3]_, and give the
        plane-stress stiffnesses `q_{ijL}` of the ply in the laminate axes of
        :meth:`.Lamina.rebuild`, Jones [3]_, Eqs. (2.85) and (2.98), pp.
        77 and 85. This is what the laminate stiffnesses need, since the
        stress resultants integrate the plane-stress `q_{ijL}` (Jones [3]_,
        Eq. (4.24), p. 198; Reddy [1]_, Eq. (3.3.38), p. 128).

        With ``plane_stress=False`` and `e_3 \ne 0`, `q_{11}`, `q_{12}`,
        `q_{22}` are the three-dimensional stiffnesses `c_{11}`, `c_{12}`,
        `c_{22}`, and the invariants give the components 11, 22, 12, 66, 16,
        26, 44, 45 and 55 of the three-dimensional stiffness matrix of
        :meth:`get_constitutive_matrix` rotated by the ply angle about the 3
        axis, whose in-plane components transform with the same formulas as
        the plane-stress ones, but are not the plane-stress stiffnesses
        (Tsai [5]_, Sec. 3.3, "Transformation in three dimensions"). They are
        stiffer than the plane-stress `q_{ijL}`, e.g. `\lambda + 2\mu = E (1 -
        \nu)/((1 + \nu)(1 - 2\nu))` instead of `E/(1 - \nu^2)` for an
        isotropic material, and are meant
        for three-dimensional properties, e.g. for continuum shell elements
        integrated with a user material (UMAT), not for the plate and shell
        stiffnesses of :class:`.Laminate`.

        References
        ----------
        .. [1] Reddy, J. N., Mechanics of Laminated Composite Plates and Shells:
           Theory and Analysis, Second Edition, CRC Press, 2003.
           https://doi.org/10.1201/b12409
        .. [2] Gürdal, Z., Haftka, R. T. and Hajela, P., Design and
           Optimization of Laminated Composite Materials, Wiley-Interscience,
           1999.
        .. [3] Jones, R. M., Mechanics of Composite Materials, Second
           Edition, Taylor & Francis, Philadelphia, 1999.
        .. [4] ABAQUS Example Problems Manual, Vol. 1, Example 1.2.2,
           Laminated composite shell: buckling of a cylindrical panel with a
           circular hole.
        .. [5] Tsai, S. W., Theory of Composites Design, 2008.

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
        if e1 != 0 and e2 != 0:
            delta = (1-nu12*nu21-nu23*nu32-nu31*nu13-2*nu21*nu32*nu13)/(e1*e2)
            self.c11 = (1    - nu23*nu32)/(delta*e2)
            self.c12 = (nu21 + nu31*nu23)/(delta*e2)
            self.c13 = (nu31 + nu21*nu32)/(delta*e2)
            self.c22 = (1    - nu13*nu31)/(delta*e1)
            self.c23 = (nu32 + nu12*nu31)/(delta*e1)
            self.c33 = e3*(1    - nu12*nu21)/(delta*e1*e2)
        else:
            # NOTE zero e1 or e2, e.g. an antiplane core: the same constants
            #      with delta e1 e2 = 1 - nu12 nu21 - ..., without dividing by
            #      the moduli, see the Notes
            den = 1-nu12*nu21-nu23*nu32-nu31*nu13-2*nu21*nu32*nu13
            if not den > 0:
                raise ValueError('MatLamina with a non-positive determinant '
                                 '1 - nu12*nu21 - nu23*nu32 - nu31*nu13 - '
                                 '2*nu21*nu32*nu13 = %r' % den)
            self.c11 = e1*(1    - nu23*nu32)/den
            self.c12 = e1*(nu21 + nu31*nu23)/den
            self.c13 = e1*(nu31 + nu21*nu32)/den
            self.c22 = e2*(1    - nu13*nu31)/den
            self.c23 = e2*(nu32 + nu12*nu31)/den
            self.c33 = e3*(1    - nu12*nu21)/den
        self.c44 = self.g23
        self.c55 = self.g13
        self.c66 = self.g12

        self.plane_stress = plane_stress
        if plane_stress:
            # NOTE plane-stress reduced stiffnesses, sigma_33 = 0, Jones
            #      (1999), Eq. (2.66), the same as in Lamina.rebuild; with
            #      e3 = 0 the same operations as the three-dimensional branch
            den = 1 - self.nu12 * self.nu21
            self.q11 = self.e1 / den
            self.q12 = self.e1*self.nu21 / den
            self.q13 = 0.
            self.q21 = self.e2*self.nu12 / den
            self.q22 = self.e2 / den
            self.q23 = 0.
            self.q31 = 0.
            self.q32 = 0.
            self.q33 = 0.
        else:
            # from references:
            #   hansen_hvejsen_2007 page 43
            #
            #   Guerdal Z., R. T. Haftka and P. Hajela (1999), Design and
            #   Optimization of Laminated Composite Materials,
            #   Wiley-Interscience.
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
            self.q31 = self.e3*(self.nu13 + self.nu12 * self.nu23) / den
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

    cpdef void trace_normalize_plane_stress(MatLamina self) except *:
        r"""Trace-normalize the lamina properties for plane stress

        Modify the original :class:`.MatLamina` object with a
        trace-normalization: the constitutive constants ``q11``, ``q12``,
        ``q22``, ``q44``, ``q55``, ``q66`` and the invariants ``u1`` to
        ``u7`` are divided by the trace ``tr = q11 + q22 + 2*q66``, such that
        they become dimensionless trace-normalized values. These can be
        accessed using the same attributes. The other attributes are not
        modified.

        Raises
        ------
        ValueError
            If ``plane_stress`` is ``False``, i.e. after
            ``rebuild(plane_stress=False)``, since the trace is defined with
            the plane-stress stiffnesses, or if the trace is not positive,
            e.g. zero for an antiplane core, whose in-plane stiffnesses
            vanish.

        Notes
        -----
        The trace is calculated according to Eq. (1), p. 45, of Melo et al.
        [1]_, `tr = Q_{11} + Q_{22} + 2Q_{66}`, the trace of the
        plane-stress stiffness matrix, with the plane-stress reduced
        stiffnesses of :meth:`rebuild`, also when ``e3`` is given, i.e. `tr
        = (e_1 + e_2)/(1 - \nu_{12}\nu_{21}) + 2 g_{12}`.  The universal
        in-plane stress stiffness components
        `Q_{11},Q_{12},Q_{22},Q_{44},Q_{55},Q_{66}` are divided by `tr`, and
        the invariants `U_1,U_2,U_3,U_4,U_5,U_6,U_7` are calculated with the
        normalized stiffnesses, such they also become trace-normalized
        invariants; since the invariants are linear in the stiffnesses, they
        are divided by `tr` as well. The trace does not depend on the ply
        angle, since `q_{11L} + q_{22L} + 2 q_{66L} = 2(u_1 + u_5)`, see
        :meth:`get_invariant_matrix` and Jones [2]_, Sec. 2.7, Eq. (2.98), p.
        85. Since the stiffnesses of
        :func:`.laminate_from_LaminationParameters` are linear in the
        invariants, those computed with the trace-normalized invariants are
        the ones of :func:`.laminated_plate` divided by `tr`.

        References
        ----------
        .. [1] Melo, J. D. D., Bi, J., and Tsai, S. W., A Novel
           Invariant-Based Design Approach to Carbon Fiber Reinforced
           Laminates, Composite Structures, 159, 44-52, 2017.
           https://doi.org/10.1016/j.compstruct.2016.09.055
        .. [2] Jones, R. M., Mechanics of Composite Materials, Second
           Edition, Taylor & Francis, Philadelphia, 1999.

        """
        cdef double tr
        if not self.plane_stress:
            raise ValueError('trace_normalize_plane_stress requires the '
                             'plane-stress stiffnesses, call '
                             'rebuild(plane_stress=True) first')
        tr = self.q11 + self.q22 + 2*self.q66
        if not tr > 0:
            raise ValueError('trace_normalize_plane_stress: the trace q11 + '
                             'q22 + 2*q66 = %r is not positive, e.g. that '
                             'of an antiplane core with zero in-plane '
                             'stiffnesses' % tr)
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

        Builds the symmetric matrix of the stiffness constants ``c11`` to
        ``c66`` computed by :meth:`rebuild`, with the rows and columns ordered
        as 11, 22, 33, 23, 13, 12 (``c44 = g23``, ``c55 = g13``, ``c66 =
        g12``).

        Returns
        -------
        C : double[:, ::1]
            Memoryview of shape ``(6, 6)``, force per unit area.

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

        Builds the matrix of the invariants ``u1`` to ``u7`` computed by
        :meth:`rebuild`. The rows correspond to ``q11``, ``q22``, ``q12``,
        ``q66``, ``q16``, ``q26``, ``q44``, ``q45``, ``q55`` and the columns to
        the factors of ``1`` and of the lamination parameters, in the same
        combinations used by :func:`.laminate_from_LaminationParameters`.

        Returns
        -------
        U : double[:, ::1]
            Memoryview of shape ``(9, 5)``, force per unit area.

        Notes
        -----
        The stiffnesses of a ply at the angle `\theta` in the laminate axes,
        see :meth:`.Lamina.rebuild`, are linear in `\cos 2\theta`, `\sin
        2\theta`, `\cos 4\theta` and `\sin 4\theta`, with the invariants of
        :meth:`rebuild` as coefficients [1]_:

        .. math::

            \{q_{11L}, q_{22L}, q_{12L}, q_{66L}, q_{16L}, q_{26L}, q_{44L},
            q_{45L}, q_{55L}\}^T = U \{1, \cos 2\theta, \sin 2\theta, \cos
            4\theta, \sin 4\theta\}^T

        where `q_{ijL}` are the plane-stress stiffnesses of
        :meth:`.Lamina.rebuild` when ``u1`` to ``u5`` are computed from the
        plane-stress reduced stiffnesses, i.e. with ``plane_stress=True``,
        the default of :meth:`rebuild`, for any ``e3`` (Jones [1]_, Sec.
        2.7, Eq. (2.98), p. 85), and the components of the
        three-dimensional stiffness matrix rotated about the 3 axis after
        ``rebuild(plane_stress=False)``. Integrated through the thickness,
        the trigonometric terms give the lamination parameters, see
        :func:`.laminate_from_LaminationParameters`.

        References
        ----------
        .. [1] Jones, R. M., Mechanics of Composite Materials, Second
           Edition, Taylor & Francis, Philadelphia, 1999.

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
    r"""Ply of a laminate

    Created with all attributes equal to zero and ``matlamina`` equal to
    ``None``. After ``matlamina``, ``h`` and ``thetadeg`` are set,
    :meth:`rebuild` computes the trigonometric terms of the ply angle and the
    ply stiffnesses ``q11L`` to ``q55L`` in the laminate axes. Usually created
    by :func:`.laminated_plate`.

    Attributes
    ----------

    plyid : int
        Identificaiton of the composite lamina
    matlamina : :class:`.MatLamina` object
        A :class:`.MatLamina` object
    h : float
        Ply thickness, length
    thetadeg : float
        Ply angle in degrees
    cost, sint : float
        Cosine and sine of the ply angle, dimensionless, computed by
        :meth:`rebuild`
    cos2t, sin2t, cos4t, sin4t : float
        Cosine and sine of two and four times the ply angle, dimensionless,
        computed by :meth:`rebuild`
    q11L, q12L, q16L, q22L, q26L, q66L : float
        Plane-stress stiffnesses of the ply in the laminate axes, force per
        unit area, computed by :meth:`rebuild`
    q44L, q45L, q55L : float
        Transverse shear stiffnesses of the ply in the laminate axes, force
        per unit area, computed by :meth:`rebuild`

    """
    def __init__(Lamina self):
        pass

    cpdef void rebuild(Lamina self):
        r"""Update constitutive matrices

        Computes ``cost``, ``sint``, ``cos2t``, ``sin2t``, ``cos4t`` and
        ``sin4t`` from ``thetadeg``, and the ply stiffnesses ``q11L`` to
        ``q55L`` in the laminate axes from the plane-stress stiffnesses of
        ``matlamina`` in the material axes, which are computed from its
        ``e1``, ``e2``, ``nu12``, ``nu21``, ``g12``, ``g13`` and ``g23``.

        Notes
        -----
        The plane-stress reduced stiffnesses in the material axes are

        .. math::

            Q_{11} = \frac{e_1}{1 - \nu_{12}\nu_{21}}, \quad
            Q_{12} = \frac{\nu_{12} e_2}{1 - \nu_{12}\nu_{21}}, \quad
            Q_{22} = \frac{e_2}{1 - \nu_{12}\nu_{21}}, \quad
            Q_{66} = g_{12}, \quad Q_{44} = g_{23}, \quad Q_{55} = g_{13}

        such that ``e3``, ``nu13`` and ``nu23`` are not used, and thermal
        expansion is not considered. With `c = \cos\theta` and `s =
        \sin\theta`, `\theta` being the ply angle, the stiffnesses in the
        laminate axes are [1]_:

        .. math::

            q_{11L} = Q_{11} c^4 + 2 (Q_{12} + 2 Q_{66}) s^2 c^2 + Q_{22} s^4 \\
            q_{12L} = (Q_{11} + Q_{22} - 4 Q_{66}) s^2 c^2 + Q_{12} (s^4 + c^4) \\
            q_{22L} = Q_{11} s^4 + 2 (Q_{12} + 2 Q_{66}) s^2 c^2 + Q_{22} c^4 \\
            q_{16L} = (Q_{11} - Q_{12} - 2 Q_{66}) s c^3 + (Q_{12} - Q_{22} +
            2 Q_{66}) s^3 c \\
            q_{26L} = (Q_{11} - Q_{12} - 2 Q_{66}) s^3 c + (Q_{12} - Q_{22} +
            2 Q_{66}) s c^3 \\
            q_{66L} = (Q_{11} + Q_{22} - 2 Q_{12} - 2 Q_{66}) s^2 c^2 + Q_{66}
            (s^4 + c^4) \\
            q_{44L} = Q_{44} c^2 + Q_{55} s^2, \quad
            q_{45L} = (Q_{55} - Q_{44}) s c, \quad
            q_{55L} = Q_{55} c^2 + Q_{44} s^2

        References
        ----------
        .. [1] Reddy, J. N., Mechanics of Laminated Composite Plates and Shells,
           CRC Press, 2003. https://doi.org/10.1201/b12409

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

    cpdef double [:, ::1] get_transf_matrix_displ_to_lamina(Lamina self):
        r"""Return displacement transformation matrix from laminate to lamina

        Builds ``[[cost, sint, 0], [-sint, cost, 0], [0, 0, 1]]`` from the
        terms computed by :meth:`rebuild`, which maps the displacement
        components in the laminate axes to those in the material axes of the
        ply.

        Returns
        -------
        T : double[:, ::1]
            Memoryview of shape ``(3, 3)``, dimensionless.

        Notes
        -----
        With the ply angle `\theta` measured from the laminate `x` axis to
        the material direction 1, as in the stiffnesses of :meth:`rebuild`,
        the material direction 1 is `\{\cos\theta, \sin\theta, 0\}` in the
        laminate axes, such that `\{u_1, u_2, u_3\}^T = T \{u_x, u_y,
        u_z\}^T`. The matrix is orthogonal, and its transpose maps the
        material components back to the laminate axes.

        """
        return np.array([[ self.cost, self.sint, 0],
                         [-self.sint, self.cost, 0],
                         [   0,     0, 1]], dtype=DOUBLE)

    cpdef double [:, ::1] get_transf_matrix_displ_to_laminate(Lamina self):
        r"""Deprecated, use :meth:`get_transf_matrix_displ_to_lamina` instead

        Issues a :class:`DeprecationWarning` and returns the matrix of
        :meth:`get_transf_matrix_displ_to_lamina`, which maps the
        displacements from the laminate axes to the material axes, contrary
        to what the old name says.

        Returns
        -------
        T : double[:, ::1]
            Memoryview of shape ``(3, 3)``, dimensionless.

        """
        _warn_deprecated('Lamina.get_transf_matrix_displ_to_laminate',
                         'Lamina.get_transf_matrix_displ_to_lamina')
        return self.get_transf_matrix_displ_to_lamina()

    cpdef double [:, ::1] get_constitutive_matrix(Lamina self):
        r"""Return the constitutive matrix

        Builds the symmetric matrix of the ply stiffnesses in the laminate axes
        computed by :meth:`rebuild`, ``[[q11L, q12L, q16L, 0, 0], [q12L, q22L,
        q26L, 0, 0], [q16L, q26L, q66L, 0, 0], [0, 0, 0, q44L, q45L], [0, 0,
        0, q45L, q55L]]``.

        Returns
        -------
        C : double[:, ::1]
            Memoryview of shape ``(5, 5)``, force per unit area.

        """
        return np.array([[self.q11L, self.q12L, self.q16L,    0,    0],
                         [self.q12L, self.q22L, self.q26L,    0,    0],
                         [self.q16L, self.q26L, self.q66L,    0,    0],
                         [   0,    0,    0, self.q44L, self.q45L],
                         [   0,    0,    0, self.q45L, self.q55L]], dtype=DOUBLE)

    cpdef double [:, ::1] get_transf_matrix_stress_to_lamina(Lamina self):
        r"""Return stress transformation matrix from laminate to lamina

        Builds the matrix from the terms ``cost``, ``sint`` and ``sin2t``
        computed by :meth:`rebuild`, such that `\{\sigma_{11}, \sigma_{22},
        \sigma_{33}, \tau_{23}, \tau_{13}, \tau_{12}\}^T = T \{\sigma_{xx},
        \sigma_{yy}, \sigma_{zz}, \tau_{yz}, \tau_{xz}, \tau_{xy}\}^T`, with
        the material axes 1, 2, 3 and the laminate axes `x`, `y`, `z`.

        Returns
        -------
        T : double[:, ::1]
            Memoryview of shape ``(6, 6)``, dimensionless.

        See Also
        --------
        get_transf_matrix_stress_to_laminate

        Notes
        -----
        With the stress components ordered as `\{\sigma_{xx}, \sigma_{yy},
        \sigma_{zz}, \tau_{yz}, \tau_{xz}, \tau_{xy}\}`, the order of
        :meth:`.MatLamina.get_constitutive_matrix`, the matrix is the
        transformation of the stresses from the laminate axes to the material
        axes for a rotation by the ply angle `\theta` about `z`, e.g.
        `\sigma_{11} = \cos^2\theta \sigma_{xx} + \sin^2\theta \sigma_{yy} +
        \sin 2\theta \tau_{xy}` [1]_. Its inverse,
        :meth:`get_transf_matrix_stress_to_laminate`, is the same matrix for
        `-\theta`.

        References
        ----------
        .. [1] Reddy, J. N., Mechanics of Laminated Composite Plates and Shells,
           CRC Press, 2003. https://doi.org/10.1201/b12409

        """
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
        r"""Return stress transformation matrix from lamina to laminate

        Builds the matrix from the terms ``cost``, ``sint`` and ``sin2t``
        computed by :meth:`rebuild`, with the opposite sign of all the terms
        in ``sint`` and ``sin2t`` with respect to
        :meth:`get_transf_matrix_stress_to_lamina`, such that
        `\{\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \tau_{yz}, \tau_{xz},
        \tau_{xy}\}^T = T \{\sigma_{11}, \sigma_{22}, \sigma_{33}, \tau_{23},
        \tau_{13}, \tau_{12}\}^T`, with the material axes 1, 2, 3 and the
        laminate axes `x`, `y`, `z`.

        Returns
        -------
        T : double[:, ::1]
            Memoryview of shape ``(6, 6)``, dimensionless.

        See Also
        --------
        get_transf_matrix_stress_to_lamina

        Notes
        -----
        It is the inverse of :meth:`get_transf_matrix_stress_to_lamina`, i.e.
        the transformation of the stresses from the material axes to the
        laminate axes, obtained by replacing the ply angle `\theta` by
        `-\theta`.

        """
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


# NOTE methods accepted by Laminate.shear_correction, besides None
_SHEAR_CORRECTIONS = ('rohwer', 'vlachoutsis', 'whitney', 'chow',
                      'birman_bert', 'thickness_shear', 'constant')


def _check_symmetric(Laminate lam, double[::1] z, double h):
    r"""Raise if the laminate is not symmetric about its mid-surface

    Checks the extension-bending coupling about the mid-surface, computed
    from the plies with the interfaces ``z`` measured from the mid-surface,
    against `h \max|A_{ij}|`.

    Notes
    -----
    A laminate is symmetric about its mid-surface when its extension-bending
    coupling about the mid-surface vanishes, `B_{ij} = \sum_k q_{ijL}^{(k)}
    (z_{k+1}^2 - z_k^2)/2 = 0`, as assumed by the factors of ``'chow'``, see
    :meth:`.Laminate.calc_transverse_shear_stiffness`.

    """
    cdef int k, N
    cdef double za, zb
    cdef Lamina ply
    N = <int>len(lam.plies)
    Amax = 0.
    Bmid = np.zeros(6, dtype=DOUBLE)
    for k in range(N):
        ply = lam.plies[k]
        za = z[k]
        zb = z[k+1]
        q = np.array([ply.q11L, ply.q12L, ply.q16L, ply.q22L, ply.q26L,
                      ply.q66L])
        Bmid += q*(zb*zb - za*za)/2.
        Amax = max(Amax, np.abs(q).max()*(zb - za))
    if np.abs(Bmid).max() > 1e-10*h*N*Amax:
        raise ValueError("'chow' shear correction factors require a laminate "
                         "symmetric about its mid-surface, B_ij about the "
                         "mid-surface = %s" % np.array2string(Bmid,
                         precision=4))


def _thickness_shear_frequency(Gs, rhos, hs, double omega_ref):
    r"""Lowest non-zero thickness-shear frequency of a laminate

    Returns the frequency of the first thickness-shear mode of an infinite
    laminate with traction-free faces, used by the ``'thickness_shear'``
    method of :meth:`.Laminate.calc_transverse_shear_stiffness`. The
    frequency is searched from ``1e-4*omega_ref`` up to ``100*omega_ref``,
    and refined by bisection to a relative interval of ``1e-15`` or for at
    most 200 bisections.

    Parameters
    ----------
    Gs, rhos, hs : array-like
        Transverse shear modulus, density and thickness of each ply, from
        bottom to top.
    omega_ref : float
        Reference frequency that scales the search, in radians per unit
        time.

    Returns
    -------
    omega : float
        The lowest non-zero thickness-shear frequency.

    Raises
    ------
    RuntimeError
        If no frequency is found below ``100*omega_ref``.

    Notes
    -----
    Solves the anti-plane motion `u(z) e^{i \omega t}` of an infinite
    laminate, `(G u_{,z})_{,z} + \rho \omega^2 u = 0` in each ply, with
    traction-free faces and continuity of `u` and `\tau = G u_{,z}` at the
    interfaces. Starting from `\{u, \tau\} = \{1, 0\}` at the bottom face, the
    exact transfer matrix of each ply gives `\tau(\omega)` at the top face,
    whose lowest positive root is returned. With `\beta = \omega
    \sqrt{\rho/G}`, the transfer matrix of a ply of thickness `h_k` is

    .. math::

        \begin{Bmatrix} u \\ \tau \end{Bmatrix}_{top} =
        \begin{bmatrix} \cos \beta h_k & \sin(\beta h_k)/(G \beta) \\
        -G \beta \sin \beta h_k & \cos \beta h_k \end{bmatrix}
        \begin{Bmatrix} u \\ \tau \end{Bmatrix}_{bottom}

    Since `\tau \approx -\omega^2 I_0` for small `\omega`, the root is
    bracketed by the first change of sign found with geometric steps from
    ``1e-4*omega_ref``, each one multiplying the frequency by 1.01, and
    refined by bisection. ``omega_ref`` is `\sqrt{\bar{A}_{55}/\bar{I}_2}`,
    i.e. the FSDT frequency with a unit shear correction factor.

    """
    Gs = np.asarray(Gs, dtype=DOUBLE)
    rhos = np.asarray(rhos, dtype=DOUBLE)
    hs = np.asarray(hs, dtype=DOUBLE)
    slowness = np.sqrt(rhos/Gs)

    def tau_top(double omega):
        cdef double u = 1., tau = 0., beta, c, s, un
        cdef int k
        for k in range(Gs.shape[0]):
            beta = omega*slowness[k]
            c = cos(beta*hs[k])
            s = sin(beta*hs[k])
            un = c*u + s/(Gs[k]*beta)*tau
            tau = -Gs[k]*beta*s*u + c*tau
            u = un
        return tau

    cdef double wa, wb, fa, fb, wm, fm
    wa = 1e-4*omega_ref
    fa = tau_top(wa)
    wb = wa
    fb = fa
    while fb < 0:
        if wb > 1e2*omega_ref:
            raise RuntimeError('No thickness-shear frequency found below '
                               '100 times the reference frequency')
        wa, fa = wb, fb
        wb = 1.01*wa
        fb = tau_top(wb)
    for _ in range(200):
        wm = (wa + wb)/2.
        if wb - wa <= 1e-15*wm:
            break
        fm = tau_top(wm)
        if fm < 0:
            wa, fa = wm, fm
        else:
            wb, fb = wm, fm
    return (wa + wb)/2.


def _equilibrium_integrals(Lamina ply, eps_x, eps_y, double za, double zb):
    r"""`\int_{z_a}^{z_b}` of the in-plane stress gradients of one ply

    Returns the integrands of `-\tau_{xz}` and `-\tau_{yz}` of
    :meth:`.Laminate.calc_equilibrium_transverse_shear`, for the derivatives
    ``eps_x``, ``eps_y`` of the generalized strains, with shape
    ``(npts, 6)``.

    Notes
    -----
    With `c_1 = \{q_{11L}, q_{12L}, q_{16L}\}`, `c_2 = \{q_{12L}, q_{22L},
    q_{26L}\}` and `c_3 = \{q_{16L}, q_{26L}, q_{66L}\}` of the ply, and the
    derivatives `\{\varepsilon^{(0)}, \varepsilon^{(1)}\}_{,x}` and
    `\{\varepsilon^{(0)}, \varepsilon^{(1)}\}_{,y}` split into their first
    and last three entries, `\varepsilon^{(0)}_{,x}`,
    `\varepsilon^{(1)}_{,x}` and `\varepsilon^{(0)}_{,y}`,
    `\varepsilon^{(1)}_{,y}`, the integrals are computed exactly:

    .. math::

        i_x = \int_{z_a}^{z_b} \left( c_1 (\varepsilon^{(0)}_{,x} + z
        \varepsilon^{(1)}_{,x}) + c_3 (\varepsilon^{(0)}_{,y} + z
        \varepsilon^{(1)}_{,y}) \right) dz \\
        i_y = \int_{z_a}^{z_b} \left( c_2 (\varepsilon^{(0)}_{,y} + z
        \varepsilon^{(1)}_{,y}) + c_3 (\varepsilon^{(0)}_{,x} + z
        \varepsilon^{(1)}_{,x}) \right) dz

    """
    c1 = np.array([ply.q11L, ply.q12L, ply.q16L])
    c2 = np.array([ply.q12L, ply.q22L, ply.q26L])
    c3 = np.array([ply.q16L, ply.q26L, ply.q66L])
    d1 = zb - za
    d2 = (zb*zb - za*za)/2.
    ix = (d1*(eps_x[:, :3] @ c1) + d2*(eps_x[:, 3:] @ c1)
          + d1*(eps_y[:, :3] @ c3) + d2*(eps_y[:, 3:] @ c3))
    iy = (d1*(eps_y[:, :3] @ c2) + d2*(eps_y[:, 3:] @ c2)
          + d1*(eps_x[:, :3] @ c3) + d2*(eps_x[:, 3:] @ c3))
    return ix, iy


def _strain_derivatives(name, eps_x, eps_y, grad_x, grad_y):
    # NOTE grad_x and grad_y are the deprecated names of eps_x and eps_y
    if grad_x is not None:
        if eps_x is not None:
            raise TypeError("%s() got both 'eps_x' and its deprecated alias "
                            "'grad_x'" % name)
        _warn_deprecated('grad_x', 'eps_x')
        eps_x = grad_x
    if grad_y is not None:
        if eps_y is not None:
            raise TypeError("%s() got both 'eps_y' and its deprecated alias "
                            "'grad_y'" % name)
        _warn_deprecated('grad_y', 'eps_y')
        eps_y = grad_y
    if eps_x is None or eps_y is None:
        raise TypeError("%s() missing the argument '%s'"
                        % (name, 'eps_x' if eps_x is None else 'eps_y'))
    return eps_x, eps_y


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
       'intrho', 'intrhoz', 'intrhoz2', 'intrhoz3', 'intrhoz4', 'intrhoz6',
       'plies', 'stack', 'shear_correction'])


cdef class Laminate:
    r"""Laminate, with its stiffnesses and mass integrals

    Created with ``offset = 0``, ``plies = []``, ``stack = []``,
    ``shear_correction = 'rohwer'``, ``scf_k13 = scf_k23 = 5/6`` and all other
    attributes equal to zero. After ``plies``, a list of :class:`.Lamina`
    objects from the bottom to the top face, and ``offset`` are set,
    :meth:`calc_constitutive_matrix` computes ``h``, the stiffnesses and the
    mass integrals, and :meth:`calc_equivalent_properties` the equivalent
    moduli. The stiffnesses are thickness integrals of the ply stiffnesses
    ``qijL`` of :class:`.Lamina`, with `z` measured from the reference
    surface. Usually created by :func:`.laminated_plate`,
    :func:`.isotropic_plate` or :func:`.laminate_from_LaminationParameters`.

    The object can be pickled: the public attributes and the instance
    dictionary are saved, and the through-thickness transverse shear
    distribution used by :meth:`calc_transverse_shear_stress` is recomputed
    when needed.

    Attributes
    ----------

    plies : list
        List of plies
    stack : list
        List of angles for each ply
    h : float
        Total thickness of the laminate, length
    offset : float
        Offset at the normal direction, length
    e1, e2 : float
        Equivalent laminate moduli in directions 1 and 2, force per unit area
    g12 : float
        Equivalent laminate shear modulus in the 12 direction, force per unit
        area
    nu12, nu21 : float
        Equivalent laminate Poisson ratios in the 12 and 21 directions,
        dimensionless
    A11, A12, A16, A22, A26, A66 : float
        Membrane stiffnesses `\int q_{ijL} dz`, force per unit length
    B11, B12, B16, B22, B26, B66 : float
        Membrane-bending coupling stiffnesses `\int q_{ijL} z dz`, force
    D11, D12, D16, D22, D26, D66 : float
        Bending stiffnesses `\int q_{ijL} z^2 dz`, force times length
    E11, E12, E16, E22, E26, E66 : float
        Stiffnesses `\int q_{ijL} z^3 dz`, force times length squared
    F11, F12, F16, F22, F26, F66 : float
        Stiffnesses `\int q_{ijL} z^4 dz`, force times length cubed
    H11, H12, H16, H22, H26, H66 : float
        Stiffnesses `\int q_{ijL} z^6 dz`, force times length to the fifth
        power
    D44, D45, D55 : float
        Transverse shear stiffnesses `\int q_{ijL} z^2 dz`, force times
        length, see ``Dts``
    F44, F45, F55 : float
        Transverse shear stiffnesses `\int q_{ijL} z^4 dz`, force times
        length cubed, see ``Fts``
    A44, A45, A55 : float
        Transverse shear stiffnesses of the first-order shear deformation
        theory (FSDT), **with the shear correction already applied**
        according to ``shear_correction``. They are ready to be used in
        the element and no shear correction factor should be applied to
        them downstream. See :meth:`.calc_transverse_shear_stiffness`.
        Force per unit length.
    Abar44, Abar45, Abar55 : float
        Constant-strain (uncorrected) transverse shear stiffnesses
        `\bar{A}_{ts} = \sum_k C_s^{(k)} h_k`. These are the terms to be used
        together with ``Dts`` and ``Fts`` in the third-order shear
        deformation theory (TSDT). Force per unit length.
    Abarbar44, Abarbar45, Abarbar55 : float
        Constant-stress transverse shear stiffnesses `\bar{\bar{A}}_{ts} = h^2
        [\sum_k (C_s^{(k)})^{-1} h_k]^{-1}`, for comparison only. Equal to
        ``nan`` when a ply has a singular `C_s^{(k)}`. Force per unit length.
    shear_correction : str or None
        Method used to obtain ``A44``, ``A45``, ``A55`` from the ply data:
        ``'rohwer'``, ``'vlachoutsis'``, ``'whitney'``, ``'chow'``,
        ``'birman_bert'``, ``'thickness_shear'``, ``'constant'`` or ``None``,
        see :meth:`.calc_transverse_shear_stiffness`. Default is
        ``'rohwer'``.
    scf_k13, scf_k23 : float
        Reported shear correction ratios ``A55/Abar55`` and ``A44/Abar44``.
        They are informative only, the correction is already inside ``A44``,
        ``A45``, ``A55``. Dimensionless.
    intrho : float
        Integral `\int_{-h/2+offset}^{+h/2+offset} \rho(z) dz`, used in
        equivalent single layer finite element mass matrices, mass per unit
        area
    intrhoz : float
        Integral `\int_{-h/2+offset}^{+h/2+offset} \rho(z)z dz`, used in
        equivalent single layer finite element mass matrices, mass per unit
        area times length
    intrhoz2 : float
        Integral `\int_{-h/2+offset}^{+h/2+offset} \rho(z)z^2 dz`, used in
        equivalent single layer finite element mass matrices, mass per unit
        area times length squared
    intrhoz3 : float
        Integral `\int_{-h/2+offset}^{+h/2+offset} \rho(z)z^3 dz`, used in
        the mass matrices of the third-order shear deformation theory (TSDT),
        mass per unit area times length cubed
    intrhoz4 : float
        Integral `\int_{-h/2+offset}^{+h/2+offset} \rho(z)z^4 dz`, used in
        the mass matrices of the TSDT, mass per unit area times length to the
        fourth power
    intrhoz6 : float
        Integral `\int_{-h/2+offset}^{+h/2+offset} \rho(z)z^6 dz`, used in
        the mass matrices of the TSDT, mass per unit area times length to the
        sixth power

    Notes
    -----
    In the notation of Reddy [1]_, the stiffnesses relate the stress
    resultants to the generalized strains of the reference surface,
    `\varepsilon^{(0)} = \{\varepsilon_{xx}^{(0)}, \varepsilon_{yy}^{(0)},
    \gamma_{xy}^{(0)}\}` and `\varepsilon^{(1)} = \{\varepsilon_{xx}^{(1)},
    \varepsilon_{yy}^{(1)}, \gamma_{xy}^{(1)}\}`, the in-plane strains being
    `\varepsilon^{(0)} + z \varepsilon^{(1)}`. In the CLPT and in the FSDT:

    .. math::

        \begin{Bmatrix} N \\ M \end{Bmatrix} =
        \begin{bmatrix} A & B \\ B & D \end{bmatrix}
        \begin{Bmatrix} \varepsilon^{(0)} \\ \varepsilon^{(1)}
        \end{Bmatrix}, \qquad
        \begin{Bmatrix} Q_y \\ Q_x \end{Bmatrix} = A_{ts}
        \begin{Bmatrix} \gamma_{yz} \\ \gamma_{xz} \end{Bmatrix}

    with `N = \{N_{xx}, N_{yy}, N_{xy}\}` and `M = \{M_{xx}, M_{yy},
    M_{xy}\}`, the last relation of the FSDT only, where the constant
    transverse shear strains through the thickness require the shear
    correction contained in ``A44``, ``A45``, ``A55``. The TSDT uses as well
    ``E``, ``F``, ``H``, ``Dts`` and ``Fts``, with the uncorrected
    ``Abar44``, ``Abar45``, ``Abar55``, since the TSDT needs no shear
    correction. The stiffnesses are computed by
    :meth:`calc_constitutive_matrix`. With ``offset`` different than zero,
    the reference surface is not the mid-surface, and ``B`` is not zero
    even for an isotropic plate.

    References
    ----------
    .. [1] Reddy, J. N., Mechanics of Laminated Composite Plates and Shells,
       CRC Press, 2003. https://doi.org/10.1201/b12409

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
        self.intrhoz3 = 0.
        self.intrhoz4 = 0.
        self.intrhoz6 = 0.
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
        # NOTE deprecated, kept for cimporting packages, same as get_Ats
        _warn_deprecated('Laminate.get_Atrans', 'Laminate.get_Ats')
        return self.get_Ats()
    cdef double [:, ::1] get_Abar_ts(Laminate self):
        return np.array([[self.Abar44, self.Abar45],
                         [self.Abar45, self.Abar55]], dtype=DOUBLE)
    cdef double [:, ::1] get_Abarbar_ts(Laminate self):
        return np.array([[self.Abarbar44, self.Abarbar45],
                         [self.Abarbar45, self.Abarbar55]], dtype=DOUBLE)
    cdef double [:, ::1] get_Dts(Laminate self):
        return np.array([[self.D44, self.D45],
                         [self.D45, self.D55]], dtype=DOUBLE)
    cdef double [:, ::1] get_Fts(Laminate self):
        return np.array([[self.F44, self.F45],
                         [self.F45, self.F55]], dtype=DOUBLE)
    cdef double [:, ::1] get_Dtrans(Laminate self):
        # NOTE deprecated, kept for cimporting packages, same as get_Dts
        _warn_deprecated('Laminate.get_Dtrans', 'Laminate.get_Dts')
        return self.get_Dts()
    cdef double [:, ::1] get_Ftrans(Laminate self):
        # NOTE deprecated, kept for cimporting packages, same as get_Fts
        _warn_deprecated('Laminate.get_Ftrans', 'Laminate.get_Fts')
        return self.get_Fts()
    cdef double [:, ::1] get_ABD(Laminate self):
        return np.array([[self.A11, self.A12, self.A16, self.B11, self.B12, self.B16],
                         [self.A12, self.A22, self.A26, self.B12, self.B22, self.B26],
                         [self.A16, self.A26, self.A66, self.B16, self.B26, self.B66],
                         [self.B11, self.B12, self.B16, self.D11, self.D12, self.D16],
                         [self.B12, self.B22, self.B26, self.D12, self.D22, self.D26],
                         [self.B16, self.B26, self.B66, self.D16, self.D26, self.D66]], dtype=DOUBLE)
    @property
    def A(self):
        r"""Membrane stiffness matrix

        Built from the attributes ``A11`` to ``A66``; a new array is returned
        at each access, such that changing it does not change the attributes.

        Returns
        -------
        A : ndarray
            Array of shape ``(3, 3)``, force per unit length.

        """
        return np.asarray(self.get_A())
    @property
    def B(self):
        r"""Membrane-bending coupling stiffness matrix

        Built from the attributes ``B11`` to ``B66``; a new array is returned
        at each access, such that changing it does not change the attributes.

        Returns
        -------
        B : ndarray
            Array of shape ``(3, 3)``, force.

        """
        return np.asarray(self.get_B())
    @property
    def D(self):
        r"""Bending stiffness matrix

        Built from the attributes ``D11`` to ``D66``; a new array is returned
        at each access, such that changing it does not change the attributes.

        Returns
        -------
        D : ndarray
            Array of shape ``(3, 3)``, force times length.

        """
        return np.asarray(self.get_D())
    @property
    def E(self):
        r"""Stiffness matrix of the third-order shear deformation theory

        Built from the attributes ``E11`` to ``E66``; a new array is returned
        at each access, such that changing it does not change the attributes.

        Returns
        -------
        E : ndarray
            Array of shape ``(3, 3)``, force times length squared.

        """
        return np.asarray(self.get_E())
    @property
    def F(self):
        r"""Stiffness matrix of the third-order shear deformation theory

        Built from the attributes ``F11`` to ``F66``; a new array is returned
        at each access, such that changing it does not change the attributes.

        Returns
        -------
        F : ndarray
            Array of shape ``(3, 3)``, force times length cubed.

        """
        return np.asarray(self.get_F())
    @property
    def H(self):
        r"""Stiffness matrix of the third-order shear deformation theory

        Built from the attributes ``H11`` to ``H66``; a new array is returned
        at each access, such that changing it does not change the attributes.

        Returns
        -------
        H : ndarray
            Array of shape ``(3, 3)``, force times length to the fifth power.

        """
        return np.asarray(self.get_H())
    @property
    def Ats(self):
        r"""Transverse shear stiffness matrix ``[[A44, A45], [A45, A55]]``

        Index 4 corresponds to `yz` and index 5 to `xz`, such that `\{Q_y,
        Q_x\}^T = A_{ts} \{\gamma_{yz}, \gamma_{xz}\}^T`. The shear correction
        is already applied, see :meth:`.calc_transverse_shear_stiffness`.
        A new array is returned at each access, such that changing it does not
        change the attributes.

        Returns
        -------
        Ats : ndarray
            Array of shape ``(2, 2)``, force per unit length.

        """
        return np.asarray(self.get_Ats())
    @property
    def Atrans(self):
        r"""Deprecated, use :attr:`.Ats` instead

        Returns the same corrected matrix as :attr:`.Ats`, issuing a
        :class:`DeprecationWarning`.

        Returns
        -------
        Atrans : ndarray
            Array of shape ``(2, 2)``, force per unit length.

        """
        warnings.warn("'Laminate.Atrans' is deprecated, use 'Laminate.Ats' "
                      "instead, which contains the shear correction",
                      DeprecationWarning, stacklevel=1)
        return np.asarray(self.get_Ats())
    @property
    def Abar_ts(self):
        r"""Constant-strain ``[[Abar44, Abar45], [Abar45, Abar55]]``

        Uncorrected transverse shear stiffness, to be used with ``Dts`` and
        ``Fts`` in the third-order shear deformation theory (TSDT).
        A new array is returned at each access, such that changing it does not
        change the attributes.

        Returns
        -------
        Abar_ts : ndarray
            Array of shape ``(2, 2)``, force per unit length.

        """
        return np.asarray(self.get_Abar_ts())
    @property
    def Abarbar_ts(self):
        r"""Constant-stress ``[[Abarbar44, Abarbar45], [Abarbar45, Abarbar55]]``

        Built from the attributes ``Abarbar44``, ``Abarbar45``, ``Abarbar55``.
        A new array is returned at each access, such that changing it does not
        change the attributes.

        Returns
        -------
        Abarbar_ts : ndarray
            Array of shape ``(2, 2)``, force per unit length.

        """
        return np.asarray(self.get_Abarbar_ts())
    @property
    def Dts(self):
        r"""TSDT transverse shear stiffness ``[[D44, D45], [D45, D55]]``

        Built from the attributes ``D44``, ``D45``, ``D55``.
        A new array is returned at each access, such that changing it does not
        change the attributes.

        Returns
        -------
        Dts : ndarray
            Array of shape ``(2, 2)``, force times length.

        """
        return np.asarray(self.get_Dts())
    @property
    def Fts(self):
        r"""TSDT transverse shear stiffness ``[[F44, F45], [F45, F55]]``

        Built from the attributes ``F44``, ``F45``, ``F55``.
        A new array is returned at each access, such that changing it does not
        change the attributes.

        Returns
        -------
        Fts : ndarray
            Array of shape ``(2, 2)``, force times length cubed.

        """
        return np.asarray(self.get_Fts())
    @property
    def Dtrans(self):
        r"""Deprecated, use :attr:`.Dts` instead

        Returns the same matrix as :attr:`.Dts`, issuing a
        :class:`DeprecationWarning`.

        Returns
        -------
        Dtrans : ndarray
            Array of shape ``(2, 2)``, force times length.

        """
        _warn_deprecated('Laminate.Dtrans', 'Laminate.Dts')
        return np.asarray(self.get_Dts())
    @property
    def Ftrans(self):
        r"""Deprecated, use :attr:`.Fts` instead

        Returns the same matrix as :attr:`.Fts`, issuing a
        :class:`DeprecationWarning`.

        Returns
        -------
        Ftrans : ndarray
            Array of shape ``(2, 2)``, force times length cubed.

        """
        _warn_deprecated('Laminate.Ftrans', 'Laminate.Fts')
        return np.asarray(self.get_Fts())
    @property
    def ABD(self):
        r"""Stiffness matrix ``[[A, B], [B, D]]``

        Built from the attributes ``A11`` to ``D66``; a new array is returned
        at each access, such that changing it does not change the attributes.

        Returns
        -------
        ABD : ndarray
            Array of shape ``(6, 6)``, with the dimensions of ``A``, ``B`` and
            ``D``.

        """
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
                      stacklevel=1)
        if self.shear_correction is None:
            self.shear_correction = 'rohwer'
        self.calc_transverse_shear_stiffness()
        return self.scf_k13, self.scf_k23


    cdef int _calc_transverse_shear_distribution(Laminate self) except -1:
        r"""Coefficients of the distribution matrix `f^{(k)}(z)` of ``'rohwer'``

        For each ply `k`, stores the 2x2 matrices `F_0`, `F_1`, `F_2` in
        ``_ts_fcoef[k, 0:3]``, such that `f^{(k)}(\bar{z}) = F_0 + \bar{z} F_1
        + \bar{z}^2 F_2` and `\{\tau_{yz}, \tau_{xz}\}^T = f^{(k)}(\bar{z})
        \{Q_y, Q_x\}^T`, with `\bar{z} = z - offset` measured from the
        mid-surface. The ply interfaces, also measured from the mid-surface,
        are stored in ``_ts_z``.

        The ABD matrix is recomputed from the plies, such that the
        distribution is consistent with them even if the laminate stiffness
        attributes were modified afterwards, e.g. by :meth:`.make_symmetric`.

        Notes
        -----
        The distribution does not depend on the reference surface, and
        measuring from the mid-surface avoids the round-off of the ABD terms
        and of the polynomials for large values of ``offset``.

        The equilibrium approach of Rohwer (1988) [1]_ assumes two
        cylindrical bending states with constant membrane forces: one in
        which the resultants depend on `x` only, with `M_{xx,x} = Q_x` and
        `M_{xy,x} = Q_y`, and one in which they depend on `y` only, with
        `M_{yy,y} = Q_y` and `M_{xy,y} = Q_x`. `\tau_{xz}` is obtained from
        the first state and `\tau_{yz}` from the second one, integrating the
        equilibrium equations `\tau_{xz,z} = -\sigma_{xx,x}` and `\tau_{yz,z}
        = -\sigma_{yy,y}` from the bottom face `z_1`, where the tractions
        vanish. With `H^* = [[A, B], [B, D]]^{-1}` about the mid-surface, the
        in-plane strains of a state are `(H^*_{1:3} + z H^*_{4:6}) \{N,
        M\}`, where `H^*_{1:3}` and `H^*_{4:6}` are the first and the last
        three rows of `H^*` and `\{N, M\} = \{N_{xx}, N_{yy}, N_{xy}, M_{xx},
        M_{yy}, M_{xy}\}`. With `c_1^{(k)} = \{q_{11L}, q_{12L}, q_{16L}\}`
        and `c_2^{(k)} = \{q_{12L}, q_{22L}, q_{26L}\}` the rows of ply `k`
        that give `\sigma_{xx}` and `\sigma_{yy}`, and `\Pi(z) = z H^*_{1:3}
        + z^2/2 H^*_{4:6}`:

        .. math::

            \tau_{xz}(z) = \left( a_1^{(k)} - c_1^{(k)} \Pi(z) \right) L_x
            \begin{Bmatrix} Q_y \\ Q_x \end{Bmatrix}, \qquad
            \tau_{yz}(z) = \left( a_2^{(k)} - c_2^{(k)} \Pi(z) \right) L_y
            \begin{Bmatrix} Q_y \\ Q_x \end{Bmatrix} \\
            a_\alpha^{(k)} = \sum_{i=1}^{k} \left( c_\alpha^{(i)} -
            c_\alpha^{(i-1)} \right) \Pi(z_i), \qquad c_\alpha^{(0)} = 0

        for `z_k \le z \le z_{k+1}`, where `L_x` places `Q_y` and `Q_x` at
        the components `M_{xy}` and `M_{xx}` of `\{N, M\}_{,x}`, and `L_y`
        places them at `M_{yy}` and `M_{xy}` of `\{N, M\}_{,y}`. The sum
        `a_\alpha^{(k)}` makes the stresses continuous at the interfaces.
        Within each ply `f^{(k)}` is quadratic in `z`, and it vanishes at the
        top face, since `\int c_1 (H^*_{1:3} + z H^*_{4:6}) dz` through the
        thickness is the row of the identity of `N_{xx}`, analogously for
        `c_2` and `N_{yy}`, which `L_x` and `L_y` do not select.

        References
        ----------
        .. [1] Rohwer, K., Improved transverse shear stiffnesses for layered
           finite elements, DFVLR-FB 88-32, 1988.

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
        z[0] = -h/2.

        # ABD = [[A, B], [B, D]] about the mid-surface, from the plies
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
        #NOTE the condition number is tested with isfinite too, since it
        #     may be inf or nan for a singular matrix, and the inverse as
        #     well, since np.linalg.inv does not raise LinAlgError for a
        #     singular matrix on every platform, e.g. Pyodide
        cond = np.linalg.cond(M) if np.all(np.isfinite(M)) else np.inf
        if not np.isfinite(cond) or cond > 1e12:
            raise ValueError('The ABD matrix of the laminate is singular or '
                             'ill-conditioned, the transverse shear '
                             'distribution cannot be computed')
        Hs = np.ascontiguousarray(scale[:, None]*np.linalg.inv(M)*scale[None, :])
        if not np.all(np.isfinite(Hs)):
            raise ValueError('The ABD matrix of the laminate is singular or '
                             'ill-conditioned, the transverse shear '
                             'distribution cannot be computed')

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
        None of the methods depends on ``offset``.

        The method is selected by the attribute ``shear_correction``:

        - ``'rohwer'`` (default): equilibrium approach, giving the full 2x2
          matrix.

        - ``'vlachoutsis'``: the scalar factors ``k13``, ``k23`` are applied
          to the constant-strain stiffness, ``A55 = k13*Abar55``, ``A44 =
          k23*Abar44``, and ``A45 = (k13 + k23)/2*Abar45``.

        - ``'whitney'``: the result is identical to ``'vlachoutsis'``.

        - ``'chow'``: for laminates symmetric about the mid-surface only,
          for which the result is identical to ``'vlachoutsis'``.

        - ``'birman_bert'``: average shear strain factors, applied as in
          ``'vlachoutsis'``. It gives `k = 1` for a homogeneous plate.

        - ``'thickness_shear'``: dynamic factors, applied as in
          ``'vlachoutsis'``. It gives `k = \pi^2/12` for a homogeneous plate.
          The ply densities must be positive.

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

        Without plies, ``A44``, ``A45``, ``A55``, ``Abar44``, ``Abar45``,
        ``Abar55``, ``Abarbar44``, ``Abarbar45`` and ``Abarbar55`` are set to
        zero and nothing else is updated.

        Raises
        ------
        ValueError
            If ``shear_correction`` is not recognized; if a ply has a
            singular transverse shear constitutive matrix `C_s^{(k)}`, e.g.
            ``g13 = 0`` or ``g23 = 0`` (all methods except ``'constant'``
            and ``None``); if the ABD matrix of the laminate is singular
            (``'rohwer'``); if the laminate is not symmetric (``'chow'``); or
            if a ply has a non-positive density (``'thickness_shear'``). Also
            if the transverse shear compliance obtained with ``'rohwer'`` is
            singular, or if the sum of ``q11L`` (direction 1) or ``q22L``
            (direction 2) times the ply thicknesses is not positive
            (``'vlachoutsis'``, ``'whitney'``, ``'chow'``, ``'birman_bert'``).

        Notes
        -----
        All of the methods are evaluated with `z` measured from the
        mid-surface, which avoids the round-off of large values of
        ``offset``. The constant-strain stiffness `\bar{A}_{ts} = \sum_k
        C_s^{(k)} h_k` corresponds to transverse shear strains constant
        through the thickness, and the constant-stress stiffness
        `\bar{\bar{A}}_{ts} = h^2 [\sum_k (C_s^{(k)})^{-1} h_k]^{-1}` to the
        constant stresses `\{\tau_{yz}, \tau_{xz}\} = \{Q_y, Q_x\}/h`, with
        the same complementary energy.

        - ``'rohwer'``: equilibrium approach of Rohwer (1988) [1]_. The
          transverse shear stresses are obtained from the equilibrium of two
          cylindrical bending states, with zero tractions at the bottom and
          top faces and continuity at every interface, `\{\tau_{yz},
          \tau_{xz}\}^T = f^{(k)}(z) \{Q_y, Q_x\}^T`, see
          :meth:`._calc_transverse_shear_distribution`. The 2x2 stiffness is
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
          inherent to the method, Rohwer (1988) [1]_, Fig. 7, p. 20, and Eq.
          (44): for his `[0, 90]_s` CFRP laminate of Table 2, computed in axes
          rotated by 45 degrees, ``A44 = A55`` are 7.0 % and ``A45`` is
          13.5 % above the tensor rotation of the result in the axes of the
          plies; ``A44`` or ``A55`` alone deviate by up to 7.2 %, at 41 and
          49 degrees.

        - ``'vlachoutsis'``: the scalar factors `k_{13}`, `k_{23}` of
          Vlachoutsis (1992) [2]_ are applied to the constant-strain
          stiffness, ``A55 = k13*Abar55``, ``A44 = k23*Abar44``, and the
          ad-hoc ``A45 = (k13 + k23)/2*Abar45``, which cannot represent the
          coupling of angle-ply laminates with ``Abar45 = 0``. The factors use
          `C_{11}` and `C_{22}` of each ply in laminate axes and the
          direction-wise neutral surfaces, being exact only for specially
          orthotropic plies. For the direction 1, with `C_{11} = q_{11L}` and
          `C_{55} = q_{55L}` of each ply:

          .. math::

              z_n = \frac{\int C_{11} z dz}{\int C_{11} dz}, \quad
              R = \int C_{11} (z - z_n)^2 dz, \quad
              g(z) = -\int_{z_1}^{z} C_{11} (\bar{z} - z_n) d\bar{z}, \quad
              k_{13} = \frac{R^2}{\bar{A}_{55} \int g^2/C_{55} dz}

          i.e. the energy of the transverse shear stress `\tau_{xz} = Q_x
          g(z)/R` of the one-dimensional equilibrium `\tau_{xz,z} =
          -\sigma_{xx,x}` is equated with that of the FSDT, `Q_x^2/(k_{13}
          \bar{A}_{55})`. Analogously for `k_{23}` with `C_{22} = q_{22L}`,
          `C_{44} = q_{44L}` and `\bar{A}_{44}`. The integrals of `g` are
          evaluated with three Gauss-Legendre points per ply, which is exact,
          `g^2` being a polynomial of degree 4 within each ply.

        - ``'whitney'``: the scalar factors of Whitney (1973) [3]_, Eqs.
          (3)-(7), extending Chow (1971) [4]_ to unsymmetric orthotropic
          laminates in cylindrical bending. The derivative of the bending
          stress, Eq. (3), `\sigma_{xx,x} = -C_{11}(B_{11} - A_{11}z)
          Q_x/(D_{11}A_{11} - B_{11}^2)`, vanishes at the neutral surface
          `z_n = B_{11}/A_{11}` and, since `R = (D_{11}A_{11} -
          B_{11}^2)/A_{11}`, equals `C_{11}(z - z_n) Q_x/R`. Integrated with
          the equilibrium `\sigma_{xx,x} + \tau_{xz,z} = 0`, Eq. (4), it gives
          the `\tau_{xz}` of ``'vlachoutsis'``, Eq. (5), and the energy
          equivalence, Eqs. (6)-(7), gives the same `k_{13}`. The result is
          identical to ``'vlachoutsis'``.

        - ``'chow'``: the scalar factors of Chow (1971) [4]_, Eqs. (11)-(12),
          with the transverse shear stress of a symmetric laminate obtained
          from the one-dimensional equilibrium `\tau_{xz} = Q_x g(z)/D_{11}`,
          `g(z) = -\int_{-h/2}^{z} C_{11} \bar{z} d\bar{z}`, and `k_{13} =
          D_{11}^2 / (\bar{A}_{55} \int g^2/C_{55} dz)`, analogously for
          `k_{23}` with `C_{22}`, `D_{22}` and `C_{44}`. The coordinate is
          measured from the mid-surface and the laminate must be symmetric
          about it, see :func:`._check_symmetric`. For symmetric laminates
          the result is identical to ``'vlachoutsis'``, whose neutral
          surfaces then coincide with the mid-surface. The factors are
          applied as in ``'vlachoutsis'``.

        - ``'birman_bert'``: the average shear strain factors of Birman and
          Bert (2002) [5]_, Eq. (13) and Appendix B, obtained by equating the
          thickness integral of the transverse shear strain of the
          equilibrium distribution with that of the FSDT, `k_{13} = R_1 h /
          (\bar{A}_{55} \int g_1/C_{55} dz)`, analogously for `k_{23}`. The
          equilibrium distribution `g_\alpha(z)` and the bending stiffness
          `R_\alpha` about the direction-wise neutral surface are the same as
          in ``'vlachoutsis'``, which generalizes Eq. (13), derived for
          symmetric cross sections, to unsymmetric laminates. It gives `k =
          1` for a homogeneous plate. The factors are applied as in
          ``'vlachoutsis'``. Note that the value recommended by Birman and
          Bert (2002) for sandwich structures, `K = 1` applied to `\bar{A}_{55}
          = \sum_k C_{55}^{(k)} h_k` (their Eqs. A5-A7), is ``None``.

        - ``'thickness_shear'``: dynamic factors of Yang, Norris and Stavsky
          (1966) [6]_, as described by Sun and Whitney (1973) [7]_, extending
          the criterion of Mindlin (1951) [8]_ to laminates, i.e. matching
          the frequency of the first thickness-shear mode of the FSDT,
          `\omega^2 = k_{13} \bar{A}_{55} / \bar{I}_2`, with the exact
          frequency of the laminate, computed with a transfer matrix through
          the plies for the anti-plane motion `u(z)`, with `\rho^{(k)}
          \ddot{u} = (C_{55}^{(k)} u_{,z})_{,z}`, traction-free faces and
          continuity of `u` and `\tau_{xz}` at every interface, see
          :func:`._thickness_shear_frequency`. The rotatory inertia about the
          centre of mass, `\bar{I}_2 = I_2 - I_1^2/I_0`, with `I_n = \int
          \rho z^n dz`, makes the result independent of ``offset``.
          Analogously for `k_{23}` with `C_{44}`. The coupling `C_{45}`
          between the two directions is neglected in the eigenproblem, and
          the factors are applied as in ``'vlachoutsis'``. It gives `k =
          \pi^2/12` for a homogeneous plate. The ply densities must be
          positive.

        - ``'constant'``: `k = 5/6`, the factor that ``'rohwer'``,
          ``'vlachoutsis'``, ``'whitney'`` and ``'chow'`` give for a
          homogeneous plate.

        References
        ----------
        .. [1] Rohwer, K., Improved transverse shear stiffnesses for layered
           finite elements, DFVLR-FB 88-32, 1988.
        .. [2] Vlachoutsis, S., Shear correction factors for plates and
           shells, International Journal for Numerical Methods in
           Engineering, 33, 1537-1552, 1992.
        .. [3] Whitney, J. M., Shear correction factors for orthotropic
           laminates under static load, Journal of Applied Mechanics, 40,
           302-304, 1973.
        .. [4] Chow, T. S., On the propagation of flexural waves in an
           orthotropic laminated plate and its response to an impulsive load,
           Journal of Composite Materials, 5, 306-319, 1971.
        .. [5] Birman, V. and Bert, C. W., On the choice of shear correction
           factor in sandwich structures, Journal of Sandwich Structures and
           Materials, 4, 83-95, 2002.
        .. [6] Yang, P. C., Norris, C. H. and Stavsky, Y., Elastic wave
           propagation in heterogeneous plates, International Journal of
           Solids and Structures, 2, 665-684, 1966.
        .. [7] Sun, C. T. and Whitney, J. M., Theories for the dynamic
           response of laminated plates, AIAA Journal, 11, 178-183, 1973.
        .. [8] Mindlin, R. D., Influence of rotatory inertia and shear on
           flexural motions of isotropic, elastic plates, Journal of Applied
           Mechanics, 18, 31-38, 1951.

        """
        cdef int k, ig, N, alpha, singular_ply
        cdef double h, det, za, zb, zm, dz, zg, zg2, wdz
        cdef double i44, i45, i55
        cdef double f00, f01, f10, f11, g00, g01, g10, g11
        cdef double S00, S01, S11, detS
        cdef double Sbb44, Sbb45, Sbb55
        cdef double Dk, Gk, num, den, zn, R, d, I, J, gz, gza, kappa
        cdef double rhok, I0, I1, I2, I2c, omega
        cdef double k13, k23
        cdef double [:, :, :, ::1] fc
        cdef double [::1] z
        cdef double [3] gp, gw
        cdef Lamina ply

        # 3-point Gauss-Legendre, exact for polynomials up to degree 5
        gp[0] = -0.7745966692414834; gp[1] = 0.; gp[2] = 0.7745966692414834
        gw[0] = 5./9.; gw[1] = 8./9.; gw[2] = 5./9.

        mode = self.shear_correction
        if not (mode is None or mode in _SHEAR_CORRECTIONS):
            raise ValueError("shear_correction must be one of %s or None, "
                             "got %r" % (', '.join([repr(m) for m in
                             _SHEAR_CORRECTIONS]), mode))

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
        # NOTE heights measured from the mid-surface, since none of the
        #      methods depends on offset, avoiding the round-off of large ones
        for k in range(N + 1):
            z[k] += -h/2.
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

        elif mode in ('vlachoutsis', 'whitney', 'chow', 'birman_bert'):
            if mode == 'chow':
                _check_symmetric(self, z, h)
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
                    raise ValueError('%r shear correction factors require '
                                     'positive in-plane stiffnesses' % mode)
                if mode == 'chow':
                    # NOTE Chow measures z from the mid-surface, which is the
                    #      neutral surface of a symmetric laminate
                    zn = 0.
                else:
                    zn = num/den
                # I = int g^2/G dz (energy), J = int g/G dz (average strain)
                R = 0; I = 0; J = 0; gza = 0
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
                        J += gw[ig]*dz*gz/Gk
                    # g at the top interface of this ply
                    gza = gza - Dk*(0.5*(zb*zb - za*za) - zn*(zb - za))
                if mode == 'birman_bert':
                    kappa = R*h/(d*J)
                else:
                    kappa = R*R/(d*I)
                if alpha == 0:
                    k13 = kappa
                else:
                    k23 = kappa
            self.A44 = k23*self.Abar44
            self.A45 = (k13 + k23)/2.*self.Abar45
            self.A55 = k13*self.Abar55

        elif mode == 'thickness_shear':
            k13 = 0; k23 = 0
            for alpha in range(2):
                Gs = np.zeros(N, dtype=DOUBLE)
                rhos = np.zeros(N, dtype=DOUBLE)
                hs = np.zeros(N, dtype=DOUBLE)
                I0 = 0; I1 = 0; I2 = 0; d = 0
                for k in range(N):
                    ply = self.plies[k]
                    Gk = ply.q55L if alpha == 0 else ply.q44L
                    if not Gk > 0:
                        raise _singular_Cs_error(k, ply)
                    rhok = ply.matlamina.rho
                    if not rhok > 0:
                        raise ValueError("Ply %d (plyid=%d) has density "
                                "rho=%g, 'thickness_shear' requires positive "
                                "densities" % (k, ply.plyid, rhok))
                    za = z[k]
                    zb = z[k+1]
                    Gs[k] = Gk
                    rhos[k] = rhok
                    hs[k] = ply.h
                    I0 += rhok*(zb - za)
                    I1 += rhok*(zb*zb - za*za)/2.
                    I2 += rhok*(zb*zb*zb - za*za*za)/3.
                    d += Gk*ply.h
                I2c = I2 - I1*I1/I0
                omega = _thickness_shear_frequency(Gs, rhos, hs,
                                                   (d/I2c)**0.5)
                kappa = omega*omega*I2c/d
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

        Evaluates, for the transverse shear forces ``Qy`` and ``Qx``, the
        through-thickness distribution of the transverse shear stresses of the
        equilibrium approach used by ``shear_correction='rohwer'``, see
        :meth:`.calc_transverse_shear_stiffness`, which is quadratic within
        each ply, vanishes at the bottom and top faces and is continuous at
        the ply interfaces, e.g. for failure criteria.

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

        Notes
        -----
        Evaluates the equilibrium distribution of Rohwer (1988) [1]_, see
        :meth:`._calc_transverse_shear_distribution`:

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

        References
        ----------
        .. [1] Rohwer, K., Improved transverse shear stiffnesses for layered
           finite elements, DFVLR-FB 88-32, 1988.

        """
        cdef int k, lo, hi, mid, N
        cdef double tol, zm
        cdef double [::1] zi
        cdef double [:, :, :, ::1] fc

        if not self._ts_ready:
            self._calc_transverse_shear_distribution()
        # NOTE the distribution is stored about the mid-surface
        zi = self._ts_z
        fc = self._ts_fcoef
        zm = z - self.offset
        N = <int>zi.shape[0] - 1
        tol = 1e-12*(zi[N] - zi[0])
        if not (zi[0] - tol <= zm <= zi[N] + tol):
            raise ValueError('z=%g is outside the laminate, [%g, %g]'
                             % (z, zi[0] + self.offset, zi[N] + self.offset))
        # last ply k with zi[k] <= zm
        lo = 0
        hi = N - 1
        while lo < hi:
            mid = (lo + hi + 1)//2
            if zi[mid] <= zm:
                lo = mid
            else:
                hi = mid - 1
        k = lo
        return ((fc[k, 0, 0, 0] + zm*fc[k, 1, 0, 0] + zm*zm*fc[k, 2, 0, 0])*Qy
              + (fc[k, 0, 0, 1] + zm*fc[k, 1, 0, 1] + zm*zm*fc[k, 2, 0, 1])*Qx,
                (fc[k, 0, 1, 0] + zm*fc[k, 1, 1, 0] + zm*zm*fc[k, 2, 1, 0])*Qy
              + (fc[k, 0, 1, 1] + zm*fc[k, 1, 1, 1] + zm*zm*fc[k, 2, 1, 1])*Qx)


    def calc_equilibrium_transverse_shear(Laminate self, eps_x=None,
            eps_y=None, z=None, grad_x=None, grad_y=None):
        r"""Transverse shear stresses recovered from the actual strain gradients

        Integrates through the thickness, from the bottom face, the
        equilibrium of the in-plane stresses of the plies given by the
        derivatives ``eps_x`` and ``eps_y`` of the generalized strains, e.g.
        those of an FSDT solution, at ``npts`` points, and returns the
        transverse shear stresses at the heights ``z``. The stresses vanish at
        the bottom face; at the top face they do not vanish in general.

        Parameters
        ----------
        eps_x, eps_y : array-like, shape ``(6,)`` or ``(npts, 6)``
            Derivatives with respect to `x` and `y` of `\{\varepsilon_{xx}^{(0)},
            \varepsilon_{yy}^{(0)}, \gamma_{xy}^{(0)}, \varepsilon_{xx}^{(1)},
            \varepsilon_{yy}^{(1)}, \gamma_{xy}^{(1)}\}`, at ``npts`` points.
            The first three entries are per unit length and the last three
            per unit length squared. Required; the default ``None`` only
            allows the deprecated aliases ``grad_x`` and ``grad_y``.
        z : array-like or None, optional
            Heights from the reference surface. If ``None``, the three
            Gauss-Legendre points of each ply are used. Default ``None``.
            Length.
        grad_x, grad_y : array-like, optional
            Deprecated, use ``eps_x`` and ``eps_y`` instead. Default
            ``None``; when not ``None``, the value is used as ``eps_x`` or
            ``eps_y`` and a :class:`DeprecationWarning` is issued.

        Returns
        -------
        z, tau_yz, tau_xz : np.ndarray
            The heights, with shape ``(nz,)``, and the stresses, with shape
            ``(npts, nz)``. With ``z=None``, ``nz`` is three times the number
            of plies. Length for ``z`` and force per unit area for the
            stresses.

        Raises
        ------
        ValueError
            If a value of ``z`` is outside the laminate.
        TypeError
            If ``eps_x`` or ``eps_y`` is missing, or given together with its
            deprecated alias ``grad_x`` or ``grad_y``.

        See Also
        --------
        calc_aposteriori_energy

        Notes
        -----
        Integrates the three-dimensional equilibrium equations from the bottom
        face, where the tractions vanish, with the in-plane stresses of each
        ply `\sigma^{(k)} = C^{(k)} (\varepsilon^{(0)} + z\varepsilon^{(1)})`:

        .. math::

            \tau_{xz}(z) = -\int_{z_1}^{z} \left( c_1 \left[
            \varepsilon^{(0)}_{,x} + \bar{z} \varepsilon^{(1)}_{,x} \right]
            + c_3 \left[ \varepsilon^{(0)}_{,y} + \bar{z}
            \varepsilon^{(1)}_{,y} \right] \right) d\bar{z} \\
            \tau_{yz}(z) = -\int_{z_1}^{z} \left( c_2 \left[
            \varepsilon^{(0)}_{,y} + \bar{z} \varepsilon^{(1)}_{,y} \right]
            + c_3 \left[ \varepsilon^{(0)}_{,x} + \bar{z}
            \varepsilon^{(1)}_{,x} \right] \right) d\bar{z}

        where `c_1`, `c_2` and `c_3` are the rows of `C^{(k)}` that give
        `\sigma_{xx}`, `\sigma_{yy}` and `\tau_{xy}`, respectively, and `z_1`
        is the bottom face. It uses the gradients of the generalized strains
        of the actual solution, as in the a posteriori approach of Noor and
        Peters (1989) [1]_, instead of the two cylindrical bending states of
        :meth:`.calc_transverse_shear_stiffness`. The top face is traction
        free only if the gradients satisfy the in-plane equilibrium of the
        stress resultants. The integrals of each ply are computed exactly by
        :func:`._equilibrium_integrals`, with `z` measured from the
        mid-surface and the gradients `\varepsilon^{(0)}` shifted by
        ``offset`` times `\varepsilon^{(1)}`, which avoids the round-off of
        large values of ``offset``.

        References
        ----------
        .. [1] Noor, A. K. and Peters, J. M., A posteriori estimates for the
           shear correction factors in multi-layered composite cylinders,
           Journal of Engineering Mechanics, 115, 1225-1244, 1989.

        """
        eps_x, eps_y = _strain_derivatives('calc_equilibrium_transverse_shear',
                                           eps_x, eps_y, grad_x, grad_y)
        eps_x = np.atleast_2d(np.asarray(eps_x, dtype=DOUBLE))
        eps_y = np.atleast_2d(np.asarray(eps_y, dtype=DOUBLE))
        # NOTE integrated about the mid-surface, avoiding the round-off of
        #      large offsets: eps0 + z eps1 = (eps0 + offset eps1) + zm eps1
        off = self.offset
        eps_x = np.hstack((eps_x[:, :3] + off*eps_x[:, 3:], eps_x[:, 3:]))
        eps_y = np.hstack((eps_y[:, :3] + off*eps_y[:, 3:], eps_y[:, 3:]))
        zi = self._ply_interfaces()
        N = len(self.plies)
        if z is None:
            xi = np.array([-0.7745966692414834, 0., 0.7745966692414834])
            zm = np.concatenate([(zi[k] + zi[k+1])/2 + (zi[k+1] - zi[k])/2*xi
                                 for k in range(N)])
            z = zm + off
        else:
            z = np.atleast_1d(np.asarray(z, dtype=DOUBLE))
            zm = z - off
        # NOTE no negative indices, the module is compiled with wraparound=False
        tol = 1e-12*(zi[N] - zi[0])
        if np.any(zm < zi[0] - tol) or np.any(zm > zi[N] + tol):
            raise ValueError('z outside the laminate')
        acc_x = np.zeros((N + 1, eps_x.shape[0]))
        acc_y = np.zeros((N + 1, eps_x.shape[0]))
        for k in range(N):
            ix, iy = _equilibrium_integrals(self.plies[k], eps_x, eps_y,
                                            zi[k], zi[k+1])
            acc_x[k+1] = acc_x[k] + ix
            acc_y[k+1] = acc_y[k] + iy
        k_of_z = np.clip(np.searchsorted(zi, zm, side='right') - 1, 0, N - 1)
        txz = np.zeros((eps_x.shape[0], z.size))
        tyz = np.zeros((eps_x.shape[0], z.size))
        for j in range(z.size):
            k = int(k_of_z[j])
            ix, iy = _equilibrium_integrals(self.plies[k], eps_x, eps_y,
                                            zi[k], zm[j])
            txz[:, j] = -(acc_x[k] + ix)
            tyz[:, j] = -(acc_y[k] + iy)
        return z, tyz, txz


    def calc_aposteriori_energy(Laminate self, eps_x=None, eps_y=None,
            grad_x=None, grad_y=None):
        r"""Transverse shear strain energy of the recovered stresses

        With the stresses of :meth:`.calc_equilibrium_transverse_shear` at
        the three Gauss-Legendre points of each ply, computes at each point
        the transverse shear strain energies per unit area of reference
        surface in the directions 13 and 23, and the resultants of the
        stresses, all integrated through the thickness with three
        Gauss-Legendre points per ply.

        Parameters
        ----------
        eps_x, eps_y : array-like
            See :meth:`.calc_equilibrium_transverse_shear`.
        grad_x, grad_y : array-like, optional
            Deprecated, use ``eps_x`` and ``eps_y`` instead, see
            :meth:`.calc_equilibrium_transverse_shear`.

        Returns
        -------
        U13, U23, Qy, Qx : np.ndarray
            The energies and the resultants `\int \bar\tau dz` of the
            recovered stresses, with shape ``(npts,)``. Energy per unit area
            for ``U13`` and ``U23``, force per unit length for ``Qy`` and
            ``Qx``.

        Raises
        ------
        ValueError
            If a ply has a singular transverse shear constitutive matrix
            `C_s^{(k)}`, e.g. ``g13 = 0`` or ``g23 = 0``.
        TypeError
            If ``eps_x`` or ``eps_y`` is missing, or given together with its
            deprecated alias ``grad_x`` or ``grad_y``.

        Notes
        -----
        With the stresses of :meth:`.calc_equilibrium_transverse_shear` and
        the strains `\{\bar\gamma_{yz}, \bar\gamma_{xz}\} = (C_s^{(k)})^{-1}
        \{\bar\tau_{yz}, \bar\tau_{xz}\}`, the energies per unit area of
        reference surface of Noor and Peters (1989) [1]_, Eq. (7):

        .. math::

            \bar{U}_{13} = \int \frac{1}{2} C_{55} \bar\gamma_{xz}^2 dz
            \qquad
            \bar{U}_{23} = \int \frac{1}{2} C_{44} \bar\gamma_{yz}^2 dz

        integrated exactly with three Gauss-Legendre points per ply, since the
        recovered stresses are quadratic in `z` within each ply. The
        a posteriori factors follow by equating their integral over the
        domain with that of the FSDT, `\frac{1}{2} Q_x^2/\bar{A}_{55}` and
        `\frac{1}{2} Q_y^2/\bar{A}_{44}`, i.e. Eqs. (8)-(9) of Noor and
        Peters (1989), with `Q = k^0 \bar{A} \gamma^0` the shear forces of
        the FSDT. The resultants are `Q_y = \int \bar\tau_{yz} dz` and `Q_x =
        \int \bar\tau_{xz} dz`.

        References
        ----------
        .. [1] Noor, A. K. and Peters, J. M., A posteriori estimates for the
           shear correction factors in multi-layered composite cylinders,
           Journal of Engineering Mechanics, 115, 1225-1244, 1989.

        """
        eps_x, eps_y = _strain_derivatives('calc_aposteriori_energy', eps_x,
                                           eps_y, grad_x, grad_y)
        z, tyz, txz = self.calc_equilibrium_transverse_shear(eps_x, eps_y)
        zi = self._ply_interfaces()
        N = len(self.plies)
        w = np.array([5./9., 8./9., 5./9.])
        U13 = np.zeros(tyz.shape[0])
        U23 = np.zeros(tyz.shape[0])
        Qx = np.zeros(tyz.shape[0])
        Qy = np.zeros(tyz.shape[0])
        for k in range(N):
            ply = self.plies[k]
            det = ply.q44L*ply.q55L - ply.q45L*ply.q45L
            if (ply.q44L <= 0 or ply.q55L <= 0
                    or det <= 1e-12*ply.q44L*ply.q55L):
                raise _singular_Cs_error(k, ply)
            half = (zi[k+1] - zi[k])/2.
            for i in range(3):
                j = 3*k + i
                gyz = (ply.q55L*tyz[:, j] - ply.q45L*txz[:, j])/det
                gxz = (-ply.q45L*tyz[:, j] + ply.q44L*txz[:, j])/det
                U13 += w[i]*half*0.5*ply.q55L*gxz**2
                U23 += w[i]*half*0.5*ply.q44L*gyz**2
                Qx += w[i]*half*txz[:, j]
                Qy += w[i]*half*tyz[:, j]
        return U13, U23, Qy, Qx


    def _ply_interfaces(Laminate self):
        r"""Heights of the ply interfaces measured from the mid-surface"""
        h = 0.
        for ply in self.plies:
            h += ply.h
        return (-h/2. + np.concatenate(([0.], np.cumsum([ply.h for ply in
                                                         self.plies]))))


    cpdef void calc_equivalent_properties(Laminate self):
        r"""Calculate the equivalent laminate properties

        The following attributes are updated:

            ``e1``, ``e2``, ``g12``, ``nu12``, ``nu21``

        They are computed from the inverse of the 6x6 ``ABD`` matrix and from
        ``h``.

        Raises
        ------
        numpy.linalg.LinAlgError
            A subclass of ``ValueError``, if the ``ABD`` matrix is singular,
            e.g. for a laminate of antiplane cores only, whose in-plane and
            bending stiffnesses vanish: the equivalent properties are then
            undefined.

        Notes
        -----
        With `a = [[A, B], [B, D]]^{-1}`, whose upper-left 3x3 block is the
        membrane compliance of the laminate with zero moments and free
        curvatures:

        .. math::

            e_1 = \frac{1}{h a_{11}}, \quad e_2 = \frac{1}{h a_{22}}, \quad
            g_{12} = \frac{1}{h a_{66}}, \quad
            \nu_{12} = -\frac{a_{12}}{a_{11}}, \quad
            \nu_{21} = -\frac{a_{12}}{a_{22}}

        where the index 6 is the third row and column of `a`.

        """
        msg = ('The ABD matrix of the laminate is singular, its equivalent '
               'properties are undefined, e.g. for a laminate of antiplane '
               'cores only, whose in-plane and bending stiffnesses vanish')
        #NOTE explicit checks, np.linalg.inv does not raise LinAlgError for
        #     a singular matrix on every platform, e.g. not with the LAPACK
        #     of Pyodide, which returns inf or nan instead
        ABD = np.asarray(self.get_ABD(), dtype=DOUBLE)
        if not np.all(np.diag(ABD) > 0):
            raise np.linalg.LinAlgError(msg)
        try:
            AI = np.linalg.inv(ABD)
        except np.linalg.LinAlgError:
            raise np.linalg.LinAlgError(msg)
        if not np.all(np.isfinite(AI)):
            raise np.linalg.LinAlgError(msg)
        a11, a12, a22, a33 = AI[0,0], AI[0,1], AI[1,1], AI[2,2]
        self.e1 = 1./(self.h*a11)
        self.e2 = 1./(self.h*a22)
        self.g12 = 1./(self.h*a33)
        self.nu12 = - a12 / a11
        self.nu21 = - a12 / a22


    cpdef void calc_constitutive_matrix(Laminate self):
        r"""Calculate the laminate constitutive terms

        This is the commonly called ``ABD`` matrix with ``shape=(6, 6)`` when
        the classical laminated plate theory is used, or the ``ABD`` matrix
        when the first-order shear deformation theory is used, containing the
        transverse shear terms.

        The transverse shear stiffnesses ``A44``, ``A45``, ``A55`` are
        calculated at the end by :meth:`.calc_transverse_shear_stiffness`,
        with the shear correction selected by ``shear_correction`` already
        applied.

        The attributes ``h``, ``intrho``, ``intrhoz``, ``intrhoz2``,
        ``intrhoz3``, ``intrhoz4``, ``intrhoz6``, ``A11`` to ``H66``,
        ``D44``, ``D45``, ``D55``, ``F44``, ``F45`` and ``F55`` are computed
        from ``plies`` and ``offset``, with the ply stiffnesses of the last
        call to :meth:`.Lamina.rebuild` and the density ``rho`` of the
        material of each ply.

        Notes
        -----
        With the plies running from `z_1 = -h/2 + offset` to `z_{N+1} = +h/2
        + offset`, the stiffnesses of the CLPT, FSDT and TSDT in the notation
        of Reddy [1]_ and the mass integrals are, integrated exactly ply by
        ply:

        .. math::

            (A_{ij}, B_{ij}, D_{ij}, E_{ij}, F_{ij}, H_{ij}) = \sum_k
            \int_{z_k}^{z_{k+1}} q_{ijL}^{(k)} (1, z, z^2, z^3, z^4, z^6) dz,
            \quad ij = 11, 12, 16, 22, 26, 66 \\
            (D_{ij}, F_{ij}) = \sum_k \int_{z_k}^{z_{k+1}} q_{ijL}^{(k)}
            (z^2, z^4) dz, \quad ij = 44, 45, 55 \\
            (intrho, intrhoz, intrhoz2, intrhoz3, intrhoz4, intrhoz6) =
            \sum_k \int_{z_k}^{z_{k+1}} \rho^{(k)} (1, z, z^2, z^3, z^4, z^6)
            dz

        and ``A44``, ``A45``, ``A55`` follow from
        :meth:`.calc_transverse_shear_stiffness`.

        The mass integrals are the moments of the density through the
        thickness, `I_n = \sum_k \int_{z_k}^{z_{k+1}} \rho^{(k)} z^n dz`,
        required by the kinetic energy `T = \frac{1}{2} \int_V \rho
        (\dot{u}^2 + \dot{v}^2 + \dot{w}^2) dV` of the displacement field of
        each theory. The in-plane displacements of the CLPT and of the FSDT
        contain the powers `z^0` and `z`, requiring `I_0`, `I_1`, `I_2`,
        while those of the TSDT [1]_ contain `z^0`, `z` and `z^3`, whose
        products give the powers `z^0` to `z^4` and `z^6`, requiring `I_0`,
        `I_1`, `I_2`, `I_3`, `I_4`, `I_6`.

        References
        ----------
        .. [1] Reddy, J. N., Mechanics of Laminated Composite Plates and Shells,
           CRC Press, 2003. https://doi.org/10.1201/b12409

        """
        cdef double h0, hk_1, hk, tmp_hk, tmp_hk_1
        self.h = 0.
        self.intrho = 0.
        self.intrhoz = 0.
        self.intrhoz2 = 0.
        self.intrhoz3 = 0.
        self.intrhoz4 = 0.
        self.intrhoz6 = 0.
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
            self.intrhoz3 += ply.matlamina.rho*(hk*hk*hk*hk/4. - hk_1*hk_1*hk_1*hk_1/4.)
            self.intrhoz4 += ply.matlamina.rho*(hk*hk*hk*hk*hk/5. - hk_1*hk_1*hk_1*hk_1*hk_1/5.)
            self.intrhoz6 += ply.matlamina.rho*(hk*hk*hk*hk*hk*hk*hk/7. - hk_1*hk_1*hk_1*hk_1*hk_1*hk_1*hk_1/7.)

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

        The attributes ``A16``, ``A26``, ``Abar45`` and ``Abarbar45``, and
        ``A45`` unless ``shear_correction='rohwer'``, are set to zero to make
        a balanced laminate. ``Abarbar45`` is kept when it is ``nan``, i.e.
        when a ply has a singular `C_s^{(k)}`. The other attributes are not
        modified, among them ``B16``, ``B26``, ``D16``, ``D26``, ``E16``,
        ``E26``, ``F16``, ``F26``, ``H16``, ``H26``, ``D45`` and ``F45``,
        which do not vanish in a balanced laminate in general, e.g. ``B16``
        and ``B26`` of an antisymmetric `[+45/-45]` laminate.

        Raises
        ------
        RuntimeError
            If ``offset`` is not zero.

        Notes
        -----
        In a balanced laminate, for every ply at `+\theta` there is a ply at
        `-\theta` of the same material and thickness, anywhere in the stack.
        `q_{16L}`, `q_{26L}` and `q_{45L}` are odd functions of the ply angle
        and the other `q_{ijL}` are even, see :meth:`.Lamina.rebuild`. The
        contributions of such a pair of plies, of thickness `t` and with
        bottom faces at `z_a` and `z_b`, to the integral of `q_{16L} z^n`
        add up to

        .. math::

            q_{16L}(\theta) \int_{z_a}^{z_a + t} z^n dz + q_{16L}(-\theta)
            \int_{z_b}^{z_b + t} z^n dz = q_{16L}(\theta) \left(
            \int_{z_a}^{z_a + t} z^n dz - \int_{z_b}^{z_b + t} z^n dz
            \right)

        which vanishes for any `z_a`, `z_b` only for `n = 0`, see
        :meth:`.calc_constitutive_matrix`. Therefore, `A_{16}`, `A_{26}` and
        `\bar{A}_{45} = \sum_k q_{45L}^{(k)} h_k` vanish in every balanced
        laminate, whereas the integrals with `n \ge 1`, `B_{16}` to `H_{26}`,
        `D_{45}` and `F_{45}`, are proportional to `z_a - z_b` and depend on
        the positions of the plies. The determinant of `C_s^{(k)}` being
        even in `\theta`, the off-diagonal term of `(C_s^{(k)})^{-1}` is odd,
        its thickness integral vanishes and so does `\bar{\bar{A}}_{45}`.
        ``A45`` is proportional to ``Abar45`` for every ``shear_correction``
        other than ``'rohwer'``, see :meth:`.calc_transverse_shear_stiffness`,
        and vanishes with it. With ``'rohwer'``, ``A45`` follows from the
        through-thickness distribution of the transverse shear stresses,
        which depends on the positions of the plies and on the whole ABD
        matrix, and does not vanish in general, e.g. for a `[+30/-30/0]_s`
        laminate, even with ``g13 = g23``, i.e. `q_{45L} = 0` in every ply.

        """
        if self.offset != 0.:
            raise RuntimeError('Laminates with offset cannot be made balanced!')
        self.A16 = 0.
        self.A26 = 0.

        # transverse shear
        self.Abar45 = 0.
        if not np.isnan(self.Abarbar45):
            self.Abarbar45 = 0.
        if self.shear_correction != 'rohwer':
            self.A45 = 0.


    cpdef void make_orthotropic(Laminate self):
        r"""Make an orthotropic laminate

        The attributes ``A16``, ``A26``, ``B16``, ``B26``, ``D16``, ``D26``,
        ``E16``, ``E26``, ``F16``, ``F26``, ``H16``, ``H26`` and the
        transverse shear couplings ``A45``, ``Abar45``, ``Abarbar45``,
        ``D45``, ``F45`` are set to zero to make an orthotropic laminate.
        ``Abarbar45`` is kept when it is ``nan``, i.e. when a ply has a
        singular `C_s^{(k)}`, as in :meth:`.make_balanced`. The other
        attributes are not modified.

        Raises
        ------
        RuntimeError
            If ``offset`` is not zero.

        Notes
        -----
        These terms vanish in a laminate whose plies are all at 0 or 90
        degrees, for which `q_{16L} = q_{26L} = q_{45L} = 0`, see
        :meth:`.Lamina.rebuild`. The terms `ij = 16, 26` of `A_{ij}` to
        `H_{ij}` and `D_{45}`, `F_{45}`, `\bar{A}_{45}` are thickness integrals
        of `q_{16L}`, `q_{26L}` and `q_{45L}`, see
        :meth:`.calc_constitutive_matrix` and
        :meth:`.calc_transverse_shear_stiffness`, and `\bar{\bar{A}}_{45}`
        vanishes with the diagonal `C_s^{(k)}` of every ply. With every ply
        diagonal and the ABD matrix decoupled between `\{N_{xx}, N_{yy},
        M_{xx}, M_{yy}\}` and `\{N_{xy}, M_{xy}\}`, the transverse shear
        distribution of ``'rohwer'`` is diagonal, and ``A45`` vanishes for
        every ``shear_correction``.

        """
        if self.offset != 0.:
            raise RuntimeError('Laminates with offset cannot be made orthotropic!')
        self.A16 = 0.
        self.A26 = 0.
        self.B16 = 0.
        self.B26 = 0.
        self.D16 = 0.
        self.D26 = 0.

        # TSDT
        self.E16 = 0.
        self.E26 = 0.
        self.F16 = 0.
        self.F26 = 0.
        self.H16 = 0.
        self.H26 = 0.

        # transverse shear
        self.A45 = 0.
        self.Abar45 = 0.
        if not np.isnan(self.Abarbar45):
            self.Abarbar45 = 0.
        self.D45 = 0.
        self.F45 = 0.


    cpdef void make_symmetric(Laminate self):
        """Make a symmetric laminate

        The `B_{ij}` terms of the constitutive matrix are set to zero.

        The attributes ``E11`` to ``E66`` of the TSDT are also set to zero.
        The other attributes are not modified.

        Raises
        ------
        RuntimeError
            If ``offset`` is not zero.

        Notes
        -----
        In a laminate symmetric about its mid-surface the integrands of the
        odd powers of `z` are odd functions of `z`, see
        :meth:`.calc_constitutive_matrix`, such that `B_{ij} = E_{ij} = 0`,
        `ij = 11, 12, 16, 22, 26, 66`, the integrals of `z` and `z^3`. The
        integrands of the even powers of `z`, i.e. of `A_{ij}`, `D_{ij}`,
        `F_{ij}`, `H_{ij}` and of the transverse shear terms `D_{ij}`,
        `F_{ij}`, `ij = 44, 45, 55`, are even functions of `z`, which the
        symmetry does not make zero.

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


    cpdef void make_smeared(Laminate self):
        r"""Make a laminated with smeared properties

        The `B_{ij}` terms of the constitutive matrix are set to zero.

        The `D_{ij}` terms are calculated from the membrane terms `A_{ij}`
        according to `D_{ij} = (h^2 A_{ij})/12`, where `h` is the
        laminate thickness.

        Raises
        ------
        NotImplementedError
            If ``offset`` is not zero.

        Notes
        -----
        `D_{ij} = h^2 A_{ij}/12` and `B_{ij} = 0` are the stiffnesses of a
        homogeneous plate with the in-plane stiffnesses `A_{ij}/h`, i.e. of a
        laminate whose plies are smeared through the thickness, the limit of
        many repetitions of the stacking sequence.

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

            ``xiA``, ``xiB``, ``xiD``, ``xiAts``

        They are the attributes of the returned object, computed from the ply
        angles and from the ply heights measured from the reference surface,
        i.e. including ``offset``, divided by the total thickness. The
        object is not modified, except that :meth:`.Lamina.rebuild` is called
        for each ply.

        Returns
        -------
        lp : :class:`.LaminationParameters`
            New object with the lamination parameters, dimensionless.

        Raises
        ------
        ValueError
            If there are no plies.

        Notes
        -----
        With `\bar{z} = z/h`, `z` measured from the reference surface, and
        `\theta(\bar{z})` the ply angle:

        .. math::

            \{\xi_{A1}, \xi_{A2}, \xi_{A3}, \xi_{A4}\} = \int \{\cos 2\theta,
            \sin 2\theta, \cos 4\theta, \sin 4\theta\} d\bar{z} \\
            \{\xi_{B1}, \xi_{B2}, \xi_{B3}, \xi_{B4}\} = 4 \int \{\cos
            2\theta, \sin 2\theta, \cos 4\theta, \sin 4\theta\} \bar{z}
            d\bar{z} \\
            \{\xi_{D1}, \xi_{D2}, \xi_{D3}, \xi_{D4}\} = 12 \int \{\cos
            2\theta, \sin 2\theta, \cos 4\theta, \sin 4\theta\} \bar{z}^2
            d\bar{z} \\
            \{\xi_{A_{ts}1}, \xi_{A_{ts}2}\} = \int \{\cos 2\theta, \sin
            2\theta\} d\bar{z}

        integrated exactly ply by ply, such that ``xiAts1 = xiA1`` and
        ``xiAts2 = xiA2``. With these definitions, the stiffnesses of the
        laminate are linear in the lamination parameters, see
        :func:`.laminate_from_LaminationParameters`, which must be given the
        same ``offset`` to reproduce the stiffnesses of the laminate. The
        `\xi_{Ai}` and `\xi_{A_{ts}i}` do not depend on ``offset``, the
        `\xi_{Bi}` and `\xi_{Di}` do.

        The `\xi_{Ai}` and `\xi_{Di}` are those of Krupa et al. [1]_, who
        integrate from `z = -h/2` to `z = h/2`, i.e. with ``offset = 0``. The
        `\xi_{Bi}` and `\xi_{A_{ts}i}`, not defined there, follow the same
        pattern, with the factor 4 of `\xi_{Bi}` corresponding to the factor
        `h^2/4` of `B_{ij}` in :func:`.laminate_from_LaminationParameters`.

        References
        ----------
        .. [1] Krupa, E. P., Cooper, J. E., Pirrera, A. and Nangia, R.,
           Improved aerostructural performance via aeroservoelastic
           tailoring of a composite wing, The Aeronautical Journal, 122
           (1255), 1442-1474, 2018. https://doi.org/10.1017/aer.2018.66

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

            lp.xiAts1 += Efac * ply.cos2t
            lp.xiAts2 += Efac * ply.sin2t

        return lp


cpdef LaminationParameters make_balanced_LP(LaminationParameters lp):
    r"""Make balanced lamination parameters

    The lamination parameters `\xi_{A2}`, `\xi_{A4}` and `\xi_{A_{ts}2}` are
    set to null to make a balanced laminate. The other lamination parameters
    are not modified, among them `\xi_{B2}`, `\xi_{B4}`, `\xi_{D2}` and
    `\xi_{D4}`, which do not vanish in a balanced laminate in general.

    The object ``lp`` is modified in place and returned.

    Parameters
    ----------
    lp : :class:`.LaminationParameters`
        Lamination parameters, modified in place.

    Returns
    -------
    lp : :class:`.LaminationParameters`
        The same object given as input.

    See Also
    --------
    make_symmetric_LP, make_orthotropic_LP

    Notes
    -----
    In a balanced laminate, with a ply at `-\theta` for every ply at
    `+\theta` of the same material and thickness, anywhere in the stack, the
    thickness integrals of the odd functions `\sin 2\theta` and `\sin
    4\theta` vanish, see :meth:`.Laminate.calc_lamination_parameters`, and
    so do `A_{16}`, `A_{26}` and `\bar{A}_{45}`, see
    :meth:`.Laminate.make_balanced`. The integrals weighted by `\bar{z}` and
    `\bar{z}^2`, `\xi_{B2}`, `\xi_{B4}`, `\xi_{D2}` and `\xi_{D4}`, depend on
    the positions of the plies, e.g. `\xi_{B2} \ne 0` for an antisymmetric
    `[+45/-45]` laminate.

    """
    lp.xiA2 = 0
    lp.xiA4 = 0
    lp.xiAts2 = 0
    return lp


cpdef LaminationParameters make_symmetric_LP(LaminationParameters lp):
    r"""Make symmetric lamination parameters

    The lamination parameters `\xi_{Bi}` are set to null to make a symmetric
    laminate.

    The object ``lp`` is modified in place and returned.

    Parameters
    ----------
    lp : :class:`.LaminationParameters`
        Lamination parameters, modified in place.

    Returns
    -------
    lp : :class:`.LaminationParameters`
        The same object given as input.

    See Also
    --------
    make_balanced_LP, make_orthotropic_LP

    Notes
    -----
    In a laminate symmetric about its mid-surface, the integrands of
    `\xi_{Bi}`, linear in `\bar{z}`, are odd functions of `\bar{z}`, see
    :meth:`.Laminate.calc_lamination_parameters`, such that `\xi_{Bi} = 0`
    and the coupling stiffnesses `B_{ij}` vanish. This holds for the
    lamination parameters about the mid-surface, i.e. with ``offset = 0``,
    as :meth:`.Laminate.make_symmetric`: about another reference surface
    `\xi_{Bi} = 4 (e/h) \xi_{Ai}` in a symmetric laminate, see
    :func:`.laminate_from_LaminationParameters`.

    """
    lp.xiB1 = 0
    lp.xiB2 = 0
    lp.xiB3 = 0
    lp.xiB4 = 0
    return lp


cpdef LaminationParameters make_orthotropic_LP(LaminationParameters lp):
    r"""Make orthotropic lamination parameters

    The lamination parameters `\xi_{A2}`, `\xi_{A4}`, `\xi_{B2}`, `\xi_{B4}`,
    `\xi_{D2}`, `\xi_{D4}` and `\xi_{A_{ts}2}` are set to null to make an
    orthotropic laminate.

    The object ``lp`` is modified in place and returned.

    Parameters
    ----------
    lp : :class:`.LaminationParameters`
        Lamination parameters, modified in place.

    Returns
    -------
    lp : :class:`.LaminationParameters`
        The same object given as input.

    See Also
    --------
    make_balanced_LP, make_symmetric_LP

    Notes
    -----
    The lamination parameters of the sine terms give the stiffnesses
    `A_{16}`, `A_{26}`, `B_{16}`, `B_{26}`, `D_{16}`, `D_{26}` and
    `\bar{A}_{45}`, see :func:`.laminate_from_LaminationParameters`, and
    they vanish when all plies are at 0 or 90 degrees, see
    :meth:`.Laminate.make_orthotropic`. The `\xi_{D2}` and `\xi_{D4}` are
    related to the bend-twist coupling and become often very small for
    balanced laminates with a large amount of plies.

    """
    lp.xiA2 = 0
    lp.xiA4 = 0
    lp.xiB2 = 0
    lp.xiB4 = 0
    lp.xiD2 = 0
    lp.xiD4 = 0
    lp.xiAts2 = 0
    return lp


cpdef Laminate laminate_from_LaminationParameters(double thickness, MatLamina
        mat, LaminationParameters lp, double offset=0.):
    r"""Return a :class:`.Laminate` object based in the thickness, material and
    lamination parameters

    In a new :class:`.Laminate` object, the stiffnesses ``A11`` to ``D66`` and
    ``Abar44``, ``Abar45``, ``Abar55`` are computed from ``thickness``,
    ``offset``, from the invariants ``u1`` to ``u7`` of ``mat`` and from the
    lamination parameters of ``lp``, with the combinations of
    :meth:`.MatLamina.get_invariant_matrix`. No plies are created, ``offset``
    is stored, and the mass integrals ``intrho``, ``intrhoz``, ``intrhoz2``,
    ``intrhoz3``, ``intrhoz4``, ``intrhoz6`` and
    the equivalent properties ``e1``, ``e2``, ``g12``, ``nu12``, ``nu21`` are
    not computed, remaining zero.
    The stiffnesses ``E11`` to ``H66``, ``D44``, ``D45``, ``D55``, ``F44``,
    ``F45`` and ``F55`` are also zero.
    The transverse shear stiffnesses ``A44``, ``A45``, ``A55`` are equal to
    the constant-strain ``Abar44``, ``Abar45``, ``Abar55``,
    ``shear_correction`` is ``None``, ``scf_k13 = scf_k23 = 1`` and the
    constant-stress ``Abarbar44``, ``Abarbar45``, ``Abarbar55`` are ``nan``.

    Parameters
    ----------
    thickness : float
        The total thickness of the laminate, length
    mat : :class:`.MatLamina` object
        Material object
    lp : :class:`.LaminationParameters` object
        The container class with all lamination parameters already defined
    offset : float, optional
        Offset of the mid-surface of the laminate from the reference surface,
        as ``offset`` of :class:`.Laminate`, with which ``lp`` was computed
        by :meth:`.Laminate.calc_lamination_parameters`. Length. Default
        ``0.``.

    Returns
    -------
    lam : :class:`.Laminate`
        laminate with the constitutive matrices already calculated

    Notes
    -----
    The stiffnesses of a ply are linear in the trigonometric functions of
    its angle, with the invariants of ``mat`` as coefficients, see
    :meth:`.MatLamina.get_invariant_matrix`, and the lamination parameters
    are their thickness integrals, see
    :meth:`.Laminate.calc_lamination_parameters`, with `\bar{z} = z/h` and
    `z` measured from the reference surface, i.e. the laminate occupies
    `\bar{e} - 1/2 \le \bar{z} \le \bar{e} + 1/2`, `\bar{e}` = ``offset``
    `/h`. For a laminate of a single material, with `h` = ``thickness``,
    `\int dz = h \int d\bar{z}`, `\int z dz = h^2/4 \cdot 4\int \bar{z}
    d\bar{z}` and `\int z^2 dz = h^3/12 \cdot 12 \int \bar{z}^2 d\bar{z}`,
    and the integrals of the constant terms `\int d\bar{z} = 1`, `4 \int
    \bar{z} d\bar{z} = 4\bar{e}` and `12 \int \bar{z}^2 d\bar{z} = 1 + 12
    \bar{e}^2`, give:

    .. math::

        A_{11} = h (u_1 + u_2 \xi_{A1} + u_3 \xi_{A3}), \quad
        A_{12} = h (u_4 - u_3 \xi_{A3}), \quad
        A_{22} = h (u_1 - u_2 \xi_{A1} + u_3 \xi_{A3}) \\
        A_{16} = h (u_2/2 \, \xi_{A2} + u_3 \xi_{A4}), \quad
        A_{26} = h (u_2/2 \, \xi_{A2} - u_3 \xi_{A4}), \quad
        A_{66} = h (u_5 - u_3 \xi_{A3}) \\
        B_{11} = \frac{h^2}{4} (4\bar{e} u_1 + u_2 \xi_{B1} + u_3 \xi_{B3}),
        \quad
        B_{12} = \frac{h^2}{4} (4\bar{e} u_4 - u_3 \xi_{B3}), \quad
        B_{22} = \frac{h^2}{4} (4\bar{e} u_1 - u_2 \xi_{B1} + u_3 \xi_{B3}) \\
        B_{16} = \frac{h^2}{4} (u_2/2 \, \xi_{B2} + u_3 \xi_{B4}), \quad
        B_{26} = \frac{h^2}{4} (u_2/2 \, \xi_{B2} - u_3 \xi_{B4}), \quad
        B_{66} = \frac{h^2}{4} (4\bar{e} u_5 - u_3 \xi_{B3}) \\
        \bar{A}_{44} = h (u_6 + u_7 \xi_{A_{ts}1}), \quad
        \bar{A}_{45} = -h u_7 \xi_{A_{ts}2}, \quad
        \bar{A}_{55} = h (u_6 - u_7 \xi_{A_{ts}1})

    and `D_{ij}` as `A_{ij}`, with the factor `h^3/12` instead of `h`,
    `\xi_{Di}` instead of `\xi_{Ai}` and `(1 + 12\bar{e}^2) u_1`, `(1 +
    12\bar{e}^2) u_4`, `(1 + 12\bar{e}^2) u_5` instead of `u_1`, `u_4`,
    `u_5`. The `\xi_{Bi}` and `\xi_{Di}` about the reference surface are
    related to those about the mid-surface, `\xi_{Bi}^0` and `\xi_{Di}^0`,
    by `\xi_{Bi} = \xi_{Bi}^0 + 4\bar{e} \xi_{Ai}` and `\xi_{Di} =
    \xi_{Di}^0 + 6\bar{e} \xi_{Bi}^0 + 12\bar{e}^2 \xi_{Ai}`, such that
    `B = B^0 + e A` and `D = D^0 + 2 e B^0 + e^2 A`, the shift of the
    stiffness integrals to another reference plane (Vasiliev and Morozov
    [2]_, Eq. (5.28), p. 235), with `B^0`, `D^0` about the mid-surface.
    ``A`` and the transverse shear stiffnesses do not depend on ``offset``,
    and ``lp`` with an offset may lie outside the feasible region of the
    lamination parameters about the mid-surface.

    With the plane-stress invariants of
    :meth:`.MatLamina.rebuild`, its default ``plane_stress=True``, these are
    the stiffnesses of the classical lamination theory, the thickness
    integrals of the plane-stress `q_{ijL}` (Jones [1]_, Eq. (4.24), p.
    198), equal to those of :func:`.laminated_plate` for the same layup and
    ``offset``, also when ``e3``, ``nu13`` and ``nu23`` are given; after
    ``mat.rebuild(plane_stress=False)`` they integrate the rotated
    three-dimensional stiffnesses instead, see :meth:`.MatLamina.rebuild`.
    Since the through-thickness distribution of the plies is not known from
    the lamination parameters, no shear correction can be computed, hence
    the uncorrected transverse shear stiffnesses.

    References
    ----------
    .. [1] Jones, R. M., Mechanics of Composite Materials, Second Edition,
       Taylor & Francis, Philadelphia, 1999.
    .. [2] Vasiliev, V. V. and Morozov, E. V., Mechanics and Analysis of
       Composite Materials, Elsevier, Oxford, 2001.

    """
    cdef double fB, fD
    lam = Laminate()
    lam.h = thickness
    lam.offset = offset
    # NOTE int z dz and int z^2 dz of the constant terms with z measured from
    #      the reference surface: h^2/4 fB = h offset and h^3/12 fD = h^3/12
    #      + h offset^2
    fB = 4*offset/thickness
    fD = 1 + 12*(offset/thickness)**2

    lam.A11 = lam.h*(mat.u1 + mat.u2*lp.xiA1 + 0*lp.xiA2 + mat.u3*lp.xiA3 + 0*lp.xiA4)
    lam.A12 = lam.h*(mat.u4 + 0*lp.xiA1 + 0*lp.xiA2 + (-1)*mat.u3*lp.xiA3 + 0*lp.xiA4)
    lam.A22 = lam.h*(mat.u1 + (-1)*mat.u2*lp.xiA1 + 0*lp.xiA2 + mat.u3*lp.xiA3 + 0*lp.xiA4)
    lam.A16 = lam.h*(0 + 0*lp.xiA1 + mat.u2/2.*lp.xiA2 + 0*lp.xiA3 + mat.u3*lp.xiA4)
    lam.A26 = lam.h*(0 + 0*lp.xiA1 + mat.u2/2.*lp.xiA2 + 0*lp.xiA3 + (-1)*mat.u3*lp.xiA4)
    lam.A66 = lam.h*(mat.u5 + 0*lp.xiA1 + 0*lp.xiA2 + (-1)*mat.u3*lp.xiA3 + 0*lp.xiA4)

    lam.B11 = lam.h*lam.h/4.*(fB*mat.u1 + mat.u2*lp.xiB1 + 0*lp.xiB2 + mat.u3*lp.xiB3 + 0*lp.xiB4)
    lam.B12 = lam.h*lam.h/4.*(fB*mat.u4 + 0*lp.xiB1 + 0*lp.xiB2 + (-1)*mat.u3*lp.xiB3 + 0*lp.xiB4)
    lam.B22 = lam.h*lam.h/4.*(fB*mat.u1 + (-1)*mat.u2*lp.xiB1 + 0*lp.xiB2 + mat.u3*lp.xiB3 + 0*lp.xiB4)
    lam.B16 = lam.h*lam.h/4.*(0*lp.xiB1 + mat.u2/2.*lp.xiB2 + 0*lp.xiB3 + mat.u3*lp.xiB4)
    lam.B26 = lam.h*lam.h/4.*(0*lp.xiB1 + mat.u2/2.*lp.xiB2 + 0*lp.xiB3 + (-1)*mat.u3*lp.xiB4)
    lam.B66 = lam.h*lam.h/4.*(fB*mat.u5 + 0*lp.xiB1 + 0*lp.xiB2 + (-1)*mat.u3*lp.xiB3 + 0*lp.xiB4)

    lam.D11 = lam.h*lam.h*lam.h/12.*(fD*mat.u1 + mat.u2*lp.xiD1 + 0*lp.xiD2 + mat.u3*lp.xiD3 + 0*lp.xiD4)
    lam.D12 = lam.h*lam.h*lam.h/12.*(fD*mat.u4 + 0*lp.xiD1 + 0*lp.xiD2 + (-1)*mat.u3*lp.xiD3 + 0*lp.xiD4)
    lam.D22 = lam.h*lam.h*lam.h/12.*(fD*mat.u1 + (-1)*mat.u2*lp.xiD1 + 0*lp.xiD2 + mat.u3*lp.xiD3 + 0*lp.xiD4)
    lam.D16 = lam.h*lam.h*lam.h/12.*(0 + 0*lp.xiD1 + mat.u2/2.*lp.xiD2 + 0*lp.xiD3 + mat.u3*lp.xiD4)
    lam.D26 = lam.h*lam.h*lam.h/12.*(0 + 0*lp.xiD1 + mat.u2/2.*lp.xiD2 + 0*lp.xiD3 + (-1)*mat.u3*lp.xiD4)
    lam.D66 = lam.h*lam.h*lam.h/12.*(fD*mat.u5 + 0*lp.xiD1 + 0*lp.xiD2 + (-1)*mat.u3*lp.xiD3 + 0*lp.xiD4)

    lam.Abar44 = lam.h*(mat.u6 + mat.u7*lp.xiAts1 + 0*lp.xiAts2)
    lam.Abar45 = lam.h*(0 + 0*lp.xiAts1 + (-1)*mat.u7*lp.xiAts2)
    lam.Abar55 = lam.h*(mat.u6 + (-1)*mat.u7*lp.xiAts1 + 0*lp.xiAts2)
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
        double xiAts1=0, double xiAts2=0, xiAtrans1=None, xiAtrans2=None,
        double offset=0.):
    r"""Return a :class:`.Laminate` object based in the thickness, material and
    lamination parameters

    The 14 lamination parameters according to the first-order shear
    deformation theory are `\xi_{A1} \cdots \xi_{A4}`, `\xi_{B1} \cdots
    \xi_{B4}`, `\xi_{D1} \cdots \xi_{D4}`, `\xi_{A_{ts}1}` and
    `\xi_{A_{ts}2}`. They are stored in a new
    :class:`.LaminationParameters` object, passed to
    :func:`.laminate_from_LaminationParameters`.

    Parameters
    ----------
    thickness : float
        The total thickness of the laminate, length
    matlamina : :class:`.MatLamina` object
        Material object
    xiA1, xiA2, xiA3, xiA4 : float
        Lamination parameters `\xi_{A1} \cdots \xi_{A4}`, dimensionless.
    xiB1, xiB2, xiB3, xiB4 : float
        Lamination parameters `\xi_{B1} \cdots \xi_{B4}`, dimensionless.
    xiD1, xiD2, xiD3, xiD4 : float
        Lamination parameters `\xi_{D1} \cdots \xi_{D4}`, dimensionless.
    xiAts1, xiAts2 : float, optional
        Lamination parameters `\xi_{A_{ts}1}` and `\xi_{A_{ts}2}`,
        dimensionless. Default ``0``.
    xiAtrans1, xiAtrans2 : float, optional
        Deprecated, use ``xiAts1`` and ``xiAts2`` instead. Default ``None``;
        when not ``None``, the value replaces ``xiAts1`` or ``xiAts2`` and a
        :class:`DeprecationWarning` is issued.
    offset : float, optional
        Offset of the mid-surface from the reference surface, about which
        the lamination parameters are defined, see
        :func:`.laminate_from_LaminationParameters`. Length. Default ``0.``.


    Returns
    -------
    lam : :class:`.Laminate`
        laminate with the constitutive matrices already calculated. See
        :func:`.laminate_from_LaminationParameters` for the transverse shear
        stiffnesses.

    Notes
    -----
    Note that `\xi_{A_{ts}1}` and `\xi_{A_{ts}2}` are optional and usually
    equal to zero, becoming important only when the transverse shear modulus
    is different in the two directions, i.e.  when `G_{13} \ne G{23}`, since
    they multiply the invariant `u_7 = (q_{44} - q_{55})/2 = (g_{23} -
    g_{13})/2`, see :func:`.laminate_from_LaminationParameters`.

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
    lp.xiAts1 = xiAts1
    lp.xiAts2 = xiAts2
    if xiAtrans1 is not None:
        _warn_deprecated('xiAtrans1', 'xiAts1')
        lp.xiAts1 = xiAtrans1
    if xiAtrans2 is not None:
        _warn_deprecated('xiAtrans2', 'xiAts2')
        lp.xiAts2 = xiAtrans2
    return laminate_from_LaminationParameters(thickness, matlamina, lp,
                                              offset)


_GRADABD_STATE = ('gradAij', 'gradBij', 'gradDij', 'gradAtsij')


cdef class GradABD:
    r"""Container to store the gradients of the ABD matrices with respect to
    the lamination parameters

    Created with all gradients equal to zero, which are updated by
    :meth:`calc_LP_grad`. The object can
    be pickled. The deprecated property ``gradAtransij`` is an alias of
    ``gradAtsij``.

    Attributes
    ----------
    gradAij, gradBij, gradDij, gradAtsij : double[:, ::1]
        The shapes of these gradient matrices are:

            gradAij: (6, 5)
            gradBij: (6, 5)
            gradDij: (6, 5)
            gradAtsij: (3, 3)

        They contain the gradients of each laminate stiffness with respect to
        the thickness and respective lamination parameters. Column 0 has the
        dimension of the respective stiffness divided by a length and the
        other columns the dimension of the stiffness. The transverse
        shear terms are the constant-strain (uncorrected) ``Abar44``,
        ``Abar45``, ``Abar55``, as given by
        :func:`.laminate_from_LaminationParameters`. The rows and
        columns correspond to::

            gradAij:

                h xiA1 xiA2 xiA3 xiA4
            A11
            A12
            A16
            A22
            A26
            A66

            gradBij:

                h xiB1 xiB2 xiB3 xiB4
            B11
            B12
            B16
            B22
            B26
            B66

            gradDij:

                h xiD1 xiD2 xiD3 xiD4
            D11
            D12
            D16
            D22
            D26
            D66

            gradAtsij:

                   h xiAts1 xiAts2
            Abar44
            Abar45
            Abar55

    """
    def __cinit__(GradABD self):
        # NOTE in __cinit__, such that the memoryviews are also initialized
        #      when the object is created by __new__, e.g. when unpickling,
        #      since the module is compiled with initializedcheck=False
        self.gradAij = np.zeros((6, 5), dtype=DOUBLE)
        self.gradBij = np.zeros((6, 5), dtype=DOUBLE)
        self.gradDij = np.zeros((6, 5), dtype=DOUBLE)
        self.gradAtsij = np.zeros((3, 3), dtype=DOUBLE)

    def __reduce__(GradABD self):
        state = {}
        for name in _GRADABD_STATE:
            state[name] = np.asarray(getattr(self, name)).copy()
        return copyreg.__newobj__, (type(self), ), state

    def __setstate__(GradABD self, dict state):
        for name, value in state.items():
            # NOTE pickles created before 'gradAtransij' became 'gradAtsij'
            if name == 'gradAtransij':
                name = 'gradAtsij'
            setattr(self, name, np.ascontiguousarray(value, dtype=DOUBLE))

    @property
    def gradAtransij(self):
        r"""Deprecated, use :attr:`.gradAtsij` instead

        Issues a :class:`DeprecationWarning` and returns ``gradAtsij``.

        Returns
        -------
        gradAtsij : double[:, ::1]
            Gradients of shape ``(3, 3)``, see :class:`.GradABD`.

        """
        _warn_deprecated('GradABD.gradAtransij', 'GradABD.gradAtsij')
        return self.gradAtsij
    @gradAtransij.setter
    def gradAtransij(self, value):
        r"""Deprecated, use :attr:`.gradAtsij` instead

        Issues a :class:`DeprecationWarning` and sets ``gradAtsij``.

        Parameters
        ----------
        value : double[:, ::1]
            2-D C-contiguous array of float64, of shape ``(3, 3)``, see
            :class:`.GradABD`.

        """
        _warn_deprecated('GradABD.gradAtransij', 'GradABD.gradAtsij')
        self.gradAtsij = value

    cpdef void calc_LP_grad(GradABD self, double thickness, MatLamina mat, LaminationParameters lp,
            double offset=0.):
        r"""Gradients of the shell stiffnesses with respect to the thickness and
        lamination parameters

        Updates ``gradAij``, ``gradBij``, ``gradDij`` and ``gradAtsij`` with
        the derivatives of the stiffnesses computed by
        :func:`.laminate_from_LaminationParameters` from ``thickness``,
        from the invariants ``u1`` to ``u7`` of ``mat``, from the
        lamination parameters of ``lp`` and from ``offset``; column 0
        contains the derivatives with respect to the thickness, at constant
        lamination parameters and ``offset``. The entries ``[0, 2]``, ``[1, 1]`` and
        ``[2, 2]`` of ``gradAtsij``, which do not depend on the lamination
        parameters, are not updated. All six rows of ``gradAij``,
        ``gradBij`` and ``gradDij`` are updated, including row 5 (``A66``,
        ``B66``, ``D66``), whose derivatives with respect to `\xi_{A3}`,
        `\xi_{B3}` and `\xi_{D3}` are `-h u_3`, `-h^2/4 \, u_3` and `-h^3/12
        \, u_3`.

        Parameters
        ----------
        thickness : float
            The total thickness of the laminate, length
        mat : :class:`.MatLamina` object
            Material object
        lp : :class:`.LaminationParameters` object
            The container class with all lamination parameters already defined
        offset : float, optional
            Offset of the mid-surface from the reference surface, about which
            the lamination parameters are defined, see
            :func:`.laminate_from_LaminationParameters`. Length. Default ``0.``.

        Returns
        -------
        None
            The attributes of the object are updated.

        Notes
        -----
        The stiffnesses of :func:`.laminate_from_LaminationParameters` are
        `A = h \, U \{1, \xi_A\}`, `B = h e \, U \{1, 0\} + h^2/4 \, U \{0,
        \xi_B\}`, `D = (h^3/12 + h e^2) \, U \{1, 0\} + h^3/12 \, U \{0,
        \xi_D\}` and `\bar{A}_{ts} = h \, U_{ts} \{1, \xi_{A_{ts}}\}`, with
        `e` = ``offset``, `U` the rows of the in-plane terms of
        :meth:`.MatLamina.get_invariant_matrix` and `U_{ts}` those of the
        transverse shear terms. At constant lamination parameters and
        offset, `\partial A/\partial h = U \{1, \xi_A\}`, `\partial
        B/\partial h = e \, U \{1, 0\} + h/2 \, U \{0, \xi_B\}`, `\partial
        D/\partial h = (h^2/4 + e^2) \, U \{1, 0\} + h^2/4 \, U \{0,
        \xi_D\}` and `\partial \bar{A}_{ts}/\partial h = U_{ts} \{1,
        \xi_{A_{ts}}\}`, and at constant thickness the derivatives with respect to the
        lamination parameters are the columns of `U` multiplied by `h`,
        `h^2/4` and `h^3/12`, respectively, and those of `U_{ts}` by `h`.

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
        for i in range(6):
            for j in range(4):
                self.gradAij[i, j+1] = h*gradinv[i, j]

        # d(B11 B12 B16 B22 B26 B66) / dh
        self.gradBij[0, 0] = offset*mat.u1 + h/2.*(mat.u2*lp.xiB1 + 0*lp.xiB2 + mat.u3*lp.xiB3 + 0*lp.xiB4)
        self.gradBij[1, 0] = offset*mat.u4 + h/2.*(0*lp.xiB1 + 0*lp.xiB2 + (-1)*mat.u3*lp.xiB3 + 0*lp.xiB4)
        self.gradBij[2, 0] = h/2.*(0*lp.xiB1 + mat.u2/2.*lp.xiB2 + 0*lp.xiB3 + mat.u3*lp.xiB4)
        self.gradBij[3, 0] = offset*mat.u1 + h/2.*((-1)*mat.u2*lp.xiB1 + 0*lp.xiB2 + mat.u3*lp.xiB3 + 0*lp.xiB4)
        self.gradBij[4, 0] = h/2.*(0*lp.xiB1 + mat.u2/2.*lp.xiB2 + 0*lp.xiB3 + (-1)*mat.u3*lp.xiB4)
        self.gradBij[5, 0] = offset*mat.u5 + h/2.*(0*lp.xiB1 + 0*lp.xiB2 + (-1)*mat.u3*lp.xiB3 + 0*lp.xiB4)

        # d(B11 B12 B16 B22 B26 B66) / d(xiB1, xiB2, xiB3, xiB4)
        for i in range(6):
            for j in range(4):
                self.gradBij[i, j+1] = h*h/4.*gradinv[i, j]

        # d(D11 D12 D16 D22 D26 D66) / dh
        self.gradDij[0, 0] = offset*offset*mat.u1 + h*h/4.*(mat.u1 + mat.u2*lp.xiD1 + 0*lp.xiD2 + mat.u3*lp.xiD3 + 0*lp.xiD4)
        self.gradDij[1, 0] = offset*offset*mat.u4 + h*h/4.*(mat.u4 + 0*lp.xiD1 + 0*lp.xiD2 + (-1)*mat.u3*lp.xiD3 + 0*lp.xiD4)
        self.gradDij[2, 0] = h*h/4.*(0 + 0*lp.xiD1 + mat.u2/2.*lp.xiD2 + 0*lp.xiD3 + mat.u3*lp.xiD4)
        self.gradDij[3, 0] = offset*offset*mat.u1 + h*h/4.*(mat.u1 + (-1)*mat.u2*lp.xiD1 + 0*lp.xiD2 + mat.u3*lp.xiD3 + 0*lp.xiD4)
        self.gradDij[4, 0] = h*h/4.*(0 + 0*lp.xiD1 + mat.u2/2.*lp.xiD2 + 0*lp.xiD3 + (-1)*mat.u3*lp.xiD4)
        self.gradDij[5, 0] = offset*offset*mat.u5 + h*h/4.*(mat.u5 + 0*lp.xiD1 + 0*lp.xiD2 + (-1)*mat.u3*lp.xiD3 + 0*lp.xiD4)

        # d(D11 D12 D16 D22 D26 D66) / d(xiD1, xiD2, xiD3, xiD4)
        for i in range(6):
            for j in range(4):
                self.gradDij[i, j+1] = h*h*h/12.*gradinv[i, j]

        # d(A44, A45, A55) / dh
        self.gradAtsij[0, 0] = (mat.u6 + mat.u7*lp.xiAts1 + 0*lp.xiAts2)
        self.gradAtsij[1, 0] = (0 + 0*lp.xiAts1 + (-1)*mat.u7*lp.xiAts2)
        self.gradAtsij[2, 0] = (mat.u6 + (-1)*mat.u7*lp.xiAts1 + 0*lp.xiAts2)

        # d(A44, A45, A55) / d(xiAts1, xiAts2)
        self.gradAtsij[0, 1] = h*mat.u7
        self.gradAtsij[1, 2] = h*(-mat.u7)
        self.gradAtsij[2, 1] = h*(-mat.u7)


cpdef Laminate n_double_laminate(double thickness, int n, double[::1] angles_deg, MatLamina matlamina):
    r"""Create a N-double laminated plate using a faster code

    This code is considerably faster than :func:`composites.utils.n_double_laminate`.

    In a new :class:`.Laminate` object, the stiffnesses ``A11`` to ``A66`` and
    ``D11`` to ``D66`` are computed from ``thickness``, from the trace and
    the trace-normalized invariants of ``matlamina`` and from the first ``n``
    angles of ``angles_deg``. The object ``matlamina`` is modified in place
    by :meth:`.MatLamina.trace_normalize_plane_stress`. No plies are created,
    ``stack`` is empty, all the other stiffnesses are zero and
    :meth:`.Laminate.calc_equivalent_properties` is called.

    An N-double laminate consists of `[\pm\phi_1,\pm\phi_2, \cdots,
    \pm\phi_n]`.

    Parameters
    ----------
    thickness : float
        Total plate thickness, length.
    n : int
        Number of angle pairs,m defines the "N" in "N-double".
    angles_deg : 1-D array of float64
        List of `\phi_n` of the N-double laminate. C-contiguous array of size
        at least ``n``, of which the first ``n`` angles are used, in degrees.
    matlamina : MatLamina
        See :class:`.MatLamina` for details. Modified in place.

    Returns
    -------
    lam : :class:`.Laminate`
        The N-double laminate.

    See Also
    --------
    composites.utils.n_double_laminate, composites.utils.double_double_laminate

    Notes
    -----
    With the principle of homogenization, at the limit where many plies are
    used we have that `B=0`. Based on the double-double laminate as described
    by Shrivastava et al. [1]_. The plies at `\pm\phi_i` have the same
    thickness, such that the stiffnesses are the averages of the ply
    stiffnesses of :meth:`.MatLamina.get_invariant_matrix` over the `2n`
    angles, written with the trace `tr` and the trace-normalized invariants
    `u_i^*` of :meth:`.MatLamina.trace_normalize_plane_stress` [2]_:

    .. math::

        A_{11}^* = u_1^* + \frac{1}{n} \sum_{i=1}^{n} \left( u_2^* \cos
        2\phi_i + u_3^* \cos 4\phi_i \right), \quad
        A_{22}^* = u_1^* + \frac{1}{n} \sum_{i=1}^{n} \left( -u_2^* \cos
        2\phi_i + u_3^* \cos 4\phi_i \right) \\
        A_{12}^* = u_4^* - \frac{1}{n} \sum_{i=1}^{n} u_3^* \cos 4\phi_i,
        \quad
        A_{66}^* = u_5^* - \frac{1}{n} \sum_{i=1}^{n} u_3^* \cos 4\phi_i,
        \quad A_{16}^* = A_{26}^* = 0 \\
        A_{ij} = tr \, A_{ij}^* h, \qquad D_{ij} = tr \, A_{ij}^* h^3/12

    where the bending stiffnesses are those of the homogenized laminate,
    see :meth:`.Laminate.make_smeared`.

    References
    ----------
    .. [1] Shrivastava, S., Sharma, N., Tsai, S. W., and Mohite, P. M., D
       and DD-Drop Layup Optimization of Aircraft Wing Panels under
       Multi-Load Case Design Environment, Composite Structures, 248, 112518,
       2020.
    .. [2] Melo, J. D. D., Bi, J., and Tsai, S. W., A Novel Invariant-Based
       Design Approach to Carbon Fiber Reinforced Laminates, Composite
       Structures, 159, 44-52, 2017.

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
