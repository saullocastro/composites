"""
Composites Core Utils Module (:mod:`composites.utils`)
==============================================================

.. currentmodule:: composites.utils

Functions that create a :class:`.Laminate` from the stacking sequence and the
ply data (:func:`.laminated_plate`), for an isotropic plate
(:func:`.isotropic_plate`) or for double-double and N-double laminates
(:func:`.double_double_laminate`, :func:`.n_double_laminate`, also available
as ``double_double_plate`` and ``n_double_plate``), and a :class:`.MatLamina`
from a tuple of material properties (:func:`.read_laminaprop`).

"""
import warnings

import numpy as np
from numpy import cos, sin, deg2rad

from .core import (MatLamina, Lamina, Laminate, LaminationParameters,
        laminate_from_lamination_parameters)


def read_laminaprop(laminaprop, rho=0):
    r"""Returns a :class:`.MatLamina` object based on an input ``laminaprop`` tuple

    The entries of ``laminaprop`` and ``rho`` are stored in a new
    :class:`.MatLamina` object, the Poisson's ratios ``nu21``, ``nu31`` and
    ``nu32`` are computed from ``nu12``, ``nu13``, ``nu23`` and the Young
    moduli, and :meth:`.MatLamina.rebuild` is called. With 2 entries, ``g12 =
    g13 = g23 = e/(2*(1+nu))``; with 2 or 6 entries, ``e3``, ``nu13`` and
    ``nu23`` are zero. Another number of entries fails an ``assert``
    statement.

    Parameters
    ----------
    laminaprop : list or tuple
        For the most general case of tri-axial stress, use a tuple containing
        the folliwing entries::

            laminaprop = (e1, e2, nu12, g12, g13, g23, e3, nu13, nu23)

        For isotropic materials aiming calculations with tri-axial stresses,
        use::

            g = e/(2*(1+nu))
            laminaprop = (e, e, nu, g, g, g, e, nu, nu)

        For othotropic materials with in-plane stresses the user can only
        supply::

            laminaprop = (e1, e2, nu12, g12, g13, g23)

        For isotropic materials with in-plane stresses the user can only
        supply::

            laminaprop = (e, nu) # new


        ======  ==============================
        symbol  value
        ======  ==============================
        e1      Young Module in direction 1
        e2      Young Module in direction 2
        nu12    12 Poisson's ratio
        g12     12 Shear Modulus
        g13     13 Shear Modulus
        g23     23 Shear Modulus
        e3      Young Module in direction 3
        nu13    13 Poisson's ratio
        nu23    23 Poisson's ratio
        ======  ==============================

        The Young and shear moduli have the dimension of force per unit area,
        the Poisson's ratios are dimensionless.


    rho : float, optional
        Material density, mass per unit volume. Default ``0``.


    Returns
    -------
    matlam : MatLamina
        A :class:`.MatLamina` object.

    Notes
    -----
    The Poisson's ratios not given follow from the reciprocity relations of
    an orthotropic material, `\nu_{ij}/E_i = \nu_{ji}/E_j`, which make its
    compliance matrix symmetric [1]_:

    .. math::

        \nu_{21} = \nu_{12} \frac{E_2}{E_1}, \quad
        \nu_{31} = \nu_{13} \frac{E_3}{E_1}, \quad
        \nu_{32} = \nu_{23} \frac{E_3}{E_2}

    and for an isotropic material `G = E/(2(1 + \nu))`. With 2 or 6 entries,
    ``e3 = 0`` makes the constitutive constants of :meth:`.MatLamina.rebuild`
    the plane-stress reduced stiffnesses.

    References
    ----------
    .. [1] Reddy, J. N., Mechanics of Laminated Composite Plates and Shells:
       Theory and Analysis, Second Edition, CRC Press, 2004.

    """
    matlam = MatLamina()

    #laminaProp = (e1, e2, nu12, g12, g13, g23, e3, nu13, nu23)
    assert len(laminaprop) in (2, 6, 9), ('Invalid entry for laminaprop: ' +
                                          str(laminaprop))
    if len(laminaprop) == 2: #ISOTROPIC in-plane stress new
        e = laminaprop[0]
        nu = laminaprop[1]
        g = e/(2*(1+nu))
        laminaprop = (e, e, nu, g, g, g, 0, 0, 0)
    elif len(laminaprop) == 6: #ORTHOTROPIC in-plane stress
        laminaprop = tuple(list(laminaprop) + [0, 0, 0])
    elif len(laminaprop) == 9: #ORTHOTROPIC 3D stress state
        laminaprop = tuple(laminaprop)

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


def _shear_correction_from_calc_scf(calc_scf, shear_correction):
    if calc_scf is None:
        return shear_correction
    warnings.warn("The argument 'calc_scf' is deprecated, use "
                  "'shear_correction' instead, with calc_scf=True mapped to "
                  "shear_correction='rohwer' and calc_scf=False mapped to "
                  "shear_correction=None", DeprecationWarning, stacklevel=3)
    return 'rohwer' if calc_scf else None


def laminated_plate(stack, plyt=None, laminaprop=None, rho=0., plyts=None,
        laminaprops=None, rhos=None, offset=0., calc_scf=None,
        shear_correction='rohwer'):
    r"""Read a laminate stacking sequence data.

    :class:`.Laminate` object is returned based on the inputs given.

    One :class:`.Lamina` is created for each entry of ``stack``, with its
    thickness, its :class:`.MatLamina` given by :func:`.read_laminaprop` with
    its ``laminaprop`` and density, and :meth:`.Lamina.rebuild` called. Then
    :meth:`.Laminate.calc_constitutive_matrix` and
    :meth:`.Laminate.calc_equivalent_properties` are called.

    Parameters
    ----------
    stack : list
        Angles of the stacking sequence in degrees, from the bottom to the top
        ply.
    plyt : float, optional
        When all plies have the same thickness, ``plyt`` can be supplied.
        ``plyt`` or ``plyts`` must be supplied. Length. Default ``None``.
    laminaprop : tuple, optional
        When all plies have the same material properties, ``laminaprop``
        can be supplied. ``laminaprop`` or ``laminaprops`` must be supplied.
        For orthotropic plies, the ``laminaprop`` should be::

            laminaprop = (E11, E22, nu12, G12, G13, G23)

        For isotropic plies, the ``laminaprop`` should be::

            laminaprop = (E, nu)

        See :func:`.read_laminaprop`. Default ``None``.
    rho : float, optional
        Uniform material density to be used for all plies, mass per unit
        volume. Default ``0.``.
    plyts : list, optional
        A list of floats with the thickness of each ply. When given,
        ``plyt`` is ignored. Length. Default ``None``.
    laminaprops : list, optional
        A list of tuples with a laminaprop for each ply. When given,
        ``laminaprop`` is ignored. Default ``None``.
    rhos : list, optional
        A list of floats with the material density of each ply. When given,
        ``rho`` is ignored. Mass per unit volume. Default ``None``.
    offset : float, optional
        Offset along the normal axis about the mid-surface, which influences
        the laminate properties. Length. Default ``0.``.
    calc_scf : bool, optional
        Deprecated, use ``shear_correction`` instead. When given, ``True`` is
        mapped to ``shear_correction='rohwer'`` and ``False`` to
        ``shear_correction=None``, overriding ``shear_correction``, and a
        :class:`DeprecationWarning` is issued. Default ``None``.
    shear_correction : str or None, optional
        Method used to compute the transverse shear stiffnesses ``A44``,
        ``A45``, ``A55``, which are returned with the correction already
        applied: ``'rohwer'`` (default, equilibrium approach),
        ``'vlachoutsis'``, ``'whitney'`` (same as ``'vlachoutsis'``),
        ``'chow'`` (symmetric laminates), ``'birman_bert'`` (average shear
        strain), ``'thickness_shear'`` (dynamic, requires ``rho``),
        ``'constant'`` (5/6) or ``None`` (no correction). See
        :meth:`.Laminate.calc_transverse_shear_stiffness`.

    Returns
    -------
    lam : :class:`.Laminate`
        The laminate, with its stiffnesses and equivalent properties
        computed.

    Raises
    ------
    ValueError
        If neither ``plyt`` nor ``plyts`` is supplied, or if neither
        ``laminaprop`` nor ``laminaprops`` is supplied.

    """
    lam = Laminate()
    lam.offset = offset
    lam.stack = list(stack)
    lam.shear_correction = _shear_correction_from_calc_scf(calc_scf,
            shear_correction)

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

    plies = []
    lam.h = 0.
    for plyt, laminaprop, thetadeg, rho in zip(plyts, laminaprops, stack, rhos):
        laminaprop = laminaprop
        ply = Lamina()
        ply.thetadeg = float(thetadeg)
        ply.h = plyt
        lam.h += ply.h
        ply.matlamina = read_laminaprop(laminaprop, rho)
        ply.rebuild()
        plies.append(ply)
    lam.plies = plies

    lam.calc_constitutive_matrix()
    lam.calc_equivalent_properties()

    return lam


def isotropic_plate(thickness, E, nu, offset=0., calc_scf=None, rho=0.,
        shear_correction='rohwer'):
    r"""Read data for an isotropic plate

    :class:`.Laminate` object is returned based on the inputs given.

    Calls :func:`.laminated_plate` with a single ply of angle ``0``,
    thickness ``thickness`` and ``laminaprop=(E, nu)``.

    Parameters
    ----------
    thickness : float
        Plate thickness, length.
    E : float
        Young modulus, force per unit area.
    nu : float
        Poisson's ratio, dimensionless.
    rho : float, optional
        Material density, mass per unit volume. Default ``0.``.
    offset : float, optional
        Offset along the normal axis about the mid-surface, which influences
        the extension-bending coupling (B matrix). Length. Default ``0.``.
    calc_scf : bool, optional
        Deprecated, use ``shear_correction`` instead. When given, ``True`` is
        mapped to ``shear_correction='rohwer'`` and ``False`` to
        ``shear_correction=None``, overriding ``shear_correction``, and a
        :class:`DeprecationWarning` is issued. Default ``None``.
    shear_correction : str or None, optional
        See :func:`.laminated_plate`. For an isotropic plate ``'rohwer'``,
        ``'vlachoutsis'``, ``'whitney'``, ``'chow'`` and ``'constant'`` give
        ``A44 = A55 = 5/6 G h``,
        ``'birman_bert'`` gives ``G h`` and ``'thickness_shear'`` gives
        ``pi^2/12 G h``, the latter requiring ``rho > 0``. Default
        ``'rohwer'``.

    Returns
    -------
    lam : :class:`.Laminate`
        The plate, with its stiffnesses and equivalent properties computed.

    See Also
    --------
    laminated_plate

    """
    shear_correction = _shear_correction_from_calc_scf(calc_scf,
            shear_correction)
    return laminated_plate(plyt=thickness, stack=[0], laminaprop=(E, nu),
            rho=rho, offset=offset, shear_correction=shear_correction)


def double_double_laminate(thickness, phideg, psideg, laminaprop=None, rho=0.):
    r"""Create a double-double laminate

    A double-double (DD) laminate consists of `[\pm\phi,\pm\psi]`, with
    ``phideg=`` `\phi`, and ``psideg=`` `\psi`.

    In a new :class:`.Laminate` object, the stiffnesses ``A11`` to ``A66`` and
    ``D11`` to ``D66`` are computed from ``thickness``, from the trace and the
    trace-normalized invariants of the material given by ``laminaprop`` and
    from the two angles. No plies are created, ``stack`` is ``[phideg,
    psideg]``, all the other stiffnesses are zero and
    :meth:`.Laminate.calc_equivalent_properties` is called.

    Parameters
    ----------
    thickness : float
        Total plate thickness, length.
    phideg : float
        Angle `\phi` of the DD laminate, in degrees.
    psideg : float
        Angle `\psi` of the DD laminate, in degrees.
    laminaprop : tuple
        See :func:`.read_laminaprop` for details. The default ``None`` is not
        accepted by :func:`.read_laminaprop`, such that a value must be given.
    rho : float, optional
        Material density, mass per unit volume. Default ``0.``.

    Returns
    -------
    lam : :class:`.Laminate`
        The double-double laminate.

    See Also
    --------
    n_double_laminate

    Notes
    -----
    With the principle of homogenization, at the limit where many plies are
    used we have that `B=0`, see Shrivastava et al. [1]_. The four plies
    `\pm\phi`, `\pm\psi` have the same thickness, such that the stiffnesses
    are the averages of the ply stiffnesses of
    :meth:`.MatLamina.get_invariant_matrix` over the four angles, written with
    the trace `tr` and the trace-normalized invariants `u_i^*` of
    :meth:`.MatLamina.trace_normalize_plane_stress` [2]_. With `\alpha =
    (\phi + \psi)/2` and `\beta = (\phi - \psi)/2`, the averages of `\cos
    2\theta` and `\cos 4\theta` are `\cos 2\alpha \cos 2\beta` and `\cos
    4\alpha \cos 4\beta`, and those of `\sin 2\theta` and `\sin 4\theta`
    vanish:

    .. math::

        A_{11}^* = u_1^* + u_2^* \cos 2\alpha \cos 2\beta + u_3^* \cos
        4\alpha \cos 4\beta, \quad
        A_{22}^* = u_1^* - u_2^* \cos 2\alpha \cos 2\beta + u_3^* \cos
        4\alpha \cos 4\beta \\
        A_{12}^* = u_4^* - u_3^* \cos 4\alpha \cos 4\beta, \quad
        A_{66}^* = u_5^* - u_3^* \cos 4\alpha \cos 4\beta, \quad
        A_{16}^* = A_{26}^* = 0 \\
        A_{ij} = tr \, A_{ij}^* h, \qquad D_{ij} = tr \, A_{ij}^* h^3/12

    where the bending stiffnesses are those of the homogenized laminate, see
    :meth:`.Laminate.make_smeared`.

    References
    ----------
    .. [1] Shrivastava, S., Sharma, N., Tsai, S. W., and Mohite, P. M., D and
       DD-Drop Layup Optimization of Aircraft Wing Panels under Multi-Load
       Case Design Environment, Composite Structures, 248, 112518, 2020.
    .. [2] Melo, J. D. D., Bi, J., and Tsai, S. W., A Novel Invariant-Based
       Design Approach to Carbon Fiber Reinforced Laminates, Composite
       Structures, 159, 44-52, 2017.

    """
    m = read_laminaprop(laminaprop, rho)
    tr = m.q11 + m.q22 + 2*m.q66
    m.trace_normalize_plane_stress()
    lam = Laminate()
    lam.h = thickness
    phi = deg2rad(phideg + psideg)/2.
    psi = deg2rad(phideg - psideg)/2.
    A11_star = m.u1 + m.u2*cos(2*phi)*cos(2*psi) + m.u3*cos(4*phi)*cos(4*psi)
    A12_star = m.u4 - m.u3*cos(4*phi)*cos(4*psi)
    A16_star = (m.u2/2.*(+sin(2*deg2rad(phideg))
                         +sin(-2*deg2rad(phideg))
                         +sin(2*deg2rad(psideg))
                         +sin(-2*deg2rad(psideg)))
                + m.u3*(+sin(4*deg2rad(phideg))
                        +sin(-4*deg2rad(phideg))
                        +sin(4*deg2rad(psideg))
                        +sin(-4*deg2rad(psideg))
                    ))
    A22_star = m.u1 - m.u2*cos(2*phi)*cos(2*psi) + m.u3*cos(4*phi)*cos(4*psi)
    A26_star = (m.u2/2.*(+sin(2*deg2rad(phideg))
                         +sin(-2*deg2rad(phideg))
                         +sin(2*deg2rad(psideg))
                         +sin(-2*deg2rad(psideg)))
                - m.u3*(+sin(4*deg2rad(phideg))
                        +sin(-4*deg2rad(phideg))
                        +sin(4*deg2rad(psideg))
                        +sin(-4*deg2rad(psideg))
                    ))
    A66_star = m.u5 - m.u3*cos(4*phi)*cos(4*psi)
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
    lam.stack = [phideg, psideg]

    return lam


def n_double_laminate(thickness, angles_deg, laminaprop=None, rho=0.):
    r"""Create a N-double laminate

    This code is considerably slower than :func:`composites.core.n_double_laminate`.

    An N-double laminate consists of `[\pm\phi_1,\pm\phi_2, \cdots,
    \pm\phi_n]`.

    In a new :class:`.Laminate` object, the stiffnesses ``A11`` to ``A66`` and
    ``D11`` to ``D66`` are computed from ``thickness``, from the trace and the
    trace-normalized invariants of the material given by ``laminaprop`` and
    from the angles. No plies are created, ``stack`` is ``angles_deg``, all
    the other stiffnesses are zero and
    :meth:`.Laminate.calc_equivalent_properties` is called.

    Parameters
    ----------
    thickness : float
        Total plate thickness, length.
    angles_deg : list of float
        List of `\phi_n` of the N-double laminate, in degrees.
    laminaprop : tuple
        See :func:`.read_laminaprop` for details. The default ``None`` is not
        accepted by :func:`.read_laminaprop`, such that a value must be given.
    rho : float, optional
        Material density, mass per unit volume. Default ``0.``.

    Returns
    -------
    lam : :class:`.Laminate`
        The N-double laminate.

    See Also
    --------
    composites.core.n_double_laminate, double_double_laminate

    Notes
    -----
    With the principle of homogenization, at the limit where many plies are
    used we have that `B=0`. Based on the double-double laminate as described
    by Shrivastava et al. [1]_. The stiffnesses are the same as those of
    :func:`composites.core.n_double_laminate`, with `n` the number of angles
    of ``angles_deg``.

    References
    ----------
    .. [1] Shrivastava, S., Sharma, N., Tsai, S. W., and Mohite, P. M., D and
       DD-Drop Layup Optimization of Aircraft Wing Panels under Multi-Load
       Case Design Environment, Composite Structures, 248, 112518, 2020.

    """
    m = read_laminaprop(laminaprop, rho)
    tr = m.q11 + m.q22 + 2*m.q66
    m.trace_normalize_plane_stress()
    lam = Laminate()
    lam.h = thickness
    A11_star = m.u1
    A12_star = m.u4
    A22_star = m.u1
    A16_star = 0
    A26_star = 0
    A66_star = m.u5
    for angle_deg in angles_deg:
        angle = deg2rad(angle_deg)
        A11_star += m.u2*(cos(2*angle) + cos(-2*angle))/(2*len(angles_deg))
        A11_star += m.u3*(cos(4*angle) + cos(-4*angle))/(2*len(angles_deg))

        A12_star -= m.u3*(cos(4*angle) + cos(-4*angle))/(2*len(angles_deg))

        A22_star -= m.u2*(cos(2*angle) + cos(-2*angle))/(2*len(angles_deg))
        A22_star += m.u3*(cos(4*angle) + cos(-4*angle))/(2*len(angles_deg))

        A66_star -= m.u3*(cos(4*angle) + cos(-4*angle))/(2*len(angles_deg))

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
    lam.stack = angles_deg

    return lam


# NOTE aliases for back-compatibility
double_double_plate = double_double_laminate
n_double_plate = n_double_laminate
