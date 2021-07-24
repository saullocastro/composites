"""
Composites Core Utils Module (:mod:`composites.utils`)
==============================================================

.. currentmodule:: composites.utils

"""
import numpy as np
from numpy import cos, sin

from .core import (MatLamina, Lamina, Laminate, LaminationParameters,
        laminate_from_lamination_parameters)

def read_laminaprop(laminaprop, rho=0):
    """Returns a :class:`.MatLamina` object based on an input ``laminaprop`` tuple

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

            laminaprop = (e1, e2, nu12) # legacy, kept for compatibility with old codes

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
    assert len(laminaprop) in (2, 3, 6, 9), ('Invalid entry for laminaprop: ' +
                                             str(laminaprop))
    if len(laminaprop) == 3: #ISOTROPIC in-plane stress legacy
        e = laminaprop[0]
        nu = laminaprop[2]
        g = e/(2*(1+nu))
        laminaprop = (e, e, nu, g, g, g, 0, 0, 0)
    elif len(laminaprop) == 2: #ISOTROPIC in-plane stress new
        e = laminaprop[0]
        nu = laminaprop[1]
        g = e/(2*(1+nu))
        laminaprop = (e, e, nu, g, g, g, 0, 0, 0)
    elif len(laminaprop) == 6: #ORTHOTROPIC in-plane stress
        laminaprop = tuple(list(laminaprop) + [0, 0, 0])
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


def laminated_plate(stack, plyt=None, laminaprop=None, rho=0., plyts=None,
        laminaprops=None, rhos=None, offset=0., calc_scf=True):
    """Read a laminate stacking sequence data.

    :class:`.Laminate` object is returned based on the inputs given.

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
        If True, use :func:`.Laminate.calc_scf` to compute shear correction
        factors, otherwise the default value of 5/6 is used

    Notes
    -----
    ``plyt`` or ``plyts`` must be supplied
    ``laminaprop`` or ``laminaprops`` must be supplied

    For orthotropic plies, the ``laminaprop`` should be::

        laminaprop = (E11, E22, nu12, G12, G13, G23)

    For isotropic plies, the ``laminaprop`` should be::

        laminaprop = (E, nu)

    """
    lam = Laminate()
    lam.offset = offset
    lam.stack = list(stack)

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
    if calc_scf:
        lam.calc_scf()

    return lam


def isotropic_plate(thickness, E, nu, offset=0., calc_scf=True, rho=0.):
    """Read data for an isotropic plate

    :class:`.Laminate` object is returned based on the inputs given.

    Parameters
    ----------
    thickness : float
        Plate thickness.
    E : float
        Young modulus.
    nu : float, optional
        Poisson's ratio.
    rho : float, optional
        Material density
    offset : float, optional
        Offset along the normal axis about the mid-surface, which influences
        the extension-bending coupling (B matrix).
    calc_scf : bool, optional
        If True, use :func:`.Laminate.calc_scf` to compute shear correction
        factors, otherwise the default value of 5/6 is used.

    """
    return laminated_plate(plyt=thickness, stack=[0], laminaprop=(E, nu),
            rho=rho, offset=offset, calc_scf=calc_scf)

def double_double_plate(thickness, phideg, psideg, laminaprop=None,
        rho=0., calc_scf=True):
    """Create a double-double laminated plate

    A double-double (DD) laminate consists of :math:`[\pm\phi,\pm\psi]`, with
    ``phideg=`` :math:`\phi`, and ``psideg=`` :math:`\psi`. With the
    principle of homogenization, at the limit where many plies are used
    we have that :math:`B=0`. Reference:

        [1] Shrivastava, S., Sharma, N., Tsai, S. W., and Mohite, P. M., 2020,
        “D and DD-Drop Layup Optimization of Aircraft Wing Panels under
        Multi-Load Case Design Environment,” Compos. Struct., 248(January), p.
        112518.

    Parameters
    ----------
    thickness : float
        Total plate thickness.
    phideg : float
        Angle :math:`\psi` of the DD laminate.
    psideg : float
        Angle :math:`\phi` of the DD laminate.
    laminaprop : tuple
        See :func:`.read_laminaprop` for details.
    rho : float, optional
        Material density
    calc_scf : bool, optional
        If True, use :func:`.Laminate.calc_scf` to compute shear correction
        factors, otherwise the default value of 5/6 is used.

    """
    m = read_laminaprop(laminaprop, rho)
    tr = m.q11 + m.q22 + 2*m.q66
    m.trace_normalize_plane_stress()
    lam = Laminate()
    lam.h = thickness
    phi = np.deg2rad(phideg)
    psi = np.deg2rad(psideg)
    A11_star = m.u1 + m.u2*cos(2*phi)*cos(2*psi) + m.u3*cos(4*phi)*cos(4*psi)
    A12_star = m.u4 - m.u3*cos(4*phi)*cos(4*psi)
    A16_star = 0 #m.u2/2.*sin(2*phi)*sin(2*psi) + m.u3*sin(4*phi)*sin(4*psi)
    A22_star = m.u1 - m.u2*cos(2*phi)*cos(2*psi) + m.u3*cos(4*phi)*cos(4*psi)
    A26_star = 0 #m.u2/2.*sin(2*phi)*sin(2*psi) - m.u3*sin(4*phi)*sin(4*psi)
    A66_star = m.u5 - m.u3*cos(4*phi)*cos(4*psi)
    lam.A11 = tr*A11_star*lam.h
    lam.A12 = tr*A12_star*lam.h
    lam.A16 = tr*A16_star*lam.h
    lam.A22 = tr*A22_star*lam.h
    lam.A26 = tr*A26_star*lam.h
    lam.A66 = tr*A66_star*lam.h
    lam.D11 = tr*A11_star*lam.h**3/12
    lam.D12 = tr*A12_star*lam.h**3/12
    lam.D16 = tr*A16_star*lam.h**3/12
    lam.D22 = tr*A22_star*lam.h**3/12
    lam.D26 = tr*A26_star*lam.h**3/12
    lam.D66 = tr*A66_star*lam.h**3/12

    return lam


