r"""
===========================================================================
Implementations based on Kassapoglou's book (:mod:`composites.kassapoglou`)
===========================================================================

.. currentmodule::composites.kassapoglou

"""
import numpy as np
from numpy import pi, tan


def calc_Nxx_crit(a, b, m, n, D11, D12, D22, D66):
    r"""Calculate uniaxial compression buckling for a composite plate

    The output of this function is the result of Eq. 6.6, section 6.2 page 129.

    Reference:

        Kassapoglou. Design and Analysis of Composite Structures. 2nd Edition. John Wiley & Sons Ltd, 2013.

    Parameters
    ----------
    a, b : float
        Plate length and width.
    m, n : int
        Number of half-waves along the plate length and width, respectively.
    D11, D12, D22, D66 : float
        Terms of the D matrix.

    Result
    ------
    Nxx_crit : float
        Critical uniaxial compression buckling load.

    """
    AR = a/b
    Nxx_crit = pi**2*(D11*m**4
                      + 2*(D12 + 2*D66)*m**2*n**2*AR**2
                      + D22*n**4*AR**4)/(a**2*m**2)
    return Nxx_crit


def calc_Nxy_crit(a, D11, D12, D16, D22, D66):
    r"""Calculate shear buckling for a composite plate

    The output of this function is the result of Eq. 6.28, section 6.3 page
    137. Variables `AR` and `\alpha` are solved using a Newton-Raphson scheme
    that finds the solution of Eqs. 6.29 and 6.30 simultaneously.

    Reference:

        Kassapoglou. Design and Analysis of Composite Structures. 2nd Edition. John Wiley & Sons Ltd, 2013.

    Parameters
    ----------
    a : float
        Plate length.
    D11, D12, D16, D22, D66 : float
        Terms of the D matrix.

    Result
    ------
    Nxy_crit : float
        Critical shear buckling load.

    """
    alpha = np.pi/6
    max_iter = 50
    for i in range(max_iter):
        AR = (D11/(D11*tan(alpha)**4 + 2*(D12 + 2*D66)*tan(alpha)**2 + D22))**(1/4)
        expr = 3*D11*AR**4*tan(alpha)**4 + (6*D11*AR**2 + 2*(D12 + 2*D66)*AR**4)*tan(alpha)**2 - (D11 + 2*(D12 + 2*D66)*AR**2 + D22*AR**4)
        if np.isclose(expr, 0):
            break
        dexpr_dalpha = 3*AR**4*D11*(4*tan(alpha)**2 + 4)*tan(alpha)**3 + (AR**4*(2*D12 + 4*D66) + 6*AR**2*D11)*(2*tan(alpha)**2 + 2)*tan(alpha)
        dalpha = -expr/dexpr_dalpha
        alpha += dalpha
    else:
        raise RuntimeError('Newton-Raphson scheme did not converge!')
    Nxy_crit = pi**2/(2*AR**2*a**2*tan(alpha))*(
        D11*(1 + 6*tan(alpha)**2*AR**2 + tan(alpha)**4*AR**4)
        + 2*(D12 + 2*D66)*(AR**2 + AR**4*tan(alpha)**2)
        + D22*AR**4)
    return Nxy_crit

