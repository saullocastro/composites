r"""
=====================================================
Kassapoglou's methods (:mod:`composites.kassapoglou`)
=====================================================

.. currentmodule::composites.kassapoglou


Methods based on Kassapoglou's book.


Reference:

    Kassapoglou. Design and Analysis of Composite Structures. 2nd Edition. John Wiley & Sons Ltd, 2013.


"""
import numpy as np
from numpy import pi, tan, sqrt


def calc_Nxx_crit(a, b, m, n, D11, D12, D22, D66):
    r"""Calculate uniaxial compression buckling for a composite plate

    The output of this function is the result of Eq. 6.6, section 6.2 page 129.
    If `m` or `n` is set to ``None``, the function searchers for the critical number of
    half-waves in the corresponding direction, up to 10 half-waves.

    Reference:

        Kassapoglou. Design and Analysis of Composite Structures. 2nd Edition. John Wiley & Sons Ltd, 2013.

    Parameters
    ----------
    a, b : float
        Plate length and width.
    m, n : int or None
        Number of half-waves along the plate length and width, respectively.
    D11, D12, D22, D66 : float
        Terms of the D matrix.

    Result
    ------
    Nxx_crit : float
        Critical uniaxial compression buckling load.

    """
    if m is None and n is None:
        raise NotImplementedError("Only m or n can be None, not both")
    if m is None:
        m = np.arange(1, 11)
    if n is None:
        n = np.arange(1, 11)
    AR = a/b
    m = np.atleast_1d(m)
    n = np.atleast_1d(n)
    Nxx_crit = pi**2*(D11*m**4
                      + 2*(D12 + 2*D66)*m**2*n**2*AR**2
                      + D22*n**4*AR**4)/(a**2*m**2)
    return Nxx_crit.min()


def calc_Nxy_crit(a, D11, D12, D16, D22, D66, rtol=1e-5, atol=1e-6, max_iter=50):
    r"""Calculate shear buckling for a composite plate

    The output of this function is the result of Eq. 6.28, section 6.4 page
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
    rtol, atol : float
        Relative and absolute tolerances used to solve Eq. 6.30.
    max_iter : int
        Maximum number of iterations used in the Newton-Raphson scheme.

    Result
    ------
    Nxy_crit : float
        Critical shear buckling load.

    """
    alpha = np.pi/6
    for i in range(max_iter):
        AR = (D11/(D11*tan(alpha)**4 + 2*(D12 + 2*D66)*tan(alpha)**2 + D22))**(1/4)
        expr = 3*D11*AR**4*tan(alpha)**4 + (6*D11*AR**2 + 2*(D12 + 2*D66)*AR**4)*tan(alpha)**2 - (D11 + 2*(D12 + 2*D66)*AR**2 + D22*AR**4)
        if np.isclose(expr, 0, atol=1e-6, rtol=1e-5):
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


def calc_beff(b, Px, Pcr, A11, A12, A22):
    r"""Calculate the effective width of a plate

    The output of this function is the result of Eq. 7.15, section 7.1 page
    164. The effective width `b_{eff}` is defined as:

    `\int_y N_{xx} dy = 2{N_x}_{max} b_{eff}`


    Reference:

        Kassapoglou. Design and Analysis of Composite Structures. 2nd Edition. John Wiley & Sons Ltd, 2013.

    Parameters
    ----------
    b : float
        Plate width.
    Px : float
        Applied compressive force.
    Pcr : float
        Critical buckling force.
    A11, A12, A22: float
        Terms of the A matrix.

    Result
    ------
    beff : float
        Effective width of the plate.

    """
    beff = b/(2*(1 + 2*(1 + A12/A11)*(1 - Pcr/Px)*A11/(A11 + 3*A22)))
    return beff


def calc_Nxx_crit_combined_shear(k, a, b, D11, D12, D22, D66):
    r"""Calculate combined uniaxial-shear buckling load for a composite plate

    The output of this function is the result of Eq. 6.34, section 6.5 page
    142. This function calculates the critical `N_{xx}` buckling load under the
    current level of shear load given by `N_{xy} = k N_{xx}`.

    Reference:

        Kassapoglou. Design and Analysis of Composite Structures. 2nd Edition. John Wiley & Sons Ltd, 2013.

    Parameters
    ----------
    k : float
        Load ratio defined as `k=N_{xy}/N_{xx}`, where `N_{xy}` is the current
        shear load, and `N_{xx}` the current compression load.
    a, b : float
        Plate length and width.
    D11, D12, D22, D66 : float
        Terms of the D matrix.

    Result
    ------
    N0 : float
        Critical `N_{xx}` buckling load under the current level of shear load
        given by `N_{xy} = k N_{xx}`.

    """
    den = 2 - 8192/81*a**2*k**2/(b**2*pi**4)
    term = 9 + 65536/81*a**2*k**2/(b**2*pi**4)
    rhs_term1 = 5 + sqrt(term)
    rhs_term2 = 5 - sqrt(term)
    Nxx_crit1 = pi**2/a**2*(D11 + 2*(D12 + 2*D66)*a**2/b**2 + D22*a**4/b**4)/den*rhs_term1
    Nxx_crit2 = pi**2/a**2*(D11 + 2*(D12 + 2*D66)*a**2/b**2 + D22*a**4/b**4)/den*rhs_term2
    return min(abs(Nxx_crit1), abs(Nxx_crit2))


def calc_Nxx_crit_combined_shear_full(Nxy, a, b, D11, D12, D16, D22, D26, D66):
    r"""Calculate combined uniaxial-shear buckling load for a composite plate

    This solution does not ignore the D16 and D26 terms. Furthermore, this
    solution takes as input `N_{xy}` instead of the `k=N_{xy}/N_{xx}`
    originally used by Kassapoglou.

    Based on section 6.5. This function calculates the critical `N_{xx}` by
    `N_{xy}`.

    Reference:

        Kassapoglou. Design and Analysis of Composite Structures. 2nd Edition. John Wiley & Sons Ltd, 2013.

    Parameters
    ----------
    Nxy : float
        Shear load `N_{xy}`.
    a, b : float
        Plate length and width.
    D11, D12, D16, D22, D26, D66 : float
        All terms of the D matrix.

    Result
    ------
    Nxx_crit : float
        Critical `N_{xx}` buckling load under the current level of shear load
        given by `N_{xy}`.

    """
    a11 = pi**4*D11*b/(4*a**3) + pi**4*D12/(2*a*b) + pi**4*D22*a/(4*b**3) + pi**4*D66/(a*b)
    a12 = -160*pi**2*D16/(9*a**2) - 160*pi**2*D26/(9*b**2)
    a21 = -160*pi**2*D16/(9*a**2) - 160*pi**2*D26/(9*b**2)
    a22 = 4*pi**4*D11*b/a**3 + 8*pi**4*D12/(a*b) + 4*pi**4*D22*a/b**3 + 16*pi**4*D66/(a*b)
    term = 16384*Nxy**2 + 4608*Nxy*a12 + 4608*Nxy*a21 + 1296*a11**2 - 648*a11*a22 + 1296*a12*a21 + 81*a22**2
    Nxx_crit1 = a*(36*a11 + 9*a22 - sqrt(term))/(18*pi**2*b)
    Nxx_crit2 = a*(36*a11 + 9*a22 + sqrt(term))/(18*pi**2*b)
    return min(abs(Nxx_crit1), abs(Nxx_crit2))

