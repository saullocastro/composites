r"""
=====================================================
Kassapoglou's methods (:mod:`composites.kassapoglou`)
=====================================================

.. currentmodule:: composites.kassapoglou


Methods based on Kassapoglou's book: buckling loads of rectangular plates
under compression, shear and combined compression and shear, and the
effective width of a buckled plate, see the functions of this module.

Notes
-----
The functions implement the equations of Kassapoglou [1]_, chapters 6 and
7, whose numbers are given in each function.

References
----------
.. [1] Kassapoglou, C., Design and Analysis of Composite Structures, Second
   Edition, John Wiley & Sons, 2013.

"""
import numpy as np
from numpy import pi, tan, sqrt


def calc_Nxx_crit(a, b, m, n, D11, D12, D22, D66):
    r"""Calculate uniaxial compression buckling for a composite plate

    If `m` or `n` is set to ``None``, the function searchers for the critical number of
    half-waves in the corresponding direction, up to 10 half-waves.

    The load is evaluated for every value of `m` and `n`, the ``None`` one
    replaced by ``1, 2, ..., 10``, and the minimum is returned. `m` and `n`
    may also be given as arrays, which are broadcast against each other.

    Parameters
    ----------
    a, b : float
        Plate length and width, dimension length.
    m, n : int or None
        Number of half-waves along the plate length and width, respectively.
    D11, D12, D22, D66 : float
        Terms of the D matrix, dimension force times length.

    Returns
    -------
    Nxx_crit : float
        Critical uniaxial compression buckling load, force per unit length.

    Raises
    ------
    NotImplementedError
        If both `m` and `n` are ``None``.

    Notes
    -----
    The output of this function is the result of Eq. 6.6, section 6.2 page
    129 of [1]_, with `AR = a/b`:

    .. math::

        N_{xx,crit} = \frac{\pi^2}{a^2 m^2} \left( D_{11} m^4 + 2 (D_{12} +
        2 D_{66}) m^2 n^2 AR^2 + D_{22} n^4 AR^4 \right)

    i.e. the buckling load of a rectangular plate with the buckling mode
    `\sin(m \pi x/a) \sin(n \pi y/b)`, exact for a specially orthotropic
    plate, `D_{16} = D_{26} = 0`, simply supported on its four edges.

    References
    ----------
    .. [1] Kassapoglou, C., Design and Analysis of Composite Structures,
       Second Edition, John Wiley & Sons, 2013.

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

    Variables `AR` and `\alpha` of the buckling pattern are solved using a
    Newton-Raphson scheme. The iterations start from `\alpha = \pi/6` and
    stop when the residual is close to zero according to
    :func:`numpy.isclose`. `D16` is not used.

    Parameters
    ----------
    a : float
        Width of the long plate, i.e. its short dimension, along `x`,
        dimension length.
    D11, D12, D16, D22, D66 : float
        Terms of the D matrix, dimension force times length.
    rtol, atol : float, optional
        Relative and absolute tolerances intended for the residual of the
        Newton-Raphson scheme. Default ``1e-5`` and ``1e-6``. They are not
        used: the convergence check uses the fixed values ``atol=1e-6`` and
        ``rtol=1e-5``.
    max_iter : int, optional
        Maximum number of iterations used in the Newton-Raphson scheme.
        Default ``50``.

    Returns
    -------
    Nxy_crit : float
        Critical shear buckling load, force per unit length.

    Raises
    ------
    RuntimeError
        If the Newton-Raphson scheme does not converge within `max_iter`
        iterations.

    Notes
    -----
    The output of this function is the result of Eq. 6.28, section 6.4 page
    137 of [1]_, for a long rectangular plate of width `a`, simply supported
    on its long edges `x = 0` and `x = a`, made of a symmetric laminate with
    `D_{16} = D_{26} = 0`, under the shear load `N_{xy}` (pp. 136-137). The
    buckling mode, Eq. 6.20,

    .. math::

        w = w_o \sin\frac{\pi x}{a} \sin\frac{\pi (y - x \tan\alpha)}{L}

    is zero along the long edges and along lines inclined by the angle
    `\alpha` to the `x` axis, separated by the distance `L` along the long
    direction `y` (Fig. 6.10), such that the length of the plate does not
    enter, and `AR = a/L` (p. 137). Eq. 6.28 follows from the minimum of the
    total potential energy with respect to `w_o`, Eqs. 6.21-6.27, and the
    critical `AR` and `\alpha` minimize it, Eqs. 6.29 and 6.30 (p. 138).
    Variables `AR` and `\alpha` are solved using a Newton-Raphson scheme
    that finds the solution of Eqs. 6.29 and 6.30 simultaneously:

    .. math::

        AR = \left( \frac{D_{11}}{D_{11} \tan^4\alpha + 2 (D_{12} + 2
        D_{66}) \tan^2\alpha + D_{22}} \right)^{1/4} \\
        r(\alpha) = 3 D_{11} AR^4 \tan^4\alpha + (6 D_{11} AR^2 + 2 (D_{12}
        + 2 D_{66}) AR^4) \tan^2\alpha - (D_{11} + 2 (D_{12} + 2 D_{66})
        AR^2 + D_{22} AR^4) = 0 \\
        N_{xy,crit} = \frac{\pi^2}{2 AR^2 a^2 \tan\alpha} \left( D_{11} (1
        + 6 \tan^2\alpha AR^2 + \tan^4\alpha AR^4) + 2 (D_{12} + 2 D_{66})
        (AR^2 + AR^4 \tan^2\alpha) + D_{22} AR^4 \right)

    The tolerances ``rtol`` and ``atol`` are intended for the residual `r`
    of Eq. 6.30. At each iteration `AR` is updated with Eq. 6.29 and
    `\alpha` is corrected by `-r/(\partial r/\partial \alpha)`, the
    derivative being taken at constant `AR`.

    References
    ----------
    .. [1] Kassapoglou, C., Design and Analysis of Composite Structures,
       Second Edition, John Wiley & Sons, 2013.

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

    The effective width `b_{eff}` is defined as:

    `\int_y N_{xx} dy = 2{N_x}_{max} b_{eff}`

    `Px` and `Pcr` enter only through the ratio ``Pcr/Px``, and the terms
    of the A matrix only through ratios between them, such that the result
    has the dimension of `b`.

    Parameters
    ----------
    b : float
        Plate width, dimension length.
    Px : float
        Applied compressive force, of the same dimension as `Pcr`.
    Pcr : float
        Critical buckling force, of the same dimension as `Px`.
    A11, A12, A22 : float
        Terms of the A matrix, dimension force per unit length.

    Returns
    -------
    beff : float
        Effective width of the plate, dimension length.

    Notes
    -----
    The output of this function is the result of Eq. 7.15, section 7.1 page
    164 of [1]_:

    .. math::

        b_{eff} = \frac{b}{2 \left( 1 + 2 \left( 1 + \frac{A_{12}}{A_{11}}
        \right) \left( 1 - \frac{P_{cr}}{P_x} \right) \frac{A_{11}}{A_{11}
        + 3 A_{22}} \right)}

    At `P_x = P_{cr}`, `b_{eff} = b/2`, i.e. the whole width is effective;
    the expression is meant for the post-buckling regime, `P_x \ge P_{cr}`.

    References
    ----------
    .. [1] Kassapoglou, C., Design and Analysis of Composite Structures,
       Second Edition, John Wiley & Sons, 2013.

    """
    beff = b/(2*(1 + 2*(1 + A12/A11)*(1 - Pcr/Px)*A11/(A11 + 3*A22)))
    return beff


def calc_Nxx_crit_combined_shear(k, a, b, D11, D12, D22, D66):
    r"""Calculate combined uniaxial-shear buckling load for a composite plate

    This function calculates the critical `N_{xx}` buckling load under the
    current level of shear load given by `N_{xy} = k N_{xx}`. `D16` and `D26`
    are neglected.

    The equation gives two roots, and the smaller of their absolute values
    is returned.

    Parameters
    ----------
    k : float
        Load ratio defined as `k=N_{xy}/N_{xx}`, where `N_{xy}` is the current
        shear load, and `N_{xx}` the current compression load. Dimensionless.
    a, b : float
        Plate length and width, dimension length.
    D11, D12, D22, D66 : float
        Terms of the D matrix, dimension force times length.

    Returns
    -------
    N0 : float
        Critical `N_{xx}` buckling load under the current level of shear load
        given by `N_{xy} = k N_{xx}`, force per unit length, non-negative.

    See Also
    --------
    calc_Nxx_crit_combined_shear_full

    Notes
    -----
    The output of this function is the result of Eq. 6.34, section 6.5 page
    142 of [1]_:

    .. math::

        N_{xx,crit} = \frac{\pi^2}{a^2} \frac{D_{11} + 2 (D_{12} + 2 D_{66})
        a^2/b^2 + D_{22} a^4/b^4}{2 - \frac{8192}{81} \frac{a^2 k^2}{b^2
        \pi^4}} \left( 5 \pm \sqrt{9 + \frac{65536}{81} \frac{a^2
        k^2}{b^2 \pi^4}} \right)

    It is the solution of :func:`calc_Nxx_crit_combined_shear_full` with
    `D_{16} = D_{26} = 0` and `N_{xy} = k N_{xx}`, and with `k = 0` it is
    the load of :func:`calc_Nxx_crit` with `m = n = 1`.

    References
    ----------
    .. [1] Kassapoglou, C., Design and Analysis of Composite Structures,
       Second Edition, John Wiley & Sons, 2013.

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

    This function calculates the critical `N_{xx}` by `N_{xy}`.

    The equation gives two roots, and the smaller of their absolute values
    is returned.

    Parameters
    ----------
    Nxy : float
        Shear load, force per unit length, the `k N_o` of Eqs. 6.32a and 6.33
        of Kassapoglou, i.e. the negative of the `N_{xy}` of the energy, Eq.
        6.32, see Notes. Its sign matters only when `D_{16}` or `D_{26}` is
        not zero.
    a, b : float
        Plate length and width, dimension length.
    D11, D12, D16, D22, D26, D66 : float
        All terms of the D matrix, dimension force times length.

    Returns
    -------
    Nxx_crit : float
        Critical compression load `N_o = -N_x` under the shear load given by
        ``Nxy``, positive in compression, force per unit length,
        non-negative.

    See Also
    --------
    calc_Nxx_crit_combined_shear

    Notes
    -----
    Based on section 6.5 of [1]_, pp. 138-142, for a plate simply supported
    on its four edges (Fig. 6.12). The Rayleigh-Ritz method used therein is
    applied with the two terms of Eq. 6.31, `w = w_1 \sin(\pi x/a) \sin(\pi
    y/b) + w_2 \sin(2 \pi x/a) \sin(2 \pi y/b)`, to the energy of Eq. 6.32
    (p. 140), here including `D_{16}` and `D_{26}`, with the resultants
    `N_x` and `N_{xy}` positive in tension:

    .. math::

        \Pi = \int_0^a \int_0^b \left[ \frac{1}{2} \left( D_{11} w_{,xx}^2
        + 2 D_{12} w_{,xx} w_{,yy} + D_{22} w_{,yy}^2 + 4 D_{66} w_{,xy}^2 + 4
        D_{16} w_{,xx} w_{,xy} + 4 D_{26} w_{,yy} w_{,xy} \right) +
        \frac{1}{2} N_x w_{,x}^2 + N_{xy} w_{,x} w_{,y} \right] dx dy

    As in Eq. 6.32a (p. 141), the compression load is `N_o = -N_x`, positive
    in compression, and the shear load enters as `k N_o`, with `k =
    N_{xy}/N_x`, such that `k N_o = -N_{xy}`. The argument ``Nxy`` of this
    function is `k N_o`, and the returned load is `N_o`. The stationarity
    of `\Pi` with respect to `w_1` and `w_2` gives `([a] - [b]) \{w_1,
    w_2\} = 0`, with

    .. math::

        a_{11} = \frac{\pi^4 D_{11} b}{4 a^3} + \frac{\pi^4 D_{12}}{2 a b}
        + \frac{\pi^4 D_{22} a}{4 b^3} + \frac{\pi^4 D_{66}}{a b}, \quad
        a_{12} = a_{21} = -\frac{160 \pi^2 D_{16}}{9 a^2} - \frac{160
        \pi^2 D_{26}}{9 b^2} \\
        a_{22} = \frac{4 \pi^4 D_{11} b}{a^3} + \frac{8 \pi^4 D_{12}}{a b}
        + \frac{4 \pi^4 D_{22} a}{b^3} + \frac{16 \pi^4 D_{66}}{a b} \\
        b_{11} = \frac{\pi^2 N_o b}{4 a}, \quad
        b_{22} = \frac{\pi^2 N_o b}{a}, \quad
        b_{12} = b_{21} = -\frac{32 k N_o}{9}

    which, with `D_{16} = D_{26} = 0`, is the generalized eigenvalue problem
    of Eqs. 6.33 (p. 141). `\det([a] - [b]) = 0` is a quadratic equation in
    `N_o`, whose roots are

    .. math::

        N_o = \frac{a \left( 36 a_{11} + 9 a_{22} \mp \sqrt{t}
        \right)}{18 \pi^2 b}, \quad
        t = 16384 (k N_o)^2 + 4608 k N_o (a_{12} + a_{21}) + 1296 a_{11}^2
        - 648 a_{11} a_{22} + 1296 a_{12} a_{21} + 81 a_{22}^2

    With `D_{16} = D_{26} = 0`, `a_{12} = a_{21} = 0` and only `(k N_o)^2`
    enters, as in Eq. 6.34 (p. 142); otherwise the sign of the shear load
    changes the result.

    References
    ----------
    .. [1] Kassapoglou, C., Design and Analysis of Composite Structures,
       Second Edition, John Wiley & Sons, 2013.

    """
    a11 = pi**4*D11*b/(4*a**3) + pi**4*D12/(2*a*b) + pi**4*D22*a/(4*b**3) + pi**4*D66/(a*b)
    a12 = -160*pi**2*D16/(9*a**2) - 160*pi**2*D26/(9*b**2)
    a21 = -160*pi**2*D16/(9*a**2) - 160*pi**2*D26/(9*b**2)
    a22 = 4*pi**4*D11*b/a**3 + 8*pi**4*D12/(a*b) + 4*pi**4*D22*a/b**3 + 16*pi**4*D66/(a*b)
    term = 16384*Nxy**2 + 4608*Nxy*a12 + 4608*Nxy*a21 + 1296*a11**2 - 648*a11*a22 + 1296*a12*a21 + 81*a22**2
    Nxx_crit1 = a*(36*a11 + 9*a22 - sqrt(term))/(18*pi**2*b)
    Nxx_crit2 = a*(36*a11 + 9*a22 + sqrt(term))/(18*pi**2*b)
    return min(abs(Nxx_crit1), abs(Nxx_crit2))

