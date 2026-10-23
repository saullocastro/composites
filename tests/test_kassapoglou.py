import sys
sys.path.append('..')

import numpy as np
from numpy import pi
import pytest

from composites import laminated_plate
from composites.kassapoglou import (calc_Nxx_crit,
                                    calc_Nxy_crit,
                                    calc_Nxx_crit_combined_shear,
                                    calc_Nxx_crit_combined_shear_full,
                                    calc_beff,
                                    )

# NOTE Kassapoglou, C., Design and Analysis of Composite Structures, Second
#      Edition, John Wiley & Sons. The reference values are the vertices of
#      the curves of the figures, read from the vector graphics of the PDF of
#      the book, in the units of the figure axes.

# Fig. 6.4, p. 129, the curves of m = 1, 2 and 3 half-waves, (length/width,
# buckling load in N/mm)
FIG_6_4 = {
    1: [(0.6001, 0.7631), (0.6482, 0.7378), (0.6974, 0.7117), (0.751, 0.6864),
        (0.7996, 0.6676), (1.1997, 0.6603), (1.307, 0.6791), (1.4143, 0.6971)],
    2: [(1.4143, 0.6971), (1.4293, 0.6937), (1.4484, 0.6898), (1.4729, 0.6864),
        (1.5025, 0.6825), (1.7998, 0.653), (1.8294, 0.6495), (1.8534, 0.6495),
        (1.8976, 0.6457), (2.073, 0.6457), (2.1507, 0.6495), (2.2289, 0.6495),
        (2.268, 0.653), (2.2971, 0.653), (2.3217, 0.653), (2.3463, 0.6568),
        (2.3804, 0.6568), (2.4144, 0.6603), (2.449, 0.6641)],
    3: [(2.449, 0.6641), (2.4681, 0.6641), (2.5022, 0.6603), (2.5267, 0.6603),
        (2.5608, 0.6568), (2.6335, 0.6568), (2.7218, 0.653), (2.7609, 0.6495),
        (2.8, 0.6495)],
    }

# Fig. 6.11, p. 138, the curve "Present", (a in mm, buckling load in N/mm)
FIG_6_11 = [(101.58, 77.66), (158.78, 35.51), (168.29, 30.06),
            (172.55, 27.72), (177.8, 25.37), (187.31, 21.86), (195.84, 19.11),
            (204.26, 17.58), (212.68, 16.0), (221.21, 15.23), (231.82, 14.42),
            (242.32, 13.66), (254.02, 12.48), (267.8, 11.31), (282.56, 10.14),
            (297.43, 9.37), (314.38, 8.2), (348.17, 6.63), (380.97, 5.45),
            (412.68, 4.69), (444.5, 3.88), (507.93, 3.11)]

# Fig. 6.14, p. 144, the curve "approx", Eq. 6.34, (Nx/Nxcrit, Nxy/Nxycrit)
FIG_6_14 = [(1.0, 0.0176), (0.9984, 0.0352), (0.9935, 0.0703),
            (0.9851, 0.1055), (0.9735, 0.1407), (0.9603, 0.1716),
            (0.9453, 0.2025), (0.9286, 0.2334), (0.9103, 0.2643),
            (0.8921, 0.2909), (0.8838, 0.2994), (0.8805, 0.3037),
            (0.8738, 0.3127), (0.8655, 0.3217), (0.8572, 0.3346),
            (0.8356, 0.3612), (0.824, 0.3745), (0.8107, 0.3921),
            (0.7957, 0.4054), (0.7792, 0.423), (0.7658, 0.4363),
            (0.7475, 0.4539), (0.7359, 0.4672), (0.7243, 0.4758),
            (0.711, 0.4891), (0.696, 0.5024), (0.6828, 0.5109),
            (0.6695, 0.5242), (0.653, 0.5375), (0.6347, 0.5508),
            (0.5963, 0.5817), (0.5566, 0.6122), (0.4054, 0.7314),
            (0.3789, 0.7533), (0.3539, 0.7709), (0.3307, 0.7842),
            (0.3091, 0.8018), (0.289, 0.8151), (0.2707, 0.8284),
            (0.2359, 0.8503), (0.128, 0.9211), (0.1064, 0.9339),
            (0.0847, 0.9472), (0.0682, 0.9562), (0.0533, 0.9648),
            (0.0448, 0.9695), (0.0383, 0.9738), (0.0265, 0.9824)]

# Fig. 7.10, p. 165, (A22/A11, A12/A11): [(Pcr/Px, beff/a), ...]
FIG_7_10 = {
    (10., 0.7): [(0.0996, 0.4522), (0.2005, 0.4571), (0.3001, 0.4604),
                 (0.3998, 0.4653), (0.4502, 0.4686), (0.5007, 0.4703),
                 (0.6003, 0.4752), (0.6999, 0.4802), (0.749, 0.4834),
                 (0.7995, 0.4867), (0.85, 0.49), (0.9004, 0.4916),
                 (1.0, 0.4966)],
    (1., 0.3): [(0.0996, 0.3125), (0.2005, 0.3256), (0.3001, 0.3404),
                (0.3998, 0.3568), (0.4502, 0.3651), (0.5007, 0.3749),
                (0.5511, 0.3831), (0.6003, 0.393), (0.6494, 0.4045),
                (0.6999, 0.416), (0.749, 0.4275), (0.7995, 0.439),
                (0.85, 0.4522), (0.9004, 0.467), (1.0, 0.4966)],
    (0.1, 0.05): [(0.0996, 0.2007), (0.2005, 0.2155), (0.3001, 0.2319),
                  (0.3493, 0.2401), (0.3998, 0.25), (0.4502, 0.2615),
                  (0.5007, 0.273), (0.5511, 0.2862), (0.6003, 0.3009),
                  (0.6494, 0.3158), (0.6999, 0.3338), (0.749, 0.3536),
                  (0.7995, 0.3749), (0.85, 0.3996), (0.8758, 0.4127),
                  (0.9004, 0.4275), (0.925, 0.444), (0.9509, 0.4604),
                  (1.0, 0.4966)],
    }

# D matrix of Fig. 6.4, p. 129, in N m (p. 143)
D11, D12, D22, D66 = 0.66, 0.47, 0.66, 0.49


def _distance_to_curve(x, y, points):
    """Distance of each point (x, y) to the polyline through ``points``"""
    P = np.array(points)
    A, B = P[:-1], P[1:]
    AB = B - A
    dist = []
    for q in np.column_stack((x, y)):
        t = np.clip(((q - A)*AB).sum(axis=1)/(AB*AB).sum(axis=1), 0, 1)
        dist.append(np.min(np.linalg.norm(A + t[:, None]*AB - q, axis=1)))
    return np.array(dist)


def test_calc_Nxx_crit():
    """Kassapoglou's Fig. 6.4, p. 129, Eq. 6.7

    The D matrix is that of the figure, in N m as stated on p. 143. The text
    on p. 129 gives a plate with a fixed length of 508 mm, but the figure is
    reproduced with a fixed width of 254 mm, used here, within 2 %.

    """
    b = 0.254
    for m, points in FIG_6_4.items():
        for AR, Nfig in points:
            a = AR*b
            # N/m to N/mm
            assert np.isclose(calc_Nxx_crit(a, b, m, 1, D11, D12, D22, D66)/1e3,
                              Nfig, rtol=0.02)
            # the envelope is the curve of the figure
            assert np.isclose(calc_Nxx_crit(a, b, None, 1, D11, D12, D22,
                                            D66)/1e3, Nfig, rtol=0.02)
    assert np.isclose(calc_Nxx_crit(b, b, 1, 1, D11, D12, D22, D66),
                      calc_Nxx_crit(b, b, 1, None, D11, D12, D22, D66))
    with pytest.raises(NotImplementedError):
        calc_Nxx_crit(b, b, None, None, D11, D12, D22, D66)


def test_calc_Nxy_crit():
    """Kassapoglou's Fig. 6.11, p. 138, Eqs. 6.28-6.30

    The (0/90)8 laminate of p. 138, eight plies with Ex = Ey = 68.9 GPa, vxy
    = 0.05, Gxy = 4.83 GPa and a ply thickness of 0.1905 mm; sixteen plies
    would give loads eight times those of the figure. With Ex = Ey the ply
    angle does not change the D matrix. The transverse shear moduli, not
    given, are not used.

    The drawn curve of the figure is reproduced within 0.5 % at the vertices
    a = 101.58, 177.80, 242.32, 348.17, 412.68 and 507.93 mm, and within 11 %
    at the other vertices of the drawn curve.

    """
    laminaprop = (68.9e9, 68.9e9, 0.05, 4.83e9, 4.83e9, 4.83e9)
    lam = laminated_plate(stack=[0]*8, plyt=0.1905e-3, laminaprop=laminaprop)
    close = (101.58, 177.8, 242.32, 348.17, 412.68, 507.93)
    for a, Nfig in FIG_6_11:
        # mm to m, N/m to N/mm
        N = calc_Nxy_crit(a/1e3, lam.D11, lam.D12, lam.D16, lam.D22,
                          lam.D66)/1e3
        assert np.isclose(N, Nfig, rtol=0.005 if a in close else 0.11)

    with pytest.raises(RuntimeError):
        calc_Nxy_crit(0.1, -lam.D11, lam.D12, lam.D16, lam.D22, lam.D66)


def test_calc_Nxx_crit_combined_shear_full():
    """Kassapoglou's Fig. 6.14, p. 144, Eqs. 6.34, 6.35 and 6.37

    The curve "approx" is Eq. 6.34 normalized by Eqs. 6.35 and 6.37, and it
    does not depend on the plate dimensions nor on the D matrix, which are
    not given. The D matrix and the plate width of Fig. 6.4 are used. The
    computed points lie within 0.005 of the curve of the figure, whose axes
    run from 0 to 1.

    """
    b = 0.254
    for a in (b, 2*b):
        K = D11 + 2*(D12 + 2*D66)*a**2/b**2 + D22*a**4/b**4
        # Eqs. 6.35 and 6.37
        Nxcrit = pi**2/a**2*K
        Nxycrit = 9*pi**4*b/(32*a**3)*K
        assert np.isclose(Nxcrit, calc_Nxx_crit(a, b, 1, 1, D11, D12, D22,
                                                D66))
        assert np.isclose(Nxcrit, calc_Nxx_crit_combined_shear(0, a, b, D11,
                                                               D12, D22, D66))
        x_values = []
        y_values = []
        # within the range of the figure
        for Nxy in np.linspace(0.02*Nxycrit, 0.98*Nxycrit):
            Nxx = calc_Nxx_crit_combined_shear_full(Nxy, a, b, D11, D12, 0,
                                                    D22, 0, D66)
            assert np.isclose(Nxx, calc_Nxx_crit_combined_shear(Nxy/Nxx, a,
                              b, D11, D12, D22, D66))
            x_values.append(Nxx/Nxcrit)
            y_values.append(Nxy/Nxycrit)
        assert _distance_to_curve(x_values, y_values, FIG_6_14).max() < 0.005


def test_calc_beff():
    """Kassapoglou's Fig. 7.10, p. 165, Eq. 7.15

    The plate width, ``b`` here and `a` in the book, and the forces enter
    only through the ratios shown in the figure. The figure is reproduced
    within 2 %.

    """
    b = 1.
    A11 = 1.
    Pcr = 1.
    for (r22, r12), points in FIG_7_10.items():
        Pcr_Px, beff_b = np.array(points).T
        beff = calc_beff(b, Pcr/Pcr_Px, Pcr, A11, r12*A11, r22*A11)
        assert np.allclose(beff/b, beff_b, rtol=0.02)
    # the whole width is effective at buckling
    assert np.isclose(calc_beff(b, Pcr, Pcr, A11, 0.3*A11, A11), b/2)


if __name__ == '__main__':
    test_calc_Nxx_crit()
    test_calc_Nxy_crit()
    test_calc_Nxx_crit_combined_shear_full()
    test_calc_beff()
