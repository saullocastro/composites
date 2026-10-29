r"""
Transverse shear stiffness of FSDT laminates

The values published by Rohwer (1988) and Vlachoutsis (1992) are collected in
``ROHWER_1988`` and ``VLACHOUTSIS_1992``, with the place in each paper where
they are printed. The cases that neither paper covers are compared with the
independent NumPy implementation of the section "Reference implementation".

Conventions: ``Ats = [[A44, A45], [A45, A55]]``, 4 <-> yz, 5 <-> xz.

Rohwer orders the transverse shear quantities as {xz, yz}, so that his ``H11``,
``H12`` and ``H22`` are ``A55``, ``A45`` and ``A44``, with ``Hbar`` and
``Hbarbar`` the constant-strain (``Abar``) and constant-stress (``Abarbar``)
stiffnesses, and his bending-twisting terms ``D13`` and ``D23`` are ``D16`` and
``D26``. Vlachoutsis' ``K13`` and ``K23`` are ``scf_k13`` and ``scf_k23``.
"""
import sys
sys.path.append('..')
from decimal import Decimal

import numpy as np
import pytest

from composites import laminated_plate, isotropic_plate
from composites.utils import read_laminaprop
from composites.core import Laminate, Lamina


# Published values
# ----------------
#
# Results computed in the papers are kept as strings, exactly as printed, so
# that each comparison uses half a unit of the last printed digit as tolerance,
# see ``assert_printed``. Input data are kept as floats.

# Rohwer, K. "Improved transverse shear stiffnesses for layered finite elements",
# DFVLR-FB 88-32, 1988. Pages as printed in the report.
ROHWER_1988 = {
    # Section 4.1, p. 17: for a homogeneous isotropic plate "the computed
    # relation between H and Hbar or Hbarbar always amounts to 5/6,
    # independent of the Poisson's ratio or the position of the reference
    # surface"
    'isotropic': dict(H_Hbar=5/6, H_Hbarbar=5/6),

    # Table 1, p. 19: "CFRP material constants", kN/mm^2
    'CFRP': dict(E_L=138.0, E_T=9.3, G_LT=4.6, G_TT=2.3, nu_LT=0.3),

    # Section 4.1, p. 18: "A laminate with a stacking sequence of [0,90]_s is
    # built from this material where each of the four layers is 0.25mm thick"
    'cross_ply': dict(stack=[0, 90, 90, 0], plyt=0.25,
        # Table 2, p. 19: "Shear stiffnesses in kN/mm"
        Hbar11='3.450', Hbar22='3.450',
        Hbarbar11='3.067', Hbarbar22='3.067',
        H11='2.313', H22='2.521',
        # p. 18: "Coupling terms are zero."
        H12=0.,
        # p. 18: "the reductions of the improved shear stiffnesses H amount to
        # 0.75 and 0.82 for the x- and y-direction, respectively, as compared
        # with the Hbarbar-values"
        H11_Hbarbar11='0.75', H22_Hbarbar22='0.82'),

    # Section 4.1, pp. 17-18, and Fig. 5, p. 18: symmetric sandwich whose
    # faces and core "shall consist of homogeneous, isotropic material; a
    # Poisson's ratio of nu = 0.3 is assumed". p. 18: "The minimum value
    # H11/Hbarbar11 is a little less than 0.45", which is reached by the
    # lowest curve of Fig. 5, E_f/E_c = 10000. "A little less" is read here as
    # less than 0.01 below.
    'sandwich': dict(nu=0.3, Ef_Ec=10000., min_H11_Hbarbar11=(0.44, 0.45)),

    # Section 4.2, pp. 20-21, and Fig. 8, p. 21: simply supported rectangular
    # plate made of the cross_ply laminate, with "a side length ratio of a/b =
    # 1.5", under a load that "follows a double cosine function with the
    # maximum at the plate center and zero at the edges". Ratios between the
    # Mindlin and the Kirchhoff center deflections, h being the plate thickness.
    'plate': dict(a_b=1.5,
        # p. 21: "At a/h = 10, for instance, the Mindlin displacement exceeds
        # the Kirchhoff one by 25% and 28% for Hbar and Hbarbar, respectively,
        # whereas with H it is 36% higher."
        a_h_10=dict(Hbar='1.25', Hbarbar='1.28', H='1.36'),
        # p. 21: "At a/h = 19, the transverse shear influence reaches 10%, and
        # for a/h = 11 the deviation from Kirchhoff is already more than 30%",
        # which refers to the H curve, the only one above 30% at a/h = 11
        a_h_19=dict(H='1.10'),
        a_h_11=dict(H_more_than=1.30)),
}

# NOTE Rohwer (1988), Sec. 4.3, pp. 22-23, gives a second example, the
#      sandwich solar panel of Fig. 9, p. 22: faces of four CFRP layers of
#      0.125 mm, [0,90,+45,-45] on top and [-45,+45,90,0] at the bottom, the
#      -45 layers adhered to the core, a foam core of 9 mm with E = 0.1
#      kN/mm^2 and nu = 0.3 (Table 3, p. 22), and the CFRP of Table 1. Its
#      stiffnesses, Table 4, p. 23, are not used as a test, because they
#      cannot be reproduced from these data and are not consistent with the
#      stated lay-up:
#      - the transverse shear stiffnesses, H11 = 0.416, H12 = 0.008 and
#        H22 = 0.403 kN/mm, depend on the foam only through its shear
#        modulus; with nu = 0.3, G = 0.0385 kN/mm^2, Eq. (43) gives 0.400,
#        0.008 and 0.387; all three round to the printed values for G
#        between 0.03999 and 0.04005 kN/mm^2, i.e. nu = 0.25 with E = 0.1,
#        which is the suspected value of the computation;
#      - D13 = D23 = 9.36 kNmm, the bending-twisting couplings D16 and D26,
#        are reproduced, 9.3573, but D11 = 1366.08, D12 = 395.09 and D22 =
#        1321.97 kNmm are not, +0.12, +0.62 and +0.47 per cent with nu = 0.3
#        and nearly the same with nu = 0.25;
#      - the printed bending entries contradict each other for this lay-up:
#        for any ply material, the +-45 layers and the isotropic core add the
#        same amount to D11 and D22 and Qbar16 = +-(Q11 - Q22)/4 at +-45
#        degrees, such that (D11 - D22)/D13 = 4 (I0 - I90)/(I+45 - I-45),
#        with I the integral of z^2 over the layers of each angle, a number
#        given by the layer positions alone: 4.2162 for Fig. 9, against
#        4.709 to 4.716 from the printed values within half a unit of their
#        last digits. No change of the material constants, of the foam or of
#        the CFRP, can reproduce Table 4; its bending entries were computed
#        for layer positions, or a lay-up, other than those of Fig. 9, or
#        are misprinted. A search over all orderings of the face layers and
#        over the material constants and thicknesses found no single
#        correction of the data that reproduces the whole table.
#      The only other publication of this panel found is by the same author,
#      Rohwer, K., "Improved transverse shear stiffnesses for anisotropic
#      sandwich plates and shells", Sandwich Constructions 2, Vol. 1, EMAS,
#      1992, pp. 465-484, not available here. The sandwich case is tested
#      instead against Vlachoutsis (1992) and Allen (1969), see
#      test_vlachoutsis_published_sandwich,
#      test_sandwich_vlachoutsis_1992_example_close_to_allen_1969 and
#      test_sandwich_antiplane_core_converges_to_allen_1969.

# Vlachoutsis, S. "Shear correction factors for plates and shells", Int.
# Journal for Numerical Methods in Engineering, Vol. 33, 1537-1552, 1992. Pages
# as printed in the journal.
VLACHOUTSIS_1992 = {
    # Eq. (23), p. 1541: "The homogeneous case is that where K1 = K2 = ... =
    # Ks = K = 5/6"
    'homogeneous': dict(K=5/6),

    # Section 4, "Laminate of n plies-One factor", pp. 1541-1542: "alternate
    # plies of 0deg, 90deg orientations relative to the x-axis of the plate
    # with the two centre layers oriented at 90deg to provide symmetry. All
    # plies have the same thickness (see Figure 3)", with 0deg at the bottom
    # in Figure 3, p. 1542. The factors are given for n = 120.
    'cross_ply': dict(n=120,
        # Eqs. (26a-d), p. 1541, and Eqs. (28a,b), p. 1542
        I=dict(E1_E2=25., G12_E2=0.5, G13_E2=0.5, G23_E2=0.2, nu12=0.25,
               K13='0.6808', K23='0.6794'),
        # Eqs. (27a-d), p. 1541, and Eqs. (30a,b), p. 1542
        II=dict(E1_E2=40., G12_E2=0.6, G13_E2=0.6, G23_E2=0.5, nu12=0.25,
                K13='0.8316', K23='0.8208')),

    # Section 4, "Symmetrical sandwich", p. 1546: "Example. If Dc/Df = Gc/Gf =
    # 10^-3 and p = c/h = 0.95 there is obtained", with p the core-to-total
    # thickness ratio of Eq. (31), rf the share of the transverse shear energy
    # of the facings of Eq. (48), and Kf, Kc the factors of the facings and of
    # the core of Eqs. (46a,b)
    'sandwich_example': dict(p=0.95, Dc_Df=1e-3, Gc_Gf=1e-3,
        rf='1.755e-5', # Eq. (49)
        K='0.01964', # Eq. (50)
        Kf='3.513e-7', # Eq. (51a)
        Kc='1.053'), # Eq. (51b)

    # Section 5, "Sandwich square plate with isotropic facings and core",
    # p. 1548, with "nu1 = nu2 = 0.3" in the caption of Fig. 11. "The index 1
    # denotes the facings and the index 2 the core."
    'sandwich_isotropic': dict(nu=0.3,
        # Table I, p. 1548: "Shear correction factors and energy ratio c/h =
        # 0.8"
        table_I=dict(c_h=0.8, rows=[
            # E1/E2, K, K1, K2, rf
            (1., '0.8333', '0.07133', '1.024', '0.01712'),
            (5., '0.5324', '0.01083', '1.184', '0.01130'),
            (10., '0.3525', '0.003442', '1.225', '0.006974'),
            (15., '0.2625', '0.001667', '1.241', '0.005012'),
            (50., '0.09385', '0.0001701', '1.265', '0.001678'),
            ]),
        # Table II, p. 1548: "Shear correction factors and energy ratio E1/E2 =
        # 15"
        table_II=dict(E1_E2=15., rows=[
            # c/h, K, K1, K2, rf
            (0.9, '0.4087', '0.001164', '1.088', '0.001779'),
            (0.8, '0.2625', '0.001667', '1.241', '0.005012'),
            (0.6, '0.1636', '0.002826', '1.772', '0.01570'),
            (0.4, '0.1400', '0.005677', '3.162', '0.03883'),
            ])),

    # Section 5, "Shallow spherical sandwich shell under local loading",
    # pp. 1550-1551, with "facings (isotropic) and core (antiplane)". The
    # geometry is given in the caption of Fig. 13, p. 1550, in m, with d = c +
    # t by Eq. (68) and t the facing thickness.
    'spherical_shell': dict(R=2., h=0.025, c=0.024, d=0.0245, nu=0.3,
        # caption of Fig. 14, p. 1551; alpha of Eq. (70) fixes the ratio
        # between the Young's modulus of the facings and the shear modulus of
        # the core through Eq. (69)
        alpha=0.1, K='0.07769'),
}


# Chow, T. S. "On the propagation of flexural waves in an orthotropic laminated
# plate and its response to an impulsive load", Journal of Composite
# Materials, Vol. 5, 306-319, 1971. Pages as printed in the journal.
CHOW_1971 = {
    # Section VII, p. 316: "a 4-layer symmetric cross-ply fiber-reinforced
    # composite in the directions (0, 90, 90, 0), and 0 degree being parallel
    # to x-axis of the body; the thickness of each layer is the same", with
    # Eq. (32): "G12/E22 = 0.6, G23/E22 = 0.5, nu12 = 0.25 together with
    # E11/E22 = 3 for glass-epoxy and E11/E22 = 40 for graphite-epoxy".
    # Table 1, p. 316: "Values of k1 and k2 Calculated from the Equations
    # (12)", printed as k1^2 and k2^2, which are scf_k13 and scf_k23. G13 is
    # not given; G13 = G12, a transversely isotropic ply, is the reading
    # under which Eq. (12) reproduces the other three values of Table 1,
    # whereas G13 = G23 gives 0.858, 0.758, 0.868 and 0.479
    'cross_ply': dict(stack=[0, 90, 90, 0], G12_E2=0.6, G23_E2=0.5,
        nu12=0.25, rows=[
            # E11/E22, k1^2, k2^2
            (3., '0.815', '0.807'),
            (40., '0.828', '0.521'),
            ]),
    # Table 1, p. 316, and p. 312: "homogeneous plate isotropic or
    # anisotropic", k1^2 = 5/6
    'homogeneous': dict(K=5/6),
}

# NOTE The value k2^2 = 0.521 of Chow (1971), Table 1, p. 316, graphite-epoxy,
#      is not reproduced by his own Eqs. (11)-(12), p. 311, with the data of
#      Eq. (32): they give 0.52165, which rounds to 0.522, whereas the three
#      other values of the table, 0.815, 0.807 and 0.828, are reproduced to
#      half a unit of their last digit with the same equations and data,
#      0.81459, 0.80700 and 0.82773. The Poisson ratio has no effect, the
#      factor 1/(1 - nu12 nu21) cancelling in Eq. (12), and the printed
#      value is not a truncation, since 0.81459 is printed as 0.815. Within
#      the stated data, 0.521 needs E11/E22 between 40.10 and 40.77, or
#      G13/E22 between 0.5971 and 0.5996 with E11/E22 = 40, i.e. no reading
#      of the printed data. The graphite-epoxy column of Table 2, p. 317, is
#      not consistent with these data either: D11/D22 = 5.98 for E11/E22 =
#      40 against the printed 6.98, and (D12 + 2 D66)/D22 = 0.247 against
#      0.063, whereas the glass-epoxy column, 2.2 and 1.13, is reproduced
#      within 1 %, 2.20 and 1.14. Hence k2^2 of graphite-epoxy is checked
#      against Eqs. (11)-(12), written out in chow_whitney_factor, and not
#      against the printed digits, see test_chow_published_cross_ply.

# Whitney, J. M. "Shear correction factors for orthotropic laminates under
# static load", Journal of Applied Mechanics, Vol. 40, 302-304, 1973. Pages as
# printed in the journal. p. 303: "each ply has the following unidirectional
# properties: EL = 25 x 10^6 psi, ET = 10^6 psi, nu_LT = 0.25, G_LT = 0.5 x
# 10^6 psi, G_TT = 0.2 x 10^6 psi", "Both a 4-layer symmetric and a 2-layer
# unsymmetric square plate are considered", [0/90/90/0] and [0/90] in his
# Fig. 1, "For the two-ply laminate k1^2 = k2^2 = 0.8212, and for the
# four-ply laminate k1^2 = 0.5952 and k2^2 = 0.7205". p. 304, a sandwich with
# faces of "unidirectional graphite/epoxy oriented at 0 deg with respect to
# the x-axis" of thickness h/10, and the core "Exx = Eyy = 0.04 x 10^6 psi,
# Ezz = 0.5 x 10^6 psi, Gxz = Gyz = 0.06 x 10^6 psi, Gxy = 0.016 x 10^6 psi,
# nu_zx = nu_zy = nu_xy = 0.25", "k1^2 = 0.4098 and k2^2 = 0.6915".
WHITNEY_1973 = {
    'ply': (25., 1., 0.25, 0.5, 0.5, 0.2),
    'laminates': [([0, 90], '0.8212', '0.8212'),
                  ([0, 90, 90, 0], '0.5952', '0.7205')],
    'sandwich': dict(core=(0.04, 0.04, 0.25, 0.016, 0.06, 0.06), p=0.8,
                     k13='0.4098', k23='0.6915'),
}

# NOTE The two factors of the sandwich of Whitney (1973), p. 304, are not
#      reproduced by his own Eqs. (5)-(7), p. 303, with the data of p. 304,
#      whereas the three factors of his laminates, p. 303, are reproduced to
#      half a unit of their last digit by the same equations, see
#      test_whitney_published_laminates:
#      - k1^2: Eq. (7) gives 0.40936, 4.4 units of the last printed digit,
#        0.11 %, below 0.4098;
#      - k2^2: Eq. (7) gives 0.67199, 2.8 % below 0.6915.
#      Eq. (7) uses the plane-stress Q11 of each layer, Eq. (3); with Young's
#      moduli instead, or with the three-dimensional stiffness for e_zz = 0,
#      including Ezz of the core, the factors are 0.40936 to 0.40937 and
#      0.6720 to 0.6735, so the deviations do not come from these readings.
#      Within the data, 0.4098 needs a face Gxz of 0.4993 instead of 0.5, or
#      faces of 0.0998 h instead of h/10, and 0.6915 needs a face Gyz of
#      0.187 instead of G_TT = 0.2, a core Gyz of 0.0642 instead of 0.06, or
#      faces of 0.0907 h; no single change of the printed data reproduces
#      both. The sandwich is therefore checked against Eqs. (5)-(7), written
#      out in chow_whitney_factor, and not against the printed digits, see
#      test_whitney_published_sandwich.

# Birman, V. and Bert, C. W. "On the choice of shear correction factor in
# sandwich structures", Journal of Sandwich Structures and Materials, Vol. 4,
# 83-95, 2002. Pages as printed in the journal.
BIRMAN_BERT_2002 = {
    # p. 89, below Eq. (13): "It is easy to show that in the case of a single
    # layer, K = 1."
    'homogeneous': dict(K=1.),
    # p. 90: "Type 2 sandwich is a typical shipbuilding configuration with
    # E-glass/vinyl ester facings and a balsa core. The properties of the
    # facings and core for Type 2 are: Ef = 24.4 GPa, Gf = 2.89 GPa and Ec =
    # 2.41 GPa, Gc = 0.103 GPa. Two geometries considered for this beam were
    # hf = 5 mm, hc = 20 mm and hf = hc = 10 mm." Table 1, p. 90, column
    # "Average shear strain", Eq. (14). The beam stiffness Q11 is read as E,
    # i.e. nu = 0. The Type 1 rows of the same table are not reproduced by
    # Eq. (14) with the printed data, which gives 0.984 and 0.955 instead of
    # 0.986 and 0.960, and they are therefore left out.
    'type_2': dict(Ef=24.4, Gf=2.89, Ec=2.41, Gc=0.103, rows=[
        # hf, hc, K
        (5., 20., '0.124'),
        (10., 10., '0.109'),
        ]),
}


def assert_printed(value, printed, rtol=None):
    r"""Assert that ``value`` agrees with the number ``printed`` in a paper

    The tolerance is half a unit of the last printed digit, i.e. ``value`` must
    round to ``printed``, unless the relative tolerance ``rtol`` is given.
    """
    ref = float(printed)
    if rtol is None:
        tol = 0.5*10.**Decimal(printed).as_tuple().exponent
    else:
        tol = rtol*abs(ref)
    assert abs(value - ref) <= tol*(1 + 1e-12), (value, printed)


_m = ROHWER_1988['CFRP']
# laminaprop = (E11, E22, nu12, G12, G13, G23) of a transversely isotropic ply
CFRP = (_m['E_L'], _m['E_T'], _m['nu_LT'], _m['G_LT'], _m['G_LT'], _m['G_TT'])


# Laminates and helpers
# ---------------------

def random_laminate(seed, offset=None):
    rng = np.random.default_rng(seed)
    nplies = int(rng.integers(2, 13))
    stack = rng.choice([0., 90., 45., -45., 30., -60., 15.], nplies)
    stack[rng.random(nplies) < 0.3] = rng.uniform(-90, 90)
    plyts = rng.uniform(0.1, 0.4, nplies)
    laminaprops = []
    for i in range(nplies):
        E1 = rng.uniform(40., 180.)
        E2 = rng.uniform(5., 12.)
        G12 = rng.uniform(3., 6.)
        G13 = G12
        G23 = rng.uniform(1.5, 4.)
        laminaprops.append((E1, E2, rng.uniform(0.2, 0.35), G12, G13, G23))
    if offset is None:
        offset = rng.uniform(-1., 1.)
    return laminated_plate(stack, plyts=plyts, laminaprops=laminaprops,
            offset=offset)


def sandwich(tf=0.1, tc=0.8, Ef_Ec=1e3, nu=0.3, offset=0., **kwargs):
    r"""Symmetric sandwich with isotropic faces of thickness ``tf`` and core
    of thickness ``tc``"""
    Ef = 70.
    laminaprops = [(Ef, nu), (Ef/Ef_Ec, nu), (Ef, nu)]
    return laminated_plate([0., 0., 0.], plyts=[tf, tc, tf],
            laminaprops=laminaprops, offset=offset, **kwargs)


def vlachoutsis_sandwich_K(p, Dc_Df, Gc_Gf):
    r"""Shear correction factors of a symmetrical sandwich, Vlachoutsis (1992)

    Closed forms of pp. 1543-1545: ``A(p)``, Eq. (32); ``T1``-``T5``, Eqs.
    (37)-(41); the single factor ``K``, Eq. (44); the factor ``K_ant`` of the
    antiplane core, ``Dc/Df = 0``, Eq. (45); and the energy ratio of the
    faces ``r_f``, Eq. (48). ``p = c/h``, Eq. (31), ``D = E/(1 - nu12 nu21)``
    of each material, Eq. (33), and ``G`` its transverse shear modulus, Eq.
    (34). ``K`` is the factor of the whole section, ``A55 = K Abar55``, Eqs.
    (24)-(25), p. 1541.
    """
    A = (1 - p)**3*(3*p**2 + 9*p + 8)/15
    T1 = (1 - p**3) + p**3*Dc_Df
    T2 = (1 - p)/Gc_Gf + p
    T3 = ((1 - p**2)**2 + 8/15*Dc_Df**2*p**4
          + 4/3*Dc_Df*p**2*(1 - p**2))
    T4 = A*Gc_Gf + p*T3
    T5 = A*Gc_Gf + p*(1 - p**2)**2
    K = 4/9*T1**2/(T2*T4)
    K_ant = 4/9*(1 - p**3)**2/(T2*T5)
    r_f = A*Gc_Gf/T4
    return K, K_ant, r_f


def chow_whitney_factor(c, G, plyts):
    r"""Shear correction factor of Whitney (1973), Eqs. (3)-(7), p. 303

    Written out here independently of composites, for orthotropic layers in
    cylindrical bending about one direction: ``c`` the plane-stress reduced
    stiffness ``Q11`` of each layer in that direction, ``G`` its transverse
    shear modulus in the same direction and ``plyts`` the thicknesses, from
    the bottom face. With ``z`` from the mid-surface, ``(A11, B11, D11) =
    int Q11 (1, z, z^2) dz`` and ``D = D11 A11 - B11^2``, Eq. (3), the shear
    stress of Eq. (5) is ``tau = [a + Q11 z (2 B11 - A11 z)] Qx/(2 D)``,
    with the constants ``a`` given by the continuity at the interfaces and
    ``tau = 0`` at the bottom face, and Eq. (7) gives ``k^2 = [A55 int
    tau^2/(G Qx^2) dz]^-1``, ``A55 = int G dz``. For a symmetric laminate
    ``B11 = 0`` and this is Eqs. (11)-(12) of Chow (1971), p. 311, as stated
    by Whitney, p. 303. The 3-point Gauss-Legendre rule per layer is exact
    for ``tau^2``, of degree 4.
    """
    plyts = np.asarray(plyts, dtype=float)
    z = -plyts.sum()/2 + np.concatenate(([0.], np.cumsum(plyts)))
    A11 = sum(ck*(z[k+1] - z[k]) for k, ck in enumerate(c))
    B11 = sum(ck*(z[k+1]**2 - z[k]**2)/2 for k, ck in enumerate(c))
    D11 = sum(ck*(z[k+1]**3 - z[k]**3)/3 for k, ck in enumerate(c))
    D = D11*A11 - B11**2
    xg, wg = np.polynomial.legendre.leggauss(3)
    a = -c[0]*z[0]*(2*B11 - A11*z[0])
    integral = 0.
    for k, (ck, Gk) in enumerate(zip(c, G)):
        if k > 0:
            a += (c[k-1] - ck)*z[k]*(2*B11 - A11*z[k])
        zg = (z[k] + z[k+1])/2 + xg*(z[k+1] - z[k])/2
        tau = (a + ck*zg*(2*B11 - A11*zg))/(2*D)
        integral += np.sum(wg*tau**2/Gk)*(z[k+1] - z[k])/2
    A55 = sum(Gk*(z[k+1] - z[k]) for k, Gk in enumerate(G))
    return 1/(A55*integral)


def _orthotropic_factors(stack, plyts, laminaprops):
    r"""``(k13, k23)`` of :func:`chow_whitney_factor` for plies at 0 or 90
    degrees, ``laminaprop = (E1, E2, nu12, G12, G13, G23)``"""
    cx, cy, Gx, Gy = [], [], [], []
    for theta, (E1, E2, nu12, G12, G13, G23) in zip(stack, laminaprops):
        delta = 1 - nu12**2*E2/E1
        Q11, Q22 = E1/delta, E2/delta
        if theta == 0:
            cx.append(Q11); cy.append(Q22); Gx.append(G13); Gy.append(G23)
        else:
            assert theta == 90
            cx.append(Q22); cy.append(Q11); Gx.append(G23); Gy.append(G13)
    return (chow_whitney_factor(cx, Gx, plyts),
            chow_whitney_factor(cy, Gy, plyts))


def _interfaces(lam):
    return lam.offset - lam.h/2 + np.concatenate(([0.],
        np.cumsum([ply.h for ply in lam.plies])))


def _gauss_integral(lam, func):
    zi = _interfaces(lam)
    xi, w = np.polynomial.legendre.leggauss(3)
    out = 0
    for k, ply in enumerate(lam.plies):
        zm, dz = (zi[k] + zi[k+1])/2, (zi[k+1] - zi[k])/2
        for x, wi in zip(xi, w):
            out = out + wi*dz*func(k, ply, zm + dz*x)
    return out


# Exactness
# ---------

@pytest.mark.parametrize('mode', ['rohwer', 'vlachoutsis', 'constant'])
@pytest.mark.parametrize('nu', [0.0, 0.3, 0.49])
@pytest.mark.parametrize('offset', [0., +0.37, -1.5])
def test_isotropic_five_sixths(mode, nu, offset):
    ref = ROHWER_1988['isotropic']
    K = VLACHOUTSIS_1992['homogeneous']['K']
    h = 2.
    E = 70e3
    G = E/(2*(1 + nu))
    lam = isotropic_plate(thickness=h, E=E, nu=nu, offset=offset,
            shear_correction=mode)
    assert np.allclose(lam.Abar_ts, G*h*np.eye(2), rtol=1e-14)
    assert np.allclose(lam.Abarbar_ts, G*h*np.eye(2), rtol=1e-14)
    for Ats, Ats_ref in [(lam.A55, lam.Abar55), (lam.A44, lam.Abar44)]:
        assert np.isclose(Ats/Ats_ref, ref['H_Hbar'], rtol=1e-12, atol=0)
    for Ats, Ats_ref in [(lam.A55, lam.Abarbar55), (lam.A44, lam.Abarbar44)]:
        assert np.isclose(Ats/Ats_ref, ref['H_Hbarbar'], rtol=1e-12, atol=0)
    assert abs(lam.A45) <= 1e-12*G*h
    assert np.isclose(lam.scf_k13, K, rtol=1e-12, atol=0)
    assert np.isclose(lam.scf_k23, K, rtol=1e-12, atol=0)


# Published values, Rohwer (1988)
# -------------------------------

def test_rohwer_published_cross_ply():
    ref = ROHWER_1988['cross_ply']
    lam = laminated_plate(ref['stack'], plyt=ref['plyt'], laminaprop=CFRP)
    assert_printed(lam.Abar55, ref['Hbar11'])
    assert_printed(lam.Abar44, ref['Hbar22'])
    assert_printed(lam.Abarbar55, ref['Hbarbar11'])
    assert_printed(lam.Abarbar44, ref['Hbarbar22'])
    # NOTE closed forms of Eqs. (22) and (25): in the yz plane G23 for the 0
    #      plies and G13 for the 90 plies, the reverse in the xz plane
    G13, G23 = CFRP[4], CFRP[5]
    t = ref['plyt']
    h = 4*t
    for Abar in (lam.Abar44, lam.Abar55):
        assert np.isclose(Abar, 2*t*G23 + 2*t*G13, rtol=1e-14, atol=0)
    for Abarbar in (lam.Abarbar44, lam.Abarbar55):
        assert np.isclose(Abarbar, h**2/(2*t/G23 + 2*t/G13), rtol=1e-14, atol=0)
    # NOTE improved stiffness, Eq. (43), against the independent
    #      implementation RohwerReference
    Ats = RohwerReference(stack_from_laminate(lam)).Ats
    assert np.allclose(lam.Ats, Ats, rtol=1e-12, atol=1e-14*Ats.max())
    assert_printed(lam.A55, ref['H11'])
    # NOTE A44 = 2.52525, and the independent RohwerReference gives the same,
    #      which Rohwer prints as H22 = 2.521, 0.17 % lower; the other three
    #      values of Table 2 round correctly, and both 2.52525/3.067 and
    #      2.521/3.067 round to the reduction 0.82 printed on p. 18, checked
    #      below to half a unit of its last digit. The relative tolerance only
    #      records this deviation of the printed value; the anchors of A44
    #      are Eq. (43), written out in RohwerReference, and the ratio 0.82
    assert_printed(lam.A44, ref['H22'], rtol=1.7e-3)
    assert abs(lam.A45 - ref['H12']) < 1e-12
    assert_printed(lam.A55/lam.Abarbar55, ref['H11_Hbarbar11'])
    assert_printed(lam.A44/lam.Abarbar44, ref['H22_Hbarbar22'])
    assert np.isclose(lam.scf_k13, lam.A55/lam.Abar55, rtol=1e-14)
    assert np.isclose(lam.scf_k23, lam.A44/lam.Abar44, rtol=1e-14)
    assert np.allclose(lam.Ats, [[lam.A44, lam.A45], [lam.A45, lam.A55]])


def test_rohwer_published_sandwich_minimum():
    r"""Fig. 5 of Rohwer (1988), p. 18, against Vlachoutsis (1992), Eq. (44)

    For isotropic faces and core with one Poisson ratio, the cylindrical
    bending of Rohwer's Eq. (34) gives the bending stress of a beam, as
    assumed by Vlachoutsis, such that every point of the curves of Fig. 5,
    ``H11/Hbarbar11 = K Abar55/Abarbar55``, follows from the closed form of
    Eq. (44), with ``Dc/Df = Gc/Gf = Ec/Ef``, and from the closed forms of
    Rohwer's Eqs. (22) and (25); they are compared to round-off: 1e-13, plus
    ``3 eps/(1 - p)`` for the closed form, which forms ``1 - p`` from ``p =
    c/h``, with the machine epsilon ``eps`` as absolute error, and raises it
    to the third power in ``A(p)``, Eq. (32). Measured: at most 8e-15 for
    ``1 - p > 1e-3`` and 1.5e-12 at ``1 - p = 6e-6``. The minimum
    of the lowest curve, ``Ef/Ec = 10000``, is 0.44580 at ``log(a_c/a_f) =
    -1.33``, measured with a bounded scalar search of the closed forms, and
    the printed "a little less than 0.45" is read as less than 0.01 below.
    """
    ref = ROHWER_1988['sandwich']
    ratios = []
    nu = ref['nu']
    Ef = 70.
    Ec = Ef/ref['Ef_Ec']
    Gf, Gc = Ef/(2*(1 + nu)), Ec/(2*(1 + nu))
    # abscissa log(a_c/a_f) of Fig. 5
    for log_ac_af in np.linspace(-6., 6., 241):
        tf, tc = 1., 10**log_ac_af
        lam = sandwich(tf=tf, tc=tc, Ef_Ec=ref['Ef_Ec'], nu=nu)
        h = 2*tf + tc
        p = tc/h
        K = vlachoutsis_sandwich_K(p, Ec/Ef, Gc/Gf)[0]
        Abar = 2*tf*Gf + tc*Gc
        Abarbar = h**2/(2*tf/Gf + tc/Gc)
        assert np.isclose(lam.Abar55, Abar, rtol=1e-13, atol=0)
        assert np.isclose(lam.Abarbar55, Abarbar, rtol=1e-13, atol=0)
        rtol = 1e-13 + 3*np.finfo(float).eps/(1 - p)
        assert np.isclose(lam.A55/lam.Abarbar55, K*Abar/Abarbar, rtol=rtol, atol=0)
        ratios.append(lam.A55/lam.Abarbar55)
    low, high = ref['min_H11_Hbarbar11']
    assert low < min(ratios) < high


def _center_deflection_ratio(lam, Ats, a_h, a_b):
    r"""Mindlin over Kirchhoff center deflection of a simply supported plate

    Navier solution for a specially orthotropic laminate, ``D16 = D26 = 0``
    and ``A45 = 0``, under the load `q_0 \sin(\pi x/a) \sin(\pi y/b)`, which
    is Rohwer's double cosine load with the origin moved to a corner.
    """
    a = a_h*lam.h
    b = a/a_b
    al, be = np.pi/a, np.pi/b
    D11, D12, D22, D66 = lam.D11, lam.D12, lam.D22, lam.D66
    A44, A55 = Ats[0, 0], Ats[1, 1]
    # unknowns: amplitudes of w, phi_x and phi_y
    K = np.array([
        [A55*al**2 + A44*be**2, A55*al, A44*be],
        [A55*al, D11*al**2 + D66*be**2 + A55, (D12 + D66)*al*be],
        [A44*be, (D12 + D66)*al*be, D66*al**2 + D22*be**2 + A44]])
    w_mindlin = np.linalg.solve(K, [1., 0., 0.])[0]
    w_kirchhoff = 1/(D11*al**4 + 2*(D12 + 2*D66)*al**2*be**2 + D22*be**4)
    return w_mindlin/w_kirchhoff


def test_rohwer_published_plate_deflection():
    ref = ROHWER_1988['plate']
    cross_ply = ROHWER_1988['cross_ply']
    lam = laminated_plate(cross_ply['stack'], plyt=cross_ply['plyt'],
            laminaprop=CFRP)
    def ratio(Ats, a_h):
        return _center_deflection_ratio(lam, Ats, a_h, ref['a_b'])
    assert_printed(ratio(lam.Abar_ts, 10), ref['a_h_10']['Hbar'])
    assert_printed(ratio(lam.Abarbar_ts, 10), ref['a_h_10']['Hbarbar'])
    assert_printed(ratio(lam.Ats, 10), ref['a_h_10']['H'])
    assert_printed(ratio(lam.Ats, 19), ref['a_h_19']['H'])
    assert ratio(lam.Ats, 11) > ref['a_h_11']['H_more_than']


# NOTE the sandwich panel of Rohwer (1988), Sec. 4.3, Table 4, p. 23, is not
#      used as a test, see the NOTE after ROHWER_1988


# Published values, Vlachoutsis (1992)
# ------------------------------------

@pytest.mark.parametrize('data_set', ['I', 'II'])
def test_vlachoutsis_published_120_plies(data_set):
    ref = VLACHOUTSIS_1992['cross_ply']
    m = ref[data_set]
    E2 = 1.
    laminaprop = (m['E1_E2']*E2, E2, m['nu12'], m['G12_E2']*E2,
            m['G13_E2']*E2, m['G23_E2']*E2)
    n = ref['n']
    stack = [0, 90]*(n//4) + [90, 0]*(n//4)
    lam = laminated_plate(stack, plyt=1/n, laminaprop=laminaprop,
            shear_correction='vlachoutsis')
    assert_printed(lam.scf_k13, m['K13'])
    assert_printed(lam.scf_k23, m['K23'])
    assert np.isclose(lam.A55, lam.scf_k13*lam.Abar55, rtol=1e-14)
    assert np.isclose(lam.A44, lam.scf_k23*lam.Abar44, rtol=1e-14)


def _vlachoutsis_ply_factors(lam):
    r"""Per-ply factors and energy shares from the recovered stresses

    For ``Qx = 1`` the share of the transverse shear energy stored in ply `k`
    is `U_k/U = I_k/I` and, since `R^2/I = A_{55}`, Vlachoutsis' Eq. (21)
    becomes `K_k = A_{55} (U_k/U)/(G_k h_k)`.
    """
    N = len(lam.plies)
    def energy(k, ply, z):
        tau_xz = lam.calc_transverse_shear_stress(z, 0., 1.)[1]
        return np.eye(N)[k]*tau_xz**2/ply.q55L
    share = _gauss_integral(lam, energy)
    share = share/share.sum()
    Gh = np.array([ply.q55L*ply.h for ply in lam.plies])
    return lam.A55*share/Gh, share


def _vlachoutsis_sandwiches():
    ref = VLACHOUTSIS_1992['sandwich_isotropic']
    params = []
    table = ref['table_I']
    for E1_E2, K, K1, K2, rf in table['rows']:
        params.append(pytest.param(E1_E2, table['c_h'], K, K1, K2, rf,
            id='Table I, E1/E2=%g' % E1_E2))
    table = ref['table_II']
    for c_h, K, K1, K2, rf in table['rows']:
        params.append(pytest.param(table['E1_E2'], c_h, K, K1, K2, rf,
            id='Table II, c/h=%g' % c_h))
    # faces and core of isotropic materials with a common Poisson's ratio
    # have Dc/Df = Gc/Gf = Ec/Ef
    ex = VLACHOUTSIS_1992['sandwich_example']
    params.append(pytest.param(1/ex['Gc_Gf'], ex['p'], ex['K'], ex['Kf'],
        ex['Kc'], ex['rf'], id='Eqs. (49)-(51)'))
    return params


@pytest.mark.parametrize('E1_E2, c_h, K, K1, K2, rf', _vlachoutsis_sandwiches())
def test_vlachoutsis_published_sandwich(E1_E2, c_h, K, K1, K2, rf):
    r"""Isotropic sandwich of Vlachoutsis (1992), Tables I and II, p. 1548,
    and Eqs. (49)-(51), p. 1546

    For faces and core of isotropic materials with one Poisson ratio,
    ``Dc/Df = Gc/Gf = E2/E1``, the factor K is the closed form of Eq. (44)
    and r_f that of Eq. (48), compared to round-off (rtol 1e-12, 1e-10 for
    r_f, a ratio of energies integrated by Gauss points). The cylindrical
    bending of Rohwer (1988), Eq. (34), gives then the bending stress of a
    beam, as assumed by Vlachoutsis, Sec. 2, such that ``'rohwer'`` gives
    the same K. Measured: K of both modes within 2e-15 of Eq. (44) and r_f
    within 3e-15 of Eq. (48) for all rows. K, r_f and the per-ply factors
    K1 and K2 must also round to the printed digits.
    """
    nu = VLACHOUTSIS_1992['sandwich_isotropic']['nu']
    K_eq44, _, rf_eq48 = vlachoutsis_sandwich_K(c_h, 1/E1_E2, 1/E1_E2)
    for mode in ['rohwer', 'vlachoutsis']:
        lam = sandwich(tf=(1 - c_h)/2, tc=c_h, Ef_Ec=E1_E2, nu=nu,
                shear_correction=mode)
        assert np.isclose(lam.scf_k13, K_eq44, rtol=1e-12, atol=0)
        assert np.isclose(lam.scf_k23, K_eq44, rtol=1e-12, atol=0)
    assert_printed(lam.scf_k13, K)
    assert_printed(lam.scf_k23, K)
    kappa, share = _vlachoutsis_ply_factors(lam)
    assert np.isclose(share[0] + share[2], rf_eq48, rtol=1e-10, atol=0)
    assert_printed(kappa[0], K1)
    assert_printed(kappa[2], K1)
    assert_printed(kappa[1], K2)
    assert_printed(share[0] + share[2], rf)


@pytest.mark.parametrize('mode', ['vlachoutsis', 'rohwer'])
def test_vlachoutsis_published_spherical_shell(mode):
    r"""Spherical sandwich shell of Vlachoutsis (1992), pp. 1550-1551

    The facings are isotropic and the core antiplane, so K is the closed form
    ``K_ant`` of Eq. (45), p. 1545, with ``Gc/Gf`` from Eqs. (69) and (70)
    and ``Gf = Ef/(2 (1 + nu))``, here 0.0776843. The core is antiplane,
    ``laminaprop = (0, 0, 0, 0, Gc, Gc)``, with zero moduli in the plane of
    the facings (Allen, 1969, p. 14), for which Eq. (44) of Vlachoutsis
    reduces to Eq. (45); measured, both modes agree with Eq. (45) within
    3e-15, checked with ``rtol = 1e-12`` for round-off. The printed ``K =
    0.07769`` of the caption of Fig. 14,
    p. 1551, which is reported from the finite element computations of his
    Reference 8, is 0.57 units of its last digit above Eq. (45), just outside
    half a unit: it would need ``alpha = 0.099992`` instead of 0.1, so the
    rounding of the printed data does not explain it. The printed value is
    therefore checked with one unit of its last digit, the anchor of the
    test being Eq. (45) of the same paper.
    """
    ref = VLACHOUTSIS_1992['spherical_shell']
    R, h, c, d, nu, alpha = (ref[key] for key in
            ['R', 'h', 'c', 'd', 'nu', 'alpha'])
    t = d - c
    assert np.isclose(c + 2*t, h)
    # Eq. (70) solved for e_s, then Eq. (69) solved for Ef/Gc
    e_s = np.sqrt(3*alpha**2/(1 + alpha**2))
    Ef_Gc = e_s*h*R*np.sqrt(12*(1 - nu**2))/(3*t*d)
    Ef = 70e3
    Gc = Ef/Ef_Gc
    # antiplane core, Dc/Df = 0
    core = (0., 0., 0., 0., Gc, Gc)
    lam = laminated_plate([0, 0, 0], plyts=[t, c, t],
            laminaprops=[(Ef, nu), core, (Ef, nu)], shear_correction=mode)
    Gf = Ef/(2*(1 + nu))
    K, K_ant, _ = vlachoutsis_sandwich_K(c/h, 0., Gc/Gf)
    assert K == K_ant
    for k in [lam.scf_k13, lam.scf_k23]:
        assert np.isclose(k, K_ant, rtol=1e-12, atol=0)
        # NOTE one unit of the last printed digit, see the docstring
        assert abs(k - float(ref['K'])) <= 1e-5


# Allen (1969)
# ------------
#
# Allen, H. G., "Analysis and Design of Structural Sandwich Panels", Pergamon
# Press, 1969, Chap. 2, pp. 8-47,
# https://doi.org/10.1016/B978-0-08-012870-2.50006-7. Pages as printed in
# the book.

def test_sandwich_vlachoutsis_1992_example_close_to_allen_1969():
    r"""Example of Vlachoutsis (1992), p. 1546, and Allen (1969), Eq. (2.15b)

    For ``Dc/Df = Gc/Gf = 1e-3`` and ``p = c/h = 0.95`` Vlachoutsis prints K
    = 0.01964, Eq. (50), and K = 0.01865, Eqs. (52)-(53), for a constant core
    shear stress Q/c ignoring the faces, "5 per cent below". The shear
    stiffness of Allen (1969), Eq. (2.15b), p. 17, ``A G = G_c d^2/c`` per
    unit width, with ``d = (h + c)/2``, Eq. (2.3), p. 9, the distance between
    the centroids of the faces, includes the lever arm d of the faces and
    rounds to the same printed K = 0.01964; composites gives A55/(G_c d^2/c)
    - 1 = 9.6e-5. The deviation from Allen is bounded with the terms that
    Allen neglects in :func:`test_sandwich_antiplane_core_converges_to_allen_1969`;
    here the core is not antiplane, so only the printed digits are compared.
    """
    nu = 0.3
    h = 1.
    p = 0.95
    t = (1 - p)*h/2
    c = p*h
    d = (h + c)/2
    Ec, Ef = 1., 1000.
    Gc, Gf = Ec/(2*(1 + nu)), Ef/(2*(1 + nu))
    lam = laminated_plate([0, 0, 0], plyts=[t, c, t],
                          laminaprops=[(Ef, nu), (Ec, nu), (Ef, nu)])
    Abar = 2*t*Gf + c*Gc
    assert np.isclose(lam.Abar55, Abar, rtol=1e-14, atol=0)
    assert_printed(lam.scf_k13, '0.01964')
    assert_printed(Gc*d**2/c/Abar, '0.01964')
    assert_printed(Gc*c/Abar, '0.01865')
    assert 9e-5 < lam.A55/(Gc*d**2/c) - 1 < 1e-4


def test_sandwich_antiplane_core_converges_to_allen_1969():
    r"""Antiplane core: closed form and convergence to Allen (1969)

    An antiplane core is "an idealized core in which the modulus of
    elasticity in planes parallel with the faces is zero but the shear
    modulus in planes perpendicular to the faces is finite", ``E_c = 0``, p.
    14, given here as ``laminaprop = (0, 0, nu, 0, G_c, G_c)``, see
    :func:`.read_laminaprop`. Then the shear stress is constant in the core,
    Eq.
    (2.12), p. 13, ``tau_c = Q E_f t d/(2 D)``, with ``D = E_f (t^3/6 + t
    d^2/2)``, Eq. (2.2), p. 9, i.e. ``tau_c = Q/(d phi)``, ``phi = 1 +
    t^2/(3 d^2)``, and parabolic in the faces, Eq. (2.9), p. 11, ``tau_f = Q
    E_f (h^2/4 - z^2)/(2 D)``. The complementary energy, Rohwer (1988), Eq.
    (43), gives::

        A = (G_c d^2/c) phi^2/(1 + eps),
        eps = (G_c/G_f) 2 t (h^2/3 - h t/2 + t^2/5)/(c d^2)

    with eps the ratio of the face to the core shear energy. This is Eq. (45)
    of Vlachoutsis (1992), p. 1545, ``A = K_ant sum_i G_i h_i``, verified
    here to round-off. Allen's shear stiffness ``G_c d^2/c``, Eq. (2.15b),
    p. 17, neglects the local bending stiffness of the faces, the first term
    of Eq. (2.2), "less than 1%" of the second when ``3 (d/t)^2 > 100``, Eq.
    (2.4), p. 10, i.e. when ``phi - 1 < 0.01``, and the shear strains of the
    faces, footnote of p. 17, i.e. eps. Hence the derived bounds, which tend
    to zero as ``d/t`` and ``G_f/G_c`` grow::

        -eps <= A/(G_c d^2/c) - 1 <= phi^2 - 1

    The family below has ``t/c = 10^-k/2`` and ``G_c/G_f = 10^-k/10``, i.e. a
    constant share of the core in ``Abar = sum_i G_i h_i``. Measured, with
    ``E_f = 70 GPa``, ``nu = 0.3``, ``c = 10 mm``:

    ===  =======  =========  =============  ===========  ======
     k    eps      phi^2-1    A/A_Allen - 1  composites   K
                                             vs Eq. (45)
    ===  =======  =========  =============  ===========  ======
     1    3.4e-4   1.51e-3    1.17e-3        7e-16        0.1003
     2    3.3e-6   1.65e-5    1.32e-5        3e-14        0.0918
     3    3.3e-8   1.67e-7    1.33e-7        2e-13        0.0910
    ===  =======  =========  =============  ===========  ======

    The deviation from Allen decreases by a factor 100 per step, as both
    bounds, while K tends to 1/11, such that the constant 5/6 overestimates
    the shear stiffness 9.2 times. The column "composites vs Eq. (45)" is
    the largest relative difference of ``A44`` and ``A55`` of ``'rohwer'``
    and ``'vlachoutsis'`` and of the independent ``RohwerReference`` from
    ``K_ant sum_i G_i h_i``. With ``E_c = 0``, Eq. (44) of Vlachoutsis
    reduces to Eq. (45), ``K = K_ant``, checked exactly, and the only
    allowance added to the derived bounds is 1e-12 for round-off, also the
    tolerance against Eq. (45).
    """
    Ef = 70e3
    nu = 0.3
    Gf = Ef/(2*(1 + nu))
    c = 10.
    lf = (Ef, Ef, nu, Gf, Gf, Gf)
    devs = []
    for k in [1, 2, 3]:
        t = c*10.**(-k)/2
        Gc = Gf*10.**(-k)/10
        # NOTE antiplane core, E_c = 0 in the plane of the faces
        lc = (0., 0., nu, 0., Gc, Gc)
        plyts = [t, c, t]
        h = c + 2*t
        p = c/h
        d = (h + c)/2
        phi = 1 + t**2/(3*d**2)
        eps = Gc/Gf*2*t*(h**2/3 - h*t/2 + t**2/5)/(c*d**2)
        A_allen = Gc*d**2/c
        Abar = 2*t*Gf + c*Gc
        K, K_ant, _ = vlachoutsis_sandwich_K(p, 0., Gc/Gf)
        # NOTE the closed form above is Eq. (45) of Vlachoutsis (1992), and
        #      Eq. (44) reduces to it for Dc/Df = 0
        assert np.isclose(A_allen*phi**2/(1 + eps), K_ant*Abar, rtol=1e-13, atol=0)
        assert K == K_ant
        slack = 1e-12
        for mode in ['rohwer', 'vlachoutsis']:
            lam = laminated_plate([0, 0, 0], plyts=plyts,
                                  laminaprops=[lf, lc, lf],
                                  shear_correction=mode)
            assert np.isclose(lam.A55, K*Abar, rtol=1e-12, atol=0)
            assert np.isclose(lam.A44, K*Abar, rtol=1e-12, atol=0)
        Ats = RohwerReference(stack_from_laminate(lam)).Ats
        assert np.allclose(Ats, K*Abar*np.eye(2), rtol=1e-12,
                           atol=1e-12*K*Abar)
        dev = lam.A55/A_allen - 1
        assert -eps - slack <= dev <= phi**2 - 1 + slack
        devs.append(dev)
        # NOTE Allen's condition (2.4) holds
        assert 3*(d/t)**2 > 100
        assert (5/6)/K > 8
    # NOTE convergence: the deviation from Allen decreases with k, as the
    #      bounds, by a factor 100 per step
    assert devs[0] > 50*devs[1] > 0 and devs[1] > 50*devs[2] > 0
    assert devs[-1] < 2e-7


# Properties
# ----------

SEEDS = list(range(12))


@pytest.mark.parametrize('seed', SEEDS)
def test_offset_invariance(seed):
    rng = np.random.default_rng(1000 + seed)
    nplies = int(rng.integers(2, 9))
    stack = rng.uniform(-90, 90, nplies)
    Ats = []
    for offset in [0., +0.4, -2.0]:
        lam = laminated_plate(stack, plyt=0.25, laminaprop=CFRP,
                offset=offset)
        Ats.append(lam.Ats)
    for other in Ats[1:]:
        assert np.allclose(other, Ats[0], rtol=1e-9,
                atol=1e-9*np.abs(Ats[0]).max())


@pytest.mark.parametrize('seed', SEEDS)
def test_upper_bound(seed):
    lam = random_laminate(seed)
    eig = np.linalg.eigvalsh(lam.Abar_ts - lam.Ats)
    assert np.all(eig > 0)
    assert lam.scf_k13 <= 1.
    assert lam.scf_k23 <= 1.


def test_upper_bound_sandwich():
    lam = sandwich()
    assert np.all(np.linalg.eigvalsh(lam.Abar_ts - lam.Ats) > 0)


@pytest.mark.parametrize('seed', SEEDS)
def test_stacking_order(seed):
    rng = np.random.default_rng(2000 + seed)
    nplies = int(rng.integers(1, 7))
    half = list(rng.uniform(-90, 90, nplies))
    # symmetric laminate
    stack = half + half[::-1]
    lam = laminated_plate(stack, plyt=0.2, laminaprop=CFRP)
    rev = laminated_plate(stack[::-1], plyt=0.2, laminaprop=CFRP)
    assert np.allclose(rev.Ats, lam.Ats, rtol=1e-12, atol=1e-12)
    # unsymmetric laminate turned upside down
    unsym = laminated_plate(half, plyt=0.2, laminaprop=CFRP)
    unsym_rev = laminated_plate(half[::-1], plyt=0.2, laminaprop=CFRP)
    assert np.allclose(unsym_rev.Ats, unsym.Ats, rtol=1e-10, atol=1e-12)
    # rotating every ply by 90 deg swaps 13 and 23
    rot = laminated_plate([t + 90 for t in stack], plyt=0.2, laminaprop=CFRP)
    assert np.isclose(rot.scf_k13, lam.scf_k23, rtol=1e-10)
    assert np.isclose(rot.scf_k23, lam.scf_k13, rtol=1e-10)
    assert np.isclose(rot.A44, lam.A55, rtol=1e-10)
    assert np.isclose(rot.A55, lam.A44, rtol=1e-10)
    assert np.isclose(rot.A45, -lam.A45, rtol=1e-10, atol=1e-12)


# Reference implementation
# ------------------------
#
# Independent NumPy implementation of the equilibrium approach of Rohwer
# (1988), used where no published value exists. Equation numbers are Rohwer's,
# but the ordering is the one of this module, tau = {tau_yz, tau_xz} and Q =
# {Qy, Qx}, instead of Rohwer's {xz, yz}.

class LayerStack:
    r"""Layer-wise description of a laminate

    Parameters
    ----------
    z : (N+1,) array
        Layer interface coordinates measured from the reference surface.
    C : (N, 3, 3) array
        In-plane constitutive matrix of each layer, in the laminate axes.
    Cs : (N, 2, 2) array
        Transverse shear constitutive matrix of each layer, ordered as
        ``[[C44, C45], [C45, C55]]``.
    """
    def __init__(self, z, C, Cs):
        self.z = np.asarray(z, dtype=float)
        self.C = np.asarray(C, dtype=float)
        self.Cs = np.asarray(Cs, dtype=float)
        self.N = self.C.shape[0]
        assert self.z.shape == (self.N + 1,)
        assert self.Cs.shape == (self.N, 2, 2)
        self.h_k = np.diff(self.z)
        self.h = self.z[-1] - self.z[0]

    @property
    def ABD(self):
        r"""``[[A, B], [B, D]]`` of Rohwer's Eq. (17)"""
        z = self.z
        A = np.einsum('k,kij->ij', z[1:] - z[:-1], self.C)
        B = np.einsum('k,kij->ij', (z[1:]**2 - z[:-1]**2)/2, self.C)
        D = np.einsum('k,kij->ij', (z[1:]**3 - z[:-1]**3)/3, self.C)
        return np.block([[A, B], [B, D]])


def stack_from_laminate(lam):
    r"""Build a :class:`LayerStack` from a :class:`.Laminate`"""
    z = [-lam.h/2 + lam.offset]
    C, Cs = [], []
    for ply in lam.plies:
        z.append(z[-1] + ply.h)
        C.append([[ply.q11L, ply.q12L, ply.q16L],
                  [ply.q12L, ply.q22L, ply.q26L],
                  [ply.q16L, ply.q26L, ply.q66L]])
        Cs.append([[ply.q44L, ply.q45L],
                   [ply.q45L, ply.q55L]])
    return LayerStack(z, C, Cs)


def Ats_constant_strain(st):
    r"""Rohwer's Eq. (22), `\bar{A}_{ts} = \sum_k C_s^{(k)} h_k`"""
    return np.einsum('k,kij->ij', st.h_k, st.Cs)


def Ats_constant_stress(st):
    r"""Rohwer's Eq. (25), `\bar{\bar{A}}_{ts} = h^2 [\sum_k (C_s^{(k)})^{-1}
    h_k]^{-1}`"""
    S = np.einsum('k,kij->ij', st.h_k, np.linalg.inv(st.Cs))
    return st.h**2*np.linalg.inv(S)


class RohwerReference:
    r"""Equilibrium approach of Rohwer (1988)

    Attributes
    ----------
    Ats : (2, 2) array
        ``[[A44, A45], [A45, A55]]`` of Rohwer's Eq. (43).
    """
    def __init__(self, st):
        self.st = st
        # Eq. (19), the lower-left block is B*^T
        Hstar = np.linalg.inv(st.ABD)
        self._Astar, self._Bstar = Hstar[:3, :3], Hstar[:3, 3:]
        self._BstarT, self._Dstar = Hstar[3:, :3], Hstar[3:, 3:]
        # Eq. (34), columns ordered as {Qy, Qx}
        self._Lx = np.zeros((6, 2))
        self._Lx[3, 1] = 1. # M_xx,x = Q_x
        self._Lx[5, 0] = 1. # M_xy,x = Q_y
        self._Ly = np.zeros((6, 2))
        self._Ly[4, 0] = 1. # M_yy,y = Q_y
        self._Ly[5, 1] = 1. # M_xy,y = Q_x
        # Eqs. (39) and (40), a^(k) = sum_i (c^(i) - c^(i-1)) Pi(z_i)
        self._a_x = np.zeros((st.N, 6))
        self._a_y = np.zeros((st.N, 6))
        acc_x, acc_y = np.zeros(6), np.zeros(6)
        c1_prev, c2_prev = np.zeros(3), np.zeros(3)
        for k in range(st.N):
            Pi_k = self._Pi(st.z[k])
            acc_x = acc_x + (st.C[k, 0] - c1_prev) @ Pi_k
            acc_y = acc_y + (st.C[k, 1] - c2_prev) @ Pi_k
            self._a_x[k], self._a_y[k] = acc_x, acc_y
            c1_prev, c2_prev = st.C[k, 0], st.C[k, 1]
        # Eq. (43), the integrand is a polynomial of fourth order in z, which
        # the 3-point Gauss-Legendre rule integrates exactly
        xi, w = np.polynomial.legendre.leggauss(3)
        Cs_inv = np.linalg.inv(st.Cs)
        compliance = np.zeros((2, 2))
        for k in range(st.N):
            zm, dz = (st.z[k] + st.z[k+1])/2, (st.z[k+1] - st.z[k])/2
            for x, wi in zip(xi, w):
                f = self.f(zm + dz*x, k)
                compliance += wi*dz*f.T @ Cs_inv[k] @ f
        Ats = np.linalg.inv(compliance)
        self.Ats = (Ats + Ats.T)/2

    def _Pi(self, z):
        r"""``[z A* + z^2/2 B*^T, z B* + z^2/2 D*]`` of Eqs. (37) and (38)"""
        return np.hstack([z*self._Astar + z**2/2*self._BstarT,
                          z*self._Bstar + z**2/2*self._Dstar])

    def f(self, z, k=None):
        r"""Distribution matrix of Eq. (41), ``tau = f(z) Q``"""
        st = self.st
        if k is None:
            k = min(max(np.searchsorted(st.z, z, side='right') - 1, 0),
                    st.N - 1)
        Pi = self._Pi(z)
        row_y = (self._a_y[k] - st.C[k, 1] @ Pi) @ self._Ly
        row_x = (self._a_x[k] - st.C[k, 0] @ Pi) @ self._Lx
        return np.vstack([row_y, row_x])

    def tau(self, z, Q):
        r"""``{tau_yz, tau_xz}`` at ``z`` for ``Q = {Qy, Qx}``"""
        return self.f(z) @ np.asarray(Q, dtype=float)


# Regression against the reference implementation
# ------------------------------------------------

@pytest.mark.parametrize('stack, sign', [
    ([45, -45, -45, 45], -1),
    ([-45, 45, 45, -45], +1),
    ])
def test_rohwer_angle_ply_coupling(stack, sign):
    r"""The ``[0,90]s`` laminate of Rohwer (1988), Table 2, rotated by 45
    degrees

    Not published as numbers: Fig. 7, p. 20, plots the stiffnesses of this
    laminate against the rotation angle, with ``H11 = H22``, about 0.85
    ``Hbarbar11``, and a coupling ``-H12/Hbarbar11`` of about 0.04 at 45
    degrees, read from the figure; only the bound ``|A45| > 0.03
    Abarbar55``, that the coupling is not negligible, is taken from it.
    The values are checked against the independent implementation of Eq.
    (43), RohwerReference, and with properties: ``Abar45 = 0``, the 0 and 90
    plies of the shear moduli cancelling; mirroring the laminate, ``y -> -y``,
    changes the sign of A45 only; and rotating every ply by 90 degrees maps
    one laminate onto the other, swapping A44 and A55 and changing the sign
    of A45, hence ``A44 = A55``. Measured: A45 = -+0.12042, scf_k13 = scf_k23
    = 0.75029, ``A55/Abarbar55 = 0.844`` and ``|A45|/Abarbar55 = 0.0393``.
    """
    lam = laminated_plate(stack, plyt=0.25, laminaprop=CFRP)
    mirror = laminated_plate([-t for t in stack], plyt=0.25, laminaprop=CFRP)
    rotated = laminated_plate([t + 90 for t in stack], plyt=0.25,
                              laminaprop=CFRP)
    ref = RohwerReference(stack_from_laminate(lam)).Ats
    assert np.allclose(lam.Ats, ref, rtol=1e-12, atol=1e-12*ref.max())
    # the ad-hoc (k13 + k23)/2*Abar45 cannot represent this coupling
    assert abs(lam.Abar45) < 1e-14*lam.Abar55
    assert sign*lam.A45 > 0.03*lam.Abarbar55
    assert np.isclose(mirror.A45, -lam.A45, rtol=1e-12, atol=0)
    assert np.allclose([mirror.A44, mirror.A55], [lam.A44, lam.A55],
                       rtol=1e-12, atol=0)
    assert np.allclose([rotated.A44, rotated.A55, rotated.A45],
                       [lam.A55, lam.A44, -lam.A45], rtol=1e-12,
                       atol=0)
    assert np.isclose(lam.A44, lam.A55, rtol=1e-12, atol=0)
    assert np.isclose(lam.scf_k13, lam.scf_k23, rtol=1e-12, atol=0)
    assert lam.scf_k13 < 1.


def _tensor_rotation(Ats, thetadeg):
    r"""Tensor rotation of ``Ats = [[A44, A45], [A45, A55]]``, Rohwer
    (1988), Eq. (44), p. 20, to the frame in which every ply angle is
    increased by ``thetadeg``"""
    c = np.cos(np.deg2rad(thetadeg))
    s = np.sin(np.deg2rad(thetadeg))
    T = np.array([[c, s], [-s, c]])
    return T @ Ats @ T.T


def test_rohwer_rotation_dependence_cross_ply():
    r"""The result of ``'rohwer'`` is not invariant to a rotation of the axes

    Rohwer (1988), Fig. 7, p. 20, plots the stiffnesses of the ``[0,90]s``
    laminate of Table 2 against the rotation angle, computed with Eq. (43) in
    the rotated axes, solid lines, and by the tensor rotation of the result
    in the axes of the plies, Eq. (44), dashed lines: "Deviations from the
    dashed lines are evident, especially around a rotation angle of 45 deg".
    The statement of :meth:`.Laminate.calc_transverse_shear_stiffness` is
    checked here, with both curves from the independent RohwerReference: at
    45 degrees ``A44 = A55`` are 7.0 % and ``A45`` is 13.5 % above the
    tensor rotation, and A44 or A55 alone deviate by up to 7.2 %, at 41 and
    49 degrees (1-degree steps; at 0.5-degree steps also at 41 degrees),
    all to half a unit of the last digit given. Read from Fig. 7 at 45
    degrees, the solid line is at about 0.85 and the dashed one at about
    0.79 of ``Hbarbar11``, consistent with 0.844 and 0.789.
    """
    stack = [0, 90, 90, 0]
    lam0 = laminated_plate(stack, plyt=0.25, laminaprop=CFRP)
    H0 = RohwerReference(stack_from_laminate(lam0)).Ats
    assert np.allclose(lam0.Ats, H0, rtol=1e-12, atol=1e-14*H0.max())
    devs = []
    for theta in range(0, 91):
        lam = laminated_plate([t + theta for t in stack], plyt=0.25,
                              laminaprop=CFRP)
        H = RohwerReference(stack_from_laminate(lam)).Ats
        assert np.allclose(lam.Ats, H, rtol=1e-12, atol=1e-14*H.max())
        Hr = _tensor_rotation(H0, theta)
        # NOTE the constant-strain and constant-stress stiffnesses follow
        #      the tensor rotation, which fixes its sign convention; the
        #      frame dependence is that of the corrected stiffness only
        for name in ['Abar_ts', 'Abarbar_ts']:
            ref = _tensor_rotation(getattr(lam0, name), theta)
            assert np.allclose(getattr(lam, name), ref, rtol=1e-12,
                               atol=1e-12*ref.max())
        devs.append(max(H[0, 0]/Hr[0, 0], H[1, 1]/Hr[1, 1]) - 1)
        if theta == 45:
            assert np.isclose(H[0, 0], H[1, 1], rtol=1e-12, atol=0)
            assert_printed(100*(H[0, 0]/Hr[0, 0] - 1), '7.0')
            assert_printed(100*(H[0, 1]/Hr[0, 1] - 1), '13.5')
            assert_printed(H[1, 1]/lam0.Abarbar55, '0.844')
            assert_printed(Hr[1, 1]/lam0.Abarbar55, '0.789')
        if theta == 0:
            assert np.allclose(H, Hr, rtol=1e-14, atol=1e-14*H.max())
    assert_printed(100*max(devs), '7.2')
    assert int(np.argmax(devs)) in (41, 49)


def _grid():
    cases = []
    for offset in [0., 0.3, -1.2]:
        for name, stack in [
                ('unidirectional', [0, 0, 0, 0]),
                ('cross-ply', [0, 90, 0, 90]),
                ('angle-ply', [30, -30, 30, -30]),
                ('quasi-isotropic', [0, 45, -45, 90, 90, -45, 45, 0]),
                ('unsymmetric', [0, 45, 90, -30, 15]),
                ]:
            cases.append((name, offset,
                laminated_plate(stack, plyt=0.125, laminaprop=CFRP,
                    offset=offset)))
        cases.append(('sandwich', offset, sandwich(offset=offset)))
        cases.append(('sandwich angle faces', offset,
            laminated_plate([30, 0, -60], plyts=[0.1, 0.8, 0.1],
                laminaprops=[CFRP, (0.14, 0.14, 0.3, 0.05, 0.05, 0.02), CFRP],
                offset=offset)))
    for seed in SEEDS:
        cases.append(('random', seed, random_laminate(seed)))
    return cases


@pytest.mark.parametrize('name, param, lam', _grid())
def test_regression_reference(name, param, lam):
    st = stack_from_laminate(lam)
    ref = RohwerReference(st).Ats
    atol = 1e-10*np.abs(ref).max()
    assert np.allclose(lam.Ats, ref, rtol=1e-10, atol=atol)
    assert np.allclose(lam.Abar_ts, Ats_constant_strain(st), rtol=1e-12)
    assert np.allclose(lam.Abarbar_ts, Ats_constant_stress(st), rtol=1e-10,
            atol=atol)


def test_sandwich_soft_core():
    r"""Isotropic sandwich with ``c/h = 0.8`` and ``Ef/Ec = 1000``

    The closed form of Vlachoutsis (1992), Eq. (44), with ``p = 0.8`` and
    ``Dc/Df = Gc/Gf = 1e-3``, gives K = 0.0050821, for ``'rohwer'`` as well,
    see :func:`test_vlachoutsis_published_sandwich`; measured within 7e-16.
    The previous implementation returned 1.2564.
    """
    K = vlachoutsis_sandwich_K(0.8, 1e-3, 1e-3)[0]
    for mode in ['rohwer', 'vlachoutsis']:
        lam = sandwich(shear_correction=mode)
        assert np.isclose(lam.scf_k13, K, rtol=1e-12, atol=0)
        assert np.isclose(lam.scf_k23, K, rtol=1e-12, atol=0)


# Other modes and errors
# ----------------------

def test_modes():
    lam = laminated_plate([45, -45, 0, 90], plyt=0.25, laminaprop=CFRP)
    Abar = lam.Abar_ts
    lam.shear_correction = None
    lam.calc_transverse_shear_stiffness()
    assert np.allclose(lam.Ats, Abar)
    assert lam.scf_k13 == 1. and lam.scf_k23 == 1.
    lam.shear_correction = 'constant'
    lam.calc_transverse_shear_stiffness()
    assert np.allclose(lam.Ats, 5/6*Abar)
    lam.shear_correction = 'vlachoutsis'
    lam.calc_transverse_shear_stiffness()
    k13, k23 = lam.scf_k13, lam.scf_k23
    assert np.isclose(lam.A45, (k13 + k23)/2*Abar[0, 1])
    lam.shear_correction = 'rohwer'
    lam.calc_transverse_shear_stiffness()
    assert np.allclose(lam.Ats, RohwerReference(stack_from_laminate(lam)).Ats)
    lam.shear_correction = 'wrong'
    with pytest.raises(ValueError, match='shear_correction'):
        lam.calc_transverse_shear_stiffness()


def test_shear_correction_argument():
    stack = [45, -45, 0, 90]
    ref = laminated_plate(stack, plyt=0.25, laminaprop=CFRP)
    assert ref.shear_correction == 'rohwer'
    for mode in ['rohwer', 'vlachoutsis', 'constant', None]:
        lam = laminated_plate(stack, plyt=0.25, laminaprop=CFRP,
                shear_correction=mode)
        assert lam.shear_correction == mode
        ref.shear_correction = mode
        ref.calc_transverse_shear_stiffness()
        assert np.allclose(lam.Ats, ref.Ats, rtol=1e-14)
    plate = isotropic_plate(thickness=2., E=70., nu=0.3,
            shear_correction=None)
    assert np.allclose(plate.Ats, plate.Abar_ts)
    with pytest.raises(ValueError, match='shear_correction'):
        laminated_plate(stack, plyt=0.25, laminaprop=CFRP,
                shear_correction='rohwer1988')


@pytest.mark.parametrize('calc_scf, mode', [(True, 'rohwer'), (False, None)])
def test_deprecated_calc_scf_argument(calc_scf, mode):
    with pytest.warns(DeprecationWarning, match='calc_scf'):
        lam = laminated_plate([0, 90, 90, 0], plyt=0.25, laminaprop=CFRP,
                calc_scf=calc_scf, shear_correction='constant')
    assert lam.shear_correction == mode
    with pytest.warns(DeprecationWarning, match='calc_scf'):
        plate = isotropic_plate(thickness=2., E=70., nu=0.3,
                calc_scf=calc_scf)
    assert plate.shear_correction == mode


def test_deprecated_calc_scf_method():
    ref = ROHWER_1988['cross_ply']
    lam = laminated_plate(ref['stack'], plyt=ref['plyt'], laminaprop=CFRP,
            shear_correction=None)
    with pytest.warns(DeprecationWarning, match='calc_scf'):
        k13, k23 = lam.calc_scf()
    assert lam.shear_correction == 'rohwer'
    assert (k13, k23) == (lam.scf_k13, lam.scf_k23)
    assert_printed(lam.A55, ref['H11'])
    assert np.isclose(k13, lam.A55/lam.Abar55, rtol=1e-14)
    assert np.isclose(k23, lam.A44/lam.Abar44, rtol=1e-14)


def test_deprecated_Atrans():
    lam = laminated_plate([45, -45, -45, 45], plyt=0.25, laminaprop=CFRP)
    with pytest.warns(DeprecationWarning, match='Ats'):
        Atrans = lam.Atrans
    assert np.array_equal(Atrans, lam.Ats)


def _manual_laminate(stack, plyt, laminaprops, shear_correction):
    lam = Laminate()
    lam.shear_correction = shear_correction
    for thetadeg, laminaprop in zip(stack, laminaprops):
        ply = Lamina()
        ply.thetadeg = thetadeg
        ply.h = plyt
        ply.matlamina = read_laminaprop(laminaprop)
        ply.rebuild()
        lam.plies.append(ply)
    lam.calc_constitutive_matrix()
    return lam


def test_singular_Cs():
    laminaprops = [CFRP, (138., 9.3, 0.3, 4.6, 0., 2.3), CFRP]
    with pytest.raises(ValueError, match='Ply 1'):
        laminated_plate([0, 30, 0], plyt=0.25, laminaprops=laminaprops)
    # Vlachoutsis only needs q55L > 0 and q44L > 0
    with pytest.raises(ValueError, match='Ply 1'):
        _manual_laminate([0, 0, 0], 0.25, laminaprops, 'vlachoutsis')
    lam = _manual_laminate([0, 30, 0], 0.25, laminaprops, 'vlachoutsis')
    assert 0 < lam.scf_k13 < 1 and 0 < lam.scf_k23 < 1
    # modes that do not need inv(Cs) do not raise
    lam = _manual_laminate([0, 30, 0], 0.25, laminaprops, None)
    assert np.allclose(lam.Ats, lam.Abar_ts)
    assert np.isnan(lam.Abarbar55)
    lam = _manual_laminate([0, 30, 0], 0.25, laminaprops, 'constant')
    assert np.allclose(lam.Ats, 5/6*lam.Abar_ts)


def test_singular_ABD():
    lam = _manual_laminate([0, 0], 0.25, [CFRP, CFRP], None)
    for ply in lam.plies:
        ply.h = 0.
    lam.shear_correction = 'rohwer'
    with pytest.raises(ValueError):
        lam.calc_constitutive_matrix()


def test_empty_laminate():
    lam = Laminate()
    lam.calc_constitutive_matrix()
    assert lam.A44 == 0. and lam.A45 == 0. and lam.A55 == 0.


# Transverse shear stress recovery
# --------------------------------

def _f(lam, z):
    r"""Distribution matrix f(z), {tau_yz, tau_xz} = f(z) {Qy, Qx}"""
    return np.column_stack([lam.calc_transverse_shear_stress(z, 1., 0.),
                            lam.calc_transverse_shear_stress(z, 0., 1.)])


@pytest.mark.parametrize('offset', [0., +0.37, -1.5])
def test_isotropic_parabola(offset):
    h = 3.
    lam = isotropic_plate(thickness=h, E=70e3, nu=0.3, offset=offset)
    Qx, Qy = 7.3, -2.1
    tol = 1e-13*3/(2*h) # machine precision relative to the peak stress
    for z in np.linspace(-h/2, h/2, 21) + offset:
        zbar = z - offset
        tau_yz, tau_xz = lam.calc_transverse_shear_stress(z, Qy, Qx)
        assert abs(tau_xz - 3*Qx/(2*h)*(1 - 4*zbar**2/h**2)) <= tol*abs(Qx)
        assert abs(tau_yz - 3*Qy/(2*h)*(1 - 4*zbar**2/h**2)) <= tol*abs(Qy)
        # no Qy contribution in tau_xz, and no Qx contribution in tau_yz
        assert abs(lam.calc_transverse_shear_stress(z, 1., 0.)[1]) <= tol
        assert abs(lam.calc_transverse_shear_stress(z, 0., 1.)[0]) <= tol


@pytest.mark.parametrize('seed', SEEDS)
def test_free_surfaces_and_continuity(seed):
    lam = random_laminate(seed)
    zi = _interfaces(lam)
    scale = 1/lam.h
    assert np.allclose(_f(lam, zi[0]), 0, atol=1e-10*scale)
    assert np.allclose(_f(lam, zi[-1]), 0, atol=1e-10*scale)
    eps = 1e-13*lam.h
    for zk in zi[1:-1]:
        below = _f(lam, zk - eps)
        above = _f(lam, zk + eps)
        assert np.allclose(below, above, atol=1e-10*scale)


@pytest.mark.parametrize('seed', SEEDS)
def test_in_plane_resultants(seed):
    lam = random_laminate(seed)
    # partition of Hstar = inv(ABD) with P(z) = dPi/dz = [A* + z B*^T, B* + z D*]
    Hstar = np.linalg.inv(lam.ABD)
    def CP(k, ply, z):
        C = np.array([[ply.q11L, ply.q12L, ply.q16L],
                      [ply.q12L, ply.q22L, ply.q26L],
                      [ply.q16L, ply.q26L, ply.q66L]])
        P = np.hstack([Hstar[:3, :3] + z*Hstar[3:, :3],
                       Hstar[:3, 3:] + z*Hstar[3:, 3:]])
        return C @ P
    assert np.allclose(_gauss_integral(lam, CP), np.hstack([np.eye(3),
        np.zeros((3, 3))]), atol=1e-10)
    # the distribution must integrate to the shear forces, which exercises
    # the integration constants of the implementation
    assert np.allclose(_gauss_integral(lam, lambda k, ply, z: _f(lam, z)),
            np.eye(2), atol=1e-10)


@pytest.mark.parametrize('name, param, lam', _grid())
def test_stress_regression_reference(name, param, lam):
    ref = RohwerReference(stack_from_laminate(lam))
    zi = _interfaces(lam)
    rng = np.random.default_rng(3)
    Q = np.array([0.7, -1.3])
    for z in np.concatenate((zi, rng.uniform(zi[0], zi[-1], 10))):
        tau = lam.calc_transverse_shear_stress(z, Q[0], Q[1])
        tau_ref = ref.tau(z, Q)
        assert np.allclose(tau, tau_ref, rtol=1e-10,
                atol=1e-10*np.abs(Q).max()/lam.h)


def test_stress_modes_and_errors():
    stack = [30, -45, 90]
    ref = laminated_plate(stack, plyt=0.25, laminaprop=CFRP)
    for mode in ['vlachoutsis', 'constant', None]:
        lam = laminated_plate(stack, plyt=0.25, laminaprop=CFRP,
                shear_correction=mode)
        assert lam.calc_transverse_shear_stress(0.1, 1., 2.) == \
            ref.calc_transverse_shear_stress(0.1, 1., 2.)
    with pytest.raises(ValueError, match='outside'):
        ref.calc_transverse_shear_stress(0.4, 1., 1.)
    with pytest.raises(ValueError, match='outside'):
        ref.calc_transverse_shear_stress(-0.38, 1., 1.)


# Chow (1971), Birman and Bert (2002) and thickness-shear frequency
# -----------------------------------------------------------------

@pytest.mark.parametrize('E1_E2, k13, k23', CHOW_1971['cross_ply']['rows'])
def test_chow_published_cross_ply(E1_E2, k13, k23):
    ref = CHOW_1971['cross_ply']
    E2 = 1.
    laminaprop = (E1_E2*E2, E2, ref['nu12'], ref['G12_E2']*E2,
                  ref['G12_E2']*E2, ref['G23_E2']*E2)
    lam = laminated_plate(ref['stack'], plyt=0.25, laminaprop=laminaprop,
            shear_correction='chow')
    # NOTE Eqs. (11)-(12) of Chow (1971), p. 311, written out in
    #      chow_whitney_factor, to round-off
    k1, k2 = _orthotropic_factors(ref['stack'], [0.25]*4, [laminaprop]*4)
    assert np.isclose(lam.scf_k13, k1, rtol=1e-12, atol=0)
    assert np.isclose(lam.scf_k23, k2, rtol=1e-12, atol=0)
    assert_printed(lam.scf_k13, k13)
    if E1_E2 == 40.:
        # NOTE k2^2 = 0.52165 by Eqs. (11)-(12), printed as 0.521, which is
        #      not consistent with the data of the paper, see the NOTE after
        #      CHOW_1971; only the closed form above is checked
        assert abs(lam.scf_k23 - float(k23)) > 0.5e-3
    else:
        assert_printed(lam.scf_k23, k23)


@pytest.mark.parametrize('seed', SEEDS)
@pytest.mark.parametrize('offset', [0., 0.8])
def test_chow_equals_vlachoutsis_symmetric(seed, offset):
    # for a symmetric laminate the neutral surfaces of Vlachoutsis coincide
    # with the mid-surface used by Chow
    rng = np.random.default_rng(2000 + seed)
    half = list(rng.uniform(-90, 90, int(rng.integers(1, 5))))
    plyts = list(rng.uniform(0.1, 0.4, len(half)))
    stack = half + half[::-1]
    plyts = plyts + plyts[::-1]
    lams = [laminated_plate(stack, plyts=plyts, laminaprop=CFRP,
                offset=offset, shear_correction=mode)
            for mode in ['chow', 'vlachoutsis']]
    assert np.allclose(lams[0].Ats, lams[1].Ats, rtol=1e-12, atol=0)


@pytest.mark.parametrize('stack', [[0, 90], [0, 90, 45], [30, -30]])
def test_chow_unsymmetric(stack):
    with pytest.raises(ValueError, match='symmetric'):
        laminated_plate(stack, plyt=0.25, laminaprop=CFRP,
                shear_correction='chow')


@pytest.mark.parametrize('mode, K', [
    ('chow', CHOW_1971['homogeneous']['K']),
    ('birman_bert', BIRMAN_BERT_2002['homogeneous']['K']),
    ('thickness_shear', np.pi**2/12)])
@pytest.mark.parametrize('nu', [0.0, 0.3, 0.49])
@pytest.mark.parametrize('offset', [0., +0.37, -1.5])
def test_homogeneous_limits(mode, K, nu, offset):
    # pi^2/12 is the value of Mindlin (1951), independent of nu
    lam = isotropic_plate(thickness=2., E=70e3, nu=nu, rho=2.7e-6,
            offset=offset, shear_correction=mode)
    assert np.isclose(lam.scf_k13, K, rtol=1e-12, atol=0)
    assert np.isclose(lam.scf_k23, K, rtol=1e-12, atol=0)
    assert np.allclose(lam.Ats, K*lam.Abar_ts, rtol=1e-12)


def _birman_bert_eq14(Qf, Gf, Qc, Gc, hf, hc):
    r"""Eq. (14) of Birman and Bert (2002), symmetric sandwich beam"""
    x = hc/(2*hf + hc)
    num = 2/3*(Qf*(1 - x**3) + Qc*x**3)
    den = ((Gf*(1 - x) + Gc*x)*(Qf*(2/3 - x + x**3/3)/Gf
           + Qf*x*(1 - x**2)/Gc + 2/3*Qc*x**3/Gc))
    return num/den


def _beam_sandwich(Ef, Gf, Ec, Gc, hf, hc, **kwargs):
    faces = (Ef, Ef, 0., Gf, Gf, Gf)
    core = (Ec, Ec, 0., Gc, Gc, Gc)
    return laminated_plate([0, 0, 0], plyts=[hf, hc, hf],
            laminaprops=[faces, core, faces], **kwargs)


@pytest.mark.parametrize('hf, hc, K', BIRMAN_BERT_2002['type_2']['rows'])
def test_birman_bert_published_type_2(hf, hc, K):
    ref = BIRMAN_BERT_2002['type_2']
    lam = _beam_sandwich(ref['Ef'], ref['Gf'], ref['Ec'], ref['Gc'], hf, hc,
            shear_correction='birman_bert')
    assert_printed(lam.scf_k13, K)
    assert_printed(lam.scf_k23, K)


@pytest.mark.parametrize('Ef_Ec, Gf_Gc, hf_hc', [(1e1, 1e1, 0.25),
    (1e3, 30., 0.05), (5e2, 1e3, 1.), (1., 1., 0.5)])
def test_birman_bert_closed_form(Ef_Ec, Gf_Gc, hf_hc):
    Ef, Gf, hc = 70., 26., 1.
    lam = _beam_sandwich(Ef, Gf, Ef/Ef_Ec, Gf/Gf_Gc, hf_hc*hc, hc,
            shear_correction='birman_bert')
    K = _birman_bert_eq14(Ef, Gf, Ef/Ef_Ec, Gf/Gf_Gc, hf_hc*hc, hc)
    assert np.isclose(lam.scf_k13, K, rtol=1e-12)


def _thickness_shear_fe(lam, alpha, nelem=4000):
    r"""Lowest non-zero thickness-shear frequency with linear finite
    elements through the thickness, as an independent reference"""
    import scipy.linalg
    zi = _interfaces(lam)
    nodes = [zi[0]]
    G = []
    rho = []
    for k, ply in enumerate(lam.plies):
        n = max(4, int(nelem*ply.h/lam.h))
        nodes += list(np.linspace(zi[k], zi[k+1], n + 1)[1:])
        G += [ply.q55L if alpha == 0 else ply.q44L]*n
        rho += [ply.matlamina.rho]*n
    nodes = np.array(nodes)
    K = np.zeros((nodes.size, nodes.size))
    M = np.zeros((nodes.size, nodes.size))
    for e in range(nodes.size - 1):
        le = nodes[e+1] - nodes[e]
        K[e:e+2, e:e+2] += G[e]/le*np.array([[1, -1], [-1, 1]])
        M[e:e+2, e:e+2] += rho[e]*le/6*np.array([[2, 1], [1, 2]])
    return np.sqrt(np.sort(scipy.linalg.eigh(K, M, eigvals_only=True))[1])


@pytest.mark.parametrize('kwargs', [
    dict(stack=[0, 0, 0], plyts=[0.5, 8., 0.5],
         laminaprops=[(70e3, 0.3), (70., 0.3), (70e3, 0.3)],
         rhos=[2.7e-6, 1e-7, 2.7e-6], offset=0.7),
    dict(stack=[0, 90, 0, 90, 90, 0, 90, 0], plyt=0.125, laminaprop=CFRP,
         rho=1.6e-6),
    dict(stack=[0, 90, 45], plyts=[0.2, 0.5, 0.3], laminaprop=CFRP,
         rhos=[1.6e-6, 1.0e-6, 2.0e-6]),
    ])
def test_thickness_shear_finite_elements(kwargs):
    pytest.importorskip('scipy')
    lam = laminated_plate(shear_correction='thickness_shear', **kwargs)
    I2c = lam.intrhoz2 - lam.intrhoz**2/lam.intrho
    for alpha, kappa, Abar in [(0, lam.scf_k13, lam.Abar55),
                               (1, lam.scf_k23, lam.Abar44)]:
        omega = _thickness_shear_fe(lam, alpha)
        # FSDT: omega^2 = kappa*Abar/I2c, and O(1/nelem^2) for the FE
        assert np.isclose(kappa, omega**2*I2c/Abar, rtol=1e-6)


@pytest.mark.parametrize('mode', ['birman_bert', 'thickness_shear'])
@pytest.mark.parametrize('seed', SEEDS[:4])
def test_new_modes_offset_invariance(mode, seed):
    rng = np.random.default_rng(3000 + seed)
    nplies = int(rng.integers(2, 9))
    stack = rng.uniform(-90, 90, nplies)
    rhos = rng.uniform(1e-6, 3e-6, nplies)
    Ats = [laminated_plate(stack, plyt=0.25, laminaprop=CFRP, rhos=rhos,
               offset=offset, shear_correction=mode).Ats
           for offset in [0., +0.4, -2.0]]
    for other in Ats[1:]:
        assert np.allclose(other, Ats[0], rtol=1e-9,
                atol=1e-9*np.abs(Ats[0]).max())


def test_thickness_shear_density():
    with pytest.raises(ValueError, match='densit'):
        laminated_plate([0, 90], plyt=0.25, laminaprop=CFRP,
                shear_correction='thickness_shear')
    with pytest.raises(ValueError, match='densit'):
        laminated_plate([0, 90], plyt=0.25, laminaprop=CFRP,
                rhos=[1.6e-6, 0.], shear_correction='thickness_shear')


# Whitney (1973)
# --------------

@pytest.mark.parametrize('stack, k13, k23', WHITNEY_1973['laminates'])
@pytest.mark.parametrize('mode', ['whitney', 'vlachoutsis'])
def test_whitney_published_laminates(stack, k13, k23, mode):
    r"""Laminates of Whitney (1973), p. 303, to half a unit of the last
    printed digit, and against Eqs. (3)-(7), written out in
    chow_whitney_factor, to round-off; ``[0, 90]`` is unsymmetric, ``B11 !=
    0``. G13 = G_LT and G23 = G_TT of the transversely isotropic ply."""
    plyts = [0.25]*len(stack)
    lam = laminated_plate(stack, plyts=plyts, laminaprop=WHITNEY_1973['ply'],
            shear_correction=mode)
    k1, k2 = _orthotropic_factors(stack, plyts,
                                  [WHITNEY_1973['ply']]*len(stack))
    assert np.isclose(lam.scf_k13, k1, rtol=1e-12, atol=0)
    assert np.isclose(lam.scf_k23, k2, rtol=1e-12, atol=0)
    assert_printed(lam.scf_k13, k13)
    assert_printed(lam.scf_k23, k23)


def test_whitney_published_sandwich():
    r"""Sandwich of Whitney (1973), p. 304, against his Eqs. (3)-(7)

    The faces of thickness h/10 are the 0-degree ply of p. 303 and the core
    is given by its in-plane plane-stress data, ``Exx = Eyy``, ``nu_xy`` and
    ``Gxy``, and ``Gxz = Gyz``; Ezz does not enter Eq. (3). Eqs. (3)-(7),
    written out in chow_whitney_factor, give 0.40936 and 0.67199, compared
    to round-off. The printed k1^2 = 0.4098 and k2^2 = 0.6915 are not
    consistent with these data, see the NOTE after WHITNEY_1973, and they
    are not used beyond checking that they are not reproduced.
    """
    ref = WHITNEY_1973['sandwich']
    p = ref['p']
    plyts = [(1 - p)/2, p, (1 - p)/2]
    laminaprops = [WHITNEY_1973['ply'], ref['core'], WHITNEY_1973['ply']]
    for mode in ['whitney', 'vlachoutsis']:
        lam = laminated_plate([0, 0, 0], plyts=plyts,
                              laminaprops=laminaprops, shear_correction=mode)
        k1, k2 = _orthotropic_factors([0, 0, 0], plyts, laminaprops)
        assert np.isclose(lam.scf_k13, k1, rtol=1e-12, atol=0)
        assert np.isclose(lam.scf_k23, k2, rtol=1e-12, atol=0)
    # NOTE 4.4 units of the last digit below 0.4098 and 2.8 % below 0.6915
    assert abs(lam.scf_k13 - float(ref['k13'])) > 0.5e-4
    assert abs(lam.scf_k23 - float(ref['k23'])) > 0.5e-4


@pytest.mark.parametrize('seed', SEEDS)
def test_whitney_equals_vlachoutsis(seed):
    lams = [random_laminate(seed) for _ in range(2)]
    for lam, mode in zip(lams, ['whitney', 'vlachoutsis']):
        lam.shear_correction = mode
        lam.calc_transverse_shear_stiffness()
    assert np.allclose(lams[0].Ats, lams[1].Ats, rtol=1e-14, atol=0)


# A posteriori recovery (Noor and Peters, 1989)
# ---------------------------------------------

@pytest.mark.parametrize('seed', SEEDS[:6])
def test_equilibrium_recovery_cylindrical_bending(seed):
    r"""With the gradients of the two cylindrical bending states of Rohwer
    (1988), the recovered stresses are the distribution f(z) Q of
    calc_transverse_shear_stress, and the top face is traction free"""
    lam = random_laminate(seed)
    Q = np.array([0.7, -1.3]) # Qy, Qx
    Hs = np.linalg.inv(lam.ABD)
    Lx = np.zeros((6, 2)); Lx[3, 1] = 1.; Lx[5, 0] = 1.
    Ly = np.zeros((6, 2)); Ly[4, 0] = 1.; Ly[5, 1] = 1.
    # tau_xz from the x state only, tau_yz from the y state only
    zi = _interfaces(lam)
    zs = np.linspace(zi[0], zi[-1], 17)
    _, _, txz = lam.calc_equilibrium_transverse_shear(Hs @ Lx @ Q,
            np.zeros(6), z=zs)
    _, tyz, _ = lam.calc_equilibrium_transverse_shear(np.zeros(6),
            Hs @ Ly @ Q, z=zs)
    for j, z in enumerate(zs):
        ref = lam.calc_transverse_shear_stress(z, Q[0], Q[1])
        assert np.isclose(tyz[0, j], ref[0], atol=1e-9*np.abs(Q).max()/lam.h)
        assert np.isclose(txz[0, j], ref[1], atol=1e-9*np.abs(Q).max()/lam.h)
    assert abs(txz[0, -1]) < 1e-9*np.abs(Q).max()/lam.h


def test_aposteriori_energy_homogeneous():
    r"""Homogeneous plate in cylindrical bending: parabola, energy
    3/5 Qx^2/(G h), i.e. k = 5/6 with Eq. (8) of Noor and Peters (1989)"""
    h, E, nu = 2., 70e3, 0.3
    G = E/(2*(1 + nu))
    lam = isotropic_plate(thickness=h, E=E, nu=nu)
    Qx = 1.7
    D = E*h**3/(12*(1 - nu**2))
    eps_x = np.zeros(6)
    # M_xx,x = Qx with the plane strain curvature of cylindrical bending
    eps_x[3] = Qx/D
    U13, U23, Qy_r, Qx_r = lam.calc_aposteriori_energy(eps_x, np.zeros(6))
    assert np.isclose(Qx_r[0], Qx, rtol=1e-12)
    assert np.isclose(U13[0], 0.6*Qx**2/(G*h), rtol=1e-12)
    assert np.isclose(0.5*Qx**2/(G*h)/U13[0], 5/6, rtol=1e-12)
    assert abs(U23[0]) < 1e-20


# Noor, A. K., "Stability of multilayered composite plates", Fibre Sci.
# Technol., 8, 81-89, 1975, Table 1, p. 83: "composite shear correction
# factors" of cross-ply laminates in cylindrical bending, following Chow and
# Whitney, with EL/ET = 30, GLT/ET = 0.6, GTT/ET = 0.5, nu_LT = 0.25, the total
# thickness of the 0 and 90 deg layers being the same (p. 84)
NOOR_1975 = {
    'ply': (30., 1., 0.25, 0.6, 0.6, 0.5),
    'rows': [
        # NL, kappa(1), kappa(2); even NL skew-symmetric, odd NL symmetric
        (2, '0.6421', '0.6421'), (4, '0.6523', '0.6523'),
        (6, '0.7422', '0.7422'), (10, '0.7947', '0.7947'),
        (3, '0.8274', '0.5412'), (5, '0.8732', '0.5914'),
        (7, '0.8769', '0.6749'), (9, '0.8736', '0.7170'),
    ],
}


@pytest.mark.parametrize('NL, k13, k23', NOOR_1975['rows'])
def test_noor_1975_published_factors(NL, k13, k23):
    stack = [0 if i % 2 == 0 else 90 for i in range(NL)]
    n0 = sum(1 for t in stack if t == 0)
    n90 = NL - n0
    plyts = [0.5/n0 if t == 0 else 0.5/n90 for t in stack]
    lam = laminated_plate(stack, plyts=plyts, laminaprop=NOOR_1975['ply'],
            shear_correction='whitney')
    assert_printed(lam.scf_k13, k13)
    assert_printed(lam.scf_k23, k23)


# Error branches
# --------------

def test_thickness_shear_frequency():
    from composites.core import _thickness_shear_frequency
    # homogeneous layer with traction-free faces: omega = pi/h sqrt(G/rho)
    omega = _thickness_shear_frequency([4.], [1.], [2.], 1.)
    assert np.isclose(omega, np.pi/2.*np.sqrt(4./1.), rtol=1e-14)
    # the search stops at 100 times the reference frequency
    with pytest.raises(RuntimeError, match='No thickness-shear frequency'):
        _thickness_shear_frequency([1.], [1.], [1.], 1e-3)


def test_empty_laminate_stress():
    with pytest.raises(ValueError, match='0 plies'):
        Laminate().calc_transverse_shear_stress(0., 1., 0.)


def test_ill_conditioned_ABD():
    # A = h E/(1 - nu^2) [[1, nu], [nu, 1]], singular for nu -> 1
    nu = 1 - 1e-12
    with pytest.raises(ValueError, match='ill-conditioned'):
        isotropic_plate(thickness=1., E=70e3, nu=nu)
    lam = isotropic_plate(thickness=1., E=70e3, nu=nu,
                          shear_correction='constant')
    assert np.isclose(lam.A55, 5/6*70e3/(2*(1 + nu)), rtol=1e-12)


def test_non_finite_Cs():
    lam = laminated_plate([0, 90, 0], plyt=0.25, laminaprop=CFRP,
                          shear_correction=None)
    lam.plies[1].q44L = np.nan
    lam.shear_correction = 'rohwer'
    with pytest.raises(ValueError, match='Singular transverse shear'):
        lam.calc_transverse_shear_stiffness()


@pytest.mark.parametrize('mode', ['vlachoutsis', 'whitney', 'chow',
                                  'birman_bert'])
def test_non_positive_in_plane_stiffness(mode):
    lam = laminated_plate([0, 90, 0], plyt=0.25, laminaprop=CFRP,
                          shear_correction=None)
    for ply in lam.plies:
        ply.q11L = 0.
        ply.q22L = 0.
    lam.shear_correction = mode
    with pytest.raises(ValueError, match='positive in-plane stiffnesses'):
        lam.calc_transverse_shear_stiffness()


@pytest.mark.parametrize('mode', ['whitney', 'chow', 'birman_bert',
                                  'thickness_shear'])
def test_singular_Cs_scalar_modes(mode):
    laminaprops = [CFRP, (138., 9.3, 0.3, 4.6, 0., 2.3), CFRP]
    # positive densities, such that 'thickness_shear' reaches the Cs check
    lam = laminated_plate([0, 0, 0], plyt=0.25, laminaprops=laminaprops,
                          rho=1.6e-6, shear_correction=None)
    lam.shear_correction = mode
    with pytest.raises(ValueError, match='Ply 1'):
        lam.calc_transverse_shear_stiffness()


def test_equilibrium_z_outside():
    lam = laminated_plate([0, 90], plyt=0.25, laminaprop=CFRP, offset=0.1)
    zero = np.zeros(6)
    # the faces themselves are accepted
    z, _, _ = lam.calc_equilibrium_transverse_shear(zero, zero,
                                                    z=[-0.15, 0.35])
    assert np.allclose(z, [-0.15, 0.35])
    for z in ([0.36], [-0.16]):
        with pytest.raises(ValueError, match='outside'):
            lam.calc_equilibrium_transverse_shear(zero, zero, z=z)


def test_aposteriori_energy_singular_Cs():
    laminaprops = [CFRP, (138., 9.3, 0.3, 4.6, 0., 2.3), CFRP]
    lam = _manual_laminate([0, 0, 0], 0.25, laminaprops, None)
    eps_x = np.zeros(6)
    eps_x[3] = 1.
    with pytest.raises(ValueError, match='Ply 1'):
        lam.calc_aposteriori_energy(eps_x, np.zeros(6))


def test_deprecated_grad_x_grad_y():
    lam = laminated_plate([0, 90, 45, -45, 30], plyt=0.2, laminaprop=CFRP,
                          offset=0.1)
    rng = np.random.default_rng(1)
    eps_x = rng.normal(size=(3, 6))
    eps_y = rng.normal(size=(3, 6))
    zs = np.linspace(-0.4, 0.6, 7)
    ref = lam.calc_equilibrium_transverse_shear(eps_x, eps_y, zs)
    ref_kw = lam.calc_equilibrium_transverse_shear(eps_x=eps_x, eps_y=eps_y,
                                                   z=zs)
    for r, a in zip(ref, ref_kw):
        assert np.array_equal(r, a)
    with pytest.warns(DeprecationWarning) as record:
        old = lam.calc_equilibrium_transverse_shear(grad_x=eps_x,
                                                    grad_y=eps_y, z=zs)
    assert len(record) == 2
    assert "'grad_x' is deprecated, use 'eps_x'" in str(record[0].message)
    assert "'grad_y' is deprecated, use 'eps_y'" in str(record[1].message)
    for r, a in zip(ref, old):
        assert np.array_equal(r, a)
    with pytest.warns(DeprecationWarning, match="'grad_y'"):
        mixed = lam.calc_equilibrium_transverse_shear(eps_x, grad_y=eps_y,
                                                      z=zs)
    for r, a in zip(ref, mixed):
        assert np.array_equal(r, a)
    U = lam.calc_aposteriori_energy(eps_x, eps_y)
    with pytest.warns(DeprecationWarning) as record:
        U_old = lam.calc_aposteriori_energy(grad_x=eps_x, grad_y=eps_y)
    assert len(record) == 2
    for r, a in zip(U, U_old):
        assert np.array_equal(r, a)
    for f in (lam.calc_equilibrium_transverse_shear,
              lam.calc_aposteriori_energy):
        with pytest.raises(TypeError, match="both 'eps_x'"):
            f(eps_x, eps_y, grad_x=eps_x)
        with pytest.raises(TypeError, match="both 'eps_y'"):
            f(eps_x, eps_y=eps_y, grad_y=eps_y)
        with pytest.raises(TypeError, match="missing the argument 'eps_x'"):
            f(eps_y=eps_y)
        with pytest.raises(TypeError, match="missing the argument 'eps_y'"):
            f(eps_x)


# Large offsets
# -------------

@pytest.mark.parametrize('mode', ['rohwer', 'vlachoutsis', 'whitney', 'chow',
                                  'birman_bert', 'thickness_shear',
                                  'constant', None])
def test_large_offset(mode):
    # computed about the mid-surface, hence identical for any offset, while
    # the round-off of the heights from the reference surface used to grow
    # with (offset/h)^2, e.g. 2e-3 for 'rohwer' at offset/h = 1e4
    stack = [0, 90, 45, 90, 0] if mode == 'chow' else [0, 90, 45, -45, 30]
    Ats = [laminated_plate(stack, plyt=0.2, laminaprop=CFRP, rho=1.6e-6,
                           offset=offset, shear_correction=mode).Ats
           for offset in [0., 1e2, -1e4, 1e6]]
    for other in Ats[1:]:
        assert np.array_equal(other, Ats[0])


def test_large_offset_stresses():
    stack = [0, 90, 45, -45, 30]
    lam0 = laminated_plate(stack, plyt=0.2, laminaprop=CFRP)
    zs = np.linspace(-0.5, 0.5, 11)
    tau0 = np.array([lam0.calc_transverse_shear_stress(z, 0.7, -1.3)
                     for z in zs])
    rng = np.random.default_rng(0)
    eps_x = rng.normal(size=6)
    eps_y = rng.normal(size=6)
    _, tyz0, txz0 = lam0.calc_equilibrium_transverse_shear(eps_x, eps_y,
                                                           z=zs)
    U0 = np.array(lam0.calc_aposteriori_energy(eps_x, eps_y))
    for offset in [1e2, -1e4, 1e6]:
        lam = laminated_plate(stack, plyt=0.2, laminaprop=CFRP,
                              offset=offset)
        # only the representation of z + offset is affected by round-off
        tol = 1e-14*abs(offset)/lam.h
        tau = np.array([lam.calc_transverse_shear_stress(z + offset, 0.7,
                                                         -1.3)
                        for z in zs])
        assert np.allclose(tau, tau0, rtol=0, atol=tol*np.abs(tau0).max())
        # same strain field, eps0 + z eps1, about the new reference surface
        eps_x0 = eps_x.copy()
        eps_y0 = eps_y.copy()
        eps_x0[:3] -= offset*eps_x0[3:]
        eps_y0[:3] -= offset*eps_y0[3:]
        z, tyz, txz = lam.calc_equilibrium_transverse_shear(eps_x0, eps_y0,
                                                            z=zs + offset)
        assert np.array_equal(z, zs + offset)
        assert np.allclose(txz, txz0, rtol=0, atol=tol*np.abs(txz0).max())
        assert np.allclose(tyz, tyz0, rtol=0, atol=tol*np.abs(tyz0).max())
        U = np.array(lam.calc_aposteriori_energy(eps_x0, eps_y0))
        assert np.allclose(U, U0, rtol=0, atol=tol*np.abs(U0).max())


def test_example_shear_stress_cross_ply():
    # [0/90/0] CFRP plate, SI units, under the shear force Qx only
    laminaprop = (142e9, 8.7e9, 0.28, 5.1e9, 5.1e9, 3.2e9)
    lam = laminated_plate([0, 90, 0], plyt=0.5e-3, laminaprop=laminaprop)
    h = lam.h
    Qy, Qx = 0., 1e3 # N/m
    zs = np.linspace(-h/2, h/2, 31)
    tau = np.array([lam.calc_transverse_shear_stress(z, Qy, Qx) for z in zs])
    tau_yz, tau_xz = tau[:, 0], tau[:, 1]
    tol = 1e-10*Qx/h
    # traction-free faces
    assert abs(tau_xz[0]) < tol and abs(tau_xz[-1]) < tol
    # orthotropic plies in the laminate axes: Qx gives no tau_yz
    assert np.allclose(tau_yz, 0, atol=tol)
    # the stresses integrate to the shear forces, and their complementary
    # energy defines Ats: 1/2 Q^T Ats^-1 Q = int 1/2 tau^T Cs^-1 tau dz
    Q = np.array([Qy, Qx])
    xi, w = np.polynomial.legendre.leggauss(3)
    zk = -h/2
    resultant = np.zeros(2)
    Ustar = 0.
    for ply in lam.plies:
        Cs = np.asarray(ply.get_constitutive_matrix())[3:, 3:]
        for x, wi in zip(xi, w):
            z = zk + ply.h*(x + 1)/2
            t = np.array(lam.calc_transverse_shear_stress(z, Qy, Qx))
            resultant += wi*ply.h/2*t
            Ustar += wi*ply.h/2*0.5*t @ np.linalg.solve(Cs, t)
        zk += ply.h
    assert np.allclose(resultant, Q, rtol=1e-12, atol=1e-12*Qx)
    assert np.isclose(Ustar, 0.5*Q @ np.linalg.solve(lam.Ats, Q), rtol=1e-12)
    # the constant-strain stiffness overestimates the shear stiffness
    assert lam.scf_k13 < 1
