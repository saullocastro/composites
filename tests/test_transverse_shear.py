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

# Rohwer, K. "Improved transverse shear stiffness for layered finite elements",
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

    # Section 4.3, pp. 22-23: solar panel of Fig. 9, p. 22, with "A foam core
    # of 9mm thickness is covered with faces from CFRP, layered [0,90,+45,-45],
    # where the -45deg-layers are adhered to the core. A layer thickness of
    # 0.125mm leads to a total thickness of 10mm." The stack is symmetric, the
    # angle of the isotropic core is irrelevant.
    'solar_panel': dict(stack=[0, 90, 45, -45, 0, -45, 45, 90, 0],
        plyts=[0.125]*4 + [9.] + [0.125]*4,
        # Table 3, p. 22: "Foam material constants", kN/mm^2
        E_foam=0.1, nu_foam=0.3,
        # Table 4, p. 23: "Sandwich stiffnesses", kNmm and kN/mm
        D11='1366.08', D12='395.09', D22='1321.97', D13='9.36', D23='9.36',
        H11='0.416', H12='0.008', H22='0.403'),
}

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

@pytest.mark.parametrize('mode', ['rohwer', 'vlachoutsis'])
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
    assert_printed(lam.A55, ref['H11'])
    # A44 = 2.5252, 0.17 % above the printed value
    assert_printed(lam.A44, ref['H22'], rtol=2e-3)
    assert abs(lam.A45 - ref['H12']) < 1e-12
    assert_printed(lam.A55/lam.Abarbar55, ref['H11_Hbarbar11'])
    assert_printed(lam.A44/lam.Abarbar44, ref['H22_Hbarbar22'])
    assert np.isclose(lam.scf_k13, lam.A55/lam.Abar55, rtol=1e-14)
    assert np.isclose(lam.scf_k23, lam.A44/lam.Abar44, rtol=1e-14)
    assert np.allclose(lam.Ats, [[lam.A44, lam.A45], [lam.A45, lam.A55]])


def test_rohwer_published_sandwich_minimum():
    ref = ROHWER_1988['sandwich']
    ratios = []
    # abscissa log(a_c/a_f) of Fig. 5
    for log_ac_af in np.linspace(-6., 6., 241):
        lam = sandwich(tf=1., tc=10**log_ac_af, Ef_Ec=ref['Ef_Ec'],
                nu=ref['nu'])
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


def test_rohwer_published_solar_panel():
    ref = ROHWER_1988['solar_panel']
    foam = (ref['E_foam'], ref['nu_foam'])
    laminaprops = [CFRP]*4 + [foam] + [CFRP]*4
    lam = laminated_plate(ref['stack'], plyts=ref['plyts'],
            laminaprops=laminaprops)
    assert_printed(lam.D16, ref['D13'])
    assert_printed(lam.D26, ref['D23'])
    assert_printed(lam.A45, ref['H12'])
    # up to 0.62 % above the printed values
    assert_printed(lam.D11, ref['D11'], rtol=0.01)
    assert_printed(lam.D12, ref['D12'], rtol=0.01)
    assert_printed(lam.D22, ref['D22'], rtol=0.01)
    # A55 = 0.3999 and A44 = 0.3872, 3.9 % below the printed values. Both are
    # governed by the foam shear modulus, and they round to the printed values
    # when the foam is given G = 0.04 kN/mm^2 instead of the E/(2(1 + nu)) =
    # 0.0385 kN/mm^2 that follows from Table 3
    assert_printed(lam.A55, ref['H11'], rtol=0.05)
    assert_printed(lam.A44, ref['H22'], rtol=0.05)


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
    nu = VLACHOUTSIS_1992['sandwich_isotropic']['nu']
    lam = sandwich(tf=(1 - c_h)/2, tc=c_h, Ef_Ec=E1_E2, nu=nu,
            shear_correction='vlachoutsis')
    assert_printed(lam.scf_k13, K)
    assert_printed(lam.scf_k23, K)
    kappa, share = _vlachoutsis_ply_factors(lam)
    assert_printed(kappa[0], K1)
    assert_printed(kappa[2], K1)
    assert_printed(kappa[1], K2)
    assert_printed(share[0] + share[2], rf)


@pytest.mark.parametrize('mode', ['vlachoutsis', 'rohwer'])
def test_vlachoutsis_published_spherical_shell(mode):
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
    # antiplane core, Dc/Df -> 0
    Ec = 1e-9*Ef
    core = (Ec, Ec, 0., Gc, Gc, Gc)
    lam = laminated_plate([0, 0, 0], plyts=[t, c, t],
            laminaprops=[(Ef, nu), core, (Ef, nu)], shear_correction=mode)
    # 0.077684, 0.01 % below the printed value
    assert_printed(lam.scf_k13, ref['K'], rtol=1e-4)
    assert_printed(lam.scf_k23, ref['K'], rtol=1e-4)


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

@pytest.mark.parametrize('stack, A45', [
    ([45, -45, -45, 45], -0.1204),
    ([-45, 45, 45, -45], +0.1204),
    ])
def test_rohwer_angle_ply_coupling(stack, A45):
    # not published, values from the reference implementation
    lam = laminated_plate(stack, plyt=0.25, laminaprop=CFRP)
    # the ad-hoc (k13 + k23)/2*Abar45 cannot represent this coupling
    assert abs(lam.Abar45) < 1e-12
    assert abs(lam.A45 - A45) < 5e-5
    assert abs(lam.scf_k13 - 0.7503) < 5e-5
    assert abs(lam.scf_k23 - 0.7503) < 5e-5


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
    lam = sandwich()
    # the previous implementation returned 1.2564
    assert abs(lam.scf_k13 - 0.0051) < 5e-5
    assert abs(lam.scf_k23 - 0.0051) < 5e-5


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
