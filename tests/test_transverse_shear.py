r"""
Transverse shear stiffness of FSDT laminates

The reference values were produced with ``transverse_shear_reference.py``,
which is an independent NumPy implementation of the same equations.

Conventions: ``Ats = [[A44, A45], [A45, A55]]``, 4 <-> yz, 5 <-> xz.
"""
import sys
sys.path.append('..')

import numpy as np
import pytest

from composites import laminated_plate, isotropic_plate
from composites.utils import read_laminaprop
from composites.core import Laminate, Lamina

from transverse_shear_reference import (stack_from_laminate, rohwer_Ats,
        Ats_constant_strain, Ats_constant_stress)


CFRP = (138., 9.3, 0.3, 4.6, 4.6, 2.3) # kN/mm^2, Rohwer (1988)


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


def sandwich(offset=0., stack=(0., 0., 0.)):
    Ef, nu = 70., 0.3
    laminaprops = [(Ef, nu), (Ef*1e-3, nu), (Ef, nu)]
    return laminated_plate(list(stack), plyts=[0.1, 0.8, 0.1],
            laminaprops=laminaprops, offset=offset)


# Exactness
# ---------

@pytest.mark.parametrize('nu', [0.0, 0.3, 0.49])
@pytest.mark.parametrize('offset', [0., +0.37, -1.5])
def test_isotropic_five_sixths(nu, offset):
    h = 2.
    E = 70e3
    G = E/(2*(1 + nu))
    lam = isotropic_plate(thickness=h, E=E, nu=nu, offset=offset)
    assert np.isclose(lam.A55, 5/6*G*h, rtol=1e-12, atol=0)
    assert np.isclose(lam.A44, 5/6*G*h, rtol=1e-12, atol=0)
    assert abs(lam.A45) <= 1e-12*G*h
    assert np.isclose(lam.scf_k13, 5/6, rtol=1e-12, atol=0)
    assert np.isclose(lam.scf_k23, 5/6, rtol=1e-12, atol=0)
    assert np.allclose(lam.Abar_ts, G*h*np.eye(2), rtol=1e-14)
    assert np.allclose(lam.Abarbar_ts, G*h*np.eye(2), rtol=1e-14)


# Published values
# ----------------

def test_rohwer_published_cross_ply():
    lam = laminated_plate([0, 90, 90, 0], plyt=0.25, laminaprop=CFRP)
    tol = 5e-5
    assert abs(lam.Abar55 - 3.4500) < tol
    assert abs(lam.Abar44 - 3.4500) < tol
    assert abs(lam.Abarbar55 - 3.0667) < tol
    # Rohwer (1988) publishes 2.313
    assert abs(lam.A55 - 2.3131) < tol
    # Rohwer (1988) publishes 2.521, the 0.17 % gap is the paper's rounding
    # of the input data
    assert abs(lam.A44 - 2.5252) < tol
    assert abs(lam.A45) < 1e-12
    assert abs(lam.scf_k13 - 0.6705) < tol
    assert abs(lam.scf_k23 - 0.7320) < tol
    assert np.allclose(lam.Ats, [[lam.A44, lam.A45], [lam.A45, lam.A55]])


@pytest.mark.parametrize('ratios, k13, k23', [
    ((25., 0.5, 0.5, 0.2, 0.25), 0.6808, 0.6794),
    ((40., 0.6, 0.6, 0.5, 0.25), 0.8316, 0.8208),
    ])
def test_vlachoutsis_published_120_plies(ratios, k13, k23):
    E1, G12, G13, G23, nu12 = ratios
    laminaprop = (E1, 1., nu12, G12, G13, G23)
    stack = [0, 90]*30 + [90, 0]*30
    lam = laminated_plate(stack, plyt=1/120, laminaprop=laminaprop,
            shear_correction='vlachoutsis')
    assert abs(lam.scf_k13 - k13) < 5e-5
    assert abs(lam.scf_k23 - k23) < 5e-5
    assert np.isclose(lam.A55, k13*lam.Abar55, atol=5e-5*lam.Abar55)
    assert np.isclose(lam.A44, k23*lam.Abar44, atol=5e-5*lam.Abar44)


@pytest.mark.parametrize('stack, A45', [
    ([45, -45, -45, 45], -0.1204),
    ([-45, 45, 45, -45], +0.1204),
    ])
def test_rohwer_angle_ply_coupling(stack, A45):
    lam = laminated_plate(stack, plyt=0.25, laminaprop=CFRP)
    # the ad-hoc (k13 + k23)/2*Abar45 cannot represent this coupling
    assert abs(lam.Abar45) < 1e-12
    assert abs(lam.A45 - A45) < 5e-5
    assert abs(lam.scf_k13 - 0.7503) < 5e-5
    assert abs(lam.scf_k23 - 0.7503) < 5e-5


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


# Regression against the reference implementation
# ------------------------------------------------

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
        cases.append(('sandwich', offset, sandwich(offset)))
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
    ref = rohwer_Ats(st).Ats
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
    assert np.allclose(lam.Ats, rohwer_Ats(stack_from_laminate(lam)).Ats)
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
    lam = laminated_plate([0, 90, 90, 0], plyt=0.25, laminaprop=CFRP,
            shear_correction=None)
    with pytest.warns(DeprecationWarning, match='calc_scf'):
        k13, k23 = lam.calc_scf()
    assert lam.shear_correction == 'rohwer'
    assert (k13, k23) == (lam.scf_k13, lam.scf_k23)
    assert abs(k13 - 0.6705) < 5e-5
    assert abs(k23 - 0.7320) < 5e-5
    assert abs(lam.A55 - 2.3131) < 5e-5


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


def _interfaces(lam):
    return lam.offset - lam.h/2 + np.concatenate(([0.],
        np.cumsum([ply.h for ply in lam.plies])))


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


def _gauss_integral(lam, func):
    zi = _interfaces(lam)
    xi, w = np.polynomial.legendre.leggauss(3)
    out = 0
    for k, ply in enumerate(lam.plies):
        zm, dz = (zi[k] + zi[k+1])/2, (zi[k+1] - zi[k])/2
        for x, wi in zip(xi, w):
            out = out + wi*dz*func(k, ply, zm + dz*x)
    return out


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
    ref = rohwer_Ats(stack_from_laminate(lam))
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
