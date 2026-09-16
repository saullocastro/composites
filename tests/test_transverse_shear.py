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
    lam = laminated_plate(stack, plyt=1/120, laminaprop=laminaprop)
    lam.shear_correction = 'vlachoutsis'
    lam.calc_transverse_shear_stiffness()
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


def test_calc_scf_returns_ratios():
    lam = laminated_plate([0, 90, 90, 0], plyt=0.25, laminaprop=CFRP)
    k13, k23 = lam.calc_scf()
    assert (k13, k23) == (lam.scf_k13, lam.scf_k23)
    assert abs(k13 - 0.6705) < 5e-5
    assert abs(k23 - 0.7320) < 5e-5


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
