import sys
sys.path.append('..')

import copy
import pickle
import types
import warnings

import numpy as np
import pytest

from composites.utils import laminated_plate
from composites.core import (MatLamina, Lamina, Laminate,
                             LaminationParameters, GradABD,
                             _LAMINATE_STATE)


laminaprop = (142e9, 8.7e9, 0.28, 5.1e9, 5.1e9, 3.2e9)
LAMINATE_MATRICES = ('A', 'B', 'D', 'E', 'F', 'H', 'ABD', 'Ats', 'Abar_ts',
                     'Abarbar_ts', 'Dtrans', 'Ftrans')


def roundtrip(obj):
    return pickle.loads(pickle.dumps(obj, protocol=pickle.HIGHEST_PROTOCOL))


def public_attrs(obj):
    """Names of the writable attributes, i.e. the "cdef public" members"""
    names = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', DeprecationWarning)
        for name, d in vars(type(obj)).items():
            if not isinstance(d, types.GetSetDescriptorType):
                continue
            try:
                setattr(obj, name, getattr(obj, name))
            except AttributeError:
                continue
            names.append(name)
    return names


def assert_same(a, b):
    assert type(a) is type(b)
    if isinstance(a, (MatLamina, Lamina, Laminate, LaminationParameters)):
        assert a is not b
        names = public_attrs(a)
        assert names
        for name in names:
            assert_same(getattr(a, name), getattr(b, name))
        if isinstance(a, Laminate):
            for name in LAMINATE_MATRICES:
                np.testing.assert_array_equal(getattr(a, name),
                                              getattr(b, name))
    elif isinstance(a, list):
        assert len(a) == len(b)
        for ai, bi in zip(a, b):
            assert_same(ai, bi)
    else:
        np.testing.assert_array_equal(a, b)


def make_laminate(shear_correction='rohwer', offset=0.):
    lam = laminated_plate([0, 45, -45, 90, 90, 30], plyt=1e-3,
            laminaprop=laminaprop, rho=1600., offset=offset,
            shear_correction=shear_correction)
    lam.calc_equivalent_properties()
    return lam


def test_laminate_state_names_complete():
    # every "cdef public" attribute must be saved when pickling
    assert sorted(public_attrs(Laminate())) == sorted(_LAMINATE_STATE)


def test_matlamina():
    lam = make_laminate()
    mat = lam.plies[0].matlamina
    assert_same(mat, roundtrip(mat))
    assert_same(MatLamina(), roundtrip(MatLamina()))


def test_lamina():
    lam = make_laminate()
    ply = lam.plies[1]
    ply2 = roundtrip(ply)
    assert_same(ply, ply2)
    np.testing.assert_array_equal(ply.get_constitutive_matrix(),
                                  ply2.get_constitutive_matrix())
    assert_same(Lamina(), roundtrip(Lamina()))


def test_lamination_parameters():
    lp = make_laminate().calc_lamination_parameters()
    assert_same(lp, roundtrip(lp))
    assert_same(LaminationParameters(), roundtrip(LaminationParameters()))


@pytest.mark.parametrize('shear_correction, offset', [
    ('rohwer', 0.), ('rohwer', 1.5e-3), ('constant', 0.), (None, 0.)])
def test_laminate(shear_correction, offset):
    lam = make_laminate(shear_correction, offset)
    lam2 = roundtrip(lam)
    assert_same(lam, lam2)
    assert lam2.shear_correction == shear_correction
    # plies sharing the same material still share it after unpickling
    lam.plies[-1].matlamina = lam.plies[0].matlamina
    lam2 = roundtrip(lam)
    assert lam2.plies[0].matlamina is lam2.plies[-1].matlamina


def test_laminate_default():
    lam = Laminate()
    lam2 = roundtrip(lam)
    assert_same(lam, lam2)
    assert lam2.plies == [] and lam2.stack == []
    assert lam2.shear_correction == 'rohwer'


@pytest.mark.parametrize('shear_correction', ['rohwer', 'constant'])
def test_laminate_transverse_shear_cache(shear_correction):
    lam = make_laminate(shear_correction, offset=0.5e-3)
    zs = np.linspace(-lam.h/2 + lam.offset, lam.h/2 + lam.offset, 23)
    # fills the transverse shear distribution cache before pickling
    ref_stress = [lam.calc_transverse_shear_stress(z, 3., -2.) for z in zs]
    lam2 = roundtrip(lam)
    stress = [lam2.calc_transverse_shear_stress(z, 3., -2.) for z in zs]
    np.testing.assert_array_equal(stress, ref_stress)

    lam2 = roundtrip(lam)
    lam.calc_transverse_shear_stiffness()
    lam2.calc_transverse_shear_stiffness()
    assert_same(lam, lam2)

    lam2 = roundtrip(lam)
    with pytest.warns(DeprecationWarning):
        ref_scf = lam.calc_scf()
    with pytest.warns(DeprecationWarning):
        scf = lam2.calc_scf()
    np.testing.assert_array_equal(scf, ref_scf)
    assert_same(lam, lam2)


def test_laminate_deepcopy():
    lam = make_laminate()
    lam.calc_transverse_shear_stress(0., 1., 1.)
    lam2 = copy.deepcopy(lam)
    assert_same(lam, lam2)
    assert lam2.plies[0] is not lam.plies[0]
    np.testing.assert_array_equal(lam2.calc_transverse_shear_stress(0., 1., 1.),
                                  lam.calc_transverse_shear_stress(0., 1., 1.))
    assert_same(Laminate(), copy.deepcopy(Laminate()))


class LaminateSubclass(Laminate):
    pass


def test_laminate_subclass():
    lam = LaminateSubclass()
    lam.h = 2.
    lam.name = 'skin'
    lam2 = roundtrip(lam)
    assert type(lam2) is LaminateSubclass
    assert lam2.h == 2.
    assert lam2.name == 'skin'


def assert_same_gradabd(a, b):
    for name in ('gradAij', 'gradBij', 'gradDij', 'gradAtransij'):
        np.testing.assert_array_equal(getattr(a, name), getattr(b, name))


def test_gradabd():
    lam = make_laminate()
    grad = GradABD()
    grad.calc_LP_grad(lam.h, lam.plies[0].matlamina,
                      lam.calc_lamination_parameters())
    assert np.any(np.asarray(grad.gradDij) != 0)
    grad2 = roundtrip(grad)
    assert_same_gradabd(grad, grad2)
    assert_same_gradabd(grad, copy.deepcopy(grad))
    assert_same_gradabd(GradABD(), roundtrip(GradABD()))
    # restored arrays are independent and usable by calc_LP_grad
    grad2.calc_LP_grad(lam.h, lam.plies[0].matlamina,
                       lam.calc_lamination_parameters())
    assert_same_gradabd(grad, grad2)
