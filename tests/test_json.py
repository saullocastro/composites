import sys
sys.path.append('..')

import io
import json

import numpy as np
import pytest

import composites
from composites import (laminated_plate, to_dict, from_dict, to_json,
                        from_json, save_json, load_json)
from composites.core import (MatLamina, Lamina, Laminate,
                             LaminationParameters, GradABD,
                             laminate_from_lamination_parameters)
from composites.json_io import (_MATLAMINA_FIELDS, _LAMINA_FIELDS,
                                _LP_FIELDS, FORMAT_VERSION)

from test_pickle import (public_attrs, assert_same, make_laminate,
                         laminaprop)


def roundtrip(obj):
    s = to_json(obj)
    # strict JSON, e.g. no NaN, readable by JSON.parse in JavaScript
    def reject(name):
        raise AssertionError('non-standard JSON constant %s' % name)
    json.loads(s, parse_constant=reject)
    return from_json(s)


@pytest.mark.parametrize('cls, fields', [
    (MatLamina, _MATLAMINA_FIELDS),
    (Lamina, _LAMINA_FIELDS + ('plyid', 'matlamina')),
    (LaminationParameters, _LP_FIELDS)])
def test_field_names_complete(cls, fields):
    # every "cdef public" attribute must be saved
    assert sorted(public_attrs(cls())) == sorted(fields)


def test_matlamina():
    mat = make_laminate().plies[0].matlamina
    assert_same(mat, roundtrip(mat))
    assert_same(MatLamina(), roundtrip(MatLamina()))


def test_lamina():
    ply = make_laminate().plies[1]
    ply.plyid = 7
    ply2 = roundtrip(ply)
    assert_same(ply, ply2)
    np.testing.assert_array_equal(ply.get_constitutive_matrix(),
                                  ply2.get_constitutive_matrix())
    ply3 = roundtrip(Lamina())
    assert_same(Lamina(), ply3)
    assert ply3.matlamina is None


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
    zs = np.linspace(-lam.h/2 + lam.offset, lam.h/2 + lam.offset, 11)
    for z in zs:
        assert (lam2.calc_transverse_shear_stress(z, 3., -2.)
                == lam.calc_transverse_shear_stress(z, 3., -2.))


def test_laminate_materials():
    lam = make_laminate()
    # laminated_plate creates one material per ply
    d = to_dict(lam)
    assert len(d['data']['materials']) == len(lam.plies)
    lam2 = from_dict(d)
    assert lam2.plies[0].matlamina is not lam2.plies[1].matlamina
    # a shared material is stored once and still shared after loading
    for ply in lam.plies:
        ply.matlamina = lam.plies[0].matlamina
    d = to_dict(lam)
    assert len(d['data']['materials']) == 1
    lam2 = from_dict(d)
    assert_same(lam, lam2)
    assert all(ply.matlamina is lam2.plies[0].matlamina for ply in lam2.plies)
    # a ply without material
    lam.plies[2].matlamina = None
    lam2 = roundtrip(lam)
    assert lam2.plies[2].matlamina is None
    assert lam2.plies[3].matlamina is lam2.plies[0].matlamina


def test_laminate_default():
    lam2 = roundtrip(Laminate())
    assert_same(Laminate(), lam2)
    assert lam2.plies == [] and lam2.stack == []
    assert lam2.shear_correction == 'rohwer'


def test_laminate_from_lamination_parameters():
    mat = make_laminate().plies[0].matlamina
    lam = laminate_from_lamination_parameters(6e-3, mat, 0.1, 0.2, 0.,
                                              0., 0., 0., 0., 0., 0.3, -0.1,
                                              0., 0.)
    assert lam.plies == []
    assert_same(lam, roundtrip(lam))


def test_laminate_stack_numpy():
    lam = laminated_plate(np.array([0, 90, 90, 0]), plyt=1e-3,
                          laminaprop=laminaprop)
    lam.stack = list(np.array([0, 90, 90, 0]))
    lam2 = roundtrip(lam)
    assert lam2.stack == [0, 90, 90, 0]
    assert all(type(angle) is int for angle in lam2.stack)
    lam.stack = [0., 45.5]
    assert roundtrip(lam).stack == [0., 45.5]


def test_laminate_non_finite():
    lam = make_laminate()
    lam.Abarbar44 = np.nan
    lam.Abarbar55 = np.inf
    lam.Abarbar45 = -np.inf
    s = to_json(lam)
    assert '"NaN"' in s and '"Infinity"' in s and '"-Infinity"' in s
    lam2 = roundtrip(lam)
    assert np.isnan(lam2.Abarbar44)
    assert lam2.Abarbar55 == np.inf
    assert lam2.Abarbar45 == -np.inf


def test_gradabd():
    lam = make_laminate()
    grad = GradABD()
    grad.calc_LP_grad(lam.h, lam.plies[0].matlamina,
                      lam.calc_lamination_parameters())
    assert np.any(np.asarray(grad.gradDij) != 0)
    grad2 = roundtrip(grad)
    for name in ('gradAij', 'gradBij', 'gradDij', 'gradAtsij'):
        np.testing.assert_array_equal(getattr(grad, name),
                                      getattr(grad2, name))
    # restored arrays are usable by calc_LP_grad
    grad2.calc_LP_grad(lam.h, lam.plies[0].matlamina,
                       lam.calc_lamination_parameters())


def test_dict_is_json_types():
    def check(value):
        if isinstance(value, dict):
            for k, v in value.items():
                assert type(k) is str
                check(v)
        elif isinstance(value, list):
            for v in value:
                check(v)
        else:
            assert value is None or type(value) in (str, int, float)
    for obj in (make_laminate(), GradABD(), Lamina(), MatLamina(),
                LaminationParameters()):
        d = to_dict(obj)
        check(d)
        assert d['type'] == type(obj).__name__
        assert d['format_version'] == FORMAT_VERSION
        assert d['composites_version'] == composites.__version__


def test_save_load(tmp_path):
    lam = make_laminate()
    fname = tmp_path / 'laminate.json'
    save_json(lam, fname)
    assert_same(lam, load_json(fname))
    assert_same(lam, load_json(str(fname)))
    with open(fname) as f:
        assert json.load(f)['type'] == 'Laminate'
    buf = io.StringIO()
    save_json(lam, buf, indent=None)
    assert '\n' not in buf.getvalue()
    buf.seek(0)
    assert_same(lam, load_json(buf))


def test_errors():
    with pytest.raises(TypeError):
        to_dict([1, 2])
    d = to_dict(make_laminate())
    with pytest.raises(ValueError, match='newer'):
        from_dict(dict(d, format_version=FORMAT_VERSION + 1))
    with pytest.raises(ValueError, match='format_version'):
        from_dict(dict(d, format_version='1'))
    with pytest.raises(ValueError, match='Unknown type'):
        from_dict(dict(d, type='Laminate2'))
    with pytest.raises(ValueError, match='Unknown keys'):
        from_dict(dict(d, extra=1))
    with pytest.raises(ValueError, match='Unknown keys for Laminate: A99'):
        from_dict(dict(d, data=dict(d['data'], A99=1.)))
    with pytest.raises(ValueError, match='Expected a number'):
        from_dict(dict(d, data=dict(d['data'], A11=None)))
    # only the strings of the non-finite values are accepted
    for value in ('abc', '1.0', 'nan', 'inf', '-inf', 'infinity', ''):
        with pytest.raises(ValueError, match='"NaN", "Infinity" or'):
            from_dict(dict(d, data=dict(d['data'], A11=value)))
    with pytest.raises(ValueError, match='"NaN", "Infinity" or'):
        from_dict(dict(d, data=dict(d['data'], stack=[0, '90'])))
    with pytest.raises(ValueError, match='Expected a number'):
        from_dict(dict(d, data=dict(d['data'], A11=[1.])))
    with pytest.raises(ValueError, match='material index'):
        plies = [dict(d['data']['plies'][0], matlamina=99)]
        from_dict(dict(d, data=dict(d['data'], plies=plies)))
    with pytest.raises(ValueError, match='shear_correction'):
        from_dict(dict(d, data=dict(d['data'], shear_correction=1)))
    with pytest.raises(ValueError, match='Expected a JSON object'):
        from_dict([])
    g = to_dict(GradABD())
    with pytest.raises(ValueError, match='shape'):
        from_dict(dict(g, data=dict(g['data'], gradAij=[[1., 2.]])))
    with pytest.raises(ValueError, match='Expected a matrix for gradAij'):
        from_dict(dict(g, data=dict(g['data'], gradAij=5.)))
    with pytest.raises(ValueError, match='JSON object for MatLamina'):
        from_dict(dict(to_dict(MatLamina()), data=[]))
    ply = to_dict(make_laminate().plies[0])
    with pytest.raises(ValueError, match='JSON object for Lamina'):
        from_dict(dict(ply, data=[]))
    for plyid in (1.5, True, '1'):
        with pytest.raises(ValueError, match='integer for plyid'):
            from_dict(dict(ply, data=dict(ply['data'], plyid=plyid)))
    with pytest.raises(ValueError, match='JSON object for Lamina'):
        from_dict(dict(d, data=dict(d['data'], plies=[1])))


@pytest.mark.parametrize('cls', [MatLamina, Lamina, Laminate,
                                 LaminationParameters, GradABD])
def test_missing_keys_keep_defaults(cls):
    d = to_dict(cls())
    d['data'] = {}
    obj = from_dict(d)
    assert type(obj) is cls
    if cls is GradABD:
        for name in ('gradAij', 'gradBij', 'gradDij', 'gradAtsij'):
            np.testing.assert_array_equal(getattr(obj, name),
                                          getattr(cls(), name))
    else:
        assert_same(cls(), obj)
    del d['data']
    assert type(from_dict(d)) is cls
