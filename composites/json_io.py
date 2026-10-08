r"""
Saving and loading in JSON (:mod:`composites.json_io`)
==============================================================

.. currentmodule:: composites.json_io

The objects :class:`.MatLamina`, :class:`.Lamina`, :class:`.Laminate`,
:class:`.LaminationParameters` and :class:`.GradABD` are saved to and loaded
from JSON with::

    from composites import laminated_plate, save_json, load_json

    lam = laminated_plate([0, 45, -45, 90], plyt=1e-3, laminaprop=laminaprop)
    save_json(lam, 'laminate.json')
    lam2 = load_json('laminate.json')

The JSON is plain text, independent of the Python and composites versions,
and it is safe to load from untrusted sources, unlike pickle. It can also be
exchanged with JavaScript when running in the browser through Pyodide.

All public attributes are stored, including the derived ones, such as the
``ABD`` matrix and the constitutive terms of each ply, such that the loaded
object is identical to the saved one. A :class:`.MatLamina` shared by several
plies of a :class:`.Laminate` is stored once and is still shared after
loading. The transverse shear distribution used by
:meth:`.Laminate.calc_transverse_shear_stress` is not stored and is
recomputed on demand. Attributes of Python subclasses that are not public
attributes of the base class are not stored, and the loaded object has the
type of the base class.

The values ``nan``, ``inf`` and ``-inf``, e.g. ``Abarbar44`` of a ply with a
singular transverse shear matrix, are stored as the strings ``"NaN"``,
``"Infinity"`` and ``"-Infinity"``, such that the output is strict JSON,
readable by ``JSON.parse`` in JavaScript.

The format of the saved data is::

    {
      "type": "Laminate",
      "format_version": 1,
      "composites_version": "0.9.20",
      "data": {...}
    }

"""
import json
import math
import numbers
import os

import numpy as np

from .version import __version__
from .core import (MatLamina, Lamina, Laminate, LaminationParameters,
                   GradABD, _LAMINATE_STATE, _GRADABD_STATE)


FORMAT_VERSION = 1

# NOTE names of the "cdef public" attributes of each class, they must follow
#      core.pxd (checked in tests/test_json.py)
_LP_FIELDS = tuple(['xi%s%d' % (m, i) for m in 'ABD' for i in (1, 2, 3, 4)]
                   + ['xiAts1', 'xiAts2'])
_MATLAMINA_FIELDS = (
    'e1', 'e2', 'e3', 'g12', 'g13', 'g23',
    'nu12', 'nu21', 'nu13', 'nu31', 'nu23', 'nu32',
    'rho', 'a1', 'a2', 'a3', 'tref',
    'st1', 'st2', 'sc1', 'sc2', 'ss12',
    'q11', 'q12', 'q13', 'q21', 'q22', 'q23', 'q31', 'q32', 'q33',
    'q44', 'q55', 'q66',
    'c11', 'c12', 'c13', 'c22', 'c23', 'c33', 'c44', 'c55', 'c66',
    'u1', 'u2', 'u3', 'u4', 'u5', 'u6', 'u7')
_LAMINA_FIELDS = (
    'h', 'thetadeg', 'cost', 'cos2t', 'cos4t', 'sint', 'sin2t', 'sin4t',
    'q11L', 'q12L', 'q22L', 'q16L', 'q26L', 'q66L', 'q44L', 'q45L', 'q55L')
_LAMINATE_FIELDS = tuple(name for name in _LAMINATE_STATE
                         if name not in ('plies', 'stack', 'shear_correction'))
_GRADABD_SHAPES = {'gradAij': (6, 5), 'gradBij': (6, 5), 'gradDij': (6, 5),
                   'gradAtsij': (3, 3)}


def _encode_float(value):
    value = float(value)
    if math.isfinite(value):
        return value
    if math.isnan(value):
        return 'NaN'
    return 'Infinity' if value > 0 else '-Infinity'


_NON_FINITE = {'NaN': math.nan, 'Infinity': math.inf, '-Infinity': -math.inf}


def _decode_float(value):
    # NOTE only the strings written by _encode_float are accepted, such that
    #      e.g. "1.0" or "nan" are rejected as malformed
    if isinstance(value, str):
        if value in _NON_FINITE:
            return _NON_FINITE[value]
        raise ValueError('Expected a number, or "NaN", "Infinity" or '
                         '"-Infinity", got %r' % (value, ))
    if not isinstance(value, numbers.Real) or isinstance(value, bool):
        raise ValueError('Expected a number, got %r' % (value, ))
    return float(value)


def _encode_number(value):
    if isinstance(value, numbers.Integral) and not isinstance(value, bool):
        return int(value)
    return _encode_float(value)


def _decode_number(value):
    if isinstance(value, numbers.Integral) and not isinstance(value, bool):
        return int(value)
    return _decode_float(value)


def _encode_fields(obj, names):
    return {name: _encode_float(getattr(obj, name)) for name in names}


def _check_keys(data, allowed, typename):
    if not isinstance(data, dict):
        raise ValueError('Expected a JSON object for %s, got %r'
                         % (typename, type(data).__name__))
    unknown = set(data) - set(allowed)
    if unknown:
        raise ValueError('Unknown keys for %s: %s'
                         % (typename, ', '.join(sorted(unknown))))


def _decode_fields(obj, data, names):
    for name in names:
        if name in data:
            setattr(obj, name, _decode_float(data[name]))


def _matlamina_to_data(mat):
    return _encode_fields(mat, _MATLAMINA_FIELDS)


def _matlamina_from_data(data):
    _check_keys(data, _MATLAMINA_FIELDS, 'MatLamina')
    mat = MatLamina()
    _decode_fields(mat, data, _MATLAMINA_FIELDS)
    return mat


def _lamina_to_data(ply, matlamina):
    data = {'plyid': int(ply.plyid)}
    data.update(_encode_fields(ply, _LAMINA_FIELDS))
    data['matlamina'] = matlamina
    return data


def _lamina_from_data(data, matlamina):
    _check_keys(data, ('plyid', 'matlamina') + _LAMINA_FIELDS, 'Lamina')
    ply = Lamina()
    if 'plyid' in data:
        plyid = data['plyid']
        if not isinstance(plyid, numbers.Integral) or isinstance(plyid, bool):
            raise ValueError('Expected an integer for plyid, got %r'
                             % (plyid, ))
        ply.plyid = plyid
    _decode_fields(ply, data, _LAMINA_FIELDS)
    ply.matlamina = matlamina
    return ply


def _laminate_to_data(lam):
    # NOTE each MatLamina is stored once, plies refer to it by index
    materials = []
    index = {}
    plies = []
    for ply in lam.plies:
        mat = ply.matlamina
        if mat is None:
            plies.append(_lamina_to_data(ply, None))
            continue
        if id(mat) not in index:
            index[id(mat)] = len(materials)
            materials.append(_matlamina_to_data(mat))
        plies.append(_lamina_to_data(ply, index[id(mat)]))
    data = _encode_fields(lam, _LAMINATE_FIELDS)
    data['shear_correction'] = lam.shear_correction
    data['stack'] = [_encode_number(angle) for angle in lam.stack]
    data['materials'] = materials
    data['plies'] = plies
    return data


def _laminate_from_data(data):
    _check_keys(data, _LAMINATE_FIELDS
                + ('shear_correction', 'stack', 'materials', 'plies'),
                'Laminate')
    lam = Laminate()
    _decode_fields(lam, data, _LAMINATE_FIELDS)
    if 'shear_correction' in data:
        shear_correction = data['shear_correction']
        if shear_correction is not None and not isinstance(shear_correction,
                                                           str):
            raise ValueError('Expected a string or null for '
                             'shear_correction, got %r' % (shear_correction, ))
        lam.shear_correction = shear_correction
    lam.stack = [_decode_number(angle) for angle in data.get('stack', [])]
    materials = [_matlamina_from_data(mat)
                 for mat in data.get('materials', [])]
    plies = []
    for ply in data.get('plies', []):
        if not isinstance(ply, dict):
            raise ValueError('Expected a JSON object for Lamina, got %r'
                             % type(ply).__name__)
        i = ply.get('matlamina')
        if i is None:
            mat = None
        elif (isinstance(i, numbers.Integral) and not isinstance(i, bool)
              and 0 <= i < len(materials)):
            mat = materials[i]
        else:
            raise ValueError('Invalid material index of a ply: %r' % (i, ))
        plies.append(_lamina_from_data(ply, mat))
    lam.plies = plies
    return lam


def _lp_to_data(lp):
    return _encode_fields(lp, _LP_FIELDS)


def _lp_from_data(data):
    _check_keys(data, _LP_FIELDS, 'LaminationParameters')
    lp = LaminationParameters()
    _decode_fields(lp, data, _LP_FIELDS)
    return lp


def _gradabd_to_data(grad):
    return {name: [[_encode_float(v) for v in row]
                   for row in np.asarray(getattr(grad, name))]
            for name in _GRADABD_STATE}


def _gradabd_from_data(data):
    _check_keys(data, _GRADABD_STATE, 'GradABD')
    grad = GradABD()
    for name in _GRADABD_STATE:
        if name not in data:
            continue
        try:
            value = np.array([[_decode_float(v) for v in row]
                              for row in data[name]], dtype=np.float64)
        except TypeError:
            raise ValueError('Expected a matrix for %s' % name)
        if value.shape != _GRADABD_SHAPES[name]:
            raise ValueError('Expected shape %s for %s, got %s'
                             % (_GRADABD_SHAPES[name], name, value.shape))
        setattr(grad, name, np.ascontiguousarray(value))
    return grad


def _lamina_alone_to_data(ply):
    mat = ply.matlamina
    return _lamina_to_data(ply, None if mat is None
                           else _matlamina_to_data(mat))


def _lamina_alone_from_data(data):
    if not isinstance(data, dict):
        raise ValueError('Expected a JSON object for Lamina, got %r'
                         % type(data).__name__)
    mat = data.get('matlamina')
    return _lamina_from_data(data, None if mat is None
                             else _matlamina_from_data(mat))


# NOTE Laminate before Lamina and the others, since isinstance is used
_TYPES = (
    ('Laminate', Laminate, _laminate_to_data, _laminate_from_data),
    ('Lamina', Lamina, _lamina_alone_to_data, _lamina_alone_from_data),
    ('MatLamina', MatLamina, _matlamina_to_data, _matlamina_from_data),
    ('LaminationParameters', LaminationParameters, _lp_to_data,
     _lp_from_data),
    ('GradABD', GradABD, _gradabd_to_data, _gradabd_from_data),
)


def to_dict(obj):
    r"""Convert an object into a dictionary of JSON-compatible types

    Parameters
    ----------
    obj : :class:`.Laminate`, :class:`.Lamina`, :class:`.MatLamina`, :class:`.LaminationParameters` or :class:`.GradABD`
        The object to be converted.

    Returns
    -------
    d : dict
        Dictionary with the keys ``'type'``, ``'format_version'``,
        ``'composites_version'`` and ``'data'``, containing only ``dict``,
        ``list``, ``str``, ``int``, ``float`` and ``None``.

    """
    for typename, cls, encode, _ in _TYPES:
        if isinstance(obj, cls):
            return {'type': typename,
                    'format_version': FORMAT_VERSION,
                    'composites_version': __version__,
                    'data': encode(obj)}
    raise TypeError('Cannot convert an object of type %r to JSON'
                    % type(obj).__name__)


def from_dict(d):
    r"""Create an object from a dictionary created by :func:`.to_dict`

    Parameters
    ----------
    d : dict
        Dictionary as returned by :func:`.to_dict`.

    Returns
    -------
    obj : :class:`.Laminate`, :class:`.Lamina`, :class:`.MatLamina`, :class:`.LaminationParameters` or :class:`.GradABD`
        The object. Attributes missing in ``d`` keep their default values.

    Raises
    ------
    ValueError
        If ``d`` is not valid, has unknown keys, or was created by a newer
        format version.

    """
    if not isinstance(d, dict):
        raise ValueError('Expected a JSON object, got %r' % type(d).__name__)
    _check_keys(d, ('type', 'format_version', 'composites_version', 'data'),
                'the saved object')
    format_version = d.get('format_version')
    if (not isinstance(format_version, numbers.Integral)
            or isinstance(format_version, bool) or format_version < 1):
        raise ValueError('Invalid format_version: %r' % (format_version, ))
    if format_version > FORMAT_VERSION:
        raise ValueError('format_version %d was saved by a newer version of '
                         'composites (%s), this version reads up to %d'
                         % (format_version, d.get('composites_version'),
                            FORMAT_VERSION))
    typename = d.get('type')
    for name, _, _, decode in _TYPES:
        if typename == name:
            return decode(d.get('data', {}))
    raise ValueError('Unknown type: %r' % (typename, ))


def to_json(obj, **kwargs):
    r"""Convert an object into a JSON string

    Parameters
    ----------
    obj : :class:`.Laminate`, :class:`.Lamina`, :class:`.MatLamina`, :class:`.LaminationParameters` or :class:`.GradABD`
        The object to be converted.
    kwargs : dict, optional
        Passed to :func:`json.dumps`, e.g. ``indent=2``.

    Returns
    -------
    s : str
        The JSON string.

    """
    return json.dumps(to_dict(obj), allow_nan=False, **kwargs)


def from_json(s):
    r"""Create an object from a JSON string created by :func:`.to_json`

    Parameters
    ----------
    s : str or bytes
        The JSON string.

    Returns
    -------
    obj : :class:`.Laminate`, :class:`.Lamina`, :class:`.MatLamina`, :class:`.LaminationParameters` or :class:`.GradABD`
        The object.

    """
    return from_dict(json.loads(s))


def save_json(obj, fname, indent=2):
    r"""Save an object to a JSON file

    Parameters
    ----------
    obj : :class:`.Laminate`, :class:`.Lamina`, :class:`.MatLamina`, :class:`.LaminationParameters` or :class:`.GradABD`
        The object to be saved.
    fname : str, path-like or file object
        Name of the file, or a file object opened in text mode for writing.
    indent : int or None, optional
        Indentation of the JSON output, ``None`` for the most compact form.

    """
    s = to_json(obj, indent=indent)
    if hasattr(fname, 'write'):
        fname.write(s)
    else:
        with open(os.fspath(fname), 'w', encoding='utf-8') as f:
            f.write(s)


def load_json(fname):
    r"""Load an object from a JSON file created by :func:`.save_json`

    Parameters
    ----------
    fname : str, path-like or file object
        Name of the file, or a file object opened for reading.

    Returns
    -------
    obj : :class:`.Laminate`, :class:`.Lamina`, :class:`.MatLamina`, :class:`.LaminationParameters` or :class:`.GradABD`
        The object.

    """
    if hasattr(fname, 'read'):
        return from_json(fname.read())
    with open(os.fspath(fname), 'r', encoding='utf-8') as f:
        return from_json(f.read())
