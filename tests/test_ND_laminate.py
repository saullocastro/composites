r"""N-double laminate of :mod:`composites.utils`

Verifies that the ``A11``, ``A12``, ``A22``, ``A66`` and ``D11``, ``D12``,
``D22``, ``D66`` of :func:`composites.utils.n_double_laminate` with the
angles ``[15, 30]`` are those of :func:`.laminated_plate` with the stack
``[15/-15/30/-30]`` repeated 32 times, for the IM6/epoxy material of
``data``, the material data of ``tests/test_DD_laminate.py``.

Notes
-----
The stiffnesses of the N-double laminate are those of the homogenized
limit of many plies, see the Notes of :func:`composites.utils.n_double_laminate`;
the ``A`` and ``D`` of the stack of 128 plies agree with them within the
default tolerance of :func:`numpy.allclose`.

"""
import sys
sys.path.append('..')

import numpy as np

from composites.utils import (read_laminaprop, n_double_laminate,
                              laminated_plate)

data = {
'IM6/epoxy': dict(Ex=203e9, Ey=11.20e9, vx=0.32, Es=8.40e9, tr=232e9),
'IM7/977-3': dict(Ex=191e9, Ey=9.94e9, vx=0.35, Es=7.79e9, tr=218e9),
'T4708/MR60H': dict(Ex=142e9, Ey=7.72e9, vx=0.34, Es=3.80e9, tr=158e9),
}

def test_ABD():
    d = data['IM6/epoxy']
    laminaprop = (d['Ex'], d['Ey'], d['vx'], d['Es'], d['Es'], d['Es'])
    thickness = 2.3
    phideg = 15
    psideg = 30
    lam = n_double_laminate(thickness, [phideg, psideg], laminaprop)

    num_repeats = 32
    stack = [phideg, -phideg, psideg, -psideg]*num_repeats
    plyt = thickness/len(stack)
    lam_ref = laminated_plate(stack, plyt, laminaprop)

    assert np.allclose(lam.A11, lam_ref.A11)
    assert np.allclose(lam.A12, lam_ref.A12)
    assert np.allclose(lam.A22, lam_ref.A22)
    assert np.allclose(lam.A66, lam_ref.A66)
    assert np.allclose(lam.D11, lam_ref.D11)
    assert np.allclose(lam.D12, lam_ref.D12)
    assert np.allclose(lam.D22, lam_ref.D22)
    assert np.allclose(lam.D66, lam_ref.D66)

if __name__ == '__main__':
    test_ABD()
