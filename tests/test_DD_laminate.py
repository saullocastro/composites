import sys
sys.path.append('..')

import numpy as np

from composites.utils import (read_laminaprop, double_double_laminate,
                              laminated_plate)

# Melo, J. D. D., Bi, J. and Tsai, S. W., "A novel invariant-based design
# approach to carbon fiber reinforced laminates", Composite Structures, 159,
# 44-52, 2017, https://doi.org/10.1016/j.compstruct.2016.09.055, Table 1,
# p. 45: Ex, Ey, nux, Es and the trace Tr of Eq. (1), p. 45, of the
# plane-stress stiffness matrix, in GPa, printed as integers
data = {
'IM6/epoxy': dict(Ex=203e9, Ey=11.20e9, vx=0.32, Es=8.40e9, tr=232e9),
'IM7/977-3': dict(Ex=191e9, Ey=9.94e9, vx=0.35, Es=7.79e9, tr=218e9),
'T4708/MR60H': dict(Ex=142e9, Ey=7.72e9, vx=0.34, Es=3.80e9, tr=158e9),
}

def test_trace_normalized():
    r"""
    Reference:

        Melo, J. D. D., Bi, J., and Tsai, S. W., 2017, “A Novel Invariant-Based
        Design Approach to Carbon Fiber Reinforced Laminates,” Compos. Struct.,
        159, pp. 44–52. https://doi.org/10.1016/j.compstruct.2016.09.055

    The trace is ``tr = Q11 + Q22 + 2 Q66``, Eq. (1), p. 45, of the
    plane-stress stiffnesses, with the materials of Table 1, p. 45, see
    ``data``. The invariants are written with ``Q66 = G12`` and the
    transverse shear ones with ``Q44 = G23`` and ``Q55 = G13``.

    """
    for material, d in data.items():
        m = read_laminaprop((d['Ex'], d['Ey'], d['vx'], d['Es'], d['Es'], d['Es']))
        tr = m.q11 + m.q22 + 2*m.q66
        # NOTE Tr of Melo et al. (2017), Table 1, p. 45, is printed as an
        #      integer in GPa: the traces of Eq. (1) computed from Ex, Ey,
        #      nux and Es differ from it by +0.09 %, -0.09 % and +0.17 %,
        #      within half a unit of the last printed digit
        assert np.isclose(tr, d['tr'], atol=0.5e9, rtol=0)
        q11 = m.q11
        q12 = m.q12
        q22 = m.q22
        q44 = m.q44
        q55 = m.q55
        q66 = m.q66
        q11 /= tr
        q12 /= tr
        q22 /= tr
        q44 /= tr
        q55 /= tr
        q66 /= tr
        u1 = (3*q11 + 3*q22 + 2*q12 + 4*q66) / 8.
        u2 = (q11 - q22) / 2.
        u3 = (q11 + q22 - 2*q12 - 4*q66) / 8.
        u4 = (q11 + q22 + 6*q12 - 4*q66) / 8.
        u5 = (u1 - u4) / 2.
        u6 = (q44 + q55) / 2.
        u7 = (q44 - q55) / 2.
        tr_norm_inv = (u1, u2, u3, u4, u5, u6, u7)
        m.trace_normalize_plane_stress()
        tr_norm_inv2 = (m.u1, m.u2, m.u3, m.u4, m.u5, m.u6, m.u7)
        assert np.allclose(tr_norm_inv, tr_norm_inv2)

def test_ABD():
    d = data['IM6/epoxy']
    laminaprop = (d['Ex'], d['Ey'], d['vx'], d['Es'], d['Es'], d['Es'])
    thickness = 2.3
    phideg = 15
    psideg = 30
    lam = double_double_laminate(thickness, phideg, psideg, laminaprop)

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
    test_trace_normalized()
    test_ABD()
