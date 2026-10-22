import sys
sys.path.append('..')

import numpy as np
import pytest

from composites.utils import (read_laminaprop, laminated_plate,
        isotropic_plate)
from composites.core import (laminate_from_LaminationParameters,
                             laminate_from_lamination_parameters,
                             make_balanced_LP, make_orthotropic_LP,
                             make_symmetric_LP, Lamina, GradABD,
                             LaminationParameters)


def test_lampar_tri_axial():
    E = 71e9
    nu = 0.33
    G = E/(2*(1+nu))
    lamprop = (E, E, nu, G, G, G, E, nu, nu)
    rho = 0
    thickness = 1
    matlamina = read_laminaprop(lamprop, rho)
    matlamina.get_constitutive_matrix()
    matlamina.get_invariant_matrix()
    ply = Lamina()
    ply.thetadeg = 45.
    ply.h = 3.
    ply.matlamina = matlamina
    ply.get_transf_matrix_displ_to_laminate()
    ply.get_constitutive_matrix()
    ply.get_transf_matrix_stress_to_lamina()
    ply.get_transf_matrix_stress_to_laminate()

    lam = laminate_from_lamination_parameters(thickness, matlamina,
        0.5, 0.4, -0.3, -0.6,
        0.5, 0.4, -0.3, -0.6,
        0.5, 0.4, -0.3, -0.6,
        0.5, 0.4)
    A = np.array([[1.05196816e+11, 5.18133569e+10, 0.00000000e+00],
                  [5.18133569e+10, 1.05196816e+11, 0.00000000e+00],
                  [0.00000000e+00, 0.00000000e+00, 2.66917293e+10]])
    B = np.array([[0, 0, 0],
                  [0, 0, 0],
                  [0, 0, 0]])
    D = np.array([[8.76640130e+09, 4.31777974e+09, 0.00000000e+00],
                  [4.31777974e+09, 8.76640130e+09, 0.00000000e+00],
                  [0.00000000e+00, 0.00000000e+00, 2.22431078e+09]])
    Ats = np.array([[2.66917293e+10, 0.00000000e+00],
                    [0.00000000e+00, 2.66917293e+10]])
    assert np.allclose(lam.A, A)
    lam.make_symmetric()
    assert np.allclose(lam.B, B)
    assert np.allclose(lam.D, D)
    assert np.allclose(lam.Ats, Ats)
    assert np.allclose(lam.Abar_ts, Ats)
    ABD = lam.ABD
    assert np.allclose(ABD[:3, :3], A)
    assert np.allclose(ABD[3:, 3:], D)


def test_lampar_plane_stress():
    E = 71e9
    nu = 0.33
    lamprop = (E, nu)
    rho = 0
    thickness = 1
    matlamina = read_laminaprop(lamprop, rho)
    matlamina.get_constitutive_matrix()
    matlamina.get_invariant_matrix()
    ply = Lamina()
    ply.thetadeg = 45.
    ply.h = 3.
    ply.matlamina = matlamina
    ply.get_transf_matrix_displ_to_laminate()
    ply.get_constitutive_matrix()
    ply.get_transf_matrix_stress_to_lamina()
    ply.get_transf_matrix_stress_to_laminate()

    lam = laminate_from_lamination_parameters(thickness, matlamina,
        0.5, 0.4, -0.3, -0.6,
        0.5, 0.4, -0.3, -0.6,
        0.5, 0.4, -0.3, -0.6,
        0.5, 0.4)
    A = np.array([[7.96768040e+10, 2.62933453e+10, 1.14440918e-06],
                  [2.62933453e+10, 7.96768040e+10, -1.14440918e-06],
                  [1.14440918e-06, -1.14440918e-06, 2.66917293e+10]])
    B = np.array([[0., 0., 0.],
                  [0., 0., 0.],
                  [0., 0., 0.]])
    D = np.array([[6.63973366e+09, 2.19111211e+09, 9.53674316e-08],
                  [2.19111211e+09, 6.63973366e+09, -9.53674316e-08],
                  [9.53674316e-08, -9.53674316e-08, 2.22431078e+09]])
    Ats = np.array([[2.66917293e+10, 0.00000000e+00],
                    [0.00000000e+00, 2.66917293e+10]])
    assert np.allclose(lam.A, A)
    assert np.allclose(lam.E, lam.E.T)
    assert np.allclose(lam.F, lam.F.T)
    assert np.allclose(lam.H, lam.H.T)
    lam.make_symmetric()
    assert np.allclose(lam.B, B)
    assert np.allclose(lam.D, D)
    assert np.allclose(lam.Ats, Ats)
    assert np.allclose(lam.Abar_ts, Ats)
    ABD = lam.ABD
    assert np.allclose(ABD[:3, :3], A)
    assert np.allclose(ABD[3:, 3:], D)


def test_laminated_plate_tri_axial():
    lamprop = (71e9, 7e9, 0.28, 7e9, 7e9, 7e9, 7e9, 0.28, 0.28)
    stack = [0, 45, 90]
    plyt = 0.000125
    lam = laminated_plate(stack, plyt, lamprop)
    A = np.array([[ 13280892.30559593, 2198758.85719477, 2015579.57848837],
                  [  2198758.85719477,13280892.30559593, 2015579.57848837],
                  [  2015579.57848837, 2015579.57848837, 4083033.36210029]])
    B = np.array([[ -1.00778979e+03, 0.00000000e+00, 8.53496487e-15],
                  [  0.00000000e+00, 1.00778979e+03, 5.31743621e-14],
                  [  8.53496487e-15, 5.31743621e-14, 0.00000000e+00]])
    D = np.array([[ 0.1708233 , 0.01057886, 0.00262445],
                  [ 0.01057886, 0.1708233 , 0.00262445],
                  [ 0.00262445, 0.00262445, 0.0326602 ]])
    Abar_ts = np.array([[ 2625000.,       0.],
                        [       0., 2625000.]])
    Ats = np.array([[2014905.27266639, -27210.29377157],
                    [ -27210.29377157, 2014905.27266639]])
    assert np.allclose(lam.A, A)
    assert np.allclose(lam.B, B)
    assert np.allclose(lam.D, D)
    assert np.allclose(lam.Abar_ts, Abar_ts)
    assert np.allclose(lam.Ats, Ats)


def test_laminated_plate_plane_stress():
    lamprop = (71e9, 7e9, 0.28, 7e9, 7e9, 7e9)
    stack = [0, 45, 90]
    plyt = 0.000125
    lam = laminated_plate(stack, plyt, lamprop)
    A = np.array([[ 13280892.30559593, 2198758.85719477, 2015579.57848837],
                  [  2198758.85719477,13280892.30559593, 2015579.57848837],
                  [  2015579.57848837, 2015579.57848837, 4083033.36210029]])
    B = np.array([[ -1.00778979e+03, 0.00000000e+00, 8.53496487e-15],
                  [  0.00000000e+00, 1.00778979e+03, 5.31743621e-14],
                  [  8.53496487e-15, 5.31743621e-14, 0.00000000e+00]])
    D = np.array([[ 0.1708233 , 0.01057886, 0.00262445],
                  [ 0.01057886, 0.1708233 , 0.00262445],
                  [ 0.00262445, 0.00262445, 0.0326602 ]])
    E = np.array([[-1.96833943e-05, 0.00000000e+00, 1.66698533e-22],
                  [0.00000000e+00, 1.96833943e-05, 1.03856176e-21],
                  [1.66698533e-22, 1.03856176e-21, 0.00000000e+00]])
    F = np.array([[3.63890059e-09, 1.87551265e-10, 6.15106073e-12],
                  [1.87551265e-10, 3.63890059e-09, 6.15106073e-12],
                  [6.15106073e-12, 6.15106073e-12, 6.53329570e-10]])
    H = np.array([[9.14779627e-17, 4.61039304e-18, 1.71625578e-20],
                  [4.61039304e-18, 9.14779627e-17, 1.71625578e-20],
                  [1.71625578e-20, 1.71625578e-20, 1.63068348e-17]])
    Abar_ts = np.array([[ 2625000.,       0.],
                        [       0., 2625000.]])
    Ats = np.array([[2014905.27266639, -27210.29377157],
                    [ -27210.29377157, 2014905.27266639]])
    Dts = np.array([[0.03076172, 0.],
                    [0., 0.03076172]])
    Fts = np.array([[6.48880005e-10, 0.00000000e+00],
                    [0.00000000e+00, 6.48880005e-10]])

    assert np.allclose(lam.A, A)
    assert np.allclose(lam.B, B)
    assert np.allclose(lam.D, D)
    assert np.allclose(lam.E, E)
    assert np.allclose(lam.F, F)
    assert np.allclose(lam.H, H)
    assert np.allclose(lam.Abar_ts, Abar_ts)
    assert np.allclose(lam.Ats, Ats)
    assert np.allclose(lam.Dts, Dts)
    assert np.allclose(lam.Fts, Fts)
    with pytest.warns(DeprecationWarning):
        lam.calc_scf()
    lam.calc_equivalent_properties()
    lp = lam.calc_lamination_parameters()
    matlamina = lam.plies[0].matlamina
    thickness = lam.h
    lam_2 = laminate_from_LaminationParameters(thickness, matlamina, lp)
    assert np.allclose(lam_2.A, lam.A)
    assert np.allclose(lam_2.B, lam.B)
    assert np.allclose(lam_2.D, lam.D)
    assert np.allclose(lam_2.Ats, lam.Abar_ts)

    lam.make_balanced()
    make_balanced_LP(lp)
    lam_2 = laminate_from_LaminationParameters(thickness, matlamina, lp)
    assert np.allclose(lam_2.A, lam.A)
    assert np.allclose(lam_2.B, lam.B)
    assert np.allclose(lam_2.D, lam.D)
    assert np.allclose(lam_2.Ats, lam.Abar_ts)

    lam.make_orthotropic()
    make_orthotropic_LP(lp)
    lam_2 = laminate_from_LaminationParameters(thickness, matlamina, lp)
    assert np.allclose(lam_2.A, lam.A)
    assert np.allclose(lam_2.B, lam.B)
    assert np.allclose(lam_2.D, lam.D)
    assert np.allclose(lam_2.Ats, lam.Abar_ts)

    lam = laminated_plate(stack, plyt, lamprop)
    lp = lam.calc_lamination_parameters()
    lam.make_symmetric()
    make_symmetric_LP(lp)
    lam_2 = laminate_from_LaminationParameters(thickness, matlamina, lp)
    assert np.allclose(lam_2.A, lam.A)
    assert np.allclose(lam_2.B, lam.B)
    assert np.allclose(lam_2.D, lam.D)
    assert np.allclose(lam_2.Ats, lam.Abar_ts)

    lam = laminated_plate(stack, plyt, lamprop)
    lp = lam.calc_lamination_parameters()
    lam.make_smeared()
    assert np.allclose(lam.D, lam.h**2/12*lam.A)
    assert np.allclose(lam.B, 0)


def test_isotropic_plate():
    E = 71e9
    nu = 0.28
    thick = 0.000125
    lam = isotropic_plate(thickness=thick, E=E, nu=nu)
    A = np.array([[9629991.31944444, 2696397.56944444,       0.   ],
                  [2696397.56944444, 9629991.31944444,       0.   ],
                  [      0.        ,       0.        , 3466796.875]])
    D = np.array([[0.01253905, 0.00351093, 0.        ],
                  [0.00351093, 0.01253905, 0.        ],
                  [0.        , 0.        , 0.00451406]])
    Abar_ts = np.array([[3466796.875,       0.   ],
                        [      0.   , 3466796.875]])
    assert np.allclose(lam.A, A)
    assert np.allclose(lam.B, 0)

    assert np.allclose(lam.D, D)
    assert np.allclose(lam.Abar_ts, Abar_ts)
    assert np.allclose(lam.Ats, 5/6*Abar_ts)



def test_errors():
    E = 71e9
    nu = 0.28
    thick = 0.000125
    lam = isotropic_plate(thickness=thick, E=E, nu=nu)
    lam.offset = 1.
    try:
        lam.make_balanced()
    except RuntimeError:
        pass
    try:
        lam.make_orthotropic()
    except RuntimeError:
        pass
    try:
        lam.make_symmetric()
    except RuntimeError:
        pass
    try:
        lam.make_smeared()
    except RuntimeError:
        pass
    try:
        lam.plies = []
        lam.calc_lamination_parameters()
    except ValueError:
        pass


def test_laminate_LP_gradients():
    E = 71e9
    nu = 0.33
    lamprop = (E, nu)
    rho = 0
    thickness = 1
    matlamina = read_laminaprop(lamprop, rho)
    lp = LaminationParameters()
    lp.xiA1 = 0.5
    lp.xiA2 = 0.4
    lp.xiA3 = -0.3
    lp.xiA4 = -0.6
    lp.xiB1 = 0.5
    lp.xiB2 = 0.4
    lp.xiB3 = -0.3
    lp.xiB4 = -0.6
    lp.xiD1 = 0.5
    lp.xiD2 = 0.4
    lp.xiD3 = -0.3
    lp.xiD4 = -0.6
    gradABD = GradABD()
    gradABD.calc_LP_grad(thickness, matlamina, lp)
    print(gradABD.gradAij)


def test_laminated_plate_missing_arguments():
    laminaprop = (71e3, 71e3, 0.33)
    with pytest.raises(ValueError, match='plyt or plyts'):
        laminated_plate([0, 90], laminaprop=laminaprop)
    with pytest.raises(ValueError, match='laminaprop or laminaprops'):
        laminated_plate([0, 90], plyt=0.125)


def test_deprecated_xiAtrans():
    lp = LaminationParameters()
    with pytest.warns(DeprecationWarning, match='xiAts1'):
        lp.xiAtrans1 = 0.3
    with pytest.warns(DeprecationWarning, match='xiAts2'):
        lp.xiAtrans2 = -0.2
    assert (lp.xiAts1, lp.xiAts2) == (0.3, -0.2)
    with pytest.warns(DeprecationWarning, match='xiAts1'):
        assert lp.xiAtrans1 == 0.3
    with pytest.warns(DeprecationWarning, match='xiAts2'):
        assert lp.xiAtrans2 == -0.2

    matlamina = read_laminaprop((142e9, 8.7e9, 0.28, 5.1e9, 5.1e9, 3.2e9))
    args = (1e-3, matlamina) + (0.,)*12
    ref = laminate_from_lamination_parameters(*args, xiAts1=0.3, xiAts2=-0.2)
    with pytest.warns(DeprecationWarning, match='xiAts'):
        lam = laminate_from_lamination_parameters(*args, xiAtrans1=0.3,
                                                  xiAtrans2=-0.2)
    assert np.array_equal(lam.Abar_ts, ref.Abar_ts)


def test_deprecated_gradAtransij():
    lam = laminated_plate([0, 45, 90], plyt=0.125e-3,
            laminaprop=(142e9, 8.7e9, 0.28, 5.1e9, 5.1e9, 3.2e9))
    grad = GradABD()
    grad.calc_LP_grad(lam.h, lam.plies[0].matlamina,
                      lam.calc_lamination_parameters())
    with pytest.warns(DeprecationWarning, match='gradAtsij'):
        gradAtransij = grad.gradAtransij
    assert np.array_equal(gradAtransij, grad.gradAtsij)
    with pytest.warns(DeprecationWarning, match='gradAtsij'):
        grad.gradAtransij = np.ones((3, 3))
    assert np.array_equal(grad.gradAtsij, np.ones((3, 3)))


def test_deprecated_Dtrans_Ftrans():
    lam = laminated_plate([0, 45, 90], plyt=0.125e-3,
            laminaprop=(142e9, 8.7e9, 0.28, 5.1e9, 5.1e9, 3.2e9))
    with pytest.warns(DeprecationWarning, match='Dts'):
        Dtrans = lam.Dtrans
    with pytest.warns(DeprecationWarning, match='Fts'):
        Ftrans = lam.Ftrans
    assert np.array_equal(Dtrans, lam.Dts)
    assert np.array_equal(Ftrans, lam.Fts)


def _abd_lp(thickness, matlamina, xiA, xiB, xiD, xiAts):
    lam = laminate_from_lamination_parameters(thickness, matlamina, *xiA,
            *xiB, *xiD, xiAts1=xiAts[0], xiAts2=xiAts[1])
    A = [lam.A11, lam.A12, lam.A16, lam.A22, lam.A26, lam.A66]
    B = [lam.B11, lam.B12, lam.B16, lam.B22, lam.B26, lam.B66]
    D = [lam.D11, lam.D12, lam.D16, lam.D22, lam.D26, lam.D66]
    Ats = [lam.Abar44, lam.Abar45, lam.Abar55]
    return np.array(A), np.array(B), np.array(D), np.array(Ats)


def test_laminate_LP_gradients_finite_differences():
    # NOTE orthotropic material, since u2 = u3 = u7 = 0 for an isotropic one
    matlamina = read_laminaprop((142e9, 8.7e9, 0.28, 5.1e9, 5.1e9, 3.2e9))
    thickness = 1.3e-3
    xi = {'A': [0.5, 0.4, -0.3, -0.6],
          'B': [0.2, -0.1, 0.35, 0.15],
          'D': [-0.4, 0.3, 0.25, -0.2],
          'Ats': [0.3, -0.2]}
    lp = LaminationParameters()
    for key, values in xi.items():
        for i, value in enumerate(values):
            setattr(lp, 'xi%s%d' % (key, i + 1), value)
    grad = GradABD()
    grad.calc_LP_grad(thickness, matlamina, lp)
    grads = {'A': np.asarray(grad.gradAij), 'B': np.asarray(grad.gradBij),
             'D': np.asarray(grad.gradDij), 'Ats': np.asarray(grad.gradAtsij)}
    keys = ['A', 'B', 'D', 'Ats']

    def stiffnesses(h, xi):
        return dict(zip(keys, _abd_lp(h, matlamina, xi['A'], xi['B'],
                                      xi['D'], xi['Ats'])))

    # column 0: central differences with respect to the thickness
    dh = 1e-6*thickness
    plus = stiffnesses(thickness + dh, xi)
    minus = stiffnesses(thickness - dh, xi)
    for key in keys:
        fd = (plus[key] - minus[key])/(2*dh)
        scale = np.abs(fd).max()
        assert np.allclose(grads[key][:, 0], fd, rtol=1e-6, atol=1e-6*scale)

    # other columns: central differences with respect to each lamination
    # parameter (the stiffnesses are linear in them)
    dxi = 1e-4
    for key in keys:
        for j in range(len(xi[key])):
            xi_plus = {k: list(v) for k, v in xi.items()}
            xi_minus = {k: list(v) for k, v in xi.items()}
            xi_plus[key][j] += dxi
            xi_minus[key][j] -= dxi
            fd = (stiffnesses(thickness, xi_plus)[key]
                - stiffnesses(thickness, xi_minus)[key])/(2*dxi)
            scale = np.abs(grads[key]).max()
            assert np.allclose(grads[key][:, j + 1], fd, rtol=1e-6,
                               atol=1e-9*scale), (key, j)


def test_matlamina_constitutive_matrix_inverse_compliance():
    # NOTE e2 != e3 and nu23 != 0, otherwise nu23*nu23 = nu23*nu32
    e1, e2, e3 = 142e9, 8.7e9, 10.1e9
    nu12, nu13, nu23 = 0.28, 0.31, 0.45
    g12, g13, g23 = 5.1e9, 5.3e9, 3.2e9
    matlamina = read_laminaprop((e1, e2, nu12, g12, g13, g23, e3, nu13,
                                 nu23))
    # compliance matrix, rows and columns ordered as 11, 22, 33, 23, 13, 12,
    # with nu21/e2 = nu12/e1, nu31/e3 = nu13/e1 and nu32/e3 = nu23/e2
    S = np.array([[1/e1, -nu12/e1, -nu13/e1, 0, 0, 0],
                  [-nu12/e1, 1/e2, -nu23/e2, 0, 0, 0],
                  [-nu13/e1, -nu23/e2, 1/e3, 0, 0, 0],
                  [0, 0, 0, 1/g23, 0, 0],
                  [0, 0, 0, 0, 1/g13, 0],
                  [0, 0, 0, 0, 0, 1/g12]])
    C_ref = np.linalg.inv(S)
    C = np.asarray(matlamina.get_constitutive_matrix())
    assert np.allclose(C, C_ref, rtol=1e-12, atol=1e-12*C_ref.max())
    # the 3D constitutive constants q_ij are the same inverse
    Q = np.array([[matlamina.q11, matlamina.q12, matlamina.q13],
                  [matlamina.q21, matlamina.q22, matlamina.q23],
                  [matlamina.q31, matlamina.q32, matlamina.q33]])
    assert np.allclose(Q, C_ref[:3, :3], rtol=1e-12, atol=0)
    assert np.allclose([matlamina.q44, matlamina.q55, matlamina.q66],
                       np.diag(C_ref)[3:], rtol=1e-12, atol=0)


def test_density_moments():
    # unsymmetric laminate, ply densities and an offset, against an
    # independent Gauss-Legendre integration of rho*z**n ply by ply
    plyts = [0.2e-3, 0.5e-3, 0.3e-3, 0.4e-3]
    rhos = [1600., 1200., 2700., 900.]
    offset = 0.37e-3
    lam = laminated_plate([0, 45, 90, -30], plyts=plyts,
            laminaprop=(142e9, 8.7e9, 0.28, 5.1e9, 5.1e9, 3.2e9), rhos=rhos,
            offset=offset)
    xi, w = np.polynomial.legendre.leggauss(8)
    zbot = -sum(plyts)/2 + offset
    ref = dict.fromkeys((0, 1, 2, 3, 4, 6), 0.)
    for plyt, rho in zip(plyts, rhos):
        z = zbot + plyt*(xi + 1)/2
        for n in ref:
            ref[n] += rho*np.sum(w*z**n)*plyt/2
        zbot += plyt
    names = {0: 'intrho', 1: 'intrhoz', 2: 'intrhoz2', 3: 'intrhoz3',
             4: 'intrhoz4', 6: 'intrhoz6'}
    for n, name in names.items():
        assert np.isclose(getattr(lam, name), ref[n], rtol=1e-12, atol=0), name
    # recomputed from zero by every call
    lam.calc_constitutive_matrix()
    for n, name in names.items():
        assert np.isclose(getattr(lam, name), ref[n], rtol=1e-12, atol=0), name
    # not computed from lamination parameters, see the docstring
    lam_lp = laminate_from_LaminationParameters(sum(plyts),
            lam.plies[0].matlamina, lam.calc_lamination_parameters())
    for name in names.values():
        assert getattr(lam_lp, name) == 0.


def test_example_laminated_plate():
    # quasi-isotropic CFRP laminate [0/45/-45/90]s, SI units
    laminaprop = (142e9, 8.7e9, 0.28, 5.1e9, 5.1e9, 3.2e9) # E11, E22, nu12, G12, G13, G23
    stack = [0, 45, -45, 90, 90, -45, 45, 0] # from the bottom to the top
    lam = laminated_plate(stack, plyt=0.125e-3, laminaprop=laminaprop,
                          rho=1600.)
    # ABD relates {Nxx, Nyy, Nxy, Mxx, Myy, Mxy} to the generalized strains
    # {exx0, eyy0, gxy0, kxx, kyy, kxy}
    ABD = lam.ABD
    assert ABD.shape == (6, 6)
    assert np.allclose(ABD[:3, :3], lam.A) and np.allclose(ABD[:3, 3:], lam.B)
    assert np.allclose(ABD[3:, 3:], lam.D)
    # symmetric laminate: no membrane-bending coupling
    assert np.allclose(lam.B, 0, atol=1e-12*lam.A11*lam.h)
    # quasi-isotropic membrane stiffness
    assert np.isclose(lam.A11, lam.A22, rtol=1e-12)
    assert np.isclose(lam.A66, (lam.A11 - lam.A12)/2, rtol=1e-12)
    assert np.allclose([lam.A16, lam.A26], 0, atol=1e-12*lam.A11)
    # but not quasi-isotropic in bending: the +-45 plies give D16, D26
    assert abs(lam.D16) > 1e-3*lam.D11
    # equivalent membrane properties
    assert np.isclose(lam.e1, lam.e2, rtol=1e-12)
    assert np.isclose(lam.g12, lam.e1/(2*(1 + lam.nu12)), rtol=1e-12)
    # FSDT: transverse shear stiffness with the shear correction applied
    assert lam.shear_correction == 'rohwer'
    assert np.allclose(lam.Ats, lam.Ats.T)
    assert np.isclose(lam.scf_k13, lam.A55/lam.Abar55, rtol=1e-15)
    assert np.isclose(lam.scf_k23, lam.A44/lam.Abar44, rtol=1e-15)
    # mass per unit area
    assert np.isclose(lam.intrho, 1600.*lam.h, rtol=1e-12)


def test_example_isotropic_plate_offset():
    # closed form of an isotropic plate whose mid-surface is at z = d
    E, nu, h, d, rho = 70e9, 0.3, 2e-3, 0.6e-3, 2700.
    lam = isotropic_plate(thickness=h, E=E, nu=nu, offset=d, rho=rho)
    C = E/(1 - nu**2)*np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu)/2]])
    G = E/(2*(1 + nu))
    atol = 1e-12*E
    assert np.allclose(lam.A, h*C, rtol=1e-12, atol=atol*h)
    assert np.allclose(lam.B, d*h*C, rtol=1e-12, atol=atol*h*d)
    assert np.allclose(lam.D, h*(d**2 + h**2/12)*C, rtol=1e-12, atol=atol*h**3)
    # the shear correction of an isotropic plate is 5/6
    assert np.allclose(lam.Ats, 5/6*G*h*np.eye(2), rtol=1e-12, atol=atol*h)
    # with s uniform in [-h/2, h/2]: <s^2> = h^2/12, <s^4> = h^4/80,
    # <s^6> = h^6/448, and int (d + s)^n ds = h <(d + s)^n>
    m = [1., d, d**2 + h**2/12, d**3 + d*h**2/4,
         d**4 + d**2*h**2/2 + h**4/80, None,
         d**6 + 15/12*d**4*h**2 + 15/80*d**2*h**4 + h**6/448]
    assert np.allclose(lam.E, h*m[3]*C, rtol=1e-12, atol=atol*h**4)
    assert np.allclose(lam.F, h*m[4]*C, rtol=1e-12, atol=atol*h**5)
    assert np.allclose(lam.H, h*m[6]*C, rtol=1e-12, atol=atol*h**7)
    for name, n in (('intrho', 0), ('intrhoz', 1), ('intrhoz2', 2),
                    ('intrhoz3', 3), ('intrhoz4', 4), ('intrhoz6', 6)):
        assert np.isclose(getattr(lam, name), rho*h*m[n], rtol=1e-12, atol=0)
    # without offset the odd moments vanish
    lam0 = isotropic_plate(thickness=h, E=E, nu=nu, rho=rho)
    assert abs(lam0.intrhoz) < 1e-15*rho*h**2
    assert abs(lam0.intrhoz3) < 1e-15*rho*h**4
    assert np.allclose(lam0.B, 0, atol=atol*h**2)


def test_example_tsdt_matrices():
    # unsymmetric laminate with an offset
    laminaprop = (142e9, 8.7e9, 0.28, 5.1e9, 5.1e9, 3.2e9)
    lam = laminated_plate([0, 45, 90, -30], plyts=[0.2e-3, 0.5e-3, 0.3e-3,
        0.4e-3], laminaprop=laminaprop, offset=0.37e-3)
    # 13 x 13 constitutive matrix of the TSDT, the generalized strains being
    # {eps0, eps1, eps3, gamma0, gamma2}, with eps = eps0 + z eps1 + z^3 eps3
    # and gamma = {gamma_yz, gamma_xz} = gamma0 + z^2 gamma2
    Z3 = np.zeros((3, 3))
    Z32 = np.zeros((3, 2))
    K = np.block([[lam.A, lam.B, lam.E, Z32, Z32],
                  [lam.B, lam.D, lam.F, Z32, Z32],
                  [lam.E, lam.F, lam.H, Z32, Z32],
                  [Z32.T, Z32.T, Z32.T, lam.Abar_ts, lam.Dts],
                  [Z32.T, Z32.T, Z32.T, lam.Dts, lam.Fts]])
    assert np.allclose(K, K.T, rtol=0, atol=0)
    # strain energy per unit area integrated ply by ply, against e^T K e / 2
    rng = np.random.default_rng(1)
    h = lam.h
    scale = np.array(3*[1.] + 3*[1/h] + 3*[1/h**3] + 2*[1.] + 2*[1/h**2])
    xi, w = np.polynomial.legendre.leggauss(6)
    zbot = -h/2 + lam.offset
    for _ in range(5):
        e = 1e-3*rng.uniform(-1, 1, 13)*scale
        U = 0.
        zk = zbot
        for ply in lam.plies:
            C = np.asarray(ply.get_constitutive_matrix())
            for x, wi in zip(xi, w):
                z = zk + ply.h*(x + 1)/2
                eps = e[0:3] + z*e[3:6] + z**3*e[6:9]
                gamma = e[9:11] + z**2*e[11:13]
                U += wi*ply.h/2*0.5*(eps @ C[:3, :3] @ eps
                                     + gamma @ C[3:, 3:] @ gamma)
            zk += ply.h
        assert np.isclose(0.5*e @ K @ e, U, rtol=1e-12)
