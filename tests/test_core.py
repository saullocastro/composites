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
