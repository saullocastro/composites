import sys
sys.path.append('..')

import numpy as np

from composites.utils import (read_laminaprop, laminated_plate,
        isotropic_plate)
from composites.core import (laminate_from_LaminationParameters,
        laminate_from_lamination_parameters, force_balanced_LP,
        force_orthotropic_LP, force_symmetric_LP, Lamina,
        GradABDE, LaminationParameters)


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
    E = np.array([[2.66917293e+10,  0.00000000e+00],
                  [0.00000000e+00,  2.66917293e+10]])
    assert np.allclose(lam.A, A)
    lam.force_symmetric()
    assert np.allclose(lam.B, B)
    assert np.allclose(lam.D, D)
    assert np.allclose(lam.E, E)
    ABD = lam.ABD
    assert np.allclose(ABD[:3, :3], A)
    assert np.allclose(ABD[3:, 3:], D)
    ABDE = lam.ABDE
    assert np.allclose(ABDE[:3, :3], A)
    assert np.allclose(ABDE[3:6, 3:6], D)


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
    E = np.array([[2.66917293e+10, 0.00000000e+00],
                  [0.00000000e+00, 2.66917293e+10]])
    assert np.allclose(lam.A, A)
    lam.force_symmetric()
    assert np.allclose(lam.B, B)
    assert np.allclose(lam.D, D)
    assert np.allclose(lam.E, E)
    ABD = lam.ABD
    assert np.allclose(ABD[:3, :3], A)
    assert np.allclose(ABD[3:, 3:], D)
    ABDE = lam.ABDE
    assert np.allclose(ABDE[:3, :3], A)
    assert np.allclose(ABDE[3:6, 3:6], D)


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
    E = np.array([[ 2625000.,       0.],
                  [       0., 2625000.]])
    assert np.allclose(lam.A, A)
    assert np.allclose(lam.B, B)
    assert np.allclose(lam.D, D)
    assert np.allclose(lam.E, E)


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
    E = np.array([[ 2625000.,       0.],
                  [       0., 2625000.]])
    assert np.allclose(lam.A, A)
    assert np.allclose(lam.B, B)
    assert np.allclose(lam.D, D)
    assert np.allclose(lam.E, E)
    lam.calc_scf()
    lam.calc_equivalent_properties()
    lp = lam.calc_lamination_parameters()
    matlamina = lam.plies[0].matlamina
    thickness = lam.h
    lam_2 = laminate_from_LaminationParameters(thickness, matlamina, lp)
    assert np.allclose(lam_2.A, lam.A)
    assert np.allclose(lam_2.B, lam.B)
    assert np.allclose(lam_2.D, lam.D)
    assert np.allclose(lam_2.E, lam.E)

    lam.force_balanced()
    force_balanced_LP(lp)
    lam_2 = laminate_from_LaminationParameters(thickness, matlamina, lp)
    assert np.allclose(lam_2.A, lam.A)
    assert np.allclose(lam_2.B, lam.B)
    assert np.allclose(lam_2.D, lam.D)
    assert np.allclose(lam_2.E, lam.E)

    lam.force_orthotropic()
    force_orthotropic_LP(lp)
    lam_2 = laminate_from_LaminationParameters(thickness, matlamina, lp)
    assert np.allclose(lam_2.A, lam.A)
    assert np.allclose(lam_2.B, lam.B)
    assert np.allclose(lam_2.D, lam.D)
    assert np.allclose(lam_2.E, lam.E)

    lam = laminated_plate(stack, plyt, lamprop)
    lp = lam.calc_lamination_parameters()
    lam.force_symmetric()
    force_symmetric_LP(lp)
    lam_2 = laminate_from_LaminationParameters(thickness, matlamina, lp)
    assert np.allclose(lam_2.A, lam.A)
    assert np.allclose(lam_2.B, lam.B)
    assert np.allclose(lam_2.D, lam.D)
    assert np.allclose(lam_2.E, lam.E)


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
    E = np.array([[3466796.875,       0.   ],
                  [      0.   , 3466796.875]])
    assert np.allclose(lam.A, A)
    assert np.allclose(lam.B, 0)

    assert np.allclose(lam.D, D)
    assert np.allclose(lam.E, E)



def test_errors():
    E = 71e9
    nu = 0.28
    thick = 0.000125
    lam = isotropic_plate(thickness=thick, E=E, nu=nu)
    lam.offset = 1.
    try:
        lam.force_balanced()
    except RuntimeError:
        pass
    try:
        lam.force_orthotropic()
    except RuntimeError:
        pass
    try:
        lam.force_symmetric()
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
    gradABDE = GradABDE()
    gradABDE.calc_LP_grad(thickness, matlamina, lp)
    print(gradABDE.gradAij)
