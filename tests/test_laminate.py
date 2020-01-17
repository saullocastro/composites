import sys
sys.path.append('..')

import numpy as np

from composites.laminate import (read_lamination_parameters, read_stack,
        read_isotropic)


def test_lampar():
    lamprop = (71e9, 71e9, 0.33)
    lam = read_lamination_parameters(1, lamprop, None,
        0.5, 0.4, -0.3, -0.6,
        0.5, 0.4, -0.3, -0.6,
        0.5, 0.4, -0.3, -0.6,
        0.5, 0.4, -0.3, -0.6)
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
    lam.calc_lamination_parameters()
    lam.calc_ABDE_from_lamination_parameters()
    assert np.allclose(lam.A, A)
    assert np.allclose(lam.B, B)
    assert np.allclose(lam.D, D)
    assert np.allclose(lam.E, E)


def test_read_stack():
    lamprop = (71e9, 7e9, 0.28, 7e9, 7e9, 7e9)
    stack = [0, 45, 90]
    plyt = 0.000125
    lam = read_stack(stack, plyt, lamprop)
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
    lam.calc_equivalent_modulus()
    lam.calc_lamination_parameters()
    lam.calc_ABDE_from_lamination_parameters()
    #TODO A, B and D are changing from the original, check!
    A = np.array([[ 13589503.90225179,  2502486.88587513,  2026742.01957523],
                  [  2502486.88587513, 13589503.90225179,  2026742.01957523],
                  [  2026742.01957523,  2026742.01957523,  4084254.25409417]])
    B = np.array([[ -1.01337101e+03,  0.00000000e+00,  8.68715094e-15],
                  [  0.00000000e+00,  1.01337101e+03,  5.33639272e-14],
                  [  8.68715094e-15,  5.33639272e-14,  0.00000000e+00]])
    D = np.array([[ 0.17445256, 0.01412545, 0.00263899],
                  [ 0.01412545, 0.17445256, 0.00263899],
                  [ 0.00263899, 0.00263899, 0.03266179]])
    E = np.array([[ 2625000.,       0.],
                  [       0., 2625000.]])
    assert np.allclose(lam.A, A)
    assert np.allclose(lam.B, B)
    assert np.allclose(lam.D, D)
    assert np.allclose(lam.E, E)

    lam.force_symmetric()
    assert np.allclose(lam.B, 0*B)

    lam.calc_ABDE_from_lamination_parameters()
    lam.force_orthotropic()
    A = np.array([[ 13589503.90225179,  2502486.88587513,  0],
                  [  2502486.88587513, 13589503.90225179,  0],
                  [  0,  0,  4084254.25409417]])
    B = np.array([[ -1.01337101e+03,  0.00000000e+00,  0],
                  [  0.00000000e+00,  1.01337101e+03,  0],
                  [  0,  0,  0.00000000e+00]])
    D = np.array([[ 0.17445256, 0.01412545, 0],
                  [ 0.01412545, 0.17445256, 0],
                  [ 0, 0, 0.03266179]])
    assert np.allclose(lam.A, A)
    assert np.allclose(lam.B, B)
    assert np.allclose(lam.D, D)
    lam.force_balanced_LP()
    lam.force_symmetric_LP()

def test_read_isotropic():
    E = 71e9
    nu = 0.28
    thick = 0.000125
    lam = read_isotropic(thickness=thick, E=E, nu=nu)
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

if __name__ == '__main__':
    test_lampar()
    test_read_stack()
    test_read_isotropic()

