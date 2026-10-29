r"""Plane-stress invariants of MatLamina when e3, nu13 and nu23 are given

The stiffnesses of a :class:`.Laminate` are those of plane stress,
`\sigma_{33} = 0`: the reduced stiffnesses depend on E1, E2, nu12, G12, G13
and G23 only (Reddy, J. N., Mechanics of Laminated Composite Plates and
Shells: Theory and Analysis, 2nd ed., CRC Press, 2003, Eq. (2.4.4), p. 101;
Jones, R. M., Mechanics of Composite Materials, 2nd ed., Taylor & Francis,
1999, Eqs. (2.64) and (2.66), p. 72), and the invariants of Tsai and Pagano
are those of the reduced stiffnesses (Jones, 1999, Eq. (2.99), p. 86).

Reference stiffnesses computed in the tests, independently of composites,
with the classical lamination theory

    (A_ij, B_ij, D_ij) = sum_k int_{z_k}^{z_k+1} Qbar_ij^(k) (1, z, z^2) dz

(Reddy, 2003, Eq. (3.3.38), p. 128; Jones, 1999, Eq. (4.24), p. 198), where
Qbar^(k) is the plane-stress reduced stiffness of ply k in the laminate
axes, see ``qbar``.
"""
import sys
sys.path.append('..')

import numpy as np
import pytest

from composites.utils import read_laminaprop, laminated_plate
from composites.core import (laminate_from_LaminationParameters,
                             laminate_from_lamination_parameters, GradABD,
                             LaminationParameters)


def qbar(E1, E2, nu12, G12, thetadeg):
    r"""Plane-stress reduced stiffness of a ply in the laminate axes

    ``Q = [[E1, nu12 E2, 0], [nu12 E2, E2, 0], [0, 0, G12 (1 - nu12
    nu21)]]/(1 - nu12 nu21)`` in the material axes, with ``nu21 = nu12
    E2/E1`` (Jones, 1999, Eq. (2.66), p. 72), and ``Qbar = Ts^-1 Q Te``
    (Jones, 1999, Eq. (2.83), p. 77, with Reuter's matrix), where ``Ts`` and
    ``Te`` transform the stresses and the engineering strains ``{xx, yy,
    xy}`` from the laminate axes to the material axes, whose direction 1 is
    at ``thetadeg`` from the laminate x axis, counterclockwise. With ``E1 =
    0``, e.g. an antiplane core, ``nu12 E2 = 0`` and ``nu12 nu21 = 0``, the
    limit of a positive-definite material, ``|nu12| < (E1/E2)**0.5``
    (Jones, 1999, Eq. (2.51), p. 69).
    """
    if E1 != 0:
        nu21 = nu12*E2/E1
    else:
        assert nu12*E2 == 0
        nu21 = 0.
    den = 1 - nu12*nu21
    Q = np.array([[E1/den, nu12*E2/den, 0.],
                  [nu12*E2/den, E2/den, 0.],
                  [0., 0., G12]])
    c = np.cos(np.deg2rad(thetadeg))
    s = np.sin(np.deg2rad(thetadeg))
    Ts = np.array([[c*c, s*s, 2*c*s], [s*s, c*c, -2*c*s],
                   [-c*s, c*s, c*c - s*s]])
    Te = np.array([[c*c, s*s, c*s], [s*s, c*c, -c*s],
                   [-2*c*s, 2*c*s, c*c - s*s]])
    return np.linalg.solve(Ts, Q @ Te)


def clt_abd(stack, plyts, laminaprops, offset=0.):
    r"""A, B, D of the classical lamination theory, plies from the bottom,
    ``laminaprops`` with the entries ``(E1, E2, nu12, G12)`` first"""
    h = sum(plyts)
    z = offset - h/2 + np.concatenate(([0.], np.cumsum(plyts)))
    A = np.zeros((3, 3))
    B = np.zeros((3, 3))
    D = np.zeros((3, 3))
    for k, (thetadeg, lp) in enumerate(zip(stack, laminaprops)):
        Qb = qbar(lp[0], lp[1], lp[2], lp[3], thetadeg)
        A += Qb*(z[k+1] - z[k])
        B += Qb*(z[k+1]**2 - z[k]**2)/2
        D += Qb*(z[k+1]**3 - z[k]**3)/3
    return A, B, D


def abar_ts(stack, plyts, G13, G23):
    r"""Constant-strain transverse shear stiffness ``sum_k Qbar_s^(k) h_k``,
    ``[[Qbar44, Qbar45], [Qbar45, Qbar55]]`` of Reddy (2003), Eq. (3.4.17b),
    p. 137, with ``Q44 = G23`` and ``Q55 = G13``, Eq. (2.4.4b), p. 101"""
    Abar = np.zeros((2, 2))
    for thetadeg, hk in zip(stack, plyts):
        c = np.cos(np.deg2rad(thetadeg))
        s = np.sin(np.deg2rad(thetadeg))
        Abar += hk*np.array([[G23*c*c + G13*s*s, (G13 - G23)*c*s],
                             [(G13 - G23)*c*s, G23*s*s + G13*c*c]])
    return Abar


def assert_abd(A, B, D, A_ref, B_ref, D_ref, h):
    r"""Round-off tolerance, where the terms that vanish are compared with
    an absolute tolerance relative to the scale of the respective matrix"""
    scale = np.abs(A_ref).max()
    assert np.allclose(A, A_ref, rtol=1e-12, atol=1e-14*scale)
    assert np.allclose(B, B_ref, rtol=1e-12, atol=1e-14*scale*h)
    assert np.allclose(D, D_ref, rtol=1e-12, atol=1e-14*scale*h**2)


# Materials given with the nine constants of read_laminaprop, e3 != 0: an
# isotropic one and an orthotropic one with E2 != E3, nu13 != nu23 and G13 !=
# G23
G_ISO = 71e9/(2*(1 + 0.33))
LAMINAPROPS_9 = {
    'isotropic': (71e9, 71e9, 0.33, G_ISO, G_ISO, G_ISO, 71e9, 0.33, 0.33),
    'orthotropic': (150e9, 8.7e9, 0.28, 5.0e9, 5.5e9, 3.1e9, 10.1e9, 0.3,
                    0.45),
}
# stack from the bottom ply, ply thicknesses
LAYUPS = {
    'cross-ply symmetric': ([0, 90, 90, 0], [0.125e-3]*4),
    'unsymmetric': ([0, 45, 90], [0.125e-3]*3),
    'antisymmetric': ([45, -45], [0.125e-3]*2),
    'unbalanced unsymmetric': ([15, -70, 35, 90, 0],
                               [0.1e-3, 0.2e-3, 0.125e-3, 0.15e-3, 0.3e-3]),
}


@pytest.mark.parametrize('layup', LAYUPS.keys())
@pytest.mark.parametrize('material', LAMINAPROPS_9.keys())
def test_lamination_parameters_three_dimensional_constants(material, layup):
    r"""Both routes give the plane-stress stiffnesses when e3 is given

    :func:`.laminated_plate` and :func:`.laminate_from_LaminationParameters`,
    with the lamination parameters of the same layup, must give A, B and D
    of the classical lamination theory, ``clt_abd``, and the constant-strain
    transverse shear stiffness ``abar_ts``, to round-off, also for e3 != 0,
    and the same values as for the six-entry ``laminaprop`` without e3, nu13
    and nu23.
    """
    lamprop = LAMINAPROPS_9[material]
    stack, plyts = LAYUPS[layup]
    h = sum(plyts)
    # NOTE no shear correction, such that Ats = Abar_ts in both routes
    lam = laminated_plate(stack, plyts=plyts, laminaprop=lamprop,
                          shear_correction=None)
    mat = read_laminaprop(lamprop)
    assert mat.e3 != 0 and mat.plane_stress
    lp = lam.calc_lamination_parameters()
    lam_lp = laminate_from_LaminationParameters(lam.h, mat, lp)

    A, B, D = clt_abd(stack, plyts, [lamprop]*len(stack))
    Abar = abar_ts(stack, plyts, lamprop[4], lamprop[5])
    for p in (lam, lam_lp):
        assert_abd(p.A, p.B, p.D, A, B, D, h)
        assert np.allclose(p.Abar_ts, Abar, rtol=1e-12,
                           atol=1e-14*np.abs(Abar).max())
    assert np.allclose(lam_lp.Ats, lam.Ats, rtol=1e-12,
                       atol=1e-14*np.abs(Abar).max())

    # NOTE e3, nu13 and nu23 do not change the plane-stress invariants
    mat6 = read_laminaprop(lamprop[:6])
    for name in ('q11', 'q12', 'q21', 'q22', 'q44', 'q55', 'q66', 'u1', 'u2',
                 'u3', 'u4', 'u5', 'u6', 'u7'):
        assert getattr(mat, name) == getattr(mat6, name), name

    # NOTE the gradients use the same invariants: at constant lamination
    #      parameters A is linear in h and D cubic, dA/dh = A/h and dD/dh =
    #      3 D/h
    grad = GradABD()
    grad.calc_LP_grad(lam.h, mat, lp)
    ij = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
    gradA = np.array([A[i, j] for i, j in ij])/h
    gradD = np.array([D[i, j] for i, j in ij])*3/h
    assert np.allclose(np.asarray(grad.gradAij)[:, 0], gradA, rtol=1e-12,
                       atol=1e-14*np.abs(gradA).max())
    assert np.allclose(np.asarray(grad.gradDij)[:, 0], gradD, rtol=1e-12,
                       atol=1e-14*np.abs(gradD).max())

    # NOTE three-dimensional option: stiffer than plane stress, the
    #      lamination-parameter route then no longer gives the plane-stress
    #      stiffnesses
    mat.rebuild(plane_stress=False)
    lam_3d = laminate_from_LaminationParameters(lam.h, mat, lp)
    assert lam_3d.A11 > (1 + 1e-3)*lam.A11
    assert lam_3d.A22 > (1 + 1e-3)*lam.A22


@pytest.mark.parametrize('offset', ['+0.2 mm', '-0.2 mm', 'one ply'])
@pytest.mark.parametrize('layup', LAYUPS.keys())
@pytest.mark.parametrize('material', LAMINAPROPS_9.keys())
def test_lamination_parameters_offset(material, layup, offset):
    r"""Lamination parameters of a laminate with an offset

    With ``offset`` `= e`, the laminate occupies `-h/2 + e \le z \le h/2 +
    e`, and :meth:`.Laminate.calc_lamination_parameters` integrates about
    the reference surface. :func:`.laminate_from_LaminationParameters`
    given the same offset must give the A, B and D of ``clt_abd`` with the
    offset, which also satisfy the shift of the stiffness integrals to
    another reference plane, `B = B^0 + e A` and `D = D^0 + 2 e B^0 + e^2
    A`, with `B^0`, `D^0` about the mid-surface (Vasiliev and Morozov,
    Mechanics and Analysis of Composite Materials, Elsevier, 2001, Eq.
    (5.28), p. 235). The constant-strain transverse shear stiffness does not
    depend on the offset. The derivatives of GradABD with respect to the
    thickness, at constant lamination parameters and offset, are compared
    with central differences of the stiffnesses, polynomials of degree 3 in
    h, whose truncation error is `(\Delta h/h)^2 = 10^{-12}` and whose
    round-off error is about `10^{-16} h/\Delta h = 10^{-10}`, relative.
    """
    lamprop = LAMINAPROPS_9[material]
    stack, plyts = LAYUPS[layup]
    e = {'+0.2 mm': 0.2e-3, '-0.2 mm': -0.2e-3, 'one ply': plyts[0]}[offset]
    h = sum(plyts)
    lam = laminated_plate(stack, plyts=plyts, laminaprop=lamprop, offset=e,
                          shear_correction=None)
    mat = read_laminaprop(lamprop)
    lp = lam.calc_lamination_parameters()
    lam_lp = laminate_from_LaminationParameters(lam.h, mat, lp, e)
    assert lam_lp.offset == e

    A, B, D = clt_abd(stack, plyts, [lamprop]*len(stack), offset=e)
    A0, B0, D0 = clt_abd(stack, plyts, [lamprop]*len(stack))
    assert_abd(A, B, D, A0, B0 + e*A0, D0 + 2*e*B0 + e**2*A0, h)
    Abar = abar_ts(stack, plyts, lamprop[4], lamprop[5])
    for p in (lam, lam_lp):
        assert_abd(p.A, p.B, p.D, A, B, D, h)
        assert np.allclose(p.Abar_ts, Abar, rtol=1e-12,
                           atol=1e-14*np.abs(Abar).max())
    assert np.allclose(lam_lp.Ats, lam.Ats, rtol=1e-12,
                       atol=1e-14*np.abs(Abar).max())

    # NOTE the same through the function with the lamination parameters as
    #      arguments
    names = ['xi%s%d' % (m, i) for m in 'ABD' for i in (1, 2, 3, 4)]
    lam_lp2 = laminate_from_lamination_parameters(lam.h, mat,
            *[getattr(lp, name) for name in names], lp.xiAts1, lp.xiAts2,
            offset=e)
    assert_abd(lam_lp2.A, lam_lp2.B, lam_lp2.D, A, B, D, h)

    grad = GradABD()
    grad.calc_LP_grad(lam.h, mat, lp, e)
    ij = ('11', '12', '16', '22', '26', '66')
    def stiffnesses(hh):
        p = laminate_from_LaminationParameters(hh, mat, lp, e)
        return {m: np.array([getattr(p, m + s) for s in ij]) for m in 'ABD'}
    dh = 1e-6*h
    plus, minus = stiffnesses(h + dh), stiffnesses(h - dh)
    for m, g in (('A', grad.gradAij), ('B', grad.gradBij),
                 ('D', grad.gradDij)):
        fd = (plus[m] - minus[m])/(2*dh)
        assert np.allclose(np.asarray(g)[:, 0], fd, rtol=1e-8,
                           atol=1e-8*np.abs(fd).max()), m


@pytest.mark.parametrize('material', LAMINAPROPS_9.keys())
def test_trace_normalized_three_dimensional_constants(material):
    r"""Trace normalization with e3 given

    The trace of Melo et al. (2017), Eq. (1), p. 45,
    https://doi.org/10.1016/j.compstruct.2016.09.055, is that of the
    plane-stress stiffness matrix, ``tr = Q11 + Q22 + 2 Q66`` with ``Q11 =
    E1/(1 - nu12 nu21)``, ``Q22 = E2/(1 - nu12 nu21)``, ``Q66 = G12`` (Jones,
    1999, Eq. (2.66), p. 72), independent of e3. The stiffnesses of the
    lamination-parameter route are linear in the invariants, such that with
    the trace-normalized invariants they are those of laminated_plate
    divided by tr. After ``rebuild(plane_stress=False)`` the trace
    normalization is refused.
    """
    lamprop = LAMINAPROPS_9[material]
    E1, E2, nu12, G12 = lamprop[:4]
    nu21 = nu12*E2/E1
    tr = (E1 + E2)/(1 - nu12*nu21) + 2*G12
    stack, plyts = LAYUPS['unbalanced unsymmetric']
    lam = laminated_plate(stack, plyts=plyts, laminaprop=lamprop,
                          shear_correction=None)
    m = read_laminaprop(lamprop)
    assert np.isclose(m.q11 + m.q22 + 2*m.q66, tr, rtol=1e-14, atol=0)
    m.trace_normalize_plane_stress()
    assert np.isclose(2*(m.u1 + m.u5), 1., rtol=1e-14, atol=0)
    lam_lp = laminate_from_LaminationParameters(lam.h, m,
            lam.calc_lamination_parameters())
    assert_abd(lam_lp.A*tr, lam_lp.B*tr, lam_lp.D*tr, lam.A, lam.B, lam.D,
               lam.h)
    assert np.allclose(lam_lp.Abar_ts*tr, lam.Abar_ts, rtol=1e-12,
                       atol=1e-14*np.abs(lam.Abar_ts).max())

    m3d = read_laminaprop(lamprop)
    m3d.rebuild(plane_stress=False)
    with pytest.raises(ValueError, match='plane-stress stiffnesses'):
        m3d.trace_normalize_plane_stress()


@pytest.mark.parametrize('material', LAMINAPROPS_9.keys())
def test_plane_stress_attribute_saved(material):
    r"""``plane_stress`` survives pickling and JSON

    It is read-only, pickled with the other attributes, and restored from
    ``q33`` when loading from JSON, see :mod:`composites.json_io`.
    """
    import pickle
    from composites.json_io import to_json, from_json
    m = read_laminaprop(LAMINAPROPS_9[material])
    with pytest.raises(AttributeError):
        m.plane_stress = False
    for plane_stress in (True, False):
        m.rebuild(plane_stress=plane_stress)
        for m2 in (pickle.loads(pickle.dumps(m)), from_json(to_json(m))):
            assert m2.plane_stress is plane_stress
            assert m2.q11 == m.q11 and m2.u1 == m.u1


# Antiplane cores: "an idealized core in which the modulus of elasticity in
# planes parallel with the faces is zero but the shear modulus in planes
# perpendicular to the faces is finite" (Allen, H. G., Analysis and Design of
# Structural Sandwich Panels, Pergamon Press, 1969, p. 14), with G13 != G23,
# given with 6 and with 9 entries, the latter with a through-thickness
# modulus E3, as a honeycomb
G13_CORE, G23_CORE = 50e6, 30e6
ANTIPLANE_CORES = {
    '6 entries': (0., 0., 0.3, 0., G13_CORE, G23_CORE),
    '9 entries': (0., 0., 0.3, 0., G13_CORE, G23_CORE, 1.2e9, 0., 0.),
}


@pytest.mark.parametrize('core', ANTIPLANE_CORES.keys())
def test_antiplane_core_matlamina(core):
    r"""MatLamina of an antiplane core, e1 = e2 = g12 = 0

    The plane-stress reduced stiffnesses of Jones (1999), Eq. (2.66), p. 72,
    vanish, Q11 = Q12 = Q22 = Q66 = 0, and Q44 = G23, Q55 = G13, such that
    u1 to u5 vanish and u6 = (G23 + G13)/2, u7 = (G23 - G13)/2. The
    three-dimensional stiffness matrix, the limit e1, e2 -> 0 of the inverse
    compliance, has C33 = E3, C44 = G23, C55 = G13 and all other terms zero.
    The reciprocal Poisson's ratios are stored as zero, and the trace of
    Melo et al. (2017), Eq. (1), being zero, the trace normalization is
    refused.
    """
    lamprop = ANTIPLANE_CORES[core]
    e3 = lamprop[6] if len(lamprop) == 9 else 0.
    m = read_laminaprop(lamprop)
    assert m.nu21 == 0 and m.nu31 == 0 and m.nu32 == 0
    for name in ('q11', 'q12', 'q13', 'q21', 'q22', 'q23', 'q31', 'q32',
                 'q33', 'q66', 'u1', 'u2', 'u3', 'u4', 'u5'):
        assert getattr(m, name) == 0, name
    assert (m.q44, m.q55) == (G23_CORE, G13_CORE)
    assert m.u6 == (G23_CORE + G13_CORE)/2
    assert m.u7 == (G23_CORE - G13_CORE)/2
    C_ref = np.zeros((6, 6))
    C_ref[2, 2] = e3
    C_ref[3, 3] = G23_CORE
    C_ref[4, 4] = G13_CORE
    C = np.asarray(m.get_constitutive_matrix())
    assert np.all(np.isfinite(C))
    assert np.array_equal(C, C_ref)
    # NOTE the limit of a positive-definite material, e1 = e2 = 1e-8 E3
    if e3 > 0:
        eps = 1e-8*e3
        m_eps = read_laminaprop((eps, eps) + lamprop[2:])
        C_eps = np.asarray(m_eps.get_constitutive_matrix())
        assert np.allclose(C, C_eps, rtol=0, atol=1e-7*e3)
    with pytest.raises(ValueError, match='not positive'):
        m.trace_normalize_plane_stress()
    # three-dimensional option
    m.rebuild(plane_stress=False)
    assert m.q33 == e3
    for name in ('q11', 'q12', 'q13', 'q22', 'q23', 'u1', 'u2', 'u3', 'u4',
                 'u5'):
        assert getattr(m, name) == 0, name


def test_zero_young_modulus_errors():
    r"""Undefined Poisson's ratios with a zero Young modulus

    A positive-definite material has ``|nu_ij| < (E_i/E_j)**0.5`` (Jones,
    1999, Eq. (2.51), p. 69), such that ``E_i = 0`` with ``E_j > 0``
    requires ``nu_ij = 0``; otherwise ``nu_ji = nu_ij E_j/E_i`` is undefined
    and ``read_laminaprop`` raises ``ValueError``. With ``E_i = E_j = 0``
    any ``nu_ij`` is accepted. A laminate of antiplane cores only has no
    in-plane and bending stiffness, and its equivalent properties and its
    shear correction factors are undefined.
    """
    for lamprop, match in [
            ((0., 8.7e9, 0.28, 5e9, 5e9, 3e9), 'requires nu12 = 0'),
            ((0., 0., 0.3, 0., 5e7, 3e7, 1e9, 0.2, 0.), 'requires nu13 = 0'),
            ((142e9, 0., 0.28, 0., 5e7, 3e7, 1e9, 0.2, 0.3),
             'requires nu23 = 0')]:
        with pytest.raises(ValueError, match=match):
            read_laminaprop(lamprop)
    # accepted: e1 = 0 with nu12 = 0, and e1 = e2 = 0 with any nu12
    m = read_laminaprop((0., 8.7e9, 0., 5e9, 5e9, 3e9))
    assert (m.q11, m.q12, m.q22, m.q66) == (0., 0., 8.7e9, 5e9)
    read_laminaprop((0., 0., 0.45, 0., 5e7, 3e7))
    core = ANTIPLANE_CORES['6 entries']
    with pytest.raises(ValueError, match='singular'):
        laminated_plate([0, 30], plyts=[1e-3, 2e-3], laminaprop=core,
                        shear_correction=None)
    with pytest.raises(ValueError, match='singular'):
        laminated_plate([0, 30], plyts=[1e-3, 2e-3], laminaprop=core)


@pytest.mark.parametrize('core', ANTIPLANE_CORES.keys())
@pytest.mark.parametrize('stack', [[0, 0, 0], [0, 30, 90], [45, -20, 0]])
def test_antiplane_core_sandwich(core, stack):
    r"""Sandwich with an antiplane core, unequal faces and an offset

    The A, B and D of :func:`.laminated_plate` are those of the classical
    lamination theory, ``clt_abd``, in which the core contributes nothing,
    and the constant-strain transverse shear stiffness is ``abar_ts``, where
    only G13 and G23 of the core enter, rotated by its angle. All are
    finite. For a core of a single material, the lamination-parameter route
    gives A = B = D = 0 and the same constant-strain transverse shear
    stiffness, for any offset.
    """
    face = (71e9, 7e9, 0.28, 4e9, 4e9, 2.5e9)
    lamprop = ANTIPLANE_CORES[core]
    plyts = [0.5e-3, 10e-3, 0.8e-3]
    offset = 0.3e-3
    laminaprops = [face, lamprop, face]
    lam = laminated_plate(stack, plyts=plyts, laminaprops=laminaprops,
                          offset=offset)
    A, B, D = clt_abd(stack, plyts, laminaprops, offset=offset)
    assert_abd(lam.A, lam.B, lam.D, A, B, D, sum(plyts))
    Abar = np.zeros((2, 2))
    for thetadeg, hk, lp in zip(stack, plyts, laminaprops):
        Abar += abar_ts([thetadeg], [hk], lp[4], lp[5])
    assert np.allclose(lam.Abar_ts, Abar, rtol=1e-12,
                       atol=1e-14*np.abs(Abar).max())
    for name in ('ABD', 'Ats', 'Abarbar_ts'):
        assert np.all(np.isfinite(getattr(lam, name))), name
    # NOTE no shear correction factor exceeds one, and the core carries the
    #      transverse shear: the corrected stiffness is far below Abar
    assert 0 < lam.scf_k13 < 1 and 0 < lam.scf_k23 < 1
    # NOTE lamination parameters of the core alone
    hc = plyts[1]
    for e in (0., offset):
        lp = LaminationParameters()
        c2, s2 = (np.cos(np.deg2rad(2*stack[1])),
                  np.sin(np.deg2rad(2*stack[1])))
        lp.xiA1, lp.xiA2, lp.xiAts1, lp.xiAts2 = c2, s2, c2, s2
        lp.xiB1, lp.xiB2 = 4*e/hc*c2, 4*e/hc*s2
        lp.xiD1, lp.xiD2 = (1 + 12*(e/hc)**2)*c2, (1 + 12*(e/hc)**2)*s2
        lp.xiA3 = lp.xiB3 = lp.xiD3 = 1.
        lam_lp = laminate_from_LaminationParameters(hc,
                read_laminaprop(lamprop), lp, e)
        for name in ('A', 'B', 'D'):
            assert np.array_equal(getattr(lam_lp, name), np.zeros((3, 3)))
        Abar_c = abar_ts([stack[1]], [hc], lamprop[4], lamprop[5])
        assert np.allclose(lam_lp.Abar_ts, Abar_c, rtol=1e-12,
                           atol=1e-14*np.abs(Abar_c).max())
