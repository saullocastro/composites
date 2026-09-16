r"""
Transverse shear stiffness / shear correction factors for FSDT laminated plates.

Implements, in the notation of the section
"Transverse shear correction factors" (\label{sec:plate-shear-correction}):

  * Eq. (plate-shear-three-stiffnesses)
        Abar   : constant-strain stiffness   Abar   = sum_k Cs_k h_k
        Abarbar: constant-stress stiffness   Abarbar= h^2 [sum_k inv(Cs_k) h_k]^-1
  * Eq. (plate-shear-Ats-final)   Rohwer (1988) equilibrium-based 2x2 stiffness
  * Eq. (plate-shear-vlach-K)     Vlachoutsis (1992) scalar factors k13, k23
  * Eq. (plate-shear-vlach-perply) per-ply factors
  * Eqs. (plate-shear-sandwich-*) closed-form symmetric-sandwich results

Ordering conventions of the text are respected throughout:

    tau^T = {tau_yz, tau_xz}     Q^T = {Q_y, Q_x}
    Ats   = [[A44, A45], [A45, A55]]        index 4 <-> yz,  index 5 <-> xz

Layer geometry: z is measured from the *reference* surface; the stack goes
from z_1 = -h/2 + d to z_{N+1} = +h/2 + d, where d is the offset of the
reference surface (positive when the mid-surface sits above the reference
surface), matching `composites.Laminate.offset`.

The laminate input is a `composites.Laminate` object
(https://github.com/saullocastro/composites) but any object exposing
`plies` (with `h`, `q11L`, ..., `q55L`), `h` and `offset` works.
"""
import numpy as np

__all__ = [
    "LayerStack", "stack_from_laminate",
    "Ats_constant_strain", "Ats_constant_stress",
    "RohwerResult", "rohwer_Ats",
    "VlachoutsisResult", "vlachoutsis_scf",
    "sandwich_A", "sandwich_closed_form", "sandwich_antiplane",
]

# 3-point Gauss-Legendre on [-1, 1]: exact up to degree 5.
# The integrands below are polynomials of degree 4 in z within each layer,
# so this quadrature is exact (not an approximation).
_GP = np.array([-np.sqrt(3.0 / 5.0), 0.0, np.sqrt(3.0 / 5.0)])
_GW = np.array([5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0])


class LayerStack:
    """Minimal layer-wise description of a laminate.

    Parameters
    ----------
    z : (N+1,) array
        Layer interface coordinates measured from the reference surface.
    C : (N, 3, 3) array
        In-plane constitutive matrix of each layer, rotated to the element
        axes -- Eq. (plate-shear-layer-C).
    Cs : (N, 2, 2) array
        Transverse shear constitutive matrix of each layer, ordered
        [[C44, C45], [C45, C55]] -- Eq. (plate-shear-layer-Cs).
    """

    def __init__(self, z, C, Cs):
        self.z = np.asarray(z, dtype=float)
        self.C = np.asarray(C, dtype=float)
        self.Cs = np.asarray(Cs, dtype=float)
        self.N = self.C.shape[0]
        assert self.z.shape == (self.N + 1,)
        assert self.Cs.shape == (self.N, 2, 2)
        self.h_k = np.diff(self.z)
        self.h = self.z[-1] - self.z[0]

    # -- laminate stiffnesses -------------------------------------------------
    @property
    def ABD(self):
        """6x6 [[A, B], [B, D]] of Eq. (plate-ABD)."""
        z = self.z
        A = np.einsum("k,kij->ij", z[1:] - z[:-1], self.C)
        B = np.einsum("k,kij->ij", (z[1:] ** 2 - z[:-1] ** 2) / 2.0, self.C)
        D = np.einsum("k,kij->ij", (z[1:] ** 3 - z[:-1] ** 3) / 3.0, self.C)
        return np.block([[A, B], [B, D]])


def stack_from_laminate(lam):
    """Build a :class:`LayerStack` from a ``composites.Laminate`` object."""
    z = [-lam.h / 2.0 + lam.offset]
    C, Cs = [], []
    for ply in lam.plies:
        z.append(z[-1] + ply.h)
        C.append([[ply.q11L, ply.q12L, ply.q16L],
                  [ply.q12L, ply.q22L, ply.q26L],
                  [ply.q16L, ply.q26L, ply.q66L]])
        Cs.append([[ply.q44L, ply.q45L],
                   [ply.q45L, ply.q55L]])
    return LayerStack(z, C, Cs)


# ---------------------------------------------------------------------------
# Eq. (plate-shear-three-stiffnesses)
# ---------------------------------------------------------------------------
def Ats_constant_strain(st):
    r"""$\bar{A}_{ts} = \sum_k C_s^{(k)} h_k$ (uncorrected stiffness)."""
    return np.einsum("k,kij->ij", st.h_k, st.Cs)


def Ats_constant_stress(st):
    r"""$\bar{\bar{A}}_{ts} = h^2 [\sum_k (C_s^{(k)})^{-1} h_k]^{-1}$."""
    S = np.einsum("k,kij->ij", st.h_k, np.linalg.inv(st.Cs))
    return st.h ** 2 * np.linalg.inv(S)


# ---------------------------------------------------------------------------
# Rohwer (1988) -- Section \ref{sec:plate-shear-rohwer}
# ---------------------------------------------------------------------------
class RohwerResult:
    """Container for the output of :func:`rohwer_Ats`."""

    def __init__(self, Ats, compliance, a_x, a_y, stack, Pi, Lx, Ly):
        self.Ats = Ats
        self.compliance = compliance
        self._a_x, self._a_y = a_x, a_y
        self._st = stack
        self._Pi, self._Lx, self._Ly = Pi, Lx, Ly

    @property
    def A44(self):
        return self.Ats[0, 0]

    @property
    def A45(self):
        return self.Ats[0, 1]

    @property
    def A55(self):
        return self.Ats[1, 1]

    def f(self, z, k=None):
        """Distribution matrix $f^{(k)}(z)$ of Eq. (plate-shear-f), 2x2.

        ``tau = {tau_yz, tau_xz}^T = f(z) Q`` with ``Q = {Q_y, Q_x}^T``.
        """
        st = self._st
        if k is None:
            k = min(max(np.searchsorted(st.z, z, side="right") - 1, 0), st.N - 1)
        Pi = self._Pi(z)
        row_y = (self._a_y[k] - st.C[k, 1] @ Pi) @ self._Ly
        row_x = (self._a_x[k] - st.C[k, 0] @ Pi) @ self._Lx
        return np.vstack([row_y, row_x])

    def tau(self, z, Q):
        """Transverse shear stresses {tau_yz, tau_xz} at height ``z``."""
        return self.f(z) @ np.asarray(Q, dtype=float)


def rohwer_Ats(stack):
    r"""Equilibrium-based transverse shear stiffness, Eq. (plate-shear-Ats-final).

    Returns
    -------
    RohwerResult
        ``.Ats`` is the 2x2 matrix $[[A_{44}, A_{45}], [A_{45}, A_{55}]]$.
    """
    st = stack
    Hstar = np.linalg.inv(st.ABD)
    Astar = Hstar[:3, :3]
    Bstar = Hstar[:3, 3:]
    BstarT = Hstar[3:, :3]          # = Bstar.T, see remark in the text
    Dstar = Hstar[3:, 3:]

    def Pi(z):
        """Eq. (plate-shear-Pi), 3x6."""
        return np.hstack([z * Astar + 0.5 * z ** 2 * BstarT,
                          z * Bstar + 0.5 * z ** 2 * Dstar])

    # Eq. (plate-shear-relation-factors)
    Lx = np.zeros((6, 2))
    Lx[3, 1] = 1.0       # M_xx,x = Q_x
    Lx[5, 0] = 1.0       # M_xy,x = Q_y
    Ly = np.zeros((6, 2))
    Ly[4, 0] = 1.0       # M_yy,y = Q_y
    Ly[5, 1] = 1.0       # M_xy,y = Q_x

    # Eq. (plate-shear-constants): a^(k) = sum_i (c^(i) - c^(i-1)) Pi(z_i)
    a_x = np.zeros((st.N, 6))
    a_y = np.zeros((st.N, 6))
    acc_x = np.zeros(6)
    acc_y = np.zeros(6)
    c1_prev = np.zeros(3)
    c2_prev = np.zeros(3)
    for k in range(st.N):
        Pi_k = Pi(st.z[k])
        acc_x = acc_x + (st.C[k, 0] - c1_prev) @ Pi_k
        acc_y = acc_y + (st.C[k, 1] - c2_prev) @ Pi_k
        a_x[k] = acc_x
        a_y[k] = acc_y
        c1_prev = st.C[k, 0]
        c2_prev = st.C[k, 1]

    # Eq. (plate-shear-complementary): sum_k int f^T inv(Cs) f dz
    compliance = np.zeros((2, 2))
    Cs_inv = np.linalg.inv(st.Cs)
    for k in range(st.N):
        za, zb = st.z[k], st.z[k + 1]
        zm, dz = 0.5 * (za + zb), 0.5 * (zb - za)
        for xi, w in zip(_GP, _GW):
            z = zm + dz * xi
            Pi_z = Pi(z)
            row_y = (a_y[k] - st.C[k, 1] @ Pi_z) @ Ly
            row_x = (a_x[k] - st.C[k, 0] @ Pi_z) @ Lx
            f = np.vstack([row_y, row_x])
            compliance += w * dz * f.T @ Cs_inv[k] @ f

    Ats = np.linalg.inv(compliance)
    Ats = 0.5 * (Ats + Ats.T)
    return RohwerResult(Ats, compliance, a_x, a_y, st, Pi, Lx, Ly)


# ---------------------------------------------------------------------------
# Vlachoutsis (1992) -- Section \ref{sec:plate-shear-vlachoutsis}
# ---------------------------------------------------------------------------
class VlachoutsisResult:
    def __init__(self, k13, k23, zn, R, d, I, kappa_k, g):
        self.k13, self.k23 = k13, k23
        self.zn = zn            # (zn1, zn2)
        self.R = R              # (R1, R2)
        self.d = d              # (d1, d2) = (Abar55, Abar44)
        self.I = I              # (I1, I2)
        self.kappa_k = kappa_k  # (N, 2) per-ply factors, columns (13, 23)
        self.g = g              # g[alpha](z) callables


def _g_alpha(z_if, Dk, zn):
    r"""Return ``g(z)`` of Eq. (plate-shear-vlach-tau) as a callable."""
    # running value of g at each interface
    g_if = np.zeros(len(z_if))
    for k in range(len(Dk)):
        za, zb = z_if[k], z_if[k + 1]
        g_if[k + 1] = g_if[k] - Dk[k] * (
            0.5 * (zb ** 2 - za ** 2) - zn * (zb - za))

    def g(z):
        k = min(max(np.searchsorted(z_if, z, side="right") - 1, 0), len(Dk) - 1)
        za = z_if[k]
        return g_if[k] - Dk[k] * (0.5 * (z ** 2 - za ** 2) - zn * (z - za))

    return g, g_if


def vlachoutsis_scf(stack):
    r"""Shear correction factors of Eq. (plate-shear-vlach-K).

    ``D_1(z) = C_{11}^{(k)}`` and ``D_2(z) = C_{22}^{(k)}`` are taken from the
    layer matrix already rotated to the element axes, so that the restriction
    to specially orthotropic layers is the user's responsibility.
    """
    st = stack
    z = st.z
    out = {}
    kappa_k = np.zeros((st.N, 2))
    zn_l, R_l, d_l, I_l, g_l = [], [], [], [], []

    for col, (alpha, Gidx) in enumerate([(0, (1, 1)), (1, (0, 0))]):
        # alpha = 0 -> direction 1 (x), uses C11 and G13 = Cs[1,1]
        # alpha = 1 -> direction 2 (y), uses C22 and G23 = Cs[0,0]
        Dk = st.C[:, alpha, alpha].copy()
        Gk = st.Cs[:, Gidx[0], Gidx[1]].copy()

        # Eq. (plate-shear-zn)
        num = np.sum(Dk * (z[1:] ** 2 - z[:-1] ** 2) / 2.0)
        den = np.sum(Dk * (z[1:] - z[:-1]))
        zn = num / den
        R = np.sum(Dk * ((z[1:] - zn) ** 3 - (z[:-1] - zn) ** 3) / 3.0)

        g, _ = _g_alpha(z, Dk, zn)

        d = np.sum(Gk * st.h_k)                      # = Abar55 / Abar44
        I_k = np.zeros(st.N)
        for k in range(st.N):
            za, zb = z[k], z[k + 1]
            zm, dz = 0.5 * (za + zb), 0.5 * (zb - za)
            acc = 0.0
            for xi, w in zip(_GP, _GW):
                acc += w * dz * g(zm + dz * xi) ** 2
            I_k[k] = acc / Gk[k]
        I = I_k.sum()

        kappa = R ** 2 / (d * I)                     # Eq. (plate-shear-vlach-K)
        # Eq. (plate-shear-vlach-perply)
        kappa_k[:, col] = R ** 2 * I_k / ((Gk * st.h_k) * I ** 2)

        out[col] = kappa
        zn_l.append(zn); R_l.append(R); d_l.append(d); I_l.append(I); g_l.append(g)

    return VlachoutsisResult(out[0], out[1], tuple(zn_l), tuple(R_l),
                             tuple(d_l), tuple(I_l), kappa_k, tuple(g_l))


# ---------------------------------------------------------------------------
# Closed-form symmetric sandwich -- Eqs. (plate-shear-sandwich-*)
# ---------------------------------------------------------------------------
def sandwich_A(p):
    r"""$A(p)$ of Eq. (plate-shear-sandwich-p), left-hand expression."""
    return 2.0 * ((1 - p) / 2.0 - (1 - p ** 3) / 3.0 + (1 - p ** 5) / 10.0)


def sandwich_A_factored(p):
    r"""$A(p)$, right-hand (factored) expression of Eq. (plate-shear-sandwich-p)."""
    return (1 - p) ** 3 * (3 * p ** 2 + 9 * p + 8) / 15.0


def sandwich_closed_form(p, Dc_Dfs, Gc_Gfs):
    r"""Eqs. (plate-shear-sandwich-T), (plate-shear-sandwich-K), (-rf).

    Parameters
    ----------
    p : float
        Core-to-total thickness ratio ``h_c/h``.
    Dc_Dfs : float
        Ratio ``D_c / D_fs`` of Eq. (plate-shear-Dalpha).
    Gc_Gfs : float
        Ratio ``G_c / G_fs`` of transverse shear moduli.

    Returns
    -------
    dict with ``kappa``, ``kappa_fs``, ``kappa_c``, ``r_fs`` and the ``T_i``.
    """
    A = sandwich_A(p)
    Gfs_Gc = 1.0 / Gc_Gfs
    T1 = (1 - p ** 3) + p ** 3 * Dc_Dfs
    T2 = Gfs_Gc * (1 - p) + p
    T3 = ((1 - p ** 2) ** 2 + 8.0 / 15.0 * Dc_Dfs ** 2 * p ** 4
          + 4.0 / 3.0 * Dc_Dfs * p ** 2 * (1 - p ** 2))
    T4 = A * Gc_Gfs + p * T3
    T5 = A * Gc_Gfs + p * (1 - p ** 2) ** 2
    T6 = A + Gfs_Gc * p * (1 - p ** 2) ** 2
    T7 = A + Gfs_Gc * p * T3
    res = dict(A=A, T1=T1, T2=T2, T3=T3, T4=T4, T5=T5, T6=T6, T7=T7)
    res["kappa"] = 4.0 / 9.0 * T1 ** 2 / (T2 * T4)
    res["kappa_fs"] = (4.0 / 9.0 * A * T1 ** 2 / ((1 - p) * T7 ** 2)
                       if p != 1 else np.nan)
    res["kappa_c"] = 4.0 / 9.0 * T1 ** 2 * T3 / T4 ** 2
    res["r_fs"] = A * Gc_Gfs / T4
    return res


def sandwich_antiplane(p, Gc_Gfs):
    r"""Eq. (plate-shear-sandwich-antiplane): the ``D_c/D_fs -> 0`` limit."""
    A = sandwich_A(p)
    Gfs_Gc = 1.0 / Gc_Gfs
    T2 = Gfs_Gc * (1 - p) + p
    T5 = A * Gc_Gfs + p * (1 - p ** 2) ** 2
    T6 = A + Gfs_Gc * p * (1 - p ** 2) ** 2
    return dict(
        kappa=4.0 / 9.0 * (1 - p ** 3) ** 2 / (T2 * T5),
        kappa_fs=(4.0 / 9.0 * A * (1 - p ** 3) ** 2 / ((1 - p) * T6 ** 2)
                  if p != 1 else np.nan),
        kappa_c=4.0 / 9.0 * (1 - p ** 3) ** 2 * (1 - p ** 2) ** 2 / T5 ** 2,
        T2=T2, T5=T5, T6=T6, A=A)


def sandwich_simple(p, Gc_Gfs):
    r"""Engineering estimate of Eq. (plate-shear-sandwich-simple)."""
    return 1.0 / (1.0 + ((1 - p) / p) / Gc_Gfs)


# ---------------------------------------------------------------------------
# Convenience wrappers taking a ``composites.Laminate`` directly
# ---------------------------------------------------------------------------
def scf_from_laminate(lam):
    """``(k13, k23)`` of Eq. (plate-shear-vlach-K).

    Drop-in replacement for ``composites.Laminate.calc_scf``::

        lam.scf_k13, lam.scf_k23 = scf_from_laminate(lam)
    """
    v = vlachoutsis_scf(stack_from_laminate(lam))
    return v.k13, v.k23


def Ats_from_laminate(lam):
    """2x2 equilibrium-based ``A_ts`` of Eq. (plate-shear-Ats-final)."""
    return rohwer_Ats(stack_from_laminate(lam)).Ats
