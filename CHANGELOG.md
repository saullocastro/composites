# Changelog

## 0.9.0 (unreleased)

### Breaking: transverse shear stiffness now includes the shear correction

`Laminate.A44`, `A45` and `A55` used to be the uncorrected constant-strain
stiffness `sum_k Cs_k h_k`, and a separate factor (`scf_k13`, `scf_k23`) was
meant to be applied downstream. They now come **with the correction already
applied**, computed by default with the equilibrium approach of

> Rohwer, K. "Improved transverse shear stiffness for layered finite
> elements", DFVLR-FB 88-32, 1988.

which returns the full 2x2 matrix directly, including the `A45` coupling of
angle-ply laminates. It is valid for unsymmetric laminates and does not depend
on `offset`. For a homogeneous plate it gives exactly `5/6 G h`.

- `A44`, `A45`, `A55`: corrected stiffness, ready for FSDT elements.
- `Laminate.Ats` (new): `[[A44, A45], [A45, A55]]`, index 4 <-> yz,
  5 <-> xz.
- `Abar44`, `Abar45`, `Abar55` and `Laminate.Abar_ts` (new): the
  constant-strain, uncorrected stiffness (what `A44`, `A45`, `A55` used to
  be).
- `Abarbar44`, `Abarbar45`, `Abarbar55` and `Laminate.Abarbar_ts` (new): the
  constant-stress stiffness `h^2 [sum_k inv(Cs_k) h_k]^-1`, for comparison.
- `scf_k13`, `scf_k23`: now the **reported ratios** `A55/Abar55` and
  `A44/Abar44`. They are informative only and must not be applied again.
- `Laminate.shear_correction` and the `shear_correction` argument of
  `laminated_plate` and `isotropic_plate` (new): `'rohwer'` (default),
  `'vlachoutsis'`, `'constant'` (5/6) or `None` (no correction).
- `Laminate.calc_transverse_shear_stiffness()` (new): called by
  `calc_constitutive_matrix()`. It raises `ValueError` if a ply has
  `g13 = 0` or `g23 = 0` (singular `Cs`), or if the ABD matrix is singular.
  Use `shear_correction=None` for plies without transverse shear properties.
- `Laminate.calc_transverse_shear_stress(z, Qy, Qx)` (new): recovers
  `(tau_yz, tau_xz)` through the thickness from the same equilibrium
  distribution, zero at the free surfaces and continuous across plies, for
  use in failure criteria.
- Laminates created from lamination parameters
  (`laminate_from_lamination_parameters`,
  `laminate_from_LaminationParameters`) have no ply distribution. They keep
  `A44 = Abar44`, etc., with `shear_correction = None` and
  `scf_k13 = scf_k23 = 1`.
- `Laminate` has new C attributes, so packages that `cimport composites` must
  be recompiled against this version.

### Fixed

The previous `Laminate.calc_scf`, which claimed to implement Vlachoutsis
(1992), returned wrong factors, often above 1, which a shear correction factor
can never be. Specifically, it:

1. accumulated the ply bending stiffness over the plies, so the result
   depended on the stacking order (`[0,90]s` gave `(1.2343, 0.4278)`);
2. cancelled the transverse shear moduli, so the factors did not depend on
   `g13`, `g23` (a soft-core sandwich gave `1.2564` instead of `0.0051`);
3. rotated `E1`, `E2` with a non-tensorial formula, wrong except at 0 and 90
   degrees;
4. used `offset` instead of the neutral surface of each direction.

The `'vlachoutsis'` mode is a corrected implementation. It is exact for
specially orthotropic plies, and its `A45 = (k13 + k23)/2 Abar45` is ad hoc.

### Deprecated

- The `calc_scf` argument of `laminated_plate` and `isotropic_plate`:
  `calc_scf=True` maps to `shear_correction='rohwer'` and `calc_scf=False` to
  `shear_correction=None`.
- `Laminate.calc_scf()`: recomputes the stiffness, switching `None` to
  `'rohwer'`, and returns `(scf_k13, scf_k23)`.
- `Laminate.Atrans`: use `Laminate.Ats`, which returns the same corrected
  matrix.

### Migration

- FSDT elements: use `lam.Ats` (or `A44`, `A45`, `A55`) directly and **remove
  any multiplication by `scf_k13`, `scf_k23` or 5/6**, otherwise the
  correction is applied twice.
- To reproduce the previous values of `A44`, `A45`, `A55`, read `Abar44`,
  `Abar45`, `Abar55` (or `lam.Abar_ts`), or pass `shear_correction=None`.
- TSDT: the transverse shear terms that go with `Dtrans` and `Ftrans` are
  the uncorrected `lam.Abar_ts`, because the third-order theory needs no
  shear correction.
- Replace `calc_scf=True` with the default and `calc_scf=False` with
  `shear_correction=None`. Replace `lam.Atrans` with `lam.Ats`.
- Through-thickness transverse shear stresses: use
  `lam.calc_transverse_shear_stress(z, Qy, Qx)` rather than `Cs @ gamma`,
  which is constant within each ply and non-zero at the free surfaces.
