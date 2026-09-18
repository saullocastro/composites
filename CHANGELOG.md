# Changelog

## 0.9.1 (2026-09-18)

### Breaking: transverse shear stiffness now includes the shear correction

`Laminate.A44`, `A45` and `A55` used to be the uncorrected constant-strain
stiffness `sum_k Cs_k h_k`, and a separate factor (`scf_k13`, `scf_k23`) was
meant to be applied downstream. They now come **with the correction already
applied**, computed by default with the equilibrium approach of

> Rohwer, K. "Improved transverse shear stiffness for layered finite
> elements", DFVLR-FB 88-32, 1988.

which returns the full 2x2 matrix directly, including the `A45` coupling of
angle-ply laminates. It is valid for unsymmetric laminates and does not depend
on `offset`. For a homogeneous plate it gives exactly `5/6 G h`. For completeness,
the approach of Vlachoutsis is also available:

> Vlachoutsis, S. "Shear correction factors for plates and shells", Int.
> Journal for Numerical Methods in Engineering, Vol. 33, 1537-1552, 1992.

The following important changes are worth mentioning:

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
- `Laminate` and `GradABD` can be pickled and deep-copied. `Laminate` stores
  its public attributes only; the transverse shear distribution used by
  `calc_transverse_shear_stress` is not stored and is recomputed on demand.
  Python subclasses of `Laminate` are supported.

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


## 0.8.6 (2026-04-10)

- Python 3.14 support; the coverage job runs on a newer Python.
- New optimisation tutorials after Le Riche and Haftka: smeared laminate with
  finite differences, smeared laminate with gradients, and lamination
  parameters with gradients.
- Improved verifications in the ghost-layer tutorial.

## 0.8.4 (2025-09-05)

- Build fixes for the Linux distributions, using the latest `cibuildwheel` on
  GitHub Actions.
- Support for the Third-order Shear Deformation Theory (TSDT): new `E`, `F`,
  `H` plane-stress matrices (`E11` ... `H66`) and the transverse terms
  `A44`, `A45`, `A55`, `D44`, `D45`, `D55`, `F44`, `F45`, `F55`, with
  `Atrans`, `Dtrans` and `Ftrans`.
- **Breaking:** the previous `E44`, `E45`, `E55` became `A44`, `A45`, `A55`,
  `get_ABDE` was removed, and `GradABDE` was renamed to `GradABD`.
  Packages that `cimport composites` must be recompiled.

## 0.7.3 (2025-02-12)

- Faster n-double laminate, implemented in Cython (`core.n_double_laminate`).
- **Breaking:** `double_double_plate` and `n_double_plate` were renamed to
  `double_double_laminate` and `n_double_laminate`.
- New tutorial on polar plots of the laminate stiffness.
- Python 3.13 support.

## 0.7.1 (2024-09-02)

- New tutorial on lightweight discrete optimisation with the ghost-layer
  method.
- Separate pytest and coverage workflows.

## 0.7.0 (2024-03-15)

- New `composites.kassapoglou` module, after Kassapoglou's *Design and
  Analysis of Composite Structures*: buckling under uniaxial compression
  (`calc_Nxx_crit`), under shear (`calc_Nxy_crit`), under combined
  compression and shear (`calc_Nxx_crit_combined_shear`,
  `calc_Nxx_crit_combined_shear_full`), and post-buckling effective width
  (`calc_beff`).
- `Laminate.make_smeared()`, for a laminate with smeared properties.

## 0.6.5 (2024-02-29)

- Fixed the gradient of the `E` matrix with respect to the lamination
  parameters.
- Fixed the Linux wheels, now built with `cibuildwheel`, and the
  `python_requires` flag.
- Python 3.12 support.
- Fixed an array assignment error in the lamination parameter gradients.
- Improved documentation, with links to the source code.

## 0.5.25 (2023-03-24)

- Removed the dependency on the NumPy C API, so the compiled package works
  with all NumPy versions; relaxed the NumPy version constraint.
- Improved OpenMP build, fixed dtypes.

## 0.5.17 (2022-12-11)

- Installation with `pip install`; Linux wheels fixed.
- Python 3.11 support; Python 3.7 dropped.
- Travis CI removed in favour of GitHub Actions.

## 0.5.4 (2022-06-24)

- macOS wheels and a single source distribution per release.
- Analytical gradients of `Aij`, `Bij`, `Dij` with respect to the lamination
  parameters (`laminate_LP_gradients`).
- Fixed the lamination parameters, which were not read correctly.
- Fixed the 3D stress state input of `read_laminaprop`.
- Tests on macOS; small performance improvements and docstrings.

## 0.4.22 (2022-03-10)

- The equivalent properties are calculated by default in the `utils`
  functions.
- Double-double and n-double laminates return a reference stacking sequence.
- Python 3.10 support; Python 3.6 dropped.
- Citation file and Zenodo DOI.
- License changed from the 2-clause to the 3-clause BSD.

## 0.4.12 (2021-08-04)

- n-double laminate (`n_double_plate`).
- Double-double (DD) laminate (`double_double_plate`).
- Plane stress is the default behavior.
- Simplified `Laminate`, which no longer needs `rebuild()`, and simplified
  trace normalization.
- **Breaking:** removed support for the legacy input.
- Documentation renders equations with MathJax; improved release
  descriptions.
- All classes and functions can be `cimport`ed by other Cython packages, with
  declarations in `core.pxd`.
- Full compatibility with finite element mass matrices of plates and shells,
  supporting laminated plates with plies of different densities.
- `thetadeg` replaces `theta`.
- Coverage with Codecov.

## 0.3.7 (2021-06-26)

- `calc_equivalent_properties`.
- `rho` argument in `isotropic_plate`.
- New mass integral `intrhoz`, needed for plates with `B != 0`.
- **Breaking:** rewritten in Cython (`composites.core`), with large
  simplifications and speed improvements; not compatible with previous
  versions.
- Integrated mass properties.
- Documentation deployed automatically; automatic releases.


## 0.2.5 (2020-01-20)

- `read_isotropic`, clean-up and documentation.
- CI and publishing with GitHub Actions.
- Independent of previous packages.
- From the sub-module of meshless 0.1.19.
- Syntax improvements.


## 0.1.1 (2018-06-14)

- Initial release, from the sub-module of compmech 0.7.2.
- Removed the pyNastran dependency.
