# Changelog

## Unreleased

### Behaviour change: plane-stress invariants of `MatLamina` for any `laminaprop`

- `MatLamina.rebuild()` computes `q11`, `q12`, `q21`, `q22` and the
  invariants `u1` to `u5` with the plane-stress reduced stiffnesses,
  `sigma33 = 0` (Jones, 1999, Eq. (2.66); Reddy, 2003, Eq. (2.4.4)), also
  when `e3`, `nu13` and `nu23` are given, i.e. for a 9-entry `laminaprop`,
  and sets `q13`, `q23`, `q31`, `q32` and `q33` to zero. Before, with
  `e3 != 0` they were the three-dimensional stiffnesses, such that
  `laminate_from_LaminationParameters()`,
  `laminate_from_lamination_parameters()`, `GradABD.calc_LP_grad()`,
  `trace_normalize_plane_stress()`, `n_double_laminate()` and
  `double_double_laminate()` gave stiffer `A`, `B`, `D` than
  `laminated_plate()` for the same material and layup, e.g. `(lambda +
  2 mu) h` instead of `E h/(1 - nu^2)` for an isotropic material, 32% higher
  for `nu = 0.33`. Both routes now agree to round-off for any layup. The
  results with a 2- or 6-entry `laminaprop` (`e3 = 0`) and those of
  `laminated_plate()` are unchanged, bit for bit.
- New option `MatLamina.rebuild(plane_stress=False)`, which keeps the
  three-dimensional constants and invariants, e.g. for continuum shell
  elements; the choice is stored in the new read-only attribute
  `MatLamina.plane_stress`, `True` by default, which is pickled and, in
  JSON, restored from `q33` without a change of the format version.
  `trace_normalize_plane_stress()` raises `ValueError` for
  `plane_stress=False`, since the trace of Melo et al. (2017), Eq. (1), is
  that of the plane-stress stiffness matrix. The stiffness constants `c11`
  to `c66` and `get_constitutive_matrix()` remain three-dimensional.
- Tests: new `tests/test_plane_stress_invariants.py`, where
  `laminated_plate()` and `laminate_from_LaminationParameters()` are
  compared, for an isotropic and an orthotropic material with nine
  constants and four layups, including unbalanced and unsymmetric ones,
  with the classical lamination theory to `rtol = 1e-12`, together with the
  trace normalization and `GradABD`. `test_lampar_tri_axial` and
  `test_matlamina_constitutive_matrix_inverse_compliance` check the
  plane-stress default and the three-dimensional option. In
  `tests/test_DD_laminate.py`, the trace is checked against Table 1 of Melo
  et al. (2017) to half a unit of its last printed digit, and the reference
  invariants use `Q66` and `Q44`, `Q55` in their places (the values are
  unchanged, all shear moduli of the test being equal).

### Fix: lamination parameters with an offset

- `laminate_from_LaminationParameters()` and
  `laminate_from_lamination_parameters()` take a new argument `offset`,
  default `0.`, and `GradABD.calc_LP_grad()` too. Since 0.3.0,
  `Laminate.calc_lamination_parameters()` integrates `xiB` and `xiD` about
  the reference surface, i.e. including `offset`, while
  `laminate_from_LaminationParameters()` assumed the mid-surface: it left
  out the terms `offset*A` of `B` and `offset**2*A` of `D` of the invariant
  `u1`, `u4`, `u5` parts, such that the round trip of a laminate with an
  offset gave wrong `B` and `D`, e.g. `B11` 4.8 times too small for a
  `[0/90]s` laminate with `offset = 0.2 mm` and `0.125 mm` plies. Given the
  same `offset`, the stiffnesses now satisfy `B = B0 + e A` and `D = D0 +
  2 e B0 + e^2 A` (Vasiliev and Morozov, 2001, Eq. (5.28)), and reproduce
  `laminated_plate()` to round-off; the created laminate stores `offset`.
  `A` and `Abar44`, `Abar45`, `Abar55` do not depend on the offset. The
  derivatives of `GradABD` with respect to the thickness are taken at
  constant lamination parameters and offset. Results with `offset = 0` are
  unchanged.
- Tests: `test_lamination_parameters_offset` in
  `tests/test_plane_stress_invariants.py`, with offsets of +0.2 mm,
  -0.2 mm and one ply thickness, for the materials and layups above,
  against the classical lamination theory to `rtol = 1e-12`, the shift of
  Vasiliev and Morozov, the transverse shear stiffness and central
  differences of `GradABD`.

## 0.10.1 (2026-10-26)

### Python versions

- Python 3.9 is no longer supported, it reached its end of life in October 2025: Python 3.10 or newer is required, and the tests and wheels cover Python 3.10 to 3.14.

### Packaging

- Packaging: the metadata, with the version, in `pyproject.toml` (PEP 621), built with `setuptools>=77` and Cython; `setup.py` keeps only the Cython extension with the compiler flags, including those of the WebAssembly build, and writes `composites/version.py` from `pyproject.toml`. The license is given as the SPDX expression `BSD-3-Clause` instead of a classifier, wheel is no longer a build requirement, and README.md, LICENSE and version.py are no longer installed in the root of the Python environment (`data_files`).

### Mass integrals of the TSDT

- `Laminate.intrhoz3`, `intrhoz4` and `intrhoz6`: the moments of the density
  `int rho z^n dz` for `n = 3, 4, 6`, computed by `calc_constitutive_matrix`
  next to `intrho`, `intrhoz`, `intrhoz2`, with `z` measured from the
  reference surface. Together they are the moments `I_0` to `I_4` and `I_6`
  required by the mass matrices of the third-order shear deformation theory.
  They are zero for laminates without plies, e.g. created from lamination
  parameters.

### Breaking: JSON format version 2

The JSON files now have `format_version` 2, which stores `intrhoz3`,
`intrhoz4` and `intrhoz6`. `from_dict`, `from_json` and `load_json` read only
the current format version and raise `ValueError` for an older or a newer
one, such that the files of format 1, saved by version 0.9.21, must be
created again. Every key of the format is required: a key missing from
`'data'`, at any level, raises `ValueError` naming the missing keys and the
class, instead of leaving the attribute at its default value.

### Breaking: sign of the shear load of the Kassapoglou combined buckling

`composites.kassapoglou.calc_Nxx_crit_combined_shear_full(Nxy, ...)` now
takes `Nxy` positive in tension, the `N_xy` of the energy of Kassapoglou's
Eq. 6.32, instead of its negative, the `k*N_o` of Eqs. 6.32a and 6.33. The
result changes only when `D16` or `D26` is not zero: a caller that passed
`k*N_o` must now pass `-k*N_o` to keep its results. The returned load is
unchanged, the compression load `N_o = -N_x`, positive in compression.

### Fixes

- `GradABD.calc_LP_grad` left the derivatives `dA66/dxiA3`, `dB66/dxiB3` and
  `dD66/dxiD3` at zero, instead of `-h*u3`, `-h^2/4*u3` and `-h^3/12*u3`.
- `MatLamina.rebuild`: the stiffness constant `c11` used `nu23*nu23` instead
  of `nu23*nu32`, differing by 4.3% from the inverse of the compliance matrix
  for `e2 != e3`, and `q31` used `nu12*nu32` instead of `nu12*nu23`,
  differing by 4.65%. `get_constitutive_matrix()` is now the inverse of the
  orthotropic compliance matrix, and `q31 = q13`.
- `Laminate.make_symmetric()` set to zero `F16`, `F26`, `H16`, `H26`, `D45`
  and `F45`, which are integrals of even powers of `z` and do not vanish in a
  symmetric laminate, e.g. `F16`, `H16` and `D45` of a `[30/-45/60]s`
  laminate. Only the integrals of odd powers, `B11` to `B66` and `E11` to
  `E66`, are now set to zero; the other attributes keep the values of
  `calc_constitutive_matrix()`. This changes the TSDT stiffnesses and the
  transverse shear terms `D45`, `F45` after `make_symmetric()`.
- `Laminate.make_orthotropic()` set to zero only `A16`, `A26`, `B16`, `B26`,
  `D16` and `D26`. It now also sets to zero the other terms that vanish in a
  laminate with all plies at 0 or 90 degrees, the integrals of `q16L`,
  `q26L` and `q45L`: `E16`, `E26`, `F16`, `F26`, `H16`, `H26`, `D45`, `F45`,
  `Abar45`, `Abarbar45` and `A45`. This changes the TSDT stiffnesses and the
  transverse shear coupling terms after `make_orthotropic()`. `Abarbar45`
  is kept when it is `nan`, i.e. when a ply has a singular transverse shear
  stiffness matrix, as by `make_balanced()`.
- `make_orthotropic_LP()` left `xiAts2` unchanged, which gives `Abar45` in
  `laminate_from_LaminationParameters()` and vanishes when all plies are at
  0 or 90 degrees; it is now set to zero, as `Abar45` by
  `Laminate.make_orthotropic()`.
- `Laminate.make_balanced()` set to zero `B16` and `B26`, which do not
  vanish in a balanced laminate in general, e.g. in an antisymmetric
  `[+45/-45]` laminate, and left `Abar45`, `Abarbar45` and `A45`, which
  vanish in every balanced laminate. It now sets to zero only the terms that
  vanish for any position of the `+theta` and `-theta` plies, the integrals
  of `z^0` of `q16L`, `q26L` and `q45L`: `A16`, `A26`, `Abar45`,
  `Abarbar45` (unless `nan`), and `A45` except with
  `shear_correction='rohwer'`, whose `A45` does not vanish in general, e.g.
  for a `[+30/-30/0]s` laminate. `B16`, `B26` and the other integrals of
  `z^n`, `n >= 1` (`D16` to `H26`, `D45`, `F45`) keep the values of
  `calc_constitutive_matrix()`.
- `make_balanced_LP()` left `xiAts2` unchanged, which gives `Abar45` in
  `laminate_from_LaminationParameters()` and vanishes in every balanced
  laminate; it is now set to zero, as `Abar45` by `Laminate.make_balanced()`.

### Deprecated

- `Lamina.get_transf_matrix_displ_to_laminate()`: use
  `Lamina.get_transf_matrix_displ_to_lamina()`. The matrix maps the
  displacements from the laminate axes to the material axes of the ply,
  contrary to what the old name said; the matrix itself is unchanged. The old
  name still works, but it emits a `DeprecationWarning`.
- The `grad_x` and `grad_y` arguments of
  `Laminate.calc_equilibrium_transverse_shear()` and
  `Laminate.calc_aposteriori_energy()`: use `eps_x` and `eps_y`, the
  derivatives of the generalized strains with respect to `x` and `y`. The
  positional order is unchanged, `(eps_x, eps_y, z=None)` and
  `(eps_x, eps_y)`, such that positional calls are not affected. The old
  keywords still work, but they emit a `DeprecationWarning`; giving both
  `eps_x` and `grad_x`, or `eps_y` and `grad_y`, raises `TypeError`.

### Clean-up

- Partial derivatives are named with an underscore before the variable, e.g.
  `eps_x` for the derivative of `eps` with respect to `x`: the internal
  derivatives of the generalized strains of
  `calc_equilibrium_transverse_shear` (`gx`, `gy`), and the derivative
  `expr_alpha` (`dexpr_dalpha`) of the Newton-Raphson iterations of
  `kassapoglou.calc_Nxy_crit` and of its generator in `theory/`. The results
  are unchanged.

### Documentation

- Pages on the laminate constitutive matrices of the CLPT, FSDT and TSDT, the
  thickness integrals for mass matrices, the conventions for defining a
  laminate, the choice of the transverse shear correction, saving and
  loading in JSON, installing and building, and the changelog.
- Usage examples extracted from the tests, one API page per module, and two
  tutorials: plate properties for the CLPT, FSDT and TSDT, and transverse
  shear stresses through the thickness.
- Reddy's book is cited as Reddy (2003), https://doi.org/10.1201/b12409.
- The Notes of the integrals of the equilibrium shear stresses use the
  derivatives `\varepsilon^{(0)}_{,x}`, `\varepsilon^{(1)}_{,x}`, as those of
  `calc_equilibrium_transverse_shear`, instead of the symbols `g_x^{(0)}`,
  `g_x^{(1)}`.
- The transverse shear tutorial uses `eps_x` and `eps_y`, and the
  sensitivities of the gradient-based optimisation tutorials follow the same
  naming, e.g. `lambda_b_h0`, `D_h0`, `D11_xiD1` instead of `dlambda_b_dh0`,
  `dD_dh0`, `dD11_dxiD1`.

## 0.9.21 (2026-10-08)

### Saving and loading in JSON

JSON replaces pickle as the way to save objects to files. The new module
`composites.json_io`, also available from `composites`, provides:

- `save_json(obj, fname)` and `load_json(fname)`: save to and load from a
  JSON file, given by its name or as a file object.
- `to_json(obj)` and `from_json(s)`: the same with a JSON string.
- `to_dict(obj)` and `from_dict(d)`: the same with a dictionary of
  JSON-compatible types, e.g. to embed the objects in other JSON documents.

They support `Laminate`, `Lamina`, `MatLamina`, `LaminationParameters` and
`GradABD`. Unlike pickle, the JSON files are readable, do not depend on the
Python or composites versions, and are safe to load from untrusted sources.
All public attributes are stored, such that the loaded object is identical to
the saved one, and a `MatLamina` shared by several plies is stored once and
still shared after loading. `nan` and `inf` are stored as the strings
`"NaN"`, `"Infinity"` and `"-Infinity"`, such that the output is strict JSON.
The files carry a `format_version`, and loading a file from a newer format
raises `ValueError`.

Pickling `Laminate` and `GradABD` is still supported for `copy.deepcopy` and
for passing objects between processes, e.g. with `multiprocessing`, but it is
no longer recommended to save objects to files.

### Pyodide (WebAssembly)

- A WebAssembly wheel for Pyodide 314 (CPython 3.14), tagged
  `pyemscripten_2026_0_wasm32` (PEP 783), is built with cibuildwheel, tested
  in Pyodide on every push and published on PyPI, such that
  `micropip.install('composites')` works in the browser.
- `setup.py` detects the Emscripten cross-compilation of `pyodide-build` and
  leaves out the GCC runtime link flags.
- NumPy is no longer a build requirement, since no module uses its C API; it
  is still a runtime requirement.

### Maintenance

- Python 3.8 support dropped, `python_requires='>=3.9'`.
- GitHub Actions updated: `actions/checkout@v7`, `actions/setup-python@v7`,
  `codecov/codecov-action@v7`, `softprops/action-gh-release@v3` and
  cibuildwheel 4.3.0.
- `setup.py` imports `Extension` from `setuptools` instead of the removed
  `distutils`.

## 0.9.12 (2026-10-02)

### New shear correction methods

Four values of `shear_correction` were added to
`Laminate.calc_transverse_shear_stiffness()`, `laminated_plate` and
`isotropic_plate`. All of them return scalar factors, applied as for
`'vlachoutsis'`: `A55 = k13*Abar55`, `A44 = k23*Abar44` and
`A45 = (k13 + k23)/2*Abar45`.

- `'whitney'`: Whitney (1973), Eqs. (3)-(7), extension of Chow (1971) to
  unsymmetric orthotropic laminates, identical to `'vlachoutsis'`; it
  reproduces the factors printed by Whitney (1973).
- `'chow'`: Chow (1971), Eqs. (11)-(12), energy equivalence with the
  one-dimensional equilibrium distribution of a symmetric laminate. It raises
  `ValueError` for unsymmetric laminates, and it is identical to
  `'vlachoutsis'` for symmetric ones.
- `'birman_bert'`: average shear strain factor of Birman and Bert (2002),
  Eq. (13), generalized to unsymmetric laminates with the direction-wise
  neutral surface of `'vlachoutsis'`. It gives 1 for a homogeneous plate.
- `'thickness_shear'`: dynamic factor of Yang, Norris and Stavsky (1966),
  extending the criterion of Mindlin (1951) to laminates: the frequency of
  the first thickness-shear mode of the FSDT is matched with the exact one of
  the laminate, computed with a transfer matrix through the plies. It gives
  `pi^2/12` for a homogeneous plate and requires positive ply densities.

### A posteriori transverse shear stresses and energies

- `Laminate.calc_equilibrium_transverse_shear(grad_x, grad_y, z=None)`: the
  transverse shear stresses recovered by integrating the 3D equilibrium
  equations with given gradients of the generalized strains, e.g. those of an
  FSDT solution, as in the a posteriori approach of Noor and Peters (1989).
- `Laminate.calc_aposteriori_energy(grad_x, grad_y)`: the transverse shear
  strain energies of the recovered stresses per direction, and their
  resultants, for the energy equivalence of Noor and Peters (1989),
  Eqs. (8)-(9).

### Fixes

- All `shear_correction` methods, `calc_transverse_shear_stress` and the a
  posteriori helpers are evaluated with the heights measured from the
  mid-surface. The results, which do not depend on `offset`, are now identical
  for any `offset`; before, the round-off grew with `(offset/h)^2`, e.g. a
  relative error of 2e-3 in `Ats` of `'rohwer'` for `offset/h = 1e4`, and
  wrong values without any error for larger offsets.
- `calc_equilibrium_transverse_shear` raised `IndexError` on Linux with
  Python 3.12 (negative index with `wraparound=False`).
- Pickling a `GradABD` created without `__init__`, e.g. by `__new__`, crashed
  with a segmentation fault; its arrays are now allocated in `__cinit__`.

### Deprecated

The remaining `Atrans` names were renamed to `Ats`, and `Dtrans`, `Ftrans`
to `Dts`, `Fts`. The old names still work, but they emit a
`DeprecationWarning`:

- `Laminate.Dtrans`, `Laminate.Ftrans`: use `Laminate.Dts`, `Laminate.Fts`.
- `LaminationParameters.xiAtrans1`, `xiAtrans2`: use `xiAts1`, `xiAts2`.
- The `xiAtrans1`, `xiAtrans2` arguments of
  `laminate_from_lamination_parameters`: use `xiAts1`, `xiAts2`.
- `GradABD.gradAtransij`: use `gradAtsij`. Pickles that contain
  `gradAtransij` still load.
- `Laminate.get_Atrans()`, `get_Dtrans()`, `get_Ftrans()` (`cdef`): use
  `get_Ats()`, `get_Dts()`, `get_Fts()`.

Packages that `cimport composites` and access `xiAtrans1`, `xiAtrans2` or
`gradAtransij` from Cython must switch to the new names and be recompiled.
Pickled `LaminationParameters` objects from older versions cannot be loaded.

## 0.9.2 (2026-09-25)

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

### Performance: Cython directives and compiler flags (0.9.2)

- `composites/core.pyx` is compiled with `initializedcheck=False`. Cython no
  longer checks, at each memoryview access, that the memoryview has been
  assigned. Reading an unassigned memoryview attribute from Python still
  raises `AttributeError`, because the property getters keep their own check.
- `setup.py` passed `linetrace=True` to every build. Without `CYTHON_TRACE`
  the line tracing itself stays compiled out, but the directive also
  generates profiling hooks that are active by default and run at every call
  of a Python-visible function. The directive is now set only for a coverage
  build, see below.
- Compiler flags: `-O3` and `-fno-math-errno` for GCC and Clang, and an
  explicit `/O2` for MSVC, which setuptools already combines with `/GL`. The
  OpenMP flags are removed, since the module does not use `prange`. Flags
  that change floating-point results, such as `/fp:fast` or `-ffast-math`,
  and flags that tie a wheel to the CPU that built it, such as
  `-march=native`, are deliberately not used.

The results are bit-identical to the previous build, checked on the `ABD`,
`E` and `Ats` matrices of laminates with 1, 8 and 33 plies, and the whole
test suite passes.

Measured on an AMD x86-64 laptop CPU (family 25, model 117), Windows 11,
Python 3.13, NumPy 2.4, Cython 3.2, MSVC, with the process pinned to one core
at high priority, minimum of 3 interleaved runs of 5 repetitions each. The
GCC and Clang flags were not measured. Cost per call from Python
(`benchmarks/bench_build_flags.py`), in microseconds, for a CFRP laminate;
each cell gives the cost before and after this change, with the relative
difference:

| Function | 8 plies | 32 plies | 128 plies |
| - | - | - | - |
| `laminated_plate` | 102.26 / 94.49 (-8%) | 232.65 / 215.25 (-7%) | 750.23 / 675.84 (-10%) |
| `Laminate.calc_constitutive_matrix` | 70.71 / 67.76 (-4%) | 160.56 / 152.04 (-5%) | 507.97 / 484.29 (-5%) |
| `Laminate.calc_transverse_shear_stiffness` | 35.87 / 35.44 (-1%) | 39.10 / 39.62 (+1%) | 48.59 / 47.52 (-2%) |
| `Laminate.calc_lamination_parameters` | 8.83 / 7.44 (-16%) | 32.92 / 30.39 (-8%) | 136.57 / 118.95 (-13%) |
| `laminate_from_LaminationParameters` | 2.63 / 2.04 (-22%) | 2.63 / 1.99 (-24%) | 2.65 / 2.05 (-23%) |
| `GradABD.calc_LP_grad` | 1.95 / 1.91 (-2%) | 1.94 / 1.91 (-2%) | 1.94 / 1.91 (-1%) |
| `n_double_laminate` | 10.31 / 9.78 (-5%) | 10.37 / 9.86 (-5%) | 10.23 / 9.91 (-3%) |

Differences within about 3% are measurement noise. Creating a laminate is 7
to 10% faster, and the functions dominated by memoryview access, such as
`laminate_from_LaminationParameters` and `calc_lamination_parameters`, up to
24% faster.

### Coverage build (0.9.2)

`setup.py` now recognises a coverage build, requested either with
`CYTHON_TRACE_NOGIL` in the environment or with `--define CYTHON_TRACE_NOGIL`,
as in `.github/workflows/coverage.yml`, and compiles it unoptimized with line
tracing. Since Python 3.12 Cython traces through `sys.monitoring` by default,
which the `Cython.Coverage` plugin does not follow, so the coverage build now
requests the legacy tracing with `CYTHON_USE_SYS_MONITORING=0`. Measured with
the command of the workflow on Python 3.13: the report used to list only the
Python modules, 188 statements at 100%, without `core.pyx`; it now includes
`core.pyx`, 833 statements at 99%. Switching between a coverage build and a
normal one regenerates `core.cpp`, which cythonize would otherwise reuse from
the other mode.


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
