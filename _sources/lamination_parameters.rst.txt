Lamination parameters and double-double laminates
=================================================

The ply stiffnesses are linear in `\cos 2\theta`, `\sin 2\theta`, `\cos
4\theta` and `\sin 4\theta`, with the material invariants as coefficients,
Eq. :eq:`lam_qbar_invariants`. For a laminate of a single material, the
thickness integrals of the trigonometric terms, the lamination parameters,
therefore describe the stacking sequence completely as far as the stiffness
is concerned: a few continuous variables replace the ply angles, and the
stiffnesses are linear in them. This is the basis of the design of laminates
with lamination parameters, of their gradients for gradient-based
optimization, and of the double-double and N-double laminates, whose
stiffnesses are averages over a few angles.


Definition
----------

With `\bar{z} = z/h`, `z` measured from the reference surface, and
`\theta(\bar{z})` the ply angle, :meth:`.Laminate.calc_lamination_parameters`
computes, exactly ply by ply, the 14 parameters stored in a
:class:`.LaminationParameters`:

.. math::
    :label: lp_definition

    \{\xi_{A1}, \xi_{A2}, \xi_{A3}, \xi_{A4}\} = \int \{\cos 2\theta,
    \sin 2\theta, \cos 4\theta, \sin 4\theta\} \, d\bar{z} \\
    \{\xi_{B1}, \xi_{B2}, \xi_{B3}, \xi_{B4}\} = 4 \int \{\cos
    2\theta, \sin 2\theta, \cos 4\theta, \sin 4\theta\} \, \bar{z} \,
    d\bar{z} \\
    \{\xi_{D1}, \xi_{D2}, \xi_{D3}, \xi_{D4}\} = 12 \int \{\cos
    2\theta, \sin 2\theta, \cos 4\theta, \sin 4\theta\} \, \bar{z}^2 \,
    d\bar{z} \\
    \{\xi_{A_{ts}1}, \xi_{A_{ts}2}\} = \int \{\cos 2\theta, \sin
    2\theta\} \, d\bar{z}

the integrals running over the thickness, `\bar{e} - 1/2 \le \bar{z} \le
\bar{e} + 1/2`, `\bar{e} = d/h` with `d` = ``offset``. The `\xi_{Ai}` and
`\xi_{Di}` are those of Krupa et al. [1]_, who integrate about the
mid-surface; the `\xi_{Bi}` and the transverse shear parameters follow the
same pattern, and `\xi_{A_{ts}i} = \xi_{Ai}`. The factors 1, 4 and 12 make
the integrals of a constant equal to 1, `4\bar{e}` and `1 + 12\bar{e}^2`,
i.e. 1, 0 and 1 about the mid-surface, where `|\xi_{Ai}| \le 1` and
`|\xi_{Di}| \le 1`. The `\xi_{Ai}` do not depend on the reference surface,
whereas

.. math::
    :label: lp_shift

    \xi_{Bi} = \xi_{Bi}^0 + 4 \bar{e} \xi_{Ai}, \qquad
    \xi_{Di} = \xi_{Di}^0 + 6 \bar{e} \xi_{Bi}^0 + 12 \bar{e}^2 \xi_{Ai}

with the superscript 0 for the parameters about the mid-surface, the
counterpart of Eq. :eq:`cm_shift`. Parameters with an offset may therefore
lie outside the feasible region of the parameters about the mid-surface.


Stiffness from lamination parameters
------------------------------------

Integrating Eq. :eq:`lam_qbar_invariants` through the thickness with `\int
dz = h \int d\bar{z}`, `\int z \, dz = h^2/4 \cdot 4 \int \bar{z} \,
d\bar{z}` and `\int z^2 dz = h^3/12 \cdot 12 \int \bar{z}^2 d\bar{z}` gives
the stiffnesses of :func:`.laminate_from_LaminationParameters`:

.. math::
    :label: lp_stiffness

    A_{11} = h (u_1 + u_2 \xi_{A1} + u_3 \xi_{A3}), \quad
    A_{12} = h (u_4 - u_3 \xi_{A3}), \quad
    A_{22} = h (u_1 - u_2 \xi_{A1} + u_3 \xi_{A3}) \\
    A_{16} = h \left( \frac{u_2}{2} \xi_{A2} + u_3 \xi_{A4} \right), \quad
    A_{26} = h \left( \frac{u_2}{2} \xi_{A2} - u_3 \xi_{A4} \right), \quad
    A_{66} = h (u_5 - u_3 \xi_{A3}) \\
    B_{11} = \frac{h^2}{4} (4\bar{e} u_1 + u_2 \xi_{B1} + u_3 \xi_{B3}),
    \quad \dots, \qquad
    D_{11} = \frac{h^3}{12} ((1 + 12\bar{e}^2) u_1 + u_2 \xi_{D1} + u_3
    \xi_{D3}), \quad \dots \\
    \bar{A}_{44} = h (u_6 + u_7 \xi_{A_{ts}1}), \quad
    \bar{A}_{45} = -h u_7 \xi_{A_{ts}2}, \quad
    \bar{A}_{55} = h (u_6 - u_7 \xi_{A_{ts}1})

where `B_{ij}` and `D_{ij}` follow the pattern of `A_{ij}` with the
parameters `\xi_{Bi}` and `\xi_{Di}`, and with the constant terms `u_1`,
`u_4`, `u_5` multiplied by `4\bar{e}` and `1 + 12\bar{e}^2`, respectively.
:func:`.laminate_from_lamination_parameters` does the same from the 14
values given as arguments. With the default plane-stress invariants of
:meth:`.MatLamina.rebuild`, and the same ``offset``, these are exactly the
stiffnesses of :func:`.laminated_plate` for the layup that gave the
parameters (Jones [2]_, Eq. (4.24)). The parameters `\xi_{A_{ts}i}` multiply
`u_7 = (G_{23} - G_{13})/2`, and only matter when `G_{13} \ne G_{23}`.

The lamination parameters do not describe the through-thickness
distribution of the plies, such that a laminate created from them has no
plies and only `A`, `B`, `D` and `\bar{A}_{ts}`: `E`, `F`, `H`, `D_{ts}`,
`F_{ts}` and the mass integrals are zero, `\bar{\bar{A}}_{ts}` is ``nan``, no
shear correction can be computed, ``shear_correction`` is ``None`` and
``Ats = Abar_ts``, see :doc:`shear_correction`. Such a laminate serves the
CLPT and, with a user factor on ``Abar_ts``, the FSDT.


Special laminates
-----------------

The functions :func:`.make_balanced_LP`, :func:`.make_symmetric_LP` and
:func:`.make_orthotropic_LP` set to zero, in place, the parameters that
vanish in the respective laminates, consistently with the methods of
:doc:`constitutive_matrices`:

- balanced, a ply at `-\theta` for every ply at `+\theta`: the integrals
  of `\sin 2\theta` and `\sin 4\theta` without weight, `\xi_{A2}`,
  `\xi_{A4}`, `\xi_{A_{ts}2}`, whereas `\xi_{B2}`, `\xi_{B4}`, `\xi_{D2}`,
  `\xi_{D4}` depend on the positions of the plies;
- symmetric about the mid-surface: `\xi_{Bi} = 0`, which holds for the
  parameters about the mid-surface only; with an offset `\xi_{Bi} = 4\bar{e}
  \xi_{Ai}` in a symmetric laminate, Eq. :eq:`lp_shift`;
- orthotropic, all plies at 0 or 90 degrees: all the sine parameters. The
  `\xi_{D2}` and `\xi_{D4}`, related to the bending-twisting coupling,
  become small in balanced laminates with many plies.


Gradients
---------

At constant thickness the stiffnesses of Eq. :eq:`lp_stiffness` are linear in
the lamination parameters, with gradients that depend only on the thickness
and on the material, and at constant lamination parameters and offset they
are polynomials of `h`. :meth:`.GradABD.calc_LP_grad` stores both in the
attributes of :class:`.GradABD`: ``gradAij``, ``gradBij`` and ``gradDij``,
of shape ``(6, 5)``, rows ``11, 12, 16, 22, 26, 66``, column 0 the
derivative with respect to `h` and columns 1 to 4 those with respect to the
four parameters of the matrix, and ``gradAtsij``, of shape ``(3, 3)``, rows
``Abar44, Abar45, Abar55``, columns `h`, `\xi_{A_{ts}1}`, `\xi_{A_{ts}2}`.
With `U` the rows of the in-plane terms of
:meth:`.MatLamina.get_invariant_matrix`, `A = h \, U \{1, \xi_A\}`, `B = h d
\, U \{1, 0\} + h^2/4 \, U \{0, \xi_B\}` and `D = (h^3/12 + h d^2) \, U \{1,
0\} + h^3/12 \, U \{0, \xi_D\}`, such that

.. math::
    :label: lp_gradients

    \frac{\partial A}{\partial h} = U \{1, \xi_A\}, \quad
    \frac{\partial B}{\partial h} = d \, U \{1, 0\} + \frac{h}{2} U \{0,
    \xi_B\}, \quad
    \frac{\partial D}{\partial h} = \left( \frac{h^2}{4} + d^2 \right) U
    \{1, 0\} + \frac{h^2}{4} U \{0, \xi_D\}

and the derivatives with respect to the parameters are the columns of `U`
times `h`, `h^2/4` and `h^3/12`; analogously for `\bar{A}_{ts}`. The
transverse shear gradients are those of the uncorrected `\bar{A}_{ts}`, the
``Ats`` of the lamination-parameter route.


Double-double and N-double laminates
------------------------------------

A double-double laminate is the repetition of the sub-laminate
`[\pm\phi, \pm\psi]`, and an N-double laminate that of `[\pm\phi_1,
\pm\phi_2, \cdots, \pm\phi_n]`, all plies of the same material and
thickness (Shrivastava et al. [3]_). In the limit of many repetitions the
laminate is homogenized: `B = 0`, and the stiffnesses are the averages of
the ply stiffnesses over the `2n` angles. They are written with the trace
`\mathrm{tr}` and the trace-normalized invariants `u_i^*` of
:meth:`.MatLamina.trace_normalize_plane_stress` (Melo et al. [4]_):

.. math::
    :label: lp_nd

    A_{11}^* = u_1^* + \frac{1}{n} \sum_{i=1}^{n} \left( u_2^* \cos
    2\phi_i + u_3^* \cos 4\phi_i \right), \quad
    A_{22}^* = u_1^* + \frac{1}{n} \sum_{i=1}^{n} \left( -u_2^* \cos
    2\phi_i + u_3^* \cos 4\phi_i \right) \\
    A_{12}^* = u_4^* - \frac{1}{n} \sum_{i=1}^{n} u_3^* \cos 4\phi_i,
    \quad
    A_{66}^* = u_5^* - \frac{1}{n} \sum_{i=1}^{n} u_3^* \cos 4\phi_i,
    \quad A_{16}^* = A_{26}^* = 0 \\
    A_{ij} = \mathrm{tr} \, A_{ij}^* h, \qquad D_{ij} = \mathrm{tr} \,
    A_{ij}^* h^3/12

the bending stiffnesses being those of the smeared laminate of
:meth:`.Laminate.make_smeared`. For the double-double laminate, with
`\alpha = (\phi + \psi)/2` and `\beta = (\phi - \psi)/2`, the averages are
`\cos 2\alpha \cos 2\beta` and `\cos 4\alpha \cos 4\beta`, which
:func:`.double_double_laminate` uses. The three functions give the same
stiffnesses:

- :func:`.double_double_laminate`, from ``laminaprop`` and the two angles;
- :func:`composites.utils.n_double_laminate`, from ``laminaprop`` and a list
  of angles;
- :func:`composites.core.n_double_laminate`, the compiled version, from a
  :class:`.MatLamina`, which it trace-normalizes in place, and the first
  ``n`` angles of an array.

The laminates have no plies: only `A` and `D` are computed, all the other
stiffnesses, the transverse shear ones included, are zero, as are the mass
integrals, and the equivalent properties of
:meth:`.Laminate.calc_equivalent_properties` are computed. They serve the
CLPT, with the homogenization assumption; an FSDT or TSDT model needs the
stiffnesses of a stacking sequence, e.g. from :func:`.laminated_plate`.


Verification
------------

``tests/test_core.py`` and ``tests/test_plane_stress_invariants.py`` check
that the lamination parameters of a laminate give back its `A`, `B`, `D`
and `\bar{A}_{ts}` through :func:`.laminate_from_LaminationParameters`, also
with an offset and with three-dimensional constants, against a
classical-lamination-theory implementation of the tests and the shift of
Vasiliev and Morozov [5]_, Eq. (5.28); the special-laminate functions; and
the gradients of :class:`.GradABD` against central differences.
``tests/test_DD_laminate.py``, ``tests/test_ND_laminate.py`` and
``tests/test_ND_laminate_core.py`` check the trace against Melo et al.
[4]_, Table 1, and the `A`, `D` of the double-double and N-double laminates
against those of :func:`.laminated_plate` with 32 repetitions of the
sub-laminate.


References
----------

.. [1] Krupa, E. P., Cooper, J. E., Pirrera, A. and Nangia, R., Improved
   aerostructural performance via aeroservoelastic tailoring of a composite
   wing, The Aeronautical Journal, 122(1255), 1442-1474, 2018.
   https://doi.org/10.1017/aer.2018.66
.. [2] Jones, R. M., Mechanics of Composite Materials, Second Edition,
   Taylor & Francis, Philadelphia, 1999.
.. [3] Shrivastava, S., Sharma, N., Tsai, S. W. and Mohite, P. M., D and
   DD-drop layup optimization of aircraft wing panels under multi-load case
   design environment, Composite Structures, 248, 112518, 2020.
   https://doi.org/10.1016/j.compstruct.2020.112518
.. [4] Melo, J. D. D., Bi, J. and Tsai, S. W., A Novel Invariant-Based
   Design Approach to Carbon Fiber Reinforced Laminates, Composite
   Structures, 159, 44-52, 2017.
   https://doi.org/10.1016/j.compstruct.2016.09.055
.. [5] Vasiliev, V. V. and Morozov, E. V., Mechanics and Analysis of
   Composite Materials, Elsevier, Oxford, 2001.
