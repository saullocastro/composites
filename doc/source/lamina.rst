Ply constitutive relations and invariants
=========================================

Every stiffness of a laminate is a thickness integral of the stiffnesses of
its plies. This chapter describes how ``composites`` obtains them: from the
engineering constants of the material, :class:`.MatLamina`, to the
plane-stress stiffnesses of a ply rotated to the laminate axes,
:class:`.Lamina`, and to the material invariants on which the lamination
parameters of :doc:`lamination_parameters` are built.


Material axes and engineering constants
---------------------------------------

The material axes of a ply are 1, the fibre direction, 2, the in-plane
direction transverse to the fibres, and 3, the normal to the ply. An
orthotropic material is described by the Young moduli `E_1`, `E_2`, `E_3`,
the shear moduli `G_{12}`, `G_{13}`, `G_{23}` and the Poisson's ratios
`\nu_{12}`, `\nu_{13}`, `\nu_{23}`, the attributes ``e1`` to ``nu23`` of
:class:`.MatLamina`. The other three Poisson's ratios follow from the
reciprocity relations `\nu_{ij}/E_i = \nu_{ji}/E_j`, which make the
compliance matrix symmetric [1]_, and are computed by
:func:`.read_laminaprop`:

.. math::
    :label: lam_reciprocity

    \nu_{21} = \nu_{12} \frac{E_2}{E_1}, \qquad
    \nu_{31} = \nu_{13} \frac{E_3}{E_1}, \qquad
    \nu_{32} = \nu_{23} \frac{E_3}{E_2}

:func:`.read_laminaprop` accepts 2, 6 or 9 constants, see
:doc:`conventions`: with 2, ``(E, nu)``, the material is isotropic with `G =
E/(2(1 + \nu))`; with 6 or 9 the shear moduli are taken as given, and
:class:`.MatLamina` never recomputes a modulus from the others. With 2 or 6
constants `E_3 = \nu_{13} = \nu_{23} = 0`.

The constants must be those of a positive-definite material (Jones [2]_,
Sec. 2.4.2): `E_i \ge 0` and `G_{ij} \ge 0`, Eq. (2.45), `|\nu_{ij}| <
(E_i/E_j)^{1/2}` for `E_j > 0`, Eq. (2.51), and

.. math::
    :label: lam_delta

    \Delta = 1 - \nu_{12}\nu_{21} - \nu_{23}\nu_{32} - \nu_{31}\nu_{13} -
    2\nu_{21}\nu_{32}\nu_{13} > 0

Eq. (2.48). :meth:`.MatLamina.rebuild` raises ``ValueError`` naming the
violated conditions. The strict inequalities of Jones are relaxed to admit
zero moduli as the limit of a positive-definite material, in which the
stiffnesses stay finite while the compliances diverge. Two limits are useful
in practice:

- the antiplane core of a sandwich, `E_1 = E_2 = G_{12} = 0` with finite
  `G_{13}`, `G_{23}`, "an idealized core in which the modulus of elasticity
  in planes parallel with the faces is zero but the shear modulus in planes
  perpendicular to the faces is finite" (Allen [3]_, p. 14), which adds only
  transverse shear stiffness to a laminate. With `E_i = 0 < E_j` the bound
  above requires `\nu_{ij} = 0`, and :func:`.read_laminaprop` stores
  `\nu_{ji} = 0`;
- a ply with a singular transverse shear stiffness, `G_{13} = 0` or `G_{23}
  = 0`, which is accepted by the material but rejected by every shear
  correction except ``'constant'`` and ``None``, see :doc:`shear_correction`.


Three-dimensional stiffness
---------------------------

The stiffness constants ``c11`` to ``c66`` of :class:`.MatLamina` are those
of the three-dimensional orthotropic material, the inverse of the compliance
matrix with `S_{11} = 1/E_1`, `S_{12} = -\nu_{12}/E_1`, `S_{13} =
-\nu_{13}/E_1`, `S_{22} = 1/E_2`, `S_{23} = -\nu_{23}/E_2`, `S_{33} = 1/E_3`,
`S_{44} = 1/G_{23}`, `S_{55} = 1/G_{13}`, `S_{66} = 1/G_{12}` [1]_. With
`\Delta` of Eq. :eq:`lam_delta`:

.. math::
    :label: lam_c3d

    C_{11} = \frac{E_1 (1 - \nu_{23}\nu_{32})}{\Delta}, \quad
    C_{12} = \frac{E_1 (\nu_{21} + \nu_{31}\nu_{23})}{\Delta}, \quad
    C_{13} = \frac{E_1 (\nu_{31} + \nu_{21}\nu_{32})}{\Delta} \\
    C_{22} = \frac{E_2 (1 - \nu_{13}\nu_{31})}{\Delta}, \quad
    C_{23} = \frac{E_2 (\nu_{32} + \nu_{12}\nu_{31})}{\Delta}, \quad
    C_{33} = \frac{E_3 (1 - \nu_{12}\nu_{21})}{\Delta} \\
    C_{44} = G_{23}, \quad C_{55} = G_{13}, \quad C_{66} = G_{12}

:meth:`.MatLamina.get_constitutive_matrix` returns this ``6 x 6`` matrix in
the order `\{11, 22, 33, 23, 13, 12\}`. With a zero `E_1` or `E_2` the same
constants are computed without dividing by the moduli, which keeps them
finite, e.g. an antiplane core given with 9 constants ``(0, 0, nu12, 0, G13,
G23, E3, 0, 0)`` has `C_{33} = E_3`, `C_{44} = G_{23}`, `C_{55} = G_{13}`
and all other constants zero.


Plane stress
------------

The plate and shell theories of :doc:`constitutive_matrices` neglect the
normal stress, `\sigma_{33} = 0`. Eliminating `\varepsilon_{33}` from the
three-dimensional relation gives the reduced stiffnesses `Q_{ij} = C_{ij} -
C_{i3} C_{j3}/C_{33}`, `i, j = 1, 2` (Jones [2]_, Eq. (2.64)), which depend on
six constants only [1]_ [2]_:

.. math::
    :label: lam_q

    Q_{11} = \frac{E_1}{1 - \nu_{12}\nu_{21}}, \quad
    Q_{12} = \frac{\nu_{12} E_2}{1 - \nu_{12}\nu_{21}}, \quad
    Q_{22} = \frac{E_2}{1 - \nu_{12}\nu_{21}} \\
    Q_{66} = G_{12}, \quad Q_{44} = G_{23}, \quad Q_{55} = G_{13}

The transverse shear stiffnesses `Q_{44}`, `Q_{55}` are kept, since the
FSDT and the TSDT use them; `E_3`, `\nu_{13}` and `\nu_{23}` do not enter.
These are the stiffnesses that :meth:`.Lamina.rebuild` rotates and that
:meth:`.MatLamina.rebuild` stores in ``q11`` to ``q66`` with its default
``plane_stress=True``, also when 9 constants are given. With
``rebuild(plane_stress=False)``, ``q11`` to ``q33`` and the invariants become
the three-dimensional `C_{ij}` of Eq. :eq:`lam_c3d`, which are stiffer, e.g.
`E(1 - \nu)/((1 + \nu)(1 - 2\nu))` instead of `E/(1 - \nu^2)` for an
isotropic material: they are meant for three-dimensional models, such as
continuum shell elements, and not for the laminate stiffnesses of
:class:`.Laminate`. The read-only attribute ``plane_stress`` records the
choice. Thermal expansion is not considered in the stiffnesses; the
strength allowables ``st1``, ``sc1``, ``ss12``, etc. are stored only.


Rotation to the laminate axes
-----------------------------

The ply angle `\theta` is measured from the laminate axis `x` to the
material direction 1, positive from `x` towards `y`, see :doc:`conventions`.
With `c = \cos\theta` and `s = \sin\theta`, the plane-stress stiffnesses of
the ply in the laminate axes, the attributes ``q11L`` to ``q55L`` computed
by :meth:`.Lamina.rebuild`, are [1]_:

.. math::
    :label: lam_qbar

    q_{11L} = Q_{11} c^4 + 2 (Q_{12} + 2 Q_{66}) s^2 c^2 + Q_{22} s^4 \\
    q_{12L} = (Q_{11} + Q_{22} - 4 Q_{66}) s^2 c^2 + Q_{12} (s^4 + c^4) \\
    q_{22L} = Q_{11} s^4 + 2 (Q_{12} + 2 Q_{66}) s^2 c^2 + Q_{22} c^4 \\
    q_{16L} = (Q_{11} - Q_{12} - 2 Q_{66}) s c^3 + (Q_{12} - Q_{22} +
    2 Q_{66}) s^3 c \\
    q_{26L} = (Q_{11} - Q_{12} - 2 Q_{66}) s^3 c + (Q_{12} - Q_{22} +
    2 Q_{66}) s c^3 \\
    q_{66L} = (Q_{11} + Q_{22} - 2 Q_{12} - 2 Q_{66}) s^2 c^2 + Q_{66}
    (s^4 + c^4) \\
    q_{44L} = Q_{44} c^2 + Q_{55} s^2, \quad
    q_{45L} = (Q_{55} - Q_{44}) s c, \quad
    q_{55L} = Q_{55} c^2 + Q_{44} s^2

In matrix form, `C^{(k)}` and `C_s^{(k)}` of ply `k` are

.. math::
    :label: lam_C_Cs

    C^{(k)} = \begin{bmatrix} q_{11L} & q_{12L} & q_{16L} \\
    q_{12L} & q_{22L} & q_{26L} \\ q_{16L} & q_{26L} & q_{66L}
    \end{bmatrix}^{(k)}, \qquad
    C_s^{(k)} = \begin{bmatrix} q_{44L} & q_{45L} \\ q_{45L} & q_{55L}
    \end{bmatrix}^{(k)}

relating `\{\sigma_{xx}, \sigma_{yy}, \tau_{xy}\}` to `\{\varepsilon_{xx},
\varepsilon_{yy}, \gamma_{xy}\}` and `\{\tau_{yz}, \tau_{xz}\}` to
`\{\gamma_{yz}, \gamma_{xz}\}`; :meth:`.Lamina.get_constitutive_matrix`
returns both blocks in one ``5 x 5`` matrix. The terms `q_{16L}`, `q_{26L}`
and `q_{45L}` are odd functions of `\theta` and vanish at 0 and 90 degrees,
the other terms are even; this is what makes some laminate stiffnesses
vanish in balanced, orthotropic or symmetric laminates, see
:doc:`constitutive_matrices`.

The transformation of the stresses from the laminate to the material axes,
e.g. to evaluate a failure criterion, is
:meth:`.Lamina.get_transf_matrix_stress_to_lamina`, with
`\sigma_{11} = c^2 \sigma_{xx} + s^2 \sigma_{yy} + 2 s c \, \tau_{xy}` as
its first row, and its inverse is
:meth:`.Lamina.get_transf_matrix_stress_to_laminate`, the same matrix for
`-\theta`; both are ``6 x 6`` in the order `\{11, 22, 33, 23, 13, 12\}`. The
displacements transform with the rotation matrix of
:meth:`.Lamina.get_transf_matrix_displ_to_lamina`.


Material invariants
-------------------

Written with the multiple angles `2\theta` and `4\theta`, the stiffnesses of
Eq. :eq:`lam_qbar` are linear combinations of five constants of the
material that do not depend on the ply angle, the invariants of Tsai and
Pagano (Jones [2]_, Sec. 2.7, Eq. (2.99)), completed by two transverse
shear invariants [4]_:

.. math::
    :label: lam_invariants

    u_1 = \frac{3 Q_{11} + 3 Q_{22} + 2 Q_{12} + 4 Q_{66}}{8}, \quad
    u_2 = \frac{Q_{11} - Q_{22}}{2}, \quad
    u_3 = \frac{Q_{11} + Q_{22} - 2 Q_{12} - 4 Q_{66}}{8} \\
    u_4 = \frac{Q_{11} + Q_{22} + 6 Q_{12} - 4 Q_{66}}{8}, \quad
    u_5 = \frac{u_1 - u_4}{2}, \quad
    u_6 = \frac{Q_{44} + Q_{55}}{2}, \quad
    u_7 = \frac{Q_{44} - Q_{55}}{2}

the attributes ``u1`` to ``u7`` of :class:`.MatLamina`, such that

.. math::
    :label: lam_qbar_invariants

    q_{11L} = u_1 + u_2 \cos 2\theta + u_3 \cos 4\theta, \quad
    q_{22L} = u_1 - u_2 \cos 2\theta + u_3 \cos 4\theta \\
    q_{12L} = u_4 - u_3 \cos 4\theta, \quad
    q_{66L} = u_5 - u_3 \cos 4\theta \\
    q_{16L} = \frac{u_2}{2} \sin 2\theta + u_3 \sin 4\theta, \quad
    q_{26L} = \frac{u_2}{2} \sin 2\theta - u_3 \sin 4\theta \\
    q_{44L} = u_6 + u_7 \cos 2\theta, \quad
    q_{45L} = -u_7 \sin 2\theta, \quad
    q_{55L} = u_6 - u_7 \cos 2\theta

The coefficients of this relation are the ``9 x 5`` matrix of
:meth:`.MatLamina.get_invariant_matrix`, whose columns multiply `\{1, \cos
2\theta, \sin 2\theta, \cos 4\theta, \sin 4\theta\}`. Integrated through the
thickness, the trigonometric terms become the lamination parameters, and the
laminate stiffnesses become linear in them, see :doc:`lamination_parameters`.
For an isotropic material `u_2 = u_3 = u_7 = 0`, and the ply angle has no
effect.


Trace and trace normalization
-----------------------------

The trace of the plane-stress stiffness matrix of a ply, `\mathrm{tr} = Q_{11}
+ Q_{22} + 2 Q_{66}` (Melo et al. [5]_, Eq. (1)), does not depend on the ply
angle: from Eq. :eq:`lam_qbar_invariants`, `q_{11L} + q_{22L} + 2 q_{66L} =
2 (u_1 + u_5)`. :meth:`.MatLamina.trace_normalize_plane_stress` divides
``q11``, ``q12``, ``q22``, ``q44``, ``q55``, ``q66`` and ``u1`` to ``u7`` by
the trace, in place, which makes them dimensionless. Since the laminate
stiffnesses of :doc:`lamination_parameters` are linear in the invariants,
those computed with the trace-normalized invariants are the dimensional ones
divided by `\mathrm{tr}`; this is the basis of the trace-based design of Tsai
and of the double-double laminates. The normalization requires the
plane-stress constants and a positive trace, i.e. it is refused after
``rebuild(plane_stress=False)`` and for an antiplane core.


References
----------

.. [1] Reddy, J. N., Mechanics of Laminated Composite Plates and Shells:
   Theory and Analysis, Second Edition, CRC Press, 2003.
   https://doi.org/10.1201/b12409
.. [2] Jones, R. M., Mechanics of Composite Materials, Second Edition,
   Taylor & Francis, Philadelphia, 1999.
.. [3] Allen, H. G., Analysis and Design of Structural Sandwich Panels,
   Pergamon Press, Oxford, 1969, Chap. 2.
   https://doi.org/10.1016/B978-0-08-012870-2.50006-7
.. [4] ABAQUS Example Problems Manual, Vol. 1, Example 1.2.2, Laminated
   composite shell: buckling of a cylindrical panel with a circular hole.
.. [5] Melo, J. D. D., Bi, J. and Tsai, S. W., A Novel Invariant-Based
   Design Approach to Carbon Fiber Reinforced Laminates, Composite
   Structures, 159, 44-52, 2017.
   https://doi.org/10.1016/j.compstruct.2016.09.055
