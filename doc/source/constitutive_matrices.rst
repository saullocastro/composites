Laminate stiffness and plate theories
=====================================

A laminate is modelled as an equivalent single layer: the displacements are
assumed through the thickness in terms of functions of the reference surface
alone, and the three-dimensional strain energy, integrated through the
thickness, leaves the stiffness of the laminate as a set of thickness
integrals of the ply stiffnesses of :doc:`lamina`. A :class:`.Laminate`
contains the integrals required by three such theories: the classical
laminated plate theory (CLPT), the first-order shear deformation theory
(FSDT) and the third-order shear deformation theory (TSDT) of Reddy [1]_.
All of them are computed at once by :meth:`.Laminate.calc_constitutive_matrix`,
which is called by :func:`.laminated_plate` and :func:`.isotropic_plate`,
such that the same object serves a model of any of the three theories.


Displacement fields and generalized strains
-------------------------------------------

With `z` measured from the reference surface, the in-plane displacements of
the three theories are [1]_:

.. math::
    :label: cm_displacements

    \text{CLPT:} \quad u = u_0 - z w_{,x}, \quad v = v_0 - z w_{,y} \\
    \text{FSDT:} \quad u = u_0 + z \phi_x, \quad v = v_0 + z \phi_y \\
    \text{TSDT:} \quad u = u_0 + z \phi_x - \frac{4}{3h^2} z^3 (\phi_x +
    w_{,x}), \quad v = v_0 + z \phi_y - \frac{4}{3h^2} z^3 (\phi_y + w_{,y})

and `w(x, y)` in all of them, i.e. the transverse normal is inextensible,
`\varepsilon_{zz} = 0`. The CLPT keeps the normal straight and normal to the
deformed reference surface; the FSDT keeps it straight but lets it rotate
independently, with the rotations `\phi_x`, `\phi_y`; the TSDT lets it warp
with a cubic term chosen such that the transverse shear strains vanish at `z
= \pm h/2`. The in-plane strains `\varepsilon = \{\varepsilon_{xx},
\varepsilon_{yy}, \gamma_{xy}\}` and the transverse shear strains `\gamma =
\{\gamma_{yz}, \gamma_{xz}\}` are written with the generalized strains of the
reference surface as:

.. math::
    :label: cm_strains

    \varepsilon = \varepsilon^{(0)} + z \varepsilon^{(1)} + z^3
    \varepsilon^{(3)}, \qquad \gamma = \gamma^{(0)} + z^2 \gamma^{(2)}

where `\varepsilon^{(3)}` and `\gamma^{(2)}` exist only in the TSDT, and
`\gamma^{(0)}` vanishes in the CLPT. The FSDT gives a constant
`\gamma^{(0)}` through the thickness, which requires a shear correction, see
:doc:`shear_correction`, whereas the TSDT satisfies the traction-free
conditions at the faces and needs no shear correction.


Thickness integrals
-------------------

With `C^{(k)}` the plane-stress stiffness of ply `k` in the laminate axes
and `C_s^{(k)}` its transverse shear stiffness, Eq. :eq:`lam_C_Cs`, and the
ply occupying `z_k \le z \le z_{k+1}`, from `z_1 = -h/2 + d` to `z_{N+1} =
h/2 + d`, `d` being ``offset``, see :doc:`conventions`:

.. math::
    :label: cm_integrals

    (A, B, D, E, F, H) = \sum_k \int_{z_k}^{z_{k+1}} C^{(k)} (1, z, z^2, z^3,
    z^4, z^6) \, dz \\
    (\bar{A}_{ts}, D_{ts}, F_{ts}) = \sum_k \int_{z_k}^{z_{k+1}} C_s^{(k)}
    (1, z^2, z^4) \, dz

integrated exactly, ply by ply, e.g. `A_{ij} = \sum_k q_{ijL}^{(k)} (z_{k+1}
- z_k)`, `B_{ij} = \sum_k q_{ijL}^{(k)} (z_{k+1}^2 - z_k^2)/2` and `D_{ij} =
\sum_k q_{ijL}^{(k)} (z_{k+1}^3 - z_k^3)/3`. The FSDT does not use
`\bar{A}_{ts}` directly but the corrected `A_{ts}`, see
:meth:`.Laminate.calc_transverse_shear_stiffness` and
:doc:`shear_correction`. The moments of the density through the thickness,
computed in the same loop, are described in :doc:`mass_integrals`.


Constitutive relations
----------------------

The stress resultants `N`, `M`, `P` are the integrals of `\{\sigma_{xx},
\sigma_{yy}, \tau_{xy}\}` times `1`, `z` and `z^3`, and `Q`, `R` those of
`\{\tau_{yz}, \tau_{xz}\}` times `1` and `z^2`, such that `Q = \{Q_y,
Q_x\}`. The strain energy per unit area of the reference surface is half the
product of the resultants by the generalized strains, and the constitutive
relations become:

.. math::
    :label: cm_clpt

    \text{CLPT:} \quad
    \begin{Bmatrix} N \\ M \end{Bmatrix} =
    \begin{bmatrix} A & B \\ B & D \end{bmatrix}
    \begin{Bmatrix} \varepsilon^{(0)} \\ \varepsilon^{(1)} \end{Bmatrix}

.. math::
    :label: cm_fsdt

    \text{FSDT:} \quad
    \begin{Bmatrix} N \\ M \end{Bmatrix} =
    \begin{bmatrix} A & B \\ B & D \end{bmatrix}
    \begin{Bmatrix} \varepsilon^{(0)} \\ \varepsilon^{(1)} \end{Bmatrix},
    \qquad Q = A_{ts} \gamma^{(0)}

.. math::
    :label: cm_tsdt

    \text{TSDT:} \quad
    \begin{Bmatrix} N \\ M \\ P \end{Bmatrix} =
    \begin{bmatrix} A & B & E \\ B & D & F \\ E & F & H \end{bmatrix}
    \begin{Bmatrix} \varepsilon^{(0)} \\ \varepsilon^{(1)} \\
    \varepsilon^{(3)} \end{Bmatrix}, \qquad
    \begin{Bmatrix} Q \\ R \end{Bmatrix} =
    \begin{bmatrix} \bar{A}_{ts} & D_{ts} \\ D_{ts} & F_{ts} \end{bmatrix}
    \begin{Bmatrix} \gamma^{(0)} \\ \gamma^{(2)} \end{Bmatrix}

The transverse shear matrices are ordered as ``[[44, 45], [45, 55]]``, with
the index 4 corresponding to `yz` and the index 5 to `xz`, in the same order
as `\gamma = \{\gamma_{yz}, \gamma_{xz}\}` and `Q = \{Q_y, Q_x\}`. Some
references, e.g. Rohwer [2]_, use the opposite order `\{xz, yz\}`.


Symbols and attributes
----------------------

=======================================  ===========================================  ================================
Symbol                                   :class:`.Laminate` attribute                 Theories
=======================================  ===========================================  ================================
`A`, `B`, `D`                            :attr:`.Laminate.A`, :attr:`.Laminate.B`,    CLPT, FSDT, TSDT
                                         :attr:`.Laminate.D`, ``A11`` to ``D66``
`[[A, B], [B, D]]`                       :attr:`.Laminate.ABD`                        CLPT, FSDT
`E`, `F`, `H`                            :attr:`.Laminate.E`, :attr:`.Laminate.F`,    TSDT
                                         :attr:`.Laminate.H`, ``E11`` to ``H66``
`A_{ts}`, corrected                      :attr:`.Laminate.Ats`, ``A44``, ``A45``,     FSDT
                                         ``A55``
`\bar{A}_{ts}`, uncorrected              :attr:`.Laminate.Abar_ts`, ``Abar44``,       TSDT
                                         ``Abar45``, ``Abar55``
`D_{ts}`                                 :attr:`.Laminate.Dts`, ``D44``, ``D45``,     TSDT
                                         ``D55``
`F_{ts}`                                 :attr:`.Laminate.Fts`, ``F44``, ``F45``,     TSDT
                                         ``F55``
`\bar{\bar{A}}_{ts}`, constant stress    :attr:`.Laminate.Abarbar_ts`                 for comparison
`C^{(k)}`                                ``q11L``, ``q12L``, ``q16L``, ``q22L``,      all
                                         ``q26L``, ``q66L`` of :class:`.Lamina`
`C_s^{(k)}`                              ``q44L``, ``q45L``, ``q55L`` of              all
                                         :class:`.Lamina`
`h`, offset of the mid-surface `d`       ``h``, ``offset``                            all
=======================================  ===========================================  ================================

The subscripts 16, 26 and 66 refer to the `xy` components. In the
transverse shear blocks of the TSDT, some references denote `\bar{A}_{ts}`,
`D_{ts}` and `F_{ts}` with the same letters `A`, `D` and `F` used for the
in-plane stiffnesses; in ``composites`` they have their own names.


Reference surface
-----------------

All the integrals are taken with `z` measured from the reference surface,
the surface of the generalized strains. Shifting the reference surface by
`d` from the mid-surface changes the odd and the higher integrals, but not
`A` nor the transverse shear stiffness `\bar{A}_{ts}`: with the superscript
0 for the integrals about the mid-surface (Vasiliev and Morozov [3]_, Eq.
(5.28), p. 235),

.. math::
    :label: cm_shift

    B = B^0 + d A, \qquad D = D^0 + 2 d B^0 + d^2 A

such that `B \ne 0` for any laminate with an offset, `B = d h C` for an
isotropic plate of plane-stress stiffness `C`. The equivalent properties,
their conditioning check and the shear corrections are evaluated about the
mid-surface, which makes them independent of ``offset`` and avoids the
round-off of large offsets.


Special laminates
-----------------

Since `q_{16L}`, `q_{26L}` and `q_{45L}` are odd functions of the ply angle
and the other terms are even, see :doc:`lamina`, the stacking sequence makes
some integrals vanish. Four methods impose these simplifications on the
stored stiffnesses, all of them requiring ``offset = 0``:

- symmetric laminate, :meth:`.Laminate.make_symmetric`: the integrands of
  the odd powers of `z` are odd about the mid-surface, `B = E = 0`;
- balanced laminate, with a ply at `-\theta` of the same material and
  thickness for every ply at `+\theta`, anywhere in the stack,
  :meth:`.Laminate.make_balanced`: only the integrals of `z^0` of the odd
  terms vanish, `A_{16} = A_{26} = \bar{A}_{45} = \bar{\bar{A}}_{45} = 0`,
  whereas `B_{16}`, `D_{16}`, ..., `D_{45}`, `F_{45}` depend on the positions
  of the plies and are kept. ``A45`` is also set to zero, except with the
  shear correction ``'rohwer'``, whose `A_{45}` follows from the
  bending-twisting coupling and does not vanish in general, see
  :doc:`shear_correction`;
- orthotropic laminate, all plies at 0 or 90 degrees,
  :meth:`.Laminate.make_orthotropic`: every integral of `q_{16L}`,
  `q_{26L}` and `q_{45L}` vanishes, i.e. the terms 16 and 26 of `A` to `H`
  and the terms 45 of the transverse shear matrices;
- smeared laminate, the limit of many repetitions of the stacking sequence,
  :meth:`.Laminate.make_smeared`: `B = 0` and `D = h^2 A/12`, the
  stiffnesses of a homogeneous plate with the in-plane stiffness `A/h`. The
  double-double and N-double laminates of :doc:`lamination_parameters` are
  built with this assumption.

The same simplifications of the lamination parameters are made by
:func:`.make_symmetric_LP`, :func:`.make_balanced_LP` and
:func:`.make_orthotropic_LP`.


Equivalent engineering properties
---------------------------------

:meth:`.Laminate.calc_equivalent_properties`, called by
:func:`.laminated_plate`, gives the membrane moduli of the laminate with
zero moments and free curvatures, from `a = [[A, B], [B, D]]^{-1}`:

.. math::
    :label: cm_equivalent

    E_1 = \frac{1}{h a_{11}}, \quad E_2 = \frac{1}{h a_{22}}, \quad
    G_{12} = \frac{1}{h a_{66}}, \quad
    \nu_{12} = -\frac{a_{12}}{a_{11}}, \quad
    \nu_{21} = -\frac{a_{12}}{a_{22}}

where the index 6 is the third row and column of `a`, stored in ``e1``,
``e2``, ``g12``, ``nu12`` and ``nu21`` of the :class:`.Laminate`. For an
unsymmetric laminate these are not the moduli of `A^{-1}`, since the
curvatures are free. They are undefined when the ``ABD`` matrix is
singular, e.g. for a laminate of antiplane cores only, or for a single ply
at an angle with `G_{12} = 0`, whose ``ABD`` has a positive diagonal and is
singular nevertheless. The method raises ``numpy.linalg.LinAlgError`` when
the diagonal is not positive, when the inverse is not finite, or when the
1-norm condition number of the matrix about the mid-surface, scaled by its
diagonal such that it does not depend on the units, exceeds `10^{12}`.


Choosing the theory
-------------------

The three theories use the same in-plane integrals `A`, `B`, `D` and differ
in the treatment of the transverse shear:

- the CLPT neglects the transverse shear deformation and needs only `A`,
  `B`, `D`. It is the theory of thin laminates, and the one of the closed
  forms of :doc:`buckling_kassapoglou`;
- the FSDT adds `A_{ts}`, with the accuracy of the shear correction. The
  transverse shear compliance relative to the bending compliance grows with
  `(h/L)^2`, `L` being an in-plane length, and with the ratio of the
  in-plane to the transverse shear moduli, such that "in CFRP laminates the
  transverse shear effect must be accounted for up to higher slenderness
  rates than in homogeneous structures" (Rohwer [2]_, p. 24), and in
  sandwiches always;
- the TSDT represents the parabolic variation of the transverse shear
  strains and needs no correction, at the cost of the additional generalized
  strains `\varepsilon^{(3)}`, `\gamma^{(2)}` and of the matrices `E`, `F`,
  `H`, `D_{ts}`, `F_{ts}`. Its displacement field, Eq.
  :eq:`cm_displacements`, assumes the reference surface at the mid-surface
  for the traction-free conditions at `\pm h/2`, whereas all the matrices of
  a :class:`.Laminate` are integrated about the reference surface, including
  ``offset``.

Laminates created from lamination parameters, or as double-double or
N-double laminates, have only part of these stiffnesses, see
:doc:`lamination_parameters`.


Assembling the matrices in a model
----------------------------------

A finite element or Ritz model gathers the generalized strains in one vector
`e` and the blocks above in one symmetric matrix `K_c`, such that the strain
energy per unit area is `e^T K_c e / 2` and the stiffness matrix is
`\int B^T K_c B \, dA`, where `B` relates the degrees of freedom to `e`. The
common layouts are:

- CLPT, ``6 x 6``: `e = \{\varepsilon^{(0)}, \varepsilon^{(1)}\}`, `K_c =`
  ``lam.ABD``.
- FSDT, ``8 x 8``: `e = \{\varepsilon^{(0)}, \varepsilon^{(1)},
  \gamma^{(0)}\}`, ``K_c = [[ABD, 0], [0, Ats]]``.
- TSDT, ``13 x 13``: `e = \{\varepsilon^{(0)}, \varepsilon^{(1)},
  \varepsilon^{(3)}, \gamma^{(0)}, \gamma^{(2)}\}`, ``K_c = [[A, B, E, 0,
  0], [B, D, F, 0, 0], [E, F, H, 0, 0], [0, 0, 0, Abar_ts, Dts], [0, 0, 0,
  Dts, Fts]]``.

In an FSDT model with a shear correction factor of the user's choice, the
factor is applied to ``Abar_ts`` instead of using ``Ats``, see
:doc:`shear_correction`.


References
----------

.. [1] Reddy, J. N., Mechanics of Laminated Composite Plates and Shells, CRC
   Press, 2003. https://doi.org/10.1201/b12409
.. [2] Rohwer, K., Improved Transverse Shear Stiffnesses for Layered Finite
   Elements, DFVLR-FB 88-32, 1988.
.. [3] Vasiliev, V. V. and Morozov, E. V., Mechanics and Analysis of
   Composite Materials, Elsevier, Oxford, 2001.
