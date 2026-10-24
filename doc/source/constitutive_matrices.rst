Laminate constitutive matrices per plate theory
===============================================

A :class:`.Laminate` contains the thickness integrals of the ply stiffnesses
required by three equivalent single-layer theories: the classical laminated
plate theory (CLPT), the first-order shear deformation theory (FSDT) and the
third-order shear deformation theory (TSDT) of Reddy [1]_. All of them are
computed at once by :meth:`.Laminate.calc_constitutive_matrix`, which is
called by :func:`.laminated_plate` and :func:`.isotropic_plate`, such that
the same object serves a model of any of the three theories.


Displacement fields and generalized strains
-------------------------------------------

With `z` measured from the reference surface, the in-plane displacements of
the three theories are [1]_:

.. math::

    \text{CLPT:} \quad u = u_0 - z w_{,x}, \quad v = v_0 - z w_{,y} \\
    \text{FSDT:} \quad u = u_0 + z \phi_x, \quad v = v_0 + z \phi_y \\
    \text{TSDT:} \quad u = u_0 + z \phi_x - \frac{4}{3h^2} z^3 (\phi_x +
    w_{,x}), \quad v = v_0 + z \phi_y - \frac{4}{3h^2} z^3 (\phi_y + w_{,y})

and `w(x, y)` in all of them. The in-plane strains `\varepsilon =
\{\varepsilon_{xx}, \varepsilon_{yy}, \gamma_{xy}\}` and the transverse shear
strains `\gamma = \{\gamma_{yz}, \gamma_{xz}\}` are written with the
generalized strains of the reference surface as:

.. math::

    \varepsilon = \varepsilon^{(0)} + z \varepsilon^{(1)} + z^3
    \varepsilon^{(3)}, \qquad \gamma = \gamma^{(0)} + z^2 \gamma^{(2)}

where `\varepsilon^{(3)}` and `\gamma^{(2)}` exist only in the TSDT, and
`\gamma^{(0)}` vanishes in the CLPT. The FSDT gives a constant
`\gamma^{(0)}` through the thickness, which requires a shear correction, see
:doc:`shear_correction`, whereas the TSDT satisfies the traction-free
conditions at `z = \pm h/2` and needs no shear correction.


Thickness integrals
-------------------

With `C^{(k)}` the plane-stress stiffness of ply `k` in the laminate axes
and `C_s^{(k)}` its transverse shear stiffness, both given by
:meth:`.Lamina.rebuild`, and the ply occupying `z_k \le z \le z_{k+1}`:

.. math::

    (A, B, D, E, F, H) = \sum_k \int_{z_k}^{z_{k+1}} C^{(k)} (1, z, z^2, z^3,
    z^4, z^6) \, dz \\
    (\bar{A}_{ts}, D_{ts}, F_{ts}) = \sum_k \int_{z_k}^{z_{k+1}} C_s^{(k)}
    (1, z^2, z^4) \, dz

The FSDT does not use `\bar{A}_{ts}` directly but the corrected `A_{ts}`,
see :meth:`.Laminate.calc_transverse_shear_stiffness`.


Constitutive relations
----------------------

The stress resultants `N`, `M`, `P` are the integrals of `\{\sigma_{xx},
\sigma_{yy}, \tau_{xy}\}` times `1`, `z` and `z^3`, and `Q`, `R` those of
`\{\tau_{yz}, \tau_{xz}\}` times `1` and `z^2`, such that `Q = \{Q_y,
Q_x\}`. The strain energy per unit area of the reference surface is half the
product of the resultants by the generalized strains, and the constitutive
relations become:

.. math::

    \text{CLPT:} \quad
    \begin{Bmatrix} N \\ M \end{Bmatrix} =
    \begin{bmatrix} A & B \\ B & D \end{bmatrix}
    \begin{Bmatrix} \varepsilon^{(0)} \\ \varepsilon^{(1)} \end{Bmatrix}

.. math::

    \text{FSDT:} \quad
    \begin{Bmatrix} N \\ M \end{Bmatrix} =
    \begin{bmatrix} A & B \\ B & D \end{bmatrix}
    \begin{Bmatrix} \varepsilon^{(0)} \\ \varepsilon^{(1)} \end{Bmatrix},
    \qquad Q = A_{ts} \gamma^{(0)}

.. math::

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
`h`, offset of the mid-surface           ``h``, ``offset``                            all
=======================================  ===========================================  ================================

The subscripts 16, 26 and 66 refer to the `xy` components. In the
transverse shear blocks of the TSDT, some references denote `\bar{A}_{ts}`,
`D_{ts}` and `F_{ts}` with the same letters `A`, `D` and `F` used for the
in-plane stiffnesses; in ``composites`` they have their own names. The
moments of the density used in the mass matrices are described in
:doc:`mass_integrals`.


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

These are, for instance, the layouts of the Ritz models of `panels
<https://github.com/saullocastro/panels>`_. In an FSDT model, a shear
correction factor of the user's choice is applied to ``Abar_ts`` instead of
``Ats``, see :doc:`shear_correction`. The assembly of the TSDT matrix is
shown in :doc:`ex_tsdt`.

The TSDT displacement field above, with the traction-free conditions at `z =
\pm h/2`, takes the reference surface at the mid-surface, whereas all the
matrices of a :class:`.Laminate` are integrated with `z` measured from the
reference surface, i.e. including ``offset``, see :doc:`conventions`.


References
----------

.. [1] Reddy, J. N., Mechanics of Laminated Composite Plates and Shells, CRC
   Press, 2003. https://doi.org/10.1201/b12409
.. [2] Rohwer, K., Improved Transverse Shear Stiffnesses for Layered Finite
   Elements, DFVLR-FB 88-32, 1988.
