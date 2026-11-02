Defining a laminate: conventions
================================

A laminate is usually created with :func:`.laminated_plate`, from the ply
angles, thicknesses, material properties and densities, with
:func:`.isotropic_plate` for a single isotropic layer, from lamination
parameters with :func:`.laminate_from_LaminationParameters`, or as a
double-double or N-double laminate, see :doc:`lamination_parameters`. The
conventions below apply to every :class:`.Laminate`.


Stacking sequence and thickness coordinate
------------------------------------------

The angles of ``stack``, and the lists ``plyts``, ``laminaprops`` and
``rhos`` of :func:`.laminated_plate`, run from the bottom to the top ply.
The coordinate `z` is normal to the reference surface, pointing from the
bottom to the top face, and it is measured from the reference surface,
which is the surface of the generalized strains, see
:doc:`constitutive_matrices`.


Offset
------

``offset`` is the position of the mid-surface relative to the reference
surface, `d` in the theory pages. The plies run from

.. math::
    :label: conv_offset

    z_1 = -\frac{h}{2} + d \quad \text{to} \quad z_{N+1} = +\frac{h}{2} + d

such that ``offset = 0`` puts the reference surface at the mid-surface,
``offset = +h/2`` at the bottom face and ``offset = -h/2`` at the top face.
For example, a skin modelled on its mould line, with the laminate above it,
has ``offset = h/2``. With an offset, ``B`` is not zero even for an
isotropic plate, and the same holds for the odd moments of the density, see
:doc:`constitutive_matrices` and :doc:`mass_integrals`; ``A`` and the
transverse shear stiffnesses do not change.
:meth:`.Laminate.make_balanced`, :meth:`.Laminate.make_orthotropic`,
:meth:`.Laminate.make_symmetric` and :meth:`.Laminate.make_smeared` require
``offset = 0``. The lamination parameters of
:meth:`.Laminate.calc_lamination_parameters` are integrated about the
reference surface, i.e. they include the offset, and
:func:`.laminate_from_LaminationParameters` needs the same ``offset`` to
give back the stiffnesses, see :doc:`lamination_parameters`.


Material properties
-------------------

``laminaprop`` is read by :func:`.read_laminaprop`, with 2, 6 or 9 entries:

- ``(E, nu)``: isotropic, with ``G = E/(2*(1 + nu))`` for the three shear
  moduli;
- ``(E11, E22, nu12, G12, G13, G23)``: orthotropic ply in plane stress, the
  usual case;
- ``(E11, E22, nu12, G12, G13, G23, E33, nu13, nu23)``: orthotropic ply with
  the properties of the three-dimensional constitutive matrix of
  :meth:`.MatLamina.get_constitutive_matrix`.

The index 1 is the fibre direction, 2 the in-plane transverse direction and
3 the normal direction. The Young moduli may be zero, e.g. ``(0, 0, nu12, 0,
G13, G23)`` for an antiplane core, and a zero ``Ei`` with a non-zero ``Ej``
requires ``nuij = 0``. The stiffnesses used in the laminate integrals are
always those of plane stress, also with 9 entries; ``E33``, ``nu13`` and
``nu23`` enter only the three-dimensional constants. The constitutive
relations, the admissible constants and the invariants are described in
:doc:`lamina`.


Ply angle
---------

The ply angle `\theta`, in degrees, is measured from the laminate `x` axis to
the material direction 1, positive from `x` towards `y`, i.e. a positive
rotation about `z`: the direction 1 is `\{\cos\theta, \sin\theta\}` in the
laminate axes. Reversing the sign of all the angles changes the sign of the
16 and 26 terms and of ``A45``.


Units
-----

No units are assumed, and any consistent system can be used, e.g. N, m and
kg, with the densities in kg/m\ :sup:`3`, or N, mm and t, with the densities
in t/mm\ :sup:`3`. The dimension of each attribute is given in the
documentation of :class:`.Laminate`.


Ordering of stresses and strains
--------------------------------

- In-plane: `\{\sigma_{xx}, \sigma_{yy}, \tau_{xy}\}` and `\{\varepsilon_{xx},
  \varepsilon_{yy}, \gamma_{xy}\}`, with the engineering shear strain
  `\gamma_{xy} = 2 \varepsilon_{xy}`.
- Transverse shear: `\{\tau_{yz}, \tau_{xz}\}`, `\{\gamma_{yz},
  \gamma_{xz}\}`, with `\gamma = 2 \varepsilon`, and `\{Q_y, Q_x\}`; the
  index 4 refers to `yz` and the index 5 to `xz`.
- :attr:`.Laminate.ABD` relates `\{N_{xx}, N_{yy}, N_{xy}, M_{xx}, M_{yy},
  M_{xy}\}` to `\{\varepsilon^{(0)}_{xx}, \varepsilon^{(0)}_{yy},
  \gamma^{(0)}_{xy}, \varepsilon^{(1)}_{xx}, \varepsilon^{(1)}_{yy},
  \gamma^{(1)}_{xy}\}`, the same order used by the strain gradients of
  :meth:`.Laminate.calc_equilibrium_transverse_shear`.
- :meth:`.Lamina.get_constitutive_matrix` is ``5 x 5``, in the laminate axes,
  ordered as `\{xx, yy, xy, yz, xz\}`.
- :meth:`.MatLamina.get_constitutive_matrix` is ``6 x 6``, in the material
  axes, ordered as `\{11, 22, 33, 23, 13, 12\}`, as are the stress
  transformations of :meth:`.Lamina.get_transf_matrix_stress_to_lamina`.
- :meth:`.Laminate.calc_transverse_shear_stress` receives ``Qy, Qx`` and
  returns ``tau_yz, tau_xz``.
