Defining a laminate: conventions
================================

A laminate is usually created with :func:`.laminated_plate`, from the ply
angles, thicknesses, material properties and densities, or with
:func:`.isotropic_plate` for a single isotropic layer::

    from composites import laminated_plate

    laminaprop = (142e9, 8.7e9, 0.28, 5.1e9, 5.1e9, 3.2e9)
    lam = laminated_plate([0, 45, -45, 90], plyt=0.125e-3,
                          laminaprop=laminaprop, rho=1600., offset=0.)

The conventions below apply to every :class:`.Laminate`.


Stacking sequence and thickness coordinate
------------------------------------------

The angles of ``stack``, and the lists ``plyts``, ``laminaprops`` and
``rhos``, run from the bottom to the top ply. The coordinate `z` is normal to
the reference surface, pointing from the bottom to the top face, and it is
measured from the reference surface, which is the surface of the
generalized strains, see :doc:`constitutive_matrices`.


Offset
------

``offset`` is the position of the mid-surface relative to the reference
surface, `d` below. The plies run from

.. math::

    z_1 = -\frac{h}{2} + d \quad \text{to} \quad z_{N+1} = +\frac{h}{2} + d

such that ``offset = 0`` puts the reference surface at the mid-surface,
``offset = +h/2`` at the bottom face and ``offset = -h/2`` at the top face.
For example, a skin modelled on its mould line, with the laminate above it,
has ``offset = h/2``. With an offset, ``B`` is not zero even for an
isotropic plate, `B = d h C`, see :doc:`ex_isotropic_offset`, and the same
holds for the odd moments of the density, see :doc:`mass_integrals`.
:meth:`.Laminate.make_balanced`, :meth:`.Laminate.make_orthotropic`,
:meth:`.Laminate.make_symmetric` and :meth:`.Laminate.make_smeared` require
``offset = 0``.

For any laminate, `B = B_{d=0} + d A` and `D = D_{d=0} + 2 d B_{d=0} + d^2
A`, the shift of the stiffness integrals to another reference plane
(Vasiliev and Morozov, Mechanics and Analysis of Composite Materials,
Elsevier, 2001, Eq. (5.28), p. 235), while `A` and the transverse shear
stiffnesses ``Abar44``, ``Abar45``, ``Abar55`` do not change. The lamination
parameters of :meth:`.Laminate.calc_lamination_parameters` are integrated
about the reference surface, i.e. they include the offset, and
:func:`.laminate_from_LaminationParameters` needs the same ``offset`` to
give back the stiffnesses:

.. literalinclude:: ../../tests/test_plane_stress_invariants.py
    :pyobject: test_lamination_parameters_offset


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

The Young moduli may be zero, e.g. ``(0, 0, nu12, 0, G13, G23)``, or ``(0,
0, nu12, 0, G13, G23, E33, 0, 0)``, for an antiplane core, "an idealized
core in which the modulus of elasticity in planes parallel with the faces is
zero but the shear modulus in planes perpendicular to the faces is finite"
(Allen, Analysis and Design of Structural Sandwich Panels, Pergamon Press,
1969, p. 14), which adds only transverse shear stiffness to a sandwich. A
zero ``Ei`` with a non-zero ``Ej`` requires ``nuij = 0``, see
:func:`.read_laminaprop`.

The index 1 is the fibre direction, 2 the in-plane transverse direction and 3
the normal direction. The ply stiffnesses used in the laminate integrals are
always the plane-stress ones, see :meth:`.Lamina.rebuild`; ``G13`` and
``G23`` give the transverse shear stiffnesses ``q44L``, ``q45L``, ``q55L``.
``E33``, ``nu13`` and ``nu23`` do not change them. The constitutive
constants and the invariants of :class:`.MatLamina`, used by
:func:`.laminate_from_LaminationParameters`, :class:`.GradABD` and
:meth:`.MatLamina.trace_normalize_plane_stress`, are by default those of the
same plane-stress stiffnesses, also with 9 entries, such that the
lamination parameters of :meth:`.Laminate.calc_lamination_parameters` give
back the stiffnesses of :func:`.laminated_plate`;
``MatLamina.rebuild(plane_stress=False)`` gives the three-dimensional ones,
e.g. for continuum shell elements, see :meth:`.MatLamina.rebuild`.


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
  axes, ordered as `\{11, 22, 33, 23, 13, 12\}`.
- :meth:`.Laminate.calc_transverse_shear_stress` receives ``Qy, Qx`` and
  returns ``tau_yz, tau_xz``.
