Thickness integrals for mass matrices
=====================================

The consistent mass matrix of an equivalent single-layer model follows from
the kinetic energy

.. math::

    T = \frac{1}{2} \int_V \rho \left( \dot{u}^2 + \dot{v}^2 + \dot{w}^2
    \right) dV

where the displacements `u`, `v`, `w` are those of the theory, see
:doc:`constitutive_matrices`. Integrating through the thickness separately
from the integration over the reference surface, the only laminate data left
are the moments of the density

.. math::

    I_n = \sum_k \int_{z_k}^{z_{k+1}} \rho^{(k)} z^n \, dz

computed by :meth:`.Laminate.calc_constitutive_matrix` with `z` measured from
the reference surface, as the stiffnesses:

======  ======================  ==============================
`n`     attribute               dimension
======  ======================  ==============================
0       ``intrho``              mass per unit area
1       ``intrhoz``             mass per unit area x length
2       ``intrhoz2``            mass per unit area x length^2
3       ``intrhoz3``            mass per unit area x length^3
4       ``intrhoz4``            mass per unit area x length^4
6       ``intrhoz6``            mass per unit area x length^6
======  ======================  ==============================


Moments required by each theory
-------------------------------

In the FSDT, `u = u_0 + z \phi_x`, and through the thickness:

.. math::

    \int \rho \dot{u}^2 dz = I_0 \dot{u}_0^2 + 2 I_1 \dot{u}_0
    \dot{\phi}_x + I_2 \dot{\phi}_x^2, \qquad
    \int \rho \dot{w}^2 dz = I_0 \dot{w}^2

and the same with `v_0` and `\phi_y`. The CLPT follows with `\phi_x =
-w_{,x}` and `\phi_y = -w_{,y}`. Both theories require `I_0`, `I_1` and
`I_2`, i.e. ``intrho``, ``intrhoz`` and ``intrhoz2``.

In the TSDT, `u = u_0 + z \phi_x + c z^3 \psi_x`, with `c = -4/(3h^2)` and
`\psi_x = \phi_x + w_{,x}`:

.. math::

    \int \rho \dot{u}^2 dz = I_0 \dot{u}_0^2 + 2 I_1 \dot{u}_0 \dot{\phi}_x
    + I_2 \dot{\phi}_x^2 + 2 c I_3 \dot{u}_0 \dot{\psi}_x + 2 c I_4
    \dot{\phi}_x \dot{\psi}_x + c^2 I_6 \dot{\psi}_x^2

and the same in `y`, such that `I_0`, `I_1`, `I_2`, `I_3`, `I_4` and `I_6`
are required, all the attributes of the table above.

In shells, the rotations that enter the in-plane displacements in `T` must
be the same as those used in the strains, e.g. including the corrections of
Sanders' kinematics in a cylinder, such that a rigid-body motion gives the
exact velocity of every point of the thickness.


Effect of the offset
--------------------

For a density symmetric about the reference surface the odd moments `I_1`
and `I_3` vanish. With ``offset`` different than zero, or an unsymmetric
distribution of the ply densities, `I_1 \neq 0` couples the in-plane and the
rotary inertias in the same way as ``B`` couples the membrane and bending
stiffnesses. For a homogeneous plate whose mid-surface is at `z = d`:

.. math::

    I_0 = \rho h, \qquad I_1 = \rho h d, \qquad I_2 = \rho h \left( d^2 +
    \frac{h^2}{12} \right)

see :doc:`ex_isotropic_offset`, where the higher moments are also checked.


Densities of the plies
----------------------

The density of each ply is the ``rho`` of its :class:`.MatLamina`. With
:func:`.laminated_plate`, ``rho`` gives the same density to all the plies and
``rhos`` one density per ply, e.g. for a sandwich with a light core. The
default density is zero, giving zero mass integrals, and the shear
correction ``'thickness_shear'`` requires positive densities, see
:doc:`shear_correction`.

The mass integrals are computed only from the plies: a :class:`.Laminate`
created from lamination parameters, see
:func:`.laminate_from_LaminationParameters`, or by
:func:`.double_double_laminate` or :func:`composites.utils.n_double_laminate`
has no plies,
and its mass integrals are zero.
