Transverse shear stresses
=========================

:meth:`.Laminate.calc_transverse_shear_stress` gives the transverse shear
stresses at a height `z` for given shear forces `Q_y`, `Q_x`, with the
equilibrium distribution that defines ``Ats`` for
``shear_correction='rohwer'``, see :doc:`shear_correction`. The stresses
vanish at the faces, integrate to the shear forces, and their complementary
energy is `Q^T A_{ts}^{-1} Q / 2`:

.. literalinclude:: ../../tests/test_transverse_shear.py
    :pyobject: test_example_shear_stress_cross_ply

For a homogeneous plate the distribution is the parabola `\tau_{xz} = 3
Q_x/(2h) (1 - 4 \bar{z}^2/h^2)`, with `\bar{z}` measured from the
mid-surface, for any offset:

.. literalinclude:: ../../tests/test_transverse_shear.py
    :pyobject: test_isotropic_parabola

From the strain gradients of an actual FSDT solution,
:meth:`.Laminate.calc_equilibrium_transverse_shear` recovers the stresses by
integrating the equilibrium equations through the thickness, and
:meth:`.Laminate.calc_aposteriori_energy` gives their strain energies, here
for a homogeneous plate in cylindrical bending, recovering the factor 5/6:

.. literalinclude:: ../../tests/test_transverse_shear.py
    :pyobject: test_aposteriori_energy_homogeneous
