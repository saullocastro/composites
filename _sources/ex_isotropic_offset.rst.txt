Isotropic plate with an offset
==============================

For a homogeneous isotropic plate of thickness `h` whose mid-surface is at
`z = d` from the reference surface, i.e. ``offset = d``, the thickness
integrals have the closed form

.. math::

    A = h C, \qquad B = d h C, \qquad D = h \left( d^2 + \frac{h^2}{12}
    \right) C, \qquad A_{ts} = \frac{5}{6} G h I

where `C` is the plane-stress stiffness of the material and `I` the 2x2
identity matrix. The offset creates a membrane-bending coupling `B \neq 0`
and, in the same way, the odd moments of the density, see
:doc:`mass_integrals`. The test below checks these expressions, the
matrices `E`, `F`, `H` of the TSDT and the six mass integrals against the
integrals `\int_{d - h/2}^{d + h/2} z^n dz`:

.. literalinclude:: ../../tests/test_core.py
    :pyobject: test_example_isotropic_plate_offset
