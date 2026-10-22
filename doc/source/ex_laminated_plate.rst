Laminated plate and ABD matrix
==============================

:func:`.laminated_plate` creates a :class:`.Laminate` from the stacking
sequence, from the bottom to the top ply, see :doc:`conventions`. The
``ABD`` matrix of the CLPT and the FSDT, the corrected transverse shear
stiffness ``Ats``, the equivalent properties and the mass integrals are
computed at once. The code below is extracted from one of the ``composites``
unit tests:

.. literalinclude:: ../../tests/test_core.py
    :pyobject: test_example_laminated_plate

The matrices of all the theories of a laminate with three plies are compared
with reference values in:

.. literalinclude:: ../../tests/test_core.py
    :pyobject: test_laminated_plate_plane_stress
