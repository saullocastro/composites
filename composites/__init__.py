r"""
===================================================================
Methods to calculate composite plate properties (:mod:`composites`)
===================================================================

.. currentmodule::composites

The ``composites`` module includes functions used to calculate plate properties
for laminated composites and isotropic plates.

Classical and first-order shear deformation theories are supported. For
classical plate theories or classical laminated plate theories (CLPT), the
relevant matrices are A, B, D, whereas for the first-order shear deformation
theories (FSDT) the matrices are A, B, D, E. All these matrices are part of the
`.Laminate` object. For isotropic plates, the `.Laminate` object is also used
for convenience.

Offset is supported, resulting in extension-bending coupling (B matrix)
different than zero even for isotropic plates.


The most convenient usage is probably with the
:func:`composites.laminate.read_stack()` function::

    from composites.laminate import read_stack

    laminaprop = (E11, E22, nu12, G12, G13, G23)
    plyt = ply_thickness
    stack = [0, 90, +45, -45]
    lam = read_stack(stack, plyt=plyt, laminaprop=laminaprop)


and with the :func:`composites.laminate.read_isotropic()` function::

    from composites.laminate import read_isotropic

    lam = read_isotropic(thickness=5., E=E, nu=nu)

Where the laminate stiffness matrix, the often called ``ABD`` matrix, with
``shape=(6, 6)``, can be accessed using::

    >>> lam.ABD

and when transverse shear stiffnesses are required, the ``ABDE`` matrix, with
``shape=(8, 8)``::

    >>> lam.ABDE

.. automodule:: composites.laminate
    :members:

.. automodule:: composites.lamina
    :members:

.. automodule:: composites.matlamina
    :members:

"""
from .version import __version__
