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
:class:`.Laminate` object.

For isotropic plates, the :class:`.Laminate` object is also used for
convenience, and since offsetting the mid-surface is supported, there can be an
extension-bending coupling (B matrix) different than zero even for isotropic
plates.


The most convenient usage is probably with the
:func:`composites.utils.isotropic_plate` or the
:func:`composites.utils.laminated_plate` functions::

    from composites import laminated_plate

    laminaprop = (E11, E22, nu12, G12, G13, G23)
    plyt = ply_thickness
    stack = [0, 90, +45, -45]
    plate = laminated_plate(stack, plyt=plyt, laminaprop=laminaprop)


and with the :func:`composites.utils.isotropic_plate` function::

    from composites import isotropic_plate

    plate = isotropic_plate(thickness=5., E=E, nu=nu)

Where the laminate stiffness matrix, the often called ``ABD`` matrix, with
``shape=(6, 6)``, can be accessed using::

    >>> plate.ABD

and when transverse shear stiffnesses are required, the ``ABDE`` matrix, with
``shape=(8, 8)``::

    >>> plate.ABDE

.. automodule:: composites.core
    :members:

.. automodule:: composites.utils
    :members:

"""
import os

from .version import __version__
from .utils import isotropic_plate, laminated_plate

def get_include():
    return os.path.join(os.path.dirname(__file__))
