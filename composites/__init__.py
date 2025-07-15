r"""
==================================
composites API (:mod:`composites`)
==================================

.. currentmodule::composites

The ``composites`` module includes functions used to calculate properties and
perform analysis on laminated composites and isotropic plates.

Classical, first- and third-order shear deformation theories are supported. For
classical plate theories or classical laminated plate theories (CLPT), and for
the first-order shear deformation theory (FSDT) the relevant matrices are the
A, B, D and Atrans. For the third-order shear deformation theory (TSDT) the
relevant matrices are the A, B, D, E, F, H; and the Atrans, Dtrans and Ftrans.
The matrices indicated with "trans" are 2 by 2 matrices containing the
transverse shear stiffnesses.  All these matrices are part of the
:class:`.Laminate` object.

The implementation of the CLTP, FSDT and TSDT closely follows the notation
adopted by::

    Reddy J.N., Mechanics of laminated composite plates and shells, theory and
    analysis. Second Edition, Boca Raton: CRC Press, 2004.

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

and when transverse shear stiffnesses are required, with ``shape=(2, 2)``::

    >>> plate.Atrans

.. automodule:: composites.core
    :members:

.. automodule:: composites.utils
    :members:

.. automodule:: composites.kassapoglou
    :members:

"""
import os

from .version import __version__
from .utils import isotropic_plate, laminated_plate

def get_include():
    return os.path.join(os.path.dirname(__file__))
