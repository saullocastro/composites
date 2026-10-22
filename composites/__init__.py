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
A, B, D and Ats. For the third-order shear deformation theory (TSDT) the
relevant matrices are the A, B, D, E, F, H; and the Abar_ts, Dts and Fts.
The matrices Ats, Abar_ts, Dts and Fts are 2 by 2 matrices containing
transverse shear stiffnesses, ordered as ``[[44, 45], [45, 55]]`` with index 4
corresponding to ``yz`` and index 5 to ``xz``. All these matrices are part of
the :class:`.Laminate` object.

The FSDT transverse shear stiffness ``Ats`` already contains the shear
correction, computed by default with the method ``'rohwer'``, see
:meth:`.Laminate.calc_transverse_shear_stiffness`, such that no shear
correction factor should be applied to it. The uncorrected constant-strain
stiffness, used in the TSDT, is ``Abar_ts``. Other methods are selected with
the argument ``shear_correction`` of :func:`.laminated_plate` and
:func:`.isotropic_plate`:

- ``'rohwer'`` (default): equilibrium approach, full 2x2 matrix
- ``'vlachoutsis'``: scalar factors per direction
- ``'whitney'``: identical to ``'vlachoutsis'``
- ``'chow'``: symmetric laminates only
- ``'birman_bert'``: scalar factors per direction, average shear strain
- ``'thickness_shear'``: scalar factors per direction, dynamic, requires
  positive ply densities
- ``'constant'``: 5/6
- ``None``: no correction

The transverse shear stresses and their energies can also be recovered from
the strain gradients of an FSDT solution with
:meth:`.Laminate.calc_equilibrium_transverse_shear` and
:meth:`.Laminate.calc_aposteriori_energy`.

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
    lam = laminated_plate(stack, plyt=plyt, laminaprop=laminaprop)


and with the :func:`composites.utils.isotropic_plate` function::

    from composites import isotropic_plate

    lam = isotropic_plate(thickness=5., E=E, nu=nu)

Where the laminate stiffness matrix, the often called ``ABD`` matrix, with
``shape=(6, 6)``, can be accessed using::

    >>> lam.ABD

and when transverse shear stiffnesses are required, with ``shape=(2, 2)``::

    >>> lam.Ats

with the shear correction already applied. The transverse shear stresses
through the thickness, for given shear forces ``Qy`` and ``Qx``, are obtained
with::

    >>> tau_yz, tau_xz = lam.calc_transverse_shear_stress(z, Qy, Qx)

Objects are saved to and loaded from JSON files with::

    >>> from composites import save_json, load_json
    >>> save_json(lam, 'lam.json')
    >>> lam = load_json('lam.json')

see :mod:`composites.json_io`.

.. automodule:: composites.core
    :members:

.. automodule:: composites.utils
    :members:

.. automodule:: composites.kassapoglou
    :members:

.. automodule:: composites.json_io
    :members: to_dict, from_dict, to_json, from_json, save_json, load_json

Notes
-----
The shear correction methods are, see
:meth:`.Laminate.calc_transverse_shear_stiffness`:

- ``'rohwer'`` (default): Rohwer (1988) [1]_, equilibrium of two cylindrical
  bending states, full 2x2 matrix
- ``'vlachoutsis'``: Vlachoutsis (1992) [2]_, energy equivalence with
  direction-wise neutral surfaces
- ``'whitney'``: Whitney (1973) [3]_, unsymmetric orthotropic laminates,
  identical to ``'vlachoutsis'``
- ``'chow'``: Chow (1971) [4]_, symmetric laminates only
- ``'birman_bert'``: Birman and Bert (2002) [5]_, equivalence of the average
  transverse shear strain
- ``'thickness_shear'``: Yang, Norris and Stavsky (1966) [6]_, first
  thickness-shear frequency, requires positive ply densities

Problem-dependent, a posteriori factors in the sense of Noor and Peters (1989)
[7]_ are supported by :meth:`.Laminate.calc_equilibrium_transverse_shear` and
:meth:`.Laminate.calc_aposteriori_energy`, which recover the transverse shear
stresses and their energies from the strain gradients of an FSDT solution.

The implementation of the CLTP, FSDT and TSDT closely follows the notation
adopted by Reddy [8]_, see :class:`.Laminate`.

References
----------
.. [1] Rohwer, K., Improved transverse shear stiffness for layered finite
   elements, DFVLR-FB 88-32, 1988.
.. [2] Vlachoutsis, S., Shear correction factors for plates and shells,
   International Journal for Numerical Methods in Engineering, 33,
   1537-1552, 1992.
.. [3] Whitney, J. M., Shear correction factors for orthotropic laminates
   under static load, Journal of Applied Mechanics, 40, 302-304, 1973.
.. [4] Chow, T. S., On the propagation of flexural waves in an orthotropic
   laminated plate and its response to an impulsive load, Journal of
   Composite Materials, 5, 306-319, 1971.
.. [5] Birman, V. and Bert, C. W., On the choice of shear correction factor
   in sandwich structures, Journal of Sandwich Structures and Materials, 4,
   83-95, 2002.
.. [6] Yang, P. C., Norris, C. H. and Stavsky, Y., Elastic wave propagation
   in heterogeneous plates, International Journal of Solids and Structures,
   2, 665-684, 1966.
.. [7] Noor, A. K. and Peters, J. M., A posteriori estimates for the shear
   correction factors in multi-layered composite cylinders, Journal of
   Engineering Mechanics, 115, 1225-1244, 1989.
.. [8] Reddy, J. N., Mechanics of Laminated Composite Plates and Shells:
   Theory and Analysis, Second Edition, CRC Press, Boca Raton, 2004.

"""
import os

from .version import __version__
from .utils import isotropic_plate, laminated_plate
from .json_io import (to_dict, from_dict, to_json, from_json,
                      save_json, load_json)

def get_include():
    """Directory of the installed ``composites`` package

    Returns the directory that contains the ``__init__.py`` of this package.

    TODO(author): intended use, e.g. as an include directory for Cython code
    that cimports ``composites``.

    Returns
    -------
    path : str
        Absolute or relative path of the directory, as given by the
        ``__file__`` of this module.

    """
    return os.path.join(os.path.dirname(__file__))
