Choosing the transverse shear correction
========================================

The FSDT gives transverse shear strains constant through the thickness,
whereas equilibrium requires transverse shear stresses that vary through the
thickness and vanish at the faces. The strain energy of the constant
distribution is corrected by a transverse shear stiffness `A_{ts}` smaller
than the constant-strain one, in the simplest case

.. math::

    A_{ts} = \kappa \bar{A}_{ts}, \qquad \bar{A}_{ts} = \sum_k
    \int_{z_k}^{z_{k+1}} C_s^{(k)} dz

with `\kappa = 5/6` for a homogeneous plate. For a laminate the correction
depends on the stacking sequence and on the direction, and it is computed
from the plies by :meth:`.Laminate.calc_transverse_shear_stiffness`, called
by :meth:`.Laminate.calc_constitutive_matrix`.


The corrected stiffness ``Ats``
-------------------------------

``A44``, ``A45``, ``A55``, i.e. :attr:`.Laminate.Ats`, always contain the
shear correction, and they are the transverse shear stiffnesses to be used in
an FSDT model as they are. The uncorrected ``Abar_ts`` is the one used in the
TSDT, which needs no correction, see :doc:`constitutive_matrices`.

The method is chosen with the argument ``shear_correction`` of
:func:`.laminated_plate` and :func:`.isotropic_plate`, or by changing the
attribute ``shear_correction`` and calling
:meth:`.Laminate.calc_transverse_shear_stiffness` again::

    lam = laminated_plate(stack, plyt=plyt, laminaprop=laminaprop,
                          shear_correction='vlachoutsis')
    lam.shear_correction = 'rohwer'
    lam.calc_transverse_shear_stiffness()


Methods and their scope
-----------------------

==========================  ==================================================  ===================
``shear_correction``        scope                                               homogeneous plate
==========================  ==================================================  ===================
``'rohwer'`` (default)      any laminate, also unsymmetric, full 2x2 matrix     5/6
                            including ``A45``, from equilibrium [1]_
``'vlachoutsis'``           scalar factor per direction [2]_, exact for         5/6
                            specially orthotropic plies
``'whitney'``               unsymmetric orthotropic laminates [3]_, same        5/6
                            result as ``'vlachoutsis'``
``'chow'``                  symmetric laminates only [4]_                       5/6
``'birman_bert'``           equivalence of the average shear strain [5]_        1
``'thickness_shear'``       first thickness-shear frequency [6]_, positive      `\pi^2/12`
                            ply densities required
``'constant'``              `\kappa = 5/6` for any laminate                     5/6
``None``                    no correction, ``Ats = Abar_ts``                    1
==========================  ==================================================  ===================

The scalar methods apply their factors as ``A55 = k13*Abar55``, ``A44 =
k23*Abar44`` and ``A45 = (k13 + k23)/2*Abar45``, which gives ``A45 = 0``
whenever ``Abar45 = 0``, e.g. for laminates of ``0``, ``90`` and balanced
``+-45`` plies, whereas ``'rohwer'`` captures the coupling of the
bending-twisting terms ``D16``, ``D26``. The details of each method are in
:meth:`.Laminate.calc_transverse_shear_stiffness`.

The result of ``'rohwer'`` is tied to the laminate axes `x`, `y`, i.e. it is
not invariant to a rotation of the axes. When a laminate is used with its
material direction rotated with respect to the axes of a model, a new
:class:`.Laminate` with the ply angles shifted by the rotation gives the
transverse shear stiffness in the axes of the model, instead of a tensor
rotation of ``Ats``.


``scf_k13`` and ``scf_k23``
---------------------------

The attributes ``scf_k13 = A55/Abar55`` and ``scf_k23 = A44/Abar44`` are
informative only: the correction is already inside ``Ats``, and they must not
be applied to it again. They are useful to compare the methods, or to report
the equivalent scalar factors of ``'rohwer'``.


Using another correction factor
-------------------------------

A model that requires its own shear correction factors, e.g. from a
reference solution, applies them to the uncorrected matrix::

    k13, k23 = 0.8, 0.75
    Ats_user = np.array([[k23*lam.Abar44, (k13 + k23)/2*lam.Abar45],
                         [(k13 + k23)/2*lam.Abar45, k13*lam.Abar55]])

or ``k*lam.Abar_ts`` with a single factor. With ``shear_correction=None``,
``Ats`` is equal to ``Abar_ts``.


Problem-dependent factors
-------------------------

The methods above assume static states of the laminate alone, e.g. the
cylindrical bending states of ``'rohwer'``. From the strain gradients of an
actual FSDT solution, :meth:`.Laminate.calc_equilibrium_transverse_shear`
recovers the transverse shear stresses by integrating the equilibrium
equations through the thickness, and :meth:`.Laminate.calc_aposteriori_energy`
gives their strain energies, which lead to a posteriori factors in the sense
of Noor and Peters [7]_. See :doc:`ex_transverse_shear`.


References
----------

.. [1] Rohwer, K., Improved Transverse Shear Stiffnesses for Layered Finite
   Elements, DFVLR-FB 88-32, 1988.
.. [2] Vlachoutsis, S., Shear correction factors for plates and shells,
   International Journal for Numerical Methods in Engineering, 33(7),
   1537-1552, 1992. https://doi.org/10.1002/nme.1620330712
.. [3] Whitney, J. M., Shear Correction Factors for Orthotropic Laminates
   Under Static Load, Journal of Applied Mechanics, 40(1), 302-304, 1973.
   https://doi.org/10.1115/1.3422950
.. [4] Chow, T. S., On the Propagation of Flexural Waves in an Orthotropic
   Laminated Plate and Its Response to an Impulsive Load, Journal of
   Composite Materials, 5(3), 306-319, 1971.
   https://doi.org/10.1177/002199837100500302
.. [5] Birman, V. and Bert, C. W., On the Choice of Shear Correction Factor
   in Sandwich Structures, Journal of Sandwich Structures and Materials,
   4(1), 83-95, 2002. https://doi.org/10.1177/1099636202004001180
.. [6] Yang, P. C., Norris, C. H. and Stavsky, Y., Elastic wave propagation
   in heterogeneous plates, International Journal of Solids and Structures,
   2(4), 665-684, 1966. https://doi.org/10.1016/0020-7683(66)90045-x
.. [7] Noor, A. K. and Peters, J. M., A Posteriori Estimates for the Shear
   Correction Factors in Multi-layered Composite Cylinders, Journal of
   Engineering Mechanics, 115(6), 1225-1244, 1989.
   https://doi.org/10.1061/(asce)0733-9399(1989)115:6(1225)
