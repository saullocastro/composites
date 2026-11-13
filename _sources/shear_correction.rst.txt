Transverse shear stiffness and shear correction
===============================================

The FSDT gives transverse shear strains constant through the thickness,
whereas the three-dimensional equilibrium requires transverse shear stresses
that vary through the thickness and vanish at the traction-free faces. The
FSDT accepts the inconsistency and corrects it in an average sense: its
transverse shear stiffness `A_{ts}` is chosen smaller than the
constant-strain one, such that the constant distribution stores the right
energy, or reproduces another chosen response, for the same shear forces.
This chapter derives the corrections implemented by
:meth:`.Laminate.calc_transverse_shear_stiffness`, called at the end of
:meth:`.Laminate.calc_constitutive_matrix`, compares them and guides the
choice. The conventions are those of :doc:`constitutive_matrices`: `Q =
\{Q_y, Q_x\}`, `\gamma = \{\gamma_{yz}, \gamma_{xz}\}`, `A_{ts} = [[A_{44},
A_{45}], [A_{45}, A_{55}]]`, with `C_s^{(k)}` of Eq. :eq:`lam_C_Cs`, in
general a full matrix.


The correction factor
---------------------

In the simplest case the corrected stiffness is

.. math::
    :label: sc_kappa

    A_{ts} = \kappa \bar{A}_{ts}, \qquad \bar{A}_{ts} = \sum_k
    \int_{z_k}^{z_{k+1}} C_s^{(k)} dz

For a homogeneous plate of shear modulus `G`, the equilibrium of the
linear bending stresses gives the parabola `\tau_{xz} = 3 Q_x/(2h) (1 - 4
\bar{z}^2/h^2)`, `\bar{z}` measured from the mid-surface, whose
complementary energy per unit area is `\int \tau_{xz}^2/(2G) dz = 3
Q_x^2/(5 G h)`. Equating it with the energy `Q_x^2/(2 \kappa G h)` of the
constant distribution gives `\kappa = 5/6`. The factor is not a material
property but the answer to a chosen criterion: matching instead the
frequency of the first thickness-shear mode gives `\kappa = \pi^2/12`
(Mindlin [8]_). For a laminate the stress distribution is piecewise
parabolic, continuous at the interfaces with slope discontinuities, the
energy equivalence depends on the stacking sequence and on the direction,
and the correction becomes a pair of factors or a `2 \times 2` matrix.

Three transverse shear stiffnesses are distinguished, following Rohwer
[1]_, all computed with `z` measured from the mid-surface, such that none
depends on ``offset``:

.. math::
    :label: sc_three

    \bar{A}_{ts} = \sum_k C_s^{(k)} h_k, \qquad
    \bar{\bar{A}}_{ts} = h^2 \left[ \sum_k \left( C_s^{(k)} \right)^{-1}
    h_k \right]^{-1}, \qquad
    A_{ts} = \text{corrected}

- `\bar{A}_{ts}`, :attr:`.Laminate.Abar_ts`, constant strains through the
  thickness, the uncorrected stiffness used by the TSDT;
- `\bar{\bar{A}}_{ts}`, :attr:`.Laminate.Abarbar_ts`, constant stresses
  `Q/h` with the same complementary energy, for comparison only, ``nan``
  when a ply has a singular `C_s^{(k)}`;
- `A_{ts}`, :attr:`.Laminate.Ats`, ``A44``, ``A45``, ``A55``, with the
  correction selected by ``shear_correction`` already applied, the stiffness
  of an FSDT model as it is.

Neither `\bar{A}_{ts}` nor `\bar{\bar{A}}_{ts}` gives any reduction for a
homogeneous plate, and both overestimate the transverse shear stiffness of a
laminate.


Equilibrium approach (``'rohwer'``)
-----------------------------------

The default method, Rohwer [1]_, derives the stresses from the equilibrium
of two cylindrical bending states with constant membrane forces: one in
which the resultants depend on `x` only, with `M_{xx,x} = Q_x` and
`M_{xy,x} = Q_y`, and one in which they depend on `y` only, with `M_{yy,y} =
Q_y` and `M_{xy,y} = Q_x`. With `H^* = [[A, B], [B, D]]^{-1}` about the
mid-surface, the in-plane strains of a state are `(H^*_{1:3} + z H^*_{4:6})
\{N, M\}`, `H^*_{1:3}` and `H^*_{4:6}` being the first and the last three
rows of `H^*`. Integrating `\tau_{xz,z} = -\sigma_{xx,x}` and `\tau_{yz,z} =
-\sigma_{yy,y}` from the bottom face, with the rows `c_1^{(k)}`,
`c_2^{(k)}` of `C^{(k)}` that give `\sigma_{xx}`, `\sigma_{yy}` and `\Pi(z)
= z H^*_{1:3} + z^2/2 \, H^*_{4:6}`:

.. math::
    :label: sc_rohwer_tau

    \tau_{xz} = \left( a_1^{(k)} - c_1^{(k)} \Pi(z) \right) L_x Q, \qquad
    \tau_{yz} = \left( a_2^{(k)} - c_2^{(k)} \Pi(z) \right) L_y Q, \qquad
    a_\alpha^{(k)} = \sum_{i=1}^{k} \left( c_\alpha^{(i)} - c_\alpha^{(i-1)}
    \right) \Pi(z_i)

for `z_k \le z \le z_{k+1}`, with `c_\alpha^{(0)} = 0`, where `L_x` places
`Q_y`, `Q_x` at the components `M_{xy}`, `M_{xx}` of `\{N, M\}_{,x}` and
`L_y` places them at `M_{yy}`, `M_{xy}` of `\{N, M\}_{,y}`. The sums
`a_\alpha^{(k)}` make the stresses continuous at the interfaces, and the
stresses vanish at the top face because `\int c_1 (H^*_{1:3} + z H^*_{4:6})
dz` through the thickness is the row of `N_{xx}` of the identity, which
`L_x` does not select, and analogously for `c_2`. The result is the
distribution `\{\tau_{yz}, \tau_{xz}\} = f^{(k)}(z) Q`, quadratic within each
ply, and the stiffness follows from the complementary energy:

.. math::
    :label: sc_rohwer_Ats

    A_{ts} = \left[ \sum_k \int_{z_k}^{z_{k+1}} f^{(k)T}
    \left(C_s^{(k)}\right)^{-1} f^{(k)} dz \right]^{-1}

integrated exactly with three Gauss-Legendre points per ply, the integrand
being of degree 4. The method admits general anisotropic and unsymmetric
laminates, angle plies with full `C^{(k)}` and `C_s^{(k)}` included, and
gives the coupling `A_{45}` from the bending-twisting terms `D_{16}`,
`D_{26}` through `M_{xy,x} = Q_y`, also in balanced laminates, for which
`\bar{A}_{45} = 0`. It gives `\kappa = 5/6` for a homogeneous plate.

The two states are tied to the axes `x`, `y`, such that the result is not
invariant to a rotation of the axes, which is inherent to the method
(Rohwer [1]_, Fig. 7, p. 20): for his ``[0,90]_s`` CFRP laminate of Table 2,
computed in axes rotated by 45 degrees, ``A44 = A55`` are 7.0 % and ``A45``
is 13.5 % above the tensor rotation of the result in the axes of the plies.
A laminate whose material directions are rotated with respect to the axes
of a model is therefore created with the ply angles shifted by the rotation,
which gives the stiffness in the axes of the model, rather than rotating
``Ats``. The distribution `f^{(k)}(z)` is stored and gives the transverse
shear stresses of :doc:`transverse_shear_stresses`.


Scalar factors per direction
----------------------------

The other methods give one factor per direction, `k_{13}` for `xz` and
`k_{23}` for `yz`, applied to the constant-strain stiffness as

.. math::
    :label: sc_scalar

    A_{55} = k_{13} \bar{A}_{55}, \qquad A_{44} = k_{23} \bar{A}_{44},
    \qquad A_{45} = \frac{k_{13} + k_{23}}{2} \bar{A}_{45}

The averaged `A_{45}` is ad hoc, not the result of an energy equivalence,
and vanishes whenever `\bar{A}_{45} = 0`, e.g. in every balanced laminate.

- ``'vlachoutsis'``, Vlachoutsis [2]_: the laminate is treated in each
  direction as a beam in cylindrical bending, with `C_{11} = q_{11L}` and
  `C_{55} = q_{55L}` of each ply in the direction 1, and the
  one-dimensional equilibrium `\tau_{xz,z} = -\sigma_{xx,x}` about the
  direction-wise neutral surface `z_n`:

  .. math::
      :label: sc_vlachoutsis

      z_n = \frac{\int C_{11} z \, dz}{\int C_{11} dz}, \quad
      R = \int C_{11} (z - z_n)^2 dz, \quad
      g(z) = -\int_{z_1}^{z} C_{11} (\bar{z} - z_n) d\bar{z}, \quad
      k_{13} = \frac{R^2}{\bar{A}_{55} \int g^2/C_{55} \, dz}

  i.e. the energy of `\tau_{xz} = Q_x g(z)/R` equated with `Q_x^2/(2 k_{13}
  \bar{A}_{55})`; analogously `k_{23}` with `q_{22L}`, `q_{44L}` and
  `\bar{A}_{44}`. The scalar `C_{\alpha\alpha}` per direction ignores the
  Poisson coupling, the shear-extension coupling of angle plies and `C_{45}`,
  such that the factors are exact only for specially orthotropic plies.
  Vlachoutsis also gives closed forms for sandwiches.
- ``'whitney'``, Whitney [3]_, Eqs. (3)-(7), extending Chow to unsymmetric
  orthotropic laminates in cylindrical bending: the stress gradient
  `\sigma_{xx,x} = C_{11} (z - z_n) Q_x/R`, with `z_n = B_{11}/A_{11}` and `R
  = (D_{11} A_{11} - B_{11}^2)/A_{11}`, is that of ``'vlachoutsis'``, and the
  result is identical.
- ``'chow'``, Chow [4]_, Eqs. (11)-(12): the same with `z` measured from
  the mid-surface, `g(z) = -\int C_{11} \bar{z} d\bar{z}` and `k_{13} =
  D_{11}^2/(\bar{A}_{55} \int g^2/C_{55} dz)`, valid for laminates symmetric
  about the mid-surface only, for which it is identical to
  ``'vlachoutsis'``; an unsymmetric laminate raises ``ValueError``.
- ``'birman_bert'``, Birman and Bert [5]_, Eq. (13) and Appendix B: the
  thickness integral of the transverse shear strain of the equilibrium
  distribution, instead of its energy, is equated with that of the FSDT,
  `k_{13} = R h/(\bar{A}_{55} \int g/C_{55} dz)`, with `g` and `R` of
  ``'vlachoutsis'``. It gives `k = 1` for a homogeneous plate.
- ``'thickness_shear'``, Yang, Norris and Stavsky [6]_, as described by Sun
  and Whitney [7]_, extending the criterion of Mindlin [8]_: the frequency of
  the first thickness-shear mode of the FSDT, `\omega^2 = k_{13}
  \bar{A}_{55}/\bar{I}_2`, is matched with the exact frequency of the
  anti-plane motion of the laminate, `\rho^{(k)} \ddot{u} = (C_{55}^{(k)}
  u_{,z})_{,z}` with traction-free faces, computed with the exact transfer
  matrix of each ply. `\bar{I}_2 = I_2 - I_1^2/I_0` is the rotatory inertia
  about the centre of mass, see :doc:`mass_integrals`, and `C_{45}` is
  neglected. It gives `k = \pi^2/12` for a homogeneous plate and requires
  positive ply densities.
- ``'constant'``: `k = 5/6` for any laminate, ``A_ts = 5/6 Abar_ts``.
- ``None``: no correction, ``Ats = Abar_ts``.

The ratios ``scf_k13 = A55/Abar55`` and ``scf_k23 = A44/Abar44`` are stored
for every method. They are informative only: the correction is already
inside ``Ats`` and must not be applied again. They allow the methods to be
compared, or report the equivalent scalar factors of ``'rohwer'``.


Comparison and choice
---------------------

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

``'rohwer'`` and ``'vlachoutsis'`` integrate the same equilibrium equation in
cylindrical bending. For symmetric cross-ply laminates, where the couplings
that ``'vlachoutsis'`` ignores are small, they nearly coincide; for
unsymmetric and angle-ply laminates only ``'rohwer'`` accounts for the full
laminate compliance and gives a consistent `A_{45}`, which is why it is the
default. The scalar methods remain useful to reproduce published factors, to
compare with closed forms of sandwiches, and when one factor per direction
is required. ``'thickness_shear'`` is the dynamic criterion, appropriate when
the thickness-shear frequencies matter, and ``'birman_bert'`` the criterion
of the average strain discussed by its authors for sandwiches; the value
they recommend for sandwich structures, `K = 1` applied to `\bar{A}_{55}`,
is ``None``.

The choice matters in proportion to the transverse shear compliance: little
for thin plates, decisively for thick laminates and sandwiches with a soft
core, whose factors may be far below 5/6. For the ``[0,90]_s`` CFRP plate of
Rohwer [1]_, Sec. 4.2, simply supported, with `a/b = 1.5` and a distributed
load, the transverse shear influence on the centre deflection reaches 10 %
at `a/h = 19` and exceeds 30 % at `a/h = 11` (p. 21, Fig. 8). At `a/h = 10`,
"the Mindlin displacement exceeds the Kirchhoff one by 25% and 28% for Hbar
and Hbarbar, respectively, whereas with H it is 36% higher" (p. 21), `H`
being his equilibrium stiffness.

A model that requires its own factors, e.g. from a reference solution,
applies them to the uncorrected matrix, `A_{ts} = \kappa \bar{A}_{ts}` or
Eq. :eq:`sc_scalar` with ``lam.Abar44``, ``lam.Abar45``, ``lam.Abar55``.
Laminates created from lamination parameters have no plies, hence no
through-thickness distribution: their ``Ats`` is ``Abar_ts`` and
``shear_correction`` is ``None``, see :doc:`lamination_parameters`. The
method is chosen with the argument ``shear_correction`` of
:func:`.laminated_plate` and :func:`.isotropic_plate`, or by changing the
attribute ``shear_correction`` and calling
:meth:`.Laminate.calc_transverse_shear_stiffness` again.

All methods except ``'constant'`` and ``None`` raise ``ValueError`` for a ply
with a singular `C_s^{(k)}`, e.g. `G_{13} = 0`; ``'rohwer'`` also for a
singular ``ABD`` matrix, and the scalar methods for zero in-plane stiffness
in a direction.

The methods assume static states of the laminate alone. Factors of an actual
FSDT solution, in the sense of Noor and Peters [9]_, follow from the
recovered stresses of :doc:`transverse_shear_stresses`.


Verification
------------

``tests/test_transverse_shear.py`` compares the methods with the values
printed in the references, to half a unit of their last digit, e.g. the
``[0,90]_s`` laminate of Rohwer [1]_, Table 2, p. 19, the cross-ply laminates
of 120 plies and the isotropic sandwiches of Vlachoutsis [2]_, the laminates
of Whitney [3]_ and Chow [4]_, and the sandwich beams of Birman and Bert
[5]_; with closed forms, such as the factor of a symmetrical sandwich of
Vlachoutsis [2]_, Eqs. (44), (45) and (48), and, for an antiplane core, the
shear stiffness `G_c d^2/c` of Allen [10]_, Eq. (2.15b), p. 17; and with an
independent implementation of Rohwer's equations for random laminates, for
which the results are also checked to be independent of ``offset``.

Some printed values are not reproduced, and are checked against the closed
forms of their papers instead, as stated in the Notes of the tests. Chow
[4]_, Table 1, p. 316, gives `k_2^2 = 0.521` for the graphite-epoxy
``[0,90]_s`` plate, whereas his Eqs. (11)-(12), p. 311, give 0.52165 with his
data, Eq. (32), which reproduce the three other values of the table. Whitney
[3]_, p. 304, gives `k_1^2 = 0.4098` and `k_2^2 = 0.6915` for a sandwich,
whereas his Eqs. (3)-(7), p. 303, give 0.40936 and 0.67199 with his data.
The sandwich panel of Rohwer [1]_, Sec. 4.3, Table 4, p. 23, is not used:
its transverse shear stiffnesses are reproduced only with a foam shear
modulus of 0.04 kN/mm\ :sup:`2`, i.e. `\nu = 0.25` instead of the stated
0.3, and its bending stiffnesses contradict the stated lay-up.


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
.. [7] Sun, C. T. and Whitney, J. M., Theories for the dynamic response of
   laminated plates, AIAA Journal, 11, 178-183, 1973.
   https://doi.org/10.2514/3.50448
.. [8] Mindlin, R. D., Influence of rotatory inertia and shear on flexural
   motions of isotropic, elastic plates, Journal of Applied Mechanics, 18,
   31-38, 1951. https://doi.org/10.1115/1.4010217
.. [9] Noor, A. K. and Peters, J. M., A Posteriori Estimates for the Shear
   Correction Factors in Multi-layered Composite Cylinders, Journal of
   Engineering Mechanics, 115(6), 1225-1244, 1989.
   https://doi.org/10.1061/(asce)0733-9399(1989)115:6(1225)
.. [10] Allen, H. G., Analysis and Design of Structural Sandwich Panels,
   Pergamon Press, 1969, Chap. 2.
   https://doi.org/10.1016/B978-0-08-012870-2.50006-7
