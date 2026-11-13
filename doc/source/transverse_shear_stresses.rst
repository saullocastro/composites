Transverse shear stresses and a posteriori energies
===================================================

The shear correction of :doc:`shear_correction` makes the FSDT store the
right transverse shear energy, but not the right stresses: `C_s^{(k)}
\gamma^{(0)}` is constant within each ply, discontinuous at the interfaces
and non-zero at the free faces, where the tractions are known to vanish. When
`\tau_{xz}` and `\tau_{yz}` are needed, e.g. for an interlaminar failure
criterion, they are recovered by integrating the three-dimensional
equilibrium equations through the thickness, from the in-plane stresses of
the plies. :class:`.Laminate` offers two recoveries: from the shear forces,
with the distribution of the equilibrium approach, and from the strain
gradients of an actual solution.


From the shear forces
---------------------

:meth:`.Laminate.calc_transverse_shear_stress` evaluates, at a height `z`
from the reference surface, the distribution of Rohwer [1]_, Eq.
:eq:`sc_rohwer_tau`:

.. math::
    :label: ts_f

    \begin{Bmatrix} \tau_{yz} \\ \tau_{xz} \end{Bmatrix} =
    f^{(k)}(z) \begin{Bmatrix} Q_y \\ Q_x \end{Bmatrix}

for the shear forces of the FSDT, e.g. `Q = A_{ts} \gamma^{(0)}`. The
distribution is quadratic within each ply, continuous at the interfaces,
zero at both faces, and its resultants are `Q`. It depends only on the
in-plane stiffnesses of the plies, through the inverse of the ``ABD``
matrix, and not on ``shear_correction``: it is stored by
:meth:`.Laminate.calc_transverse_shear_stiffness` with ``'rohwer'`` and
computed on the first call otherwise, and it does not depend on ``offset``,
being computed about the mid-surface. For a homogeneous plate it is the
parabola `\tau_{xz} = 3 Q_x/(2h) (1 - 4 \bar{z}^2/h^2)`, `\bar{z}` measured
from the mid-surface.

The distribution assumes the two cylindrical bending states of Rohwer, i.e.
constant membrane forces and the moment gradients `M_{xx,x} = Q_x`,
`M_{xy,x} = Q_y`, `M_{yy,y} = Q_y`, `M_{xy,y} = Q_x`. It is the natural
recovery when only the shear forces of a model are known.


From the strain gradients
-------------------------

When the derivatives of the generalized strains are available, e.g. from the
shape functions of a finite element or of a Ritz solution,
:meth:`.Laminate.calc_equilibrium_transverse_shear` integrates the
equilibrium of the actual in-plane stresses `\sigma^{(k)} = C^{(k)}
(\varepsilon^{(0)} + z \varepsilon^{(1)})` from the bottom face `z_1`:

.. math::
    :label: ts_recovery

    \tau_{xz}(z) = -\int_{z_1}^{z} \left( c_1 \left[
    \varepsilon^{(0)}_{,x} + \bar{z} \varepsilon^{(1)}_{,x} \right]
    + c_3 \left[ \varepsilon^{(0)}_{,y} + \bar{z}
    \varepsilon^{(1)}_{,y} \right] \right) d\bar{z} \\
    \tau_{yz}(z) = -\int_{z_1}^{z} \left( c_2 \left[
    \varepsilon^{(0)}_{,y} + \bar{z} \varepsilon^{(1)}_{,y} \right]
    + c_3 \left[ \varepsilon^{(0)}_{,x} + \bar{z}
    \varepsilon^{(1)}_{,x} \right] \right) d\bar{z}

i.e. `\tau_{xz,z} = -(\sigma_{xx,x} + \tau_{xy,y})` and `\tau_{yz,z} =
-(\tau_{xy,x} + \sigma_{yy,y})`, with `c_1`, `c_2`, `c_3` the rows of
`C^{(k)}` that give `\sigma_{xx}`, `\sigma_{yy}`, `\tau_{xy}`. The arguments
``eps_x`` and ``eps_y`` are the derivatives of `\{\varepsilon^{(0)},
\varepsilon^{(1)}\}` in the order of :attr:`.Laminate.ABD`, at any number of
points, and the stresses are returned at the heights ``z``, by default the
three Gauss-Legendre points of each ply. The integrals are exact, ply by
ply. This is the a posteriori approach of Noor and Peters [2]_: it uses the
gradients of the solution instead of assumed states. The stresses vanish at
the bottom face by construction, and at the top face only if the gradients
satisfy the in-plane equilibrium of the stress resultants, `N_{xx,x} +
N_{xy,y} = 0` and `N_{xy,x} + N_{yy,y} = 0`, which an approximate solution
does only approximately; the residual at the top face is a measure of that
error. With the gradients of the two cylindrical bending states, the
recovery gives back the distribution of Eq. :eq:`ts_f`.


A posteriori energies and factors
---------------------------------

:meth:`.Laminate.calc_aposteriori_energy` evaluates, with the recovered
stresses and the strains `\{\bar\gamma_{yz}, \bar\gamma_{xz}\} =
(C_s^{(k)})^{-1} \{\bar\tau_{yz}, \bar\tau_{xz}\}`, the transverse shear
strain energies per unit area of Noor and Peters [2]_, Eq. (7):

.. math::
    :label: ts_energy

    \bar{U}_{13} = \int \frac{1}{2} C_{55} \bar\gamma_{xz}^2 \, dz, \qquad
    \bar{U}_{23} = \int \frac{1}{2} C_{44} \bar\gamma_{yz}^2 \, dz

with `C_{44} = q_{44L}`, `C_{55} = q_{55L}` of each ply, together with the
resultants `Q_y = \int \bar\tau_{yz} dz`, `Q_x = \int \bar\tau_{xz} dz` of
the recovered stresses, integrated exactly with three Gauss-Legendre points
per ply. Equating, over the domain of the model, the FSDT energy with a
factor `k_{13}`, `Q_x^2/(2 k_{13} \bar{A}_{55})`, with the energy of the
recovered stresses gives the a posteriori factor

.. math::
    :label: ts_factor

    k_{13} = \frac{\int Q_x^2/(2 \bar{A}_{55}) \, dA}{\int \bar{U}_{13} \,
    dA}

and analogously `k_{23}` with `Q_y`, `\bar{A}_{44}` and `\bar{U}_{23}`
(Noor and Peters [2]_, Eqs. (8)-(9)), where `Q` are the shear forces of the
FSDT solution, obtained with an initial factor. For a homogeneous plate in
cylindrical bending, `\bar{U}_{13} = 3 Q_x^2/(5 G h)` and `k_{13} = 5/6`.
Unlike the factors of :doc:`shear_correction`, these depend on the problem:
on the loads, the boundary conditions and the position in the structure.
A ply with a singular `C_s^{(k)}` raises ``ValueError``.


Verification
------------

``tests/test_transverse_shear.py`` checks that the stresses vanish at the
faces, are continuous at the interfaces and integrate to the shear forces
for random laminates, that they agree with an independent implementation of
Rohwer's equations, that they do not depend on ``offset`` or on
``shear_correction``, that the recovery from the gradients of the two
cylindrical bending states gives Eq. :eq:`ts_f` with a traction-free top
face, and the energy and the factor 5/6 of a homogeneous plate.


References
----------

.. [1] Rohwer, K., Improved Transverse Shear Stiffnesses for Layered Finite
   Elements, DFVLR-FB 88-32, 1988.
.. [2] Noor, A. K. and Peters, J. M., A Posteriori Estimates for the Shear
   Correction Factors in Multi-layered Composite Cylinders, Journal of
   Engineering Mechanics, 115(6), 1225-1244, 1989.
   https://doi.org/10.1061/(asce)0733-9399(1989)115:6(1225)
