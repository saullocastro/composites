Buckling of rectangular plates (Kassapoglou)
============================================

The functions of :mod:`composites.kassapoglou` give closed-form and two-term
Rayleigh-Ritz estimates of the buckling loads of flat rectangular laminated
plates and of the effective width of a post-buckled plate, following
Kassapoglou [1]_, chapters 6 and 7. They are preliminary-design tools: they
take the bending stiffnesses `D_{ij}` of the CLPT, e.g. ``lam.D11`` of a
:class:`.Laminate`, assume simply supported edges and neglect the transverse
shear deformation, see :doc:`constitutive_matrices`. The loads are
resultants per unit length, positive in compression for the returned
critical loads.


Energy and assumptions
----------------------

For a deflection `w(x, y)` of a plate of length `a` along `x` and width `b`
along `y`, under the in-plane resultants `N_x`, `N_{xy}`, positive in
tension, the total potential energy is (Kassapoglou [1]_, Eq. 6.32, here
with `D_{16}` and `D_{26}`):

.. math::
    :label: kp_energy

    \Pi = \int_0^a \int_0^b \left[ \frac{1}{2} \left( D_{11} w_{,xx}^2
    + 2 D_{12} w_{,xx} w_{,yy} + D_{22} w_{,yy}^2 + 4 D_{66} w_{,xy}^2 + 4
    D_{16} w_{,xx} w_{,xy} + 4 D_{26} w_{,yy} w_{,xy} \right) +
    \frac{1}{2} N_x w_{,x}^2 + N_{xy} w_{,x} w_{,y} \right] dx \, dy

With a symmetric laminate, `B = 0`, the bending is decoupled from the
membrane state, which is uniform; the methods below assume this. The
bending-twisting terms `D_{16}`, `D_{26}` are neglected by all functions but
:func:`.calc_Nxx_crit_combined_shear_full`.


Uniaxial compression
--------------------

With the mode `w = \sin(m \pi x/a) \sin(n \pi y/b)` of a specially
orthotropic plate, `D_{16} = D_{26} = 0`, simply supported on four edges,
the buckling load is exact (Kassapoglou [1]_, Eq. 6.6, p. 129), with `AR =
a/b`:

.. math::
    :label: kp_nxx

    N_{xx,crit} = \frac{\pi^2}{a^2 m^2} \left( D_{11} m^4 + 2 (D_{12} +
    2 D_{66}) m^2 n^2 AR^2 + D_{22} n^4 AR^4 \right)

:func:`.calc_Nxx_crit` evaluates it for the given numbers of half-waves, and
searches the critical one, from 1 to 10, in the direction given as
``None``; usually ``n = 1`` and ``m = None``.


Shear
-----

For a long plate of width `a`, simply supported on its long edges, under
the shear load `N_{xy}`, the mode of Kassapoglou [1]_, Eq. 6.20,

.. math::
    :label: kp_shear_mode

    w = w_o \sin\frac{\pi x}{a} \sin\frac{\pi (y - x \tan\alpha)}{L}

has nodal lines inclined by `\alpha` and separated by `L`, such that the
length of the plate does not enter. Minimizing the energy with respect to
`w_o` gives the load, Eq. 6.28, with `AR = a/L`:

.. math::
    :label: kp_nxy

    N_{xy,crit} = \frac{\pi^2}{2 AR^2 a^2 \tan\alpha} \left( D_{11} (1
    + 6 \tan^2\alpha AR^2 + \tan^4\alpha AR^4) + 2 (D_{12} + 2 D_{66})
    (AR^2 + AR^4 \tan^2\alpha) + D_{22} AR^4 \right)

and minimizing it with respect to `AR` and `\alpha` gives Eqs. 6.29 and
6.30:

.. math::
    :label: kp_nxy_min

    AR = \left( \frac{D_{11}}{D_{11} \tan^4\alpha + 2 (D_{12} + 2
    D_{66}) \tan^2\alpha + D_{22}} \right)^{1/4} \\
    3 D_{11} AR^4 \tan^4\alpha + (6 D_{11} AR^2 + 2 (D_{12} + 2 D_{66})
    AR^4) \tan^2\alpha - (D_{11} + 2 (D_{12} + 2 D_{66}) AR^2 + D_{22}
    AR^4) = 0

which :func:`.calc_Nxy_crit` solves by Newton-Raphson iterations in
`\alpha`, starting from `\pi/6`, with `AR` updated from the first equation
at each iteration; ``RuntimeError`` is raised without convergence in
``max_iter`` iterations. The argument ``D16`` is not used, and the sign of
`N_{xy}` does not matter for `D_{16} = D_{26} = 0`.


Combined compression and shear
------------------------------

For a plate simply supported on four edges under compression `N_o = -N_x`
and shear, Kassapoglou [1]_, Sec. 6.5, uses the two-term deflection of
Eq. 6.31, `w = w_1 \sin(\pi x/a) \sin(\pi y/b) + w_2 \sin(2 \pi x/a)
\sin(2 \pi y/b)`, in the energy of Eq. :eq:`kp_energy`. The stationarity
with respect to `w_1` and `w_2` gives the eigenvalue problem `([a] - [b])
\{w_1, w_2\} = 0`, with

.. math::
    :label: kp_combined

    a_{11} = \frac{\pi^4 D_{11} b}{4 a^3} + \frac{\pi^4 D_{12}}{2 a b}
    + \frac{\pi^4 D_{22} a}{4 b^3} + \frac{\pi^4 D_{66}}{a b}, \quad
    a_{12} = a_{21} = -\frac{160 \pi^2 D_{16}}{9 a^2} - \frac{160
    \pi^2 D_{26}}{9 b^2} \\
    a_{22} = \frac{4 \pi^4 D_{11} b}{a^3} + \frac{8 \pi^4 D_{12}}{a b}
    + \frac{4 \pi^4 D_{22} a}{b^3} + \frac{16 \pi^4 D_{66}}{a b} \\
    b_{11} = \frac{\pi^2 N_o b}{4 a}, \quad
    b_{22} = \frac{\pi^2 N_o b}{a}, \quad
    b_{12} = b_{21} = \frac{32 N_{xy}}{9}

whose determinant is a quadratic equation in `N_o`; the smaller root in
absolute value is returned. Two functions solve it:

- :func:`.calc_Nxx_crit_combined_shear`, the closed form of Kassapoglou,
  Eq. 6.34, p. 142, with `D_{16} = D_{26} = 0` and the shear given as the
  ratio `k = N_{xy}/N_x`, such that the shear load grows with the
  compression:

  .. math::
      :label: kp_combined_closed

      N_{xx,crit} = \frac{\pi^2}{a^2} \frac{D_{11} + 2 (D_{12} + 2 D_{66})
      a^2/b^2 + D_{22} a^4/b^4}{2 - \frac{8192}{81} \frac{a^2 k^2}{b^2
      \pi^4}} \left( 5 \pm \sqrt{9 + \frac{65536}{81} \frac{a^2
      k^2}{b^2 \pi^4}} \right)

  which reduces to Eq. :eq:`kp_nxx` with `m = n = 1` for `k = 0`;
- :func:`.calc_Nxx_crit_combined_shear_full`, with `D_{16}` and `D_{26}`
  and a fixed shear load ``Nxy``, positive in tension as in Eq.
  :eq:`kp_energy`. With `D_{16} = D_{26} = 0` only `N_{xy}^2` enters and
  the result is that of the closed form for `N_{xy} = k N_o`; otherwise the
  sign of the shear load changes the buckling load.


Post-buckling: effective width
------------------------------

Beyond buckling, the compression `N_{xx}` concentrates at the supported
edges. The effective width `b_{eff}`, defined by `\int_y N_{xx} \, dy = 2
N_{x,max} b_{eff}`, is given by Kassapoglou [1]_, Eq. 7.15, p. 164:

.. math::
    :label: kp_beff

    b_{eff} = \frac{b}{2 \left( 1 + 2 \left( 1 + \frac{A_{12}}{A_{11}}
    \right) \left( 1 - \frac{P_{cr}}{P_x} \right) \frac{A_{11}}{A_{11}
    + 3 A_{22}} \right)}

computed by :func:`.calc_beff` from the membrane stiffnesses and the ratio of
the buckling load to the applied load. At `P_x = P_{cr}` the whole width is
effective, `b_{eff} = b/2`, and the expression is meant for `P_x \ge
P_{cr}`.


Choice and limits
-----------------

The functions give the first estimate of a design: :func:`.calc_Nxx_crit`
for compression, :func:`.calc_Nxy_crit` for shear of long plates,
:func:`.calc_Nxx_crit_combined_shear` for proportional compression and
shear, :func:`.calc_Nxx_crit_combined_shear_full` when the bending-twisting
coupling or a fixed shear load matters, and :func:`.calc_beff` for the
post-buckling strength. They assume simply supported edges, a symmetric
laminate and thin plates, without transverse shear deformation, which
overestimates the buckling loads of thick plates and sandwiches; clamped
edges, other loads, unsymmetric laminates or shear-deformable plates require
a Ritz or finite element model with the stiffnesses of
:doc:`constitutive_matrices`.


Verification
------------

``tests/test_kassapoglou.py`` compares the functions with the curves of
Figs. 6.4, 6.11, 6.14 and 7.10 of Kassapoglou [1]_, read from the vector
graphics of the book, and the sign convention of ``Nxy`` with the energy of
Eq. :eq:`kp_energy`, whose Hessian, integrated numerically, is singular at
the returned load.


References
----------

.. [1] Kassapoglou, C., Design and Analysis of Composite Structures, Second
   Edition, John Wiley & Sons, 2013. https://doi.org/10.1002/9781118536933
