Matrices of the third-order shear deformation theory
====================================================

The TSDT uses the in-plane matrices ``A``, ``B``, ``D``, ``E``, ``F``, ``H``
and the uncorrected transverse shear matrices ``Abar_ts``, ``Dts``, ``Fts``,
see :doc:`constitutive_matrices`. The test below assembles the ``13 x 13``
constitutive matrix of the generalized strains `\{\varepsilon^{(0)},
\varepsilon^{(1)}, \varepsilon^{(3)}, \gamma^{(0)}, \gamma^{(2)}\}` and
checks that it gives the strain energy per unit area obtained by integrating
the ply strain energy density through the thickness:

.. literalinclude:: ../../tests/test_core.py
    :pyobject: test_example_tsdt_matrices
