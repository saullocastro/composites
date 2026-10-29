Documentation for the ``composites`` module
===========================================

High-performance methods for analysis and design of composites.
With the ``composites`` module, you are able to calculate:

* A, B, D, E, F, H plane-stress stiffness matrices, see
  :doc:`constitutive_matrices`

  - A, B, D, for classical plate theory (CLT, or CLPT)
  - A44, A45, A55 for first-order shear deformation theory (FSDT), with
    the shear correction already applied, by default with the
    equilibrium approach of Rohwer (1988); also available: Vlachoutsis
    (1992), Whitney (1973), Chow (1971), Birman and Bert (2002), the
    thickness-shear frequency of Yang, Norris and Stavsky (1966), and a
    constant 5/6, see :doc:`shear_correction`
  - E, F, H for third-order shear deformation theory (TSDT)

* The moments of the density through the thickness used in the mass
  matrices of the CLPT, FSDT and TSDT, see :doc:`mass_integrals`

* Transverse shear stresses through the thickness, from the shear forces or
  recovered a posteriori from the strain gradients of an FSDT solution (Noor
  and Peters, 1989), with the corresponding transverse shear strain energies

* Material invariants, trace-normalized or not

* Lamination parameters based on material invariants

* Stiffness matrices (ABD) based on lamination parameters

* Based on Kassapoglou's book, local buckling under compression, shear, and
  post-buckling

* Saving and loading of laminates, plies, materials and lamination parameters
  in JSON, see :doc:`saving_loading`

``composites`` also runs in the browser with Pyodide, see :doc:`install`.


Code repository
---------------

https://github.com/saullocastro/composites


Citing this library
-------------------

Castro, S. G. P. Methods for analysis and design of composites (Version 0.10.2) [Computer software]. 2026. https://doi.org/10.5281/zenodo.2871782

Bibtex::

    @misc{composites2026,
        author = {Castro, Saullo G. P.},
        doi = {10.5281/zenodo.2871782},
        title = {{Methods for analysis and design of composites (Version 0.10.2)}},
        year = 2026
        }


Theory and conventions
----------------------

.. toctree::
    :maxdepth: 1

    conventions.rst
    constitutive_matrices.rst
    mass_integrals.rst
    shear_correction.rst
    saving_loading.rst


Usage examples
--------------

.. toctree::
    :maxdepth: 1

    ex_laminated_plate.rst
    ex_isotropic_offset.rst
    ex_tsdt.rst
    ex_transverse_shear.rst
    ex_json.rst


Tutorials
---------

.. toctree::
    :maxdepth: 2

    tutorials.rst


composites API
--------------

.. toctree::
    :maxdepth: 2

    api.rst


Installing composites
---------------------

.. toctree::
    :maxdepth: 1

    install.rst


Changelog
---------

.. toctree::
    :maxdepth: 1

    changelog.rst


License
-------

.. literalinclude:: ../../LICENSE
    :encoding: latin-1


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
