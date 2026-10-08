Documentation for the ``composites`` module
===========================================

High-performance methods for analysis and design of composites.
With the ``composites`` module, you are able to calculate:

* A, B, D, E, F, H plane-stress stiffness matrices
    - A, B, D, for classical plate theory (CLT, or CLPT)
    - A44, A45, A55 for first-order shear deformation theory (FSDT), with
      the shear correction already applied, by default with the
      equilibrium approach of Rohwer (1988); also available: Vlachoutsis
      (1992), Whitney (1973), Chow (1971), Birman and Bert (2002), the
      thickness-shear frequency of Yang, Norris and Stavsky (1966), and a
      constant 5/6
    - E, F, H for third-order shear deformation theory (TSDT)

* Transverse shear stresses through the thickness, from the shear forces or
  recovered a posteriori from the strain gradients of an FSDT solution (Noor
  and Peters, 1989), with the corresponding transverse shear strain energies

* Material invariants, trace-normalized or not

* Lamination parameters based on material invariants

* Stiffness matrices (ABD) based on lamination parameters

* Based on Kassapoglou's book, local buckling under compression, shear, and
  post-buckling 

* Saving and loading of laminates, plies, materials and lamination parameters
  in JSON


Running in the browser with Pyodide
-----------------------------------

From version 0.9.20, a WebAssembly wheel is published on PyPI for Pyodide 314
(CPython 3.14), such that ``composites`` runs in the browser, e.g. in
JupyterLite or in a web page with Pyodide::

    import micropip
    await micropip.install('composites')

    from composites import laminated_plate, to_json
    lam = laminated_plate([0, 45, -45, 90], plyt=0.125e-3,
                          laminaprop=(142e9, 8.7e9, 0.28, 5.1e9, 5.1e9, 3.2e9))
    print(lam.ABD)
    s = to_json(lam)  # strict JSON, readable with JSON.parse in JavaScript


Code repository
---------------

https://github.com/saullocastro/composites


Citing this library
-------------------

Castro, S. G. P. Methods for analysis and design of composites (Version 0.9.20) [Computer software]. 2026. https://doi.org/10.5281/zenodo.2871782

Bibtex::
    
    @misc{composites2026,
        author = {Castro, Saullo G. P.},
        doi = {10.5281/zenodo.2871782},
        title = {{Methods for analysis and design of composites (Version 0.9.20)}},
        year = 2026
        }

Tutorials
---------

.. toctree::
    :maxdepth: 2

    tutorials.rst


composites API
--------------

.. toctree::
    :maxdepth: 1

    api.rst


License
-------

.. literalinclude:: ../../LICENSE
    :encoding: latin-1


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

