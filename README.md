Github Actions status:
[![pytest](https://github.com/saullocastro/composites/actions/workflows/pytest.yml/badge.svg)](https://github.com/saullocastro/composites/actions/workflows/pytest.yml)
[![Documentation](https://github.com/saullocastro/composites/actions/workflows/auto_doc.yml/badge.svg)](https://github.com/saullocastro/composites/actions/workflows/auto_doc.yml)
[![Deploy](https://github.com/saullocastro/composites/actions/workflows/pythonpublish.yml/badge.svg)](https://github.com/saullocastro/composites/actions/workflows/pythonpublish.yml)

Coverage status:

[![codecov](https://github.com/saullocastro/composites/actions/workflows/coverage.yml/badge.svg)](https://github.com/saullocastro/composites/actions/workflows/coverage.yml)
[![Codecov Status](https://codecov.io/gh/saullocastro/composites/branch/master/graph/badge.svg?token=KD9D8G8D2P)](https://codecov.io/gh/saullocastro/composites)


Methods for analysis and design of composites
=============================================

High-performance module to calculate properties of laminated composite
materials. Usually, this module is used to calculate:

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

* Saving and loading of laminates, plies, materials and lamination parameters
  in JSON


Documentation
-------------

The documentation is available on: https://saullocastro.github.io/composites.


Running in the browser with Pyodide
-----------------------------------

From version 0.9.21, a WebAssembly wheel is published on PyPI for Pyodide 314
(CPython 3.14), such that ``composites`` runs in the browser, e.g. in
JupyterLite or in a web page with Pyodide::

    import micropip
    await micropip.install('composites')

    from composites import laminated_plate, to_json
    lam = laminated_plate([0, 45, -45, 90], plyt=0.125e-3,
                          laminaprop=(142e9, 8.7e9, 0.28, 5.1e9, 5.1e9, 3.2e9))
    print(lam.ABD)
    s = to_json(lam)  # strict JSON, readable with JSON.parse in JavaScript


Citing this repository
----------------------

Castro, SGP. Methods for analysis and design of composites (Version 0.9.21) [Computer software]. 2026. https://doi.org/10.5281/zenodo.2871782

Bibtex :
    
    @misc{composites2026,
        author = {Castro, Saullo G. P.},
        doi = {10.5281/zenodo.2871782},
        title = {{Methods for analysis and design of composites (Version 0.9.21)}},
        year = 2026
        }


History
-------

- version 0.1.0, from sub-module of compmech 0.7.2
- version 0.2.2, from sub-module of meshless 0.1.19
- version 0.2.3 onwards: independent of previous packages
- version 0.3.0 onwards: with fast Cython version, not compatible with previous versions
- version 0.4.0 onwards: fast Cython and cimportable by other packages, full
  compatibility with finite element mass matrices of plates and shells,
  supporting laminated plates with materials of different densities
- version 0.5.4 onwards: verified lamination parameters, analytical gradients
  of Aij, Bij, Dij with respect to lamination parameters, supportting MAC-OS
- version 0.5.17 onwards: installing with pip
- version 0.6.0 onwards: cibuildwheel to distribute for Linux
- version 0.7.0 onwards: added Kassapoglou's module
- version 0.8.0 onwards: support for Third-order Shear Deformation Theory (TSDT)
- version 0.9.1 onwards: A44, A45, A55 with the shear correction already applied (Rohwer, 1988), transverse shear stress recovery, picklable Laminate and GradABD, improved Cython and compiler flags (0.9.2), see CHANGELOG.md
- version 0.9.12 onwards: shear correction methods of Whitney (1973), Chow (1971), Birman and Bert (2002) and Yang, Norris and Stavsky (1966), a posteriori transverse shear stresses and energies (Noor and Peters, 1989), Atrans, Dtrans and Ftrans deprecated in favour of Ats, Dts and Fts, see CHANGELOG.md
- version 0.9.21 onwards: saving and loading in JSON, support for Pyodide (WebAssembly) with wheels on PyPI, see CHANGELOG.md


License
-------
Distrubuted under the 3-Clause BSD license
(https://raw.github.com/saullocastro/composites/master/LICENSE).

Contact: S.G.P.Castro@tudelft.nl.

