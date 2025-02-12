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

* A, B, D, E plane-stress stiffness matrices for plates
    - A, B, D, for classical plate theory (CLT, or CLPT)
    - E for first-order shear deformation theory (FSDT)

* Material invariants, trace-normalized or not

* Lamination parameters based on material invariants

* Stiffness matrices (ABDE) based on lamination parameters


Citing this repository
----------------------

Castro, SGP. Methods for analysis and design of composites (Version
0.7.3) [Computer software]. 2025. https://doi.org/10.5281/zenodo.2871782

Bibtex :
    
    @misc{composites2025,
        author = {Castro, Saullo G. P.},
        doi = {10.5281/zenodo.2871782},
        title = {{Methods for analysis and design of composites}}
        }

Documentation
-------------

The documentation is available on: https://saullocastro.github.io/composites.


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


License
-------
Distrubuted under the 3-Clause BSD license
(https://raw.github.com/saullocastro/composites/master/LICENSE).

Contact: S.G.P.Castro@tudelft.nl.

