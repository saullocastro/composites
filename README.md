Travis-CI status:

[![Build Status](https://travis-ci.com/saullocastro/composites.svg?branch=master)](https://travis-ci.com/saullocastro/composites)

Github Actions status:

[![Actions Status](https://github.com/saullocastro/composites/workflows/pytest/badge.svg)](https://github.com/saullocastro/composites/actions)

Coverage status:

[![Coverage Status](https://coveralls.io/repos/github/saullocastro/composites/badge.png?branch=master)](https://coveralls.io/github/saullocastro/composites?branch=master)
[![Codecov Status](https://codecov.io/gh/saullocastro/composites/branch/master/graph/badge.svg?token=KD9D8G8D2P)](https://codecov.io/gh/saullocastro/composites)


Methods to calculate composite plate properties
===============================================

High-performance module to calculate properties of laminated composite
materials. Usually, this module is used to calculate:

- A, B, D, E plane-stress stiffness matrices for plates
-- A, B, D, for classical plate theory (CLT, or CLPT)
-- E for first-order shear deformation theory (FSDT)
- Material invariants, trace-normalized or not
- Lamination parameters based on material invariants
- Stiffness matrices (ABDE) based on lamination parameters


Citing this repository
----------------------

Castro, S. G. P. Methods to calculate composite plate properties (Version
0.5.1) [Computer software]. 2022. https://doi.org/10.5281/zenodo.2871782

Bibtex :
    
    @misc{composites2022,
        author = {Castro, Saullo G. P.},
        doi = {10.5281/zenodo.2871782},
        title = {{Methods to calculate composite plate properties (Version 0.5.1) [Computer software]. 2022}}
        }

Documentation
-------------

The documentation is available on: https://saullocastro.github.io/composites.


History
-------

- version 0.1.0, from sub-module of compmech 0.7.2
- version 0.2.2, from sub-module of meshless 0.1.19
- version 0.2.3 onwards, independent of previous packages
- version 0.3.0 onwards, with fast Cython version, not compatible with previous versions
- version 0.4.0 onwards, fast Cython and cimportable by other packages, full
  compatibility with finite element mass matrices of plates and shells,
  supporting laminated plates with materials of different densities
- version 0.5.0 onwards, verified lamination parameters, analytical gradients
  of Aij, Bij, Dij with respect to lamination parameters
  

License
-------
Distrubuted under the 3-Clause BSD license
(https://raw.github.com/saullocastro/composites/master/LICENSE).

Contact: S.G.P.Castro@tudelft.nl.

