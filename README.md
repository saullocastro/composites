Travis-CI status:

[![Build Status](https://travis-ci.com/saullocastro/composites.svg?branch=master)](https://travis-ci.com/saullocastro/composites)

Github Actions status:

[![Actions Status](https://github.com/saullocastro/composites/workflows/pytest/badge.svg)](https://github.com/saullocastro/composites/actions)

Coverage status:

[![Coverage Status](https://coveralls.io/repos/github/saullocastro/composites/badge.svg?branch=master)](https://coveralls.io/github/saullocastro/composites?branch=master)
[![Codecov Status](https://codecov.io/gh/saullocastro/composites/branch/master/graph/badge.svg?token=KD9D8G8D2P)](https://codecov.io/gh/saullocastro/composites)


Methods to calculate structural properties of plates
====================================================

Usually, this module is used to calculate:

- A, B, D, E stiffness matrices for plates
-- A, B, D, for classical plate theory (CLT, or CLPT)
-- E for first-order shear deformation theory (FSDT)

- Lamination parameters based on material invariants

Documentation
===

The documentation is available on: https://saullocastro.github.io/composites/


History
===

- version 0.1.0, from sub-module of compmech 0.7.2
- version 0.2.2, from sub-module of meshless 0.1.19
- version 0.2.3 onwards, independent of previous packages
- version 0.3.0 onwards, with fast Cython version, not compatible with previous versions
- version 0.4.0 onwards, fast Cython and cimportable by other packages, full
  compatibility with finite element mass matrices of plates and shells,
  supporting laminated plates with materials of different densities
  

License
-------
Distrubuted in the 2-Clause BSD license (https://raw.github.com/saullocastro/composites/master/LICENSE).

Contact: castrosaullo@gmail.com

