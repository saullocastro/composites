Travis-CI status:

[![Build Status](https://travis-ci.com/saullocastro/composites.svg?branch=master)](https://travis-ci.com/saullocastro/composites)

Github Actions status:

[![Actions Status](https://github.com/saullocastro/composites/workflows/pytest/badge.svg)](https://github.com/saullocastro/composites/actions)

Coverage status:

[![Coverage Status](https://coveralls.io/repos/github/saullocastro/composites/badge.svg?branch=master)](https://coveralls.io/github/saullocastro/composites?branch=master)


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

Repository created from sub-modules of "compmech" and "meshless" repositories.
The first adopted version number came from the following reasoning:

- from compmech 0.7.2: composites 0.1.0
- from meshless 0.1.19: composites 0.2.2;
- version 0.2.3 onwards is independent of previous packages
- version 0.3.0 onwards with fast Cython version, not compatible with previous versions


License
-------
Distrubuted in the 2-Clause BSD license (https://raw.github.com/saullocastro/composites/master/LICENSE).

Contact: castrosaullo@gmail.com

