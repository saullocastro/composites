Documentation for the ``composites`` module
===========================================

High-performance methods for analysis and design of composites.
With the ``composites`` module, you are able to calculate:

* A, B, D, E plane-stress stiffness matrices for plates
    - A, B, D, for classical plate theory (CLT, or CLPT)
    - E for first-order shear deformation theory (FSDT)

* Material invariants, trace-normalized or not

* Lamination parameters based on material invariants

* Stiffness matrices (ABDE) based on lamination parameters

* Based on Kassapoglou's book, local buckling under compression, shear, and
  post-buckling 


Code repository
---------------

https://github.com/saullocastro/composites


Citing this library
-------------------

Castro, S. G. P. Methods for analysis and design of composites (Version 0.7.1) [Computer software]. 2024. https://doi.org/10.5281/zenodo.2871782

Bibtex :
    
    @misc{composites2024,
        author = {Castro, Saullo G. P.},
        doi = {10.5281/zenodo.2871782},
        title = {{Methods for analysis and design of composites (Version 0.7.1) [Computer software]. 2024}}
        }


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

