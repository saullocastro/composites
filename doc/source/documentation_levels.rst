How the documentation is organised
==================================

The documentation of ``composites`` is written in four cumulative levels,
each with its own place, so that a reader, or a tool, can go from what a
function does to why it does it without the levels repeating each other.

.. list-table::
    :header-rows: 1

    * - Level
      - Where
      - Content
    * - 1, functional
      - docstrings of the modules, classes, functions and methods of
        ``composites/`` (``.pyx``, ``.pxd`` and ``.py`` files), except their
        ``Notes`` and ``References`` sections; comments in the code;
        docstrings and comments of the tests; the README
      - what each resource does: a one-line summary, an extended summary
        (what, not why), ``Parameters`` with type, unit, shape and
        default, ``Returns``, ``Raises``, ``See Also`` and the
        ``Attributes`` of the Cython classes; short usage ``Examples``
        only on the main entry points, such as :func:`.laminated_plate`
    * - 2, theory
      - the ``Notes`` and ``References`` sections of the docstrings, used
        only for theory
      - the theory that each resource implements, local to it: the
        equations in the notation of the code, the assumptions and limits
        of validity, the meaning of each parameter in that context, and
        the published sources
    * - 3, theory in context
      - the pages of this documentation listed under "Theory and
        conventions", with :doc:`install` and :doc:`pyodide`
      - how the pieces of level 2 fit together, like a book: the
        conventions shared by the modules, the lamina, the constitutive
        matrices, the shear correction, the transverse shear stresses, the
        lamination parameters, the density moments, Kassapoglou's buckling
        formulas and the JSON format, with short usage snippets and links
        to the API instead of repeating level 2
    * - 4, worked examples
      - the notebooks in ``notebooks/``, listed under "Worked examples" and
        "Tutorials"
      - complete applications of the theory with stored outputs: each step
        explained, results printed and plotted, and compared with their
        reference values

The API pages render levels 1 and 2 together. To find the theory of a
resource, read the ``Notes`` of its docstring; to see how a method fits in
the theory of laminated plates, read the pages of level 3, starting from
:doc:`conventions`.


The tests
---------

The tests in ``tests/`` are part of level 1 for their docstrings and
comments, and part of level 2 for their ``Notes`` and ``References``. Each
test module states in its docstring what it verifies and against what:

* a published value, table or figure of the literature, named in the
  docstring and cited in full in its ``References``, e.g. the shear
  correction factors of Rohwer (1988) or the buckling curves of
  Kassapoglou's book;
* a closed form computed in the test, independently of ``composites``,
  e.g. the stiffness matrices of an isotropic plate;
* a property that holds for any input, e.g. a JSON or pickle round trip
  that restores every public attribute, or the equality of an N-double
  laminate with the stack it homogenizes.

The equations and the sources of the reference values are in the ``Notes``
and ``References`` of each test module.
