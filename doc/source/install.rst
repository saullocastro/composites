Installing composites
=====================

From the distributed packages
-----------------------------

Wheels are published on PyPI for Windows, Linux and macOS, for Python 3.9 and
newer, and installed with::

    python -m pip install composites

NumPy is the only run-time requirement.


From the source code
--------------------

The module :mod:`composites.core` is written in Cython and compiled to C++,
such that installing from the source code requires a C++ compiler, e.g. MSVC
on Windows, GCC on Linux or Clang on macOS. Cython is installed
automatically as a build requirement, see ``pyproject.toml``::

    git clone https://github.com/saullocastro/composites
    cd composites
    python -m pip install .

For development, the extension can be built in place instead, such that the
tests and the documentation use the source tree when run from its root
directory::

    python setup.py build_ext --inplace
    python -m pytest tests

The documentation is built with the requirements of
``ci_requirements_doc.txt``, the notebooks also requiring pandoc::

    python -m sphinx -b html doc/source doc/build/html


In the browser, with Pyodide
----------------------------

From version 0.9.21, a WebAssembly wheel is published on PyPI for Pyodide 314
(CPython 3.14), such that ``composites`` runs in the browser, e.g. in
JupyterLite or in a web page with `Pyodide <https://pyodide.org>`_::

    import micropip
    await micropip.install('composites')

    from composites import laminated_plate, to_json
    lam = laminated_plate([0, 45, -45, 90], plyt=0.125e-3,
                          laminaprop=(142e9, 8.7e9, 0.28, 5.1e9, 5.1e9, 3.2e9))
    print(lam.ABD)
    s = to_json(lam)  # strict JSON, readable with JSON.parse in JavaScript


Using composites from Cython
----------------------------

The declarations of the classes and functions of :mod:`composites.core` are
distributed with the package in ``core.pxd``, such that Cython code can use
them with, e.g.::

    from composites.core cimport Laminate

The function :func:`composites.get_include` returns the directory of the
installed package, which contains ``core.pxd``.
