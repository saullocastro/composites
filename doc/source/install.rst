Installing composites
=====================

From the distributed packages
-----------------------------

Wheels are published on PyPI for Windows, Linux and macOS, for Python 3.10 and
newer, and installed with::

    python -m pip install composites

NumPy is the only run-time requirement.


From the source code
--------------------

The module :mod:`composites.core` is written in Cython and compiled to C++,
such that installing from the source code requires a C++ compiler, e.g. MSVC
on Windows, GCC on Linux or Clang on macOS. Cython is installed
automatically as a build requirement, see ``pyproject.toml``, which also
holds the package metadata and the version::

    git clone https://github.com/saullocastro/composites
    cd composites
    python -m pip install .

For development, the extension can be built in place instead, such that the
tests and the documentation use the source tree when run from its root
directory; ``setup.py`` also writes ``composites/version.py`` with the
version of ``pyproject.toml``::

    python setup.py build_ext --inplace
    python -m pytest tests

The documentation is built with the requirements of
``ci_requirements_doc.txt``, the notebooks also requiring pandoc::

    python -m sphinx -b html doc/source doc/build/html


In the browser, with Pyodide
----------------------------

A WebAssembly wheel is published on PyPI for Pyodide 314 (CPython 3.14),
installed in the browser with ``micropip``, see :doc:`pyodide`::

    import micropip
    await micropip.install('composites')


Using composites from Cython
----------------------------

The declarations of the classes and functions of :mod:`composites.core` are
distributed with the package in ``core.pxd``, such that Cython code can use
them with, e.g.::

    from composites.core cimport Laminate

The function :func:`composites.get_include` returns the directory of the
installed package, which contains ``core.pxd``.
