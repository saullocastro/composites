Running composites in the browser with Pyodide
==============================================

A WebAssembly wheel of ``composites`` is published on PyPI for Pyodide 314
(CPython 3.14), tagged ``pyemscripten_2026_0_wasm32`` (PEP 783), such that
``composites`` runs in the browser, e.g. in JupyterLite or in a web page
with `Pyodide <https://pyodide.org>`_, with the same API. The wheel is
built with cibuildwheel and the tests of ``tests/`` run on it in Pyodide.


Installing
----------

In a Pyodide console or a JupyterLite notebook, ``composites`` and its
run-time requirement NumPy are installed with ``micropip``::

    import micropip
    await micropip.install('composites')

after which ``composites`` is imported as usual::

    from composites import laminated_plate

    lam = laminated_plate([0, 45, -45, 90], plyt=0.125e-3,
                          laminaprop=(142e9, 8.7e9, 0.28, 5.1e9, 5.1e9, 3.2e9))
    print(lam.ABD)


Exchanging laminates with JavaScript
------------------------------------

The functions of :mod:`composites.json_io` give strict JSON, see
:doc:`saving_loading`, which JavaScript reads with ``JSON.parse``. In a web
page that has loaded Pyodide, the stiffnesses of a laminate are computed in
Python and read in JavaScript with:

.. code-block:: javascript

    const pyodide = await loadPyodide();
    await pyodide.loadPackage('micropip');
    const s = await pyodide.runPythonAsync(`
    import micropip
    await micropip.install('composites')
    from composites import laminated_plate, to_json
    lam = laminated_plate([0, 45, -45, 90], plyt=0.125e-3,
                          laminaprop=(142e9, 8.7e9, 0.28, 5.1e9, 5.1e9, 3.2e9))
    to_json(lam)
    `);
    const lam = JSON.parse(s);
    console.log(lam.type, lam.data.A11, lam.data.D11);

The object has the keys ``type``, ``format_version``,
``composites_version`` and ``data``, the last one holding the attributes of
the :class:`.Laminate`, e.g. ``A11`` to ``H66`` and the transverse shear
stiffnesses ``A44``, ``A45``, ``A55``. The string ``s`` is read back in
Python with :func:`.from_json`.
