Saving and loading (JSON)
=========================

The objects :class:`.MatLamina`, :class:`.Lamina`, :class:`.Laminate`,
:class:`.LaminationParameters` and :class:`.GradABD` are saved to and loaded
from JSON by the functions of :mod:`composites.json_io`, all available from
``composites``:

=============================================  ==============================================
:func:`.save_json`, :func:`.load_json`         file, given by its name or as a file object
:func:`.to_json`, :func:`.from_json`           JSON string
:func:`.to_dict`, :func:`.from_dict`           dictionary of JSON types, to embed in another
                                               JSON document
=============================================  ==============================================

::

    from composites import laminated_plate, save_json, load_json

    lam = laminated_plate([0, 45, -45, 90], plyt=0.125e-3,
                          laminaprop=(142e9, 8.7e9, 0.28, 5.1e9, 5.1e9, 3.2e9))
    save_json(lam, 'laminate.json')
    lam2 = load_json('laminate.json')

All the public attributes are stored, including the derived ones, such as the
stiffnesses, the mass integrals and the ply stiffnesses, such that the loaded
object is identical to the saved one without any recomputation. JSON files
are plain text, independent of the Python version and safe to load from
untrusted sources, unlike pickle, which remains supported for
``copy.deepcopy`` and ``multiprocessing`` but is not recommended for
files. The output is strict JSON, also readable with ``JSON.parse`` in
JavaScript, e.g. in the browser with Pyodide, see :doc:`install`.


Embedding in other documents
----------------------------

:func:`.to_dict` returns a dictionary containing only ``dict``, ``list``,
``str``, ``int``, ``float`` and ``None``, such that a laminate can be one
entry of the input file of a model, read back with :func:`.from_dict`, see
:doc:`ex_json`.


Format version
--------------

Each saved object carries ``'type'``, ``'format_version'`` and
``'composites_version'``, next to the attributes in ``'data'``. The current
``format_version`` is 2, which stores the mass integrals ``intrhoz3``,
``intrhoz4`` and ``intrhoz6`` of a :class:`.Laminate`. :func:`.from_dict`
reads only the ``format_version`` of the installed ``composites`` and raises
a ``ValueError`` for an older or a newer one, e.g. for the files of format 1
saved by version 0.9.21, which must be created again. Unknown keys also
raise a ``ValueError``.
