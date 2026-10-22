Saving and loading in JSON
==========================

A :class:`.Laminate` saved with :func:`.save_json` is loaded identical by
:func:`.load_json`, from a file name or a file object, see
:doc:`saving_loading`. The function ``make_laminate`` of the tests creates a
laminate with :func:`.laminated_plate`:

.. literalinclude:: ../../tests/test_json.py
    :pyobject: test_save_load

:func:`.to_dict` and :func:`.from_dict` embed a laminate in another JSON
document:

.. literalinclude:: ../../tests/test_json.py
    :pyobject: test_example_embed_in_json_document

A file saved before the mass integrals ``intrhoz3``, ``intrhoz4`` and
``intrhoz6`` existed still loads, the missing attributes being zero until
:meth:`.Laminate.calc_constitutive_matrix` is called:

.. literalinclude:: ../../tests/test_json.py
    :pyobject: test_laminate_without_new_density_moments
