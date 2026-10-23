rmdir /s /q dist
rmdir /s /q build
rmdir /s /q composites.egg-info
python -m build
twine upload --repository-url https://test.pypi.org/legacy/ dist/*
twine upload dist/*
