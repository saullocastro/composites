rmdir /s /q dist
rmdir /s /q build
rmdir /s /q composites.egg-info
python setup.py sdist bdist_wheel
twine upload --repository-url https://test.pypi.org/legacy/ dist/*
twine upload dist/*
