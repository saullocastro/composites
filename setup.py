from __future__ import absolute_import

import os
import inspect
import subprocess
from setuptools import setup, find_packages
from distutils.extension import Extension

import numpy as np
from Cython.Build import cythonize

def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout=subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        git_revision = out.strip().decode('ascii')
    except OSError:
        git_revision = "Unknown"

    return git_revision


def get_version_info(version, is_released):
    fullversion = version
    if not is_released:
        git_revision = git_version()
        fullversion += '.dev0+' + git_revision[:7]
    return fullversion


def write_version_py(version, is_released, filename='composites/version.py'):
    fullversion = get_version_info(version, is_released)
    with open("./composites/version.py", "wb") as f:
        f.write(('__version__ = "%s"\n' % fullversion).encode())
    return fullversion


# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    setupdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    return open(os.path.join(setupdir, fname)).read()


#_____________________________________________________________________________

install_requires = [
        "numpy",
        ]

CLASSIFIERS = """\

Development Status :: 5 - Production/Stable
Intended Audience :: Science/Research
Intended Audience :: Developers
Intended Audience :: Education
Intended Audience :: End Users/Desktop
Topic :: Scientific/Engineering
Topic :: Education
Topic :: Software Development
Topic :: Software Development :: Libraries :: Python Modules
Operating System :: POSIX :: BSD
Operating System :: Microsoft :: Windows
Operating System :: Unix
Programming Language :: Python :: 3.7
Programming Language :: Python :: 3.8
Programming Language :: Python :: 3.9
Programming Language :: Python :: 3.10
License :: OSI Approved :: BSD License

"""

is_released = True
version = '0.4.22'

fullversion = write_version_py(version, is_released)

data_files = [('', [
        'README.md',
        'LICENSE',
        'composites/version.py',
        ])]

package_data = {
        'composites': ['*.pxd'],
        '': ['tests/*.*'],
        }

if os.name == 'nt':
    compile_args = ['/openmp', '/O2']
    link_args = []
else:
    compile_args = ['-fopenmp']
    link_args = ['-fopenmp']
include_dirs = [
            np.get_include(),
            ]

extensions = [
    Extension('composites.core',
        sources=[
            './composites/core.pyx',
            ],
        include_dirs=include_dirs,
        extra_compile_args=compile_args,
        extra_link_args=link_args,
        language='c++'),

    ]
ext_modules = cythonize(extensions,
        compiler_directives={'linetrace': True},
        language_level='3',
        )

s = setup(
    name = "composites",
    version = fullversion,
    author = "Saullo G. P. Castro",
    author_email = "S.G.P.Castro@tudelft.nl",
    description = ("Methods to calculate properties of laminated composite materials"),
    long_description = read('README.md'),
    long_description_content_type = 'text/markdown',
    license = "BSD",
    keywords = "mechanics composite materials composites shell classical first-order laminated plate theory",
    url = "https://github.com/saullocastro/composites",
    package_data = package_data,
    data_files = data_files,
    classifiers = [_f for _f in CLASSIFIERS.split('\n') if _f],
    install_requires = install_requires,
    ext_modules = ext_modules,
    packages = find_packages(),
)
