from __future__ import absolute_import

import os
import inspect
import subprocess
from setuptools import setup, find_packages


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
        f.write(b'__version__ = "%s"\n' % fullversion.encode())
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
        "coveralls",
        "setuptools-git-version",
        ]

CLASSIFIERS = """\

Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
Intended Audience :: Education
Topic :: Scientific/Engineering :: Mathematics
License :: OSI Approved :: BSD License
Operating System :: Microsoft :: Windows
Programming Language :: Python :: 2.7
Programming Language :: Python :: 3.5
Programming Language :: Python :: 3.6
Operating System :: Unix

"""

is_released = True
version = '0.2.1'

fullversion = write_version_py(version, is_released)

data_files = [('', [
        'README.md',
        'LICENSE',
        'composites/version.py',
        ])]

package_data = {
        '': ['tests/*.*'],
        }

s = setup(
    name = "composites",
    version = fullversion,
    author = "Saullo G. P. Castro",
    author_email = "castrosaullo@gmail.com",
    description = ("Methods to calculate composite material properties"),
    license = "BSD",
    keywords = "mechanics composite materials composites shell classical laminated plate theory",
    url = "https://github.com/compmech/composites",
    packages=find_packages(),
    package_data=package_data,
    data_files=data_files,
    long_description=read('README.md'),
    classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
    install_requires=install_requires,
)
