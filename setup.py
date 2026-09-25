import platform
import os
import sys
import inspect
import subprocess
from setuptools import setup, find_packages
from distutils.extension import Extension

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
Operating System :: Microsoft :: Windows
Operating System :: Unix
Operating System :: POSIX :: BSD
Programming Language :: Python :: 3.8
Programming Language :: Python :: 3.9
Programming Language :: Python :: 3.10
Programming Language :: Python :: 3.11
Programming Language :: Python :: 3.12
Programming Language :: Python :: 3.13
Programming Language :: Python :: 3.14
License :: OSI Approved :: BSD License

"""

is_released = True
version = '0.9.2'

fullversion = write_version_py(version, is_released)

data_files = [('', [
        'README.md',
        'LICENSE',
        'composites/version.py',
        ])]

package_data = {
        'composites': ['*.pxd', '*.pyx'],
        '': ['tests/*.*'],
        }

# NOTE a coverage build, see .github/workflows/coverage.yml, requested with
#      CYTHON_TRACE_NOGIL in the environment or with --define CYTHON_TRACE...
trace = ('CYTHON_TRACE_NOGIL' in os.environ.keys()
         or any('CYTHON_TRACE' in arg for arg in sys.argv))

# NOTE flags for speed. No module uses prange, so OpenMP is not needed. GCC
#      and Clang get -O3 explicitly, because the level inherited from the
#      Python build is not guaranteed, and -fno-math-errno, which lets sqrt()
#      compile to a single instruction and changes no result. MSVC is
#      already at its fastest standard-conforming setting with the /O2 and
#      /GL that setuptools passes. Flags that change floating-point results,
#      such as /fp:fast or -ffast-math, and flags that tie a wheel to the CPU
#      that built it, such as -march=native, are deliberately left out
define_macros = []
if platform.system() == 'Windows':
    compile_args = ['/O2']
    link_args = []
elif platform.system() == 'Linux':
    compile_args = ['-O3', '-fno-math-errno']
    link_args = ['-static-libgcc', '-static-libstdc++']
else: # MAC-OS
    compile_args = ['-O3', '-fno-math-errno']
    link_args = []

if trace:
    # NOTE unoptimized, so that every traced line maps to code. Since Python
    #      3.12 Cython traces through sys.monitoring by default, which the
    #      Cython.Coverage plugin cannot follow, hence the legacy tracing
    if os.name == 'nt': # Windows
        compile_args = ['/Od']
    else: # MAC-OS or Linux
        compile_args = ['-O0']
    link_args = []
    define_macros = [('CYTHON_TRACE_NOGIL', '1'),
                     ('CYTHON_USE_SYS_MONITORING', '0')]

include_dirs = [
            ]


extensions = [
    Extension('composites.core',
        sources=[
            './composites/core.pyx',
            ],
        include_dirs=include_dirs,
        extra_compile_args=compile_args,
        extra_link_args=link_args,
        define_macros=define_macros,
        language='c++'),

    ]


def generated_with_other_trace_mode(ext):
    r"""Whether the C++ file generated from a Cython source of ``ext`` was
    generated with, or without, line tracing, the opposite of this build

    cythonize regenerates a C++ file only when its source is newer, so
    switching between a coverage build and a normal one would otherwise
    reuse the C++ file of the other mode without notice.
    """
    for source in ext.sources:
        if not source.endswith('.pyx'):
            continue
        cpp = os.path.splitext(source)[0] + '.cpp'
        if os.path.isfile(cpp):
            with open(cpp, encoding='utf-8', errors='ignore') as f:
                if ('__Pyx_TraceLine(' in f.read()) != trace:
                    return True
    return False

# NOTE line tracing only for a coverage build, since the profiling hooks it
#      generates otherwise stay active in every function call
ext_modules = cythonize(extensions,
        compiler_directives={'linetrace': trace},
        language_level='3',
        force=any(generated_with_other_trace_mode(ext) for ext in extensions),
        )

s = setup(
    name = "composites",
    version = fullversion,
    author = "Saullo G. P. Castro",
    author_email = "S.G.P.Castro@tudelft.nl",
    description = ("Methods for analysis and design of composites"),
    long_description = read('README.md'),
    long_description_content_type = 'text/markdown',
    license = "3-Clause BSD",
    keywords = "mechanics composite materials composites shell classical first-order laminated plate theory",
    url = "https://github.com/saullocastro/composites",
    package_data = package_data,
    data_files = data_files,
    classifiers = [_f for _f in CLASSIFIERS.split('\n') if _f],
    install_requires = install_requires,
    ext_modules = ext_modules,
    packages = find_packages(),
)
