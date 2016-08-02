#from setuptools import setup
from distutils.core import setup

ext_modules = []
include_dirs = []

try:
    from Cython.Build import cythonize
    import numpy
    ext_modules = cythonize("tde/strain_functions_cy/strain_functions_cy.pyx")
    include_dirs = [numpy.get_include()]
except ImportError:
    pass

setup(
    name = 'tde',
    version = '0.2.0',
    packages = ['tde'],
    ext_modules = ext_modules,
    include_dirs = include_dirs
)
