#!/usr/bin/env python

from distutils.core import *
from distutils import sysconfig
import os.path

# Get numpy include directory (works across versions)
import numpy
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

LIBCCL_LIBRARY_PATH = "/home/damonge/lib"

_nmtlib = Extension("_nmtlib",
                    ["pymaster/namaster.i"],
                    libraries = ['nmt','sharp','fftpack','c_utils','chealpix','cfitsio','gsl','gslcblas','m','gomp'],
                    include_dirs = [numpy_include, "../src/"],
                    extra_compile_args=['-O4', '-fopenmp',],
                    )

setup(name = "pymaster",
      description = "Library for pseudo-Cl computation",
      author = "David Alonso",
      version = "0.1",
      packages = ['pymaster'],
      ext_modules = [_nmtlib],
      )
