#!/usr/bin/env python3

import os, sysconfig
from glob import glob
from setuptools import setup, find_packages, Extension

# https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html#distributing-cython-modules
def no_cythonize(extensions, **_ignore):
  for extension in extensions:
    sources = []
    for sfile in extension.sources:
      path, ext = os.path.splitext(sfile)
      if ext in (".pyx", ".py"):
        if extension.language == "c++":
            ext = ".cpp"
        else:
            ext = ".c"
        sfile = path + ext
      sources.append(sfile)
    extension.sources[:] = sources
  return extensions

athena_dir = '../src/'
srcs = ['parameter_input.cpp', 
        'globals.cpp',
        'outputs/io_wrapper.cpp',
        'utils/utils.cpp']
srcs = [athena_dir + x for x in srcs]
srcs += glob(athena_dir + 'radiation/*.cpp')

extra_compile_args = sysconfig.get_config_var('CFLAGS').split()
extra_compile_args += ['-std=c++11']

pyharp = Extension('athenapp.pyharp',
  define_macros = [('MAJOR_VERSION', '0'),
                   ('MINOR_VERSION', '1')],
  include_dirs = [],
  libraries = ['netcdf'],
  #library_dirs = ['/usr/local/lib'],
  sources = ['athenapp/pyharp/pyharp.pyx'] + srcs,
  extra_compile_args = extra_compile_args,
  language = 'c++'
  )

pyathena = Extension('athenapp.pyathena',
  define_macros = [('MAJOR_VERSION', '0'),
                   ('MINOR_VERSION', '1')],
  include_dirs = [],
  libraries = [],
  #library_dirs = ['/usr/local/lib'],
  sources = ['athenapp/pyathena/pyathena.pyx'],
  extra_compile_args = extra_compile_args,
  language = 'c++'
  )

extensions = [pyathena, pyharp]

CYTHONIZE = bool(int(os.getenv("CYTHONIZE", 0)))

if CYTHONIZE:
  compiler_directives = {"language_level": 2, "embedsignature": True}
  extensions = Cython.Build.cythonize(extensions, compiler_directives = compiler_directives)
else:
  extensions = no_cythonize(extensions)

with open("requirements.txt") as fp:
  install_requires = fp.read().strip().split("\n")

with open("requirements-dev.txt") as fp:
  dev_requires = fp.read().strip().split("\n")
  
with open("README", "r") as fh:
  long_description = fh.read()

setup(
  name="snapy", # Replace with your own username
  version="0.0.6",
  author="Cheng Li",
  author_email="chengcli@umich.edu",
  description="python extension for SNAP model",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://snap.chengcli.io",
  #ext_modules = extensions,
  ext_modules = None,
  install_requires = install_requires,
  extras_require = {
    "dev": dev_requires,
    "docs": ["sphinx", "sphinx-rtd-theme"]
  },
  packages = ['snapy', 'snapy.athena', 'snapy.harp'],
  #packages = find_packages(),
  classifiers = [
      "Programming Language :: Python :: 3",
      "License :: OSI Approved :: MIT License",
      "Operating System :: OS Independent",
  ],
  python_requires='>=3.0',
)
