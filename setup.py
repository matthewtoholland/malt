from setuptools import setup, find_packages
from Cython.Build import cythonize
from setuptools.extension import Extension
import sys
import os

extensions = [Extension('esp_gen', ['malt/esp_gen.pyx'])]


setup(name='malt',
      packages=['malt'],
      include_package_data=True,
      ext_modules=cythonize(extensions, language_level="3", annotate=True),
      version='2.2.4',
      license='MIT',
      author='Matthew Holland',
      description='A package for converting pdb files into TRIPOS MOL2 files')
