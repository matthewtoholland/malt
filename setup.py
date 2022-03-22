from setuptools import setup, find_packages
from setuptools.extension import Extension
import sys
import os



setup(name='malt',
      packages=['malt'],
      include_package_data=True,
      package_data = {'malt':['Data/*']},
      version='2.3.0',
      license='MIT',
      author='Matthew Holland',
      description='A package for converting pdb files into TRIPOS MOL2 files')
