#!/usr/bin/env python

from setuptools import setup

setup(name='galah_tools',
      version='1.0',
      description='Tools for GALAH',
      author='Janez Kos',
      author_email='jkos@usyd.edu.au',
      url='https://www.galah-survey.org',
      py_modules=['galah_tools'],
      install_requires=['pyfits', 'numpy', 'psycopg2', 'csv', 'os', 'fnmatch', 'pickle', 'pyflann', 'scipy', 'matplotlib'],
     )