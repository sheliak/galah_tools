#!/usr/bin/env python

from setuptools import setup

setup(name='galah_tools',
      version='1.0',
      description='Tools for GALAH',
      author='Janez Kos',
      author_email='jkos@usyd.edu.au',
      url='https://www.galah-survey.org',
      py_modules=['galah_tools', 'sclip.sclip'],
      install_requires=['pyfits', 'numpy', 'scipy', 'matplotlib'],
     )

try:
	import psycopg2 as mdb
except:
	print 'You might want to install psycopg2 module. Otherwise you can\'t use sql databases.'

try:
	import  csv
except:
	print 'You might want to install csv module. Otherwise you can only use sql databases and not text-based databases.'

try:
	from pyflann import *
except:
	print 'You might want to install FLANN and pyflann bindings. Otherwise you will not be able to use nearest neighbour search algorithms.'