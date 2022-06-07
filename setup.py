#!/usr/bin/env python

from setuptools import setup

setup(name='galah_tools',
      version='2.0',
      description='Tools for GALAH',
      author='Janez Kos',
      author_email='janez.kos@gmail.com',
      url='https://www.galah-survey.org',
      py_modules=['galah_tools', 'sclip.sclip'],
      install_requires=['pyfits', 'numpy', 'scipy', 'matplotlib'],
     )

try:
	import  csv
except:
	print 'You might want to install csv module. Otherwise you can only use sql databases and not text-based databases.'
