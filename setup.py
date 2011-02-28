#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup

classifiers = """\
Development Status :: 5 - Production/Stable
Environment :: Console
Intended Audience :: Science/Research
Intended Audience :: Developers
Intended Audience :: Education
License :: OSI Approved :: MIT License
Operating System :: OS Independent
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Education
Topic :: Software Development :: Libraries :: Python Modules
"""

setup(name             = 'airse',
      version          = '0.0.1',
      author           = 'Filipe Fernandes',
      author_email     = 'ocefpaf@gmail.com',
      maintainer       = 'Filipe Fernandes',
      maintainer_email = 'ocefpaf@gmail.com',
      url              = 'http://ocefpaf.tiddlyspot.com/#python-airsea',
      description      = 'AirSea Libray for Python',
      long_description = """\
This module is a translation of the original AIRSEA-2.0 MATLAB toolbox routines for calculating the properties of airsea fluxes.
""",
      download_url     = 'http://pypi.python.org/packages/source/a/airsea/',
      packages         = ['airsea'], #FIXME
      classifiers      = filter(None, classifiers.split("\n")),
      platforms        = 'any',
      license          = 'MIT',
      keywords         = 'oceanography  airsea',
    )
