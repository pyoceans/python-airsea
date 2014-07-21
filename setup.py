#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import re
from setuptools import setup

VERSIONFILE = "airsea/__init__.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))


def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.txt', 'CHANGES.txt')
LICENSE = read('LICENSE.txt')


classifiers = """\
Development Status :: 1 - Planning
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

config = dict(name='airsea',
              version=verstr,
              packages=['airsea'],
              test_suite='tests',
              use_2to3=True,
              license=LICENSE,
              long_description=long_description,
              classifiers=filter(None, classifiers.split("\n")),
              description='AirSea Libray for Python',
              author='Filipe Fernandes',
              author_email='ocefpaf@gmail.com',
              maintainer='Filipe Fernandes',
              maintainer_email='ocefpaf@gmail.com',
              #url='http://pypi.python.org/pypi/airsea/',
              #download_url='%s/a/airsea/airsea-%s.tar.gz' % (source, verstr),
              platforms='any',
              keywords=['oceanography', 'data analysis', 'air-sea'],
              install_requires=install_requires)

setup(**config)
