#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'mnowotka'

import sys

try:
    from setuptools import setup
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()
    from setuptools import setup
setup(
    name='standardiser',
    version='0.1.7',
    author='Francis Atkinson',
    author_email='francis@ebi.ac.uk',
    description='Provides a simple way of standardising molecules as a prelude to e.g. molecular modelling exercises.',
    url='https://www.ebi.ac.uk/chembldb/index.php/ws',
    license='Apache License, Version 2.0',
    entry_points={
        'console_scripts': [
            'standardiser=standardiser.bin.standardiser:main']},
    packages=['standardiser',
              'standardiser.bin'],
    long_description=open('ReadMe.txt').read(),
    package_data={
        'standardiser': ['bin/*', 'data/*', 'docs/*'],
        },
    classifiers=['Development Status :: 2 - Pre-Alpha',
                 'Environment :: Console',
                 'Intended Audience :: Science/Research',
                 'License :: OSI Approved :: Apache Software License',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python',
                 'Topic :: Scientific/Engineering :: Chemistry'],
    zip_safe=False,
)
