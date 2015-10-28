#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function)

from setuptools import setup

setup(
    name='chxtools',
    version='0.0.0',
    author='Brookhaven National Laboratory',
    packages=['chxtools', 'pyXPCS', 'chxtools.pims_readers'],
    package_data={'chxtools/X-ray_database' : ['.dat']},
    include_package_data=True,
)
