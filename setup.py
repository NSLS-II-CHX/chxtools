#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function)

import sys
import warnings


from setuptools import setup, find_packages

setup(
    name='chxtools',
    version='0.0.0',
    author='Brookhaven National Laboratory',
    packages=find_packages(),
    package_data={'chxtools/X-ray_database' : ['.dat']},
    include_package_data=True,
)
