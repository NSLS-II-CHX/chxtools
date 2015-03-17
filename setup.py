#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function)

import sys
import warnings


from setuptools import setup

setup(
    name='chxtools',
    version='0.0.0',
    author='Brookhaven National Laboratory',
    packages=['chxtools'],
    entry_points={
        'console_scripts': [
            'replay = dataportal.replay.replay:main']},
    package_data={'chxtools': 'data/*'},
    include_package_data=True,
)
