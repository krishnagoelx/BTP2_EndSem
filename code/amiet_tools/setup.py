# -*- coding: utf-8 -*-
"""
Amiet_tools setup file
=================
"""

from setuptools import setup

settings = {
    'name': 'amiet_tools',
    'version': '0.0.2',
    'description': "Collection of Python functions for prediction of turbulence-flat plate interaction noise using Amiet's model",
    'license': 'BSD 3-Clause',
    'install_requires': [
        'numpy>=1.20.3',
        'scipy>=1.6.3'],
    'packages': ['amiet_tools', 'amiet_tools.classes', 'amiet_tools.functions']
    }

setup(**settings)