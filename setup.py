from setuptools import setup, find_packages
import os

import ccsnmultivar

# translate README.md to something pypi understands
import os
long_description = ''
if os.path.exists('README.txt'):
    long_description = open('README.txt').read()

setup(
    name = 'ccsnmultivar',
    packages = ['ccsnmultivar'],
    version = ccsnmultivar.__version__,
    description = 'Multivariate regression analysis of core-collapse simulations',
    long_description = long_description,
    author = 'Bill Engels',
    author_email = 'w.j.engels@gmail.com',
    url = 'https://github.com/bwengals/CCSNMultivar',
    download_url = 'https://github.com/bwengals/CCSNMultivar/tarball/0.1',
    keywords = ['multivariate','regression', 'core-collapse', 'supernova'],
    classifiers = [],
    platforms = 'any',
    install_requires=[
       'numpy',
       'scipy',
       'tabulate',
       'scikit-learn',
    ]
)
