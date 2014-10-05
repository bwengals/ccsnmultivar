from setuptools import setup, find_packages
import io
import codecs
import os
import sys

import ccsnmultivar

here = os.path.abspath(os.path.dirname(__file__))

def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)
long_description = read('README.rst')


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
    keywords = ['regression', 'core-collapse', 'supernova'],
    classifiers = [],
    platforms = 'any',
    install_requires=[
       'numpy',
       'scipy',
       'tabulate',
       'patsy',
       'pandas',
       'scikit-learn',
    ]
)
