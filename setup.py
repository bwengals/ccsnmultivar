from setuptools import setup, find_packages


setup(
  name = 'ccsnmultivar',
  packages = ['ccsnmultivar'], # this must be the same as the name above
  version = '0.1.23',
  description = 'Multivariate regression analysis of core-collapse simulations',
  author = 'Bill Engels',
  author_email = 'w.j.engels@gmail.com',
  url = 'https://github.com/bwengals/CCSNMultivar', # use the URL to the github repo
  download_url = 'https://github.com/bwengals/CCSNMultivar/tarball/0.1', # I'll explain this in a second
  keywords = ['regression', 'core-collapse', 'supernova'], # arbitrary keywords
  classifiers = [],
  install_requires=[
       'numpy',
       'scipy', 
       'tabulate',
        'patsy',
        'pandas',
        'scikit-learn',
    ]
)
