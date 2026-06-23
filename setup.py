from setuptools import setup, find_packages
import os

# with open('requirements.txt') as f:
#     required = f.read().splitlines()
# print(required)
from csdid._version import __version
print(__version)

# Long description for PyPI (rendered from the project README).
_here = os.path.abspath(os.path.dirname(__file__))
try:
    with open(os.path.join(_here, 'readme.md'), encoding='utf-8') as _f:
        long_description = _f.read()
except Exception:
    long_description = 'Difference in Difference in Python'

setup(
  name = 'csdid',
  version=__version,
  url='https://github.com/d2cml-ai/csdid',
  author='D2CML Team, Alexander Quispe, Carlos Guevara, Jhon Flroes',
  keywords=['Causal inference', 'Research'],
  license="MIT",
  description='Difference in Difference in Python',
  long_description=long_description,
  long_description_content_type='text/markdown',
  classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Topic :: Scientific/Engineering",
    ],
  install_requires=[
        'pandas',
        'numpy',
        'scipy',
        'patsy',
        'plotnine',
        'twine',
        'joblib',
        'statsmodels',
        'drdid'
  ],
  packages=find_packages(),
  package_data={
    'data': ['data/*'],
    'configs': ['configs/*']
  }
)
