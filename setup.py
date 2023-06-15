from setuptools import setup, find_packages

# with open('requirements.txt') as f:
#     required = f.read().splitlines()
# print(required)
from csdid._version import __version
print(__version)

setup(
  name = 'csdid',
  version=__version,
  url='https://github.com/d2cml-ai/csdid',
  author='D2CML Team, Alexander Quispe, Carlos Guevara, Jhon Flroes',
  keywords=['Causal inference', 'Research'],
  license="MIT",
  description='Difference in Difference in Python',
  classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Topic :: Scientific/Engineering",
    ],
  install_requires=[
        'pandas',
        'numpy<=1.24.3',
        'scipy',
        'patsy',
        'plotnine',
        'twine'
  ],
  packages=find_packages(),
  package_data={
    'data': ['data/*'],
    'configs': ['configs/*']
  }
)