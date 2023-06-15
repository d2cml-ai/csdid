from setuptools import setup, find_packages

with open('requirements.txt') as f:
  requi = f.read().splitlines()

from csdid._version import __version



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
  install_requires=requi,
  packages=find_packages(),
  package_data={
    'data': ['data/*'],
    'configs': ['configs/*']
  }
)