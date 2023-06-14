from setuptools import setup, find_packages

setup(
  name = 'csdid',
  version='0.1.66',
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
  packages=find_packages(),
  package_data={
    'data': ['data/*'],
    'configs': ['configs/*']
  }
  
)