from setuptools import setup, find_packages
from os import path
from likelihood_combiner.version import *

# Taken from https://github.com/me-manu/gammaALPs/blob/master/setup.py
here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='lklcom',
      version=get_version_pypi(),
      author="Tjark Miener",
      author_email="tmiener@ucm.es",
      description='LikelihoodCombiner combines DM-related likelihoods from different experiments.',
      long_description=long_description,
      long_description_content_type='text/x-rst',
      url='https://github.com/TjarkMiener/likelihood_combiner',
      license='BSD-3-Clause',
      packages=['likelihood_combiner'],
      install_requires=[
          'matplotlib',
          'numpy>=1.15.0',
          'scipy',
          'pandas',
          'pyyaml',
          'tables',
          ],
      include_package_data=True,
      dependencies=[],
      dependency_links=[],
      zip_safe=False)
