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
      license='GPL-3.0',
      packages=['likelihood_combiner'],
      install_requires=[
          'matplotlib',
          'numpy>=1.15.0',
          'scipy',
          'jupyter',
          'pandas',
          'pyyaml',
          'tables',
          ],
      entry_points = {
        'console_scripts': ['lklcom-local=likelihood_combiner.local:main',
                            'lklcom-cluster=likelihood_combiner.cluster:main',
                            'gLike_to_lklcom=likelihood_combiner.io:_gLike_to_lklcom',
                            'lklcom_to_gLike=likelihood_combiner.io:_lklcom_to_gLike',
                            'gLikeLimits_to_lklcomLimits=likelihood_combiner.io:_gLikeLimits_to_lklcomLimits',
                            'merge_to_lklcom=likelihood_combiner.io:_merge_to_lklcom']
      },
      include_package_data=True,
      dependencies=[],
      dependency_links=[],
      zip_safe=False)
