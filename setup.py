from setuptools import setup, find_packages

setup(name='lklcom',
      version='0.4.1',
      description='LikelihoodCombiner combines likelihoods from different experiments.',
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
