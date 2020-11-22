# LikelihoodCombiner: Combining likelihoods from different experiments.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3666282.svg)](https://doi.org/10.5281/zenodo.3666282)

![Gloryduck logo](images/Gloryduck_logo.png)

**LikelihoodCombiner** is a package under active development to combine likelihoods from different experiments. The main target of this package is the **Gloryduck** project. This project joint analysis of gamma-ray data from [*Fermi*-LAT](https://glast.sites.stanford.edu/), [HAWC](https://www.hawc-observatory.org/), [H.E.S.S.](https://www.mpi-hd.mpg.de/hfm/HESS/), [MAGIC](https://magic.mpp.mpg.de/) and [VERITAS](https://veritas.sao.arizona.edu/) to search for gamma-ray signals from dark matter annihilation in dwarf satellite galaxies.

## Install LikelihoodCombiner

### Clone Repository with Git

Clone the LikelihoodCombiner repository:

```bash
cd </installation/path>
git clone https://github.com/TjarkMiener/likelihood_combiner
```

### Install Package with Anaconda

Next, download and install [Anaconda](https://www.anaconda.com/download/), or, for a minimal installation, [Miniconda](https://conda.io/miniconda.html). Create a new conda environment that includes all the dependencies for LikelihoodCombiner:

```bash
conda env create -f </installation/path>/likelihood_combiner/environment.yml
```

Finally, install LikelihoodCombiner into the new conda environment with pip:

```bash
conda activate lklcom
cd </installation/path>/likelihood_combiner
pip install --upgrade .
```

NOTE for developers: If you wish to fork/clone the respository and make changes to any of the LikelihoodCombiner modules, the package must be reinstalled for the changes to take effect.

### Installing as a conda package

To install it as a conda package, first install Anaconda by following the instructions here: https://www.anaconda.com/distribution/.

Then, create and enter a new Python 3.8 environment with:

```bash
conda create -n [ENVIRONMENT_NAME] python=3.8
source activate [ENVIRONMENT_NAME]
```

From the environment, add the necessary channels for all dependencies:

```bash
conda config --add channels conda-forge
conda config --add channels menpo
```

Install the package:

```bash
conda install -c tmiener likelihood_combiner
```

This should automatically install all dependencies (NOTE: this may take some time, as by default MKL is included as a dependency of NumPy and it is very large).

If you want to import any functionality from LikelihoodCombiner into your own Python scripts, then you are all set. However, if you wish to make use of any of the scripts in likelihood_combiner/scripts (like {local/cluster}.py), you should also clone the repository locally and checkout the corresponding tag (i.e. for version v0.4.1):

```bash
git clone https://github.com/TjarkMiener/likelihood_combiner
git checkout v0.4.1
```

LikelihoodCombiner should already have been installed in your environment by Conda, so no further installation steps (i.e. with setuptools or pip) are necessary and you should be able to run scripts/{local/cluster}.py directly.

### Dependencies

- Python 3.8.X
- NumPy
- SciPy
- Pandas
- PyTables
- PyYAML
- Matplotlib
  
## Run the Combiner

Run LikelihoodCombiner from the command line:
  
```bash
LikelihoodCombiner_dir=</installation/path>/likelihood_combiner
python $LikelihoodCombiner_dir/scripts/{local|cluster}.py $LikelihoodCombiner_dir/config/example_config.yml 
```
  
### Mock data 

The data you can find in the LikelihoodCombiner, where produced with [gLike](https://github.com/javierrico/gLike/) using the [mock data](https://github.com/javierrico/gLike/tree/master/data). These txt files **don't** correspond to IACT observations of Segue 1 or Ursa Major II and are only included for testing the code framework.

## Uninstall LikelihoodCombiner

### Remove Anaconda Environment

First, remove the conda environment in which LikelihoodCombiner is installed and all its dependencies:

```bash
conda remove --name lklcom --all
```

### Remove LikelihoodCombiner

Next, completely remove LikelihoodCombiner from your system:

```bash
rm -rf </installation/path>/likelihood_combiner
```
