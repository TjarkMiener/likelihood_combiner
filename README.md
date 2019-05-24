# LikelihoodCombiner: Combining likelihoods from different experiments.

**LikelihoodCombiner** is a package under active development to combine likelihoods from different experiments. The main target of this package is the **Gloryduck** project. This project joint analysis of gamma-ray data from [MAGIC](https://magic.mpp.mpg.de/), [VERITAS](https://veritas.sao.arizona.edu/), [H.E.S.S.](https://www.mpi-hd.mpg.de/hfm/HESS/), [HAWC](https://www.hawc-observatory.org/) and [*Fermi*-LAT](https://glast.sites.stanford.edu/) to search for gamma-ray signals from dark matter annihilation in dwarf satellite galaxies.

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

### Dependencies

- Python 3.7.3
- NumPy
- PyTables
- PyYAML
- Matplotlib
  
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
