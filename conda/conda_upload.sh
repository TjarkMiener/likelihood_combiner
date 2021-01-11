# Only need to change these two variables
PKG_NAME=lklcom
USER=tmiener

OS=$TRAVIS_OS_NAME-64
mkdir ~/conda-bld
conda config --set anaconda_upload yes
export CONDA_BLD_PATH=~/conda-bld
LKLCOM_VERSION=`python -c "import likelihood_combiner; print(likelihood_combiner.__version__);"

export VERSION=$LKLCOM_VERSION
conda build .
anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER -l nightly $CONDA_BLD_PATH/$OS/$PKG_NAME-$LKLCOM_VERSION.tar.bz2 --force
