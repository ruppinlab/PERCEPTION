
if [ -z "$1" ]
then
    PYTHON_VERSION=3.6
else
    PYTHON_VERSION=$1
fi

git clone https://gitlab.com/bu_cnio/conda-recipes.git
conda mambabuild conda-recipes/r-beyondcell --output-folder ./
mamba install --use-local --update-deps linux-64/r-beyondcell-*.tar.bz2
