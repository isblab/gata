# gata
Analyzing the effect of GATA protein mutants on its ability to bind DNA.

## Setup
### Create a conda environment
```
conda env create -f environment.yml -n gata
conda activate gata
```

### Setup AF-pipeline
```
git clone --recursive https://github.com/isblab/af_pipeline.git
```
Add `$PATH_TO_AF_PIPELINE` to `~/.bash_profile`. Then run,
```
source ~/.bash_profile
```
