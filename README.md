# GATA-DNA binding
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

## Analysis
Run the `abcd.py` script to first create AF3 input files (`create` mode) and then to analyze the AF3 predictions (`dkadj` mode). 

## Results

The `./gata_output/` directory contains the following:

* af3_inputs
	* dimer.json: AF3 input JSON file. Upload this file on the AF3 webserver to obtain the predictions.

* analysis
	* rb_dir_path.json: JSON file containing the paths to the rigid_body .pdb file obtained from AF-Pipeline.
	* rb_figures: ChimeraX figures for the riggid bodies obtained from AF-Pipeline.
	* rigid_body: contains the rigid bodies extracted from the AF3 predictions using AF-Pipeline.
