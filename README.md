[Add PubMed link]: [![PubMed]

# GATA3-extension Mutants Rewire Lineage Identity in Luminal Breast cancer

This repository contains the scripts for the analysis of the effect of GATA protein mutants on its ability to bind DNA.

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
**Step 1**: To generate the AF3 input JSON files, run the following command
```
python gata_analysis.py -i
```

**Step 2**: Use the input JSON file to obtain predictions from the [AF3 webserver](https://alphafoldserver.com/).


**Step 3**: For the performing the analysis, run the following command,
```
python gata_analysis.py -a
```

## Results

The `./gata_output/` directory contains the following:

* af3_inputs
	* dimer.json: AF3 input JSON file. Upload this file on the AF3 webserver to obtain the predictions.

* af3_preds: contains the best predicted structure (based on the ranking score) obtained from AF3.

* analysis
	* rb_dir_path.json: JSON file containing the paths to the rigid_body .pdb file obtained from AF-Pipeline.
	* rb_figures: ChimeraX figures for the riggid bodies obtained from AF-Pipeline.
	* rigid_body: contains the rigid bodies extracted from the AF3 predictions using AF-Pipeline.


## Information
__Author(s):__ Kartik Majila

__Date__: Aug 19, 2026

**License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License.

__Testable:__ Yes

__Publication__ Nivedhya Venas, Mrunal Ratnaparkhi, Kartik Majila, Shruthi Viswanath, Dimple Notani, Radhakrishnan Sabarinathan, __GATA3-extension Mutants Rewire Lineage Identity in Luminal Breast cancer__, submitted.
