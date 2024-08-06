# Phenotype Propagator
Phenotype Propagator uses gene-level mcSEED functional annotations for in silico predictions of the presence or absence (denoted as binary: “1” or “0”) of 158 functional metabolic pathways using a semi-automated process based on a combination of the following three approaches:
In the Version 2 we use only Machine Learning (ML)-based phenotype predictions.

## Prerequisites

- R
  - data.table
  - openxlsx
  - caret
  - yaml

## Running Phenotype Propagator
To run the pipeline first prepare the input files by placing them into `input/annotation` (danatello outputs with mcSEED functional annotations for each target genome/MAG) and `input/MAG_taxonomy.txt` with species and genus names for each genome/MAG.
Then update `PhenotypePropagator.yaml` config file (including `root_dir` and other directories)

Finally, run `Rscript PhenotypePropagator.R` 

## Output files
`BPM.elsx` - The main output file containing phenotype predictions for each target genomes/MAGs and each of 158 phenotypes. The output is presented in a form of Binary Phenotype Matrix (BPM). 
`MAGpredictionsFull.xlsx` - Detailed output file containing phenotype predictions and the distribution of metabolic pathway genes in each target genome and each mcSEED subsystem.