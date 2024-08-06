# danatello
Danatello is a DIAMOND-based annotation pipeline. One can use it to propagate mcSEED annotations to isolate genomes and MAGs. The output is  used for phenotype propagation pipeline.

## Prerequisites

### Software
* `python=3.9`
* `diamond=2.1.4`
* `mmseqs2=14.7e284`
* `samtools=1.16`
* `pandas`
* `numpy`
* `tqdm`
* `scikit-learn`
* `scipy`

To create a conda environment with these tools, use:
```
mamba env install -f envs/danatello.yml
```

### Subsystem table
You will need two mcSEED database files that can be downloaded from https://doi.org/10.5281/zenodo.10041396
1. 'subsystems.txt` -  the subsystems table that include functional annotations and categories for all genes from the mcSEED subsystems.
2. 'all_no_func_representative_d.dmnd' - DIAMOND database Contains 2856 reference mcSEED proteomes clustrered to avoid redundancy and enriched by all functionally annotated proteins from mcSEED subsystems (version compiled on January 25, 2021).  
Put these files to `database/` folder

### Genome/MAG FAA files for annotation
1. Annotate FNA files (MAGs or genomes) using `Prokka`
2. Put resulting FAA files to `faa/`

## Directory structure
 - `code/` - scripts used in the pipeline
 - `database` - mcSEED database files used for annotation (i) `subsystems.txt`, (ii) 'all_no_func_representative_d.dmnd'
 - `faa/` - FAA files (isolate genome or MAG) one wants to annotate
 - `DIAMOND/` - DIAMOND output files (created by `run_donatello.sh`)
 - `annotation/` - danatello output files (created by `run_donatello.sh`)
 - `run_donatello.sh` - script that runs danatello

## Running danatello
Run the danatello pipeline. To do that, use the `run_danatello.sh` script.

### Steps
1. Edit `run_danatello.sh`:
  * Set the `db_dir` variable (path to the directory with the DIAMOND database)
  * Set the `ncpu` variable (number of  CPUs to use)
2. Run the script
```
bash run_danatello.sh
```

## Output
Output files are in the `annotation/` folder.
For each input FAA file, the pipeline will create a table with multiple columns. The most important are `ID` and `Winner` ones:
1. ID of the protein sequence
2. Functional role (best hit)

Note that each functional role has its own line, and there could be multiple lines for one sequence. Also, some roles could be repeated.