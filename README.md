# ChEMBL Submission Pipeline
This README markdown refers to the pipeline for converting supplementary tables in CSV format to ChEMBL-Submission-Ready TSV tables, specifically for Biotransformation Data. The pipeline requires both R and Python installations, although majority of it is written in R, the chemical compounds data handling is executed using RDKit in Python.

## Inputs
1. compounds.csv
2. assay_input.csv
3. activity_inputs.csv
These are the compulsory input CSV files.

## Outputs
1. REFERENCE.tsv
2. COMPOUND_RECORD.tsv
3. COMPOUND_CTAB.sdf
4. ASSAY.tsv
5. ACTIVITY.tsv
Other important TSV files such as INFO.txt and ASSAY_PARAM.tsv are highly user/experiment dependent.

## Directory Organization

### Script
The Script folder consists of two files:
#### 1. main_PYnotebook.ipynb (Jupyter notebook with Python kernel)
The notebook contains two options: one to use SDF as an input and one to use a csv as input to generate both COMPOUND_RECORD.tsv and COMPOUND_CTAB.sdf. This notebook calls for function ```python/generate_compound_files```. 
#### 2. main_Rscript.R (R script)
This R script calls three functions, ```write_activity_tsv.R```, ```write_assay_tsv.R```, ```write_reference_tsv.R```.
