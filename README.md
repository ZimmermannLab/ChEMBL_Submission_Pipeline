[![License](https://img.shields.io/badge/License-MIT%202.0-blue.svg)](https://opensource.org/licenses/MIt)
# ChEMBL Submission Pipeline
This README markdown refers to the pipeline for converting supplementary tables in CSV format to ChEMBL-Submission-Ready TSV tables, specifically for Biotransformation Data. The pipeline requires both R and Python installations, although majority of it is written in R, the chemical compounds data handling is executed using RDKit in Python.

## Inputs
1. compounds.csv
2. assay_input.csv
3. activity_inputs.csv <br>
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
This R script calls three functions, ```write_activity_tsv.R```, ```write_assay_tsv.R```, ```write_reference_tsv.R```. <br>
However, the script also contains the ##commented out script to convert Excel supplementary sheets from our use-case dataset to convert these Excel sheets to standardised CSV format specific gut bacteria-mediated biotransformation. 
Based on the output files, the whole script is divided into 4 sections:
#### REFERENCE.tsv ####
```write_reference_tsv``` takes various inputs to generate a reference TSV file according to ChEMBL rules. 

```r
RIDX <- c("HumanMicrobiome_DrugMetabolism")
DOI <- c("10.1038/s41586-019-1291-3")
REF_TYPE <- c("Publication")
TITLE <- c("Mapping human microbiome drug metabolism by gut bacteria and their genes") #mandatory
AUTHORS	<- c("Michael Zimmermann, Maria Zimmermann-Kogadeeva, Rebekka Wegmann & Andrew L. Goodman") #mandatory
ABSTRACT <- c("Individuals vary widely in their responses to medicinal drugs, which can be dangerous and expensive owing to treatment delays and adverse effects.")
output_dir <- "/ChEMBL_Submission_Pipeline/outputs
ref_tbl <- write_reference_tsv(output_dir, RIDX, DOI, TITLE, AUTHORS, ABSTRACT, REF_TYPE)
```

#### COMPOUND_RECORD.tsv and COMPOUND_CTAB.sdf ####
```generate_compound_files``` creates both TSV and SDF files and is based on RDKit. 

```python
output_dir = "/ChEMBL_Submission_Pipeline/ChEMBL_Submission_Pipeline/outputs"
input_file = "/ChEMBL_Submission_Pipeline/ChEMBL_Submission_Pipeline/inputs/compounds.csv"
RIDX = "HumanMicrobiome_DrugMetabolism"
prefix = "HMDM"

generate_compound_files(output_dir, input_file, RIDX, prefix)
```

#### ASSAY.tsv ####
```write_assay_tsv``` takes a supplementary CSV file as input to generate ASSAY.tsv.

```r
input_csv <- "/ChEMBL_Submission_Pipeline/inputs/assay_input.csv"
source <- "Zimmermann"

assay_tbl <- write_assay_tsv(output_dir, input_csv, ridx, category = c("bacteria", "enzyme", "community"), source)
```


#### ACTIVITY.tsv ####
```write_activity_tsv``` takes a supplementary CSV file as input to generate ACTIVITY.tsv.

```r
input_csv <- "ChEMBL_Submission_Pipeline/inputs/activity_input.csv"
compound_tsv <- "/ChEMBL_Submission_Pipeline/outputs/COMPOUND_RECORD.tsv" 
assay_tsv <- "/ChEMBL_Submission_Pipeline/outputs/ASSAY.tsv"


activity_tbl <- write_activity_tsv(output_dir, input_csv, compound_tsv, assay_tsv,
                   ridx, type = "Biotransformation")
```
