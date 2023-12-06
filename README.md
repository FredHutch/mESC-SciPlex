# Template repository for Hadland Lab Github team

### Once a new repository is created, please revise this README file based on project need.

### IMPORTANT: When creating a new repository using this template, one of the first things to do is create a "processed_data" folder within the project directory.

## Description

This is a template that can be used to create a repository for a new computational analysis project by the Hadland Lab. Below is an example of the data structure and an explanation of the components:

![Alt text](example_structure.png)

### results:
- contains raw scripts and data generated from each experiment/analysis in the project
- labarchives_path.txt lists where the project notebook can be found in the Hadland Lab notebook collection on Lab Archives
- experiment directories are formatted by date and can include a one word descriptor if desired
- each experiment directory should contain a brief summary.txt file, _runall script that contains the code and any data/results (i.e. graphs, tables, lists, etc.) that were created
- additional organizational logic and subdirectories within an experiment directory can be used at the discretion of the authors

### cellranger:
- contains scripts for [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome) processing of raw fastq data
- run_cellranger.sh is a template script for running the "cellranger count" function on the Fred Hutch servers using sbatch
- all names/options that begin with "CHANGE-" in the template script need to be modified for the specific project
  - **NOTE:** demultiplexing via a custom script or "cellranger mkfastq" needs to be done first, this is typically done by core staff when sequencing with the Hutch Genomics Core

### monocle:
- contains the compiled/finalized scripts for all [monocle3](https://cole-trapnell-lab.github.io/monocle3/) analyses performed in the project
- monocle_preprocess.R is a template script for processing cellranger output data into a useable cds file for downstream analysis
- additional scripts can be added for each experiment as they are finalized

### processed_data:
- contains the intermediate data files used during project analyses (i.e. rds files created by monocle preprocess script)
  - **NOTE:** this directory is already part of the .gitignore file and will not be tracked or uploaded to a remote repository due to potential large files stored here, it is present for ease of use when working on a local machine

## Additional Notes

This is a template repo, thus the codes, workflow and data structure are subject to change based on lab experience and individual user needs.

Repositories will begin as private and organized in this structure. Following manuscript submission/acceptance, they will be reorganized (potentially cloned and then reorganized) to suit the structure of the paper and then made publicly available.

The default permissions for repositories created from this template will give the creator and Brandon admin access and the Hadland Lab Github team write access.

Utilize branches when multiple people are working on a project simultaneously.

Local commits are recommended daily, or each time work is done on an experiment, to optimize tracking control. Pushes to a remote repository/branch can be done on a more intermittent basis at user discretion.

Questions, comments and suggestions about creating/using Github repos can be sent to Brandon, Adam or Rachel.
