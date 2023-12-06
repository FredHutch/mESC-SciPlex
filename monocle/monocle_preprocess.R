################################################################
#__Author:__ Adam Heck
#__Script:__ monocle_preprocess.R
#__Project:__ "INSERT PROJECT/REPO NAME HERE"
#__Summary:__ "IF WANTED LIST SUMMARY/DESCRIPTION OF ANALYSIS HERE"
#################################################################
rm(list = ls()) ##clears environment

#Load Packages
library(monocle3)
library(leidenbase)
library(VGAM)  ## required for negbinomial.size ???
library(viridis)###load package for viridis color themes
library(stringr)
library(tibble)
library(dplyr)##library for regression analysis/violin plot
library(magrittr)##regression analysis/violin
library(ggplot2)

#Set up project directory
projectdir <- "~/Desktop/PATH_TO_PROJECT_DIRECTORY" #This is the location of repo on local machine, typically store on Desktop, can change at user discretion
fastqdir <- "~/Desktop/PATH_TO_FASTQ_FILES" #Location of fastq files. This SHOULD NOT be in the same location as project/repo directory
resdir <- paste(projectdir, "results", sep = "/") #Subdirectory in project repo where results are stored
pddir <- paste(projectdir, "processed_data", sep = "/") #Subdirectory in project repo where processed data files are stored, I.e. monocle3 rds files

# for PCs the file paths have a different format. instead of "~/" they require "C:/"

# Import 10x data using load_cellranger_data()
setwd(fastqdir)
cds_1 <- load_cellranger_data("FILE_FOLDER_1")
cds_2 <- load_cellranger_data("FILE_FOLDER_2")
# file/folder names should be changed at user discretion for specific datasets

# If necesszry, use combine_cds() to merge 
cds <- combine_cds(list(cds_1, cds_2))

colnames(cds)[2] <- "gene_short_name"  ##helps with the violin plot function

# Preprocess the data, set correct number of dimensions and reduce the data
cds <- preprocess_cds(cds, num_dim = 10) #Number of dimensions used will vary for each dataset, use plot_pc_variance and plot_cells after UMAP reduction to decided how many dimensions
# plot_pc_variance_explained(cds)
plot_pc_variance_explained(cds) #general rule of thumb is the last (highest numbered) PCA component should explain >0.1 of the variance. See Monocle3 tutorial for more details.

# Remove batch effects, from samples,
cds <- align_cds(cds, alignment_group = "sample")

# Reduce dimensions, NOTE: to ensure repoduciblity, need to run with these options
cds <- reduce_dimension(cds,umap.fast_sgd = FALSE,cores=1,n_sgd_threads=1)

# Plot the cells
plot_cells(cds)
plot_cells(cds, color_cells_by = "sample") # Color cells by sample.

# Saving the full cds, unclustered
saveRDS(cds, file.path(pddir, "'DATE'_fullcds_UNCLUSTERED.rds")) #File name can be changed at user discretion, default is 'DATE_CDS-DETAILS.rds'

# Cluster cells, resolution automatically determine. Like doing this for now, less bias.
set.seed(17) ##to ensure reproducibility, need to have set.seed() function first, and random_seed = , in the options. Choose one number for both inputs.
cds = cluster_cells(cds, random_seed = 17, resolution = 5e-4) ##resolution will be different per sample, change it around to get a different number of clusters

# Plot cells by cluster to fine-tune resolution needed for analysis
plot_cells(cds, group_label_size = 6)

##always a good idea to click on your cds object in the "global environment" page and look into colData
##then listData, and see what types of variables you can label your cells by
##for knockout experiemnt you can label by treatment for example, or for embryos you can label by timepoint, E9 and E10 for example

#Once you have a cluster resolution you are satisfied with, add cluster to column to colData
cluster <- clusters(cds)
pData(cds)$cluster <- cluster   ##adds a "cluster" column to colData VERY HELPFUL, needed to get group_cells_by = "cluster" to work

# Saving the full cds, clustered
saveRDS(cds, file.path(pddir, "'DATE'_fullcds_CLUSTERED.rds")) #File name can be changed at user discretion, default is 'DATE_CDS-DETAILS.rds'

# This ends the initial processing, from this point on the analyses are project specific. If additional processing is needed, i.e. subsetting or labeling clusters, it can be added below
# For specific analyses, like find_gene_modules(), graph_test() or trajectory analyses, it is recommended to start a new script and load a preprocessed rds file
###############
###############

# Loading a preprocessed rds file
cds <- readRDS(file.path(pddir,"NAME_OF_RDS_FILE.rds"))

## An important note that if subsetting, some functions require the full cds be present in order to run properly. Therefore any scripts doing analysis on a subset of a cds should also load in the full cds.