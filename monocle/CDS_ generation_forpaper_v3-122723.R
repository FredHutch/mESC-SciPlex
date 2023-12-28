#Add packages necessary for Sciplex analysis

suppressPackageStartupMessages({
  library(rlang)
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(Matrix)
  library(tidyverse)
  library(monocle3)
  library(forcats)
  library(pheatmap)
  library(RColorBrewer)
  library(reshape2)
  library(tibble)
  library(Rcpp)
  library(reticulate)
  library(VGAM)
  library(viridis)
  library(magrittr)
  library(mygene)
  library(ggpubr)
  library(scales)
  library(org.Mm.eg.db)
  library(WriteXLS)
  #library(garnett) #Do not have loaded
  
  DelayedArray:::set_verbose_block_processing(TRUE)
  options(DelayedArray.block.size=1000e7)})


writeLines(capture.output(sessionInfo()), "sessionInfo-mESC-SciPlex_CDS-generation.txt")

#Set project directory
projectdir <- "/Users/adamheck/Desktop/mESC-SciPlex"
inputdir <- paste(projectdir, "processed_data/AH_CDS_122023", sep = "/")
plotdir <- paste(projectdir, "results/PLOTS", sep = "/")
outputdir <- paste(projectdir, "results", sep = "/")
setwd(projectdir)
# Load hashtable
hashTable1 = 
  read.table(file = "/Users/adamheck/Desktop/mESC-SciPlex/processed_data/Sanjay_files/hashTable.out",
             sep = "\t",
             header = F,
             col.names = c("sample", "Cell", "Oligo", "axis", "Count"))

# hash table pre processing
hash_dataframe1 = 
  hashTable1 %>%
  
  dplyr::select(Cell, 
                Oligo, 
                Count)%>%
  spread(key = Oligo, 
         value = Count,
         fill = 0)
hash_dataframe_cell1 = hash_dataframe1$Cell
hash_dataframe1 = hash_dataframe1[,2:ncol(hash_dataframe1)]
max_oligo = 
  colnames(hash_dataframe1)[apply(X = hash_dataframe1, 
                                 MARGIN = 1, 
                                 FUN = (which.max))]
max_slide_val = 
  apply(X = hash_dataframe1, 
        MARGIN = 1, 
        FUN = function(x){maximum_col = which.max(x)
          x[maximum_col]})
ratio_top_2 = 
  apply(X = hash_dataframe1, 
        MARGIN = 1, 
        FUN = function(x){
          maximum_col = which.max(x)
          max_val = x[maximum_col]
          x[maximum_col] = 0
          seoncd_max = which.max(x)
          second_max_val = x[seoncd_max]
          max_val/second_max_val
        })
hash_df = data.frame(
  Cell = hash_dataframe_cell1,
  max_val = max_slide_val,
  max_id = max_oligo,
  ratio_top_2 = ratio_top_2,
  hash_total)

# Import complete cds
SP_cds = 
  readRDS("/Users/adamheck/Desktop/mESC-SciPlex/processed_data/Sanjay_files/cds.RDS")
coldata_SP_cds = colData(SP_cds) %>% as_tibble()

# Filter cds based on hash table
coldata =
  colData(SP_cds) %>%
  as.data.frame() %>%
  left_join(hash_df, 
            by = "Cell")
rownames(coldata) = 
  coldata$Cell
rowdata =
  rowData(SP_cds)
rownames(rowdata) =
  rownames(SP_cds)
SP_cds = 
  new_cell_data_set(expression_data = counts(SP_cds),
                    gene_metadata = rowdata,
                    cell_metadata = coldata)
coldata_SP_cds = colData(SP_cds) %>% as_tibble()
cells_to_keep = 
  colData(SP_cds) %>%
  as.data.frame() %>%
  filter(!is.na(ratio_top_2),
         ratio_top_2 >= 5) %>%
  pull(Cell) %>%
  as.character()
SP5_cds = SP_cds[,cells_to_keep]
#view coldata
coldata_SP5_cds = colData(SP5_cds) %>% as_tibble()

# Add metadata to cds
colData(SP5_cds)$data_set = 
  ifelse(grepl(x = colData(SP5_cds)$max_id,
               pattern = "day"),
              "EB differentiation",
              "Embryo")
# Create EB conditions key for metadata
condition_key_EB = 
  data.frame(
    max_id = colData(SP5_cds) %>%
      as.data.frame() %>%
      filter(grepl(pattern = "day",
                   x = max_id)) %>%
      pull(max_id) %>%
      unique())
condition_key_EB$day = 
  stringr::str_split_fixed(condition_key_EB$max_id,
                           "_",
                           n = 2)[,1] %>%
  stringr::str_replace_all(pattern = "day",
                           replacement = "")
condition_key_EB$condition = 
  stringr::str_split_fixed(condition_key_EB$max_id,
                           "_",
                           n = 2)[,2] 
condition_key_EB = 
  condition_key_EB %>%
  left_join(data.frame(condition = 
                         seq(1,16,1) %>% 
                         as.character(),
                       BMP = 
                         c(0.3,0.3,0.3,0.3,
                           1.25,1.25,1.25,1.25,
                           5,5,5,5,
                           20,20,20,20),
                       Activin = 
                         c(0,1.25,5,20,
                           0,1.25,5,20,
                           0,1.25,5,20,
                           0,1.25,5,20)),
            by = "condition")
# add EB condition keys to cds
coldata =
  colData(SP5_cds) %>%
  as.data.frame() %>%
  left_join(condition_key_EB) %>%
  data.frame(row.names = colData(SP5_cds)$Cell)

SP5_cds = 
  new_cell_data_set(expression_data = counts(SP5_cds),
                    gene_metadata = rowdata,
                    cell_metadata = coldata)

# Create embryo condition key for metadata
condition_key_Embryo = 
  data.frame(
    max_id = colData(SP5_cds) %>%
      as.data.frame() %>%
      filter(grepl(pattern = "sm",
                   x = max_id)) %>%
      
      pull(max_id) %>%
      unique())
condition_key_Embryo$sm = 
  stringr::str_split_fixed(condition_key_Embryo$max_id,
                           "_",
                           n = 2)[,1] %>%
  stringr::str_replace_all(pattern = "sm",
                           replacement = "")
condition_key_Embryo$part = 
  stringr::str_split_fixed(condition_key_Embryo$max_id,
                           "_",
                           n = 2)[,2]  
condition_key_Embryo$part =   
  stringr::str_split_fixed(condition_key_Embryo$part,
                           "_",
                           n = 2)[,1]
#add Embryo condition keys to cds
coldata =
  colData(SP5_cds) %>%
  as.data.frame() %>%
  left_join(condition_key_Embryo) %>%
  data.frame(row.names = colData(SP5_cds)$Cell)

SP5_cds = 
  new_cell_data_set(expression_data = counts(SP5_cds),
                    gene_metadata = rowdata,
                    cell_metadata = coldata)

# Define E7 cells
condition_key_Embryo_E7 = 
  data.frame(
    max_id = colData(SP5_cds) %>%
      as.data.frame() %>%
      filter(grepl(pattern = "E7",
                   x = max_id)) %>%
      
      pull(max_id) %>%
      unique())
condition_key_Embryo_E7$stage = 
  stringr::str_split_fixed(condition_key_Embryo_E7$max_id,
                           "_",
                           n = 2)[,1] %>%
  stringr::str_replace_all(pattern = "stage",
                           replacement = "")
coldata =
  colData(SP5_cds) %>%
  as.data.frame() %>%
  left_join(condition_key_Embryo_E7) %>%
  data.frame(row.names = colData(SP5_cds)$Cell)

SP5_cds = 
  new_cell_data_set(expression_data = counts(SP5_cds),
                    gene_metadata = rowdata,
                    cell_metadata = coldata)
coldata_SP5_cds = colData(SP5_cds) %>% as_tibble()

####Subset cells according to data_set eg. Embryo vs EB differentiation####
# Filter for EB differentiation
EB3_cells_to_keep =
  colData(SP5_cds) %>%
  as.data.frame() %>%
  filter(!is.na(data_set),
         data_set == "EB differentiation") %>%
  pull(Cell) %>%
  as.character()

EB3_cds = SP5_cds[,EB3_cells_to_keep]
coldata_EB3_cds = colData(EB3_cds) %>% as_tibble()

# Filter for E7.5-E8, embryo proper and yolk sac
EYS_cells_to_keep =
  colData(Emb1_cds) %>%
  as.data.frame() %>%
  filter(!is.na(part), part %in% c("EP", "YS")) %>%
  pull(Cell) %>%
  as.character()
EYS_cds = Emb1_cds[, EYS_cells_to_keep]
coldata_EYS_cds = colData(EYS_cds) %>% as_tibble()

####PROCESS AND CLUSTER CDS####
#### Process EB cds with set seed ####
EB3_cds <- preprocess_cds(EB3_cds, num_dim = 98)
EB3_cds <- reduce_dimension(EB3_cds, umap.n_neighbors = 4L)
# Cluster cells
set.seed(17)
EB3_cds = cluster_cells(EB3_cds, resolution = 3.6e-4, random_seed = 17)
# Add umap data to cds
colData(EB3_cds)$umap1 = reducedDim(EB3_cds, type = "UMAP")[,1]
colData(EB3_cds)$umap2 = reducedDim(EB3_cds, type = "UMAP")[,2]
colData(EB3_cds)$Cluster = clusters(EB3_cds, reduction_method = "UMAP")

# Process embryo/yolk sac dataset
EYS_cds <- preprocess_cds(EYS_cds, num_dim =40)
EYS_cds <- reduce_dimension(EYS_cds, umap.n_neighbors = 4L)
# Cluster cells
set.seed (17)
EYS_cds <- cluster_cells(EYS_cds, resolution = 2e-3, random_seed = 17)
# Add umap data to cds
colData(EYS_cds)$umap1 = reducedDim(EYS_cds, type = "UMAP")[,1]
colData(EYS_cds)$umap2 = reducedDim(EYS_cds, type = "UMAP")[,2]
colData(EYS_cds)$Cluster = clusters(EYS_cds, reduction_method = "UMAP")

####Load and or save cds data frames####
#Save SP5_cds:
saveRDS(SP5_cds, file.path("/Users/adamheck/Desktop/mESC-SciPlex/processed_data/AH_CDS_122023/", "SP5.RDS"))
#Save EB3_cds:
saveRDS(EB3_cds, file.path("/Users/adamheck/Desktop/mESC-SciPlex/processed_data/AH_CDS_122023/", "EB3.RDS"))
#Save EYS_cds
saveRDS(EYS_cds, file.path("/Users/adamheck/Desktop/mESC-SciPlex/processed_data/AH_CDS_122023/", "EYS.RDS"))

# Turn off plotter
dev.off()
