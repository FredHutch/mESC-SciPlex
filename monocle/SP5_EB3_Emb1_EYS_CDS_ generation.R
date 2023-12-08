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
  setwd("/Users/adamheck/Desktop/mESC-SciPlex")
  
  DelayedArray:::set_verbose_block_processing(TRUE)
  options(DelayedArray.block.size=1000e7)})

#Add hashTable to Rsession and designate column names
#In this hashTable all cells are the same axis (1) and the same sample (sciChem_HSC_3)
#count is column indicating how many RNAs in that cell, oligo designates a particular group
hashTable1 = 
  read.table(file = "/Users/adamheck/Desktop/mESC-SciPlex/processed_data/Sanjay_files/hashTable.out",
             sep = "\t",
             header = F,
             col.names = c("sample", "Cell", "Oligo", "axis", "Count"))

# VERIFY CELL BELONGS TO A PARTICULAR GROUP AND A SUFFICIENT Number of RNAs are counted in that group
# Making separate dataframe with condition in columns - hashtag rows and count as variable... 
# A 0 is added if there is no count value 
hash_dataframe1 = 
  hashTable1 %>%
  
  dplyr::select(Cell, 
                Oligo, 
                Count)%>%
  spread(key = Oligo, 
         value = Count,
         fill = 0)
# Make value with just column Cell - allows adding back later to df
hash_dataframe_cell1 = hash_dataframe1$Cell

# Modifies df -removes replaces cell with number ?????
hash_dataframe1 = hash_dataframe1[,2:ncol(hash_dataframe1)]

# COMBINES HASHTABLE WITH CDS
# Get the Maximum slide oligo per Cell - Allows calculation of ratio_top_2.  Excludes cells with too many in multiple groups
#MARGIN = 1 means applies over rows
#which.max determines position of maximal value and then that value

#max_id column: max_oligo is condition - includes embryo data designated SP-somite pairs YS or EP yolk sac or embryo proper!
#Here an id is assigned based on which oligo has the max value - max_oligo is max_id
max_oligo = 
  colnames(hash_dataframe1)[apply(X = hash_dataframe1, 
                                 MARGIN = 1, 
                                 FUN = (which.max))]
#max_slide_val column - max_val in the end: Determines column (group) with highest count for cell 
max_slide_val = 
  apply(X = hash_dataframe1, 
        MARGIN = 1, 
        FUN = function(x){maximum_col = which.max(x)
          x[maximum_col]})
#ratio_top_2 column: Calculate ratio_top_2
#NOTE many cells are inf because there is no other cell with count max_val/0
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
#hash_total column: sums counts - not used in any further code but is useful to look at
hash_total = 
  rowSums(hash_dataframe1)

#Creates a df with columns from each above functions - 
#Cell = value created above - hash_dataframe_cell1
#max_val is the highest count found for the cell
#max_id is the oligo with the highest count
#ratio_top_2 insures cell has count and that it is not found in 2 groups but only in the group with max_value
hash_df = data.frame(
  Cell = hash_dataframe_cell1,
  max_val = max_slide_val,
  max_id = max_oligo,
  ratio_top_2 = ratio_top_2,
  hash_total)

#Add Sanjay's cds sent via dropbox to Rsession
SP_cds = 
  readRDS("/Users/adamheck/Desktop/mESC-SciPlex/processed_data/Sanjay_files/cds.RDS")
#view coldata - SP_cds contains 4 columns in coldata Cell (names cell), Sample (all sciChem_HSC_3), Size_Factor, n.umi
coldata_SP_cds = colData(SP_cds) %>% as_tibble()

#hash_df that lists conditions and embryo types as well as ratio_top_2 is added: 
#NOTE hash_df has 834,307 cells and Sanjay cds has 98,562 cells - Sanjay cds presumably has cells removed because of low n.umi!
#Most cells with inf in the ratio_top_2 (caused by division of 0) are gone (most likely weren't in Sanjay cds because n.umi was too low for these cells)
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

#view coldata - SP_cds contains 8 columns in coldata Cell (names cell), Sample (all sciChem_HSC_3), Size_Factor, n.umi, max_value, max_id, ratio_top_2, hash_total
coldata_SP_cds = colData(SP_cds) %>% as_tibble()

#filter requires Ratio_top_2 >=5 - creating new coldata with cells removed 
#SP5_cds is Ratio_top_2 >=5 - SP5_cds has 76,716 cells left

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

#### ID [MAX_ID] INFO SPREAD TO SEPARATE COLUMNS TO ALLOW SEPARATING BY FILTER ####
# Adds column to colData data_set: if has day in Max_id gets designated EB differentiation - if not is Embryo
colData(SP5_cds)$data_set = 
  ifelse(grepl(x = colData(SP5_cds)$max_id,
               pattern = "day"),
              "EB differentiation",
              "Embryo")
#Separates max_id EB info to Activin BMP concentrations with Activin and BMP being columns in colData along with condition number 1-16 #
# Create condition_key df that lists max_id components in separate columns
#Creates single column df with max_id only EB data because requires day to be in max_id (48 rows)
condition_key_EB = 
  data.frame(
    max_id = colData(SP5_cds) %>%
      as.data.frame() %>%
      filter(grepl(pattern = "day",
                   x = max_id)) %>%
      pull(max_id) %>%
      unique())
#Adds new column to condition_key df day and pulls number after day in max_id out
condition_key_EB$day = 
  stringr::str_split_fixed(condition_key_EB$max_id,
                           "_",
                           n = 2)[,1] %>%
  stringr::str_replace_all(pattern = "day",
                           replacement = "")
#Adds new column to condition_key df condition is 1-16
condition_key_EB$condition = 
  stringr::str_split_fixed(condition_key_EB$max_id,
                           "_",
                           n = 2)[,2] 
#Adds 2 new columns Activin and BMP based on condition # 1-16
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
#add EB condition keys to cds
coldata =
  colData(SP5_cds) %>%
  as.data.frame() %>%
  left_join(condition_key_EB) %>%
  data.frame(row.names = colData(SP5_cds)$Cell)

SP5_cds = 
  new_cell_data_set(expression_data = counts(SP5_cds),
                    gene_metadata = rowdata,
                    cell_metadata = coldata)

#Add columns to cds to allow filtering of different Embryo subsets: there are 5 - 7-8SP YS vs EP 12-13SP YS vs EP and E7 embryo
#Creates single column df with max_id only Embryo data because requires sm to be in max_id (19 rows)
condition_key_Embryo = 
  data.frame(
    max_id = colData(SP5_cds) %>%
      as.data.frame() %>%
      filter(grepl(pattern = "sm",
                   x = max_id)) %>%
      
      pull(max_id) %>%
      unique())
#Adds new column to condition_key df day and pulls number after day in max_id out
condition_key_Embryo$sm = 
  stringr::str_split_fixed(condition_key_Embryo$max_id,
                           "_",
                           n = 2)[,1] %>%
  stringr::str_replace_all(pattern = "sm",
                           replacement = "")
#Adds new column to condition_key df indicating whether embryo cell derives from yolk sac or embryo proper
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

#Creates single column df with max_id only E7 data because requires E7 to be in max_id (5 rows)
condition_key_Embryo_E7 = 
  data.frame(
    max_id = colData(SP5_cds) %>%
      as.data.frame() %>%
      filter(grepl(pattern = "E7",
                   x = max_id)) %>%
      
      pull(max_id) %>%
      unique())
#Add new column to E7 stage
condition_key_Embryo_E7$stage = 
  stringr::str_split_fixed(condition_key_Embryo_E7$max_id,
                           "_",
                           n = 2)[,1] %>%
  stringr::str_replace_all(pattern = "stage",
                           replacement = "")
#add Embryo E7 condition keys to cds
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
#Filter EB differentiation
EB3_cells_to_keep =
  colData(SP5_cds) %>%
  as.data.frame() %>%
  filter(!is.na(data_set),
         data_set == "EB differentiation") %>%
  pull(Cell) %>%
  as.character()

EB3_cds = SP5_cds[,EB3_cells_to_keep]

coldata_EB3_cds = colData(EB3_cds) %>% as_tibble()

#Filter embryo data: whole embryo
Embryo_cells_to_keep =
  colData(SP5_cds) %>%
  as.data.frame() %>%
  filter(!is.na(data_set),
         data_set == "Embryo") %>%
  pull(Cell) %>%
  as.character()
Emb1_cds = SP5_cds[,Embryo_cells_to_keep]
coldata_Emb1_cds = colData(Emb1_cds) %>% as_tibble()

#Filter embryo data: E7
E7_cells_to_keep =
  colData(Emb1_cds) %>%
  as.data.frame() %>%
  filter(!is.na(stage),
         stage == "E7") %>%
  pull(Cell) %>%
  as.character()
E7_cds = Emb1_cds[,E7_cells_to_keep]
coldata_E7_cds = colData(E7_cds) %>% as_tibble()

#Filter embryo data: sm 12-13 (presumably somite pair 12-13 about E8.5 to E9)
sm_late_cells_to_keep =
  colData(Emb1_cds) %>%
  as.data.frame() %>%
  filter(!is.na(sm),
         sm == "12-13") %>%
  pull(Cell) %>%
  as.character()
Emb_sm_late_cds = Emb1_cds[, sm_late_cells_to_keep]
coldata_Emb_sm_late_cds = colData(Emb_sm_late_cds) %>% as_tibble()

#Filter embryo data: sm 7-8 (presumably somite pair 7-8 about E8)
sm_early_cells_to_keep =
  colData(Emb1_cds) %>%
  as.data.frame() %>%
  filter(!is.na(sm),
         sm == "7-8") %>%
  pull(Cell) %>%
  as.character()
Emb_sm_early_cds = Emb1_cds[, sm_early_cells_to_keep]
coldata_Emb_sm_early_cds = colData(Emb_sm_early_cds) %>% as_tibble()

#Filter EP data from YS: sm 7-8
EP_cells_to_keep =
  colData(Emb1_cds) %>%
  as.data.frame() %>%
  filter(!is.na(part),
          part == "EP") %>%
  pull(Cell) %>%
  as.character()
EP_cds = Emb1_cds[, EP_cells_to_keep]
coldata_EP_cds = colData(EP_cds) %>% as_tibble()

#Filter sm data both embryo and yolksac
EYS_cells_to_keep =
  colData(Emb1_cds) %>%
  as.data.frame() %>%
  filter(!is.na(part), part %in% c("EP", "YS")) %>%
  pull(Cell) %>%
  as.character()
EYS_cds = Emb1_cds[, EYS_cells_to_keep]
coldata_EYS_cds = colData(EYS_cds) %>% as_tibble()

#Load Barb's previous analysis to compare
Barb_SP5_cds = 
  readRDS("/Users/adamheck/Desktop/mESC-SciPlex/processed_data/BVF_CDS/SP5.RDS")


####PROCESS AND CLUSTER CDS####
# Process for EB combined with Embryo cds
SP5_cds <- preprocess_cds(SP5_cds, num_dim = 50) #120 is num_dim for 12-2-20 subset data and most day_cds
SP5_cds <- reduce_dimension(SP5_cds, umap.n_neighbors = 4L) #10L for 12-2-20 subset and most day_cds but lower no. makes closer together. Used 4L for EB_cds
set.seed (17)
SP5_cds <- cluster_cells(SP5_cds, resolution = 3e-4, random_seed = 17) #Resolution - Lower no, fewer clusters - used 3e-3 for 12-2-20 subset data - used 3e-4 for cds
plot_cells(SP5_cds)
plot_cells(EB3_cds)
plot_cells(Barb_SP5_cds)
ggsave("results/PLOTS/SP5_cluster.tiff", height = 6, width = 6, dpi = 300, bg = "transparent")
plot_cells(SP5_cds) +
  facet_wrap(~data_set + ~part, ncol = 2)
ggsave("PLOTS/SP5_parts_cluster.tiff", height = 4, width = 4, dpi = 300, bg = "transparent")
plot_cells(SP5_cds, genes = c("Apoa2"), color_cells_by = "cluster", label_cell_groups = T, show_trajectory_graph = FALSE, label_branch_points = FALSE, label_leaves = FALSE,
           graph_label_size = 0.75, cell_size = 0.5) +
  facet_wrap(~data_set + ~part, ncol = 2)
plot_cells(SP5_cds, genes = Lateral_plate_Choi, color_cells_by = "cluster", label_cell_groups = T, show_trajectory_graph = FALSE, label_branch_points = FALSE, label_leaves = FALSE,
           graph_label_size = 0.75, cell_size = 0.5) 
  
#### Process EB cds with set seed ####
EB3_cds <- preprocess_cds(EB3_cds, num_dim = 98) #120 is num_dim for 12-2-20 subset data and most day_cds
EB3_cds <- reduce_dimension(EB3_cds, umap.n_neighbors = 4L) #10L for 12-2-20 subset and most day_cds but lower no. makes closer together. Used 4L for EB_cds
#cluster cells
set.seed(17)
EB3_cds = cluster_cells(EB3_cds, resolution = 3.6e-4, random_seed = 17)
#### Add umap data to cds ####
colData(EB3_cds)$umap1 = reducedDim(EB3_cds, type = "UMAP")[,1]
colData(EB3_cds)$umap2 = reducedDim(EB3_cds, type = "UMAP")[,2]
colData(EB3_cds)$Cluster = clusters(EB3_cds, reduction_method = "UMAP")

#Simple plot
ggsave("results/PLOTS/EB3_cluster.tiff", height = 6, width = 6, dpi = 300, bg = "transparent")
plot_cells(EB3_cds) +
  facet_wrap(~data_set + ~part, ncol = 2)

plot_cells(EB3_cds,
           color_cells_by = "cluster",
           label_cell_groups = T,
           show_trajectory_graph = FALSE,
           label_branch_points = FALSE,
           label_leaves = FALSE,
           graph_label_size = 0.75,
           cell_size = 0.2) +
  facet_grid(Activin~BMP) 
  ggsave("PLOTS/EB3_4x4.tiff", height = 6, width = 6, dpi = 300, bg = "transparent")  
plot_cells(EB3_cds, genes = c("Mixl1", "Cxcr4", "Etv2", "Dll4", "Gja5", "Lyve1", "Hoxb6", "Hoxa9"),
           color_cells_by = "cluster",
           label_cell_groups = T,
           show_trajectory_graph = FALSE,
           label_branch_points = FALSE,
           label_leaves = FALSE,
           graph_label_size = 0.75,
           cell_size = 0.2) 
  facet_grid(Activin~BMP)
ggsave("PLOTS/EB2_genes.tiff", height = 4, width = 6, dpi = 300, bg = "transparent")


# Process Embryo cds 
Emb1_cds <- preprocess_cds(Emb1_cds, num_dim = 8) #120 is num_dim for 12-2-20 subset data and most day_cds
Emb1_cds <- reduce_dimension(Emb1_cds, umap.n_neighbors = 4L) #10L for 12-2-20 subset and most day_cds but lower no. makes closer together. Used 4L for EB_cds
set.seed (17)
Emb1_cds <- cluster_cells(Emb1_cds, resolution = 5e-3, random_seed = 17) #Resolution - Lower no, fewer clusters - used 3e-3 for 12-2-20 subset data - used 3e-4 for cds
plot_cells(Emb1_cds)

plot_cells(Emb1_cds) +
  facet_wrap(~part, ncol = 3)
colData(Emb1_cds)$umap1 = reducedDim(Emb1_cds, type = "UMAP")[,1]
colData(Emb1_cds)$umap2 = reducedDim(Emb1_cds, type = "UMAP")[,2]
colData(Emb1_cds)$Cluster = clusters(Emb1_cds, reduction_method = "UMAP")
# Process EYS_cds
EYS_cds <- preprocess_cds(EYS_cds, num_dim =40) #120 is num_dim for 12-2-20 subset data and most day_cds
EYS_cds <- reduce_dimension(EYS_cds, umap.n_neighbors = 4L) #10L for 12-2-20 subset and most day_cds but lower no. makes closer together. Used 4L for EB_cds
set.seed (17)
EYS_cds <- cluster_cells(EYS_cds, resolution = 2e-3, random_seed = 17) #Resolution - Lower no, fewer clusters - used 3e-3 for 12-2-20 subset data - used 3e-4 for cds
plot_cells(EYS_cds)

plot_cells(EYS_cds) +
  facet_wrap(~part, ncol = 2)
colData(EYS_cds)$umap1 = reducedDim(EYS_cds, type = "UMAP")[,1]
colData(EYS_cds)$umap2 = reducedDim(EYS_cds, type = "UMAP")[,2]
colData(EYS_cds)$Cluster = clusters(EYS_cds, reduction_method = "UMAP")
# Process EP cds
EP_cds <- preprocess_cds(EP_cds, num_dim = 10) #120 is num_dim for 12-2-20 subset data and most day_cds
EP_cds <- reduce_dimension(EP_cds, umap.n_neighbors = 4L) #10L for 12-2-20 subset and most day_cds but lower no. makes closer together. Used 4L for EB_cds
EP_cds <- cluster_cells(EP_cds, resolution = 3e-4) #Resolution - Lower no, fewer clusters - used 3e-3 for 12-2-20 subset data - used 3e-4 for cds
plot_cells(EP_cds)
plot_cells(EP_cds, genes = c("Lyve1"), color_cells_by = "cluster",
           label_cell_groups = T,
           show_trajectory_graph = FALSE,
           label_branch_points = FALSE,
           label_leaves = FALSE,
           graph_label_size = 0.75,
           cell_size = 0.5) +
  facet_wrap(~sm, ncol = 3)


ggsave("PLOTS/E_YS_cluster.tiff", height = 4, width = 4, dpi = 300, bg = "transparent")
plot_cells(EYS_cds, genes = c("Cdh5", "Kdr", "Cxcr4", "Dll4", "CD44", "Runx1", "Procr", "Gja5", "Nr2f2", "Lyve1", "Myb", "Hlf", "Flt3", "Il7r", "Gfi1", "Ptprc", "Ahnak"), color_cells_by = "cluster",
           label_cell_groups = T,
           show_trajectory_graph = FALSE,
           label_branch_points = FALSE,
           label_leaves = FALSE,
           graph_label_size = 0.75,
           cell_size = 0.5) 
ggsave("PLOTS/EB2_HE.tiff", height = 4, width = 6, dpi = 300, bg = "transparent")  
#genes for Somitic, Somitic, Intermediate, Paraxial, Pharyngeal, Exe mesoderm, Allantois, Mesenchyme, Cardiomyocytes, Gut, HE prog 
plot_cells(SP5_cds, genes = c("Pcdh19", "Aldh1a2", "Osr1", "Meox1", "Alx1", "Hoxd1", "Tmem119", "Ahnak", "Tnnt2", "Gpx2", "Etv2", "Fev"), color_cells_by = "cluster",
           label_cell_groups = T,
           show_trajectory_graph = FALSE,
           label_branch_points = FALSE,
           label_leaves = FALSE,
           graph_label_size = 0.75,
           cell_size = 0.5)
ggsave("PLOTS/SP5_mesoderm.tiff", height = 4, width = 6, dpi = 300, bg = "transparent")
  facet_wrap(~part, ncol = 2)

# Process E7 cds
E7_cds <- preprocess_cds(E7_cds, num_dim = 20) #120 is num_dim for 12-2-20 subset data and most day_cds
E7_cds <- reduce_dimension(E7_cds, umap.n_neighbors = 4L) #10L for 12-2-20 subset and most day_cds but lower no. makes closer together. Used 4L for EB_cds
E7_cds <- cluster_cells(E7_cds, resolution = 3e-4) #Resolution - Lower no, fewer clusters - used 3e-3 for 12-2-20 subset data - used 3e-4 for cds
plot_cells(E7_cds)
# Process sm 12-13 or sm late cds
Emb_sm_late_cds <- preprocess_cds(Emb_sm_late_cds, num_dim = 20) #120 is num_dim for 12-2-20 subset data and most day_cds
Emb_sm_late_cds <- reduce_dimension(Emb_sm_late_cds, umap.n_neighbors = 4L) #10L for 12-2-20 subset and most day_cds but lower no. makes closer together. Used 4L for EB_cds
Emb_sm_late_cds <- cluster_cells(Emb_sm_late_cds, resolution = 3e-4) #Resolution - Lower no, fewer clusters - used 3e-3 for 12-2-20 subset data - used 3e-4 for cds
plot_cells(Emb_sm_late_cds)
plot_cells(Emb_sm_late_cds) +
  facet_wrap(~part, ncol = 2) 
# Process sm7-8 or sm early cds
Emb_sm_early_cds <- preprocess_cds(Emb_sm_early_cds, num_dim = 20) #120 is num_dim for 12-2-20 subset data and most day_cds
Emb_sm_early_cds <- reduce_dimension(Emb_sm_early_cds, umap.n_neighbors = 4L) #10L for 12-2-20 subset and most day_cds but lower no. makes closer together. Used 4L for EB_cds
Emb_sm_early_cds <- cluster_cells(Emb_sm_early_cds, resolution = 3e-4) #Resolution - Lower no, fewer clusters - used 3e-3 for 12-2-20 subset data - used 3e-4 for cds
plot_cells(Emb_sm_early_cds)
plot_cells(Emb_sm_early_cds) +
  facet_wrap(~part, ncol = 2) 

plot_cells(SP5_cds, genes = c("Cxcr4", "Stab2", "Kdr", "Lyve1", "Gja5", "Dll4", "Etv2", "Cdh5", "Ahnak", "Aldh1a2"), color_cells_by = "cluster",
           label_cell_groups = T,
           show_trajectory_graph = FALSE,
           label_branch_points = FALSE,
           label_leaves = FALSE,
           graph_label_size = 0.75,
           cell_size = 0.5)
plot_cells(Emb1_cds, genes = c("Cxcr4", "Lhx1", "Aldh1a2", "Hoxd1", "Hoxd9"), color_cells_by = "cluster",
           label_cell_groups = T,
           show_trajectory_graph = FALSE,
           label_branch_points = FALSE,
           label_leaves = FALSE,
           graph_label_size = 0.75,
           cell_size = 0.5)
plot_cells(EB3_cds, genes = c("Hoxb6"), color_cells_by = "cluster",
           label_cell_groups = T,
           show_trajectory_graph = FALSE,
           label_branch_points = FALSE,
           label_leaves = FALSE,
           graph_label_size = 0.75,
           cell_size = 0.5) 
      facet_wrap(~sm + part + stage, ncol = 3) 
      
 ####Lateral Plate Mesoderm lineage markers Prummel KD Mosimann C 2019 Nature Comm Zebrafish####
      Lateral_plate_intro <- c("Hand1", "Hand2", "Osr1", "Foxf1", "Prrx1", "Mesp1", "Etv2", "Cxcr4", "Foxh1")
#drivers for early LPM induction - loss will allow later LPM including Hand2 exp contain drl enhancer
      Lateral_plate_drivers <- c("Foxh1", "Mixl1", "Eomes")
      Lateral_plate_additional <- c( "Gata4", "Mef2c", "Tbx5", "Nkx2-5", "mir-133", "Foxo1", "Etv2", "Klf2", "Tal1", "Lmo2",
                         "Erg", "Gata2", "Runx1", "Cxcr4", "Bmp4", "Fgf9", "Apela", "Inhba" )
      
plot_cells(EYS_cds, genes=Lateral_plate_additional,
                 color_cells_by = "cluster",
                 label_cell_groups = TRUE,
                 label_branch_points = FALSE,
                 show_trajectory_graph = FALSE,
                 label_leaves = FALSE,
                 graph_label_size = 0.75,
                 cell_size = 0.5) +
        no_axes ()
      ggsave("PLOTS/LPM_intro_Prummel_Mosimann.tiff", height = 4, width = 6, dpi = 300, bg = "transparent")
plot_cells (Emb1_cds) +
    facet_wrap(~sm + part + stage, ncol = 3)
plot_cells(Emb1_cds, genes = Lateral_plate_Choi, color_cells_by = "cluster", label_cell_groups = T, show_trajectory_graph = FALSE,
                 label_branch_points = FALSE, label_leaves = FALSE,
                 graph_label_size = 0.75, cell_size = 0.5) 
plot_cells(Emb1_cds, genes = Lateral_plate_intro, color_cells_by = "cluster", label_cell_groups = T, show_trajectory_graph = FALSE,
           label_branch_points = FALSE, label_leaves = FALSE,
           graph_label_size = 0.75, cell_size = 0.5) 
    facet_wrap(~part, ncol = 3)
    facet_wrap(~sm + part + stage, ncol = 3)
plot_cells(Emb1_cds, genes = c("Krt18"), color_cells_by = "cluster", label_cell_groups = T, show_trajectory_graph = FALSE,
               label_branch_points = FALSE, label_leaves = FALSE,
               graph_label_size = 0.75, cell_size = 0.5) +
    facet_wrap(~part, ncol = 3)
    facet_wrap(~sm + part + stage, ncol = 3)
plot_cells(EB3_cds, genes = c("Krt18"), color_cells_by = "cluster", label_cell_groups = T, show_trajectory_graph = FALSE,
               label_branch_points = FALSE, label_leaves = FALSE,
               graph_label_size = 0.75, cell_size = 0.2) +
    facet_grid(Activin~BMP)
    facet_wrap(~part, ncol = 3)
    facet_wrap(~sm + part + stage, ncol = 3)
  ggsave("PLOTS/EB3_Krt8_4x4.tiff", height = 4, width = 6, dpi = 300, bg = "transparent")
facet_wrap(~sm + part + stage, ncol = 3)
facet_wrap(part, ncol = 2)
####Load and or save cds data frames####
#Save SP5_cds: has both EB and Embryo data
saveRDS(SP5_cds, file.path("/Users/adamheck/Desktop/mESC-SciPlex/processed_data/AH_CDS_122023/", "SP5.RDS"))
#Save EB2_cds: has only EB data
saveRDS(EB3_cds, file.path("/Users/adamheck/Desktop/mESC-SciPlex/processed_data/AH_CDS_122023/", "EB3.RDS"))
#Save EYS_cds
saveRDS(EYS_cds, file.path("/Users/adamheck/Desktop/mESC-SciPlex/processed_data/AH_CDS_122023/", "EYS.RDS"))
#Save Emb1_cds
saveRDS(Emb1_cds, file.path("/Users/adamheck/Desktop/mESC-SciPlex/processed_data/AH_CDS_122023/", "Emb1.RDS"))
#load SP5_cds
SP5_cds =
  readRDS("/Users/hadlandlab/Desktop/sciPlex_HSC3_Paper/SP5.RDS")
#load EB2_cds
EB2_cds =
  readRDS("/Users/hadlandlab/Desktop/sciPlex_HSC3_Paper/EB2.RDS")
#load EB3_cds
EB3_cds =
  readRDS("/Users/hadlandlab/Desktop/1-13-22_EB_analysis_PAPER/EB3.RDS")
#load EYS_cds
EYS_cds =
  readRDS("/Users/hadlandlab/Desktop/1-13-22_EB_analysis_PAPER/EYS.RDS")
#load Emb1_cds
Emb1_cds = readRDS("/Users/hadlandlab/Desktop/1-13-22_EB_analysis_PAPER/Emb1.RDS")
#load subsetted EB mesoderm mesenchyme cds
EB_M_M_cds = readRDS("/Users/hadlandlab/Desktop/1-13-22_EB_analysis_PAPER/EB_M_M.RDS")

#save files eg. dataframes
save(dfPaper, file = "dfPaper.Rdata")
load("/Users/hadlandlab/Desktop/1-13-22_EB_analysis_PAPER/dfPaper.Rdata")
save(dfPaper_1, file = "dfPaper_1.Rdata")
load("/Users/hadlandlab/Desktop/1-13-22_EB_analysis_PAPER/dfPaper_1.Rdata")
#If plots are no longer plotting try
dev.off()
