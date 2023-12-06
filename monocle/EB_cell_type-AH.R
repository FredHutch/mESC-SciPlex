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
#Brandon loads also  
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
  library(garnett)
  setwd("/Users/amheck/Desktop/mESC-SciPlex")
  
  DelayedArray:::set_verbose_block_processing(TRUE)
  options(DelayedArray.block.size=1000e7)})

####THEMES####
clear_theme = 
  function () {
    theme(
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent"))}

no_axes =
  function () {
    theme(axis.line.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())}

white_theme =   
  function () {
    theme(text = element_text(color = "white"), 
          axis.line.x = element_line(color = "white"),
          axis.line.y = element_line(color = "white"),
          axis.ticks = element_line(color = "white"),
          axis.text = element_text(color = "white"),
          strip.background = element_blank(),
          strip.text.x = element_text(color = "white"))}

#Add hashTable to Rsession and designate column names
#In this hashTable all cells are the same axis (1) and the same sample (sciChem_HSC_3)
#count is column indicating how many RNAs in that cell, oligo designates a particular group
hashTable1 = 
  read.table(file = "/Users/amheck/Desktop/mESC-SciPlex/Sanjay_code-files/hashTable.out",
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
  readRDS("/Users/amheck/Desktop/mESC-SciPlex/BVF_CDS/SP5.RDS")
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
###ERROR here currently for me######
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

####Save and load EB3_cds: has only EB data####
saveRDS(object = EB3_cds, file = "EB3.RDS")
EB3_cds = readRDS("/Users/amheck/Desktop/mESC-SciPlex/BVF_CDS/EB3.RDS")

####PLOT CLUSTERS EB_CDS (Figure2a)#####    
plot_cells(EB3_cds,
           group_cells_by = "cluster",
           color_cells_by = "cluster",
           show_trajectory_graph = F,
           label_cell_groups = T,
           group_label_size = 3.0,
           label_groups_by_cluster = T,
           label_branch_points = T,
           label_roots = T,
           label_leaves = T,
           graph_label_size = 0.1,
           cell_size = 0.2) +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        text = element_text(size=10)) 

ggsave("PLOTS/EB3_cds_cluster.tiff", width = 4, height = 4, dpi = 300, bg = "transparent")

#### Generate a heat maps to designate a cell types to a cluster (Supplemental Figure?)#######
#load RColorBrewer
library(RColorBrewer)
#Read in gene short names for gene anotations
genes = 
  read.table(file = "/Users/amheck/Desktop/mESC-SciPlex/Sanjay_code-files/gene.annotations",
             sep = "\t",
             header = F,
             col.names = c("id", "gene_short_name"))

#Gene markers that describe designated cell type
type_genes <- c(
  #epiblast#  
  "Dnmt3b", "Pou5f1", "Epcam", "Utf1", 
  #primitive streak#
  "Foxd3", "Nanog",
  #APS&primitive streak#
  "Eomes", "Hhex", "Otx2", "Mixl1", "Gsc", "Lhx1", "Pax7", "Foxa2",  
  #Caudal Epiblast#                
  "Nkx1-2", 
  #PGC
  "Ifitm3", "Dnd1", "Dppa3",
  #Notochord   
  "T", "Chrd", "Noto",       
  #Nascent mesoderm# 
  "Mesp1", "Mesp2", "Lefty2",
  #Mixed mesoderm#
  "Isl1",   
  #Intermediate mesoderm/caudal#
  "Osr1", "Cdx2", "Hes7","Hoxb1", "Cdx1", "Gbx2",
  #Paraxial mesoderm#               
  "Tbx1", "Tcf15",
  #Somitic mesoderm#
  "Aldh1a2", "Dll1", "Tbx6",
  #Pharyngeal mesoderm                
  "Tcf21", 
  #Exe mesoderm, allantois, mesenchyme
  "Hoxa10", "Hoxa11", "Tbx4", "Bmp4", "Ahnak", "Pmp22", "Tcf21", "Osr1", "Cdx2",
  #Cardiomyocytes#                
  "Smarcd3", "Acta2", "Tagln", "Myl7", "Nkx2-5", "Myl4", "Tnnt2",  
  #HE progenitors#               
  "Etv2", "Kdr", "Lmo2", "Vwf", "Dll4", "Cxcr4",
  #Endothelium#                
  "Igf2", "Anxa5", "Pecam1","Cdh5", 
  #Blood progenitors#               
  "Runx1","Cited4","Gata1", 
  #Erythroid1,2,3                     
  "Hbb-bh1", "Hba-a2", "Hba-a1", "Gypa",
  #NMP#
  "Epha5", "Cdx4", "Hoxb9", "Sox2", "Irx3",
  #Rostral neural ectoderm#
  "Hesx1", "Six3", "Pax2",
  #Neural crest#
  "Sox9", "Pax3", "Tfap2a", "Dlx2", "Sox10", 
  #Ectoderm#
  "Foxg1", "Trp63", "Grhl2", "Grhl3",
  #Endoderms#
  "Dkk1", "Krt19", "Amot", "Spink1", "Emb", "Cystm1", "Apoe", "Apoa2", "Ttr", "Tfap2c", "Ascl2", "Elf5", "Sparc", "Plat", "Lamb1")

#subset cds for the gene markers listed above
EB3_cds_type <- EB3_cds [rowData(EB3_cds)$gene_short_name %in% type_genes,]

#Create a df of cell groups (clusters)
EB3_type_group_df <- tibble::tibble(cell=row.names(colData(EB3_cds_type)), EB3_type_group="clusters"(EB3_cds_type)[colnames(EB3_cds_type)])

#generates matrix using aggregate gene expression for cell group [cluster] (what is the gene expression value?????)
EB3_type_ag_mat <- aggregate_gene_expression(EB3_cds_type, gene_group_df=NULL, EB3_type_group_df, norm_method = c("binary"),
                                             scale_agg_values = F)

#replaces gene id with gene short name#
row.names(EB3_type_ag_mat) = genes[match(row.names(EB3_type_ag_mat), genes$id),]$gene_short_name 

### Organise matrix and add annotation###
sub_type <- structure(EB3_type_ag_mat, ".dimnames" = list(c("Dnmt3b", "Pou5f1", "Epcam", "Utf1",
                                                            "Foxd3", "Nanog",
                                                            "Eomes", "Hhex", "Otx2", "Mixl1", "Gsc", "Lhx1", "Pax7", "Foxa2",
                                                            "Nkx1-2",
                                                            "Ifitm3", "Dnd1", "Dppa3",
                                                            "T", "Chrd", "Noto",
                                                            "Mesp1", "Mesp2", "Lefty2",
                                                            "Isl1",
                                                            "Osr1", "Cdx2", "Hes7","Hoxb1", "Cdx1", "Gbx2",
                                                            "Tbx1", "Tcf15",
                                                            "Aldh1a2", "Dll1", "Tbx6",
                                                            "Tcf21",
                                                            "Hoxa10", "Hoxa11", "Tbx4", "Bmp4", "Ahnak", "Pmp22", "Osr1", "Cdx2",
                                                            "Smarcd3", "Acta2", "Tagln", "Myl7", "Nkx2-5", "Myl4", "Tnnt2",
                                                            "Etv2", "Kdr", "Lmo2", "Vwf", "Dll4", "Cxcr4",
                                                            "Igf2", "Anxa5", "Pecam1","Cdh5",
                                                            "Runx1","Cited4","Gata1",
                                                            "Hbb-bh1", "Hba-a2", "Hba-a1", "Gypa",
                                                            "Epha5", "Cdx4", "Hoxb9", "Sox2", "Irx3",
                                                            "Hesx1", "Six3", "Pax2",
                                                            "Sox9", "Pax3", "Tfap2a", "Dlx2", "Sox10",
                                                            "Foxg1", "Trp63", "Grhl2", "Grhl3",
                                                            "Dkk1", "Krt19", "Amot", "Spink1", "EB3", "Cystm1", "Apoe", "Apoa2", "Ttr", "Tfap2c", "Ascl2", "Elf5", "Sparc", "Plat", "Lamb1")))

sub_anno_type <- structure(list("cell_type" = c("Epiblast", "Epiblast", "Epiblast", "Epiblast",
                                                "PS", "PS",
                                                "PS/APS", "PS/APS", "PS/APS", "PS/APS", "PS/APS", "PS/APS", "PS/APS", "PS/APS",
                                                "Caudal epiblast",
                                                "PGC", "PGC", "PGC",
                                                "Notochord", "Notochord", "Notochord",
                                                "Nascent mesoderm", "Nascent mesoderm", "Nascent mesoderm",
                                                "Mixed mesoderm",
                                                "Intermediate/Caudal mesoderm", "Intermediate/Caudal mesoderm", "Intermediate/Caudal mesoderm", "Intermediate/Caudal mesoderm", "Intermediate/Caudal mesoderm", "Intermediate/Caudal mesoderm",
                                                "Paraxial mesoderm", "Paraxial mesoderm",
                                                "Somitic mesoderm", "Somitic mesoderm", "Somitic mesoderm",
                                                "Pharyngeal mesoderm",
                                                "Exe_mesoderm/Mesenchyme/Allantois", "Exe_mesoderm/Mesenchyme/Allantois", "Exe_mesoderm/Mesenchyme/Allantois", "Exe_mesoderm/Mesenchyme/Allantois", "Exe_mesoderm/Mesenchyme/Allantois", "Exe_mesoderm/Mesenchyme/Allantois", "Exe_mesoderm/Mesenchyme/Allantois", "Exe_mesoderm/Mesenchyme/Allantois",
                                                "Cardiomyocytes", "Cardiomyocytes", "Cardiomyocytes", "Cardiomyocytes", "Cardiomyocytes", "Cardiomyocytes", "Cardiomyocytes", 
                                                "HE progenitors","HE progenitors","HE progenitors","HE progenitors","HE progenitors","HE progenitors",
                                                "Endothelium", "Endothelium", "Endothelium", "Endothelium",
                                                "Blood progenitors", "Blood progenitors", "Blood progenitors",
                                                "Erythroid", "Erythroid", "Erythroid", "Erythroid",
                                                "NMP", "NMP", "NMP", "NMP", "NMP",
                                                "Rostral neural ectoderm", "Rostral neural ectoderm", "Rostral neural ectoderm",
                                                "Neural crest", "Neural crest", "Neural crest", "Neural crest", "Neural crest",
                                                "Ectoderm", "Ectoderm", "Ectoderm", "Ectoderm",
                                                "Endoderm", "Endoderm", "Endoderm", "Endoderm", "Endoderm", "Endoderm", "Endoderm", "Endoderm", "Endoderm", "Endoderm", "Endoderm", "Endoderm", "Endoderm", "Endoderm", "Endoderm")),
                           .Names = "cell_type",
                           "row.names" = c("Dnmt3b", "Pou5f1", "Epcam", "Utf1",
                                           "Foxd3", "Nanog",
                                           "Eomes", "Hhex", "Otx2", "Mixl1", "Gsc", "Lhx1", "Pax7", "Foxa2",
                                           "Nkx1-2",
                                           "Ifitm3", "Dnd1", "Dppa3",
                                           "T", "Chrd", "Noto",
                                           "Mesp1", "Mesp2", "Lefty2",
                                           "Isl1",
                                           "Osr1", "Cdx2", "Hes7","Hoxb1", "Cdx1", "Gbx2",
                                           "Tbx1", "Tcf15",
                                           "Aldh1a2", "Dll1", "Tbx6",
                                           "Tcf21",
                                           "Hoxa10", "Hoxa11", "Tbx4", "Bmp4", "Ahnak", "Pmp22", "Osr1", "Cdx2",
                                           "Smarcd3", "Acta2", "Tagln", "Myl7", "Nkx2-5", "Myl4", "Tnnt2",
                                           "Etv2", "Kdr", "Lmo2", "Vwf", "Dll4", "Cxcr4",
                                           "Igf2", "Anxa5", "Pecam1","Cdh5",
                                           "Runx1","Cited4","Gata1",
                                           "Hbb-bh1", "Hba-a2", "Hba-a1", "Gypa",
                                           "Epha5", "Cdx4", "Hoxb9", "Sox2", "Irx3",
                                           "Hesx1", "Six3", "Pax2",
                                           "Sox9", "Pax3", "Tfap2a", "Dlx2", "Sox10",
                                           "Foxg1", "Trp63", "Grhl2", "Grhl3",
                                           "Dkk1", "Krt19", "Amot", "Spink1", "Emb", "Cystm1", "Apoe", "Apoa2", "Ttr", "Tfap2c", "Ascl2", "Elf5", "Sparc", "Plat", "Lamb1"),
                           class="data.frame")
sub_samp_ordered <- sub_type[row.names(sub_anno_type),]
pheatmap::pheatmap(sub_samp_ordered, annotation_row = sub_anno_type,
                   color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(25),
                   annotation_legend = F, cellwidth = 6, cellheight= 6,
                   cluster_rows=F, cluster_cols=T, fontsize = 5,
                   filename = "EB3_Heat_map_for_cell_type_1_13_22_binary_11_scale_F.tiff",
                   width = 8, height = 11, dpi = 300, bg = "transparent")

####Assign cell types using heatmap of clusters and marker expression - Figure2d####
EB3_assigned_cell_type_cds <- EB3_cds
colData(EB3_assigned_cell_type_cds)$assigned_cell_type <- as.character(clusters(EB3_assigned_cell_type_cds))
colData(EB3_assigned_cell_type_cds)$assigned_cell_type = dplyr::recode(colData(EB3_assigned_cell_type_cds)$assigned_cell_type,
        "4" = "Primitive streak", "5" = "Primitive streak", "13" = "Primitive streak", "14" = "Primitive streak", "23" = "Primitive streak",
        "17" = "Epiblast", "3" = "Epiblast",
        "10" = "Exe_mesoderm/Allantois", "8" = "Exe_mesoderm/Allantois", "6" = "Exe_mesoderm/Allantois",
        "7" = "Mesoderm", "22" = "Mesoderm", "30" = "Mesoderm", "1" = "Mesoderm",
        "21" = "HE progenitors", "35"= "HE progenitors", "15"= "HE progenitors", "18"= "HE progenitors", "26"= "HE progenitors",
        "25" = "HE progenitors", "33" = "HE progenitors", "39" = "HE progenitors", "24" = "HE progenitors", "29" = "HE progenitors",
        "38" = "Blood progenitors", "27" = "Blood progenitors", "32" = "Blood progenitors", "2" = "Blood progenitors",
        "9" = "Erythroid", "31" = "Erythroid", "36" = "Erythroid",
        "20" = "Undifferentiated ESC", "19" = "Undifferentiated ESC",
        "37" = "Visceral endoderm", "34" = "Visceral endoderm",
        "44" = "Definitive endoderm",
        "16" = "Ectoderm", "41" = "Ectoderm",
        "12" = "Mesenchyme", "40" = "Mesenchyme", "28" = "Mesenchyme",
        "42" = "PGC",
        "45" = "Cardiomyocytes",
        "11" = "Unknown", "43" = "Unknown", "46" = "Unknown")

#plot cells for celltype
plot_cells(EB3_assigned_cell_type_cds, group_cells_by = "cluster", color_cells_by = "assigned_cell_type", show_trajectory_graph = F, label_cell_groups = F,
           group_label_size = 4.5, label_groups_by_cluster =F, label_branch_points = F, label_roots = F, label_leaves = F,
           graph_label_size = 0.1, cell_size = 0.2) +
  
  scale_color_manual(values=c("Pink", "Maroon", "#877700", "Peru", "Khaki",
                              "Red", "Orange", "Turquoise", "SkyBlue",
                              "SlateBlue", "Yellow", "Orchid", "LightSteelBlue",
                              "SlateGray", "PeachPuff")) +
     
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.background = element_rect(color = "transparent", fill = "transparent", linetype=0)) +
  theme(panel.border = element_blank(),
        strip.background = element_blank() ,
        strip.text = element_blank()) 
  

ggsave("PLOTS/EB3_cds_assigned_cell_types_D6.tiff", width = 4, height = 4, dpi = 300, bg = "transparent")

facet_grid(Activin~BMP)

###Save assigned cell type CDS###
setwd("/Users/amheck/Desktop/mESC-SciPlex/BVF_CDS")
saveRDS(object = EB3_assigned_cell_type_cds, file = "EB3_celltype.RDS")

####Filter assigned cell type cds by day for figure 2#####
#Day4
EB3_assigned_cell_type_D4_cells_to_keep =
  colData(EB3_assigned_cell_type_cds) %>%
  as.data.frame() %>% 
  filter(!is.na(day),
         day %in% c("4")) %>%
  pull(Cell) %>%
  as.character() 
EB3_assigned_cell_type_D4_cds = EB3_assigned_cell_type_cds[,EB3_assigned_cell_type_D4_cells_to_keep]
#Day5
EB3_assigned_cell_type_D5_cells_to_keep =
  colData(EB3_assigned_cell_type_cds) %>%
  as.data.frame() %>% 
  filter(!is.na(day),
         day %in% c("5")) %>%
  pull(Cell) %>%
  as.character() 
EB3_assigned_cell_type_D5_cds = EB3_assigned_cell_type_cds[,EB3_assigned_cell_type_D5_cells_to_keep]
#Day6
EB3_assigned_cell_type_D6_cells_to_keep =
  colData(EB3_assigned_cell_type_cds) %>%
  as.data.frame() %>% 
  filter(!is.na(day),
         day %in% c("6")) %>%
  pull(Cell) %>%
  as.character() 
EB3_assigned_cell_type_D6_cds = EB3_assigned_cell_type_cds[,EB3_assigned_cell_type_D6_cells_to_keep]
#missing cardiomyocytes from D4 and D5 so use the following to color cell types!
plot_cells(EB3_assigned_cell_type_D5_cds, group_cells_by = "cluster", color_cells_by = "assigned_cell_type", show_trajectory_graph = F, label_cell_groups = F,
           group_label_size = 4.5, label_groups_by_cluster =F, label_branch_points = F, label_roots = F, label_leaves = F,
           graph_label_size = 0.1, cell_size = 0.2) +
  
  scale_color_manual(values=c("Pink", "#877700", "Peru", "Khaki",
                              "Red", "Orange", "Turquoise", "SkyBlue",
                              "SlateBlue", "Yellow", "Orchid", "LightSteelBlue",
                              "SlateGray", "PeachPuff")) +
  
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.background = element_rect(color = "transparent", fill = "transparent", linetype=0)) +
  theme(panel.border = element_blank(),
        strip.background = element_blank() ,
        strip.text = element_blank()) +
  no_axes()

ggsave("PLOTS/EB3_cds_assigned_cell_types_D5.tiff", width = 4, height = 4, dpi = 300, bg = "transparent")

####PLOT and FILTER based on EXPRESSION OF ONE VS ANOTHER GENE to CREATE DF with GENE EXPRESSION CHARACTERISTICS####
#Subset cds for genes expressed#
Paper_genes <- c("Cxcr4", "Cyp26a1", "Cyp26b1", "Aldh1a2", "Rbp1", "Crabp2", "Crabp1", "Rara", "Rarg", "Rarb",
                 "Cdx2", "Cdx1", "Hoxa9", "Hoxa11", "Arrb1", "Arrb2", "Clec1a", "Clec1b", "P2rx7",
                 "Pdgfra", "Kdr", "Bmp4", "Mixl1", "Eomes", "Lhx1", "Ahnak", "Etv2", "Itga2b",
                 "T", "Runx1", "Cd40", "Cdh5", "Cdh2", "Flt4", "Tek", "Flt1",
                 "Dll4", "Notch1", "Dlk1", "Gja5", "Vwf", "Sox17", "Gata4", "Foxf1", "Runx1",
                 "Cd44", "Procr", "Itgb3", "Itgb7", "Trim47", "Cd38", "Lyve1", "Stab2", "Ace")
Paper_cds <- EB3_assigned_cell_type_cds [rowData(EB3_assigned_cell_type_cds)$gene_short_name %in% Paper_genes,]
Paper_exprs <- SingleCellExperiment::counts(Paper_cds)
#normalize counts
Paper_exprs <- Matrix::t(Matrix::t(Paper_exprs)/size_factors(Paper_cds))
Paper_exprs <- reshape2::melt(as.matrix(Paper_exprs)) 
colnames(Paper_exprs) <- c("f_id", "Cell", "expression")
#make a matrix that has expression vs short_gene_name and cell and everything else in coldata - very useful!!!!!!
Paper_cds_coldata <- colData(Paper_cds)
#adds coldata info get error if do row data 1st
Paper_exprs <- merge(Paper_exprs, Paper_cds_coldata, by.x = "Cell", by.y = "row.names")
#add gene_short_name to dataframe or cds
#load gene annotations file
genes = read.table(file = "/Users/amheck/Desktop/mESC-SciPlex/Sanjay_code-files/gene.annotations", sep = "\t", header = F, col.names = c("f_id", "gene_short_name"))
#add column with gene short name
Paper_exprs <- merge(Paper_exprs, genes, by = "f_id")
#Create data frame from matrix and add umap cluster condition info  note: SPREAD puts gene_short_name in columns
dfPaper = data.frame(Paper_exprs)
dfPaper_1 = dfPaper %>% 
  dplyr::select(Cell,  expression,  gene_short_name, umap1,  umap2, Cluster, condition, day, assigned_cell_type)%>%
  spread(key = gene_short_name, value = expression, fill = 0)

####FILTER EXPRESSION LEVELS DATAFRAME to SUBSET DATAFRAME#####
#Figure2
#Primitive streak
dfPS_T =dfPaper_1 %>% filter(T>0.1)
#mesoderm
dfMixl1_Kdr = dfPaper_1 %>% filter(Mixl1>0.1 & Kdr>0.1)
#LPM
dfLPM_Kdr = dfPaper_1 %>% filter(Kdr>0.1)
dfLPM_Pdgfra_Kdr = dfPaper_1 %>% filter(Pdgfra>0.1 & Kdr>0.1)
dfLPM_Pdgfra_Kdr_neg = dfPaper_1 %>% filter(Pdgfra>0.1 & Kdr<0.1)
dfLPM_Pdgfra_Kdr_Cxcr4 = dfPaper_1 %>% filter(Pdgfra>0.1 & Kdr>0.1 & Cxcr4 > 0.1)
dfLPM_Mixl1_Cxcr4_Eomes_Lhx1 = dfPaper_1 %>% filter(Mixl1>0.1 & Cxcr4>0.1 & Eomes>0.1 & Lhx1>0.1)
dfLPM_Mixl1_Cxcr4 = dfPaper_1 %>% filter(Mixl1>0.1 & Cxcr4>0.1)
dfLPM_Kdr_Mixl1_Cxcr4_Pdgfra = dfPaper_1 %>% filter(Kdr>0.1 & Mixl1>0.1 & Cxcr4>0.1 & Pdgfra>0.1)
dfLPM_Dlk1 = dfPaper_1 %>% filter (Dlk1>0.1)
dfLPM_Dlk1_Kdr = dfPaper_1 %>% filter (Dlk1>0.1, Kdr>0.1)
#Endothelial
df_EC_Cdh5_Kdr = dfPaper_1 %>% filter(Cdh5>0.1 & Kdr>0.1)
#HE progenitors
dfHE_Kdr_Etv2 = dfPaper_1 %>% filter(Kdr>0.2 & Etv2>0.1)
dfHE_Kdr_Etv2_Cdh5 = dfPaper_1 %>% filter(Kdr>0.2 & Etv2>0.1 & Cdh5>0.1)
dfHE_Kdr_Etv2_Cdh5_Dll4 = dfPaper_1 %>% filter(Kdr>0.2 & Etv2>0.1 & Cdh5>0.1 & Dll4>0.1)
#HSC
dfHE_Cdh5_Dll4_Cxcr4 = dfPaper_1 %>% filter(Cxcr4>0.1 & Dll4>0.1 & Cdh5>0.1)
dfHE_Cdh5_Dll4_Cxcr4_Lyve1 = dfPaper_1 %>% filter(Cxcr4>0.1 & Dll4>0.1 & Cdh5>0.1 & Lyve1>0.1)
dfHE_Cdh5_Cxcr4_Lyve1 = dfPaper_1 %>% filter(Cxcr4>0.1 & Cdh5>0.1 & Lyve1>0.1)
dfHSC_Cxcr4_Gja5 = dfPaper_1 %>% filter(Cxcr4>0.1 & Gja5>0.1)
#AEC
dfAEC_Cdh5_Dll4 = dfPaper_1 %>% filter(Dll4>0.1 & Cdh5>0.1)
dfAEC_Dll4_Gja5_Cdh5 = dfPaper_1 %>% filter(Dll4>0.1 & Gja5>0.1 & Cdh5>0.1)
dfAEC_Dll4_Lyve1_Cdh5 = dfPaper_1 %>% filter(Dll4>0.1 & Lyve1>0.1 & Cdh5>0.1)
dfAEC_Dll4_Lyve1_Gja5_Cdh5 = dfPaper_1 %>% filter(Dll4>0.1 & Lyve1>0.1 & Gja5>0.1 & Cdh5>0.1)
dfAEC_Lyve1 = dfPaper_1 %>% filter(Lyve1>0.1)
dfAEC_Lyve1_Dll4neg_Cdh5 = dfPaper_1 %>% filter(Lyve1>0.1 & Dll4<0.1 & Cdh5>0.5)
dfAEC_Gja5_Dll4neg_Cdh5 = dfPaper_1 %>% filter(Dll4<0.1 & Gja5>0.1 & Cdh5>0.1)
dfAEC_Gja5_Cdh5 = dfPaper_1 %>% filter(Cdh5>0.1 & Gja5>0.1)
dfAEC_Lyve1_Cdh5 = dfPaper_1 %>% filter(Cdh5>0.1 & Lyve1>0.1)
dfAEC_Lyve1_Gja5_Cdh5 = dfPaper_1 %>% filter(Cdh5>0.1 & Lyve1>0.1 & Gja5>0.1)
dfAEC_Hoxa9_Cdh5 = dfPaper_1 %>% filter(Hoxa9>0.1 & Cdh5>0.1)
dfAEC_Hoxa9_Dll4_Cdh5 = dfPaper_1 %>% filter(Hoxa9>0.1 & Cdh5>0.1 & Dll4>0.1)
#RA
dfRA_Cxcr4_Cyp26a1 = dfPaper_1 %>% filter(Cxcr4>0.1 & Cyp26a1>0.1)
dfRA_Cyp26a1_Kdr = dfPaper_1 %>% filter(Kdr>0.1 & Cyp26a1>0.1)
dfRA_Cxcr4_Cyp26a1_Kdr = dfPaper_1 %>% filter(Cxcr4>0.1 & Cyp26a1>0.1 & Kdr>0.1)
dfRA_Cxcr4_Rarg = dfPaper_1 %>% filter(Cxcr4>0.1 & Rarg>0.1)
dfRA_Rarg = dfPaper_1 %>% filter(Rarg>0.1)
dfRA_Rarg_Kdr = dfPaper_1 %>% filter(Rarg>0.1 & Kdr>0.1)
dfRA_Rarg_Kdr_Cxcr4 = dfPaper_1 %>% filter(Rarg>0.1 & Kdr>0.1 & Cxcr4>0.1)
dfRA_Rarg_Crabp2 = dfPaper_1 %>% filter(Rarg>0.1 & Crabp2>0.1)
dfRA_Crabp2 = dfPaper_1 %>% filter(Crabp2>0.1)
dfRA_Crabp2_Kdr = dfPaper_1 %>% filter(Crabp2>0.1 & Kdr>0.1)
dfRA_Crabp2_Kdr_Cxcr4 = dfPaper_1 %>% filter(Crabp2>0.1 & Kdr>0.1 & Cxcr4>0.1)
dfRA_Rara_Kdr = dfPaper_1 %>% filter(Rara>0.1 & Kdr>0.1)
dfRA_Rara_Kdr_Cxcr4 = dfPaper_1 %>% filter(Rara>0.1 & Kdr>0.1 & Cxcr4>0.1)
dfRA_Rara = dfPaper_1 %>% filter(Rara>0.1)
dfRA_Rarg_Rara = dfPaper_1 %>% filter(Rarg>0.1 & Rara>0.1)
dfRA_Arrb1_Kdr = dfPaper_1 %>% filter(Arrb1>0.1 & Kdr>0.1)
dfRA_Arrb2_Kdr = dfPaper_1 %>% filter(Arrb2>0.1 & Kdr>0.1)
dfRA_Arrb1_Arrb2_Kdr = dfPaper_1 %>% filter(Arrb1>0.1 & Arrb2>0.1 & Kdr>0.1)

####Write dataframe to Excel (rstudio only shows 50 columns)####
install.packages("WriteXLS")
library(WriteXLS)
WriteXLS("dfAEC_Lyve1_Cdh5", ExcelFileName = "dfAEC_Lyve1_Cdh5.xlsx")

####PLOT df COLOR BY DAY of genes expressed single or 4x4 ####
#will not plot concentrations as 4x4 progression if do not do following - should be top df
dfLPM_Pdgfra_Kdr$condition <- factor(dfLPM_Pdgfra_Kdr$condition,
                                 levels=c("1","5","9","13",
                                          "2", "6","10","14",
                                          "3","7","11","15",
                                          "4","8","12","16")) 
#top
ggplot (data=dfLPM_Pdgfra_Kdr, mapping = aes(x=umap1, y=umap2, color = day, size = Pdgfra)) +
  scale_size_continuous(range = c(0.5, 0.5)) +
  scale_color_manual(values = c("4" = "#D41159",
                                "5" = "#0C7BDC",
                                "6" = "#FFC20A")) + geom_jitter() +
#bottom
  geom_point(data=dfPaper_1, aes(umap1, umap2), color = "gray", size = 0.05, alpha = 0.8) + geom_jitter() +
  
#add facet wrap here if want 4x4
  facet_wrap(~condition, nrow = 4) +  
  theme(plot.title=element_text(hjust=0.5, face='bold', color = 'black', size = 10)) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  no_axes()
ggsave("PLOTS/Fig2_Pdgfra_Kdr-neg_4x4.tiff", height = 3, width = 3, dpi = 300, bg = "transparent")

facet_wrap(~condition, nrow = 4)

####Subset Paper_1 according to cell type####
 
####PLOT BARS CELL NUMBER VS CONDITION USING DATAFRAME####
#.drop insures condition with 0 counts will be added to tibble
#na.rm = T in sum replaces NaN (results when divide by 0) with 0
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

####Following replaces NaN with 0 must run code above is.nan.data.frame etc before this can work
dfPaper_2[is.nan(dfPaper_2)] <- 0     

dfPaper_2 <- dfPaper_1 %>% 
  group_by (condition, day, .drop = F) %>% 
  dplyr::summarise (number_condition=dplyr::n(),
                    Pdgfra_Kdr_no=sum(Pdgfra>0.1 & Kdr>0.1), Pdgfra_Kdr_ratio = Pdgfra_Kdr_no/dplyr::n(), Pdgfra_Kdr_percent= Pdgfra_Kdr_ratio * 100,
                    Pdgfra_Kdr_neg_no=sum(Pdgfra>0.1 & Kdr<0.1), Pdgfra_Kdr_neg_ratio = Pdgfra_Kdr_neg_no/dplyr::n(), Pdgfra_Kdr_neg_percent= Pdgfra_Kdr_neg_ratio * 100,
                    Etv2_Kdr_Cdh5_no=sum(Cdh5>0.1 & Etv2>0.1 & Kdr>0.1), Etv2_Kdr_Cdh5_ratio = Etv2_Kdr_Cdh5_no/dplyr::n(), Etv2_Kdr_Cdh5_percent= Etv2_Kdr_Cdh5_ratio * 100,
                    Etv2_Kdr_no=sum( Etv2>0.1 & Kdr>0.1), Etv2_Kdr_ratio = Etv2_Kdr_no/dplyr::n(), Etv2_Kdr_percent= Etv2_Kdr_ratio * 100,
                    Kdr_Cxcr4_Mixl1_no=sum(Kdr>0.1 & Cxcr4 >0.1 & Mixl1 >0.1), Kdr_Cxcr4_Mixl1_ratio = Kdr_Cxcr4_Mixl1_no/dplyr::n(), EKdr_Cxcr4_Mixl1_percent= Kdr_Cxcr4_Mixl1_ratio * 100,
                    Cdh5_Dll4_Kdr_no=sum(Kdr>0.1 & Cdh5>0.1 & Dll4>0.1), Cdh5_Dll4_Kdr_ratio = Cdh5_Dll4_Kdr_no/dplyr::n(), Cdh5_Dll4_Kdr_percent= Cdh5_Dll4_Kdr_ratio * 100) %>% 
  dplyr::ungroup()

####Number  of cells / condition & assigned cell type - made heat map in prism####
dfPaper_3 <- dfPaper_1 %>% 
  group_by (condition, assigned_cell_type, day, .drop = F) %>% 
  dplyr::summarise (number_condition=dplyr::n()) %>% 
  dplyr::ungroup() 

####Write dataframe to Excel (plot with prism so all figures match####
install.packages("WriteXLS")
library(WriteXLS)
WriteXLS("dfPaper_2", ExcelFileName = "Condition_Day_dfPaper2.xlsx")

####If plots are no longer plotting try####
dev.off()

####Plot cells monocle####
#Zhao and Choi mesoderm-endoderm bifurcation
plot_cells(EB3_cds, genes = c("Sox17", "Mesp1", "Utf1", "Fgf5", "Foxa2", "Cer1", "Gsc", "Gata6", "T"), color_cells_by = "cluster", label_cell_groups = T, cell_size = 0.5)
#Zhao and Choi Flk1+ mesoderm and hemangiogenic lineage
plot_cells(EB3_cds, genes = c("Kdr", "Pdgfra", "Mesp1", "Tbx20", "Isl1", "Hand1", "Foxf1", "Etv2", "Tal1", "Gata2"), color_cells_by = "cluster", label_cell_groups = T, cell_size = 0.5) 
#Zhao and Choi LPM meta gene
plot_cells(EB3_cds, genes = c("Ahnak", "Ndufa4l2", "Gata6", "Cfc1", "Cxcr4", "Gata3", "Pdlim3", "Myl4", "Msx1", "Slc38a4", "Asb1", "Fgf3", "Gpc3", "Hmga2", "Lmo1", "Dkk1" ), color_cells_by = "cluster", label_cell_groups = T, cell_size = 0.5) 
plot_cells(EB3_cds, genes = c("Pecam1", "Cdh5", "Cdh1","Cdh2"), color_cells_by = "cluster", label_cell_groups = T, cell_size = 0.5) 
plot_cells(EB3_cds, genes = c("Rab15", "Arhgef19", "Apela", "Prss22"), color_cells_by = "cluster", label_cell_groups = T, cell_size = 0.5) 
plot_cells(EB3_cds, genes = c("Foxf1", "Hand1", "Bmp4", "Tagln", "Acta2", "Mesp1", "Kdr", "Etv1", "Gata2", "Lmo2", "Tal1"), color_cells_by = "cluster", label_cell_groups = T, cell_size = 0.5) 
plot_cells(EB3_cds, genes = c("Alox5", "Alox5ap"), color_cells_by = "cluster", label_cell_groups = T, cell_size = 0.5) 
ggsave("PLOTS/.tiff", height = 8, width = 10, dpi = 300, bg = "transparent") 
facet_grid(Activin~BMP)






