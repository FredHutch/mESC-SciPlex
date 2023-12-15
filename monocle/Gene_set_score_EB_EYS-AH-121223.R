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
  #library(garnett) #Do not have loaded
  
  library(WriteXLS)
  
  DelayedArray:::set_verbose_block_processing(TRUE)
  options(DelayedArray.block.size=1000e7)})

writeLines(capture.output(sessionInfo()), "sessionInfo-mESC-SciPlex_Genesetscore.txt")
# Set project directory
projectdir <- "/Users/adamheck/Desktop/mESC-SciPlex"
inputdir <- paste(projectdir, "processed_data/AH_CDS_122023", sep = "/")
plotdir <- paste(projectdir, "results/PLOTS", sep = "/")
outputdir <- paste(projectdir, "results", sep = "/")
setwd(projectdir)

####Load relevant cds####
#load SP5_cds - contains all sciPlex data 
SP5_cds = 
  readRDS("/Users/adamheck/Desktop/mESC-SciPlex/processed_data/AH_CDS_122023/SP5.RDS")
#load EB3_cds Filtered for EB UMAP added
EB3_cds = 
  readRDS("/Users/adamheck/Desktop/mESC-SciPlex/processed_data/AH_CDS_122023/EB3.RDS")
#load EYS_cds E8 and E9 separated into EP and YS no E7 has UMAP added
EYS_cds =
  readRDS("/Users/adamheck/Desktop/mESC-SciPlex/processed_data/AH_CDS_122023/EYS.RDS")

####THEMES####
#brandon's simple theme
simple_theme <-  theme(text = element_blank(), panel.grid = element_blank(),
                       axis.title = element_blank(),
                       axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       axis.line.x = element_blank(),
                       axis.line.y = element_blank())
my_theme <- theme(panel.border = element_blank(),
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
                  axis.ticks.y = element_blank())
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
setwd(plotdir)
####Plot single gene expression in Embryo/Yolk Sac for a given gene
#Set genes in list
gene_list <- c("Kdr", "Cdh5","Aldh1a2","Runx1","Rara","Rarb","Rarg","Hbb-bh1")  # Replace with actual gene names

for (gene in gene_list) {
  p <- plot_cells(EYS_cds, genes = gene, color_cells_by = "cluster",
                  label_cell_groups = FALSE,
                  show_trajectory_graph = FALSE,
                  label_branch_points = FALSE,
                  label_leaves = FALSE,
                  graph_label_size = 0.75,
                  cell_size = 0.75)
  # Save the plot
  ggsave(filename = paste0(gene, "_exp_full_EmbryoYolkSac.tiff"), plot = p, width = 6, height = 6, dpi = 300, bg = "transparent")
}
#Plot by part of embryo dissected (embryo proper or yolk sac)
plot_cells(EYS_cds, color_cells_by = "part",
           label_cell_groups = FALSE,
           show_trajectory_graph = FALSE,
           label_branch_points = FALSE,
           label_leaves = FALSE,
           graph_label_size = 0.75,
           cell_size = 0.5)

ggsave("EmbryoYS_fullcds_dissection_part.tiff", width = 6, height = 6, dpi = 300, bg = "transparent")

#Subset endothelial clusters
ec_cells =
  colData(EYS_cds) %>%
  as.data.frame() %>% 
  filter(!is.na(Cluster),
         Cluster %in% c("8", "10")) %>% 
  
  pull(Cell) %>%
  as.character() 
ec_cds = EYS_cds[,ec_cells]

#Plot Lyve1 exp in endothelial cells
plot_cells(ec_cds, genes = c("Lyve1"), color_cells_by = "cluster",
           label_cell_groups = FALSE,
           show_trajectory_graph = FALSE,
           label_branch_points = FALSE,
           label_leaves = FALSE,
           graph_label_size = 0.75,
           cell_size = 1)
ggsave("EmbryoYS_ECsub_Lyve1_exp.tiff", width = 6, height = 6, dpi = 300, bg = "transparent")

#Plot endothelial cells colored by YS/EP
plot_cells(ec_cds, color_cells_by = "part",
           label_cell_groups = FALSE,
           show_trajectory_graph = FALSE,
           label_branch_points = FALSE,
           label_leaves = FALSE,
           graph_label_size = 0.75,
           cell_size = 1)
ggsave("EmbryoYS_ECsub_dissection_part.tiff", width = 6, height = 6, dpi = 300, bg = "transparent")

#From endothelial cells, subset embryo
eb_cells =
  colData(ec_cds) %>%
  as.data.frame() %>% 
  filter(!is.na(part),
         part %in% c("EP")) %>% 
  
  pull(Cell) %>%
  as.character() 
ec_eb_cds = ec_cds[,eb_cells]

#Plot Lyve1 exp in embryo proper subset
plot_cells(ec_eb_cds, genes = c("Lyve1"), color_cells_by = "cluster",
           label_cell_groups = FALSE,
           show_trajectory_graph = FALSE,
           label_branch_points = FALSE,
           label_leaves = FALSE,
           graph_label_size = 0.75,
           cell_size = 1)
ggsave("EmbryoYS_ECsub_embryoproper_Lyve1_exp.tiff", width = 6, height = 6, dpi = 300, bg = "transparent")

#From endothelial cells, subset yolk sac
ys_cells =
  colData(ec_cds) %>%
  as.data.frame() %>% 
  filter(!is.na(part),
         part %in% c("YS")) %>% 
  
  pull(Cell) %>%
  as.character() 
ec_ys_cds = ec_cds[,ys_cells]

#Plot Lyve1 exp in yolk sac subset
plot_cells(ec_ys_cds, genes = c("Lyve1"), color_cells_by = "cluster",
           label_cell_groups = FALSE,
           show_trajectory_graph = FALSE,
           label_branch_points = FALSE,
           label_leaves = FALSE,
           graph_label_size = 0.75,
           cell_size = 1)
ggsave("EmbryoYS_ECsub_yolksac_Lyve1_exp.tiff", width = 6, height = 6, dpi = 300, bg = "transparent")

#########Estimate gene set scores##############
#Function for estimating gene set scores
#### Estimate gene set scores for signature gene lists
estimate_score <- function(cds, markers){
  cds_score = cds[fData(cds)$gene_short_name %in% markers,] 
  aggregate_score = exprs(cds_score)
  aggregate_score = Matrix::t(Matrix::t(aggregate_score) / pData(cds_score)$Size_Factor)
  aggregate_score = Matrix::colSums(aggregate_score)
  pData(cds)$score = log(aggregate_score +1) 
  return(cds)
}

###Specify markers for each gene set
#Primitive streak
PS <- c( "Wnt8a", "Trh", "Pim2")
#Lateral plate mesoderm
LPM <- c( "Mixl1", "Cxcr4", "Myl4", "Asb4", "Dkk1", "Lefty2")
#Extra Embryonic Cells (EEC) from Mesoderm
EEC <- c("Runx1", "Hapln1", "Fli1", "Ddah1", "Tal1", "Sox7", "Gata2", "Lmo2", "Col1a1", "Morc4", "Lgr5", "Lama4", "Dlk1")
#Intra Embryonic Cells (IEC) from Mesoderm
IEC <- c("Lefty2", "Lrig3", "Pcdh8", "Zic3", "Lhx1", "Dll1", "Pcdh19", "Sema6d", "Jag1", "Zic2", "Lfng", "Pmaip1", "Dll3", "Rimbp2")
#Mesenchyme
Mesenchyme <- c("Tdo2", "Colec11", "Gdf6", "Snai2", "Wisp1", "Ccdc80", "Sdpr", "Lum", "Postn", "Mab21l2")
#Cardiomyocytes
Cardiomyocytes <- c("Myl7", "Myl4", "Nkx2-5", "Tnnt2", "Mef2c", "Tnnc1", "Sfrp5", "Nexn", "Sh3bgr", "Csrp3", "Gm45123", "Cnn1", "Unc45b", "Ptges3l")
#Extra-embryonic mesoderm
EXE_mesoderm <- c("Hand1", "Col9a1", "Spin2c", "Cdx2", "Cdx4", "Hoxd11")
#Other mesodermal types
Nascent_mesoderm <- c("Lefty2", "Mesp2", "Pcdh8", "Prss22", "Cdh20", "Slpi", "Pmaip1", "Tmem229a", "Arl4d", "Rimbp2", "Tdgf1")
Somitic_mesoderm <- c("Aldh1a2", "Tbx6", "Rhbg", "Hes7")
Pharyngeal_mesoderm <- c("Gata5", "A330008L17Rik", "Nxf3", "Fibin", "Pnliprp1", "Tlx1")
Paraxial_mesoderm <- c( "Rspo1", "Ebf2", "Ptgfr", "Tbx1", "Meox2", "Ebf3")
### AEC genes from Kalucka et al 2020
AEC1_Kalucka <- c("8430408G22Rik", "Clu", "Crip1", "Fbln2", "Gja4", "Hey1", "Mecom", "Sat1", "Sema3g", "Sox17", "Tm4sf1", "Tsc22d1")
#AEC from Xu et al
AEC2_Xu <- c("Efnb2","Bmx","Gja5","Gja4","Cd200","Cxcl12","Kitl","Gkn3","Fbln2","Stmn2","Gper1","Tmem255b","C1qtnf9","Timp4","Tmem100","Sparcl1","Slc26a10","Rbp7","Alpl","Apq7","Hey1","Emp2","Lims2","Cav1","Ednrb","Dusp26","Ace","Car7","Pcp4l1","Edn1","Sema3g","Efhd1","Jam2","Cd34","Ramp3","Plat","Ltbp4","Ackr3","Col18a1","Serpinf1","Acta2","Podx1","Gprc5a","Lmcd1","Vegfc","Nxpe4","Igfbp3")


####Plot gene set scores in UMAP across entire EB dataset
#Creating list of genesets to be used
list_of_genesets <- list(Mesenchyme,Cardiomyocytes,EXE_mesoderm,AEC1_Kalucka,AEC2_Xu)
names(list_of_genesets) <- c("Mesenchyme_geneset","Cardiomyocytes_geneset", "EXE_mesoderm_geneset","AEC1_Kalucka_geneset","AEC2_Xu_geneset")
#Set colors for plotting gene set scores
mycol <- c("gray80", "gray80", "gray80", "gray80", "gray80", "red", "red4")
#For loop to generate gene set heat maps
for (geneset_name in names(list_of_genesets)) {
  # Extract the geneset from the list
  geneset <- list_of_genesets[[geneset_name]]
  # Calculate gene set score
  eb_cds <- estimate_score(EB3_cds, markers = geneset)
  #Create data frame for ordering the cells
  marker_coldata = colData(eb_cds) %>% as_tibble()
  #Reorder data frame so that cells with highest gene set score are on the top.  Reorder from lowest expression to highest.  Has to be a dataframe
  ordered_marker_coldata = marker_coldata[order(marker_coldata$score),]
  #Plot in UMAP space
  p <- ggplot (data=ordered_marker_coldata, mapping = aes(x=umap1, y=umap2, color = score, size = 0.25)) +
    scale_size_continuous(range = c(0.25, 0.25)) +
    #
    scale_color_gradientn(colours = mycol) +
    #bottom is entire data set
    geom_point(data=ordered_marker_coldata, aes(umap1, umap2), color = "gray", size = 0.25, alpha = 0.8) + geom_jitter() +
    #add facet wrap here by days
    #facet_wrap(~day, ncol = 4) +
    #Themes
    theme(legend.position = "left") +
    my_theme
  #Save the plot
  ggsave(filename = paste0(geneset_name, ".tiff"), plot = p, width = 4, height = 4, dpi = 300, bg = "white")
}

#Subset the mesoderm differentiation cells
EB_Mesoderm_Mesenchyme_cells_to_keep =
  colData(EB3_cds) %>%
  as.data.frame() %>% 
  filter(!is.na(Cluster),
         Cluster %in% c("32", "28","3","16","5","12","20","27","9","8","38","15","34","37","35","42","31","19","6")) %>% 
  
  pull(Cell) %>%
  as.character() 
meso_sub_cds = EB3_cds[,EB_Mesoderm_Mesenchyme_cells_to_keep]
#Plot meso subset
plot_cells(meso_sub_cds,
           group_cells_by = "cluster",
           color_cells_by = "cluster",
           show_trajectory_graph = F,
           label_cell_groups = T,
           group_label_size = 5.0,
           label_groups_by_cluster = T,
           label_branch_points = T,
           label_roots = T,
           label_leaves = T,
           graph_label_size = 0.1,
           cell_size = 0.5) +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        text = element_text(size=10)) 

ggsave("FIG2B_MesoSub_cds_cluster.tiff", width = 6, height = 6, dpi = 300, bg = "white")

#####Plot gene set scores in UMAP for meso subset, each day separately
#Creating list of genesets to be used
list_of_genesets <- list(PS,LPM,EEC,IEC)
names(list_of_genesets) <- c("Primitive_Streak_geneset_bydays","LateralPlateMesoderm_geneset_bydays", "ExtraEmbryonic_geneset_bydays","IntraEmbryonic_geneset_bydays")
#Set colors for plotting gene set scores
mycol <- c("gray80", "gray80", "gray80", "gray80", "gray80", "red","red2", "red4")
#For loop to generate gene set heat maps (facet wrap with day 4-6)
for (geneset_name in names(list_of_genesets)) {
  # Extract the geneset from the list
  geneset <- list_of_genesets[[geneset_name]]
  #Calculate gene set score
  gss_cds <- estimate_score(meso_sub_cds, markers = geneset) #Replace cds in estimate_scores function based on if you want all cells or a subset
  #Create data frame for ordering the cells
  coldata_gss_cds = colData(gss_cds) %>% as_tibble()
  #Reorder data frame so that cells with highest gene set score are on the top.  Reorder from lowest expression to highest.  Has to be a dataframe
  Ordered_gss_cds = coldata_gss_cds[order(coldata_gss_cds$score),]
  #Plot the gene set in UMAP space
  p <- ggplot (data=Ordered_gss_cds, mapping = aes(x=umap1, y=umap2, color = score)) +
    scale_size_continuous(range = c(0.5, 1)) +
    #
    scale_color_gradientn(colours = mycol) +
    #bottom is entire data set
    geom_point(data=Ordered_gss_cds, aes(umap1, umap2), color = "gray", size = 0.25, alpha = 0.8) + geom_jitter() +
    #add facet wrap here by days
    facet_wrap(~day, ncol = 4) +
    #Themes
    theme(legend.position = "left") +
    my_theme
  #Save the plot
  ggsave(filename = paste0(geneset_name, ".tiff"), plot = p, width = 8, height = 4, dpi = 300, bg = "white")
}

######Creating heatmaps for different geneset scores###########
#Generate PS score data
gss_cds <- estimate_score(meso_sub_cds, markers = PS)
gss_coldata = colData(gss_cds) %>% as_tibble()
colnames(gss_coldata)[colnames(gss_coldata) == "score"] <- "PS_score"
#Generate LPM score data
gss_cds <- estimate_score(meso_sub_cds, markers = LPM)
gss_coldata$LPM_score <- colData(gss_cds)$score
#Generate EEC score data
gss_cds <- estimate_score(meso_sub_cds, markers = EEC)
gss_coldata$EEC_score <- colData(gss_cds)$score
#Generate IEC score data
gss_cds <- estimate_score(meso_sub_cds, markers = IEC)
gss_coldata$IEC_score <- colData(gss_cds)$score

#Generate mid point scores for each gene set  
gss_MP <- gss_coldata %>% 
  dplyr::summarize (PS_MP=(((max(PS_score)-min(PS_score))/2) + min(PS_score)),
                    LPM_MP=(((max(LPM_score)-min(LPM_score))/2) + min(LPM_score)),
                    EEC_MP=(((max(EEC_score)-min(EEC_score))/2) + min(EEC_score)),
                    IEC_MP=(((max(IEC_score)-min(IEC_score))/2) + min(IEC_score)))

#Split into single day dataframes
day4_gss <- gss_coldata[gss_coldata$day == '4', ]
day5_gss <- gss_coldata[gss_coldata$day == '5', ]
day6_gss <- gss_coldata[gss_coldata$day == '6', ]


# Define the modified function to calculate the sum for cells greater than the threshold value, make the condition the row names, and transpose the dataframe
gss_cell_counts <- function(data, threshold_values) {
  result_df <- data.frame(condition = unique(data$condition))
  for (col in names(threshold_values)) {
    sum_results <- sapply(unique(data$condition), function(condition) {
      subset_data <- data[data$condition == condition, col]
      sum(subset_data > threshold_values[col])
    })
    result_df[[col]] <- sum_results
  }
  rownames(result_df) <- result_df$condition
  result_df$condition <- NULL
  result_df <- t(result_df)
  return(result_df)
}
# Define the threshold values for each column
threshold_values <- c(PS_score = gss_MP$PS_MP, LPM_score = gss_MP$LPM_MP, EEC_score = gss_MP$EEC_MP, IEC_score = gss_MP$IEC_MP)  # Replace with the actual threshold values for each column
#Run new function to calcualt cell counts for each day
result_df_day4 <- gss_cell_counts(day4_gss, threshold_values)
result_df_day5 <- gss_cell_counts(day5_gss, threshold_values)
result_df_day6 <- gss_cell_counts(day6_gss, threshold_values)

# Manually define the desired order for rows and columns
desired_row_order <- c("PS_score", "LPM_score", "IEC_score","EEC_score") # Replace with actual row names
desired_col_order <- c("1","5","9","13","2","6","10","14","3","7","11","15","4","8","12","16") # Replace with actual column names

#Manual set the order so it is easier to read through the table
result_df_day4 <- result_df_day4[desired_row_order, desired_col_order]
result_df_day5 <- result_df_day5[desired_row_order, desired_col_order]
result_df_day6 <- result_df_day6[desired_row_order, desired_col_order]

#Manually set the scale of the heatmap
breaks <- seq(0, 400, by = 0.1)  # Adjust the 'by' value as needed for finer or coarser color transitions
#Creat the color palette you want
color_palette <- colorRampPalette(c("gray80", "blue2", "green2","yellow2"))(length(breaks) - 1)
####For loop to genereate heatmaps
list_of_datasets <- list(result_df_day4,result_df_day5,result_df_day6)
names(list_of_datasets) <- c("Day4_GeneSet_heatmap","Day5_GeneSet_heatmap", "Day6_GeneSet_heatmap")
# Loop through each dataset
for (dataset_name in names(list_of_datasets)) {
  # Extract the dataset from the list
  dataset <- list_of_datasets[[dataset_name]]
  dataset <- dataset[desired_row_order, desired_col_order]
  #Create heatmaps
  p <- pheatmap(dataset,
                scale = "none", # No scaling, use raw cell numbers
                cluster_cols = FALSE,
                cluster_rows = FALSE,
                clustering_method = "complete",
                color = color_palette,
                breaks = breaks,
                annotation_legend = TRUE
  )
  #Save image
  ggsave(filename = paste0(dataset_name, ".tiff"), plot = p, width = 8, height = 4, dpi = 300)
}

























####################Barb's previous work###################################

plot_cells(SP5_cds)
####Subset HE progenitors and Blood progenitors#####
EB_Mesoderm_Mesenchyme_cells_to_keep =
  colData(EB3_assigned_cell_type_cds) %>%
  as.data.frame() %>% 
  filter(!is.na(Cluster),
         Cluster %in% c("7", "22", "28","30", "40", "1","39", "8", "12", "6", "10", "29", "45")) %>% 
  
  pull(Cell) %>%
  as.character() 
EB_M_M_cds = EB3_assigned_cell_type_cds[,EB_Mesoderm_Mesenchyme_cells_to_keep]
plot_cells(EB_M_M_cds)
coldata_SP5_cds = colData(SP5_cds) %>% as_tibble()
rowdata_EB_M_M_cds = rowData(EB_M_M_cds) %>% as_tibble()

#Process subsetted cells (if you want)
EB_M_M_repo_cds <- preprocess_cds(EB_M_M_cds, num_dim = 50) #120 is num_dim for 12-2-20 subset data and most day_cds
EB_M_M_repo_cds <- reduce_dimension(EB_M_M_repo_cds, umap.n_neighbors = 5L) #10L for 12-2-20 subset and most day_cds but lower no. makes closer together. Used 4L for EB_cds
set.seed (17)
EB_M_M_repo_cds <- cluster_cells(EB_M_M_repo_cds, resolution = 8e-4, random_seed = 17) #Resolution - Lower no, fewer clusters - used 3e-3 for 12-2-20 subset data - used 3e-4 for cds
colData(EB_M_M_repo_cds)$umap1 = reducedDim(EB_M_M_repo_cds, type = "UMAP")[,1]
colData(EB_M_M_repo_cds)$umap2 = reducedDim(EB_M_M_repo_cds, type = "UMAP")[,2]
colData(EB_M_M_repo_cds)$Cluster = clusters(EB_M_M_repo_cds, reduction_method = "UMAP")

plot_cells(EB_M_M_repo_cds, color_cells_by = "cluster", group_label_size = 5)
plot_cells(EB_M_M_repo_cds) +
  facet_grid(Activin~BMP)

####Subset Embryo into YS and Emb proper removing E7 since not separated
#Filter sm data both embryo and yolksac
EYS_cells_to_keep =
  colData(Emb1_cds) %>%
  as.data.frame() %>%
  filter(!is.na(part), part %in% c("EP", "YS")) %>%
  pull(Cell) %>%
  as.character()
EYS_cds = Emb1_cds[, EYS_cells_to_keep]
coldata_EYS_cds = colData(EYS_cds) %>% as_tibble()
# Process EYS_cds
EYS_cds <- preprocess_cds(EYS_cds, num_dim =6) #120 is num_dim for 12-2-20 subset data and most day_cds
EYS_cds <- reduce_dimension(EYS_cds, umap.n_neighbors = 4L) #10L for 12-2-20 subset and most day_cds but lower no. makes closer together. Used 4L for EB_cds
set.seed (17)
EYS_cds <- cluster_cells(EYS_cds, resolution = 4e-3, random_seed = 17) #Resolution - Lower no, fewer clusters - used 3e-3 for 12-2-20 subset data - used 3e-4 for cds
plot_cells(EYS_cds)
plot_cells(EYS_cds) +
  facet_wrap(~part, ncol = 2)
ggsave("PLOTS/EYS_clusters_part.tiff", height = 4, width = 6, dpi = 300, bg = "transparent")
#write EYS umap data to cds
colData(EYS_cds)$umap1 = reducedDim(EYS_cds, type = "UMAP")[,1]
colData(EYS_cds)$umap2 = reducedDim(EYS_cds, type = "UMAP")[,2]
colData(EYS_cds)$Cluster = clusters(EYS_cds, reduction_method = "UMAP")

####MARKERS####
####MESODERM MARKERS####

####Lateral Plate Mesoderm lineage markers Prummel KD Mosimann C 2019 Nature Comm Zebrafish####
Lateral_plate_intro <- c("Hand1", "Hand2", "Osr1", "Foxf1", "Prrx1", "Mesp1", "Etv2", "Cxcr4", "Foxh1")
#drivers for early LPM induction - loss will allow later LPM including Hand2 exp contain drl enhancer
Lateral_plate_drivers <- c("Foxh1", "Mixl1", "Eomes")
#do not know where these come from in the Prummel et al paper but from there....
Lateral_plate_additional <- c( "Gata4", "Mef2c", "Tbx5", "Nkx2-5", "mir-133", "Foxo1", "Etv2", "Klf2", "Tal1", "Lmo2",
                               "Erg", "Gata2", "Runx1", "Cxcr4", "Bmp4", "Fgf9", "Apela", "Inhba" )
#MY LATERAL PLATE MARKERS Eomes, Foxh1, Mixl1 (Prummel, 2019) - Gata4 (Rojas 2005) - Cxcr4 (Mcgrath Palis 1999) - Bmp4 (Chandler KJ 2009) - Foxf1 (Funayama N et al)
LPM_BVF <- c("Eomes", "Foxh1", "Mixl1", "Gata4", "Cxcr4", "Bmp4", "Foxf1")


#These are worthless unfortunately: AGM vs Yolk sac PK44 expressed genes: Li et al (Lan ) Cell and Develop Biol 2021 PK44 cells are CD41,CD43, CD45- and CD31, CD201, Kit and CD44+
Yolk_sac_PK44 <- c("Hspa1b","Hspa1a", "Lyve1", "Lum", "Stab2", "GM6654", "Ubb", "Ptbp1", "Rpl5", "Rpl35", "Maf", "Aplnr", "8430408G22Rik", "Malat1", "GM6251", "Col1a2", "Eef2", "Hsp90abl", "Col3a1")
AGM_PK44 <- c("Adarb2", "Adamts13", "Gm12657", "Gm2070", "Mrgpra6", "E030024N20Rik", "Rps15a-ps6", "Gm10012", "Cxcr2", "Gm6402", "Gm15645", "Gm15421", "H3f3a", "Actg1", "Idi1", "Rbm3", "Rn45s", "Ppia", "Gm11517")

#Zhao and Choi LPM to mesenchyme(SM) or hematopoietic mESC Development '19 from gene set score list for LPM
Lateral_plate_Choi <- c("Ahnak", "Ndufa4l2", "Gata6", "Cfc1","Cxcr4", "Gata3", "Pdlim3", "Myl4", "Msx1", "Slc38a4",
                        "Asb4", "Fgf3", "Gpc3", "Hmga2", "Lmo2",  "Dkk1", "Zfp703", "Msx2", "Foxf1", "Hand1",
                        "Bambi", "Ppic", "Amot", "Vstm2b" , "Fgf10", "Bmp4", "Nrp1", "Kdr", "Tbx3", "Pmp22",
                        "Rgs5", "Tgfb2", "Tbx20", "Tgfb3", "Bmp2", "Dusp2")
#####MESODERM#####
#Above LPM culled for YS expressors (Bmp4 and Foxf1) with added Choi markers expressed Nascent/Mixed higher than EXE
LPM <- c("Eomes", "Mixl1", "Gata4", "Cxcr4", "Pdlim3", "Myl4", "Asb4", "Dkk1", "Vstm2b", "Lefty2", "Socs2")
#Other mesoderm
Nascent_mesoderm <- c("Lefty2", "Mesp2", "Dll1", "Dll3", "Frzb", "Pcdh8", "Prss22", "Cdh20", "Slpi", "Pmaip1", "Tmem229a", "Arl4d", "Rimbp2", "Tdgf1")
Somitic_mesoderm <- c("Aldh1a2", "Hoxb3os", "Tbx6", "Rhbg", "Meox1", "Arg", "Pcdh19", "Dll1", "Dll3", "Pcdh8", "Fgf9", "Nkd2", "Aldoc", "Hes7")
Pharyngeal_mesoderm <- c("Mef2c", "Rgs5", "Tcf21", "Gata5", "Meis1", "Lancl3", "Adcyap1r1", "Tnc","A330008L17Rik", "Nxf3", "Fibin", "Pnliprp1", "Tlx1")
Paraxial_mesoderm <- c( "Adrb3", "Col23a1", "Rspo1", "Col26a1", "Prrx1", "Ebf2", "Ptgfr", "Meox1", "Tbx1", "Meox2", "Tbx18", "Tgfbi", "9130410C08Rik", "Itih5", "Ebf3")
EXE_mesoderm <- c("Hand1", "Col9a1", "Radil", "Pde4a", "Fgf15", "Bmp4", "Spin2c", "Dlk1", "Cdx2", "Cdx4")
Mesenchyme <- c("Tdo2", "Acta2", "Tnnt2", "Ahnak", "Tmem108", "Gdf6", "Snai2", "Daam1", "Wisp1", "Stard8", "Lypd6", "Ccdc80", "Sdpr", "Col3a1", "Lum", "Col5a2", "Postn", "Colec11", "Mab21l2")
Cardiomyocytes <- c("Myl7", "Myl4", "Nkx2-5", "Tnnt2", "Mef2c", "Tnnc1", "Sfrp5", "Nexn", "Sh3bgr", "Csrp3", "Gm45123", "Cnn1", "Unc45b", "Ptges3l" )

plot_cells(EYS_cds, genes = Cardiomyocytes_BVF, color_cells_by = "cluster", label_cell_groups = T, show_trajectory_graph = FALSE,
                  label_branch_points = FALSE, label_leaves = FALSE,
                  graph_label_size = 0.75, cell_size = 1) 
       facet_wrap(~part, ncol = 2)
       facet_wrap(~sm + part + stage, ncol = 3)
plot_cells(EYS_cds, genes = c("Socs2"), color_cells_by = "cluster", label_cell_groups = T, show_trajectory_graph = FALSE,
                  label_branch_points = FALSE, label_leaves = FALSE,
                  graph_label_size = 0.75, cell_size = 1) +
         facet_wrap(~part, ncol = 2)
       facet_wrap(~sm + part + stage, ncol = 3)
       
ggsave("PLOTS/EB_M_M_EXE_gt_Mixed.tiff", height = 4, width = 6, dpi = 300, bg = "transparent") 


#########Estimate gene set scores##############
#Function for estimating gene set scores
#### Estimate gene set scores for signature gene lists
estimate_score <- function(cds, markers){
  cds_score = cds[fData(cds)$gene_short_name %in% markers,] 
  aggregate_score = exprs(cds_score)
  aggregate_score = Matrix::t(Matrix::t(aggregate_score) / pData(cds_score)$Size_Factor)
  aggregate_score = Matrix::colSums(aggregate_score)
  pData(cds)$score = log(aggregate_score +1) 
  return(cds)
}

###Specify markers for each gene set



#Generate UMAP plot
plot_cells(cds, color_cells_by = 'score', cell_size = 1, show_trajectory_graph = F) + 
  scale_color_gradient2(low="gray80",mid="gray80",high="red3",midpoint=(((max(pData(cds)$score)-min(pData(cds)$score))/2)+min(pData(cds)$score)), space = "Lab") + simple_theme

#Alternate way to plot gene set scores
mycol <- c("gray80", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
plot_cells(cds, color_cells_by = 'score', cell_size = 0.5, show_trajectory_graph = F)  +
  scale_color_gradientn(colours = mycol)

#### Estimate gene set scores#####
#Brandon code for plotting of gene set scores...

EB_M_M_1_repo1_cds <- estimate_score(EB_M_M_1_cds, markers = Lateral_plate_Choi)
plot_cells(EYS_1_cds, color_cells_by = 'score4', cell_size = 1, show_trajectory_graph = F, alpha = .5) +
  scale_color_gradient2(low="gray80", mid="gray80", high="red",
                        midpoint=(((max(pData(EYS_1_cds)$score4)-min(pData(EYS_1_cds)$score4))/2) +
                                    min(pData(EYS_1_cds)$score4)), space = "Lab") +
  facet_wrap(~part, ncol = 2)
  facet_grid(Activin~BMP)
  

####EB-Embryo combined####
#Function for calculating score
SP6_cds <- SP5_cds       
estimate_score <- function(SP6_cds, markers)
  {SP6_cds_score = SP6_cds[fData(SP6_cds)$gene_short_name %in% markers,] 
  SP6_aggregate_score = exprs(SP6_cds_score)
  SP6_aggregate_score = Matrix::t(Matrix::t(SP6_aggregate_score) / pData(SP6_cds_score)$Size_Factor)
  SP6_aggregate_score = Matrix::colSums(SP6_aggregate_score)
  pData(SP6_cds)$score = log(SP6_aggregate_score +1) 
  return(SP6_cds)}
#Estimating score gets added to column score in designated cds
SP6_cds <- estimate_score(SP6_cds, markers = Lateral_plate_Choi)
####EB #### 
#EB_Mesoderm subset- reclustered
#EB_M_M_repo1_cds has Choi LPM as score1 and Gottgens EXE as score2 and Gottgens mesenchyme as score3 and Gottgens cardiomyocytes as score 4
EB_M_M_1_cds <- EB_M_M_cds       
estimate_LPM_score <- function(EB_M_M_1_cds, markers)
  {EB_M_M_1_cds_LPM_score = EB_M_M_1_cds[fData(EB_M_M_1_cds)$gene_short_name %in% markers,] 
  EB_M_M_1_cds_aggregate_LPM_score = exprs(EB_M_M_1_cds_LPM_score)
  EB_M_M_1_cds_aggregate_LPM_score = Matrix::t(Matrix::t(EB_M_M_1_cds_aggregate_LPM_score) / pData(EB_M_M_1_cds_LPM_score)$Size_Factor)
  EB_M_M_1_cds_aggregate_LPM_score = Matrix::colSums(EB_M_M_1_cds_aggregate_LPM_score)
  pData(EB_M_M_1_cds)$LPM_score = log(EB_M_M_1_cds_aggregate_LPM_score +1) 
  return(EB_M_M_1_cds)} 
EB_M_M_1_cds <- estimate_LPM_score(EB_M_M_1_cds, markers = LPM_BVF1)

estimate_score2 <- function(EYS_1_cds, markers)
  {EYS_1_cds_score2 = EYS_1_cds[fData(EYS_1_cds)$gene_short_name %in% markers,] 
  EYS_1_cds_aggregate_score2 = exprs(EYS_1_cds_score2)
  EYS_1_cds_aggregate_score2 = Matrix::t(Matrix::t(EYS_1_cds_aggregate_score2) / pData(EYS_1_cds_score2)$Size_Factor)
  EYS_1_cds_aggregate_score2 = Matrix::colSums(EYS_1_cds_aggregate_score2)
  pData(EYS_1_cds)$score2 = log(EYS_1_cds_aggregate_score2 +1) 
  return(EYS_1_cds)}
EYS_1_cds <- estimate_score2(EYS_1_cds, markers = EXE_mesoderm_BVF)

estimate_score3 <- function(EYS_1_cds, markers)
  {EYS_1_cds_score3 = EYS_1_cds[fData(EYS_1_cds)$gene_short_name %in% markers,] 
  EYS_1_cds_aggregate_score3 = exprs(EYS_1_cds_score3)
  EYS_1_cds_aggregate_score3 = Matrix::t(Matrix::t(EYS_1_cds_aggregate_score3) / pData(EYS_1_cds_score3)$Size_Factor)
  EYS_1_cds_aggregate_score3 = Matrix::colSums(EYS_1_cds_aggregate_score3)
  pData(EYS_1_cds)$score3 = log(EYS_1_cds_aggregate_score3 +1) 
  return(EYS_1_cds)}
EYS_1_cds <- estimate_score3(EYS_1_cds, markers = Somitic_mesoderm)

estimate_score4 <- function(EYS_1_cds, markers)
  {EYS_1_cds_score4 = EYS_1_cds[fData(EYS_1_cds)$gene_short_name %in% markers,] 
  EYS_1_cds_aggregate_score4 = exprs(EYS_1_cds_score4)
  EYS_1_cds_aggregate_score4 = Matrix::t(Matrix::t(EYS_1_cds_aggregate_score4) / pData(EYS_1_cds_score4)$Size_Factor)
  EYS_1_cds_aggregate_score4 = Matrix::colSums(EYS_1_cds_aggregate_score4)
  pData(EYS_1_cds)$score4 = log(EYS_1_cds_aggregate_score4 +1) 
  return(EYS_1_cds)}
EYS_1_cds <- estimate_score4(EYS_1_cds, markers = Paraxial_mesoderm)

coldata_EB_M_M_1_cds = colData(EB_M_M_1_cds) %>% as_tibble()



####CREATE DF SUMMARIZE gene_set_scores per condition####
#calculate midpoint (MP) same as Brandon code
df_EYS_1 <- coldata_EYS_1_cds %>% 
  dplyr::summarize (Mean1 = mean(score1, na.rm = T), SD1=sd(score1, na.rm=T), MSD1=(Mean1+(2*SD1)), MP1=(((max(score1)-min(score1))/2) + min(score1)),
                    Mean2 = mean(score2, na.rm = T), SD2=sd(score2, na.rm=T), MSD2=(Mean2+(2*SD2)), MP2=(((max(score2)-min(score2))/2) + min(score2)),
                    Mean3 = mean(score3, na.rm = T), SD3=sd(score3, na.rm=T), MSD3=(Mean3+(2*SD3)), MP3=(((max(score3)-min(score3))/2) + min(score3)),
                    Mean4 = mean(score4, na.rm = T), SD4=sd(score4, na.rm=T), MSD4=(Mean4+(2*SD4)), MP4=(((max(score4)-min(score4))/2) + min(score4)))

#determine cluster and condition with highest expression (MP is Midpoint is value MP from above calculation - must be entered)
df_EYS_2 <- coldata_EYS_1_cds %>% 
  group_by (condition, Cluster, .drop = F) %>% 
  dplyr::summarise (number_condition=dplyr::n(), Median1 = median(score1, na.rm=T), Mean1 = mean(score1, na.rm = T), SD1=sd(score1, na.rm=T), Min1 = min(score1), Max1 = max(score1), MSD1=(Mean1+SD1),  
                    MP1 = (1.17), Highest_no1=sum(score1>MP1), Highest_ratio1 = Highest_no1/dplyr::n(), Highest_percent1= Highest_ratio1 * 100,
                    Median2 = median(score2, na.rm=T), Mean2 = mean(score2, na.rm = T), SD2=sd(score2, na.rm=T), Min2 = min(score2), Max2 = max(score2), MSD2=(Mean2+SD2),  
                    MP2 = (.91), Highest_no2=sum(score2>MP2), Highest_ratio2 = Highest_no2/dplyr::n(), Highest_percent2= Highest_ratio2 * 100,
                    Median3 = median(score3, na.rm=T), Mean3 = mean(score3, na.rm = T), SD3=sd(score3, na.rm=T), Min3 = min(score3), Max2 = max(score3), MSD3=(Mean3+SD3),  
                    MP3 = (1.6), Highest_no3=sum(score3>MP3), Highest_ratio3 = Highest_no3/dplyr::n(), Highest_percent3= Highest_ratio3 * 100,
                    Median4 = median(score4, na.rm=T), Mean4 = mean(score4, na.rm = T), SD4=sd(score4, na.rm=T), Min4 = min(score4), Max4 = max(score4), MSD4=(Mean4+SD4),  
                    MP4 = (1.24), Highest_no4=sum(score4>MP4), Highest_ratio4 = Highest_no4/dplyr::n(), Highest_percent4= Highest_ratio4 * 100,) %>% 
  dplyr::ungroup()


####PLOT WITHOUT GRADIENT####
#Subset high for plotting on top is MP from midpoint calculation
df_EYS_LPM_High1 = coldata_EYS_1_cds %>% filter(score1>1.17 & score3<1.6)
df_EYS_EXE_High1 = coldata_EYS_1_cds %>% filter(score2>0.91)
df_EYS_SM_High1 = coldata_EYS_1_cds %>% filter(score3>1.6 )
df_EYS_PM_High1 = coldata_EYS_1_cds %>% filter(score4> 1.24)

df_EB_M_M_1_High2 <- df_LPM_BVF_High1 %>% 
  group_by (condition, Cluster, .drop = F) %>% 
  dplyr::summarise(number_condition=dplyr::n(), Mean1 = mean(score1, na.rm = T), Mean2 = mean(score2, na.rm = T), Mean3 = mean(score3, na.rm = T), Mean4 = mean(score4, na.rm = T))
dplyr::ungroup()

#puts CONDITION 1-16 in my order! Must do for all data frames being plotted
df_EYS_LPM_High1$condition <- factor(df_EYS_LPM_High1$condition, levels=c("1","5","9","13","2", "6","10","14","3","7","11","15","4","8","12","16"))
df_EYS_EXE_High1$condition <- factor(df_EYS_EXE_High1$condition, levels=c("1","5","9","13","2", "6","10","14","3","7","11","15","4","8","12","16"))
coldata_EYS_1_cds$condition <- factor(coldata_EYS_1_cds$condition, levels=c("1","5","9","13","2", "6","10","14","3","7","11","15","4","8","12","16"))
df_EYS_SM_High1$condition <- factor(df_EYS_SM_High1$condition, levels=c("1","5","9","13","2", "6","10","14","3","7","11","15","4","8","12","16"))
df_EYS_PM_High1$condition <- factor(df_EYS_PM_High1$condition, levels=c("1","5","9","13","2", "6","10","14","3","7","11","15","4","8","12","16"))

#Plot embryo data
ggplot () +
  #bottom
  geom_point(data=coldata_EYS_1_cds, aes(umap1, umap2), color = "gray80", size = 0.4, alpha = 0.5) +
  #middle layer 
  geom_point(data=df_EYS_SM_High1, aes(umap1, umap2), color = "gold", size = 0.5, alpha = 1) +
  #middle layer
  geom_point(data = df_EYS_PM_High1, aes (umap1, umap2), color = "mediumpurple", size = 0.2, alpha=1) +
  #middle layer 
  geom_point(data=df_EYS_EXE_High1, aes(umap1, umap2), color = "blue", size = 0.4, alpha = 1) +
  #Top
  geom_point(data=df_EYS_LPM_High1, aes(umap1, umap2), color = "red", size = 0.1, alpha = 1) + 
  #Add facet wrap here
  facet_wrap(~part, nrow = 1) +
  #Add themes 
  theme(plot.title=element_text(hjust=0.5, face='bold', color = 'black', size = 10)) +
  theme(legend.position = "left") +
  my_theme 
  #Save
  ggsave("PLOTS/EXE.tiff", height = 6, width = 6, dpi = 300,bg = "transparent") 
#Plot EB data
  #Top
  geom_point(data=df_Cardio_High1, aes(umap1, umap2), color = "red", size = 0.1, alpha = 1) 
  #middle layer 
  geom_point(data=df_Mesen_High1, aes(umap1, umap2), color = "blue", size = 0.5, alpha = 1) 
  #middle layer 
  geom_point(data=df_EXE_High1, aes(umap1, umap2), color = "gold", size = 0.4, alpha = 1) 
  #middle layer
  geom_point(data = df_LPM_BVF_High1, aes (umap1, umap2), color = "mediumpurple", size = 0.2, alpha=1) 
  #bottom
  geom_point(data=coldata_EB_M_M_repo1_cds, aes(umap1, umap2), color = "gray", size = 0.05, alpha = 0.8) + geom_jitter() +
  #add facet wrap here if want 4x4
  facet_wrap(~condition, nrow = 4) +
  #themes
  theme(plot.title=element_text(hjust=0.5, face='bold', color = 'black', size = 10)) +
  theme(legend.position = "left") +
  my_theme 
  #Save
  ggsave("PLOTS/LPM_BVF_4x4.tiff", height = 4, width = 6, dpi = 300, bg = "transparent") 

facet_wrap(~condition, nrow = 4) 

####PLOT WITH GRADIENT single score with highest GSS on top####
#Reorder data frame so that cells with highest gene set score are on the top.  Reorder from lowest expression to highest.  Has to be a dataframe
Ordered_EB_M_M_1_cds = coldata_EB_M_M_1_cds[order(coldata_EB_M_M_1_cds$LPM_score),]
#Reorder condition or will be plotted 1, 10, 11, 12
Ordered_EB_M_M_1_cds$condition <- factor(Ordered_EB_M_M_1_cds$condition, levels=c("1","5","9","13","2", "6","10","14","3","7","11","15","4","8","12","16"))
coldata_EB_M_M_1_cds$condition <- factor(coldata_EB_M_M_1_cds$condition, levels=c("1","5","9","13","2", "6","10","14","3","7","11","15","4","8","12","16"))
  #top is gradient do not think 0's are plotted???? do not know why you need bottom layer
ggplot (data=Ordered_EB_M_M_1_cds, mapping = aes(x=umap1, y=umap2, color = LPM_score, size = LPM_score)) +
  scale_size_continuous(range = c(0.5, 1)) +
  scale_color_gradient2(low = "gray80", mid = "gray80",  high = "red",
                        midpoint=(((max(Ordered_EB_M_M_1_cds$LPM_score)-min(Ordered_EB_M_M_1_cds$LPM_score))/2) + min(Ordered_EB_M_M_1_cds$LPM_score)), space = "Lab") +
  #bottom is entire data set
  geom_point(data=coldata_EB_M_M_1_cds, aes(umap1, umap2), color = "gray", size = 0.5, alpha = 0.8) + geom_jitter() +
  #add facet wrap here
  facet_wrap(~condition, ncol = 4) +
  #Themes
  theme(plot.title=element_text(hjust=0.5, face='bold', color = 'black', size = 10)) +
  theme(legend.position = "left") +
  my_theme
  #Save
  ggsave("PLOTS/EYS_LPM_Choi_mesoderm.tiff", height = 4, width = 6, dpi = 300, bg = "transparent")
  
####Write dataframe to Excel (rstudio only shows 50 columns)####

  WriteXLS("df_EYS_2", ExcelFileName = "EYS_mesoderm.xlsx")

#tried to reorder cds based on score - cannot figure out how to return reordered coldata to cds
Ordered_coldata =
  colData(EYS_1_cds) %>%
  as.data.frame() 
Ordered_coldata = Ordered_coldata[order(Ordered_coldata$score1),]
rowdata = rowData (EYS_1_cds)
rownames(rowdata) = rownames(EYS_1_cds)
Order1_cds = 
  new_cell_data_set(expression_data = counts(EYS_1_cds),
                    gene_metadata = rowdata,
                    cell_metadata = Ordered_coldata)












