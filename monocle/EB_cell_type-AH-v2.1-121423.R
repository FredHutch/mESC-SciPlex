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

writeLines(capture.output(sessionInfo()), "sessionInfo-mESC-SciPlex_Celltyping.txt")
# Set project directory
projectdir <- "/Users/adamheck/Desktop/mESC-SciPlex"
inputdir <- paste(projectdir, "processed_data/AH_CDS_122023", sep = "/")
plotdir <- paste(projectdir, "results/PLOTS", sep = "/")
outputdir <- paste(projectdir, "results", sep = "/")
setwd(projectdir)

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

##Load in EB3_cds
EB3_cds = readRDS("/Users/adamheck/Desktop/mESC-SciPlex/processed_data/AH_CDS_122023/EB3.RDS")

#Set output directory for plots
setwd(plotdir)
####PLOT CLUSTERS EB_CDS (Figure1B)#####    
plot_cells(EB3_cds,
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
           cell_size = 0.2) +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        text = element_text(size=10)) 

ggsave("FIG1B_EB_cds_cluster.tiff", width = 6, height = 6, dpi = 300, bg = "transparent")

####Plot single gene expression for a given gene
#Set genes in list
gene_list <- c("Hbb-bh1", "Aldh1a2", "Rara","Rarb","Rarg")  # Replace with actual gene names
#gene_list <- c("Cxcr4", "Gfi1","Hlf","Myb","Spi1","Neurl3","Phgdh","Sfrp2","Nupr1","Mycn","Gck","Ift57","Eya2")
for (gene in gene_list) {
  p <- plot_cells(EB3_cds, genes = gene, color_cells_by = "cluster",
             label_cell_groups = FALSE,
             show_trajectory_graph = FALSE,
             label_branch_points = FALSE,
             label_leaves = FALSE,
             graph_label_size = 0.75,
             cell_size = 0.5)
  # Save the plot
  ggsave(filename = paste0(gene, "_exp_full_EBcds.tiff"), plot = p, width = 6, height = 6, dpi = 300, bg = "transparent")
}

#### Generate a heat maps to designate a cell types to a cluster (Supplemental Figure?)#######
#Read in gene short names for gene anotations
genes = 
  read.table(file = "/Users/adamheck/Desktop/mESC-SciPlex/processed_data/Sanjay_files/gene.annotations",
             sep = "\t",
             header = F,
             col.names = c("id", "gene_short_name"))

#Diagnostic markers that describe designated cell type
type_genes <- c(
  #epiblast#  
  "Dnmt3b", "Pou5f1", "Epcam", "Utf1",
  #PGC#
  "Ifitm3", "Dnd1", "Dppa3",
  #primitive streak#
  "Nanog", "Eomes", "Hhex", "Otx2", "Mixl1", "Gsc", "Lhx1","Fgf8", "Pax7",
  #Mesoderm# 
  "Mesp1", "Lefty2", "Mesp2","Osr1", "Cdx2", "Hes7","Hoxb1", "Cdx1", "Gbx2","Tbx1","Meox1", "Tcf15","Aldh1a2", "Dll1", "Tbx6","Isl1","Tcf21", 
  #Exe mesoderm, allantois#
  "Hoxa10", "Hoxa11", "Tbx4", "Bmp4",
  #Mesenchyme#
  "Ahnak", "Pmp22","Krt18","Krt8","Igf2",
  #Cardiomyocytes#                
  "Smarcd3", "Acta2", "Tagln", "Myl7", "Myl4", "Tnnt2",  
  #HE/Endothelial#               
  "Etv2", "Kdr", "Anxa5","Pecam1","Cdh5", 
  #Blood progenitors#               
  "Lmo2","Runx1","Gata1", 
  #Erythroid#                    
  "Hbb-bh1", "Hba-a2", "Hba-a1", "Gypa",
  #Neuroectoderm#
  "Irx3","Hesx1", "Six3","Sox9", "Pax3", "Tfap2a","Foxd3","Dlx2","Sox10","En1", "Pax2",
  #Endoderm#
  "Dkk1", "Krt19", "Amot", "Spink1", "Emb", "Cystm1", "Apoe", "Apoa2", "Ttr")

#subset cds for the gene markers listed above
EB3_cds_type <- EB3_cds [rowData(EB3_cds)$gene_short_name %in% type_genes,]

#Create a df of cell groups (clusters)
EB3_type_group_df <- tibble::tibble(cell=row.names(colData(EB3_cds_type)), EB3_type_group="clusters"(EB3_cds_type)[colnames(EB3_cds_type)])

#generates matrix using aggregate gene expression for cell group [cluster]
EB3_type_ag_mat <- aggregate_gene_expression(EB3_cds_type, gene_group_df=NULL, EB3_type_group_df, norm_method = c("binary"),
                                             scale_agg_values = F)

#replaces gene id with gene short name#
row.names(EB3_type_ag_mat) = genes[match(row.names(EB3_type_ag_mat), genes$id),]$gene_short_name 

### Organize matrix and add annotation###
sub_type <- structure(EB3_type_ag_mat, ".dimnames" = list(c(
  #epiblast#  
  "Dnmt3b", "Pou5f1", "Epcam", "Utf1",
  #PGC#
  "Ifitm3", "Dnd1", "Dppa3",
  #primitive streak#
  "Nanog", "Eomes", "Hhex", "Otx2", "Mixl1", "Gsc", "Lhx1","Fgf8", "Pax7",
  #Mesoderm# 
  "Mesp1", "Lefty2", "Mesp2","Osr1", "Cdx2", "Hes7","Hoxb1", "Cdx1", "Gbx2","Tbx1","Meox1", "Tcf15","Aldh1a2", "Dll1", "Tbx6","Isl1","Tcf21", 
  #Exe mesoderm, allantois#
  "Hoxa10", "Hoxa11", "Tbx4", "Bmp4",
  #Mesenchyme#
  "Ahnak", "Pmp22","Krt18","Krt8","Igf2",
  #Cardiomyocytes#                
  "Smarcd3", "Acta2", "Tagln", "Myl7", "Myl4", "Tnnt2",  
  #HE/Endothelial#               
  "Etv2", "Kdr", "Anxa5","Pecam1","Cdh5", 
  #Blood progenitors#               
  "Lmo2","Runx1","Gata1", 
  #Erythroid#                    
  "Hbb-bh1", "Hba-a2", "Hba-a1", "Gypa",
  #Neuroectoderm#
  "Irx3","Hesx1", "Six3","Sox9", "Pax3", "Tfap2a","Foxd3","Dlx2","Sox10","En1", "Pax2",
  #Endoderm#
  "Dkk1", "Krt19", "Amot", "Spink1", "Emb", "Cystm1", "Apoe", "Apoa2", "Ttr")))

sub_anno_type <- structure(list("cell_type" = c("Epiblast", "Epiblast", "Epiblast", "Epiblast",
                                                "PGC", "PGC", "PGC",
                                                "PS", "PS", "PS", "PS", "PS", "PS", "PS", "PS", "PS",
                                                "Mesoderm", "Mesoderm", "Mesoderm","Mesoderm", "Mesoderm", "Mesoderm","Mesoderm", "Mesoderm", "Mesoderm","Mesoderm", "Mesoderm", "Mesoderm","Mesoderm", "Mesoderm", "Mesoderm","Mesoderm","Mesoderm",
                                                "ExMeso_Alla", "ExMeso_Alla", "ExMeso_Alla", "ExMeso_Alla", 
                                                "Mesenchyme", "Mesenchyme", "Mesenchyme", "Mesenchyme","Mesenchyme",
                                                "Cardiomyocytes", "Cardiomyocytes", "Cardiomyocytes", "Cardiomyocytes", "Cardiomyocytes", "Cardiomyocytes", 
                                                "HE/Endothelial", "HE/Endothelial", "HE/Endothelial", "HE/Endothelial","HE/Endothelial",
                                                "Bloodprog", "Bloodprog", "Bloodprog",
                                                "Erythroid", "Erythroid", "Erythroid", "Erythroid",
                                                "Neuroectoderm", "Neuroectoderm", "Neuroectoderm", "Neuroectoderm", "Neuroectoderm","Neuroectoderm", "Neuroectoderm", "Neuroectoderm", "Neuroectoderm", "Neuroectoderm","Neuroectoderm",
                                                "Endoderm", "Endoderm", "Endoderm","Endoderm", "Endoderm", "Endoderm","Endoderm", "Endoderm", "Endoderm")),
                           .Names = "cell_type",
                           "row.names" = c(
                             #epiblast#  
                             "Dnmt3b", "Pou5f1", "Epcam", "Utf1",
                             #PGC#
                             "Ifitm3", "Dnd1", "Dppa3",
                             #primitive streak#
                             "Nanog", "Eomes", "Hhex", "Otx2", "Mixl1", "Gsc", "Lhx1","Fgf8", "Pax7",
                             #Mesoderm# 
                             "Mesp1", "Lefty2", "Mesp2","Osr1", "Cdx2", "Hes7","Hoxb1", "Cdx1", "Gbx2","Tbx1","Meox1", "Tcf15","Aldh1a2", "Dll1", "Tbx6","Isl1","Tcf21", 
                             #Exe mesoderm, allantois#
                             "Hoxa10", "Hoxa11", "Tbx4", "Bmp4",
                             #Mesenchyme#
                             "Ahnak", "Pmp22","Krt18","Krt8","Igf2",
                             #Cardiomyocytes#                
                             "Smarcd3", "Acta2", "Tagln", "Myl7", "Myl4", "Tnnt2",  
                             #HE/Endothelial#               
                             "Etv2", "Kdr", "Anxa5","Pecam1","Cdh5", 
                             #Blood progenitors#               
                             "Lmo2","Runx1","Gata1", 
                             #Erythroid#                    
                             "Hbb-bh1", "Hba-a2", "Hba-a1", "Gypa",
                             #Neuroectoderm#
                             "Irx3","Hesx1", "Six3","Sox9", "Pax3", "Tfap2a","Foxd3","Dlx2","Sox10","En1", "Pax2",
                             #Endoderm#
                             "Dkk1", "Krt19", "Amot", "Spink1", "Emb", "Cystm1", "Apoe", "Apoa2", "Ttr"),
                           class="data.frame")
#Reorder based on diagnostic marker cell typing annotation
sub_samp_ordered <- sub_type[row.names(sub_anno_type),]
# Create a heatmap of diagnostic markers and clusters
breaks <- seq(0, 3, by = 0.1)
color_palette <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaks) - 1)

pheatmap::pheatmap(sub_samp_ordered, annotation_row = sub_anno_type,color = color_palette,
                   breaks = breaks, clustering_method = "ward.D",
                   annotation_legend = F, cellwidth = 6, cellheight= 6,
                   cluster_rows=F, cluster_cols=T, fontsize = 5, scale = "row",
                   filename = "FIGS2_EB_Heat_map_for_cell_type.tiff",
                   width = 8, height = 11, dpi = 300, bg = "transparent")
dev.off()#Use here b/c heatmap seems to throw off plotter

####Assign cell types using heatmap of clusters and marker expression - Figure1C####
EB3_assigned_cell_type_cds <- EB3_cds
colData(EB3_assigned_cell_type_cds)$assigned_cell_type <- as.character(clusters(EB3_assigned_cell_type_cds))
colData(EB3_assigned_cell_type_cds)$assigned_cell_type = dplyr::recode(colData(EB3_assigned_cell_type_cds)$assigned_cell_type,
        "32" = "Primitive streak", "28" = "Primitive streak", "3" = "Primitive streak", "16" = "Primitive streak", "5" = "Primitive streak", "12" = "Primitive streak",
        "26" = "Epiblast", "33" = "Epiblast", "17" = "Epiblast",
        "20" = "Exe_mesoderm/Allantois", "27" = "Exe_mesoderm/Allantois", "9" = "Exe_mesoderm/Allantois",
        "8" = "Mesoderm", "38" = "Mesoderm", "15" = "Mesoderm", "34" = "Mesoderm", "37" = "Mesoderm",
        "42" = "HE/Endothelial", "7" = "HE/Endothelial", "40" = "HE/Endothelial", "39" = "HE/Endothelial",
        "43" = "HE/Endothelial", "11" = "HE/Endothelial", "2" = "HE/Endothelial", "1" = "HE/Endothelial",
        "14" = "Blood progenitors", "18" = "Blood progenitors", "25" = "Blood progenitors", "10" = "Blood progenitors",
        "23" = "Erythroid", "24" = "Erythroid", "13" = "Erythroid",
        "22" = "Undifferentiated ESC", "4" = "Undifferentiated ESC",
        "36" = "Endoderm", "47" = "Endoderm", "29" = "Endoderm",
        "30" = "Neuroectoderm", "41" = "Neuroectoderm",
        "31" = "Mesenchyme", "19" = "Mesenchyme", "6" = "Mesenchyme",
        "44" = "PGC",
        "35" = "Cardiomyocytes",
        "21" = "Unknown", "45" = "Unknown", "46" = "Unknown")

#plot cells for celltype, full dataset
plot_cells(EB3_assigned_cell_type_cds, group_cells_by = "cluster", color_cells_by = "assigned_cell_type", show_trajectory_graph = F, label_cell_groups = F,
           group_label_size = 4.5, label_groups_by_cluster =F, label_branch_points = F, label_roots = F, label_leaves = F,
           graph_label_size = 0.1, cell_size = 0.2) +
  
  scale_color_manual(values=c("Pink", "Maroon", "PeachPuff", "Khaki",
                              "Red", "Orange", "Turquoise", "SkyBlue",
                              "SlateBlue","Peru", "Yellow", "Orchid", "LightSteelBlue",
                              "SlateGray")) +
  #facet_wrap(~condition, nrow = 4) +    
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.background = element_rect(color = "transparent", fill = "transparent", linetype=0)) +
  theme(panel.border = element_blank(),
        strip.background = element_blank() ,
        strip.text = element_blank()) 
  

ggsave("FIG1C-EB_cds_assigned_cell_types_FULL2.tiff", width = 6, height = 6, dpi = 300, bg = "transparent")

# Plot cells for cell type by Activin/BMP conditions
# Extract colData
cds_colData <- colData(EB3_assigned_cell_type_cds)
# Reorder the 'condition' factor levels
cds_colData$condition <- factor(cds_colData$condition, levels = as.character(1:16))
# Replace the colData in the cds
colData(EB3_assigned_cell_type_cds) <- cds_colData
# Plot cells
plot_cells(EB3_assigned_cell_type_cds, group_cells_by = "cluster", color_cells_by = "assigned_cell_type", show_trajectory_graph = F, label_cell_groups = F,
           group_label_size = 4.5, label_groups_by_cluster =F, label_branch_points = F, label_roots = F, label_leaves = F,
           graph_label_size = 0.1, cell_size = 0.2) +
  
  scale_color_manual(values=c("Pink", "Maroon", "PeachPuff", "Khaki",
                              "Red", "Orange", "Turquoise", "SkyBlue",
                              "SlateBlue","Peru", "Yellow", "Orchid", "LightSteelBlue",
                              "SlateGray")) +
  facet_wrap(~condition, nrow = 4, dir = "v") +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.background = element_rect(color = "transparent", fill = "transparent", linetype=0)) +
  theme(panel.border = element_blank(),
        strip.background = element_blank() ,
        strip.text = element_blank())

ggsave("FIGS4-EB_cds_assigned_cell_types_ActivinBMP-gradient.tiff", width = 6, height = 6, dpi = 300, bg = "transparent")


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
###Graph cell type for individual days in UMAP space
#For loop to run graph generation
list_of_datasets <- list(EB3_assigned_cell_type_D4_cds,EB3_assigned_cell_type_D5_cds,EB3_assigned_cell_type_D6_cds)
names(list_of_datasets) <- c("Day4_EB_cds_assigned_cell_types","Day5_EB_cds_assigned_cell_types","Day6_EB_cds_assigned_cell_types")
# Loop through each dataset
for (dataset_name in names(list_of_datasets)) {
  # Extract the dataset from the list
  dataset <- list_of_datasets[[dataset_name]]
  #plot dataset in UMAP space
  p <- plot_cells(dataset, group_cells_by = "cluster", color_cells_by = "assigned_cell_type", show_trajectory_graph = F, label_cell_groups = F,
                  group_label_size = 4.5, label_groups_by_cluster =F, label_branch_points = F, label_roots = F, label_leaves = F,
                  graph_label_size = 0.1, cell_size = 0.3) +
    
    scale_color_manual(values=c("Pink", "Maroon", "PeachPuff", "Khaki",
                                "Red", "Orange", "Turquoise", "SkyBlue",
                                "SlateBlue","Peru", "Yellow", "Orchid", "LightSteelBlue",
                                "SlateGray")) +
    
    theme(legend.position = "none",
          legend.direction = "vertical",
          legend.title = element_blank(),
          legend.text = element_blank(),
          legend.background = element_rect(color = "transparent", fill = "transparent", linetype=0)) +
    theme(panel.border = element_blank(),
          strip.background = element_blank() ,
          strip.text = element_blank()) +
    no_axes()
  # Save the plot
  ggsave(filename = paste0(dataset_name, ".tiff"), plot = p, width = 6, height = 6, dpi = 300, bg = "transparent")
}

## Create dataframe with genes of interest to use for cell counts and specific gene expresssion
# Need to filter otherwise too big to compute with full cds gene list
Paper_genes <- c("Cxcr4", "Cyp26a1", "Cyp26b1", "Aldh1a2", "Rbp1", "Crabp2", "Crabp1", "Rara", "Rarg", "Rarb",
                 "Cdx2", "Cdx1", "Hoxa9", "Hoxa11", "Arrb1", "Arrb2", "Clec1a", "Clec1b", "P2rx7",
                 "Pdgfra", "Kdr", "Bmp4", "Mixl1", "Eomes", "Lhx1", "Ahnak", "Etv2", "Itga2b",
                 "T", "Runx1", "Cd40", "Cdh5", "Cdh2", "Flt4", "Tek", "Flt1",
                 "Dll4", "Notch1", "Dlk1", "Gja5", "Vwf", "Sox17", "Gata4", "Foxf1", "Runx1",
                 "Cd44", "Procr", "Itgb3", "Itgb7", "Trim47", "Cd38", "Lyve1", "Stab2", "Ace")
Paper_cds <- EB3_assigned_cell_type_cds [rowData(EB3_assigned_cell_type_cds)$gene_short_name %in% Paper_genes,]
Paper_exprs <- SingleCellExperiment::counts(Paper_cds)
Paper_exprs <- reshape2::melt(as.matrix(Paper_exprs)) 
colnames(Paper_exprs) <- c("f_id", "Cell", "expression")
Paper_cds_coldata <- colData(Paper_cds)
Paper_exprs <- merge(Paper_exprs, Paper_cds_coldata, by.x = "Cell", by.y = "row.names")
genes = read.table(file = "/Users/adamheck/Desktop/mESC-SciPlex/processed_data/Sanjay_files/gene.annotations", sep = "\t", header = F, col.names = c("f_id", "gene_short_name"))
Paper_exprs <- merge(Paper_exprs, genes, by = "f_id")
dfPaper = data.frame(Paper_exprs)
dfPaper_1 = dfPaper %>% 
  dplyr::select(Cell,  expression,  gene_short_name, umap1,  umap2, Cluster, condition, day, assigned_cell_type)%>%
  spread(key = gene_short_name, value = expression, fill = 0)

####Generate heatmaps based on cell numbers for cell type or specific cell expression####
####Cell type
dfPaper_celltype <- dfPaper_1 %>% 
  group_by (condition, assigned_cell_type, day, .drop = F) %>% 
  dplyr::summarise (number_condition=dplyr::n()) %>% 
  dplyr::ungroup()

##############FIGURE 1D Making a heatmap for each day looking at cell type and condition#############
df <- dfPaper_celltype
#Split into single day dataframes
day4_df <- df[df$day == '4', ]
day5_df <- df[df$day == '5', ]
day6_df <- df[df$day == '6', ]

####For loop to genereate heatmaps
list_of_datasets <- list(day4_df,day5_df,day6_df)
names(list_of_datasets) <- c("Day4_celltype_heatmap","Day5_celltype_heatmap", "Day6_celltype_heatmap")
# Manually define the desired order for rows and columns
desired_row_order <- c("Epiblast", "Primitive streak","PGC","Endoderm","Neuroectoderm","Mesoderm","HE/Endothelial","Blood progenitors","Erythroid","Mesenchyme","Exe_mesoderm/Allantois","Cardiomyocytes") # Replace with actual row names
desired_col_order <- c("1","5","9","13","2","6","10","14","3","7","11","15","4","8","12","16") # Replace with actual column names
#Manually set the scale of the heatmap
breaks <- seq(0, 1400, by = 1)  # Adjust the 'by' value as needed for finer or coarser color transitions
#Creat the color palette you want
color_palette <- colorRampPalette(c("grey80", "blue", "green","yellow"))(length(breaks) - 1)
# Loop through each dataset
for (dataset_name in names(list_of_datasets)) {
  # Extract the dataset from the list
  dataset <- list_of_datasets[[dataset_name]]
  
  # Remove day column
  dataset <- dataset[, !colnames(dataset) %in% 'day']
  #Reorganize the data frame
  data_matrix <- dcast(dataset, assigned_cell_type ~ condition, value.var = "number_condition")
  #Remove N/A values and make the cell types rownames
  data_matrix[is.na(data_matrix)] <- 0
  rownames(data_matrix) <- data_matrix[,1]
  data_matrix <- data_matrix[,-1]
  #Reorder rows and columns
  data_matrix <- data_matrix[desired_row_order, desired_col_order]
  #Create heatmaps
  p <- pheatmap(data_matrix,
                scale = "none", # No scaling, use raw cell numbers
                cluster_cols = FALSE,
                cluster_rows = FALSE,
                clustering_method = "complete",
                color = color_palette,
                breaks = breaks,
                annotation_legend = TRUE
  )
  #Save image
  ggsave(filename = paste0(dataset_name, ".tiff"), plot = p, width = 8, height = 6, dpi = 300)
}


####### Filter cells based on marker expression#########
#Primitive streak
dfPS_T = dfPaper_1 %>% filter(T>0.1)
#LPM
dfLPM_Pdgfra_Kdr = dfPaper_1 %>% filter(Pdgfra>0.1 & Kdr>0.1)
dfLPM_Pdgfra_Kdr_neg = dfPaper_1 %>% filter(Pdgfra>0.1 & Kdr<0.1)
dfLPM_Pdgfra_neg_Kdr = dfPaper_1 %>% filter(Pdgfra<0.1 & Kdr>0.1)
#Vascular mesoderm
dfHE_Kdr_Etv2_dp = dfPaper_1 %>% filter(Kdr>0.1 & Etv2>0.1)
#Immature HE
dfHE_Kdr_Etv2_Cdh5 = dfPaper_1 %>% filter(Kdr>0.1 & Etv2>0.1 & Cdh5>0.1)
#Mature HE
dfHE_Kdr_Cdh5_Etv2neg = dfPaper_1 %>% filter(Kdr>0.1 & Etv2<0.1 & Cdh5>0.1)
#Mature arterial HE
dfHE_Kdr_Etv2_Cdh5_Dll4 = dfPaper_1 %>% filter(Kdr>0.1 & Etv2<0.1 & Cdh5>0.1 & Dll4>0.1)
####PLOT df COLOR BY DAY of genes expressed
####Plots for Figure 3 and Figure S5
#For loop to run graph generation
list_of_datasets <- list(dfPS_T,dfLPM_Pdgfra_Kdr,dfLPM_Pdgfra_Kdr_neg,dfHE_Kdr_Etv2_dp,dfHE_Kdr_Etv2_Cdh5,dfHE_Kdr_Cdh5_Etv2neg,dfHE_Kdr_Etv2_Cdh5_Dll4)
names(list_of_datasets) <- c("Brachyury+","Pdgfra+_Kdr+", "Pgfra+_Kdr-","Kdr+_Etv2+","Kdr+_Etv2+_Cdh5+","Kdr+_Cdh5+_Etv2-","Kdr+_Cdh5+_Dll4+_Etv2-")
# Loop through each dataset
for (dataset_name in names(list_of_datasets)) {
  # Extract the dataset from the list
  dataset <- list_of_datasets[[dataset_name]]
  
  # Create the plot using ggplot2
  p <- ggplot(data=dataset, mapping = aes(x=umap1, y=umap2, color = day)) +
    scale_size_continuous(range = c(0.2, 0.2)) +
    scale_color_manual(values = c("4" = "#D41159",
                                  "5" = "#0C7BDC",
                                  "6" = "#FFC20A")) + geom_jitter() +
    #bottom
    geom_point(data=dfPaper_1, aes(umap1, umap2), color = "gray", size = 0.05, alpha = 0.8) + geom_jitter() +
    
    #add facet wrap here if want 4x4
    #facet_wrap(~condition, nrow = 4) +  
    theme(plot.title=element_text(hjust=0.5, face='bold', color = 'black', size = 10)) +
    theme(legend.position = "right") +
    theme(panel.border = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "white", color = "white"),
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
          axis.ticks.y = element_blank(),
          legend.key = element_rect(fill = "white")) +
    ggtitle(paste("Plot for", dataset_name)) +
    no_axes()
  
  # Save the plot
  ggsave(filename = paste0(dataset_name, ".tiff"), plot = p, width = 6, height = 6, dpi = 300, bg = "transparent")
}

##############FIGURE XXXX Making a heatmap for each day looking at specific gene exp and condition#############
##KDR,Flk,T,Meso celltypes
df <- dfPaper_1 %>% 
  group_by(condition, day, .drop = F) %>% 
  dplyr::summarise(
    number_condition = dplyr::n(),
    Brach = sum(T > 0.1 & Pdgfra < 0.1 & Kdr < 0.1),
    Pdgfra_Kdr_dp = sum(Pdgfra > 0.1 & Kdr > 0.1),
    Pdgfra_Kdr_neg = sum(Pdgfra > 0.1 & Kdr < 0.1)
  ) %>%
  dplyr::ungroup()

#Split into single day dataframes
day4_df <- df[df$day == '4', ]
day5_df <- df[df$day == '5', ]
day6_df <- df[df$day == '6', ]
####For loop to genereate MESODERM heatmaps
list_of_datasets <- list(day4_df,day5_df,day6_df)
names(list_of_datasets) <- c("Day4_Mesomarker_heatmap","Day5_Mesomarker_heatmap", "Day6_Mesomarker_heatmap")
#Vector for renaming the rows
row_name_mapping <- c("1"="C1","2"="C2","3"="C3", "4"="C4","5"="C5","6"="C6","7"="C7","8"="C8","9"="C9","10"="C10","11"="C11","12"="C12","13"="C13","14"="C14","15"="C15","16"="C16")
# Manually define the desired order for rows and columns
desired_row_order <- c("Brach", "Pdgfra_Kdr_neg","Pdgfra_Kdr_dp") # Replace with actual row names
desired_col_order <- c("C1","C5","C9","C13","C2","C6","C10","C14","C3","C7","C11","C15","C4","C8","C12","C16") # Replace with actual column names
#Manually set the scale of the heatmap
breaks <- seq(0, 177, by = 0.1)  # Adjust the 'by' value as needed for finer or coarser color transitions
#Creat the color palette you want
color_palette <- colorRampPalette(c("grey80", "blue", "green","yellow"))(length(breaks) - 1)
# Loop through each dataset
for (dataset_name in names(list_of_datasets)) {
  # Extract the dataset from the list
  dataset <- list_of_datasets[[dataset_name]]
  # Change condition to the rowname
  dataset <- dataset[, -which(names(dataset) == "condition")]
  # Remove day column
  dataset <- dataset[, !colnames(dataset) %in% 'day']
  # Remove total number column
  dataset <- dataset[, !colnames(dataset) %in% 'number_condition']
  # Rename rownames
  rownames(dataset) <- row_name_mapping[as.character(rownames(dataset))]
  #Transpose the data
  dataset <- t(dataset)
  #Re-order for heatmap
  dataset <- dataset[desired_row_order, desired_col_order]
  #Plot heatmap
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

##Kdr,Etv2,Cdh5 Endothelial celltypes
df <- dfPaper_1 %>% 
  group_by(condition, day, .drop = F) %>% 
  dplyr::summarise(
    number_condition = dplyr::n(),
    VascMeso = sum(Etv2 > 0.1 & Kdr > 0.1 & Cdh5 < 0.1),
    ImmatHE = sum(Etv2 > 0.1 & Kdr > 0.1 & Cdh5 > 0.1),
    MatureHE = sum(Etv2 < 0.1 & Kdr > 0.1 & Cdh5 > 0.1),
    ArterHE = sum(Etv2 < 0.1 & Kdr > 0.1 & Cdh5 > 0.1 & Dll4 > 0.1)
  ) %>%
  dplyr::ungroup()

#Split into single day dataframes
day4_df <- df[df$day == '4', ]
day5_df <- df[df$day == '5', ]
day6_df <- df[df$day == '6', ]
####For loop to genereate MESODERM heatmaps
list_of_datasets <- list(day4_df,day5_df,day6_df)
names(list_of_datasets) <- c("Day4_Endothelial_heatmap","Day5_Endothelial_heatmap", "Day6_Endothelial_heatmap")
#Vector for renaming the rows
row_name_mapping <- c("1"="C1","2"="C2","3"="C3", "4"="C4","5"="C5","6"="C6","7"="C7","8"="C8","9"="C9","10"="C10","11"="C11","12"="C12","13"="C13","14"="C14","15"="C15","16"="C16")
# Manually define the desired order for rows and columns
desired_row_order <- c("VascMeso","ImmatHE", "MatureHE","ArterHE") # Replace with actual row names
desired_col_order <- c("C1","C5","C9","C13","C2","C6","C10","C14","C3","C7","C11","C15","C4","C8","C12","C16") # Replace with actual column names
#Manually set the scale of the heatmap
breaks <- seq(0, 575, by = 0.1)  # Adjust the 'by' value as needed for finer or coarser color transitions
#Creat the color palette you want
color_palette <- colorRampPalette(c("grey80", "blue", "green","yellow"))(length(breaks) - 1)
# Loop through each dataset
for (dataset_name in names(list_of_datasets)) {
  # Extract the dataset from the list
  dataset <- list_of_datasets[[dataset_name]]
  # Change condition to the rowname
  dataset <- dataset[, -which(names(dataset) == "condition")]
  # Remove day column
  dataset <- dataset[, !colnames(dataset) %in% 'day']
  # Remove total number column
  dataset <- dataset[, !colnames(dataset) %in% 'number_condition']
  # Rename rownames
  rownames(dataset) <- row_name_mapping[as.character(rownames(dataset))]
  #Transpose the data
  dataset <- t(dataset)
  #Re-order for heatmap
  dataset <- dataset[desired_row_order, desired_col_order]
  #Plot heatmap
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

####FIGURE 4 HE/Endothelial subset and marker expression
#Subset the mesoderm differentiation cells
endothelial_cells_keep =
  colData(EB3_assigned_cell_type_cds) %>%
  as.data.frame() %>% 
  filter(!is.na(assigned_cell_type),
         assigned_cell_type %in% c("HE/Endothelial")) %>% 
  
  pull(Cell) %>%
  as.character() 
endo_sub_cds = EB3_assigned_cell_type_cds[,endothelial_cells_keep]
#Plot to check the correct cells were subsetted
plot_cells(endo_sub_cds)
#Create expression matrix from endo_sub_cds
endo_cds <- endo_sub_cds [rowData(endo_sub_cds)$gene_short_name %in% Paper_genes,]
#Create expression matrix
endo_exprs <- SingleCellExperiment::counts(endo_cds)
endo_exprs <- reshape2::melt(as.matrix(endo_exprs)) 
colnames(endo_exprs) <- c("f_id", "Cell", "expression")
endo_cds_coldata <- colData(endo_cds)
endo_exprs <- merge(endo_exprs, endo_cds_coldata, by.x = "Cell", by.y = "row.names")
genes = read.table(file = "/Users/adamheck/Desktop/mESC-SciPlex/processed_data/Sanjay_files/gene.annotations", sep = "\t", header = F, col.names = c("f_id", "gene_short_name"))
endo_exprs <- merge(endo_exprs, genes, by = "f_id")
dfendo = data.frame(endo_exprs)
dfendo_1 = dfendo %>% 
  dplyr::select(Cell,  expression,  gene_short_name, umap1,  umap2, Cluster, condition, day, assigned_cell_type)%>%
  spread(key = gene_short_name, value = expression, fill = 0)

#Pull out cells based on specific expression
dfendo_express = dfendo_1 %>% filter(Cdh5>0.1 & Dll4 >0.1)
HE_dataset <- dfendo_express

# Calculate expression patterns
HE_dataset$expression_pattern <- ifelse(HE_dataset$Lyve1 > 0.1, "Lyve1 Positive", "Lyve1 Negative")
#Calculate how many Lyve1 + or - cells there are
sum_of_levels <- sum(expression_counts["Lyve1 Positive"])
print(sum_of_levels)
#214
sum_of_levels <- sum(expression_counts["Lyve1 Negative"])
print(sum_of_levels)
#485

# Define color mappings
color_mapping <- c("Lyve1 Positive" = "darkorange", "Lyve1 Negative" = "darkorchid4")

# Plot the cells in a UMAP space, coloring them based on the expression patterns
ggplot(data=HE_dataset, mapping = aes(x=umap1, y=umap2, color = expression_pattern, size = 0.5)) +
  scale_size_continuous(range = c(0.5, 0.5)) +
  scale_color_manual(values = color_mapping) + 
  geom_jitter() +
  #bottom
  geom_point(data=dfendo, aes(umap1, umap2), color = "gray95", size = 0.1, alpha = 0.8) + geom_jitter() +
  
  #add facet wrap here if want 4x4
  #facet_wrap(~condition, nrow = 4) +  
  theme(plot.title=element_text(hjust=0.5, face='bold', color = 'black', size = 10)) +
  theme(legend.position = "right") +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "white", color = "white"),
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
        axis.ticks.y = element_blank(),
        legend.key = element_rect(fill = "white")) +  # Set the background color of the legend keys to white
  ggtitle(paste("VECad+Dll4+", dataset_name)) +
  no_axes()

ggsave("FIG4B-VECad_Dll4_cells_Lyve1_exp.tiff", width = 6, height = 4, dpi = 300, bg = "white")

###Check expression of HE genes in the context of HE/EC subset UMAP space
gene_list <- c("Cxcr4","Lyve1", "Gfi1","Hlf","Myb","Spi1","Neurl3","Phgdh","Sfrp2","Nupr1","Mycn","Gck","Ift57","Eya2")
for (gene in gene_list) {
  p <- plot_cells(endo_sub_cds, genes = gene, color_cells_by = "cluster",
                  label_cell_groups = FALSE,
                  show_trajectory_graph = FALSE,
                  label_branch_points = FALSE,
                  label_leaves = FALSE,
                  graph_label_size = 0.75,
                  cell_size = 0.5)
  # Save the plot
  ggsave(filename = paste0(gene, "_exp_HE-EC_cds.tiff"), plot = p, width = 6, height = 6, dpi = 300, bg = "transparent")
}

###Check Tessa's gene modules within the context of HE subset
estimate_score <- function(cds, markers){
  cds_score = cds[fData(cds)$gene_short_name %in% markers,] 
  aggregate_score = exprs(cds_score)
  aggregate_score = Matrix::t(Matrix::t(aggregate_score) / pData(cds_score)$Size_Factor)
  aggregate_score = Matrix::colSums(aggregate_score)
  pData(cds)$score = log(aggregate_score +1) 
  return(cds)
}
##Create gene sets
Tessmod9 <- c("Sox17", "Dll4", "Enfb2", "Hey1", "Nrp1", "Nos3","Vwf","Cd44","Cxcr4","Cdkn1c","H19","Mecom","Txnip","Igf2")
Tessmod15 <- c("Runx1","Spi1","Myb","Kit","Gfi1","Npr1","Mycn","Gck","Sfrp","Neurl3","Phgdh","Hlf","Ift57","Lmo1")

list_of_genesets <- list(Tessmod9,Tessmod15)
names(list_of_genesets) <- c("Tessamod9_geneset","Tessmod15_geneset")
#Set colors for plotting gene set scores
mycol <- c("gray80", "gray80", "gray80", "gray80", "gray80", "red1", "red4") # change as needed to highlight different populations
#For loop to generate gene set heat maps
for (geneset_name in names(list_of_genesets)) {
  # Extract the geneset from the list
  geneset <- list_of_genesets[[geneset_name]]
  # Calculate gene set score
  eb_cds <- estimate_score(endo_sub_cds, markers = geneset)
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



####If plots are no longer plotting try####
dev.off()








#### Full cell type classification for generating Fig2A heatmap#######
###### Code Barb initially used ###########
#load RColorBrewer
library(RColorBrewer)
#Read in gene short names for gene anotations
genes = 
  read.table(file = "/Users/adamheck/Desktop/mESC-SciPlex/processed_data/Sanjay_files/gene.annotations",
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
#FIGURE S2A heatmap
pheatmap::pheatmap(sub_samp_ordered, annotation_row = sub_anno_type,
                   color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(25),
                   annotation_legend = F, cellwidth = 6, cellheight= 6,
                   cluster_rows=F, cluster_cols=T, fontsize = 5,
                   filename = "FIGS2_EB_Heat_map_for_cell_type.tiff",
                   width = 8, height = 11, dpi = 300, bg = "transparent")
dev.off()#Use here b/c heatmap seems to throw off plotter


















#############################Extra Barb Code#################################

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
library(WriteXLS)
WriteXLS("dfAEC_Lyve1_Cdh5", ExcelFileName = "dfAEC_Lyve1_Cdh5.xlsx")

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






