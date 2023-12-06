
#following required to install seurat after March 2021 - do first
remove.packages(grep("spatstat", installed.packages(), value = T))
.rs.restartR()
devtools::install_version("spatstat", version = "1.64-1")
suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(Matrix)
  library(tidyverse)
  library(monocle3)
  library(forcats)
  library(pheatmap)
  library(scater)
  library(Seurat)
  library(ggrepel)
  library(spatstat.core)
  setwd("/Users/hadlandlab/Desktop/sciPlex_HSC3_8")
  
  DelayedArray:::set_verbose_block_processing(TRUE)
  options(DelayedArray.block.size=1000e7)})


####SSE (single cell experiment) from online text files####
BiocManager::install(c('rtracklayer'))
library(SingleCellExperiment)
library(scuttle)
library(scran)
library(scater)
library(uwot)
library(rtracklayer)


#####Installing ATLAS Embryo+yolk sac but not separated Gottgens 2019 data####
#There are very clear instructions for loading on github page/website
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MouseGastrulationData")

#Read about Gottgens data
browseVignettes("MouseGastrulationData")

#load Gottgens data
library(MouseGastrulationData)

#tells you sample numbers and Embryonic days - there are 37 samples
head(AtlasSampleMetadata, n = 50)

#generates sce from defined samples
sce7to8 <- EmbryoAtlasData(samples = c(16, 33, 34, 35,
                                      10, 14, 15, 30, 31, 32,
                                      2, 3, 4, 6, 19,20,
                                      23, 26, 27,
                                      8, 12, 13))
sce7.25to7.75 <- EmbryoAtlasData(samples = c(2, 3, 4, 6, 19,20,
                                             23, 26, 27,
                                             8, 9, 12, 13))
sce7.25 <- EmbryoAtlasData(samples = c(23, 26, 27))
sce7.5 <- EmbryoAtlasData(samples = c(2, 3, 4, 6, 19,20))
sce7 <- EmbryoAtlasData(samples = c(2, 3, 4, 6, 19,20,
                                    23, 26, 27,
                                    8, 9, 12, 13,
                                    10, 14, 15, 30, 31, 32))
sce7_7.5_8 <- EmbryoAtlasData(samples = c(16, 33, 34, 35,
                                          10, 14, 15, 30, 31, 32,
                                          2, 3, 4, 6, 19,20))
sce8 <- EmbryoAtlasData(samples = c(16, 33, 34, 35))

#####USED FOR MOST ANALYSES: one sample from E7 and E8 used for most comparisons any more samples and too big for further analysis####
sce78 <- EmbryoAtlasData(samples = c(16, 13, 19, 26, 32))

#Stage used for dissection of Yolk sac Embryo etc
sce8.25 <- EmbryoAtlasData(samples = c(24,25,28))
#View data
sce8.25
counts(sce)[6:9, 1:3]
head(sizeFactors(sce))
head(rowData(sce))
head(colData(sce))

#exclude technical artefacts - usually not done because not necessary
singlets <- which(!(colData(sce8.25)$doublet | colData(sce8.25)$stripped))
#plot gottgens ps way
plot(x = reducedDim(sce8.25, "umap")[singlets, 1],
     y = reducedDim(sce8.25, "umap")[singlets, 2],
  col = EmbryoCelltypeColours[colData(sce78)$celltype[singlets]],
  pch = 19,
  xaxt = "n", yaxt = "n",
  xlab = "DUMB", ylab = "DUMBER") 


#DID NOT DO: get spliced data??? did not do this....
sce8.0 <- EmbryoAtlasData(samples=c(16, 33, 34, 35), get.spliced=TRUE)
names(assays(sce))
#ALSO DID NOT DO: unfilt is all data raw nothing called
unfilt <- EmbryoAtlasData(type="raw", samples=c(1:2))
sapply(unfilt, dim)
sessionInfo()

####CONVERT TO CDS####

#row data has no columns just gene names???
#does not work
sce78 <- as(as.matrix(sce78), "sparseMatrix")
# data.matrix from r does not work no method for coercing this S4 class to a vector
matsce78 <- data.matrix(sce78, rownames.force=NA)
sce78 <- as.matrix(sce78)
#add column with rownames
sce78 <- cbind(sce78, replicate(1, rownames))

coldata_ps78 =
  colData(sce78) %>%
  as.data.frame() 

rowdata =
  rowData(sce78)

rownames(rowdata) =
  rownames(sce78) 

psYEA_cds = 
  new_cell_data_set(expression_data = counts(sce78),
                    gene_metadata = rowdata,
                    cell_metadata = coldata_psYEA) 

# Not sure what next 2 are for
coldata_psYEA_cds =
  colData(psYEA_cds) %>%
  as.data.frame() %>%
  left_join(hash_df, 
            by = "Cell")

coldata_sce8.25 =
  colData(sce8.25) %>%
  as.data.frame() 

#change SYMBOL to gene_short_name
  names(rowData(ps78_cds))[2] <- "gene_short_name"
  
  plot_cells(ps78_cds)

####SUBSET cds#####
ps_mesoderm_cells_to_keep =
    colData(ps_7.25to7.75_cds) %>%
    as.data.frame() %>% 
    filter(!is.na(celltype),
           celltype %in% c("Nascent mesoderm",
                           "Primitive Streak",
                           "Anterior Primitive Streak",
                           "Epiblast",
                           "Haematoendothelial progenitors",
                           "Endothelium",
                           "Erythroid 1",
                           "Erythroid 2",
                           "Blood progenitors 1",
                           "Blood progenitors 2",
                           "Erythroid1",
                           "Erythroid2",
                           "Erythroid3",
                           "Mesenchyme",
                           "Allantois",
                           "Exe mesoderm",
                           "Cardiomyocytes",
                           "Mixed mesoderm",
                           "Pharyngeal mesoderm",
                           "Intermediate mesoderm",
                           "Somitic mesoderm",
                           "Paraxial mesoderm",
                           "PGC")) %>%
    pull(cell) %>%
    as.character() 
ps_nascent_mesoderm_cds = ps_7.25to7.75_cds[,ps_nascent_mesoderm_cells_to_keep] 
#Nascent mesoderm subset
ps_nascent_mesoderm_cells_to_keep =
  colData(ps_7.25to7.75_cds) %>%
  as.data.frame() %>% 
  filter(!is.na(celltype),
         celltype %in% c("Nascent mesoderm",
                         "Mixed mesoderm",
                         "Primitive Streak",
                         "Anterior Primitive Streak",
                         "Epiblast",
                         "Haematoendothelial progenitors"
                         )) %>%
  pull(cell) %>%
  as.character() 
ps_nascent_mesoderm_cds = ps_7.25to7.75_cds[,ps_nascent_mesoderm_cells_to_keep] 



####JOIN EB and GOTTGENS (ps)####

#Change rowname column ENSEMBLE to id in ps_cds
names(rowData(ps78_cds))[1] <- "id"

#Intersect rowdata
rowData(EB_cds)$id_2 = 
  stringr::str_split_fixed(rowData(EB_cds)$id,
                           "\\.",
                           2)[,1]
#intersecting genes should be genes in both cds
intersecting_genes_EB_ps78 = 
  intersect(rowData(ps78_cds)$id %>%
              as.character(),
            rowData(EB_cds)$id_2)

rownames(EB_cds) =rowData(EB_cds)$id_2

EB_ps78_cds_intersected = 
  EB_cds[intersecting_genes_EB_ps78,]

ps78_cds = ps78_cds[intersecting_genes_EB_ps78,]

#following comes back as [1] TRUE - guess they are identical ?????
identical(rownames(EB_ps78_cds_intersected) %>%
            as.character(),
          rownames(ps78_cds) %>%
            as.character())
#makes joint coldata with 3 columns cell(cellname)
#annotation(stage or max_id [dayn_n(concentration)eg day6_16])
#sample(Embryo Atlas or EB Differentiation)
joint_EB_ps78_coldata = 
  rbind(colData(ps78_cds) %>%
          as.data.frame() %>%
          mutate(Cell = 
                   colnames(ps78_cds) %>%
                   as.character()) %>%
          dplyr::select(Cell,
                        annotation = stage) %>%
          mutate(sample = "Embryo Atlas"),
        colData(EB_ps78_cds_intersected) %>%
          as.data.frame() %>%
          dplyr::select(Cell,
                        annotation = max_id) %>%
          mutate(sample = "EB Differentiation"))

rownames(joint_EB_ps78_coldata) = 
  joint_EB_ps78_coldata$Cell

rowdata =
  rowData(ps78_cds) %>% 
  as.data.frame()

joint_EB_ps78_cds =
  new_cell_data_set(expression_data = cbind(counts(ps78_cds), 
                                            counts(EB_ps78_cds_intersected)),
                    cell_metadata = joint_EB_ps78_coldata,
                    gene_metadata = rowdata)
# Run Seurat CCA + Data Integration ---------------------------------------
library(scater)
library(Seurat)
library(spatstat.core)

count_joint_EB_ps78_mat = assay(joint_EB_ps78_cds)
coldata_joint_EB_ps78_cds_df = colData(joint_EB_ps78_cds) %>% as.data.frame()
coldata_joint_EB_ps78_cds_df$Cell = rownames(coldata_joint_EB_ps78_cds_df)

cds_EB_ps78_seurat = 
  CreateSeuratObject(counts = count_joint_EB_ps78_mat,
                     project = "hsc",
                     assay = "RNA",
                     meta.data = coldata_joint_EB_ps78_cds_df)

# split dataset by experiment
cds_EB_ps78.list <- SplitObject(cds_EB_ps78_seurat, split.by = "sample")

# normalize data and find variable genes
for (i in 1:length(cds_EB_ps78.list)) {
  cds_EB_ps78.list[[i]] <- NormalizeData(cds_EB_ps78.list[[i]], verbose = FALSE)
  cds_EB_ps78.list[[i]] <- FindVariableFeatures(cds_EB_ps78.list[[i]], selection.method = "vst", 
                                        nfeatures = 2000, verbose = FALSE)}
#####FIND ANCHORS: This function (FindIntegrationAnchors) requires a lot of memory exhausts vector memory if files are too big 1.5GB seurat worked####
#Takes alot of time about 30 min for 1.5GB
cds_EB_ps78.anchors <- 
  FindIntegrationAnchors(object.list = cds_EB_ps78.list, 
                         dims = 1:100, 
                         anchor.features = 1000,
                         reduction = "cca",
                         verbose = T)


# integrate data 
cds_EB_ps78.integrated <- IntegrateData(cds_EB_ps78.anchors, 
                                dims = 1:50,
                                normalization.method = c("LogNormalize"),
                                k.weight = 10)
                                
cds_EB_ps78.integrated <- ScaleData(cds_EB_ps78.integrated, verbose = FALSE)
cds_EB_ps78.integrated <- RunPCA(cds_EB_ps78.integrated, verbose = FALSE)
cds_EB_ps78.integrated <- RunUMAP(cds_EB_ps78.integrated, dims = 1:50,n.components = 2L)



DimPlot(cds_EB_ps78.integrated, group.by = "sample",split.by = "annotation",ncol = 6) +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("plots/Seurat_EB_ps78_integration_dims50_kweight10_atlas_annotation_50_2L.tiff",
         height = 20,
         width = 20,
         limitsize = F)

DimPlot(cds_EB_ps78.integrated, group.by = "sample") +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("plots/Seurat_EB_ps78_integration_dims50_kweight10_atlas_annotation_color_50_2L.tiff",
         height = 5,
         width = 5,
         limitsize = F)


seurat_umap_embeddings = 
  as.data.frame(Embeddings(cds_EB_ps78.integrated, reduction = "umap"))

coldata_with_seurat_umap_coordinates = 
  left_join(coldata_joint_EB_ps78_cds_df,
            seurat_umap_embeddings %>% 
              rownames_to_column(var = "Cell"),
            by = "Cell")



#### Transfer seurat integrated coordinates to monocle ####-----------------------

joint_cds_EB_ps78_2 = 
  Joint_cds_EBd5_ps78 %>%
  detect_genes() %>%
  estimate_size_factors(method = c("mean-geometric-mean-log-total")) %>%
  preprocess_cds() %>%
  reduce_dimension()


identical(colnames(joint_cds_EB_ps78_2) %>% as.character(),
          coldata_with_seurat_umap_coordinates$Cell %>% as.character())

reducedDim(x = joint_cds_EB_ps78_2,
           type = "UMAP") <-
  matrix(cbind(coldata_with_seurat_umap_coordinates$UMAP_1,
               coldata_with_seurat_umap_coordinates$UMAP_2), 
         ncol=2)

joint_cds_EB_ps78_2 =
  joint_cds_EB_ps78_2 %>%
  cluster_cells(resolution = 3e-4)

saveRDS(object = joint_cds_EB_ps78_2, file = "joint_cds_EB_ps78_2.RDS")
#Writing umap data to cds
colData(joint_cds_EB_ps78_2)$umap1 =
  reducedDim(joint_cds_EB_ps78_2,
             type = "UMAP")[,1]

colData(joint_cds_EB_ps78_2)$umap2 =
  reducedDim(joint_cds_EB_ps78_2,
             type = "UMAP")[,2]


colData(joint_cds_EB_ps78_2)$cluster =
  clusters(joint_cds_EB_ps78_2,
           reduction_method=  "UMAP")

coldata_joint_cds_EB_ps78_2 = 
  colData(joint_cds_EB_ps78_2) %>%
  as.data.frame() 

#### ANNOTATION ####--------------------------------------------------------------

coldata_atlas = 
  read.table(file = "/Users/hadlandlab/Desktop/meta.tab.csv", 
             header = T, 
             sep = ",")

#didnot do the following of Sanjay code because I was unable to add ps info to jointcoldata in next section
coldata_atlas$study = 
  "Pijuan"

coldata_atlas$super_sample = 
  "wt.atlas"


coldata_atlas$Cell = 
  paste(coldata_atlas$study,
        coldata_atlas$super_sample,
        coldata_atlas$sample, 
        coldata_atlas$cell, 
        coldata_atlas$barcode, 
        sep = "_")

#rename cell in coldata_atlas to Cell so pjs data can be added to coldata_joint
colnames(coldata_atlas) [colnames(coldata_atlas) == "cell"] = "Cell"

#adding EB data and pjs data including umap coordinates to coldata_joint
coldata_joint_cds_EB_ps78_3 =
  coldata_joint_cds_EB_ps78_2 %>%
  left_join(coldata_atlas %>%
              dplyr::select(Cell,
                            umap_1_pjs = umapX,
                            umap_2_pjs = umapY,
                            celltype_pjs = celltype,
                            cluster_pjs = cluster,
                            theiler_pjs = theiler),
            by = "cell")

coldata_joint_cds_EB_ps78_2 = 
  coldata_joint_cds_EB_ps78_2 %>%
  left_join(colData(EB_cds) %>%
              as.data.frame() %>%
              dplyr::select(Cell,
                            cluster_eb = Cluster,
                            condition_eb = condition,
                            BMP,
                            Activin,
                            dataset_eb = data_set,
                            day_eb = day),
            by = "Cell")

ggplot() +
  geom_point(data = 
               coldata_joint_cds_EB_ps78_2 %>%
               sample_n(size = nrow(coldata_joint_cds_EB_ps78_2),
                        replace = F) %>%
               mutate(day_alt = ifelse(!is.na(dataset_eb),
                                       day_eb,
                                       "NA")),
             aes(x = umap1,
                 y = umap2,
                 color = day_alt),
             size = 0.55,
             stroke = 0) +
  monocle3:::monocle_theme_opts() +
  clear_theme() +
  no_axes() +
  scale_color_manual(values = c("4" = "#ffeda0",
                                "5" = "#feb24c",
                                "6" = "#f03b20",
                                "NA" = "grey80")) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank()) +
  ggsave("plots/joint_EB_day_umap.png",
         height = 5,
         width = 5,
         dpi = 600,
         bg = "transparent")

#creates a data frame for cell type with median umap coordinates - allows addition of cell gottgens cell type to plot
coldata_joint_cds_EB_ps78_2 = colData(joint_cds_EB_ps78_2) %>% as.data.frame()
label_position_df = 
  coldata_joint_cds_EB_ps78_2 %>%
  group_by(celltype_pjs) %>%
  summarise(med_umap1 = median(umap1),
            med_umap2 = median(umap2)) %>%
  drop_na()

#do not know what this is for: get error in filter cannot coerce type 'closure' to vector of type 'character'
coldata_joint_cds_EB_ps78_1 %>%
  mutate()
filter(grepl(x = annotation,
             pattern = "E"))

#plots cell type with Gottgens data (dataset_eb is filtered)
####INCREASE MAX OVeRLAPS####
options(ggrepel.max.overlaps = Inf)

theme(legend.position = "top",
      strip.background = element_blank(),
      strip.text = element_blank(),
      legend.title = "none",
      legend.key = element_rect(fill = NULL, size = 10)) +
  ggsave("Plots/joint_atlas_celltype_psonly_umap.tiff",
         height = 8,
         width = 8,
         dpi = 600,
         bg = "transparent")

#plot cell types with all data including eb (gray on top) shows labels genes left panel not overlaid
plot_cells(Joint_cds_EB_ps78_2, 
           color_cells_by = "cluster",
           label_cell_groups = TRUE,
           label_branch_points = FALSE,
           show_trajectory_graph = FALSE,
           label_leaves = FALSE,
           graph_label_size = 0.75,
           cell_size = 0.4) 
ggplot() +
  geom_point(data = 
               coldata_joint_cds_EB_ps78_2,
             
             aes(x = umap1,
                 y = umap2,
                 color = celltype_pjs),
             size = 0.55,
             stroke = 0) +
  theme(strip.text = element_text(size = 6, color = "black", margin = margin(0.1,0.1,0.1,0.1, "cm")),
        strip.background = element_rect(colour = "blue", fill = "transparent"),
        panel.background = element_rect(colour = "blue", fill = "transparent"),
        panel.grid.major.x = element_line(color = "grey50"),
        panel.grid.major.y = element_line(color = "grey50"),
        panel.grid.minor.x = element_line(colour = "grey50"),
        panel.grid.minor.y = element_line(colour = "grey50"),
        panel.ontop = T) + 
  geom_label_repel(data = label_position_df,
                   aes(x = med_umap1,
                       y = med_umap2,
                       label = celltype_pjs),
                   size = 1.5,
                   label.size = .15,
                   alpha = 0.65,
                   seed = 42,
                   box.padding = 0.35,
                   point.padding = 0.25) +
  geom_label_repel(data = label_position_df,
                   aes(x = med_umap1,
                       y = med_umap2,
                       label = celltype_pjs),
                   size = 1.5,
                   label.size = .15,
                   fill = NA,
                   seed = 42,
                   box.padding = 0.35,
                   point.padding = 0.25) +
  
theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank()) +
  ggsave("Plots/joint_atlas_celltype_umap.tiff",
         height = 10,
         width = 10,
         dpi = 600,
         bg="transparent")



Lateral_plate_markers <- c("Mixl1",
                           "Eomes",
                           "Foxh1",
                           "Cxcr4",
                           "Gata4",
                           "Bmp4")
####PLOT CELLS JOINT CDS#####
#clusters with celltype
plot_cells(Joint_cds_EB_ps78_2, 
           color_cells_by = "cluster",
           label_cell_groups = T,
           label_branch_points = FALSE,
           show_trajectory_graph = FALSE,
           label_leaves = FALSE,
           graph_label_size = 0.75,
           cell_size = 0.3) +
geom_label_repel(data = label_position_df,
                   aes(x = med_umap1,
                       y = med_umap2,
                       label = celltype_pjs),
                   size = 1.5,
                   label.size = .1,
                   fill = alpha(c("white"), 0.5),
                   
                   seed = 42,
                   box.padding = 0.2,
                   point.padding = 0.2) +
theme(strip.text = element_text(size = 6, color = "black", margin = margin(0.1,0.1,0.1,0.1, "cm")),
      strip.background = element_rect(colour = "blue", fill = "transparent"),
      panel.background = element_rect(colour = "blue", fill = "transparent"),
      panel.grid.major.x = element_line(color = "grey50"),
      panel.grid.major.y = element_line(color = "grey50"),
      panel.grid.minor.x = element_line(colour = "grey50"),
      panel.grid.minor.y = element_line(colour = "grey50"),
      panel.ontop = T) + 

  no_axes()+
  
ggsave("PLOTS/EB_ps78_2_cds.tiff",
       height = 8,
       width = 8,
       dpi = 600,
       bg = "transparent")

#plot 6 genes
plot_cells(joint_cds_EB_ps78_2, genes = c("Cdh5", "Cxcr4", "Kdr", "Pdgfra", "Mixl1", "Gata4"),
           color_cells_by = "cluster",
           label_cell_groups = T,
           label_branch_points = FALSE,
           show_trajectory_graph = FALSE,
           label_leaves = FALSE,
           graph_label_size = 0.75,
           cell_size = 0.3) +
  theme(strip.text = element_text(size = 6, color = "black", margin = margin(0.1,0.1,0.1,0.1, "cm")),
        strip.background = element_rect(colour = "blue", fill = "transparent"),
        panel.background = element_rect(colour = "blue", fill = "transparent"),
        panel.grid.major.x = element_line(color = "grey50"),
        panel.grid.major.y = element_line(color = "grey50"),
        panel.grid.minor.x = element_line(colour = "grey50"),
        panel.grid.minor.y = element_line(colour = "grey50"),
        panel.ontop = T) + 
  no_axes()+

ggsave("PLOTS/joint_cds_EB_ps78_1_cds_LPMgenes.tiff",
       height = 8,
       width = 8,
       dpi = 600,
       bg = "transparent")

#plot 2x2 sample
plot_cells(joint_cds_EB_ps78_2, genes=c("Ahnak"),
           color_cells_by = "cluster",
           label_cell_groups = F,
           label_branch_points = FALSE,
           show_trajectory_graph = FALSE,
           label_leaves = FALSE,
           graph_label_size = 0.75,
           cell_size = 0.3) +
geom_label_repel(data = label_position_df,
                   aes(x = med_umap1,
                       y = med_umap2,
                       label = celltype_pjs),
                   size = 1.5,
                   label.size = .1,
                   fill = alpha(c("white"), 0.3),
                   
                   seed = 42,
                   box.padding = 0.2,
                   point.padding = 0.2) +
  
theme(strip.background = element_rect(colour = "blue", fill = "transparent"),
      panel.background = element_rect(colour = "blue", fill = "transparent"),
      panel.grid.major.x = element_line(color = "grey50"),
      panel.grid.major.y = element_line(color = "grey50"),
      panel.grid.minor.x = element_line(colour = "grey50"),
      panel.grid.minor.y = element_line(colour = "grey50"),
        panel.ontop = T) +
  
facet_grid ("sample") +


ggsave("PLOTS/Joint_EB_ps78_cds_2_Ahnak_2x2.tiff",
         height = 8,
         width = 8,
         dpi = 600,
         bg = "transparent")

# plot 5x4+1 concentration
#must run the following to get the correct order in facet_wrap - unfortunately gradient is reversed 1-4 accross top instead of left side
joint_cds_EB_ps78_2$annotation <- factor(joint_cds_EB_ps78_2$annotation,
                                            levels=c("day5_1","day5_5","day5_9","day5_13",
                                                     "day5_2", "day5_6","day5_10","day5_14",
                                                     "day5_3","day5_7","day5_11","day5_15",
                                                     "day5_4","day5_8","day5_12","day5_16",
                                                     "E7.0", "E7.25", "E7.5", "E7.75", "E8"))
plot_cells(EB_cds, genes=c("Rarg"),
           color_cells_by = "cluster",
           label_cell_groups = T,
           label_branch_points = FALSE,
           show_trajectory_graph = FALSE,
           label_leaves = FALSE,
           graph_label_size = 0.25,
           cell_size = 0.5) +

facet_wrap(~condition, ncol = 4) +
 theme(strip.text = element_text(size = 6, color = "black", margin = margin(0.1,0.1,0.1,0.1, "cm")),
      strip.background = element_rect(colour = "blue", fill = "transparent"),
        panel.background = element_rect(colour = "blue", fill = "transparent"),
        panel.grid.major.x = element_line(color = "grey50"),
        panel.grid.major.y = element_line(color = "grey50"),
        panel.grid.minor.x = element_line(colour = "grey50"),
        panel.grid.minor.y = element_line(colour = "grey50"),
        panel.ontop = T) 


ggsave("PLOTS/joint_cds_EB_ps78_2_Hoxd1_Exe_mesoderm_4x4.tiff",
       height = 8,
       width = 8,
       dpi = 600,
       bg = "transparent")

####plot cells#####
plot_cells(ps78_1_cds, genes=c("Cxcr4","Dll4","Mecom", "Gata3", "H19", "Cdkn1c", "Txnip", "Gja5", "Gja4"),
           color_cells_by = "cluster",
           label_cell_groups = TRUE,
           label_branch_points = FALSE,
           show_trajectory_graph = FALSE,
           label_leaves = FALSE,
           graph_label_size = 0.75,
           cell_size = 0.4) 

ggsave("PLOTS/ps78_1_HEHSC_diff_genes.tiff",
       height = 8,
       width = 8,
       dpi = 600,
       bg = "transparent")

plot_cells(EB_cds, genes=c("Cxcr4","Cdh4","Cdh8", "Cdh11", "Lrp1", "Mixl1"),
           color_cells_by = "cluster",
           label_cell_groups = TRUE,
           label_branch_points = FALSE,
           show_trajectory_graph = FALSE,
           label_leaves = FALSE,
           graph_label_size = 0.75,
           cell_size = 0.4) 

ggsave("PLOTS/EB_cds_sort_genes.tiff",
       height = 8,
       width = 8,
       dpi = 600,
       bg = "transparent")



#plot genes in joint_cds
 p1 = plot_cells(joint_cds_EB_ps78_1,
             genes = c("Mixl1",
                       "Cxcr4",
                       "Gata4",
                       "Eomes",
                       "Tdo2",
                       "Bmp4",
                       "Acta2",
                       "Foxf1",
                       "Irx3"),
             cell_size = 0.25,
             cell_stroke = 0,
             label_cell_groups = F,
             label_branch_points = F,
             label_groups_by_cluster = F)+
  
p1$facet$params$nrow = 3 +
p1+ ggsave("plots/joint_umap_mesoderm_genes.tiff",
         height = 8,
         width = 8,
         dpi = 600,
         bg = "transparent")
 
 #Plot ps without EB
 ggplot() +
   geom_point(data = 
                coldata_joint_cds_EB_ps78_1 %>%
                filter(is.na(dataset_eb)),
              aes(x = umap1,
                  y = umap2,
                  color = celltype_pjs),
              size = .5,
              stroke = 0) +
   scale_color_manual (values=c("pink", "gray", "#415A65", "red",
                                "#7FDFE1", "#B3F80B", "#AAAFCA", "#0BF6F8",
                                "#396F4E", "#BBDD7F", "#45BCF3", "#4A6B4B",
                                "pink", "gray", "#415A65", "red",
                                "#7FDFE1", "#B3F80B", "#AAAFCA", "#0BF6F8",
                                "#396F4E", "#BBDD7F", "#45BCF3", "#4A6B4B",
                                "pink", "gray", "#415A65", "red",
                                "#7FDFE1", "#B3F80B", "#AAAFCA", "#0BF6F8",
                                "#396F4E", "#BBDD7F", "#45BCF3", "#4A6B4B", "#45BCF3", "#4A6B4B", "#4A6B4B")) +
   guides(color = guide_legend(override.aes = list(size = 3))) +
   geom_label_repel(data = label_position_df,
                    aes(x = med_umap1,
                        y = med_umap2,
                        label = celltype_pjs),
                    size = 1.5,
                    label.size = .1,
                    fill = alpha(c("white"), 0.5),
                    
                    seed = 42,
                    box.padding = 0.2,
                    point.padding = 0.2) +
   
   monocle3:::monocle_theme_opts() +
   clear_theme() +
   no_axes() +
   theme(legend.position = "top",
         strip.background = element_blank(),
         strip.text = element_blank())+
   theme(panel.border = element_rect(fill = NA, color = "black"),
         legend.text=element_text(size=rel(.6)))
   
   ggsave("plots/joint_atlas_celltype_umap1.tiff")
   
#Plot ps with EB

plot_cells(joint_cds_EB_ps78_1,
           color_cells_by = "cluster",
           cell_size = 0.25,
           cell_stroke = 0,
           label_cell_groups = T,
           label_branch_points = F,
           label_groups_by_cluster = F) +
  geom_point(data = 
               coldata_joint_cds_EB_ps78_1 %>%
               filter(is.na(cluster_pjs)))
            
     guides(color = guide_legend(override.aes = list(size = 3))) +
     geom_label_repel(data = label_position_df,
                      aes(x = med_umap1,
                          y = med_umap2,
                          label = celltype_pjs),
                      size = 1.5,
                      label.size = .1,
                      fill = alpha(c("white"), 0.5),
                      
                      seed = 42,
                      box.padding = 0.2,
                      point.padding = 0.2) +
     
     monocle3:::monocle_theme_opts() +
     clear_theme() +
     no_axes() +
     theme(legend.position = "top",
           strip.background = element_blank(),
           strip.text = element_blank())+
     theme(panel.border = element_rect(fill = NA, color = "black"),
           legend.text=element_text(size=rel(.6)))
ggsave("plots/joint_atlas_celltype_umap_with EB.tiff") 

#saved for above for some reason removes 5000 cells and cannot overlay EB
scale_color_manual (values=c("pink", "gray", "#415A65", "red",
                             "#7FDFE1", "#B3F80B", "#AAAFCA", "#0BF6F8",
                             "#396F4E", "#BBDD7F", "#45BCF3", "#4A6B4B",
                             "pink", "gray", "#415A65", "red",
                             "#7FDFE1", "#B3F80B", "#AAAFCA", "#0BF6F8",
                             "#396F4E", "#BBDD7F", "#45BCF3", "#4A6B4B",
                             "pink", "gray", "#415A65", "red",
                             "#7FDFE1", "#B3F80B", "#AAAFCA", "#0BF6F8",
                             "#396F4E", "#BBDD7F", "#45BCF3", "#4A6B4B" ))
  
####Coldata for cds####
coldata_EB_PS_mesoderm_day45_cds = colData(EB_PS_mesoderm_day45_cds) %>% as_tibble()
####Save the current working directory####
save.image(file = 'HSC3_10_LPM_df.RData')

####My attempt to load yolk sac vs embryo data from PS-Gottgens github DIDNT WORK#####
#load raw data from pijuan-sala website originally
mat1 <- read.delim(("/Users/hadlandlab/Desktop/sciPlex_HSC3_8/E-MTAB-6970/counts_raw.txt"), header = TRUE, row.names=1, check.names=FALSE)
#load lognorm data - used this 11-9-21
mat <- read.delim(("/Users/hadlandlab/Desktop/sciPlex_HSC3_8/E-MTAB-6970/counts_logNorm.txt"), header = TRUE, row.names=1, check.names=FALSE)
#only endogenous genes ????? (didnot remove any genes)
spike.mat1 <- mat[grepl("^ERCC-", rownames(mat1)),]
mat1 <- mat[grepl("^ENSMUSG", rownames(mat1)),]
#splitting off the gene length column ???? (also didnot remove any genes - maybe already done)
gene.length <- mat1[,1]
mat <- as.matrix(mat1[,-1])
dim(mat)

#Construct SCE object
sce <- SingleCellExperiment(assays = list(counts = mat))
sce
sce1 <- SingleCellExperiment(assays = list(counts = mat1))
sce1

# Compute a log-transformed normalized expression matrix of raw data 
# Adding logcounts to the assay slot?????
sce1 <- scuttle::logNormCounts(sce1)
# Following works assess count data ?????
assay(sce, "counts")

#generating coldata matrix
#inputting data downloaded from pijuan-sala website 2019 paper
coldata <- read.delim("/Users/hadlandlab/Desktop/sciPlex_HSC3_8/E-MTAB-6970/E-MTAB-6970.txt", check.names=FALSE)
#removing unneeded columns from coldata and keeping what I hope is important
coldata <- data.frame(
  row.names=coldata[,"Factor Value[single cell identifier]"],
  embryo.part=coldata[,"Factor Value[organism part]"],
  assay.name=coldata[,"Assay Name"],
  technical.replicates=coldata[,"Comment[technical replicate group]"]
)
coldata

####CONVERT TO CDS####
#error for YEA needs to be a  matrix?
#row data has no columns just gene names???
#does not work
sce <- as(as.matrix(sce), "sparseMatrix")
# data.matrix from r does not work no method for coercing this S4 class to a vector
matsce <- data.matrix(sce, rownames.force=NA)
sce <- as.matrix(sce)
#add column with rownames
sce <- cbind(sce, replicate(1, rownames))

coldata_psYEA =
  colData(sce) %>%
  as.data.frame() 

rowdata =
  rowData(sce)

rownames(rowdata) =
  rownames(sce) 

psYEA_cds = 
  new_cell_data_set(expression_data = counts(sce),
                    gene_metadata = rowdata,
                    cell_metadata = coldata_psYEA) 
