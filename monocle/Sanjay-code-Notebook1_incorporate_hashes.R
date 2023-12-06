
suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(Matrix)
  library(tidyverse)
  library(monocle3)
  setwd("/net/trapnell/vol1/home/sanjays/projects/hsc/sciChem_HSC_3/analysis")
  source("chiSq_test_functions.R")
  DelayedArray:::set_verbose_block_processing(TRUE)
  options(DelayedArray.block.size=1000e7)
})

# Plotting function -------------------------------------------------------

clear_theme = 
  function () {
    theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")
    )}

no_axes =
  function () {
    theme(axis.line.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank()
    )}
      
white_theme =   
  function () {
  theme(text = element_text(color = "white"), 
        axis.line.x = element_line(color = "white"),
        axis.line.y = element_line(color = "white"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(color = "white"),
        strip.background = element_blank(),
        strip.text.x = element_text(color = "white")
  )}
  

cds = 
  readRDS("/net/trapnell/vol1/home/sanjays/projects/hsc/sciChem_HSC_3/final-output/cds.RDS")

hashTable = 
  read.table(file = "/net/trapnell/vol1/home/sanjays/projects/hsc/sciChem_HSC_3/hashRDS/hashTable.out",
             sep = "\t",
             header = F,
             col.names = c("sample", "Cell", "oligo", "axis", "count"))



# Read in the number of RNA UMI molecules per cell 
rna_umis = 
  read.table("/net/trapnell/vol1/home/sanjays/projects/hsc/sciChem_HSC_3/UMIs.per.cell.barcode_2",
             col.names = c("sample","Cell","n.umi")) %>%
  dplyr::select(-"sample")


# Tabulate per cell information from hash table
hashTable = 
  hashTable %>%
  dplyr::group_by(axis, Cell) %>%
  dplyr::mutate(total_within_axis = sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Cell) %>%
  dplyr::mutate(total_hash = sum(count)) %>%
  dplyr::ungroup() %>%
  left_join(rna_umis,
            by = "Cell") %>%
  mutate(in_cds = Cell %in% colnames(cds))


# Make RNA knee plot to set background UMI treshold for estimating background hashes

# The RNA cutoff was set at 300 
# The background count treshold was set at 10

breaks_for_plot = c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000, 500000)

rna_umis %>%
  dplyr::select(n.umi) %>%
  arrange(desc(n.umi)) %>%
  mutate(rank = dplyr::row_number()) %>%
  ggplot() +
  geom_line(aes(x = rank,
                y = n.umi)) +
  scale_x_log10(breaks = breaks_for_plot) +
  scale_y_log10(breaks = breaks_for_plot) +
  geom_hline(yintercept = 300, 
             color = "red", 
             size = 0.5) +
  geom_hline(yintercept = 15, 
             color = "dodgerblue", 
             size = 0.5)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 8))+
  clear_theme() +
  white_theme() +
  theme(panel.background = 
          element_blank(), 
        panel.border = element_rect(fill = NA, 
                                    colour = "white"))+
  xlab("Number of Cells") +
  ylab("RNA UMIs") +
  ggsave(filename = "plots/Notebook1_RNA_kneeplot.pdf",
         height = 2,
         width = 3) 


# Correct the  hash table by removing background estimated from debris --------

# The hashTable is in a "long" form (as opposed to wide) where every row is a Cell/hashOligo
# combination with the number of UMIs recovered for that hash Oligo. Additionally, based on the
# design of the experiment hashes can denote labeling of different things.

corrected_hash_table = 
  hashTable %>%
  # perfor background correction separately for each axis
  group_by(axis) %>% 
  dplyr::rename(Oligo = oligo,
                Count = count) %>%
  nest() %>% 
  mutate(corrected_df = purrr::map(data, .f = function(hash_table_subset){
    hash_table_subset$Oligo = 
      factor(hash_table_subset$Oligo, 
             levels = unique(hash_table_subset$Oligo))
    
    # Background Cells fall below RNA treshold
    background_cells = 
      hash_table_subset %>%
      filter(n.umi < 10) %>%
      pull(Cell) %>%
      as.character()
    
    # Test Cells that are in the CDS object
    test_cells =
      hash_table_subset %>%
      filter(in_cds) %>%
      pull(Cell) %>%
      as.character()
    
    background_sample_hashes = 
      hash_table_subset %>% 
      filter(Cell %in% background_cells) 
    
    test_cell_hashes = 
      hash_table_subset %>% 
      filter(Cell %in% test_cells) 
    
    hash_df = assign_hash_labels_return_all(test_cell_hashes, 
                                            background_sample_hashes, 
                                            downsample_rate=1)
    
    left_join(hash_table_subset,
              hash_df,
              by = c("Cell","Oligo")) 
  })) %>%
  dplyr::select(-data) %>%
  unnest()


# Recover the experimental information from hashing data -------------

# Make a matrix representation of the data
hash_dataframe = 
  corrected_hash_table %>%
  ungroup() %>%
  filter(axis == 1,
         in_cds) %>%
  dplyr::select(Cell, 
                Oligo, 
                adjusted_count)%>%
  spread(key = Oligo, 
         value = adjusted_count,
         fill = 0)

hash_dataframe_cell = hash_dataframe$Cell
hash_dataframe = hash_dataframe[,2:ncol(hash_dataframe)]

# Get the Maximum slide oligo per Cell
max_oligo = 
  colnames(hash_dataframe)[apply(X = hash_dataframe, 
                              MARGIN = 1, 
                              FUN = (which.max))]

max_slide_val = 
  apply(X = hash_dataframe, 
        MARGIN = 1, 
        FUN = function(x){
          maximum_col = which.max(x)
          
          x[maximum_col]
        })

ratio_top_2 = 
  apply(X = hash_dataframe, 
        MARGIN = 1, 
        FUN = function(x){
          maximum_col = which.max(x)
          max_val = x[maximum_col]
          x[maximum_col] = 0
          seoncd_max = which.max(x)
          second_max_val = x[seoncd_max]
          max_val/second_max_val
        })

hash_total = 
  rowSums(hash_dataframe)

hash_df = data.frame(
  Cell = hash_dataframe_cell,
  max_val = max_slide_val,
  max_id = max_oligo,
  ratio_top_2 = ratio_top_2,
  hash_total)


# Incorporate Hash Data into the CDS --------------------------------------

coldata =
  colData(cds) %>%
  as.data.frame() %>%
  left_join(hash_df, 
            by = "Cell")

rownames(coldata) = 
  coldata$Cell

rowdata =
  rowData(cds)

rownames(rowdata) =
  rownames(cds)

cds = 
  new_cell_data_set(expression_data = counts(cds),
                    gene_metadata = rowdata,
                    cell_metadata = coldata)


# Look at the distribution of hash enrichment ratios

pdf("plots/enrichment_ratio_histogram.pdf", 
    width = 1.5, height = 1.5)
colData(cds) %>%
  as.data.frame() %>%
  filter(!is.na(ratio_top_2)) %>%
  ggplot() + 
  geom_histogram(aes(x = log10(ratio_top_2)),
                 color = "white",
                 fill = "grey50",
                 size = 0.05) + 
  monocle3:::monocle_theme_opts() +
  geom_vline(xintercept = log10(5), color = "#FF0000", size = .5) +
  theme(legend.position = "none", 
        text = element_text(size = 6),  
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.2,"line"),
        legend.text.align = 0) +
  clear_theme() +
  white_theme() +
  xlab("log10(Enrichment Ratio)") +
  ylab("Cells")
dev.off()



pdf("plots/doublets_vs_singlets_barplot.pdf", 
    width = 1.5, height = 1.5)
colData(cds) %>%
  as.data.frame() %>%
  filter(!is.na(ratio_top_2)) %>%
  mutate(singlet = 
           ifelse(ratio_top_2 >= 5,
                  "Singlet",
                  "Doublet")) %>%
  ggplot() + 
  geom_bar(aes(x = singlet,
               fill = singlet),
           size = 0.25,
           color = "white") + 
  monocle3:::monocle_theme_opts() +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none", 
        text = element_text(size = 6),  
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.2,"line"),
        legend.text.align = 0,
        axis.title.x = element_blank()) +
  clear_theme()+
  white_theme() +
  xlab("") +
  ylab("Cells")
dev.off()

cells_to_keep = 
  colData(cds) %>%
  as.data.frame() %>%
  filter(!is.na(ratio_top_2),
         ratio_top_2 >= 5) %>%
  pull(Cell) %>%
  as.character()

cds = cds[,cells_to_keep]


# Add coding variables into coldata ---------------------------------------

colData(cds)$data_set = 
  ifelse(grepl(x = colData(cds)$max_id,
               pattern = "day"),
         "EB differentiation",
         "Embryo")


condition_key = 
  data.frame(
    max_id = colData(cds) %>%
    as.data.frame() %>%
    filter(grepl(pattern = "day",
                 x = max_id)) %>%
    pull(max_id) %>%
    unique())

condition_key$day = 
  stringr::str_split_fixed(condition_key$max_id,
                           "_",
                           n = 2)[,1] %>%
  stringr::str_replace_all(pattern = "day",
                           replacement = "")

condition_key$condition = 
  stringr::str_split_fixed(condition_key$max_id,
                           "_",
                           n = 2)[,2] 

condition_key = 
  condition_key %>%
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

condition_matrix =
  condition_key %>%
  mutate(new_val = BMP+Activin,
         condition = as.numeric(condition)) %>%
  dplyr::select(new_val,condition) %>%
  distinct() %>%
  arrange(condition) %>%
  pull(new_val) %>%
  matrix(nrow = 4)

library(pheatmap)
library(viridis)
pheatmap(mat = log(condition_matrix),
         color = viridis(n = 16),
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         legend = F,
         cellwidth = 10,
         cellheight = 10,
         border_color = "black",
         filename = "plots/heatmap.pdf")

  
  
coldata =
  colData(cds) %>%
  as.data.frame() %>%
  left_join(condition_key) %>%
  data.frame(row.names = colData(cds)$Cell)

cds = 
  new_cell_data_set(expression_data = counts(cds),
                    gene_metadata = rowdata,
                    cell_metadata = coldata)

cds = 
  cds %>%
  estimate_size_factors() %>%
  preprocess_cds(num_dim = 100) %>%
  reduce_dimension(umap.n_neighbors = 10L) %>%
  cluster_cells()

colData(cds)$umap1 = reducedDim(cds,
                                type = "UMAP")[,1]
colData(cds)$umap2 = reducedDim(cds,
                                type = "UMAP")[,2]

colData(cds)$Cluster = clusters(cds,
                                reduction_method = "UMAP")
# Make summary plots based on UMAP embedding ------------------------------


ggplot() +
  geom_point(data = 
               colData(cds) %>%
               as.data.frame() %>%
               arrange(data_set),
             aes(x = umap1,
                 y = umap2,
                 color = data_set),
             size = 0.35,
             stroke = 0) +
  monocle3:::monocle_theme_opts() +
  clear_theme() +
  no_axes() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("EB differentiation" = "#B35EA5",
                                "Embryo" = "#71BF44"))+
  ggsave("plots/dataset_umap.pdf",
         height = 4,
         width = 4)
  

pdf("plots/dataset_barplot.pdf", 
    width = 1.5, height = 1.5)
colData(cds) %>%
  as.data.frame() %>%
  ggplot() + 
  geom_bar(aes(x = data_set,
               fill = data_set),
           size = 0.25,
           color = "white") + 
  monocle3:::monocle_theme_opts() +
  scale_fill_manual(values = c("EB differentiation" = "#B35EA5",
                                "Embryo" = "#71BF44"))+
  theme(legend.position = "none", 
        text = element_text(size = 6),  
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.2,"line"),
        legend.text.align = 0,
        axis.title.x = element_blank()) +
  clear_theme() +
  white_theme() +
  ylab("Cells")
dev.off()



ggplot() +
  geom_point(data = 
               colData(cds) %>%
               as.data.frame() %>%
               filter(data_set == "EB differentiation"),
             aes(x = umap1,
                 y = umap2,
                 color = day),
             size = 0.25,
             stroke = 0) +
  monocle3:::monocle_theme_opts() +
  clear_theme() +
  no_axes() +
  scale_color_manual(values = c("4" = "#ffeda0",
                                "5" = "#feb24c",
                                "6" = "#f03b20")) +
  theme(legend.position = "none") +
  ggsave("plots/EB_day_umap.pdf",
         height = 4,
         width = 4)


ggplot() +
  geom_point(data = 
               colData(cds) %>%
               as.data.frame() %>%
               filter(data_set == "EB differentiation") %>%
               dplyr::select(-Activin,
                             -BMP),
             aes(x = umap1,
                 y = umap2),
             color = "grey80",
             size = 0.25,
             stroke = 0) +
  geom_point(data = 
               colData(cds) %>%
               as.data.frame() %>%
               filter(data_set == "EB differentiation"),
             aes(x = umap1,
                 y = umap2,
                 color = day),
             size = 0.25,
             stroke = 0) +
  monocle3:::monocle_theme_opts() +
  clear_theme() +
  no_axes() +
  scale_color_manual(values = c("4" = "#ffeda0",
                                "5" = "#feb24c",
                                "6" = "#f03b20")) +
  facet_grid(Activin~BMP) +
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  ggsave("plots/EB_activin_bmp_umap.png",
         height = 8,
         width = 8,
         dpi = 600,
         bg = "transparent")



ggplot() +
  geom_point(data = 
               colData(cds) %>%
               as.data.frame() %>%
               dplyr::select(-Activin,
                             -BMP),
             aes(x = umap1,
                 y = umap2),
             color = "grey80",
             size = 0.55,
             stroke = 0) +
  geom_point(data = 
               colData(cds) %>%
               as.data.frame() %>%
               filter(data_set == "EB differentiation",
                      Activin == 0),
             aes(x = umap1,
                 y = umap2,
                 color = day),
             size = 0.55,
             stroke = 0) +
  monocle3:::monocle_theme_opts() +
  clear_theme() +
  no_axes() +
  scale_color_manual(values = c("4" = "#ffeda0",
                                "5" = "#feb24c",
                                "6" = "#f03b20")) +
  facet_wrap(~BMP,
             nrow = 1) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank()) +
  ggsave("plots/EB_activin0_bmp_umap.png",
         height = 2.5,
         width = 8,
         dpi = 600,
         bg = "transparent")


p_1 = 
  plot_cells(cds,genes = c("Kdr",
                         "Hoxb6",
                         "Runx1",
                         "Mecom",
                         "Lefty2",
                         "Mesp1"),
           cell_size = 0.25,
           cell_stroke = 0,
           label_cell_groups = F,
           label_branch_points = F,
           label_groups_by_cluster = F) +
  no_axes() +
  clear_theme() +
  theme(strip.background = element_blank(),
        legend.position = "none")

p_1$facet$params$nrow = 1
p_1 +
  ggsave("plots/umap_hsc2_genes.png",
         height = 2.5,
         width = 8,
         dpi = 600,
         bg = "transparent")

p_2 = 
  plot_cells(cds,genes = c("Dll4",
                        "Hey1",
                        "Nrp1",
                        "Efnb2"),
           cell_size = 0.25,
           cell_stroke = 0,
           label_cell_groups = F,
           label_branch_points = F,
           label_groups_by_cluster = F) +
  no_axes() +
  clear_theme() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

p_2$facet$params$nrow = 1
p_2 +
  ggsave("plots/umap_art_genes.png",
         height = 2.5,
         width = 8,
         dpi = 600,
         bg = "transparent")
  


# Perform Alignment -------------------------------------------------------


# Load in data from Pijuan-Sala 2019 ----------------------------
embryo_cds = 
  readRDS(file = "~/projects/hsc/published_data/Pijuan_Sala_2019/data/wild_tpye_atlas/atlas.monocle3.RDS")


# Join the two datasets ---------------------------------------------------

rowData(cds)$id_2 = 
  stringr::str_split_fixed(rowData(cds)$id,
                           "\\.",
                           2)[,1]

intersecting_genes = 
  intersect(rowData(embryo_cds)$id %>%
              as.character(),
            rowData(cds)$id_2)

rownames(cds) =rowData(cds)$id_2


cds_intersected = 
  cds[intersecting_genes,]

embryo_cds = embryo_cds[intersecting_genes,]

identical(rownames(cds_intersected) %>%
            as.character(),
          rownames(embryo_cds) %>%
            as.character())


joint_coldata = 
  rbind(colData(embryo_cds) %>%
          as.data.frame() %>%
          mutate(Cell = 
                   colnames(embryo_cds) %>%
                   as.character()) %>%
          dplyr::select(Cell,
                        annotation = stage) %>%
          mutate(sample = "Embryo Atlas"),
        colData(cds_intersected) %>%
          as.data.frame() %>%
          dplyr::select(Cell,
                        annotation = max_id) %>%
          mutate(sample = "EB Differentiation")
  )


rownames(joint_coldata) = 
  joint_coldata$Cell

rowdata =
  rowData(embryo_cds) %>% 
  as.data.frame()


joint_cds =
  new_cell_data_set(expression_data = cbind(counts(embryo_cds), 
                                            counts(cds_intersected)),
                    cell_metadata = joint_coldata,
                    gene_metadata = rowdata)


# Run Seurat CCA + Data Integration ---------------------------------------
library(scater)
library(Seurat)

count_mat = assay(joint_cds)
coldata_df = colData(joint_cds) %>% as.data.frame()
coldata_df$Cell = rownames(coldata_df)

cds_seurat = 
  CreateSeuratObject(counts = count_mat,
                     project = "hsc",
                     assay = "RNA",
                     meta.data = coldata_df)


# split dataset by experiment
cds.list <- SplitObject(cds_seurat, split.by = "sample")

# normalize data and find variable genes
for (i in 1:length(cds.list)) {
  cds.list[[i]] <- NormalizeData(cds.list[[i]], verbose = FALSE)
  cds.list[[i]] <- FindVariableFeatures(cds.list[[i]], selection.method = "vst", 
                                        nfeatures = 2000, verbose = FALSE)
}

cds.anchors <- 
  FindIntegrationAnchors(object.list = cds.list, 
                         dims = 1:100, 
                         anchor.features = 1000,
                         reduction = "cca",
                         verbose = T)

# integrate data 
cds.integrated <- IntegrateData(cds.anchors, 
                                dims = 1:100)
cds.integrated <- ScaleData(cds.integrated, verbose = FALSE)
cds.integrated <- RunPCA(cds.integrated, verbose = FALSE)
cds.integrated <- RunUMAP(cds.integrated, dims = 1:50,n.components = 2L)


DimPlot(cds.integrated, group.by = "sample",split.by = "annotation",ncol = 6) +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("plots/Seurat_integration_atlas_annotation.pdf",
         height = 20,
         width = 20,
         limitsize = F)

DimPlot(cds.integrated, group.by = "sample") +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("plots/Seurat_integration_atlas_annotation_color.pdf",
         height = 5,
         width = 5,
         limitsize = F)


seurat_umap_embeddings = 
  as.data.frame(Embeddings(cds.integrated, reduction = "umap"))

coldata_with_seurat_umap_coordinates = 
  left_join(coldata_df,
            seurat_umap_embeddings %>% 
              rownames_to_column(var = "Cell"),
            by = "Cell")



# Transfer seurat integrated coordinates to monocle -----------------------


joint_cds = 
  joint_cds %>%
  detect_genes() %>%
  estimate_size_factors() %>%
  preprocess_cds() %>%
  reduce_dimension()

identical(colnames(joint_cds) %>% as.character(),
          coldata_with_seurat_umap_coordinates$Cell %>% as.character())

reducedDim(x = joint_cds,
           type = "UMAP") <-
  matrix(cbind(coldata_with_seurat_umap_coordinates$UMAP_1,
               coldata_with_seurat_umap_coordinates$UMAP_2), 
         ncol=2)

joint_cds =
  joint_cds %>%
  cluster_cells()

saveRDS(object = joint_cds, file = "joint_cds_atlas_EB.RDS")


colData(joint_cds)$umap1 =
  reducedDim(joint_cds,
             type = "UMAP")[,1]

colData(joint_cds)$umap2 =
  reducedDim(joint_cds,
             type = "UMAP")[,2]


colData(joint_cds)$Cluster =
  clusters(joint_cds,
           reduction_method=  "UMAP")
# annotation --------------------------------------------------------------

coldata_joint = 
  colData(joint_cds) %>%
  as.data.frame()


coldata_atlas = 
  read.table(file = "/net/trapnell/vol1/home/sanjays/projects/hsc/published_data/Pijuan_Sala_2019/data/wild_tpye_atlas/meta.tab", 
             header = T, 
             sep = "\t")


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


coldata_joint =
  coldata_joint %>%
  left_join(coldata_atlas %>%
              dplyr::select(Cell,
                            umap_1_pjs = umapX,
                            umap_2_pjs = umapY,
                            celltype_pjs = celltype,
                            cluster_pjs = cluster,
                            theiler_pjs = theiler),
            by = "Cell")

coldata_joint = 
  coldata_joint %>%
  left_join(colData(cds) %>%
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
               coldata_joint %>%
               sample_n(size = nrow(coldata_joint),
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


label_position_df = 
  coldata_joint %>%
  group_by(celltype_pjs) %>%
  summarise(med_umap1 = median(umap1),
            med_umap2 = median(umap2)) %>%
  drop_na()


coldata_joint %>%
  mutate()
  filter(grepl(x = annotation,
               pattern = "E"))

ggplot() +
  geom_point(data = 
               coldata_joint %>%
               filter(is.na(dataset_eb)),
             aes(x = umap1,
                 y = umap2,
                 color = celltype_pjs),
             size = 0.55,
             stroke = 0) +
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
  monocle3:::monocle_theme_opts() +
  clear_theme() +
  no_axes() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank()) +
  ggsave("plots/joint_atlas_celltype_umap.png",
         height = 5,
         width = 5,
         dpi = 600,
         bg = "transparent")




pp_1 = 
  plot_cells(joint_cds,
             genes = c("Kdr",
                       "Hbb-bh1",
                       "Runx1",
                       "Mecom"),
             cell_size = 0.25,
             cell_stroke = 0,
             label_cell_groups = F,
             label_branch_points = F,
             label_groups_by_cluster = F) +
  no_axes() +
  clear_theme() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

pp_1$facet$params$nrow = 1
pp_1 +
  ggsave("plots/joint_umap_hsc_genes.png",
         height = 2.5,
         width = 8,
         dpi = 600,
         bg = "transparent")

pp_2 = 
  plot_cells(joint_cds,
             genes = c("Dll4",
                       "Hey1",
                       "Nrp1",
                       "Efnb2"),
             cell_size = 0.25,
             cell_stroke = 0,
             label_cell_groups = F,
             label_branch_points = F,
             label_groups_by_cluster = F) +
  no_axes() +
  clear_theme() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

pp_2$facet$params$nrow = 1
pp_2 +
  ggsave("plots/joint_umap_art_genes.png",
         height = 2.5,
         width = 8,
         dpi = 600,
         bg = "transparent")


coldata_joint %>%
  filter(!is.na(celltype_pjs)) %>%
  mutate(annotation = factor(x = annotation,
                           levels = c("mixed_gastrulation",
                                      "E6.5","E6.75","E7.0","E7.25",
                                      "E7.5","E7.75","E8.0","E8.25","E8.5"))) %>%
  sample_n(size = sum(!is.na(celltype_pjs))) %>%
  ggplot() +
  geom_point(aes(x = umap1,
                 y = umap2,
                 color = annotation),
             size = 0.55,
             stroke = 0) +
  scale_color_viridis_d() +
  monocle3:::monocle_theme_opts() +
  no_axes() +
  clear_theme() +
  theme(legend.position = "none") +
  ggsave("plots/joint_umap_stage.png",
         height = 5,
         width = 5,
         dpi = 600,
         bg = "transparent")


# pull in other mouse hematopoiesis data ----------------------------------

all_mouse_data = readRDS("/net/trapnell/vol1/home/sanjays/projects/hsc/RDS_output/hadland_atlas_joint_zhu.RDS")
  

cells_to_analyze = 
  read.table("/net/trapnell/vol1/home/sanjays/projects/hsc/RDS_output/hemato_cells.txt",
             sep = "\t",
             header = T)[,1] %>%
  as.character()

hemato.cds = all_mouse_data[,cells_to_analyze]


all_mouse_data_no_tal1_ko = 
  all_mouse_data[,colData(all_mouse_data)$tal1]

all_mouse_data_no_tal1_ko =
  all_mouse_data_no_tal1_ko[,!(colData(all_mouse_data_no_tal1_ko)$celltype %in% c("Doublet","Stripped", "mESC_Flk1+VEC+")) ]


