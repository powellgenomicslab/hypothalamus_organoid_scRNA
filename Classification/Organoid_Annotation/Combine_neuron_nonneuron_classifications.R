#############################
### Author: Drew  Neavin
### Reason: Combine the nueonr subtype annotations to the main seurat object that has all cells
#############################

library(tidyverse)
library(data.table)
library(Seurat)
library(colorspace)
library(RColorBrewer)
library(ggeasy)
library(ggridges)



##### Set up Directories #####
dir <- "/path/to/base/directory/"
neuron_dir <- paste0(dir,"output/Classification/Neuron_Clustering/DN_subclustering/")
outdir <- paste0(dir,"output/Combined_Final_Classification/")
datadir <- paste0(dir, "output/Classification/scPred/")

dir.create(outdir)


##### Read in the Data #####
seurat_neuron <- readRDS(paste0(neuron_dir, "seurat_DNs_classified.rds"))
seurat <- readRDS(paste0(datadir, "seurat_broad_classifications.rds"))



### Add the metadata from the neuron classifications to the seurat object ###
subtypes <- seurat_neuron@meta.data
neuron_subtypes <- data.frame(subtypes[,c("updated_names", "Broad_Neurons")])
rownames(neuron_subtypes) <- rownames(subtypes)
colnames(neuron_subtypes) <- c("Neuron_Subtypes_Specific","Neuron_Subtypes_Broad")

seurat <- AddMetaData(seurat, neuron_subtypes)



### Update Broad Cell Type Categories ###
seurat$CellType <- ifelse((seurat$CellType == "Immature Oligodendrocytes" | seurat$CellType == "Oligodendrocytes"), "Oligodendrocytes", seurat$CellType) ## since have very few Oligodendrocytes and Immature Oligodendrocytes, consider them together

saveRDS(seurat, paste0(outdir,"seurat_all_singlets_classified.rds"))




##### Make UMAPS #####
umap_orig <- DimPlot(subset(seurat, subset = CellType != "unassigned"), group.by = "CellType", label = TRUE, repel = TRUE) +
				scale_color_manual(values = broad_celltype_colors, name = "Cell\nType")
ggsave(umap_orig, filename = paste0(outdir,"All_Cells_Broad_Cell_Types_UMAP.png"), width = 10)
ggsave(umap_orig, filename = paste0(outdir,"All_Cells_Broad_Cell_Types_UMAP.pdf"), width = 10)


umap_orig_no_lab <- DimPlot(subset(seurat, subset = CellType != "unassigned"), group.by = "CellType", label = FALSE, repel = TRUE) +
						scale_color_manual(values = broad_celltype_colors, name = "Cell\nType")
ggsave(umap_orig_no_lab, filename = paste0(outdir,"All_Cells_Broad_Cell_Types_UMAP_no_label.png"), width = 10)
ggsave(umap_orig_no_lab, filename = paste0(outdir,"All_Cells_Broad_Cell_Types_UMAP_no_label.pdf"), width = 10)


umap_orig_neuron_subtype <- DimPlot(subset(seurat, subset = CellType != "unassigned"), group.by = "Combined_Cell_Types", label = TRUE, repel = TRUE) +
				scale_color_manual(values = c(broad_celltype_colors,neuron_colors), name = "Cell\nType")
ggsave(umap_orig_neuron_subtype, filename = paste0(outdir,"All_Cells_Broad_Neuron_Cell_Types_UMAP.png"), width = 13)
ggsave(umap_orig_neuron_subtype, filename = paste0(outdir,"All_Cells_Broad_Neuron_Cell_Types_UMAP.pdf"), width = 13)


umap_orig_neuron_subtype_no_lab <- DimPlot(subset(seurat, subset = CellType != "unassigned"), group.by = "Combined_Cell_Types", label = FALSE, repel = TRUE) +
						scale_color_manual(values = c(broad_celltype_colors,neuron_colors), name = "Cell\nType")
ggsave(umap_orig_neuron_subtype_no_lab, filename = paste0(outdir,"All_Cells_Broad_Neuron_Cell_Types_UMAP_no_label.png"), width = 13)
ggsave(umap_orig_neuron_subtype_no_lab, filename = paste0(outdir,"All_Cells_Broad_Neuron_Cell_Types_UMAP_no_label.pdf"), width = 13)



