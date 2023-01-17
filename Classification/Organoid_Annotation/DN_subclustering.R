### Date: 1 March, 2022
### Author: Drew Neavin
### Purpose: Test for sybytes of neurons with clustering and plotting markers of interest

library(tidyverse)
library(data.table)
library(Seurat)
library(colorspace)
library(RColorBrewer)
library(ggeasy)
library(Nebulosa)



##### Set up Directories #####
dir <- "/path/to/base/directory/"
neuron_dir <- paste0(dir,"/output/Classification/Neuron_Clustering/Train_Test_Reannotation/")
outdir <- paste0(dir,"output/Classification/Neuron_Clustering/DN_subclustering/")

dir.create(outdir)


##### Read in the Data #####
seurat_neuron <- readRDS(paste0(neuron_dir, "neurons_classified_updated_names_sct_map.rds"))



##### Subset just the DNs #####
seurat_DNs <- subset(seurat_neuron, subset = Broad_Neurons == "DNs")


##### Normalize etc #####
seurat_DNs <- SCTransform(seurat_DNs, verbose = TRUE)
seurat_DNs <- RunPCA(seurat_DNs, npcs = 100)
seurat_DNs <- RunUMAP(seurat_DNs, reduction = "pca",dims = 1:100)
seurat_DNs <- FindNeighbors(seurat_DNs, reduction = "pca", dims = 1:100)


### Try multiple clustering categories
seurat_DNs_list <- list()
plot_list <- list()


clust_dir <- paste0(outdir,"Clustering/")
dir.create(clust_dir, recursive = TRUE)

DefaultAssay(seurat_DNs) <- "SCT"


for (resolution in c(seq(0.01,0.09,0.01),seq(0.1,2, by = 0.1))){
	seurat_DNs_list[[as.character(resolution)]] <- FindClusters(seurat_DNs, resolution = as.numeric(resolution))

	plot_list[[as.character(resolution)]] <- DimPlot(seurat_DNs_list[[as.character(resolution)]], reduction = "umap", label = TRUE, repel = TRUE)

	ggsave(plot_list[[as.character(resolution)]], filename = paste0(clust_dir, "Unassigned_clusters_resolution_",resolution, ".png"))
}




##### Make figures with markers #####
### General Markers ###
markers <- fread("genes_of_interest.tsv") ## provided as supplementary table in the manuscript

neuron_marker_dir <- paste0(outdir, "markers/")
dir.create(neuron_marker_dir)

DefaultAssay(seurat_DNs) <- "RNA"

for (gene in unique(markers$ENSG_ID)){
    gene_id <- markers[ENSG_ID == gene]$Gene
	if (gene %in% rownames(seurat_DNs)){
		if (rowSums(seurat_DNs[gene,]) > 0){
			if(!file.exists(paste0(neuron_marker_dir,gene_id,"_umap_seurat_sub.png"))){
				print(gene_id)
				umap <- FeaturePlot(seurat_DNs, features = gene) + labs(title = gene_id)
				ggsave(umap, filename = paste0(neuron_marker_dir, gene_id,"_umap_seurat_sub.png"))
			}
			if(!file.exists(paste0(neuron_marker_dir,gene_id,"_umap_nebulosa.png"))){
				density_plot <- plot_density(seurat_DNs, gene, pal = "plasma", reduction = "umap") + labs(title = gene_id)
				ggsave(density_plot, filename = paste0(neuron_marker_dir,gene_id,"_umap_nebulosa.png"))
			# }
		}
	}
}


### Will use resolution 0.8 and classified the cells ###
cluster_classifications <- fread("DN_subtypes_clusters.tsv", sep = "\t") ## provided on github
cluster_classifications$Idents <- as.character(cluster_classifications$Idents)


## Add DN subtypes to seurat object ###
seurat_DNs_selected <- seurat_DNs_list[["0.8"]]

idents <- data.table(Barcode = colnames(seurat_DNs_selected), Idents = as.character(Idents(seurat_DNs_selected)))


DN_subtype_table <- cluster_classifications[,c("Idents", "DN_CellType")][idents, on = "Idents"]

DN_subtype_table <- data.frame(DN_subtype_table)
rownames(DN_subtype_table) <- DN_subtype_table$Barcode
DN_subtype_table$Idents <- NULL
DN_subtype_table$Barcode <- NULL

seurat_DNs_selected <- AddMetaData(seurat_DNs_selected, DN_subtype_table)

saveRDS(seurat_DNs_selected, paste0(outdir, "seurat_DNs_classified.rds"))

DN_colors <- cluster_classifications$Color
names(DN_colors) <- cluster_classifications$DN_CellType


##### Make UMAP #####
umap_DNs <- DimPlot(seurat_DNs_selected, group.by = "DN_CellType", label = TRUE, repel = TRUE) +
				ggtitle("Dopaminergic Neuron\nSubtypes") +
				xlab("UMAP 1") +
				ylab("UMAP 2") +
				scale_color_manual(values = DN_colors)
					
						
ggsave(umap_DNs, filename = paste0(outdir,"umap_DNs_", Sys.Date(), ".png"), width = 5, height = 3.5)
ggsave(umap_DNs, filename = paste0(outdir,"umap_DNs_", Sys.Date(), ".pdf"), width = 5, height = 3.5)



umap_DNs_no_lab <- DimPlot(seurat_DNs_selected, group.by = "DN_CellType", label = FALSE, repel = TRUE) +
				ggtitle("Dopaminergic Neuron\nSubtypes") +
				xlab("UMAP 1") +
				ylab("UMAP 2") +
				scale_color_manual(values = DN_colors)


ggsave(umap_DNs_no_lab, filename = paste0(outdir,"umap_DNs_no_label_", Sys.Date(), ".png"), width = 5, height = 3.5)
ggsave(umap_DNs_no_lab, filename = paste0(outdir,"umap_DNs_no_label_", Sys.Date(), ".pdf"), width = 5, height = 3.5)




#### Test for DEGs between cell tyoes
deg_list_dir <- paste0(outdir, "DEGs/")
dir.create(deg_list_dir)

gene_conversions <- fread("/path/to/10x/matrix/features.tsv.gz", sep = "\t", col.names = c("Gene", "Gene_ID", "Assay")) ### So can convert ENSG IDs to gene IDs for easy interpretation of DEGs

Idents(seurat_DNs_selected) <- "DN_CellType"
deg_list <- list()

for (ident in unique(Idents(seurat_DNs_selected))){
	deg_list[[as.character(ident)]] <- FindMarkers(object = seurat_DNs_selected, ident.1 = ident, test.use = "MAST", latent.vars = c("percent.rb", "nCount_RNA"))
	deg_list[[as.character(ident)]]$Gene <- rownames(deg_list[[as.character(ident)]])

	### Add in Gene IDs ###
	deg_list[[as.character(ident)]] <- gene_conversions[deg_list[[as.character(ident)]], on = "Gene"]
	deg_list[[as.character(ident)]]$Assay <- NULL

	fwrite(deg_list[[as.character(ident)]], paste0(deg_list_dir, ident, "_DEGs.tsv"), sep = "\t")
}