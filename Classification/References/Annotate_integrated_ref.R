#############################
### Author: Drew  Neavin
### Date: 16 August, 2021
### Reason: Reannotate the cells based on the integrated reference
#############################



library(Seurat)
library(data.table)
library(ggplot2)
library(ggridges)
library(tidyverse)
library(Nebulosa)



##### Set up Directories #####
dir <- "/path/to/base/directory/"
int_dir <- paste0(dir,"output/Classification/Integrate_refs/")
outdir <- paste0(dir, "output/Classification/annotation/")
dir.create(outdir, recursive = TRUE)



##### Read in Data #####
seurat_integrated <- readRDS(paste0(int_dir, "integrated_human_mouse_csc_nature.rds"))



##### Cluster to identify neurons and non-neuons to separate #####
clust_dir <- paste0(outdir, "Cluster_together/")
dir.create(clust_dir)

seurat_integrated_sct_subset_list <- list()
plot_list <- list()

DefaultAssay(seurat_integrated_new_sct) <- "integrated"


for (resolution in c(seq(0.01,0.09,0.01),seq(0.1,2, by = 0.1))){
	seurat_integrated_sct_subset_list[[as.character(resolution)]] <- FindClusters(seurat_integrated_new_sct, resolution = resolution)

	plot_list[[as.character(resolution)]] <- DimPlot(seurat_integrated_sct_subset_list[[as.character(resolution)]], reduction = "umap", label = TRUE, repel = TRUE)

	ggsave(plot_list[[as.character(resolution)]], filename = paste0(clust_dir, "Integrated_clusters_resolution_",resolution, ".png"))
}



##### Make figures with markers to identify cell types #####
### General Markers ###
markers <- fread("genes_of_interest.tsv") ### List of markers tested can be found in supplementar tables for this manuscript

marker_dir <- paste0(outdir, "markers/")
dir.create(marker_dir)

DefaultAssay(seurat_integrated_new_sct) <- "RNA"

for (gene in unique(markers$Gene)){
	if (gene %in% rownames(seurat_integrated_new_sct)){
		if(!file.exists(paste0(marker_dir,gene,"_umap.png"))){
			umap <- FeaturePlot(seurat_integrated_new_sct, features = gene) + labs(title = gene)
			ggsave(umap, filename = paste0(marker_dir, gene,"_umap_seurat_sub.png"))
		}
		if(!file.exists(paste0(marker_dir,gene,"_umap_nebulosa.png"))){
			density_plot <- plot_density(seurat_integrated_new_sct, gene, pal = "plasma") + labs(title = gene)
			ggsave(density_plot, filename = paste0(marker_dir,gene,"_umap_nebulosa.png"))
		}
	}
}



##### Separate neurons and non-neurons and remake umap + cluster #####
seurat_integrated_new_sct_list <- list()

seurat_integrated_new_sct_list[["neurons"]] <- subset(seurat_integrated_sct_subset_list[["0.02"]], idents = 0)
seurat_integrated_new_sct_list[["non-neurons"]] <- subset(seurat_integrated_sct_subset_list[["0.02"]], idents = c(1,2))



### Integrate the neurons and non-neurons groups separately for clustering ###
## Subset the studies ##
seurat_integrated_new_sct_study_list <- list()
seurat_integrated_new_sct_study_list <- lapply(names(seurat_integrated_new_sct_list), function(x){
	for (ref in unique(seurat_integrated_new_sct_list[[x]]$Origin)){
		seurat_integrated_new_sct_study_list[[x]][[ref]] <- subset(seurat_integrated_new_sct_list[[x]], subset = Origin == ref)
	}
	return(seurat_integrated_new_sct_study_list[[x]])
})
names(seurat_integrated_new_sct_study_list) <- names(seurat_integrated_new_sct_list)

## Normalize data per study ##
seurat_integrated_new_sct_study_list <- lapply(seurat_integrated_new_sct_study_list, function(x){
	lapply(x, FUN = SCTransform)
})

## Identify features ##
integration_features_study <- lapply(seurat_integrated_new_sct_study_list, function(x){
	SelectIntegrationFeatures(object.list = x, nfeatures = 3000)
})
names(integration_features_study) <- names(seurat_integrated_new_sct_study_list)

## prep for integration ##
seurat_integrated_new_sct_study_list <- lapply(names(integration_features_study), function(x){
	seurat_integrated_new_sct_study_list[[x]] <- PrepSCTIntegration(object.list = seurat_integrated_new_sct_study_list[[x]], anchor.features = integration_features_study[[x]])
	return(seurat_integrated_new_sct_study_list[[x]])
})
names(seurat_integrated_new_sct_study_list) <- names(integration_features_study)

## identify anchors ##
anchors_study_list <- lapply(names(integration_features_study), function(x){
	FindIntegrationAnchors(object.list = seurat_integrated_new_sct_study_list[[x]], normalization.method = "SCT",anchor.features = integration_features_study[[x]])
})
names(anchors_study_list) <- names(integration_features_study)

## Integrate data ##
seurat_integrated_neuron_nonneuron_list <- lapply(anchors_study_list, function(x){
	IntegrateData(anchorset = x, normalization.method = "SCT")
})

## Run PCA, UMAP and find neighbors for integrated neurons and non-neurons ##
seurat_integrated_neuron_nonneuron_list <- lapply(seurat_integrated_neuron_nonneuron_list, function(x){
	RunPCA(x, npcs = 100)
})

seurat_integrated_neuron_nonneuron_list <- lapply(seurat_integrated_neuron_nonneuron_list, function(x){
	RunUMAP(x, reduction = "pca", dims = 1:100, min.dist = 0.2, spread = 1.5)
})

seurat_integrated_neuron_nonneuron_list <- lapply(seurat_integrated_neuron_nonneuron_list, function(x){
	FindNeighbors(x, reduction = "pca", dims = 1:100)
})



##### Cluster integrated datasets to identify common cell types - test multiple resolutions to identify the best for annotating cell types #####
seurat_integrated_sct_neuron_nonneuron_subset_list <- list()
plot_list <- list()

	### General Markers ###
	markers <- fread("genes_of_interest.tsv") ### List of markers tested can be found in supplementar tables for this manuscript



for( x in names(seurat_integrated_neuron_nonneuron_list)){
	clust_dir <- paste0(outdir,x,"/Clustering/")
	dir.create(clust_dir, recursive = TRUE)


	DefaultAssay(seurat_integrated_neuron_nonneuron_list[[x]]) <- "integrated"
	

	for (resolution in c(seq(0.01,0.09,0.01),seq(0.1,2, by = 0.1))){
		seurat_integrated_sct_neuron_nonneuron_subset_list[[x]][[as.character(resolution)]] <- FindClusters(seurat_integrated_neuron_nonneuron_list[[x]], resolution = resolution)

		plot_list[[x]][[as.character(resolution)]] <- DimPlot(seurat_integrated_sct_neuron_nonneuron_subset_list[[x]][[as.character(resolution)]], reduction = "umap", label = TRUE, repel = TRUE)

		ggsave(plot_list[[x]][[as.character(resolution)]], filename = paste0(clust_dir, x, "_Integrated_clusters_resolution_",resolution, ".png"))
	}

	marker_dir_neuron_nonneuron <- paste0(outdir,x,"markers/")
	dir.create(marker_dir_neuron_nonneuron)

	##### Make figures with markers to identify cell types #####
	DefaultAssay(seurat_integrated_neuron_nonneuron_list[[x]]) <-  "RNA"

	for (gene in unique(markers$Gene)){
		if (gene %in% rownames(seurat_integrated_neuron_nonneuron_list[[x]])){
			if(!file.exists(paste0(marker_dir_neuron_nonneuron,gene,"_umap.png"))){
				umap <- FeaturePlot(seurat_integrated_neuron_nonneuron_list[[x]], features = gene) + labs(title = gene)
				ggsave(umap, filename = paste0(marker_dir_neuron_nonneuron, gene,"_umap_seurat_sub.png"))
			}
			if(!file.exists(paste0(marker_dir_neuron_nonneuron,gene,"_umap_nebulosa.png"))){
				density_plot <- plot_density(seurat_integrated_neuron_nonneuron_list[[x]], gene, pal = "plasma") + labs(title = gene)
				ggsave(density_plot, filename = paste0(marker_dir_neuron_nonneuron,gene,"_umap_nebulosa.png"))
			}
		}
	}
}



### Assign some cell types, pull out some clusters and recluster
Idents(seurat_integrated_neuron_nonneuron_list[["neurons"]]) <- Idents(seurat_integrated_sct_neuron_nonneuron_subset_list[["neurons"]][["1.9"]])
seurat_integrated_neuron_nonneuron_list[["neurons"]]$CellType <- ifelse((Idents(seurat_integrated_neuron_nonneuron_list[["neurons"]]) %in% c(21,40,58)), "IPCs", "Neurons")

non_neuron_class <- data.frame(ident = seq(0,11), CellType = c("Astrocytes & Radial Glia", "OPC",  "Oligodendrocytes", "Radial Glial Cells","Tanycytes & Radial Glia", "Astrocytes", "Ependymal Cells", "Immature Oligodendrocytes & Oligodendrocytes", "Microglia", "Endothelial", "VLMCs & Unknown", "Astrocytes"))
non_neuron_class <- data.frame(left_join(data.frame(ident = as.numeric(as.character(Idents(seurat_integrated_sct_neuron_nonneuron_subset_list[["non-neurons"]][["0.2"]])))), non_neuron_class))
rownames(non_neuron_class) <- colnames(seurat_integrated_sct_neuron_nonneuron_subset_list[["non-neurons"]][["0.2"]])

seurat_integrated_neuron_nonneuron_list[["non-neurons"]] <- AddMetaData(seurat_integrated_neuron_nonneuron_list[["non-neurons"]], non_neuron_class)

### Recluster for idents 0,4,7 and 10 ###
seurat_integrated_nonneuron_list <- list()

for (group in unique(seurat_integrated_neuron_nonneuron_list[["non-neurons"]]$CellType)[grep("&", unique(seurat_integrated_neuron_nonneuron_list[["non-neurons"]]$CellType))]){
	seurat_integrated_nonneuron_list[[group]] <- subset(seurat_integrated_neuron_nonneuron_list[["non-neurons"]], subset = CellType == group)
}


seurat_integrated_nonneuron_list <- lapply(seurat_integrated_nonneuron_list, function(x){
	RunPCA(x, npcs = 100)
})

seurat_integrated_nonneuron_list[["Astrocytes & Radial Glia"]] <- RunUMAP(seurat_integrated_nonneuron_list[["Astrocytes & Radial Glia"]], reduction = "pca", dims = 1:100)
seurat_integrated_nonneuron_list[["Immature Oligodendrocytes & Oligodendrocytes"]] <- RunUMAP(seurat_integrated_nonneuron_list[["Immature Oligodendrocytes & Oligodendrocytes"]], reduction = "pca", dims = 1:30)
seurat_integrated_nonneuron_list[["Tanycytes & Radial Glia"]] <- RunUMAP(seurat_integrated_nonneuron_list[["Tanycytes & Radial Glia"]], reduction = "pca", dims = 1:30)
seurat_integrated_nonneuron_list[["VLMCs & Unknown"]] <- RunUMAP(seurat_integrated_nonneuron_list[["VLMCs & Unknown"]], reduction = "pca", dims = 1:30)




seurat_integrated_nonneuron_list <- lapply(seurat_integrated_nonneuron_list, function(x){
	FindNeighbors(x, reduction = "pca", dims = 1:100)
})



### Test multiple differet resolutions to find the best that has clusters of independent cell types that can be annotated ###
for( x in names(seurat_integrated_nonneuron_list)){
	clust_dir <- paste0(outdir,"non-neurons/",x,"/Clustering/")
	dir.create(clust_dir, recursive = TRUE)


	DefaultAssay(seurat_integrated_nonneuron_list[[x]]) <- "integrated"
	
	for (resolution in c(seq(0.01,0.09,0.01),seq(0.1,2, by = 0.1))){
		seurat_integrated_sct_nonneuron_subset_list[[x]][[as.character(resolution)]] <- FindClusters(seurat_integrated_nonneuron_list[[x]], resolution = resolution)

		plot_list[[x]][[as.character(resolution)]] <- DimPlot(seurat_integrated_sct_nonneuron_subset_list[[x]][[as.character(resolution)]], reduction = "umap", label = TRUE, repel = TRUE)

		ggsave(plot_list[[x]][[as.character(resolution)]], filename = paste0(clust_dir, x, "_Integrated_clusters_resolution_",resolution, ".png"))
	}


	marker_dir_nonneurons <- paste0(outdir,"non-neurons/",x,"markers/")
	dir.create(marker_dir_neuron_nonneuron)

	##### Make figures with markers to identify cell types #####
	DefaultAssay(seurat_integrated_nonneuron_list[[x]]) <-  "RNA"

	for (gene in unique(markers$Gene)){
		if (gene %in% rownames(seurat_integrated_nonneuron_list[[x]])){
			if(!file.exists(paste0(marker_dir_nonneurons,gene,"_umap.png"))){
				umap <- FeaturePlot(seurat_integrated_nonneuron_list[[x]], features = gene) + labs(title = gene)
				ggsave(umap, filename = paste0(marker_dir_nonneurons, gene,"_umap_seurat_sub.png"))
			}
			if(!file.exists(paste0(marker_dir_nonneurons,gene,"_umap_nebulosa.png"))){
				density_plot <- plot_density(seurat_integrated_nonneuron_list[[x]], gene, pal = "plasma") + labs(title = gene)
				ggsave(density_plot, filename = paste0(marker_dir_nonneurons,gene,"_umap_nebulosa.png"))
			}
		}
	}
}





### Finalize cell type annotations
seurat_integrated_nonneuron_list[["Astrocytes & Radial Glia"]]$CellType <- ifelse(Idents(seurat_integrated_sct_nonneuron_subset_list[["Astrocytes & Radial Glia"]][["0.2"]]) == 3, "Radial Glial Cells", "Astrocytes")
seurat_integrated_nonneuron_list[["Immature Oligodendrocytes & Oligodendrocytes"]]$CellType <- ifelse(Idents(seurat_integrated_sct_nonneuron_subset_list[["Immature Oligodendrocytes & Oligodendrocytes"]][["1.5"]]) %in% c(1,6,8), "Oligodendrocytes", "Immature Oligodendrocytes")
seurat_integrated_nonneuron_list[["Tanycytes & Radial Glia"]]$CellType <- ifelse(Idents(seurat_integrated_sct_nonneuron_subset_list[["Tanycytes & Radial Glia"]][["0.7"]]) %in% c(0,4), "Tanycytes", 
																	ifelse(Idents(seurat_integrated_sct_nonneuron_subset_list[["Tanycytes & Radial Glia"]][["0.7"]]) == 2, "Vascular & Leptomeningeninal Cells (VLMCs)", "Radial Glial Cells"))
seurat_integrated_nonneuron_list[["VLMCs & Unknown"]]$CellType <- ifelse(Idents(seurat_integrated_sct_nonneuron_subset_list[["VLMCs & Unknown"]][["0.02"]]) == 0, "Vascular & Leptomeningeninal Cells (VLMCs)", "Radial Glial Cells")



### Make table of all the classifications
CellTypes_df <- data.table(rbind(data.frame(Barcode = colnames(seurat_integrated_neuron_nonneuron_list[["neurons"]]), CellType = seurat_integrated_neuron_nonneuron_list[["neurons"]]@meta.data[,c("CellType")]), data.frame(Barcode = colnames(seurat_integrated_neuron_nonneuron_list[["non-neurons"]]), CellType = seurat_integrated_neuron_nonneuron_list[["non-neurons"]]@meta.data[,c("CellType")])))
CellTypes_df <- CellTypes_df[!grepl("&", CellType)]

CellTypes_df_nonneuron_update <- data.table(rbind(data.frame(Barcode = colnames(seurat_integrated_nonneuron_list[["Astrocytes & Radial Glia"]]), CellType = seurat_integrated_nonneuron_list[["Astrocytes & Radial Glia"]]@meta.data[,c("CellType")]), 
													data.frame(Barcode = colnames(seurat_integrated_nonneuron_list[["Immature Oligodendrocytes & Oligodendrocytes"]]), CellType = seurat_integrated_nonneuron_list[["Immature Oligodendrocytes & Oligodendrocytes"]]@meta.data[,c("CellType")]),
														data.frame(Barcode = colnames(seurat_integrated_nonneuron_list[["Tanycytes & Radial Glia"]]), CellType = seurat_integrated_nonneuron_list[["Tanycytes & Radial Glia"]]@meta.data[,c("CellType")]),
															data.frame(Barcode = colnames(seurat_integrated_nonneuron_list[["VLMCs & Unknown"]]), CellType = seurat_integrated_nonneuron_list[["VLMCs & Unknown"]]@meta.data[,c("CellType")])))

CellTypes_df <- rbind(CellTypes_df, CellTypes_df_nonneuron_update)
rownames(CellTypes_df) <- CellTypes_df$Barcode

seurat_integrated_new_sct <- AddMetaData(seurat_integrated_new_sct, CellTypes_df)


saveRDS(seurat_integrated_new_sct, paste0(outdir, "integrated_refs_sct_reannotated.rds"))





### Read in doublet classifications froom DoubletDetecting.R script and remove doublets ###
doublets <- fread(paste0(outdir, "singlets_doublet.tsv"))
doublets <- data.frame(doublets)
rownames(doublets) <- doublets$Barcode
doublets$Barcode <- NULL

seurat_integrated_new_sct <- AddMetaData(seurat_integrated_new_sct, meta = doublets)
p_doublets <- DimPlot(seurat_integrated_new_sct, reduction = "umap", group.by = "Combined_DropletType")

ggsave(p_doublets, filename = paste0(outdir, "integrated_figures_sct_reannotated_doublets.png"), width = 11)


seurat_integrated_new_sct_sing <- subset(seurat_integrated_new_sct, subset = Combined_DropletType == "singlet")

DefaultAssay(seurat_integrated_new_sct_sing) <- "integrated"

seurat_integrated_new_sct_sing <- RunPCA(seurat_integrated_new_sct_sing, npcs = 100)
seurat_integrated_new_sct_sing <- RunUMAP(seurat_integrated_new_sct_sing, reduction = "pca", dims = 1:100, min.dist = 0.2, spread = 1.5)
seurat_integrated_new_sct_sing <- RunUMAP(seurat_integrated_new_sct_sing, reduction = "pca", dims = 1:100, min.dist = 0.15, spread = 1.2)
seurat_integrated_new_sct_sing <- FindNeighbors(seurat_integrated_new_sct_sing, reduction = "pca", dims = 1:100)


# Visualization
p1_sing <- DimPlot(seurat_integrated_new_sct_sing, reduction = "umap", group.by = "Origin")
p2_sing <- DimPlot(seurat_integrated_new_sct_sing, reduction = "umap", group.by = "age")
p3_sing <- DimPlot(seurat_integrated_new_sct_sing, reduction = "umap", group.by = "Cell_Types_Merged_Simple")
p4_sing <- DimPlot(seurat_integrated_new_sct_sing, reduction = "umap", group.by = "CellType")
ggsave((p1_sing + p2_sing + p3_sing + p4_sing), filename = paste0(outdir,"integrated_figures_reannotated_sing.png"), width = 15, height = 10)

ggsave(p4_sing, filename = paste0(outdir,"CellType_integrated_figures_sct_reannotated_sing.png"), width = 11)



saveRDS(seurat_integrated_new_sct_sing, paste0(outdir, "integrated_refs_sct_reannotated_singlets.rds"))

