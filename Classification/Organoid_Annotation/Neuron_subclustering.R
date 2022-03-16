#############################
### Author: Drew  Neavin
### Reason: Use scPred to predict cell types in the hypothalamus organoids
#############################
library(data.table)
library(Seurat)
library(tidyverse)
library(ggridges)
library(Nebulosa)



##### Set up Directories #####
dir <- "/path/to/base/directory/"
datadir <- paste0(dir, "output/Classification/annotation/")
outdir <- paste0(dir,"output/Classification/Neuron_Clustering/")

dir.create(outdir, recursive = TRUE)




##### Read in the data #####
##### Subset unassigned and cluster for cluster-based annotation #####
seurat <- readRDS(paste0(datadir, "integrated_refs_sct_reannotated_singlets.rds"))



seurat_updated <- SCTransform(seurat_updated, verbose = TRUE)
seurat_updated <- RunPCA(seurat_updated, npcs = 100)
seurat_updated <- RunUMAP(seurat_updated, reduction = "pca",dims = 1:100)
seurat_updated <- FindNeighbors(seurat_updated, reduction = "pca", dims = 1:100)

### Try multiple clustering categories
seurat_updated_list <- list()
plot_list <- list()


clust_dir <- paste0(outdir,"Clustering/")
dir.create(clust_dir, recursive = TRUE)


for (resolution in c(seq(0.01,0.09,0.01),seq(0.1,2, by = 0.1))){
	seurat_updated_list[[as.character(resolution)]] <- FindClusters(seurat_updated, resolution = resolution)

	plot_list[[as.character(resolution)]] <- DimPlot(seurat_updated_list[[as.character(resolution)]], reduction = "umap", label = TRUE, repel = TRUE)

	ggsave(plot_list[[as.character(resolution)]], filename = paste0(clust_dir, "Unassigned_clusters_resolution_",resolution, ".png"))
}




##### Check DEGs for each clusters:
resolution <- 1.7
deg_list <- list()


deg_list_dir <- paste0(outdir, "DEGs/")
dir.create(deg_list_dir)
gene_conversions <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/Parkinsons_TH_Collab/data/Remapping/remapped/1205_GEX_WT_THpos_GEX_0_1_HC3GJDSXY/outs/filtered_feature_bc_matrix/features.tsv.gz", sep = "\t", col.names = c("Gene", "Gene_ID", "Assay"))


for (ident in unique(c(0,1,2,3,7,10,12,15,23,26,29,17,Idents(seurat_updated_list[[as.character(resolution)]])))){
	deg_list[[as.character(ident)]] <- FindMarkers(object = seurat_updated_list[[resolution]], ident.1 = ident)
	deg_list[[as.character(ident)]]$Gene <- rownames(deg_list[[as.character(ident)]])

	### Add in Gene IDs ###
	deg_list[[as.character(ident)]] <- gene_conversions[deg_list[[as.character(ident)]], on = "Gene"]
	deg_list[[as.character(ident)]]$Assay <- NULL

	fwrite(deg_list[[as.character(ident)]], paste0(deg_list_dir, "cluster_", ident, "_resolution_",resolution,"_DEGs.tsv"), sep = "\t")
}

saveRDS(seurat_updated_list[[as.character(resolution)]], paste0(outdir,"seurat_resolution_", resolution, ".rds"))
