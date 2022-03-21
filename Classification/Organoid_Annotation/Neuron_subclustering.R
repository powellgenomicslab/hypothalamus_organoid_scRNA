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

### Cluster at multiple different resolutions to find the one that best separtes unique cell types ###
seurat_updated_list <- list()
plot_list <- list()


clust_dir <- paste0(outdir,"Clustering/")
dir.create(clust_dir, recursive = TRUE)


for (resolution in c(seq(0.01,0.09,0.01),seq(0.1,2, by = 0.1))){
	seurat_updated_list[[as.character(resolution)]] <- FindClusters(seurat_updated, resolution = resolution)

	plot_list[[as.character(resolution)]] <- DimPlot(seurat_updated_list[[as.character(resolution)]], reduction = "umap", label = TRUE, repel = TRUE)

	ggsave(plot_list[[as.character(resolution)]], filename = paste0(clust_dir, "Unassigned_clusters_resolution_",resolution, ".png"))
}

### Resolution 1.7 was the most appropriate for the neurons ###
resolution <- 1.7


##### Check DEGs for each cluster to help with nueron annotation #####
deg_list <- list()


deg_list_dir <- paste0(outdir, "DEGs/")
dir.create(deg_list_dir)
gene_conversions <- fread("/path/to/10x/matrix/features.tsv.gz", sep = "\t", col.names = c("Gene", "Gene_ID", "Assay")) ### So can convert ENSG IDs to gene IDs for easy interpretation of DEGs


for (ident in unique(Idents(seurat_updated_list[[as.character(resolution)]]))){
	deg_list[[as.character(ident)]] <- FindMarkers(object = seurat_updated_list[[resolution]], ident.1 = ident)
	deg_list[[as.character(ident)]]$Gene <- rownames(deg_list[[as.character(ident)]])

	### Add in Gene IDs ###
	deg_list[[as.character(ident)]] <- gene_conversions[deg_list[[as.character(ident)]], on = "Gene"]
	deg_list[[as.character(ident)]]$Assay <- NULL

	fwrite(deg_list[[as.character(ident)]], paste0(deg_list_dir, "cluster_", ident, "_resolution_",resolution,"_DEGs.tsv"), sep = "\t")
}

saveRDS(seurat_updated_list[[as.character(resolution)]], paste0(outdir,"seurat_resolution_", resolution, ".rds"))



##### Subcluster clusters with "/" in them (because likely two different cell types in cluster) #####
seurat_updated_list_sub <- list()
plot_list <- list()


## Set up directory ##
neuron_marker_dir <- paste0(outdir, "Annotated_Subclustering/Markers/")
dir.create(neuron_marker_dir)


## Read in necessary files
markers <- fread("genes_of_interest.tsv") ### List of markers tested can be found in supplementar tables for this manuscript
cluster_anno_dt <- fread("cluster_annotations.tsv", sep = "\t") ### can be found on Github in Classification/Organoid_Annotation
cluster_anno_dt$Cluster_Name <- gsub(" ", "_", cluster_anno_dt$Cluster_Name) %>% gsub("/", ".", .)




for (clust_anno in grep("\\.", cluster_anno_dt$Cluster_Name, value = TRUE)){
    cluster <- cluster_anno_dt[Cluster_Name == clust_anno]$Cluster

    print(clust_anno)
    print(cluster)

    print("Normalizing")
    seurat_updated_list_sub[[resolution]][[clust_anno]] <- subset(seurat_updated_list[[resolution]], idents = cluster)
    seurat_updated_list_sub[[resolution]][[clust_anno]] <- SCTransform(seurat_updated_list_sub[[clust_anno]], verbose = TRUE)
    seurat_updated_list_sub[[resolution]][[clust_anno]] <- RunPCA(seurat_updated_list_sub[[clust_anno]], npcs = 100)
    seurat_updated_list_sub[[clust_anno]] <- RunUMAP(seurat_updated_list_sub[[clust_anno]], reduction = "pca",dims = 1:100)
    seurat_updated_list_sub[[clust_anno]] <- FindNeighbors(seurat_updated_list_sub[[clust_anno]], reduction = "pca", dims = 1:100)

    ### Cluster multiple resolutions to identify the best that will effectively separate these cell types ###
    print("Subclustering")
    dir.create(paste0(clust_dir,clust_anno, "/"))

    for (resolution in c(seq(0.01,0.09,0.01),seq(0.1,2, by = 0.1))){
        seurat_updated_list[[as.character(resolution)]] <- FindClusters(seurat_updated_list_sub[[clust_anno]], resolution = resolution)

        plot_list[[as.character(resolution)]] <- DimPlot(seurat_updated_list[[as.character(resolution)]], reduction = "umap", label = TRUE, repel = TRUE)

        ggsave(plot_list[[as.character(resolution)]], filename = paste0(clust_dir, clust_anno, "/clusters_resolution_",resolution, ".png"))
    }

    ##### Make figures with markers #####
    ### General Markers ###
    print("Making Marker Expression Plots")
    DefaultAssay(seurat_updated_list_sub[[clust_anno]]) <- "RNA"
    dir.create(paste0(neuron_marker_dir,clust_anno, "/"))

    for (gene in unique(markers$ENSG_ID)){
        gene_id <- markers[ENSG_ID == gene]$Gene
        if (gene %in% rownames(seurat_updated_list_sub[[clust_anno]])){
            if (rowSums(seurat_updated_list_sub[[clust_anno]][gene,]) > 0){
                if(!file.exists(paste0(neuron_marker_dir,clust_anno, "/",gene_id,"_umap_seurat_sub.png"))){
                    print(gene_id)
                    umap <- FeaturePlot(seurat_updated_list_sub[[clust_anno]], features = gene) + labs(title = gene_id)
                    ggsave(umap, filename = paste0(neuron_marker_dir, clust_anno, "/", gene_id,"_umap_seurat_sub.png"))
                }
                if(!file.exists(paste0(neuron_marker_dir,clust_anno, "/", gene_id,"_umap_nebulosa.png"))){
                    density_plot <- plot_density(seurat_updated_list_sub[[clust_anno]], gene, pal = "plasma", reduction = "umap") + labs(title = gene_id)
                    ggsave(density_plot, filename = paste0(neuron_marker_dir, clust_anno, "/", gene_id,"_umap_nebulosa.png"))
                }
            }
        }
    }
}



##### Add in subclustering #####
subclusters <- fread(paste0(outdir,"cluster_annotations_subcluster.tsv"))
subclusters$Cluster <- factor(subclusters$Cluster)

updated_names <- list()

for (clust in unique(subclusters$Original_Cluster)){
    clust_anno <- cluster_anno_dt[Cluster == clust]$Cluster_Name

    DefaultAssay(seurat_updated_list_sub[[clust_anno]]) <- "SCT"

    seurat_updated_list_sub[[clust_anno]] <- FindClusters(seurat_updated_list_sub[[clust_anno]], resolution = unique(subclusters[Original_Cluster == clust]$resolution))

    updated_names[[clust_anno]] <- left_join(data.frame(Idents(seurat_updated_list_sub[[clust_anno]])), subclusters[Original_Cluster == clust, c("Cluster", "Annotation")], by = c("Idents.seurat_sub..clust_anno..." = "Cluster"))
    updated_names[[clust_anno]]$Idents.seurat_updated_list_sub..clust_anno... <- NULL

    rownames(updated_names[[clust_anno]]) <- colnames(seurat_updated_list_sub[[clust_anno]])
}

updated_names_df <- do.call(rbind, updated_names)
rownames(updated_names_df) <- gsub(".+\\..+\\.", "", rownames(updated_names_df))
updated_names_df$Barcode <- rownames(updated_names_df)


### Add Cell Types to metadata ###
cluster_anno_dt$Cluster <- as.factor(cluster_anno_dt$Cluster)
seurat_idents <- data.frame(Idents(seurat_updated_list[[as.character(resolution)]]))
colnames(seurat_idents) <- "Cluster"
cell_type_dt <- left_join(seurat_idents, cluster_anno_dt, by = c("Cluster"))
rownames(cell_type_dt) <- colnames(seurat_updated_list[[as.character(resolution)]])
cell_type_dt$Barcode <- rownames(cell_type_dt)

combined_updated_cell_type <- left_join(cell_type_dt, updated_names_df, by = "Barcode")
combined_updated_cell_type$Cluster_Name <- ifelse(!is.na(combined_updated_cell_type$Annotation), combined_updated_cell_type$Annotation, combined_updated_cell_type$Cluster_Name)
combined_updated_cell_type$Cluster_Name  <- gsub(" ", "_",combined_updated_cell_type$Cluster_Name )
rownames(combined_updated_cell_type) <- combined_updated_cell_type$Barcode

combined_updated_cell_type$Barcode <- NULL
combined_updated_cell_type$Annotation <- NULL
colnames(combined_updated_cell_type)[1] <- c("Cluster")

rownames(combined_updated_cell_type) <- rownames(cell_type_dt)

seurat_annotated <- AddMetaData(seurat_updated_list[[as.character(resolution)]], combined_updated_cell_type)


saveRDS(seurat_annotated, paste0(outdir,"seurat_neurons_annotated.rds"))

