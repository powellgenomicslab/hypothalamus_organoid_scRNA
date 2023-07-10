library(data.table)
library(Seurat)
library(scPred)
library(tidyverse)
library(ggridges)
library(Nebulosa)



##### Set Up Directories #####
dir <- "/path/to/base/directory/"
outdir <- paste0(dir,"output/Classification/Neuron_Clustering/")



##### Read in Files #####
seurat <- readRDS(paste0(outdir,"seurat_neurons_annotated.rds"))
cluster_anno_dt$Cluster_Name <- gsub(" ", "_", cluster_anno_dt$Cluster_Name) %>% gsub("/", ".", .)





##### Use bootstrapping to reclassify droplets that don't annotate well - possibly because clustered together but actually shouldn't be same cell type #####
train_test_dir <- paste0(outdir, "Train_Test_Reannotation/")
dir.create(train_test_dir)

DefaultAssay(seurat) <- "RNA"


### Process the data ###
seurat4train <- seurat %>% 
	NormalizeData() %>%
	FindVariableFeatures() %>%
	ScaleData() %>%
	RunPCA() %>%
	RunUMAP(dims = 1:30)



seurat4train <- AddMetaData(seurat4train, combined_updated_cell_type)

seurat4train <- getFeatureSpace(seurat4train, "Cluster_Name")

### Train the models ###
seurat4train <- trainModel(seurat4train, model = "mda")

### Look at the probabilities ###
get_probabilities(seurat4train) %>% head()

get_scpred(seurat4train)

pTop <- plot_probabilities(seurat4train)

get_probabilities(seurat4train)[which(seurat4train@meta.data$Cluster_Name == "OXT_Neurons" & get_probabilities(seurat4train)$OXT_Neurons < 0.5),]


max_dt <- data.table(Barcode = colnames(seurat4train),
						max_prob = apply(get_probabilities(seurat4train), 1, max), 
						cell_type = colnames(get_probabilities(seurat4train))[apply(get_probabilities(seurat4train), 1, which.max)],
						original_cell_type = seurat4train@meta.data$Cluster_Name,
                        second_prob = apply(get_probabilities(seurat4train), 1, Rfast::nth,k = 2, descending = T),
                        second_assignment = gsub("scpred_","",colnames(get_probabilities(seurat4train))[apply(get_probabilities(seurat4train), 1, Rfast::nth,k = 2, descending = T, index.return = T)]))



max_dt$diff1_2 <- max_dt$max_prob - max_dt$second_prob


diff_plot <- ggplot(max_dt, aes(diff1_2)) +
    geom_density()




max_dt$updated_names <- ifelse(max_dt$diff1_2 < 0.1, "Unclassifiable", ifelse(max_dt$max_prob < 0.5, "Unassigned", max_dt$cell_type))

updated_df <- data.frame(max_dt[,c("updated_names")])
rownames(updated_df) <- max_dt$Barcode

seurat4train <- AddMetaData(seurat4train, updated_df)


### Retest ###
seurat4train2 <- getFeatureSpace(seurat4train, "updated_names")

### Train the models ###
seurat4train2 <- trainModel(seurat4train2, model = "mda")

### Look at the probabilities ###
get_probabilities(seurat4train2) %>% head()

get_scpred(seurat4train)
get_scpred(seurat4train2)

pTop <- plot_probabilities(seurat4train2)



seurat_annotated <- AddMetaData(seurat4train, updated_df)
seurat_annotated@meta.data$updated_names <- ifelse(seurat_annotated@meta.data$updated_names == "Unknown_Neuron_Subtype", seurat_annotated@meta.data$Cluster_Name, seurat_annotated@meta.data$updated_names)
seurat_annotated@meta.data$updated_names <- gsub("Unassigned", "Unknown_Neuron_Subtype", seurat_annotated@meta.data$updated_names)



saveRDS(seurat_annotated, paste0(train_test_dir, "neurons_classified.rds"))

## Remove the cells that were high probabilities of multiple different cell types
seurat_annotated_subset <- subset(seurat_annotated, subset = updated_names != "Unclassifiable")




## Normalize data ##
seurat_annotated_subset <- SCTransform(seurat_annotated_subset, verbose = TRUE)
seurat_annotated_subset <- RunPCA(seurat_annotated_subset, npcs = 100)
seurat_annotated_subset <- RunUMAP(seurat_annotated_subset, reduction = "pca",dims = 1:100)
seurat_annotated_subset <- FindNeighbors(seurat_annotated_subset, reduction = "pca", dims = 1:100)

saveRDS(seurat_annotated_subset, paste0(train_test_dir, "neurons_classified_remove_unclassifiable.rds"))


## Visualize results ##
umap_ref_reclass <- DimPlot(seurat_annotated_subset, group.by = "updated_names", label = TRUE, repel = TRUE) 

umap_ref_reclass_no_lab <- DimPlot(seurat_annotated_subset, group.by = "updated_names", label = FALSE, repel = TRUE) 

umap_orig <- DimPlot(seurat_annotated_subset, group.by = "Cluster_Name", label = TRUE, repel = TRUE) 

umap_orig_no_lab <- DimPlot(seurat_annotated_subset, group.by = "Cluster_Name", label = FALSE, repel = TRUE) 

umap_ref_reclass_assign <- DimPlot(subset(seurat_annotated, subset = updated_names != "Unassigned"), group.by = "updated_names", label = TRUE, repel = TRUE)
 





##### Make popout clusters for better visualization #####
highlight_dir <- paste0(train_test_dir, "Highlight_clusters/")
dir.create(highlight_dir, recursive = TRUE)

Idents(seurat_annotated) <- seurat_annotated$updated_names
plot_highlight_list <- list()

for (ident in unique(seurat_annotated$updated_names)){
	plot_highlight_list[[as.character(ident)]] <- DimPlot(seurat_annotated, reduction = "umap", cells.highlight = WhichCells(seurat_annotated, idents = ident), cols.highlight = c("#8B0269")) + ggtitle(as.character(ident))
}



##### Update with new names from Lily (combining some groups and renaming some) #####
seurat_annotated <- readRDS(paste0(train_test_dir, "neurons_classified_sct_map.rds"))


### Update names with Lily new IDs ###
seurat_annotated@meta.data$Broad_Neurons <- gsub("_", " ", seurat_annotated@meta.data$updated_names)


saveRDS(seurat_annotated, paste0(train_test_dir, "neurons_classified_updated_names_sct_map.rds"))



## Visualize ##
umap_orig <- DimPlot(seurat_annotated, group.by = "Broad_Neurons", label = TRUE, repel = TRUE) 

umap_orig_no_lab <- DimPlot(seurat_annotated, group.by = "Broad_Neurons", label = FALSE, repel = TRUE) 



