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
seurat <- readRDS(paste0(outdir,"seurat_resolution_1.7.rds"))
cluster_anno_dt <- fread(paste0(outdir,"cluster_annotations.tsv"), sep = "\t")
cluster_anno_dt$Cluster_Name <- gsub(" ", "_", cluster_anno_dt$Cluster_Name) %>% gsub("/", ".", .)


##### Subcluster clusters with "/" in them (because likely two different cell types in cluster) #####
seurat_sub <- list()

seurat_updated_list <- list()
plot_list <- list()


clust_dir <- paste0(outdir,"Annotated_Subclustering/Clusters/")
dir.create(clust_dir, recursive = TRUE)

markers <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/Parkinsons_TH_Collab/data/genes_of_interest.tsv")

neuron_marker_dir <- paste0(outdir, "Annotated_Subclustering/Markers/")
dir.create(neuron_marker_dir)


for (clust_anno in grep("\\.", cluster_anno_dt$Cluster_Name, value = TRUE)){
    cluster <- cluster_anno_dt[Cluster_Name == clust_anno]$Cluster

    print(clust_anno)
    print(cluster)

    print("Normalizing")
    seurat_sub[[clust_anno]] <- subset(seurat, idents = cluster)
    seurat_sub[[clust_anno]] <- SCTransform(seurat_sub[[clust_anno]], verbose = TRUE)
    seurat_sub[[clust_anno]] <- RunPCA(seurat_sub[[clust_anno]], npcs = 100)
    seurat_sub[[clust_anno]] <- RunUMAP(seurat_sub[[clust_anno]], reduction = "pca",dims = 1:100)
    seurat_sub[[clust_anno]] <- FindNeighbors(seurat_sub[[clust_anno]], reduction = "pca", dims = 1:100)

    ### Try multiple clustering categories
    print("Subclustering")
    dir.create(paste0(clust_dir,clust_anno, "/"))

    for (resolution in c(seq(0.01,0.09,0.01),seq(0.1,2, by = 0.1))){
        seurat_updated_list[[as.character(resolution)]] <- FindClusters(seurat_sub[[clust_anno]], resolution = resolution)

        plot_list[[as.character(resolution)]] <- DimPlot(seurat_updated_list[[as.character(resolution)]], reduction = "umap", label = TRUE, repel = TRUE)

        ggsave(plot_list[[as.character(resolution)]], filename = paste0(clust_dir, clust_anno, "/clusters_resolution_",resolution, ".png"))
    }

    ##### Make figures with markers #####
    ### General Markers ###
    print("Making Marker Expression Plots")
    DefaultAssay(seurat_sub[[clust_anno]]) <- "RNA"
    dir.create(paste0(neuron_marker_dir,clust_anno, "/"))

    for (gene in unique(markers$ENSG_ID)){
        gene_id <- markers[ENSG_ID == gene]$Gene
        if (gene %in% rownames(seurat_sub[[clust_anno]])){
            if (rowSums(seurat_sub[[clust_anno]][gene,]) > 0){
                if(!file.exists(paste0(neuron_marker_dir,clust_anno, "/",gene_id,"_umap_seurat_sub.png"))){
                    print(gene_id)
                    umap <- FeaturePlot(seurat_sub[[clust_anno]], features = gene) + labs(title = gene_id)
                    ggsave(umap, filename = paste0(neuron_marker_dir, clust_anno, "/", gene_id,"_umap_seurat_sub.png"))
                }
                if(!file.exists(paste0(neuron_marker_dir,clust_anno, "/", gene_id,"_umap_nebulosa.png"))){
                    density_plot <- plot_density(seurat_sub[[clust_anno]], gene, pal = "plasma", reduction = "umap") + labs(title = gene_id)
                    ggsave(density_plot, filename = paste0(neuron_marker_dir, clust_anno, "/", gene_id,"_umap_nebulosa.png"))
                }
            }
        }
    }
}


##### Add in Lily subclustering #####
subclusters <- fread(paste0(outdir,"lily_cluster_annotations_subcluster.tsv"))
subclusters$Cluster <- factor(subclusters$Cluster)

updated_names <- list()

for (clust in unique(subclusters$Original_Cluster)){
    clust_anno <- cluster_anno_dt[Cluster == clust]$Cluster_Name

    DefaultAssay(seurat_sub[[clust_anno]]) <- "SCT"

    seurat_sub[[clust_anno]] <- FindClusters(seurat_sub[[clust_anno]], resolution = unique(subclusters[Original_Cluster == clust]$resolution))

    updated_names[[clust_anno]] <- left_join(data.frame(Idents(seurat_sub[[clust_anno]])), subclusters[Original_Cluster == clust, c("Cluster", "Annotation")], by = c("Idents.seurat_sub..clust_anno..." = "Cluster"))
    updated_names[[clust_anno]]$Idents.seurat_sub..clust_anno... <- NULL

    rownames(updated_names[[clust_anno]]) <- colnames(seurat_sub[[clust_anno]])
}

updated_names_df <- do.call(rbind, updated_names)
rownames(updated_names_df) <- gsub(".+\\..+\\.", "", rownames(updated_names_df))
updated_names_df$Barcode <- rownames(updated_names_df)


### Add Cell Types to metadata ###
cluster_anno_dt$Cluster <- as.factor(cluster_anno_dt$Cluster)
cell_type_dt <- left_join(data.frame(Idents(seurat4train)), cluster_anno_dt, by = c("Idents.seurat4train." = "Cluster"))
rownames(cell_type_dt) <- colnames(seurat4train)
cell_type_dt$Barcode <- rownames(cell_type_dt)

combined_updated_cell_type <- left_join(cell_type_dt, updated_names_df, by = "Barcode")
combined_updated_cell_type$Cluster_Name <- ifelse(!is.na(combined_updated_cell_type$Annotation), combined_updated_cell_type$Annotation, combined_updated_cell_type$Cluster_Name)
combined_updated_cell_type$Cluster_Name  <- gsub(" ", "_",combined_updated_cell_type$Cluster_Name )
rownames(combined_updated_cell_type) <- combined_updated_cell_type$Barcode

combined_updated_cell_type$Barcode <- NULL
combined_updated_cell_type$Annotation <- NULL
colnames(combined_updated_cell_type)[1] <- c("Cluster")

rownames(combined_updated_cell_type) <- rownames(cell_type_dt)

seurat <- AddMetaData(seurat, combined_updated_cell_type)



##### Hack scPred train + classification to try and reclassify droplets that don't annotate well - possibly because clustered together but actually shouldn't be same cell type #####
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
seurat4train <- trainModel (seurat4train, model = "mda")

### Look at the probabilities ###
get_probabilities(seurat4train) %>% head()

get_scpred(seurat4train)

pTop <- plot_probabilities(seurat4train)
ggsave(pTop, filename = paste0(train_test_dir, "top_classification_probabilities.png"))

get_probabilities(seurat4train)[which(seurat4train@meta.data$Cluster_Name == "OXT_Neurons" & get_probabilities(seurat4train)$OXT_Neurons < 0.5),]


max_dt <- data.table(Barcode = colnames(seurat4train),
						max_prob = apply(get_probabilities(seurat4train), 1, max), 
						cell_type = colnames(get_probabilities(seurat4train))[apply(get_probabilities(seurat4train), 1, which.max)],
						original_cell_type = seurat4train@meta.data$Cluster_Name)

max_dt$updated_names <- ifelse(max_dt$max_prob < 0.5, "Unassigned", max_dt$cell_type)

table(max_dt[max_dt$max_prob < 0.5]$cell_type, max_dt[max_dt$max_prob < 0.5]$original_cell_type)

table(max_dt[max_prob < 0.5]$original_cell_type)

updated_df <- data.frame(max_dt[,c("updated_names")])
rownames(updated_df) <- max_dt$Barcode

seurat4train <- AddMetaData(seurat4train, updated_df)


### Retest ###
seurat4train2 <- getFeatureSpace(seurat4train, "updated_names")

### Train the models ###
seurat4train2 <- trainModel (seurat4train2, model = "mda")

### Look at the probabilities ###
get_probabilities(seurat4train2) %>% head()

get_scpred(seurat4train)
get_scpred(seurat4train2)

pTop <- plot_probabilities(seurat4train2)
ggsave(pTop, filename = paste0(train_test_dir, "top_classification_probabilities_round2.png"))



### Retrain ###
max_dt2 <- data.table(Barcode = colnames(seurat4train2),
						max_prob = apply(get_probabilities(seurat4train2), 1, max), 
						cell_type = colnames(get_probabilities(seurat4train2))[apply(get_probabilities(seurat4train2), 1, which.max)],
						original_cell_type = seurat4train2@meta.data$Cluster_Name)

max_dt2$updated_names <- ifelse(max_dt2$max_prob < 0.5, "Unassigned", max_dt2$cell_type)

table(max_dt2[max_dt2$max_prob < 0.5]$cell_type, max_dt2[max_dt2$max_prob < 0.5]$original_cell_type)

table(max_dt2[max_prob < 0.5]$original_cell_type)

updated_df2 <- data.frame(max_dt2[,c("updated_names")])
rownames(updated_df2) <- max_dt2$Barcode

seurat4train2 <- AddMetaData(seurat4train2, updated_df2)


### Retest ###
seurat4train3 <- getFeatureSpace(seurat4train2, "updated_names")

### Train the models ###
seurat4train3 <- trainModel (seurat4train3, model = "mda")

### Look at the probabilities ###
get_probabilities(seurat4train3) %>% head()

get_scpred(seurat4train)
get_scpred(seurat4train2)
get_scpred(seurat4train3)

pTop <- plot_probabilities(seurat4train3)
ggsave(pTop, filename = paste0(train_test_dir, "top_classification_probabilities_round3.png"))


### Not a huge improvement and don't want to overfit - will stick with one round of train-test-reclassify ###
seurat_annotated <- AddMetaData(seurat4train, updated_df)
saveRDS(seurat_annotated, paste0(train_test_dir, "neurons_classified.rds"))
seurat_annotated <- readRDS(paste0(train_test_dir, "neurons_classified.rds"))


seurat_annotated <- SCTransform(seurat_annotated, verbose = TRUE)
seurat_annotated <- RunPCA(seurat_annotated, npcs = 100)
seurat_annotated <- RunUMAP(seurat_annotated, reduction = "pca",dims = 1:100)
seurat_annotated <- FindNeighbors(seurat_annotated, reduction = "pca", dims = 1:100)

umap_ref_reclass <- DimPlot(seurat_annotated, group.by = "updated_names", label = TRUE, repel = TRUE) +
						scale_color_manual(values = c("#ef1a1e", "#ee0d80", "#ba51a0", "#953d99", "#3653a7", "#3e79bd", "#6abe43", "#a7cf37",  "#fcd502", "#faa412", "#f4801a", "#f48667", "#f7b0c7", "#d095c2", "#cbabd0", "#8187c3", "#afbde3", "#a7d48a", "#c8e18b", "#ffecac", "#ffc67a", "#fbb172"))
ggsave(umap_ref_reclass, filename = paste0(train_test_dir,"umap_reclassified.png"), width = 10)


umap_ref_reclass_no_lab <- DimPlot(seurat_annotated, group.by = "updated_names", label = FALSE, repel = TRUE) +
						scale_color_manual(values = c("#ef1a1e", "#ee0d80", "#ba51a0", "#953d99", "#3653a7", "#3e79bd", "#6abe43", "#a7cf37",  "#fcd502", "#faa412", "#f4801a", "#f48667", "#f7b0c7", "#d095c2", "#cbabd0", "#8187c3", "#afbde3", "#a7d48a", "#c8e18b", "#ffecac", "#ffc67a", "#fbb172"))
ggsave(umap_ref_reclass_no_lab, filename = paste0(train_test_dir,"umap_reclassified_no_label.png"), width = 10)





umap_orig <- DimPlot(seurat_annotated, group.by = "Cluster_Name", label = TRUE, repel = TRUE) +
						scale_color_manual(values = c("#ef1a1e", "#ee0d80", "#ba51a0", "#953d99", "#3653a7", "#3e79bd", "#6abe43", "#a7cf37",  "#fcd502", "#faa412", "#f4801a", "#f48667", "#f7b0c7", "#d095c2", "#cbabd0", "#8187c3", "#afbde3", "#a7d48a", "#c8e18b", "#ffecac", "#ffc67a", "#fbb172"))
ggsave(umap_orig, filename = paste0(train_test_dir,"umap_clust_label.png"), width = 10)


umap_orig_no_lab <- DimPlot(seurat_annotated, group.by = "Cluster_Name", label = FALSE, repel = TRUE) +
						scale_color_manual(values = c("#ef1a1e", "#ee0d80", "#ba51a0", "#953d99", "#3653a7", "#3e79bd", "#6abe43", "#a7cf37",  "#fcd502", "#faa412", "#f4801a", "#f48667", "#f7b0c7", "#d095c2", "#cbabd0", "#8187c3", "#afbde3", "#a7d48a", "#c8e18b", "#ffecac", "#ffc67a", "#fbb172"))
ggsave(umap_orig_no_lab, filename = paste0(train_test_dir,"umap_clust_label_no_label.png"), width = 10)






#### If not unassigned first, then add back original label
seurat_annotated@meta.data$updated_names <- ifelse(seurat_annotated@meta.data$updated_names == "Unassigned", seurat_annotated@meta.data$Cluster_Name, seurat_annotated@meta.data$updated_names)




umap_ref_reclass_assign <- DimPlot(subset(seurat_annotated, subset = updated_names != "Unassigned"), group.by = "updated_names", label = TRUE, repel = TRUE)
ggsave(umap_ref_reclass_assign, filename = paste0(train_test_dir,"umap_reclassified_no_unassigned.png"), width = 10)
 

saveRDS(seurat_annotated, paste0(train_test_dir, "neurons_classified_sct_map.rds"))


##### Make popout clusters for Lily #####
highlight_dir <- paste0(train_test_dir, "Highlight_clusters/")
dir.create(highlight_dir, recursive = TRUE)

Idents(seurat_annotated) <- seurat_annotated$updated_names
plot_highlight_list <- list()

for (ident in unique(seurat_annotated$updated_names)){
	plot_highlight_list[[as.character(ident)]] <- DimPlot(seurat_annotated, reduction = "umap", cells.highlight = WhichCells(seurat_annotated, idents = ident), cols.highlight = c("#8B0269")) + ggtitle(as.character(ident))
	ggsave(plot_highlight_list[[as.character(ident)]], filename = paste0(highlight_dir,as.character(ident),"_integrated.png"), width = 8)
}



##### Update with new names from Lily (combining some groups and renaming some) #####
seurat_annotated <- readRDS(paste0(train_test_dir, "neurons_classified_sct_map.rds"))


### Update names with Lily new IDs ###
seurat_annotated@meta.data$Broad_Neurons <- gsub("_", " ", ifelse(seurat_annotated@meta.data$updated_names == "KISS1_Neurons", "DNs", 
    ifelse(seurat_annotated@meta.data$updated_names == "Midbrain", "Midbrain Neurons", 
        ifelse(seurat_annotated@meta.data$updated_names == "OTP+_DNs", "DNs", 
            ifelse(seurat_annotated@meta.data$updated_names == "Pituitary_Cells", "Neuroendocrine Cells", 
                ifelse(seurat_annotated@meta.data$updated_names == "RAX+_DNs", "DNs", 
                    ifelse(seurat_annotated@meta.data$updated_names == "Telencephalon", "Preoptic Area Neurons", 
                        ifelse(seurat_annotated@meta.data$updated_names == "Unassigned", "Unknown Neuron Subtype", 
                            ifelse(seurat_annotated@meta.data$updated_names == "Mesenchymal_Cells" , "Neuroendocrine_Cells", seurat_annotated@meta.data$updated_names)))))))))


saveRDS(seurat_annotated, paste0(train_test_dir, "neurons_classified_updated_names_sct_map.rds"))
seurat_annotated <- readRDS(paste0(train_test_dir, "neurons_classified_updated_names_sct_map.rds"))

colors_dt <- fread(paste0(dir, "data/cell_type_colors.tsv"), sep = "\t")

colors_vec <- colors_dt$Color
names(colors_vec) <- colors_dt$Group
colors_vec <- colors_vec[names(colors_vec) != "Mesenchymal Cells"]


umap_orig <- DimPlot(seurat_annotated, group.by = "Broad_Neurons", label = TRUE, repel = TRUE) +
						scale_color_manual(values = colors_vec)
ggsave(umap_orig, filename = paste0(train_test_dir,"umap_clust_label_updated_", Sys.Date(), ".png"), width = 10)
ggsave(umap_orig, filename = paste0(train_test_dir,"umap_clust_label_updated_", Sys.Date(), ".pdf"), width = 10)


umap_orig_no_lab <- DimPlot(seurat_annotated, group.by = "Broad_Neurons", label = FALSE, repel = TRUE) +
						scale_color_manual(values = colors_vec) 
ggsave(umap_orig_no_lab, filename = paste0(train_test_dir,"umap_clust_label_no_label_updated_", Sys.Date(), ".png"), width = 10)
ggsave(umap_orig_no_lab, filename = paste0(train_test_dir,"umap_clust_label_no_label_updated_", Sys.Date(), ".pdf"), width = 10)





##### Test for DEG for each cell type #####
Idents(seurat_annotated) <- seurat_annotated@meta.data$Broad_Neurons

CellType_Markers <- list()

for (ct in unique(seurat_annotated@meta.data$Broad_Neurons)){
	CellType_Markers[[ct]] <- FindMarkers(seurat_annotated, ident.1 = ct,ident.2 = NULL, only.pos = TRUE)
	CellType_Markers[[ct]]$Gene <- rownames(CellType_Markers[[ct]])
}



### Add in Gene IDs to tables ###
gene_conversions <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/Parkinsons_TH_Collab/data/Remapping/remapped/1205_GEX_WT_THpos_GEX_0_1_HC3GJDSXY/outs/filtered_feature_bc_matrix/features.tsv.gz", header = FALSE)
gene_conversions$V3 <- NULL
colnames(gene_conversions) <- c("Gene", "Gene_ID")


CellType_Markers <- lapply(CellType_Markers, function(x) {
    gene_conversions[x, on = "Gene"]
})


lapply(CellType_Markers, head)


lapply(names(CellType_Markers), function(celltype){
    fwrite(CellType_Markers[[celltype]], gsub(" ", "_", paste0(outdir,celltype,"_neuron_subtypes.tsv")), sep = "\t")
})



##### Make pop out figures of neurons for Lily for paper #####
##### Figures of each cell type poppe4d #####
colors_dt <- fread(paste0(dir, "data/cell_type_colors.tsv"), sep = "\t")

neuron_colors <- colors_dt$Color
names(neuron_colors) <- colors_dt$Group
neuron_colors <- neuron_colors[names(neuron_colors) != "Mesenchymal Cells"]


cells_neurons <- list()
plot_neurons <- list()

Idents(seurat_annotated) <- "Broad_Neurons"

for (type in unique(seurat_annotated@meta.data$Broad_Neurons)[!unique(seurat_annotated@meta.data$Broad_Neurons) %in% NA]){
	cells_neurons[[type]] <- WhichCells(seurat_annotated, idents = type)

	plot_neurons[[type]] <- DimPlot(seurat_annotated, reduction = "umap",cells.highlight = list(cells_neurons[[type]]), cols.highlight = neuron_colors[type]) + 
                        ggtitle(type) 
	ggsave(plot_neurons[[type]], filename = paste0(outdir,type,"_umap.png"), width = 8)
}


