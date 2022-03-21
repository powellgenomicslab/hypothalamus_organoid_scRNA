##### Intro #####
### Author: Drew Neavin
### Date: 15 August, 2021
### Reason: Combine various QC metrics and create additional metrics after remapping with tdtomato

##### Load in Libraries #####
library(data.table)
library(Seurat)
library(ggplot2)
library(Nebulosa)
library(tidyverse)

##### Set up directories
dir <- "/path/to/base/directory/"
outdir <- paste0(dir,"output/QC/seurat_QC/")
dir.create(outdir, recursive = TRUE)



##### Read in files #####
### seurat objects ###
seurat_list <- readRDS(paste0(dir,"output/QC/Hash_Demultiplexing/seurat_list.rds"))

### make a variable that houses the list of pools
pools <- names(seurat_list)

### scds and dropletqc ###
scds_dropqc <- fread(paste0(dir,"output/QC/DropletQC_scds/scds_dropletQC_results.tsv"))

## Remove the -1 from the end of each barcode and add Barcode as rowname ##
scds_dropqc <- data.frame(scds_dropqc)
rownames(scds_dropqc) <- gsub("-1","",scds_dropqc$Barcode)


### Combine the seurat objects together
seurat <- merge(seurat_list[[1]], y = seurat_list[[2]])


## Add the scds results to the seurat objects
seurat <- AddMetaData(seurat, scds_dropqc)



### Classify final doublet calls from scds, DropletQC and hashtag demultiplexing
## scds doublet => doublet
## MULTI_ID doublet => doublet
## DropletQC cell_status empty_droplet => empty_droplet
seurat@meta.data$combined_classification <- ifelse(seurat@meta.data$MULTI_ID == "Doublet" | seurat@meta.data$scds_DropletType == "doublet", "doublet", ifelse(seurat@meta.data$cell_status == "empty_droplet", "empty_droplet", seurat@meta.data$MULTI_ID))


### Add QC features (mt and rb %)
RbGeneList <- read.delim(file = "RibosomalGeneList_GeneID_ENSG.txt") ## files provided in github rb_mt_gene_list directory
MtGeneList <- read.delim(file = "MtGeneList_GeneID_ENSG.txt") ## files provided in github rb_mt_gene_list directory

print("Calculating Mt %")
seurat <- PercentageFeatureSet(seurat, features = MtGeneList$ENSG, col.name = "percent.mt")

print("Calculating Rb %")
RbGeneList <- RbGeneList[which(RbGeneList$ENSG %in% rownames(seurat)),]
seurat <- PercentageFeatureSet(seurat, features = RbGeneList$ENSG, col.name = "percent.rb")


### Add tdtomato and TH together since tdtomato is attached to TH on one allele
seurat@meta.data$th_tdtomato <- colSums(seurat[["RNA"]]@counts[c("ENSG00000180176","tdtomato"),])


saveRDS(seurat, paste0(outdir,"seurat_all_cells.rds"))


print("Adding Cell Cycle")
cell_cycle <- fread(paste0(dir,"output/QC/cyclone_cell_cycle/CellCycleProportions.txt"), sep = "\t")
rownames(cell_cycle) <- cell_cycle$V1
cell_cycle$V1 <- NULL
seurat <- AddMetaData(seurat, cell_cycle)


##### Remove the doublets and negative assignments
seurat_sub <- subset(seurat, subset = combined_classification != "Negative")
seurat_sub <- subset(seurat_sub, subset = combined_classification != "doublet")
seurat_sub <- subset(seurat_sub, subset = cell_status == "cell")

##### Remove cells with more than 12.5% mt genes and 
seurat_sub <- subset(seurat_sub, subset = percent.mt <= 12.5)



##### Normalize and Plot results #####
### Make pre-QC figures ###
seurat_sub <- SCTransform(seurat_sub, method = "glmGamPoi", verbose = FALSE) 
seurat_sub <- RunPCA(seurat_sub, features = VariableFeatures(object = seurat_sub))
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:10)
seurat_sub <- FindClusters(seurat_sub, resolution = 0.7)
seurat_sub <- RunUMAP(seurat_sub, dims = 1:10, seed.use = 31, n.neighbors = 25)


##### Take a preliminary look at the clusters and features definining them
umap <- DimPlot(seurat_sub)
ggsave(umap, filename = paste0(outdir,"umap.png"))


saveRDS(seurat_sub, paste0(outdir,"seurat_singlets.rds"))
seurat_sub <- readRDS(paste0(outdir,"seurat_singlets.rds"))


outdir_filtered <- paste0(outdir,"singlets/")
dir.create(outdir_filtered)

### QC Figures ###
plot_mt_pct <- VlnPlot(seurat_sub, features = c( "percent.mt"), group.by = "combined_classification", pt.size = 0)

plot_rb_pct <- VlnPlot(seurat_sub, features = c( "percent.rb"), group.by = "combined_classification", pt.size = 0)

plot_n_count <- VlnPlot(seurat_sub, features = c( "nCount_RNA"), group.by = "combined_classification", pt.size = 0)

plot_nFeature_RNA <- VlnPlot(seurat_sub, features = c( "nFeature_RNA"), group.by = "combined_classification", pt.size = 0)

lib_mt <- FeatureScatter(seurat_sub, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "combined_classification")

lib_genes <- FeatureScatter(seurat_sub, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "combined_classification")

nebulosa_mt_umap <- plot_density(seurat_sub, "percent.mt", pal = "plasma")

nebulosa_rb_umap <- plot_density(seurat_sub, "percent.rb", pal = "plasma")

nebulosa_nFeature_RNA_umap <- plot_density(seurat_sub, "nFeature_RNA", pal = "plasma")

nebulosa_nCount_RNA_umap <- plot_density(seurat_sub, "nCount_RNA", pal = "plasma")

nFeature_RNA_umap <- FeaturePlot(seurat_sub, features = "nFeature_RNA")

nCount_RNA_umap <- FeaturePlot(seurat_sub, features = "nCount_RNA")

umap_location <- DimPlot(seurat_sub, group.by = "combined_classification")

umap_location_split <- DimPlot(seurat_sub, group.by = "combined_classification", split.by = "combined_classification",  ncol = 3)

mt_umap <- FeaturePlot(seurat_sub, features = "percent.mt")

rb_umap <- FeaturePlot(seurat_sub, features = "percent.rb")

umap_condition <- DimPlot(seurat_sub, group.by = "pool")

umap_condition_split <- DimPlot(seurat_sub, group.by = "pool", split.by = "pool",  ncol = 3)

umap_cell_cycle <- DimPlot(seurat_sub, group.by = "phases")

### Make a reference from ENSG IDs to Gene IDs ###
GeneConversion <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/Parkinsons_TH_Collab/data/GEX/GEX_WT_TH/features.tsv.gz", header = F, sep = "\t")
GeneConversion$V3 <- NULL
colnames(GeneConversion) <- c("ENSG_ID", "Gene_ID")



##### Check for features of interest #####

### Make a reference from ENSG IDs to Gene IDs ###
GeneConversion <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/Parkinsons_TH_Collab/data/GEX/GEX_WT_TH/features.tsv.gz", header = F, sep = "\t")
GeneConversion$V3 <- NULL
colnames(GeneConversion) <- c("ENSG_ID", "Gene_ID")



##### Check for features of interest #####
markers <- fread("genes_of_interest.tsv") ### List of markers tested can be found in supplementar tables for this manuscript

genes_dir <- paste0(outdir_filtered, "genes_of_interest/")
dir.create(genes_dir)


for (gene in unique(markers$ENSG)){
	if (gene %in% rownames(seurat_sub)){
		gene_name <- GeneConversion[ENSG_ID == gene]$Gene_ID
		if(!file.exists(paste0(genes_dir,gene_name,"_umap.png"))){
			plot <- FeaturePlot(seurat_sub, features = gene) + labs(title = gene_name)
			ggsave(plot, filename = paste0(genes_dir,gene_name,"_umap.png"))
		}
		if(!file.exists(paste0(genes_dir,gene_name,"_umap_nebulosa.png"))){
			density_plot <- plot_density(seurat_sub, gene, pal = "plasma") + labs(title = gene_name)
			ggsave(density_plot, filename = paste0(genes_dir,gene_name,"_umap_nebulosa.png"))
		}
	}
}


th_tdtomato <- FeaturePlot(object = seurat_sub,
	order = TRUE,
    features = c("ENSG00000180176", "tdtomato"),
    cols = c("grey90", "#FF0033", "blue"),
	blend = TRUE,
	blend.threshold = 0)
ggsave(th_tdtomato, filename = paste0(genes_dir,"TH_tdtomato.png"), height = 5, width = 15)

