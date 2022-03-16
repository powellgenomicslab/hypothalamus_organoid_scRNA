##### Intro #####
### Author: Drew Neavin
### Date: 28 June, 2021
### Reason: Preprocess Neonatal hypothalamus snRNA-seq dataset
### Original dataset from Huang et al, "Generation of hypothalamic arcuate organoids from human induced pluripotent stem cells"
### Data downloaded from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164101

##### Load in Libraries #####
library(data.table)
library(Seurat)
library(tidyverse)
library(Nebulosa)


##### Set up directories #####
dir <- "/path/to/base/directory/"
datadir <- paste0(dir,"data/Human_Neonatal_Hypothalamus/")
outdir <- paste0(dir,"output/Classification/Neonatal_snRNAseq/Preprocessing/")
dir.create(outdir, recursive = TRUE)


##### Read in the data #####
file_list <- list.files(datadir, pattern = ".gz")

counts_list <- lapply(file_list, function(x){
	fread(paste0(datadir,x))
})
names(counts_list) <- gsub("_DGE.+", "", file_list)

### Add sample name to barcodes to help make unique in case of same barcode in different pools ###
counts_list <- lapply(names(counts_list), function(x){
	colnames(counts_list[[x]])[2:ncol(counts_list[[x]])] <- paste0(x, "_", colnames(counts_list[[x]])[2:ncol(counts_list[[x]])])
	return(counts_list[[x]])
})

##### Merge the data together before loading into Seurat (because different number of genes causes merging problems after) #####
mymerge <- function(x, y) merge(x, y, on = "GENE", all = TRUE)
counts_combined <- Reduce(mymerge,counts_list)

### Prepare dataframe to make in to seurat object ###
genes <- counts_combined$GENE
counts_combined$GENE <- NULL

### Replace NA with 0 ###
counts_combined[is.na(counts_combined)] <- 0
counts_combined <- data.frame(counts_combined)
rownames(counts_combined) <- genes


### Use provided samples.txt to add metadata to seurat object
samples <- fread(paste0(datadir,"samples.txt"), header = FALSE)
colnames(samples) <- c("GSM", "Description")
samples$Patient <- gsub("Human hypothalamus cells_patient", "", samples$Description) %>% gsub("_run\\d+", "", .)



###### Make a Seurat object #####
seurat <- CreateSeuratObject(data.frame(counts_combined))


seurat@meta.data$GSM <- gsub("_[ACTG]+","",rownames(seurat@meta.data)) %>% gsub("N", "", .)
seurat@meta.data$Individual <- samples[data.table(GSM = seurat@meta.data[,c("GSM")]), on = "GSM"]$Patient



### Add QC features (mt and rb %)
RbGeneList <- read.delim(file = "RibosomalGeneList_GeneID_ENSG.txt") ## files provided in github rb_mt_gene_list directory
MtGeneList <- read.delim(file = "MtGeneList_GeneID_ENSG.txt") ## files provided in github rb_mt_gene_list directory

print("Calculating Mt %")
seurat <- PercentageFeatureSet(seurat, features = MtGeneList$GeneID, col.name = "percent.mt")

print("Calculating Rb %")
RbGeneList <- RbGeneList[which(RbGeneList$GeneID %in% rownames(seurat)),]
seurat <- PercentageFeatureSet(seurat, features = RbGeneList$GeneID, col.name = "percent.rb")



##### Normalize and Plot results #####
### Make pre-QC figures ###
seurat <- SCTransform(seurat, method = "glmGamPoi", verbose = FALSE)
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
seurat <- FindNeighbors(seurat, dims = 1:10)
seurat <- FindClusters(seurat, resolution = 0.5)
seurat <- RunUMAP(seurat, dims = 1:10)


saveRDS(seurat, paste0(outdir,"seurat_all_cells.rds"))


seurat@meta.data$GSM_Individual <- paste0(seurat@meta.data$GSM, "_", seurat@meta.data$Individual)


### QC Figures pre filtering ###
plot_mt_pct <- VlnPlot(seurat, features = c( "percent.mt"), group.by = "Individual", pt.size = 0)

plot_rb_pct <- VlnPlot(seurat, features = c( "percent.rb"), group.by = "Individual", pt.size = 0)

plot_n_count <- VlnPlot(seurat, features = c( "nCount_RNA"), group.by = "Individual", pt.size = 0)

plot_nFeature_RNA <- VlnPlot(seurat, features = c( "nFeature_RNA"), group.by = "Individual", pt.size = 0)

plot_mt_pct <- VlnPlot(seurat, features = c( "percent.mt"), group.by = "GSM_Individual", pt.size = 0)

plot_rb_pct <- VlnPlot(seurat, features = c( "percent.rb"), group.by = "GSM_Individual", pt.size = 0)

plot_n_count <- VlnPlot(seurat, features = c( "nCount_RNA"), group.by = "GSM_Individual", pt.size = 0)

plot_nFeature_RNA <- VlnPlot(seurat, features = c( "nFeature_RNA"), group.by = "GSM_Individual", pt.size = 0)

lib_mt <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Individual")

lib_genes <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Individual")

nebulosa_mt_umap <- plot_density(seurat, "percent.mt", pal = "plasma")

nebulosa_rb_umap <- plot_density(seurat, "percent.rb", pal = "plasma")

nebulosa_nFeature_RNA_umap <- plot_density(seurat, "nFeature_RNA", pal = "plasma")

nebulosa_nCount_RNA_umap <- plot_density(seurat, "nCount_RNA", pal = "plasma")

mt_umap <- FeaturePlot(seurat, features = "percent.mt")

rb_umap <- FeaturePlot(seurat, features = "percent.rb")

nFeature_RNA_umap <- FeaturePlot(seurat, features = "nFeature_RNA")

nCount_RNA_umap <- FeaturePlot(seurat, features = "nCount_RNA")

Individual_umap <- DimPlot(seurat, group.by = "Individual")

Individual_umap_ <- DimPlot(seurat, group.by = "Individual", .by = "Individual")




##### FILTER #####
### Use 5% mt% as filtering threshold (used for original manuscript: "Generation of hypothalamic arcuate organoids from human induced pluripotent stem cells", Huang et al)
seurat_sub <- subset(seurat, subset = percent.mt < 5 & nFeature_RNA > 400)


##### Normalize and Plot results #####
### Make pre-QC figures ###
seurat_sub <- SCTransform(seurat_sub, method = "glmGamPoi", vars.to.regress = "percent.mt") 
seurat_sub <- RunPCA(seurat_sub, features = VariableFeatures(object = seurat_sub))
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:30)
seurat_sub <- FindClusters(seurat_sub, resolution = 0.2, dims = 1:30) 
seurat_sub <- RunUMAP(seurat_sub, dims = 1:30)

outdir_filtered <- paste0(outdir,"filtered/")
dir.create(outdir_filtered)

##### Take a preliminary look at the clusters and features definining them
umap <- DimPlot(seurat_sub, label = TRUE, repel = TRUE)


saveRDS(seurat_sub, paste0(outdir_filtered,"seurat_filtered.rds"))
seurat_sub <- readRDS(paste0(outdir_filtered,"seurat_filtered.rds"))


### QC Figures pre filtering ###
plot_mt_pct <- VlnPlot(seurat_sub, features = c( "percent.mt"), group.by = "Individual", pt.size = 0)

plot_rb_pct <- VlnPlot(seurat_sub, features = c( "percent.rb"), group.by = "Individual", pt.size = 0)

plot_n_count <- VlnPlot(seurat_sub, features = c( "nCount_RNA"), group.by = "Individual", pt.size = 0)

plot_nFeature_RNA <- VlnPlot(seurat_sub, features = c( "nFeature_RNA"), group.by = "Individual", pt.size = 0)

plot_mt_pct <- VlnPlot(seurat_sub, features = c( "percent.mt"), group.by = "GSM_Individual", pt.size = 0)

plot_rb_pct <- VlnPlot(seurat_sub, features = c( "percent.rb"), group.by = "GSM_Individual", pt.size = 0)

plot_n_count <- VlnPlot(seurat_sub, features = c( "nCount_RNA"), group.by = "GSM_Individual", pt.size = 0)

plot_nFeature_RNA <- VlnPlot(seurat_sub, features = c( "nFeature_RNA"), group.by = "GSM_Individual", pt.size = 0)

lib_mt <- FeatureScatter(seurat_sub, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Individual")

lib_genes <- FeatureScatter(seurat_sub, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Individual")

nebulosa_mt_umap <- plot_density(seurat_sub, "percent.mt", pal = "plasma")

nebulosa_rb_umap <- plot_density(seurat_sub, "percent.rb", pal = "plasma")

nebulosa_nFeature_RNA_umap <- plot_density(seurat_sub, "nFeature_RNA", pal = "plasma")

nebulosa_nCount_RNA_umap <- plot_density(seurat_sub, "nCount_RNA", pal = "plasma")

mt_umap <- FeaturePlot(seurat_sub, features = "percent.mt")

rb_umap <- FeaturePlot(seurat_sub, features = "percent.rb")

nFeature_RNA_umap <- FeaturePlot(seurat_sub, features = "nFeature_RNA")

nCount_RNA_umap <- FeaturePlot(seurat_sub, features = "nCount_RNA")

Individual_umap <- DimPlot(seurat_sub, group.by = "Individual")

Individual_umap_split <- DimPlot(seurat_sub, group.by = "Individual", split.by = "Individual")



##### Assign Classifications to Clusters #####
classifications <- c("0" = "Neurons", "1" = "Astrocytes", "2" = "Oligodendrocytes", "3" = "Neurons", "4" = "Neurons", "5" = "OPC", "6" = "Microglia", "7" = "Neurons", "8" = "Vascular & Leptomeningeninal Cells (VLMCs)", "9" = "Neurons", "10" = "Ependymal Cells", "11" = "Neurons", "12" = "Immature Oligodendrocytes", "13" = "Endothelial Cells", "14" = "Neurons", "15" = "Neurons", "16" = "Unknown", "17" = "Tanycytes", "18" = "Neurons", "19" = "Neurons", "20" = "Neurons")
classification_dt <- data.table(cluster = names(classifications), Classification = classifications)

markers_all[classification_dt, on = "cluster"]

fwrite(markers_all,  paste0(outdir_filtered,"markers_all_clusters.tsv"), sep = "\t")


seurat_sub@meta.data$Classification <- classification_dt[data.table(Cluster = seurat_sub@meta.data[,"SCT_snn_res.0.2"]), on = "Cluster"]$Classification

##### 1. Remove unknown population #####
seurat_sub <- subset(seurat_sub, subset = Cell_Types_Merged_Simple != "Unknown")


saveRDS(seurat_sub, paste0(outdir_filtered, "seurat_filtered_classified.rds"))


