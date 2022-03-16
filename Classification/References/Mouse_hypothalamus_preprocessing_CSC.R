##### Intro #####
### Author: Drew Neavin
### Date: 12 July, 2021
### Reason: Plot markers for Cell Stem Cell mouse 
### Original dataset from Zhang et al, "Cascade diversification directs generation of neuronal diversity in the hypothalamus"
### Data downloaded from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151060


##### Load in Libraries #####
library(data.table)
library(Seurat)
library(ggplot2)
library(scPred)
library(tidyverse)
library(Nebulosa)
library(biomaRt)



##### Set up Directories #####
dir <- "/path/to/base/directory/"
data_dir <- paste0(dir, "data/mouse_hypothalamus_development_CSC/")
outdir <- paste0(dir,"output/Classification/Mouse_hypothalamus_CSC/time_markers/")
dir.create(outdir, recursive = TRUE)



##### Read in Data #####
counts <- fread(paste0(data_dir,"GSE151060_HY_integrated_43261cells_counts.txt.gz"))
counts <- as.data.frame(counts)
colnames(counts) <- gsub("\\.", "-", colnames(counts))
rownames(counts) <- counts$V1
counts$V1 <- NULL


meta <- fread(paste0(data_dir,"GSE151060_HY_integrated_43261cells_metadata.csv.gz"))
rownames(meta) <- meta$V1
meta$V1 <- NULL

seurat <- CreateSeuratObject(counts = counts, meta.data = meta)

seurat@meta.data$time.sum <- gsub(" ", "", seurat@meta.data$time.sum)
colnames(seurat@meta.data) <- gsub("time.sum", "time_sum", colnames(seurat@meta.data))



### Humanize data ###
##### Update the gene names
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")


mouse_genes <- data.table(ENSM = rownames(seurat))
genesV2 = data.table(getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(seurat), mart = mouse, attributesL = c("ensembl_gene_id"), martL = human, uniqueRows=TRUE))


colnames(mouse_genes) <- c("MGI.symbol")
updated_conversion <- genesV2[match(mouse_genes$MGI.symbol, genesV2$MGI.symbol)]
conversion <- data.table(mgi_symbol = mouse_genes$MGI.symbol, ENSG = updated_conversion$Gene.stable.ID)
conversion$ENSG <- ifelse(is.na(conversion$ENSG), conversion$mgi_symbol, conversion$ENSG)


##### Make humanized seurat object
if (all(rownames(counts) == conversion$mgi_symbol)){
	rownames(counts) <- make.unique(conversion$ENSG)
} else {
	print("Rownames don't match the conversion dataframe")
}

### Only keep the rows that are not NA
seurat_humanized <- CreateSeuratObject(counts, meta.data = meta)

feature_meta <- data.table(ENSG = conversion$ENSG)
rownames(feature_meta) <- rownames(seurat_humanized)


seurat_humanized[["RNA"]] <- AddMetaData(seurat_humanized[["RNA"]], feature_meta)

### Remove cells with >10% mt reads ###
seurat_humanized[["percent.mt"]] <- PercentageFeatureSet(object = seurat_humanized, pattern = "^MT-")


## Only keep cells with < 10% mt reads
seurat_humanized <- subset(seurat_humanized, subset = percent.mt < 10)


## Normalize Data ##
seurat_humanized <- SCTransform(seurat_humanized, method = "glmGamPoi")
seurat_humanized <- RunPCA(seurat_humanized)
seurat_humanized <- RunUMAP(seurat_humanized, dims = 1:30)
seurat_humanized <- FindNeighbors(seurat_humanized, dims = 1:30)
seurat_humanized <- FindClusters(seurat_humanized)

saveRDS(seurat_humanized, paste0(outdir,"mouse_humanized_seurat.rds"))

