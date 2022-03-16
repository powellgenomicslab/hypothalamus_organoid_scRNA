#############################
### Author: Drew  Neavin
### Date: 16 August, 2021
### Reason: integrate the snRNA-seq neonatal human dataset with the mouse cell stem cell and nature references to check for potential reclassifications and possibility of using them together as ref for scpred
#############################


library(data.table)
library(Seurat)
library(tidyverse)


##### Set up Directories #####
dir <- "/path/to/base/directory/"
neonatal_dir <- datadir <- paste0(dir,"output/Classification/Neonatal_snRNAseq/Preprocessing/filtered/")
mouse_csc_dir <- datadir <- paste0(dir,"output/Classification/Mouse_hypothalamus_CSC/time_markers/")
mouse_nature_dir <- paste0(dir,"output/Classification/Mouse_hypothalamus_nature/Preprocessing/")
outdir <- datadir <- paste0(dir,"output/Classification/Integrate_refs/")
dir.create(outdir, recursive = TRUE)


### Read in prepared datasets
seurat_mouse_csc <- readRDS(paste0(mouse_csc_dir,"mouse_humanized_seurat.rds"))
seurat_mouse_nature <- readRDS(paste0(mouse_nature_dir,"mouse_humanized_seurat_unique_genes.rds"))
seurat_neonatal <- readRDS(paste0(neonatal_dir,"seurat_filtered_classified.rds"))

DefaultAssay(seurat_neonatal) <- "RNA"
seurat_neonatal[["SCT"]] <- NULL




## Need to update the rownames for the neonatal object to ENSG IDs to match mouse csc ##
### Make a reference from Gene IDs to ENSG IDs ###
GeneConversion <- fread("/directflow/SCCGGroupShare/projects/DrewNeavin/Parkinsons_TH_Collab/data/GEX/GEX_WT_TH/features.tsv.gz", header = F, sep = "\t")
GeneConversion$V3 <- NULL
colnames(GeneConversion) <- c("ENSG_ID", "Gene_ID")

row_gene <- data.table(ENSG_ID = rownames(seurat_mouse_csc[["RNA"]]@counts))
row_gene_nature <- data.table(ENSG_ID = rownames(seurat_mouse_nature[["RNA"]]@counts))

updated_names <- GeneConversion[row_gene, on = "ENSG_ID"]
updated_names_nature <- GeneConversion[row_gene_nature, on = "ENSG_ID"]

updated_names$Gene_ID <- make.unique(ifelse(is.na(updated_names$Gene_ID), "unknown", updated_names$Gene_ID))
updated_names_nature$Gene_ID <- make.unique(ifelse(is.na(updated_names_nature$Gene_ID), "unknown", updated_names_nature$Gene_ID))

counts <- seurat_mouse_csc[["RNA"]]@counts
rownames(counts) <- gsub("_", "-",updated_names$Gene_ID)
rownames(counts) <- gsub("\\.", "-",updated_names$Gene_ID)

seurat_mouse_csc_new <- CreateSeuratObject(counts = counts, meta.data = seurat_mouse_csc@meta.data)


counts2 <- seurat_mouse_nature[["RNA"]]@counts
rownames(counts2) <- gsub("_", "-",updated_names_nature$Gene_ID)
rownames(counts2) <- gsub("\\.", "-",updated_names_nature$Gene_ID)

seurat_mouse_nature_new <- CreateSeuratObject(counts = counts2, meta.data = seurat_mouse_nature@meta.data)


ref_list <- list(seurat_mouse_csc_new, seurat_mouse_nature_new, seurat_neonatal)



# normalize and identify variable features for each dataset independently
ref_list <- lapply(ref_list, FUN = function(x) {
	DefaultAssay(x) <- "RNA"
    x <- SCTransform(x)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = ref_list, nfeatures = 3000)


ref_list <- PrepSCTIntegration(object.list = ref_list, anchor.features = features)


anchors <- FindIntegrationAnchors(object.list = ref_list, anchor.features = features, normalization.method = "SCT")


# this command creates an 'integrated' data assay
ref_combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")


# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(ref_combined) <- "integrated"


### Add classification columns ###
ref_combined@meta.data$Origin <- ifelse(!is.na(ref_combined@meta.data$Classification), "Human Neonatal", 
									ifelse(!is.na(ref_combined@meta.data$CellType.sum), "Mouse CSC", "Mouse Nature"))

ref_combined@meta.data$age <- ifelse(ref_combined@meta.data$Individual == "5900", "H0",
								ifelse(ref_combined@meta.data$Individual == "4438", "H42",
									ifelse(ref_combined@meta.data$Individual == "4396", "H54",
										ifelse(ref_combined@meta.data$Individual == "4389", "H79",
											ifelse(ref_combined@meta.data$Individual == "4368", "H141",
												ifelse(ref_combined@meta.data$Individual == "4458", "H225", NA))))))

ref_combined@meta.data$age <- ifelse(!is.na(ref_combined@meta.data$time), gsub("\\.R\\d", "", ref_combined@meta.data$time), ref_combined@meta.data$age)

ref_combined@meta.data$age <- ifelse(!is.na(ref_combined@meta.data$Age), ref_combined@meta.data$Age, ref_combined@meta.data$age)

ref_combined@meta.data$age <- factor(ref_combined@meta.data$age, levels = c("E11", "E14", "E15", "E17", "P0", "P2", "P7", "P10", "P23", "H0", "H42", "H54", "H79", "H141", "H225"))


ref_combined@meta.data$Original_Cell_Types <- ifelse(!is.na(ref_combined@meta.data$Classification),ref_combined@meta.data$Classification,
												ifelse(!is.na(ref_combined@meta.data$CellType.sum), ref_combined@meta.data$CellType.sum, ref_combined@meta.data$Annotation))




# Run the standard workflow for visualization and clustering

ref_combined <- RunPCA(ref_combined, npcs = 100)
ref_combined <- RunUMAP(ref_combined, reduction = "pca", dims = 1:100, min.dist = 0.2, spread = 1.5)
ref_combined <- FindNeighbors(ref_combined, reduction = "pca", dims = 1:100)



saveRDS(ref_combined, paste0(outdir, "integrated_human_mouse_csc_nature.rds"))


