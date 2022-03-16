##### Intro #####
### Author: Drew Neavin
### Date: 5 August, 2021
### Reason: Plot markers for Cell Stem Cell mouse 
### Original dataset from Romanov et al, "Molecular design of hypothalamus development"
### Data downloaded from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132730


library(Seurat)
library(tidyverse)
library(data.table)
library(biomaRt)
library(Nebulosa)


dir <- "/path/to/base/directory/"
datadir <- paste0(dir,"data/mouse_hypothalamus_development_nature/")
outdir <- paste0(dir,"output/Classification/Mouse_hypothalamus_nature/Preprocessing/")
dir.create(outdir, recursive = TRUE)


##### Read in Files #####
seurat <- readRDS(paste0(datadir,"GSE132730_TractNAE_integrated.rds"))


##### Update the gene names
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")


mouse_genes <- data.table(ENSM = rownames(seurat[["RNA"]]))
genesV2 = data.table(getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(seurat[["RNA"]]), mart = mouse, attributesL = c("ensembl_gene_id"), martL = human, uniqueRows=TRUE))


updated_conversion <- genesV2[match(mouse_genes$ENSM, genesV2$MGI.symbol)]
conversion <- data.table(ENSM = mouse_genes$ENSM, ENSG = updated_conversion$Gene.stable.ID)


### For integrated ###
mouse_genes_integrated <- data.table(ENSM = rownames(seurat))
genesV2_integrated = data.table(getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(seurat), mart = mouse, attributesL = c("ensembl_gene_id"), martL = human, uniqueRows=TRUE))

updated_conversion_integrated <- genesV2_integrated[match(mouse_genes_integrated$ENSM, genesV2_integrated$MGI.symbol)]
conversion_integrated <- data.table(ENSM = mouse_genes_integrated$ENSM, ENSG = updated_conversion_integrated$Gene.stable.ID)




##### Make humanized seurat object
counts_merged <- GetAssayData(seurat, assay = "RNA", slot = "counts")

if (all(rownames(counts_merged) == conversion$ENSM)){
	rownames(counts_merged) <- conversion$ENSG
} else {
	print("Rownames don't match the conversion dataframe")
}


### Only keep the rows that are not NA
counts_merged <- counts_merged[which(!is.na(rownames(counts_merged))),]

seurat_humanized <- CreateSeuratObject(counts_merged, meta.data = seurat@meta.data)

### Remove cells with >10% mt reads ###
seurat_humanized[["percent.mt"]] <- PercentageFeatureSet(object = seurat_humanized, pattern = "^MT-")

## Only keep cells with < 10% mt reads
seurat_humanized <- subset(seurat_humanized, subset = percent.mt < 10)

## Remove cells with unknown or non-classical cell type annotations ##
seurat_humanized <- subset(seurat_humanized, subset = Cell_Types_Merged_Simple != "Mixed mytotic cells")
seurat_humanized <- subset(seurat_humanized, subset = Cell_Types_Merged_Simple != "Mixed cells")
seurat_humanized <- subset(seurat_humanized, subset = Cell_Types_Merged_Simple != "Astroependymal cells")

saveRDS(seurat_humanized, paste0(outdir,"mouse_humanized_seurat.rds"))


