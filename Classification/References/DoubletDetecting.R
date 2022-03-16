######################################
### Reason: Some concerns about possibility of doublets in references; Will use scds and scDblFinder to test for doublets
### Author: Drew Neavin
### Date: 11 October, 2021
### run from demuxafy using the following command: singularity exec --bind /directflow Demuxafy.sif R
######################################


library(scds)
library(scDblFinder)
library(data.table)
library(tidyverse)
library(SingleCellExperiment)
library(Seurat)


### Set up Directories ###
dir <- "/path/to/base/directory/"
data_dir <- paste0(dir,"output/Classification/Integrate_multiple_refs/neonatal_mouse_csc_nature/")
outdir <- paste0(dir, "output/Classification/Integrate_multiple_refs/neonatal_mouse_csc_nature/reannotation/")
dir.create(outdir, recursive = TRUE)


### Read in data ###
##### Read in Data #####
seurat_integrated <- readRDS(paste0(data_dir, "integrated_human_mouse_csc_nature.rds"))

DefaultAssay(seurat_integrated_subset) <- "RNA"
seurat_integrated_subset[["percent.mt"]] <- PercentageFeatureSet(object = seurat_integrated_subset, pattern = "^MT-")


## Only keep cells with < 10% mt reads
seurat_integrated_subset <- subset(seurat_integrated_subset, subset = percent.mt < 10)

seurat_integrated_subset_list <- list()
for (ref in unique(seurat_integrated_subset@meta.data$Origin)){
	tmp <- subset(seurat_integrated_subset, subset = Origin == ref)
	if (ref == "Mouse CSC"){
		for (samp in unique(tmp@meta.data$time))
			if (nrow(tmp@meta.data[which(tmp@meta.data$time == samp),]) > 0){
				seurat_integrated_subset_list[[ref]][[samp]] <- subset(tmp, subset = time == samp)
			}
	}
	if (ref == "Mouse Nature"){
		for (samp in unique(paste0(tmp@meta.data$Age, "_", tmp@meta.data$SampleID))){
			age <- gsub("_.+", "", samp)
			print(age)
			id <- gsub(".+_", "", samp)
			print(id)
			print(samp)
			if (nrow(tmp@meta.data[which(tmp@meta.data$Age == age & tmp@meta.data$SampleID == id),]) > 0){
				seurat_integrated_subset_list[[ref]][[samp]] <- subset(subset(tmp, subset = (Age == age)), subset = SampleID == id)
			}
		}
	}
	if (ref == "Human Neonatal"){
		for (samp in unique(tmp@meta.data$GSM_Individual)){
			if (nrow(tmp@meta.data[which(tmp@meta.data$GSM_Individual == samp),]) > 0){
				seurat_integrated_subset_list[[ref]][[samp]] <- subset(tmp, subset = GSM_Individual == samp)
			}
		}
	}
}

### Convert to SCE Object ###
sce_list <- list()

for (ref in names(seurat_integrated_subset_list)){
	for (time in names(seurat_integrated_subset_list[[ref]])){
		sce_list[[ref]][[time]] <- SingleCellExperiment(list(counts=seurat_integrated_subset_list[[ref]][[time]][["RNA"]]@counts))
	}
}

### scds ###
## Annotate doublet using binary classification based doublet scoring:
sce_list <- lapply(sce_list, function(x){
	lapply(x, function(y){
		bcds(y, retRes = TRUE, estNdbl=TRUE)
	})
})

## Annotate doublet using co-expression based doublet scoring:
sce_list <- lapply(sce_list, function(x){
	lapply(x, function(y){
		cxds(y, retRes = TRUE, estNdbl=TRUE)
	})
})


##### scDblFinder #####
### Calculate Singlets and Doublets ###
sce_list <- lapply(sce_list, function(x){
	lapply(x, function(y){
		scDblFinder(y)
	})
})
	sce_list[[ref]] <- lapply(sce_list[[ref]], function(y){

saveRDS(sce_list, paste0(outdir,"scds_scDblFinder_predicted_sce.rds"))
sce_list <- readRDS(paste0(outdir,"scds_scDblFinder_predicted_sce.rds"))


### If cxds worked, run hybrid, otherwise use bcds annotations
Doublets <- list()

for (x in names(sce_list)){
	for (y in names(sce_list[[x]])){
		if ("cxds_score" %in% colnames(colData(sce_list[[x]][[y]]))) {
			## Combine both annotations into a hybrid annotation
			sce_list[[x]][[y]] = cxds_bcds_hybrid(sce_list[[x]][[y]], estNdbl=TRUE)
			Doublets[[x]][[y]] <- as.data.frame(cbind(rownames(colData(sce_list[[x]][[y]])), colData(sce_list[[x]][[y]])$hybrid_score, colData(sce_list[[x]][[y]])$hybrid_call))
		} else {
			print("this pool failed cxds so results are just the bcds calls")
			Doublets[[x]][[y]] <- as.data.frame(cbind(rownames(colData(sce_list[[x]][[y]])), colData(sce_list[[x]][[y]])$bcds_score, colData(sce_list[[x]][[y]])$bcds_call))
		}
		Doublets[[x]][[y]]$scDblFinder_DropletType <-  sce_list[[x]][[y]]$scDblFinder.class
		Doublets[[x]][[y]]$scDblFinder_Score <- sce_list[[x]][[y]]$scDblFinder.score
	}
}

lapply(Doublets, function(x) lapply(x, function(y) table(y$V3, y$scDblFinder_DropletType)))

Doublets_df_list <- lapply(Doublets, function(x) do.call(rbind, x))

Doublets_df <- do.call(rbind, Doublets_df_list)
colnames(Doublets_df)[1:3] <- c("Barcode", "scds_Score", "scds_DropletType")

Doublets_df$scds_DropletType <- gsub("FALSE", "singlet", Doublets_df$scds_DropletType) %>% gsub("TRUE", "doublet", .)

Doublets_df$Combined_DropletType <- ifelse((Doublets_df$scds_DropletType == "doublet" & Doublets_df$scDblFinder_DropletType == "doublet"), "doublet", "singlet")


fwrite(Doublets_df, paste0(outdir, "singlets_doublet.tsv"), sep = "\t")





