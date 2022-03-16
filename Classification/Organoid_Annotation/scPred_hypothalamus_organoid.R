#############################
### Author: Drew  Neavin
### Reason: Use scPred to predict cell types in the hypothalamus organoids
#############################

library(data.table)
library(Seurat)
library(scPred)
library(Rfast)
library(tidyverse)
library(ggridges)



##### Set up Directories #####
dir <- "/path/to/base/directory/"
QCdir <- paste0(dir,"output/QC/seurat_QC/")
datadir <- paste0(dir,"output/Classification/annotation/")
outdir <- paste0(dir,"output/Classification/scPred/")

barcodes <- "/path/to/10x/matrix/features.tsv.gz" ### use any of the 10x features.tsv.gz for this just to convert ENSG to Gene IDs for different objects

dir.create(outdir, recursive = TRUE)




##### Read in the data #####
seurat <- readRDS(paste0(QCdir,"seurat_singlets.rds")) ## hypothalamus organoid data
ref_seurat <- readRDS(paste0(datadir,"integrated_refs_sct_reannotated_singlets.rds")) ## hypothalamus reference datasets
DefaultAssay(ref_seurat) <- "RNA"



### Process the ref data ###
ref_seurat <- ref_seurat %>% 
	NormalizeData() %>%
	FindVariableFeatures() %>%
	ScaleData() %>%
	RunPCA() %>%
	RunUMAP(dims = 1:30)


### Make label of cell types ###
ref_seurat <- getFeatureSpace(ref_seurat, "CellType")

### Train the models ###
ref_seurat <- trainModel(ref_seurat, model = "mda")


########## Best predictions are with combined model - svmRadial for Radial glia and mda for the rest; tested with different models (not shown) ##########
ref_seurat_combined <- trainModel(ref_seurat, model = "svmRadial", reclassify = c("Radial Glial Cells"))

saveRDS(ref_seurat_combined, paste0(outdir, "classification_ref_seurat_combined.rds"))


### Compare the old and new models
get_scpred(ref_seurat)
get_scpred(ref_seurat_nnet)
get_scpred(ref_seurat_glm)
get_scpred(ref_seurat_combined)

pTop_combined <- plot_probabilities(ref_seurat_combined)
ggsave(pTop_combined, filename = paste0(outdir, "top_classification_probabilities_combined.png"))


### Normalize the hypothalamus data ###
## Need to update the rownames for the seurat object to match those from Allen ##
### Make a reference from ENSG IDs to Gene IDs ###
GeneConversion <- fread("/path/to/10/matrix/features.tsv.gz", header = F, sep = "\t")
GeneConversion$V3 <- NULL
colnames(GeneConversion) <- c("ENSG_ID", "Gene_ID")

row_ensg <- data.table(ENSG_ID = rownames(seurat[["RNA"]]@counts))

updated_names <- GeneConversion[row_ensg, on = "ENSG_ID"]

rownames(seurat[["RNA"]]@counts) <- updated_names$Gene_ID
rownames(seurat[["RNA"]]@data) <- updated_names$Gene_ID



seurat_combined <- scPredict(seurat, ref_seurat_combined, recompute_alignment = FALSE)

### Save predictions ###
saveRDS(seurat_combined, paste0(outdir,"seurat_neonatal_human_predicted_combined.rds"))



### Plot new umap based on predictions
umap_original_top <- DimPlot(seurat_combined, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
ggsave(umap_original_top, filename = paste0(outdir,"original_umap_predictions_top_combined.png"), width = 10)


umap_original_max_prob <- FeaturePlot(seurat_combined, features = "scpred_max", label = TRUE, repel = TRUE)
ggsave(umap_original_max_prob, filename = paste0(outdir,"original_umap_predictions_top_prob_combined.png"))


table(seurat$scpred_prediction)
table(seurat_combined$scpred_prediction)




########## Best Model was MDA + svmRadial ##########
### Check on probability distributions for cells to remove any cells that have relatively high probabilities for another cell type ###
scpred_combined_results <- data.table(seurat_combined@meta.data[,grepl("scpred", colnames(seurat_combined@meta.data))])
rownames(scpred_combined_results) <- colnames(seurat_combined)
scpred_combined_results_raw <- scpred_combined_results
scpred_combined_results_raw[,c("scpred_max","scpred_prediction","scpred_no_rejection"):=NULL]

probability_combined_df <- data.table(top_prob = apply(scpred_combined_results_raw, 1, max),
								top_assignment = gsub("scpred_","",colnames(scpred_combined_results_raw)[apply(scpred_combined_results_raw, 1, which.max)]),
								second_prob = apply(scpred_combined_results_raw, 1, Rfast::nth,k = 2, descending = T),
								second_assignment = gsub("scpred_","",colnames(scpred_combined_results_raw)[apply(scpred_combined_results_raw, 1, Rfast::nth,k = 2, descending = T, index.return = T)]))



probability_combined_df$diff1_2 <- probability_combined_df$top_prob - probability_combined_df$second_prob



### Assign those with less than 0.15 prob difference between top two probabilities as unassigned ###
seurat_combined$CellType <- ifelse(probability_combined_df$diff1_2 < 0.15, "unassigned", seurat_combined$scpred_prediction)

saveRDS(seurat_combined, paste0(outdir,"seurat_broad_classifications.rds"))



