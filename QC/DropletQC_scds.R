##### Intro #####
### Author: Drew Neavin
### Date: 15 August, 2021
### Reason: Assess for doublets and empty droplets based on transcriptome with DropletQC and scds resepectively

##### Load in Libraries #####
library(scds)
library(DropletQC)
library(data.table)
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)


##### Set up directories
dir <- "/path/to/base/directory/"
datadir <- paste0(dir,"data/")
outdir <- paste0(dir,"output/QC/DropletQC_scds/")
dir.create(outdir, recursive = TRUE)


### Make a list of the 10x directories to be analyzed
dir_list <- dir(datadir, pattern = "GEX")
pools <- gsub("_GEX_0_1_HC3GJDSXY", "", dir_list) %>% gsub("\\d+_GEX_", "", .)



##### scds #####
## Read in data
counts_list <- lapply(dir_list, function(x){
	Read10X(paste0(datadir,x, "/outs/filtered_feature_bc_matrix/"), gene.column = 1)
})
names(counts_list) <- pools


## Add poolnames to cell names so can easily match back if there are any droplets with same barcode in different pools
counts_list <- lapply(names(counts_list), function(x){
	colnames(counts_list[[x]]) <- paste0(x, "_", colnames(counts_list[[x]]))
	return(counts_list[[x]])
})
names(counts_list) <- pools


sce <- lapply(counts_list, function(x){
	SingleCellExperiment(list(counts=x))
})


## Annotate doublet using binary classification based doublet scoring:
sce = lapply(sce, function(x){
	bcds(x, retRes = TRUE, estNdbl=TRUE)
})


## Annotate doublet using co-expression based doublet scoring:
sce = lapply(sce, function(x){
	cxds(x, retRes = TRUE, estNdbl=TRUE)
})


### If cxds worked, run hybrid, otherwise use bcds annotations
	## Combine both annotations into a hybrid annotation
sce = lapply(sce, function(x){
	if ("cxds_score" %in% colnames(colData(x))) { ### sometimes cxds fails so need to account for that possiblity
		cxds_bcds_hybrid(x, estNdbl=TRUE)
	}
})

Doublets <- lapply(sce, function(x){
	if ("cxds_score" %in% colnames(colData(x))){
		as.data.frame(cbind(rownames(colData(x)), colData(x)$hybrid_score, colData(x)$hybrid_call))
	} else {
		print("this pool failed cxds so results are just the bcds calls")
		as.data.frame(cbind(rownames(colData(x)), colData(x)$bcds_score, colData(x)$bcds_call))
		}
})


## Doublet scores are now available via colData:
Doublets <- lapply(Doublets, function(x){
	colnames(x) <- c("Barcode","scds_score","scds_DropletType")
	x$scds_DropletType <- gsub("FALSE","singlet",x$scds_DropletType) 
	x$scds_DropletType <- gsub("TRUE","doublet",x$scds_DropletType)
	return(x)
})


message("writing output")
lapply(names(Doublets), function(x){
	fwrite(Doublets[[x]], paste0(outdir,x,"_scds_doublets.tsv"), sep = "\t")
})


Doublets_dt <- data.table(do.call(rbind, Doublets))


##### DropletQC #####
##### Calculate nuclear fraction #####
if (!file.exists(paste0(outdir, "nuclear_fraction_predictions.rds"))){
	nf1 <- lapply(dir_list, function(x){
		nuclear_fraction_tags(outs = paste0(datadir, x, "/outs"))
	})
	names(nf1) <- pools
	saveRDS(nf1, paste0(outdir, "nuclear_fraction_predictions.rds"))
} else {
	print("reading in file already produced")
	nf1 <- readRDS(paste0(outdir, "nuclear_fraction_predictions.rds"))
}

### Add UMIs to sce object (will need to identify empty droplets) ###
sce <- lapply(sce, function(x){
	colData(x)$umi <- colSums(counts(x))
	return(x)
})

### Add UMIs to nf1 object ###
nf1 <- lapply(names(nf1), function(x){
	nf1[[x]]$umi <- colData(sce[[x]])$umi
	return(nf1[[x]])
})
names(nf1) <- pools


### Identify empty droplets
nf1 <- lapply(nf1, function(x){
	identify_empty_drops(x)
})



nf1 <- lapply(names(nf1), function(pool){
	nf1[[pool]]$pool <- pool
	nf1[[pool]]$Barcode <- paste0(pool, "_", rownames(nf1[[pool]]))
	return(nf1[[pool]])
})

nuclear_fractions <- data.table(do.call(rbind, nf1))

scds_nuclear_fraction <- Doublets_dt[nuclear_fractions, on = "Barcode"]


p_nuclear_fractions <- ggplot(scds_nuclear_fraction, aes(x = nuclear_fraction)) +
	geom_density() +
	theme_classic() +
	facet_wrap(vars(pool), nrow = 1) + 
	theme(legend.title = element_blank()) 

ggsave(p_nuclear_fractions, filename = paste0(outdir,"nuclear_fractions.png"), width = 4, height = 3)


p_nf_scds <- ggplot(scds_nuclear_fraction, aes(x = scds_DropletType, y = nuclear_fraction)) +
		geom_violin() +
		theme_classic() +
		facet_wrap(vars(pool), nrow = 1) + 
		stat_summary(fun = "mean", size = 1,
		geom="point", color="firebrick3")

ggsave(p_nf_scds, filename = paste0(outdir,"nuclear_fractions_scds.png"), width = 4, height = 3)


fwrite(scds_nuclear_fraction, paste0(outdir,"scds_dropletQC_results.tsv"), sep = "\t")