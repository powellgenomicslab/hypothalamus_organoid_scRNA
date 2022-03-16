##### Intro #####
### Author: Drew Neavin
### Date: 5 June, 2021
### Reason: Demultiplex hashed data using seurat pipeline functions

##### Load in Libraries #####
library(data.table)
library(Seurat)



##### Set up directories
dir <- "/path/to/base/directory/"
datadir <- paste0(dir,"data/")
hashdir <- paste0(dir,"data/HTO/")
outdir <- paste0(dir,"output/QC/Hash_Demultiplexing/")
dir.create(outdir, recursive = TRUE)


### Make a list of the 10x directories to be analyzed
dir_list <- dir(datadir, pattern = "GEX")
hashdir_list <- dir(hashdir)

pools <- gsub("_GEX_0_1_HC3GJDSXY", "", dir_list) %>% gsub("\\d+_GEX_", "", .)



## Read in data
counts_list <- lapply(dir_list, function(x){
	Read10X(paste0(datadir,x,"/outs/filtered_feature_bc_matrix/"), gene.column = 1)
})
names(counts_list) <- pools


## Add poolnames to cell names so can easily match back if there are any droplets with same barcode in different pools
counts_list <- lapply(names(counts_list), function(x){
	colnames(counts_list[[x]]) <- paste0(x, "_", colnames(counts_list[[x]]))
	return(counts_list[[x]])
})
names(counts_list) <- pools


### Read in the hashing data ###
hash_list <- lapply(hashdir_list, function(x){
	Read10X(paste0(hashdir,x), gene.column =1)
})
names(hash_list) <- rev(pools)


### Add pool names to barcodes so all are unique ###
hash_list <- lapply(names(hash_list), function(x){
	colnames(hash_list[[x]]) <- paste0(x, "_", colnames(hash_list[[x]]))
	return(hash_list[[x]])
})
names(hash_list) <- rev(pools)


### Remove the unmapped row (creates noise when identifying doublets/singlets)
hash_list <- lapply(hash_list, function(x){
	x[-c(which(rownames(x) == "unmapped")),]
})

### Check sizes of both the counts and hashing matrices to see if they are the same
lapply(hash_list, dim)
lapply(counts_list, dim)
## The counts dataframe has more so subset to the barcodes in the hashing data


## Get just the cells that are also in the hash list ##
counts_list_new <- lapply(names(counts_list), function(x){
	print(x)
	colnames(counts_list[[x]]) <- gsub("-1", "", colnames(counts_list[[x]]))
	counts_list[[x]] <- counts_list[[x]][,which(colnames(counts_list[[x]]) %in% colnames(hash_list[[x]]))]
	return(counts_list[[x]])
})
names(counts_list_new) <- names(counts_list)

### Check sizes of both the counts and hashing matrices to see if they are the same
lapply(hash_list, dim)
lapply(counts_list_new, dim)


## Get just the cells that are also in the counts list ##
hash_list_new <- lapply(names(counts_list_new), function(x){
	hash_list[[x]][,which(colnames(hash_list[[x]]) %in% colnames(counts_list_new[[x]]))]
})
names(hash_list_new) <- names(counts_list_new)

### Check sizes of both the counts and hashing matrices to see if they are the same
lapply(hash_list_new, dim)
lapply(counts_list_new, dim)
## They are the same now, indicating that the same barcodes are in each group


### Create seurat object from counts data
seurat_list_norm <- lapply(counts_list_new, function(x){
        CreateSeuratObject(counts = x)
    })
names(seurat_list_norm) <- pools

## Normalize, find variable features and scale counts data
seurat_list_norm <- lapply(seurat_list_norm, function(x){
	tmp <- NormalizeData(x, verbose = TRUE)
	tmp <- FindVariableFeatures(tmp, selection.method = "mean.var.plot")
	tmp <- ScaleData(tmp, features = VariableFeatures(tmp))
	return(tmp)
})
names(seurat_list_norm) <- pools

### Add hashtag data ###
seurat_list_norm <- lapply(names(seurat_list_norm), function(x){
	seurat_list_norm[[x]][["HTO"]] <- CreateAssayObject(counts = hash_list_new[[x]])
	return(seurat_list_norm[[x]])
})
names(seurat_list_norm) <- pools



### Normalize hashtag data
seurat_list_norm <- lapply(seurat_list_norm, function(x){
	NormalizeData(x, assay = "HTO", normalization.method = "RC")
})


seurat_list_norm <- lapply(seurat_list_norm, function(x){
	MULTIseqDemux(x, assay = "HTO", autoThresh = TRUE)
})


### Save the seurat objects to file ###
saveRDS(seurat_list_norm, paste0(outdir,"seurat_list.rds"))

