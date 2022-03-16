library(scran)
library(data.table)
library(Seurat)


##### Set up directories #####
dir <- "/path/to/base/directory/"
datadir <-paste0(dir,"output/QC/seurat_QC/")
outdir <- paste0(dir,"output/QC/cyclone_cell_cycle/")
dir.create(outdir, recursive = TRUE)


##### Read in the data #####
seurat <- readRDS(paste0(datadir,"seurat_all_cells.rds")) ### From seurat_QC.R script



##### Compute Cell Cycle Using scran package and compare the results from seurat and scran #####
print("Loading in data for cell cycle determination")
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))


assigned <- cyclone(seurat[["RNA"]]@counts, pairs=hs.pairs)  #Note, this takes hours
table(assigned$phases)

assigned_dt <- data.table(phases = assigned$phases, scores = assigned$scores, normalized_scores = assigned$normalized.scores)
rownames(assigned_dt) <- colnames(seurat)

fwrite(assigned_dt, paste0(outdir,"CellCycleProportions.txt"), sep = "\t", row.names = TRUE) #Save so that can read in and don't have to wait to recompute again
