### Reason: Test differentially expressed genes in DNs vs other neuron subtypes
### Date: 11 March, 2022
### Author: Drew Neavin

##### Load in Libraries #####
library(Seurat)
library(data.table)


##### Set up Directories #####
dir <- "/path/to/base/directory/"
neuron_dir <- paste0(dir,"/output/Classification/Neuron_Clustering/Train_Test_Reannotation/")
outdir <- paste0(dir, "output/Classification/Neuron_Clustering/DN_DE/")
dir.create(outdir, recursive = TRUE)



##### Read in Data #####
seurat_neurons <- readRDS(paste0(neuron_dir, "neurons_classified_updated_names_sct_map.rds"))


##### Classify if cells are double positive (tdtomato high and detectable TH or tdtomato RNA) - pulling the most extremes of these groups (DN double positive vs non-DN double negative) #####
seurat_neurons@meta.data$`FACS+/TH+` <- paste0(ifelse(seurat_neurons@meta.data$pool == "WT_THpos", "+", "-"), "/", ifelse(colSums(seurat_neurons[["RNA"]][c("ENSG00000180176", "tdtomato"),]) > 0, "+", "-"))

seurat_neurons@meta.data$`Neuron_Subtypes_Broad_FACS+/TH+` <- paste0(seurat_neurons$Neuron_Subtypes_Broad, "_", seurat_neurons@meta.data$`FACS+/TH+`)

seurat_neurons_sub <- subset(seurat_neurons, subset = `Neuron_Subtypes_Broad_FACS+/TH+` %in% c("DNs_+/+", paste0(unique(seurat_neurons$Neuron_Subtypes_Broad)[!unique(seurat_neurons$Neuron_Subtypes_Broad) %in% "DNs"], "_", "-/-")))

unique(seurat_neurons_sub$`Neuron_Subtypes_Broad_FACS+/TH+`)


##### Test DEGs #####
Idents(seurat_neurons_sub) <- "Neuron_Subtypes_Broad_FACS+/TH+"

dn_degs <- FindMarkers(seurat_neurons_sub, test.use = "MAST", ident.1 = "DNs_+/+", logfc.threshold = 0, latent.vars = c("combined_classification", "percent.rb", "nCount_RNA")) ## correct for ribosomal %, number of UMIs and replicate

dn_degs$ENSG <- rownames(dn_degs)


### Add in Gene IDs to tables ###
gene_conversions <- fread("/path/to/10x/matrix/features.tsv.gz", sep = "\t", col.names = c("ENSG", "Gene_ID", "Assay")) ### So can convert ENSG IDs to gene IDs for easy interpretation of DEGs
gene_conversions$Assay <- NULL


dn_degs <- gene_conversions[dn_degs, on = "ENSG"]

fwrite(dn_degs, paste0(outdir, "DN_TH+_vs_neurons_TH-_DEGs.tsv.gz"), sep = "\t")


