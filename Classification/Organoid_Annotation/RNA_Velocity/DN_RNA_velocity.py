import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import os
import anndata as ad
import matplotlib.pyplot as plt
import numba
import umap
import itertools
import math
from scipy.sparse import csr_matrix
import skmisc
from plotnine import *
import re

### Set up directories
outdir = '/path/to/output/dir/'
meta_file = '/path/to/DN_metadata.csv' ## Available on the hypothalamus_organoid_scRNA github (https://github.com/powellgenomicslab/hypothalamus_organoid_scRNA)

os.makedirs(outdir)
os.chdir(outdir)

### Load in sample data info
samples = pd.read_table('/path/to/velocyto/outdir/velocyto_files.tsv', delim_whitespace=True, header=0)


### Read in data
adata = list()

for i in range(0,2):
	adata.append(scv.read(samples['Directory'][i], cache=True))


### Update droplet names to match seurat names
for i in range(0,2):
	adata[i].obs_names = [row.replace("x", "") for row in [row.replace(":", "_") for row in [row.replace(samples.Pool[i], samples.Name[i]) for row in adata[i].obs_names]]]
	adata[i].var_names_make_unique()


#### Combine the data together
adata_concat = ad.concat(adata)


### Read in metadata
meta = pd.read_csv(meta_file, index_col=False)


### Filter anndata object for droplet ids
adata_concat = adata_concat[np.isin(adata_concat.obs.index,meta["Barcode"])]


### Order the umap droplets in same orderas adata
adata_barcodes = pd.DataFrame(adata_concat.obs.index)
adata_barcodes = adata_barcodes.rename(columns = {0:'Barcode'})
meta_ordered = adata_barcodes.merge(meta, on = "Barcode")

### Update metadata columns
adata_concat.obs['DN_CellType'] = meta_ordered['DN_CellType'].values



### Filter
# Also filter cells with too few unspliced counts
adata_sep.obs['n_unspliced_counts'] = adata_sep.layers['unspliced'].sum(axis=1).A1
adata_sep.obs['n_spliced_counts'] = adata_sep.layers['spliced'].sum(axis=1).A1

## Check the spliced and unspliced counts numbers
plot = (ggplot(adata_sep.obs, aes(x='n_unspliced_counts')) + geom_density(alpha=0.1) + geom_vline(xintercept = 1500))
plot.save(filename = "unspliced_counts_density_plot.png")

plot_spliced_unspliced = (ggplot(adata_sep.obs, aes(x='n_spliced_counts', y = 'n_unspliced_counts')) + geom_point(alpha=0.5))
plot_spliced_unspliced.save(filename = "spliced_unspliced_counts_scatter_plot.png")

adata_sep = adata_sep[adata_sep.obs['n_unspliced_counts'] >= 1500]


# Filter genes
sc.pp.filter_genes(adata_sep, min_cells=20)

# Filter genes that don't have unspliced counts detected in at least 10 cells
adata_sep.var['n_unspliced_counts'] = np.count_nonzero(adata_sep.layers['unspliced'].toarray(), axis=0)
adata_sep.var['log_n_unspliced_counts'] = np.log(adata_sep.var['n_unspliced_counts'])
adata_sep.var['n_spliced_counts'] = np.count_nonzero(adata_sep.layers['spliced'].toarray(), axis=0)

plot = (ggplot(adata_sep.var, aes(x='log_n_unspliced_counts')) + geom_density(alpha=0.1) + geom_vline(xintercept = np.log(10)))
plot.save(filename = "gene_unspliced_counts_density_plot.png")

plot_spliced_unspliced = (ggplot(adata_sep.var, aes(x='n_spliced_counts', y = 'n_unspliced_counts')) + geom_point(alpha=0.5))
plot_spliced_unspliced.save(filename = "gene_spliced_unspliced_counts_scatter_plot.png")


adata_sep = adata_sep[:,adata_sep.var['n_unspliced_counts']>10]




##### RNA Velocity #####
scv.pp.filter_and_normalize(adata_sep)
scv.pp.moments(adata_sep)


### Dynamical model to get latent time
scv.tl.recover_dynamics(adata_sep, max_iter=20)
scv.tl.velocity(adata_sep, mode='dynamical')
scv.tl.velocity_graph(adata_sep)
scv.tl.recover_latent_time(adata_sep)


## UMAP + clustering
# Identify highly-variable genes.
sc.pp.highly_variable_genes(adata_sep, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata_sep, save= "variable_genes.png")

# Scale each gene to unit variance. Clip values exceeding standard deviation 10.
sc.pp.scale(adata_sep, max_value=10)

# Reduce the dimensionality of the data by running principal component analysis (PCA), which reveals the main axes of variation and denoises the data.
sc.tl.pca(adata_sep, svd_solver='arpack')
sc.pp.neighbors(adata_sep, n_neighbors=10, n_pcs=40)

# Embedding the neighborhood graph 
sc.tl.umap(adata_sep)


### Plotting
## Plot Latent Time
scv.pl.scatter(adata_sep, color='latent_time', fontsize=24, size=100, color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0,1], save= "latent_time.png")

## Set up color palette for DN cell types
celltype_palette = {'OTP+ RAX+ DNs 1': '#00b050', 'OTP+ RAX+ DNs 2': '#00b1f0', 'OTP+ RAX+ DNs 3': '#0145cd', 'OTP+ RAX+ DNs 4': '#6f70fe', 'OTP+ RAX- DNs': '#fe85e9', 'OTP- RAX- DNs': '#fed867', 'KISS1-like Neurons': '#ff0000'}

## Plot top latent genes in heatmaps
top_genes = adata_sep.var['fit_likelihood'].sort_values(ascending=False).index[:100]
scv.pl.heatmap(adata_sep, var_names = top_genes[1:100], sortby='velocity_pseudotime', col_color='DN_CellType', palette = celltype_palette, yticklabels=True, n_convolve=100, figsize=(40, 20), save= "pseudotime_heatmap.png")
scv.pl.heatmap(adata_sep, var_names = top_genes[1:100], sortby='velocity_pseudotime', col_color='DN_CellType', palette = celltype_palette, yticklabels=True, n_convolve=100, figsize=(40, 20), save= "pseudotime_heatmap.pdf")

## Plot selected markers
scv.pl.heatmap(adata_sep, var_names = ['RAX','SSH','ASCL1', 'PRDM12', 'OTP','DDC', 'NR5A2', 'TENT5A', 'KCNJ6', 'NR2F2', 'KCNJ3'], sortby='velocity_pseudotime', col_color='DN_CellType', palette = celltype_palette, yticklabels=True, n_convolve=100, save= "pseudotime_heatmap_marker_genes.png", figsize=(4, 2))
scv.pl.heatmap(adata_sep, var_names = ['RAX','SSH','ASCL1', 'PRDM12', 'OTP','DDC', 'NR5A2', 'TENT5A', 'KCNJ6', 'NR2F2', 'KCNJ3'], sortby='velocity_pseudotime', col_color='DN_CellType', palette = celltype_palette, yticklabels=True, n_convolve=100, save= "pseudotime_heatmap_marker_genes.pdf", figsize=(4, 2))

## Plot Cell Types by cluster
scv.pl.scatter(adata_sep, color='DN_CellType', palette = celltype_palette, fontsize=24, size=100, color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0,1], title='DN Subtypes', legend_loc='right margin', save= "DN_CellType_clusters_nolabels.png")
scv.pl.scatter(adata_sep, color='DN_CellType', palette = celltype_palette, fontsize=24, size=100, color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0,1], title='DN Subtypes', legend_loc='right margin', save= "DN_CellType_clusters_nolabels.pdf")

## Plot velocity on umap with various colors of cells
scv.pl.scatter(adata_sep, color='velocity_pseudotime', fontsize=24, size=100, color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0,1], save= "velocity_pseudotime_umap.png")

scv.pl.velocity_embedding_stream(adata_sep, basis='umap', color=['latent_time'], save= "stream_latent_colors.png", size=150)
scv.pl.velocity_embedding_stream(adata_sep, basis='umap', color=['latent_time'], save= "stream_latent_colors.pdf", size=150)

scv.pl.velocity_embedding_stream(adata_sep, basis='umap', color=['velocity_pseudotime'], save= "stream_latent_colors.png", size=150)

scv.pl.velocity_embedding_stream(adata_sep, basis='umap', color=['DN_CellType'], palette = celltype_palette, save= "stream_seurat_DN_CellType.png", size=150, alpha = 0.6, legend_loc='right margin')
scv.pl.velocity_embedding_stream(adata_sep, basis='umap', color=['DN_CellType'], palette = celltype_palette, save= "stream_seurat_DN_CellType.pdf", size=150, alpha = 0.6, legend_loc='right margin')

scv.pl.proportions(adata_sep, groupby='DN_CellType')
scv.pl.velocity_embedding_grid(adata_sep, basis='umap', color='DN_CellType', palette = celltype_palette, save='embedding_grid.png', title='', scale=0.25)


## Save umap axes if want to do some plots with seurat
umap_axes = pd.DataFrame(adata_sep.obsm['X_umap'])
umap_axes = umap_axes.rename({0: 'UMAP1', 1: 'UMAP2'}, axis='columns')
umap_axes['Barcode'] = adata_sep.obs.index

umap_axes.to_csv(outdir + "umap_axes.tsv", sep = "\t")


