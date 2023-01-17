import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import os
import anndata as ad
import matplotlib as plt
import numba
import umap
import itertools
import math
from scipy.sparse import csr_matrix
import skmisc
from plotnine import *

datadir = '/directflow/SCCGGroupShare/projects/DrewNeavin/Parkinsons_TH_Collab/output/QC/seurat_QC/'
basedir = '/directflow/SCCGGroupShare/projects/DrewNeavin/Parkinsons_TH_Collab/output/RNAveocity/'
outdir = '/directflow/SCCGGroupShare/projects/DrewNeavin/Parkinsons_TH_Collab/output/RNAveocity/neurons'

os.makedirs(outdir)
os.chdir(outdir)


# scvelo settings
scv.logging.print_version()
# Running scvelo 0.2.3 (python 3.8.5) on 2021-05-31 10:45.
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.set_figure_params('scvelo', figsize=(40,40))  # for beautified visualization
​
# scanpy settings
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
# scanpy==1.7.2 anndata==0.7.6 numpy==1.20.3 scipy==1.6.3 pandas==1.2.4 scikit-learn==0.24.2
sc.settings.set_figure_params(dpi=200)
​

### Load in sample data info
samples = pd.read_table('/directflow/SCCGGroupShare/projects/DrewNeavin/Parkinsons_TH_Collab/scripts/RNAvelocity/velocyto_files.tsv', delim_whitespace=True, header=0)


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
meta = pd.read_csv("/directflow/SCCGGroupShare/projects/DrewNeavin/Parkinsons_TH_Collab/output/Classification/Neonatal_snRNAseq/scPred/metadata_nnet_broad_classification.csv")


### Filter anndata object for droplet ids
adata_concat = adata_concat[np.isin(adata_concat.obs.index,meta["Unnamed: 0"])]


### Order the umap droplets in same orderas adata
adata_barcodes = pd.DataFrame(adata_concat.obs.index)
adata_barcodes = adata_barcodes.rename(columns = {0:'Cell ID'})


meta = meta.rename(columns = {'Unnamed: 0':'Cell ID'})
meta_ordered = adata_barcodes.merge(meta, on = "Cell ID")

### Update metadata columns
adata_concat.obs['scpred_prediction'] = meta_ordered['scpred_prediction'].values


adata_sep = adata_concat[meta_ordered['scpred_prediction'].isin(["Neurons"])]


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
scv.tl.velocity(adata_sep, mode = "stochastic")
scv.tl.velocity_graph(adata_sep)
scv.tl.velocity_pseudotime(adata_sep)



### Try with scvelo umap
os.chdir(outdir)

## UMAP + clustering
# Identify highly-variable genes.
sc.pp.highly_variable_genes(adata_sep, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata_sep, save= "variable_genes.png")

# Scale each gene to unit variance. Clip values exceeding standard deviation 10.
sc.pp.scale(adata_sep, max_value=10)

# Reduce the dimensionality of the data by running principal component analysis (PCA), which reveals the main axes of variation and denoises the data.
sc.tl.pca(adata_sep, svd_solver='arpack')

# compute the neighborhood graph of cells using the PCA representation of the data matrix. 
# You might simply use default values here. For the sake of reproducing Seurat’s results, let’s take the following values.
sc.pp.neighbors(adata_sep, n_neighbors=10, n_pcs=40)

# Embedding the neighborhood graph - reccomend embedding the graph in 2 dimensions using UMAP (McInnes et al., 2018).
#  It is potentially more faithful to the global connectivity of the manifold than tSNE, i.e., it better preservers trajectories.
#  In some ocassions, you might still observe disconnected clusters and similar connectivity violations. They can usually be remedied by running:
#tl.paga(adata)
#pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
#tl.umap(adata, init_pos='paga')
sc.tl.umap(adata_sep)

# Clustering the neighborhood graph - recommend the Leiden graph-clustering method (community detection based on optimizing modularity) by Traag *et al.* (2018).
# Note that Leiden clustering directly clusters the neighborhood graph of cells, which we already computed in the previous section.
sc.tl.leiden(adata_sep)
sc.pl.umap(adata_sep, color=['leiden'], save= "umap.pdf")
scv.pl.velocity_embedding_stream(adata_sep, basis='umap', color=['leiden'], save= "stream.png", size=10)
scv.pl.velocity_embedding_stream(adata_sep, basis='umap', color=['scpred_prediction'], save= "stream_broad_cell_type.png", size=10)

# Creates a plot AND adds embedded velocity vectors 'velocity_umap' to the obsm slot
scv.pl.velocity_embedding(adata_sep, basis='umap', arrow_length=2, arrow_size=1.5, dpi=150, color=['leiden'], save= "arrow.png")
scv.pl.velocity_embedding(adata_sep, basis='umap', arrow_length=2, arrow_size=1.5, dpi=150, color=['scpred_prediction'], save= "arrow_broad_cell_type.png")
# Just a different representation, nothing new is added to adata



### Dynamical model to get latent time
scv.tl.recover_dynamics(adata_sep, max_iter=20)
scv.tl.velocity(adata_sep, mode='dynamical')
scv.tl.velocity_graph(adata_sep)
scv.tl.recover_latent_time(adata_sep)

scv.pl.scatter(adata_sep, color='latent_time', fontsize=24, size=100, color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0,1], save= "latent_time.png")

top_genes = adata_sep.var['fit_likelihood'].sort_values(ascending=False).index[:100]

scv.pl.heatmap(adata_sep, var_names = top_genes[1:100], sortby='latent_time', col_color='scpred_prediction', yticklabels=True, n_convolve=100, figsize=(40, 20), save= "heatmap.png")
scv.pl.heatmap(adata_sep, var_names = top_genes[1:100], sortby='velocity_pseudotime', col_color='scpred_prediction', yticklabels=True, n_convolve=100, figsize=(40, 20), save= "pseudotime_heatmap.png")


scv.pl.scatter(adata_sep, color='scpred_prediction', fontsize=24, size=100, color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0,1], save= "scpred_predicted_clusters.png")
scv.pl.scatter(adata_sep, color='velocity_pseudotime', fontsize=24, size=100, color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0,1], save= "velocity_pseudotime_umap.png")
scv.pl.velocity_embedding_stream(adata_sep, basis='umap', color=['latent_time'], save= "stream_latent_colors.png", size=150)
scv.pl.velocity_embedding_stream(adata_sep, basis='umap', color=['scpred_prediction'], save= "stream_seurat_scpred_prediction.png", size=50)



scv.pl.scatter(adata_sep, basis=top_genes[:15], ncols=5, frameon=False, save = "top_dynamic_genes")


var_names = ['TH','DDC']
scv.pl.scatter(adata_sep, var_names, color=['scpred_prediction'], frameon=False, save = "TH_dynamic_model")
scv.pl.scatter(adata_sep, x='scpred_prediction', y=['TH'], color=['scpred_prediction'], frameon=False, save = "TH_latent")





plot = (ggplot(adata_sep.obs, aes(x='latent_time', color='Time', fill='Time')) + geom_density(alpha=0.1))
plot.save(filename = "density_plot.png")


adata_sep.write(filename=outdir +  "_adata_sep.h5ad")
adata_sep.obs.to_csv(outdir + "_metadata.csv")




adata_sep.write(filename=outdir + "adata_concat.h5ad")
adata_sep = ad.read_h5ad(filename=outdir + "adata_concat.h5ad")





scv.pl.scatter(adata_sep, color='latent_time', fontsize=24, size=100, color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0,1], dpi = 600, save="latent_time.png")
scv.pl.scatter(adata_sep, color='latent_time', fontsize=24, size=100, color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0,1], dpi = 600, save="latent_time_norescale.png")

top_genes = adata_sep.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata_sep, var_names=top_genes, sortby='latent_time', col_color='clusters', n_convolve=100, save="heatmap.png")




