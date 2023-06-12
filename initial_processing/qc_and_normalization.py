import pandas as pd
import os
import scanpy as sc
import scanpy.external as sce
from anndata import AnnData
import seaborn as sns
import matplotlib.pylab as plt
import doubletdetection
import numpy as np

figures = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/figures/qc'
sc_file = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/single_cell_files'
os.makedirs(figures, exist_ok = True)
sc.set_figure_params(dpi = 300, format = 'png')
sc.settings.figdir = figures

if __name__ == '__main__':
    ##Convert gene_id to gene_name and add metadata from gtf file
    adata = sc.read(f'{sc_file}/multiome_gex_all_cells_raw.gz.h5ad')
    print(adata)
    adata.var_names_make_unique()
    adata = adata[:,adata.var['seqname']!='chrM']
    print(adata)
    sc.pp.calculate_qc_metrics(adata, expr_type='umis', percent_top=None, log1p=False, inplace=True)
    sc.pp.filter_cells(adata, min_counts=1)
    print(adata)
    adata.obs['log10_total_umis'] = np.log10(adata.obs['total_umis'])
    adata.obs['log10_n_genes_by_umis'] = np.log10(adata.obs['n_genes_by_umis'])
    sc.pl.scatter(adata, x='log10_total_umis', y='n_genes_by_umis',color='mouse', show=False, save='genes_by_umis_pretrim_log')
    sc.pl.scatter(adata, x='total_umis', y='n_genes_by_umis', color='mouse',show=False, save='genes_by_umis_pretrim_log')
    sc.pl.highest_expr_genes(adata, n_top=20, show=False, save=f'_pretrim')
    sc.pl.violin(adata, ['log10_total_umis'],groupby='mouse',show=False, save = 'umis_pretrim')
    sc.pl.violin(adata, ['n_genes_by_umis'], groupby='mouse',show=False,save = 'genes_pretrim')
    sc.pp.filter_cells(adata, min_counts=500)
    sc.pp.filter_cells(adata, min_genes=200)
    print(adata)
    sc.pp.filter_genes(adata, min_cells=10)
    print(adata)
    sc.pl.scatter(adata, x='log10_total_umis', y='n_genes_by_umis',color ='mouse', show = False, save = 'genes_by_umis_posttrim_log')
    sc.pl.scatter(adata, x='total_umis', y='n_genes_by_umis', color ='mouse',show = False, save = 'genes_by_umis_posttrim')
    sc.pl.highest_expr_genes(adata, n_top=20, show = False, save = f'_posttrim')
    sc.pl.violin(adata, ['log10_total_umis'],groupby='mouse',show=False, save='umis_posttrim')
    sc.pl.violin(adata, ['n_genes_by_umis'],groupby='mouse',show=False, save='genes_posttrim')
    # scrublet
    sce.pp.scrublet(adata, batch_key='mouse')
    sce.pl.scrublet_score_distribution(adata, show=False,save='scrublet_scores')
    print(adata)
    adata.write(f'{sc_file}/multiome_gex_qc_trimmed_raw.gz.h5ad', compression='gzip')
    sns.jointplot(
        data=adata.obs,
        x="log10_total_umis",
        y="log10_n_genes_by_umis",
        kind="hex",
    )
    plt.savefig(os.path.join(figures, 'genes_by_umis.png'))
