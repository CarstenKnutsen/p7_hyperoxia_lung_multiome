import os
import scanpy as sc
import scanpy.external as sce
import muon as mu
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
    mudata = mu.read(f'{sc_file}/multi_all_cells_raw.h5mu')
    print(mudata)
    adata = mudata.mod['rna']
    adata.var['seqname'] = [x[0] for x in adata.var['interval'].str.split(':')]
    adata = adata[:,adata.var['seqname']!='chrM']
    print(adata)
    sc.pp.calculate_qc_metrics(adata, expr_type='umis', percent_top=None, log1p=False, inplace=True)
    mu.pp.filter_obs(adata, 'total_umis', lambda x: x >= 1)
    print(adata)
    adata.obs['log10_total_umis'] = np.log10(adata.obs['total_umis'])
    adata.obs['log10_n_genes_by_umis'] = np.log10(adata.obs['n_genes_by_umis'])
    sc.pl.scatter(adata, x='log10_total_umis', y='n_genes_by_umis',color='mouse', show=False, save='genes_by_umis_pretrim_log')
    sc.pl.scatter(adata, x='total_umis', y='n_genes_by_umis', color='mouse',show=False, save='genes_by_umis_pretrim_log')
    sc.pl.highest_expr_genes(adata, n_top=20, show=False, save=f'_pretrim')
    sc.pl.violin(adata, ['log10_total_umis'],groupby='mouse',show=False, save = 'umis_pretrim')
    sc.pl.violin(adata, ['n_genes_by_umis'], groupby='mouse',show=False,save = 'genes_pretrim')
    mu.pp.filter_obs(adata, 'total_umis', lambda x: x >= 500)
    mu.pp.filter_obs(adata, 'n_genes_by_umis', lambda x: x >= 200)
    mu.pp.filter_var(adata, 'n_cells_by_umis', lambda x: x >= 10)
    sc.pl.scatter(adata, x='log10_total_umis', y='n_genes_by_umis',color ='mouse', show = False, save = 'genes_by_umis_posttrim_log')
    sc.pl.scatter(adata, x='total_umis', y='n_genes_by_umis', color ='mouse',show = False, save = 'genes_by_umis_posttrim')
    sc.pl.highest_expr_genes(adata, n_top=20, show = False, save = f'_posttrim')
    sc.pl.violin(adata, ['log10_total_umis'],groupby='mouse',show=False, save='umis_posttrim')
    sc.pl.violin(adata, ['n_genes_by_umis'],groupby='mouse',show=False, save='genes_posttrim')
    # scrublet
    sce.pp.scrublet(adata, batch_key='mouse')
    sce.pl.scrublet_score_distribution(adata, show=False,save='scrublet_scores')
    print(adata)
    print(mudata)
    adata = adata[adata.obs['predicted_doublet']==False]
    mudata.mod['rna'] = adata
    mudata.update()
    print(mudata)

    mudata.write(f'{sc_file}/multi_all_cells_raw.h5mu')
    sns.jointplot(
        data=adata.obs,
        x="log10_total_umis",
        y="log10_n_genes_by_umis",
        kind="hex",
    )
    plt.savefig(os.path.join(figures, 'genes_by_umis.png'))
