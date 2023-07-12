import pandas as pd
import os
import scanpy as sc
import muon as mu
import itertools
import scanpy.external as sce

figures = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/figures/tissue_embedding'
sc_file = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/single_cell_files'
os.makedirs(figures, exist_ok = True)
sc.set_figure_params(dpi = 300, dpi_save = 300, format = 'png')
sc.settings.figdir = figures
if __name__ == '__main__':
    mudata = mu.read(f'{sc_file}/multi_all_cells_raw.h5mu')
    adata = mudata.mod['rna']
    adata.layers['soupx'] = adata.X.copy()
    sc.pp.normalize_total(adata, key_added=None, target_sum=1e6)
    sc.pp.log1p(adata, base=10)
    sc.pp.highly_variable_genes(adata,
                                n_top_genes = 2000,
                                batch_key='mouse'
                                )
    sc.pp.pca(adata, use_highly_variable=True)
    sce.pp.harmony_integrate(adata, key='mouse', max_iter_harmony=20)
    sc.pp.neighbors(adata, use_rep='X_pca_harmony')
    sc.tl.umap(adata, min_dist=0.5)
    sc.tl.leiden(adata, key_added='leiden')
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups_dotplot(adata, groupby='leiden',
                                    n_genes=int(150 / len(adata.obs['leiden'].unique())), show=False,
                                    save=f'leiden_markers.png')
    # adata = adata[~adata.obs['leiden'].isin(['14']),:] # Cluster 14 looks low quality
    sc.tl.umap(adata, min_dist=0.5)
    sc.pl.pca_overview(adata, color='leiden', show=False, save = True)
    sc.pl.pca_variance_ratio(adata, show = False, save ='variance_ratio')
    sc.pl.pca_loadings(adata, components = ','.join([str(x) for x in range(1,10)]), show = False, save = True)
    ## Assign lineage from gene expression
    genes = ['Col1a1', 'Cdh5', 'Ptprc', 'Epcam', 'leiden']
    genedf = sc.get.obs_df(adata, keys=genes)
    grouped = genedf.groupby("leiden")
    mean = grouped.mean()
    mean_t = mean.T
    mean_t.to_csv(f'{figures}/leiden_lineage_scores.csv')
    lineage_dict = {}
    for cluster in mean_t.columns:
        gene = mean_t[cluster].idxmax()
        if gene == 'Cdh5':
            lineage_dict[cluster] = 'endothelial'
        elif gene == 'Ptprc':
            lineage_dict[cluster] = 'immune'
        elif gene == 'Epcam':
            lineage_dict[cluster] = 'epithelial'
        elif gene == 'Col1a1':
            lineage_dict[cluster] = 'mesenchymal'
    adata.obs['lineage'] = [lineage_dict[x] for x in adata.obs['leiden']]
    mudata.write(f'{sc_file}/multi_all_cells_raw.h5mu')
    adata.obs['predicted_doublet'] = adata.obs['predicted_doublet'].astype('str')
    print(adata)
    gene_dict = {
        "mesenchymal": [
            "Col3a1",
            "G0s2",
            "Limch1",
            "Col13a1",
            "Col14a1",
            "Serpinf1",
            "Pdgfra",
            "Scara5",
            "Acta2",
            "Hhip",
            "Fgf18",
            "Wif1",
            "Tgfbi",
            "Tagln",
            "Mustn1",
            "Aard",
            "Pdgfrb",
            "Cox4i2",
            "Higd1b",
            "Wt1",
            "Lrrn4",
            "Upk3b",
            "Mki67",
            "Acta1",
            'Lgals3',
            'Tubb3',
            'Aqp3',
            "Epcam",
            "Ptprc",
            "Pecam1",
        ],
        "endothelial": [
            "Gja5",
            "Bmx",
            "Fn1",
            "Ctsh",
            "Kcne3",
            "Cdh13",
            "Car8",
            "Mmp16",
            "Slc6a2",
            "Thy1",
            "Mmrn1",
            "Ccl21a",
            "Reln",
            "Neil3",
            "Mki67",
            "Aurkb",
            "Depp1",
            "Ankrd37",
            "Peg3",
            "Mest",
            "Hpgd",
            "Cd36",
            "Car4",
            "Sirpa",
            "Fibin",
            "Col1a1",
            "Epcam",
            "Ptprc",
            "Pecam1",
        ],
        "immune": [
            "Cd68",
            "Gal",
            "Itgax",
            "Car4",
            "C1qa",
            "Plac8",
            "Batf3",
            "Itgae",
            # "Cd209a",
            "Mreg",
            "Mcpt8",
            "Retnlg",
            "Ms4a1",
            "Gzma",
            "Cd3e",
            "Areg",
            "Mki67",
            "Col1a1",
            "Epcam",
            "Ptprc",
            "Pecam1",
        ],
        "epithelial": [
            'Muc1',
            "Scg5",
            "Ascl1",
            "Lyz1",
            "Lyz2",
            "Sftpc",
            "Slc34a2",
            "S100g",
            "Sftpa1",
            "Akap5",
            "Hopx",
            "Col4a4",
            "Vegfa",
            "Lyz1",
            "Tmem212",
            "Dynlrb2",
            "Cdkn1c",
            "Tppp3",
            "Scgb3a2",
            "Cyp2f2",
            "Scgb1a1",
            "Reg3g",
            "Scgb3a1",
            "Mki67",
            "Col1a1",
            "Epcam",
            "Ptprc",
            "Pecam1",
        ],
    }
    gene_ls = []
    for k in gene_dict.keys():
        gene_ls = gene_ls + gene_dict[k]
    gene_ls = gene_ls + [
                 'leiden',
                 'lineage',
                 'n_genes_by_umis',
                 'mouse',
                 'treatment',
                 'total_umis',
                 'timepoint',
                 'predicted_doublet',
                 'doublet_score',
                 ]
    for gene in gene_ls :
        if gene == 'leiden':
            sc.pl.umap(adata, legend_loc='on data', color=gene, alpha=0.5, show = False,save=f'_{gene}')
        else:
            sc.pl.umap(adata, color=gene, alpha=0.5, show = False, save=f'_{gene}')
    pd.DataFrame(index=adata.obs_names, data = adata.obsm['X_umap']).to_csv(f'{figures}/allcells_umap.csv')
    with pd.ExcelWriter(f'{figures}/metadata_counts.xlsx', engine='xlsxwriter') as writer:
        obs_list = ['lineage', 'mouse',]
        num_obs = len(obs_list) + 1
        for ind in range(0, num_obs):
            for subset in itertools.combinations(obs_list, ind):
                if len(subset) != 0:
                    subset = list(subset)
                    if len(subset) == 1:
                        key = subset[0]
                        adata.obs[key].value_counts().to_excel(writer, sheet_name=key)
                    else:
                        key = "_".join(subset)
                        adata.obs.groupby(subset[:-1])[subset[-1]].value_counts().to_excel(writer, sheet_name=key[:31])


