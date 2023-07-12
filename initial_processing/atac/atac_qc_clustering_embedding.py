import numpy as np
import pandas as pd
import os
import scanpy as sc
import scanpy.external as sce
import muon as mu
from muon import atac as ac
import anndata
import matplotlib.pyplot as plt

figures = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/figures/atac'
sc_file = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/single_cell_files'
sc.set_figure_params(dpi=300, format="png")
sc.settings.figdir = figures
gene_dict = {
    "mesenchymal": [
        "Col3a1",
        "G0s2",
        "Limch1",
        "Col13a1",
        "Col14a1",
        "Serpinf1",
        'Dcn',
        "Pdgfra",
        "Scara5",
        "Acta2",
        "Hhip",
        "Fgf18",
        "Wif1",
        "Tgfbi",
        "Tagln",
        "Mustn1",
        'Msln',
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
        # "Klra5",
        'Klrc2',
        'Il5',
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
        'Dmbt1',
        'Aqp3',
        'Krt5',
        'Scgb1a1',
        'Krt15',
        'Foxj1',
        'Ccdc153',
        # 'Dcllk1',
        # 'Ascl3',
        'Gp2',
        "Mki67",
        "Col1a1",
        "Epcam",
        "Ptprc",
        "Pecam1",
    ],
}


def make_dir(directory):
    os.makedirs(directory, exist_ok=True)
    sc.settings.figdir = directory
    # ac.settings.figdir = directory
    return


if __name__ == "__main__":
    ## read in data and add obs to atac data

    mdata = mu.read(f'{sc_file}/multi_all_cells_raw.h5mu')
    ac.tl.add_peak_annotation_gene_names(mdata)
    mu.pp.intersect_obs(mdata)
    mdata.mod['atac'].obs[['lineage', 'celltype']] = mdata.mod['rna'].obs[['lineage', 'celltype']]
    atac = mdata.mod['atac']

    #add tf_annotations
    tf_bed = pd.read_csv(
        "/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/cellranger_output/230609_aggregate/outs/analysis/tf_analysis/peak_motif_mapping.bed",
        sep='\t',
        index_col=None,
        header=None
        )
    tf_bed['peak_name'] = tf_bed[0] + ":" + tf_bed[1].astype('string') + '-' + tf_bed[2].astype('string')
    df = tf_bed.groupby('peak_name')[3].apply(lambda x: '&'.join(x.astype(str))).reset_index()
    peak_tf_dt = pd.Series(df[3].values, index=df.peak_name)
    atac.var['tfs'] = pd.Series(peak_tf_dt)

    # add gene to var
    df = atac.uns['atac']['peak_annotation'].copy()
    df['gene'] = df.index.copy()
    df2 = df.groupby('peak')['gene'].apply(lambda x: ','.join(x.astype(str))).reset_index()
    peak_gene_dt = pd.Series(df2['gene'].values, index=df2.peak)
    atac.var['annotated_gene'] = pd.Series(peak_gene_dt)
    atac.var['annotated_gene']

    ## temp to add peak_type
    df = atac.uns['atac']['peak_annotation'].copy()
    print(df.columns)
    df['gene'] = df.index.copy()
    df2 = df.groupby('peak')['peak_type'].apply(lambda x: ','.join(x.astype(str))).reset_index()
    atac.var['peak_type'] = pd.Series(df2['peak_type'].values, index=df2.peak)
    df2 = df.groupby('peak')['distance'].apply(lambda x: ','.join(x.astype(str))).reset_index()
    atac.var['distance'] = pd.Series(df2['distance'].values, index=df2.peak)




    ## QC

    sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True)
    make_dir(f'{figures}/qc')
    sc.pl.violin(atac, ['total_counts', 'n_genes_by_counts'], jitter=0.4, multi_panel=True,show=False, save='_counts_pretrim.png')
    mu.pp.filter_obs(atac, 'total_counts', lambda x: (x >= 1000))
    sc.pl.violin(atac, ['total_counts', 'n_genes_by_counts'], jitter=0.4, multi_panel=True,show=False, save='_counts_posttrim.png')
    atac.layers["counts"] = atac.X.copy()

    ## normalize
    # mu.pl.histogram(atac, ['n_genes_by_counts', 'total_counts'],save='counts_psttrim.png')
    ac.tl.nucleosome_signal(atac, n=1e6)
    # mu.pl.histogram(atac, "nucleosome_signal", kde=True,save='nucleosome_signal.png')
    ac.tl.get_gene_annotation_from_rna(mdata['rna'])
    tss = ac.tl.tss_enrichment(mdata, n_tss=1000)  # by default, features=ac.tl.get_gene_annotation_from_rna(mdata)
    # ac.pl.tss_enrichment(tss, save='tss_enrichment.png')
    ac.pp.binarize(atac)
    ac.pp.tfidf(atac, scale_factor=1e4)
    sc.tl.pca(atac)
    sc.pl.pca_variance_ratio(atac, show=False, save="variance_ratio")
    sce.pp.harmony_integrate(atac, 'mouse')
    sc.pp.neighbors(atac, use_rep='X_pca_harmony',n_pcs=10)
    sc.tl.leiden(atac)
    sc.tl.umap(atac,min_dist=0.1)
    print(atac)
    make_dir(f'{figures}/embedding/tissue_embedding')
    color_ls = ['leiden','treatment','mouse','celltype','lineage']
    for color in color_ls:
        sc.pl.umap(atac, color=color,alpha=0.5,show=False,save=f'_{color}.png')
    for lineage in atac.obs['lineage'].cat.categories:
        figures_lin = f'{figures}/embedding/{lineage}'
        genes=gene_dict[lineage]
        make_dir(figures_lin)
        lin_atac = atac[atac.obs['lineage']==lineage].copy()
        # lin_atac.X = lin_atac.layers['counts'].copy()
        # ac.pp.binarize(lin_atac)
        # ac.pp.tfidf(lin_atac, scale_factor=1e4)
        sc.pp.pca(lin_atac)
        sce.pp.harmony_integrate(lin_atac, key='mouse', max_iter_harmony=20)
        sc.pp.neighbors(lin_atac,use_rep='X_pca_harmony',n_pcs=10)
        sc.tl.leiden(
            lin_atac,
            key_added=f"leiden_{lineage}",
        )

        sc.tl.dendrogram(lin_atac, groupby=f"leiden_{lineage}")
        sc.tl.rank_genes_groups(lin_atac, f"leiden_{lineage}", method="wilcoxon")
        print(lin_atac.obs[f"leiden_{lineage}"].cat.categories)
        sc.pl.rank_genes_groups_dotplot(
            lin_atac,
            groupby=f"leiden_{lineage}",
            n_genes=int(150 / len(lin_atac.obs[f"leiden_{lineage}"].unique())),
            show=False,
            save=f"{lineage}_leiden_markers.png",
        )
        # Visualize cell type connections, generate umap for lineage, save lineage umap to main object
        sc.tl.paga(lin_atac, groups="celltype")
        sc.pl.paga(
            lin_atac,
            color="celltype",
            show=False,
            save=f"_{lineage}_celltype.png",
        )
        sc.tl.umap(lin_atac, min_dist=0.1)
        # Add lineage umaps and leiden clusters to top level
        atac.obs[f"umap_{lineage}_1"] = np.nan
        atac.obs[f"umap_{lineage}_2"] = np.nan
        lin_atac.obs[f"umap_{lineage}_1"] = [x[0] for x in lin_atac.obsm["X_umap"]]
        lin_atac.obs[f"umap_{lineage}_2"] = [x[1] for x in lin_atac.obsm["X_umap"]]
        atac.obs[f"umap_{lineage}_1"].loc[lin_atac.obs.index] = lin_atac.obs[
            f"umap_{lineage}_1"
        ]
        atac.obs[f"umap_{lineage}_2"].loc[lin_atac.obs.index] = lin_atac.obs[
            f"umap_{lineage}_2"
        ]
        atac.obs[f"leiden_{lineage}"] = np.nan
        atac.obs[f"leiden_{lineage}"].loc[lin_atac.obs.index] = lin_atac.obs[
            f"leiden_{lineage}"
        ]
        atac.obsm[f"X_umap_{lineage}"] = atac.obs[
            [f"umap_{lineage}_1", f"umap_{lineage}_2"]
        ].to_numpy()
        del atac.obs[f"umap_{lineage}_1"]
        del atac.obs[f"umap_{lineage}_2"]
        sc.pl.pca_overview(lin_atac, color=f"leiden_{lineage}", show=False, save=True)
        sc.pl.pca_variance_ratio(lin_atac, show=False, save="variance_ratio")
        sc.pl.pca_loadings(
            lin_atac, components="1,2,3,4,5", show=False, save="loadings"
        )
        plt.close()
        plt.clf()
        metadata = f"{figures_lin}/metadata"
        os.makedirs(metadata, exist_ok=True)
        sc.settings.figdir = metadata
        print(lineage)

        print(lin_atac)
        for color in [
            # "leiden",
            f"leiden_{lineage}",
            "treatment",
            "mouse",
            "n_genes_by_counts",
            "celltype",
        ]:
            if color in ["leiden", f"leiden_{lineage}"]:
                sc.pl.umap(
                    lin_atac,
                    color=color,
                    alpha=0.5,
                    legend_loc="on data",
                    show=False,
                    save=f"{color}",
                )
            else:
                sc.pl.umap(
                    lin_atac, color=color, alpha=0.5, show=False, save=f"{color}"
                )
        gene_fn = f"{figures_lin}/genes"
        os.makedirs(gene_fn, exist_ok=True)
        sc.settings.figdir = gene_fn
        print(lin_atac.uns['atac']['peak_annotation'])
        genes = [gene for gene in genes if gene in lin_atac.uns['atac']['peak_annotation'].index.unique()]
        for gene in genes:
            ac.pl.umap(lin_atac, color=[gene], average='peak_type',use_raw=False,show=False, save=f"{gene}")
        sc.settings.figdir = figures_lin
        ac.pl.dotplot(lin_atac, genes, groupby='celltype',show=False,save=f'{lineage}_genes.png')
        sc.tl.rank_genes_groups(lin_atac, "celltype", method="wilcoxon")
        sc.pl.rank_genes_groups_dotplot(
            lin_atac,
            groupby="celltype",
            n_genes=int(100 / len(lin_atac.obs["leiden"].unique())),
            show=False,
            save=f"_{lineage}_celltype_markers.png",
        )
    print('Done witn lineages')
    mdata.write(f'{sc_file}/multi_all_cells_processed.h5mu')
    print(mdata)
        











