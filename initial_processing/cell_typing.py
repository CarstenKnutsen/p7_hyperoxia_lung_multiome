import pandas as pd
import os
import scanpy as sc
import scanpy.external as sce
import seaborn as sns
import numpy as np
import matplotlib.pylab as plt
import leidenalg as la
import northstar
import itertools
import anndata

figures = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/figures/tissue_embedding'
sc_file = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/single_cell_files'
os.makedirs(figures, exist_ok=True)
sc.set_figure_params(dpi=300, dpi_save=300, format="png")
sc.settings.figdir = figures

if __name__ == "__main__":
    # read in processed data
    adata = sc.read(f"{sc_file}/multiome_gex_processed_log_no_doublet.gz.h5ad")
    print(adata)
    # marker gene dicts for each lineage
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
            "Mki67",
            "Col1a1",
            "Epcam",
            "Ptprc",
            "Pecam1",
        ],
    }
    # create empty cell type columns for both northstar and cluster identities
    adata.obs["celltype"] = pd.Series(index=adata.obs.index, data=None, dtype="str")
    # do this for every lineage
    for lineage in adata.obs["lineage"].unique():
        print(lineage)
        genes = gene_dict[lineage]
        figures_lin = f"{figures}/{lineage}"
        os.makedirs(figures_lin, exist_ok=True)
        sc.settings.figdir = figures_lin
        lin_adata = anndata.AnnData(
            adata[adata.obs["lineage"] == lineage].X,
            obs=adata[adata.obs["lineage"] == lineage].obs,
            var=adata[adata.obs["lineage"] == lineage].var,
        )
        del lin_adata.var["highly_variable"]

        sc.pp.highly_variable_genes(lin_adata, n_top_genes=2000, batch_key="mouse")
        sc.pl.highly_variable_genes(lin_adata, save="highly_variable", show=False)
        sc.pp.pca(lin_adata, use_highly_variable=True)
        sc.pp.neighbors(lin_adata)
        # print(lineage)
        # print(lin_adata.obs["Individual"].value_counts())
        # sce.pp.bbknn(
        #     lin_adata,
        #     batch_key="Run",
        #     n_pcs=10,
        #     neighbors_within_batch=2,
        #     set_op_mix_ratio=0.5,
        # )
        # lin_adata.uns["neighbors"]["connectivities"] = lin_adata.obsp["connectivities"]
        leiden_res_dict = {
            "mesenchymal": 0.001,
            "endothelial": 0.001,
            "epithelial": 0.001,
            "immune": 0.002,
        }
        sc.tl.leiden(
            lin_adata,
            resolution=leiden_res_dict[lineage],
            partition_type=la.CPMVertexPartition,
            key_added=f"leiden_{lineage}",
        )
        # # drop low quality clusters
        # if lineage == "endothelial":
        #     lin_adata = lin_adata[
        #         lin_adata.obs[f"leiden_{lineage}"] != '8'
        #     ]  # low quality cells with no endo markers, lots of Rpl/Rps and some epi markers
        # elif lineage == 'mesenchymal':
        #     lin_adata = lin_adata[
        #         lin_adata.obs[f"leiden"] != '52'
        #         ]  # Low quality E20 cells after iteration, fall in center of umap and have no markers outside of ribosomal genes
        # # leiden cluster identity based off of cell type markers
        leiden_ct_dict = {
            "mesenchymal": {
                "0": "Alveolar fibroblast",
                "1": "Alveolar fibroblast",
                "2": "Alveolar fibroblast",
                "3": "Myofibroblast",
                "4": "Alveolar fibroblast",
                "5": "Mesothelial",
                "6": "Airway smooth muscle",
                "7": "Pericyte",
                "8": "Alveolar fibroblast",
                "9": "Adventitial fibroblast",
                "10": "Alveolar fibroblast",
                "11": "Proliferating fibroblast",
                "12": "Proliferating myofibroblast",
                "13": "Pericyte",
                "14": "Alveolar fibroblast",
                "15": "Alveolar fibroblast",
                "16": "Proliferating pericyte",
                "17": "Acta1+ weird cells",
                "18": "Vascular smooth muscle",
                "19": "Chrm2+ hyperoxia cells",
                "20": "Myofibroblast",
                "21": "Alveolar fibroblast",
                "22": "Alveolar fibroblast",
                "23": "Alveolar fibroblast",
                "24": "Alveolar fibroblast",
                "25": "Alveolar fibroblast",
                "26": "Alveolar fibroblast",
            },
            "endothelial": {
                "0": "Cap1",
                "1": "Cap1",
                "2": "Cap2",
                "3": "Cap2",
                "4": "Proliferating EC",
                "5": "Cap1",
                "6": "Arterial EC",
                "7": "Venous EC",
                "8": "Lymphatic EC",
                "9": "Cap1",
            },
            "immune": {
                "0": "Alveolar macrophage",
                "1": "Alveolar macrophage",
                "2": "Alveolar macrophage",
                "3": "Monocyte",
                "4": "Proliferating macrophage",
                "5": "T cell",
                "6": "B cell",
                "7": "Interstitial macrophage",
                "8": "Intermediate monocyte",
                "9": "Dendritic cell",
                "10": "NK cell",
                "11": "IL cell",
                "12": "Alveolar macrophage",
                "13": "Alveolar macrophage",

            },
            "epithelial": {
                "0": "AT2",
                "1": "AT2",
                "2": "AT1",
                "3": "AT2",
                "4": "AT1",
                "5": "Club",
                "6": "Club",
                "7": "Ciliated",
                "8": "Proliferating AT2",
                "9": "Basal",
                "10": "AT2",
                "11": "Club",
                "12": "Club",
                "13": "Neuroendocrine",
                "14": "Ciliated",
                "15": "AT1_AT2",

            },
        }
        lin_adata.obs["celltype"] = (
            lin_adata.obs[f"leiden_{lineage}"].map(leiden_ct_dict[lineage]).astype("str")
        )
        sc.tl.dendrogram(lin_adata, groupby=f"leiden_{lineage}")
        sc.tl.rank_genes_groups(lin_adata, f"leiden_{lineage}", method="wilcoxon")
        print(lin_adata.obs[f"leiden_{lineage}"].cat.categories)
        sc.pl.rank_genes_groups_dotplot(
            lin_adata,
            groupby=f"leiden_{lineage}",
            n_genes=int(150 / len(lin_adata.obs[f"leiden_{lineage}"].unique())),
            show=False,
            save=f"{lineage}_leiden_markers.png",
        )
        sc.pl.DotPlot(lin_adata, genes, groupby=f"leiden_{lineage}",).style(
            cmap="Reds"
        ).add_totals().savefig(
            os.path.join(figures_lin, f"dotplot_{lineage}_leiden.png"), dpi=300
        )
        plt.clf()
        adata.obs["celltype"].loc[lin_adata.obs.index] = lin_adata.obs[
            "celltype"
        ]
        # Visualize cell type connections, generate umap for lineage, save lineage umap to main object
        sc.tl.paga(lin_adata, groups="celltype")
        sc.pl.paga(
            lin_adata,
            color="celltype",
            show=False,
            save=f"_{lineage}_celltype.png",
        )
        b_dict = {"mesenchymal": 1.4, "endothelial": 1.1, "epithelial": 1.6, "immune": 0.9}
        sc.tl.umap(lin_adata, a=2, b=b_dict[lineage])
        # Add lineage umaps and leiden clusters to top level
        adata.obs[f"umap_{lineage}_1"] = np.nan
        adata.obs[f"umap_{lineage}_2"] = np.nan
        lin_adata.obs[f"umap_{lineage}_1"] = [x[0] for x in lin_adata.obsm["X_umap"]]
        lin_adata.obs[f"umap_{lineage}_2"] = [x[1] for x in lin_adata.obsm["X_umap"]]
        adata.obs[f"umap_{lineage}_1"].loc[lin_adata.obs.index] = lin_adata.obs[
            f"umap_{lineage}_1"
        ]
        adata.obs[f"umap_{lineage}_2"].loc[lin_adata.obs.index] = lin_adata.obs[
            f"umap_{lineage}_2"
        ]
        adata.obs[f"leiden_{lineage}"] = np.nan
        adata.obs[f"leiden_{lineage}"].loc[lin_adata.obs.index] = lin_adata.obs[
            f"leiden_{lineage}"
        ]
        adata.obsm[f"X_umap_{lineage}"] = adata.obs[
            [f"umap_{lineage}_1", f"umap_{lineage}_2"]
        ].to_numpy()
        del adata.obs[f"umap_{lineage}_1"]
        del adata.obs[f"umap_{lineage}_2"]
        sc.pl.pca_overview(lin_adata, color=f"leiden_{lineage}", show=False, save=True)
        sc.pl.pca_variance_ratio(lin_adata, show=False, save="variance_ratio")
        sc.pl.pca_loadings(
            lin_adata, components="1,2,3,4,5", show=False, save="loadings"
        )
        plt.close()
        plt.clf()
        metadata = f"{figures_lin}/metadata"
        os.makedirs(metadata, exist_ok=True)
        sc.settings.figdir = metadata
        for color in [
            "leiden",
            f"leiden_{lineage}",
            "lineage",
            "treatment",
            "mouse",
            "n_genes_by_umis",
            "celltype",
        ]:
            if color in ["leiden",f"leiden_{lineage}"]:
                sc.pl.umap(
                    lin_adata,
                    color=color,
                    alpha=0.5,
                    legend_loc="on data",
                    show=False,
                    save=f"{color}",
                )
            else:
                sc.pl.umap(
                    lin_adata, color=color, alpha=0.5, show=False, save=f"{color}"
                )
        gene_fn = f"{figures_lin}/genes"
        os.makedirs(gene_fn, exist_ok=True)
        sc.settings.figdir = gene_fn
        for gene in genes:
            sc.pl.umap(lin_adata, alpha=0.5, color=gene, show=False, save=f"{gene}")
        sc.settings.figdir = figures_lin
        sc.pl.DotPlot(lin_adata, genes, groupby="celltype",).add_dendrogram(
            show=True
        ).style(cmap="Reds").add_totals().savefig(
            os.path.join(figures_lin, f"dotplot_{lineage}_celltype.png"), dpi=300
        )
        plt.clf()

        pd.DataFrame(index=lin_adata.obs_names, data=lin_adata.obsm["X_umap"]).to_csv(
            f"{figures_lin}/{lineage}_umap.csv"
        )
        # lin_adata = lin_adata[~lin_adata.obs["celltype"].isna(), :]
        # lin_adata.obs = lin_adata.obs.dropna(axis=1, how="all")
        # lin_adata.write(f"{sc_file}/{lineage}_processed_log.gz.h5ad", compression="gzip")
        with pd.ExcelWriter(
            f"{figures_lin}/{lineage}_metadata_counts.xlsx", engine="xlsxwriter"
        ) as writer:
            obs_list = ["celltype", "treatment", "mouse"]
            num_obs = len(obs_list) + 1
            for ind in range(0, num_obs):
                for subset in itertools.combinations(obs_list, ind):
                    if len(subset) != 0:
                        subset = list(subset)
                        if len(subset) == 1:
                            key = subset[0]
                            lin_adata.obs[key].value_counts().to_excel(
                                writer, sheet_name=key
                            )
                        else:
                            key = "_".join(subset)
                            lin_adata.obs.groupby(subset[:-1])[
                                subset[-1]
                            ].value_counts().to_excel(writer, sheet_name=key[:31])
        sc.tl.rank_genes_groups(lin_adata, "celltype", method="wilcoxon")
        sc.pl.rank_genes_groups_dotplot(
            lin_adata,
            groupby="celltype",
            n_genes=int(100 / len(lin_adata.obs["leiden"].unique())),
            show=False,
            save=f"_{lineage}_celltype_markers.png",
        )
    print('Done witn lineages')
    adata.X = adata.layers["raw"].copy()
    del adata.layers
    # adata = adata[~adata.obs['celltype'].isna()] # drop the low quality clusters from the lineage portion
    print('trimmed')
    print(adata)
    adata.write(f"{sc_file}/multiome_gex_processed_cell_typed_raw.gz.h5ad", compression="gzip")
    sc.pp.normalize_total(adata, target_sum = 1e4)
    sc.pp.log1p(adata)
    sc.settings.figdir = f"{figures}/tissue_embedding"


    gene_ls = []
    for k in gene_dict.keys():
        gene_ls = gene_ls + gene_dict[k]
    gene_ls = gene_ls + [
                 'leiden',
                 'lineage',
                 'log10_n_genes_by_umis',
                 'mouse',
                 'treatment',
                 'log10_total_umis',

        'celltype'
                 ]
    for gene in gene_ls:
        if gene == 'leiden':
            sc.pl.umap(adata, legend_loc='on data', color=gene, alpha=0.5, show = False,save=f'_{gene}')
        else:
            sc.pl.umap(adata, color=gene, alpha=0.3, show = False, save=f'_{gene}')
    with pd.ExcelWriter(f'{figures}/tissue_embedding/metadata_counts.xlsx', engine='xlsxwriter') as writer:
        obs_list = ['lineage', 'treatment','mouse']
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


