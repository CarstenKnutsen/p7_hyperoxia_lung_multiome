import numpy as np
import pandas as pd
import os
import scanpy as sc
import sys
import muon as mu
import muon.atac as ac
sys.path.insert(1, "/home/carsten/alvira_bioinformatics/lungsc_ck")
from ck_functions import anndataks_ck

figures = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/figures/atac'
sc_file = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/single_cell_files'
deg_dir = f"{figures}/datf"
os.makedirs(deg_dir, exist_ok=True)
sc.set_figure_params(dpi=200, format="png")
sc.settings.figdir = figures

if __name__ == "__main__":
    mudata = mu.read(f'{sc_file}/multi_all_cells_processed.h5mu')
    adata = mudata.mod['atac']
    adata.X = adata.layers['counts'].copy()
    # tf-bc matrix
    matrix_dir = "/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/cellranger_output/230609_aggregate/outs/analysis/tf_analysis/filtered_tf_bc_matrix"
    adata_tf = sc.read_mtx(os.path.join(matrix_dir, "matrix.mtx"))

    motifs_path = os.path.join(matrix_dir, "motifs.tsv")
    var = pd.read_csv(motifs_path, sep='\t', index_col=0, header=None)
    barcodes_path = os.path.join(matrix_dir, "barcodes.tsv")
    obs = pd.read_csv(barcodes_path, sep='\t', index_col=0, header=None)
    adata_tf = adata_tf.T
    adata_tf.layers['raw'] = adata_tf.X.copy()
    adata_tf.obs_names = obs.index.values
    adata_tf.var_names = var.index.values
    adata_tf = adata_tf[adata.obs_names,:]
    sc.pp.normalize_total(adata_tf, target_sum=1e6)
    sc.pp.log1p(adata_tf,base=10)
    adata_tf.obs[['lineage','celltype','treatment']] = adata.obs[['lineage','celltype','treatment']]
    adata_tf_norm = adata_tf[adata_tf.obs['treatment']=='Normoxia']
    sc.tl.rank_genes_groups(
        adata_tf_norm,
        "lineage",
        method="wilcoxon",
        pts=True,
        key_added="rank_peaks_groups_lineage",
    )
    with pd.ExcelWriter(
            f"{deg_dir}/lineage_markers.xlsx", engine="xlsxwriter"
    ) as writer:
        for ct in adata_tf.obs["lineage"].cat.categories:
            df = sc.get.rank_genes_groups_df(
                adata_tf_norm, key="rank_peaks_groups_lineage", group=ct
            )
            df.index = df['names']
            df.to_excel(writer, sheet_name=f"{ct} v rest"[:31])


    print(adata_tf)
    for lineage in adata_tf.obs["lineage"].unique():
        print(lineage)
        figures_lin = f"{deg_dir}/{lineage}"
        os.makedirs(figures_lin, exist_ok=True)
        sc.settings.figdir = figures_lin
        lin_adata_tf = adata_tf[adata_tf.obs["lineage"] == lineage]
        lin_adata_tf_norm = lin_adata_tf[lin_adata_tf.obs["treatment"] == 'Normoxia']
        ## broad ct markers
        sc.tl.rank_genes_groups(
            lin_adata_tf_norm,
            "celltype",
            method="wilcoxon",
            pts=True,
            key_added="rank_peaks_groups_celltype",
        )
        with pd.ExcelWriter(
                f"{figures_lin}/{lineage}_normoxia_celltype_markers.xlsx", engine="xlsxwriter"
        ) as writer:
            for ct in lin_adata_tf.obs["celltype"].cat.categories:
                df= sc.get.rank_genes_groups_df(
                    lin_adata_tf_norm, key="rank_peaks_groups_celltype", group=ct
                )
                df.index = df['names']
                df.to_excel(writer, sheet_name=f"{ct} v rest"[:31])
        ##  ct v every other ct
        figures_lin_comp = f"{figures_lin}/cell_type_comparisons"
        os.makedirs(figures_lin_comp, exist_ok=True)
        for ct in lin_adata_tf_norm.obs["celltype"].cat.categories:
            print(ct)
            with pd.ExcelWriter(
                    f"{figures_lin_comp}/{ct}.xlsx", engine="xlsxwriter"
            ) as writer:
                for ct2 in lin_adata_tf_norm.obs["celltype"].cat.categories:
                    cts_adata_tf = lin_adata_tf_norm[lin_adata_tf_norm.obs["celltype"].isin([ct, ct2])]
                    try:
                        sc.tl.rank_genes_groups(
                            lin_adata_tf_norm,
                            "celltype",
                            groups=[ct, ct2],
                            method="wilcoxon",
                            pts=True,
                            key_added="rank_peaks_groups_celltype",
                        )
                        df = sc.get.rank_genes_groups_df(
                            cts_adata_tf, key="rank_peaks_groups_celltype", group=ct
                        )
                        df.index = df['names']
                        df.to_excel(writer, sheet_name=f"{ct} v {ct2}"[:31])
                    except:
                        print(ct)
                        print(ct2)
                        print('no comp')

        ## treatment
        with pd.ExcelWriter(
                f"{figures_lin}/hyperoxia_datf.xlsx", engine="xlsxwriter"
        ) as writer:
            for ct in lin_adata_tf.obs["celltype"].cat.categories:
                ct_adata_tf = lin_adata_tf[lin_adata_tf.obs["celltype"] == ct]
                try:
                    sc.tl.rank_genes_groups(
                        ct_adata_tf,
                        "treatment",
                        method="wilcoxon",
                        pts=True,
                        key_added="rank_peaks_groups_treatment",
                    )
                    df = sc.get.rank_genes_groups_df(
                        ct_adata_tf,
                        key="rank_peaks_groups_treatment",
                        group="Hyperoxia",
                    )
                    df.index = df['names']
                    df.to_excel(writer, sheet_name=f"{ct}")
                except:
                    print(ct)
                    print('no hyperoxia comparison')
                    continue

    print('DONE!')
