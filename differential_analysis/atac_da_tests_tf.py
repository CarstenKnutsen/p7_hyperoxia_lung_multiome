import pandas as pd
import os
import scanpy as sc
import numpy as np


# https://github.com/scverse/scanpy/issues/987
def obs_key_wise_subsampling(adata, obs_key, N):
    '''
    Subsample each class to same cell numbers (N). Classes are given by obs_key pointing to categorical in adata.obs.
    '''
    counts = adata.obs[obs_key].value_counts()
    # subsample indices per group defined by obs_key
    indices = [np.random.choice(adata.obs_names[adata.obs[obs_key]==group], size=N, replace=True) for group in counts.index]
    selection = np.hstack(np.array(indices))
    return adata[selection].copy()

figures = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/figures/atac'
sc_file = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/single_cell_files/share'
deg_dir = f"{figures}/datf"
os.makedirs(deg_dir, exist_ok=True)
sc.set_figure_params(dpi=200, format="png")
sc.settings.figdir = figures

if __name__ == "__main__":
    adata_tf = sc.read(f'{sc_file}/p7_multiome_tf_processed.gz.h5ad')
    sc.pp.normalize_total(adata_tf, target_sum=1e6)
    sc.pp.log1p(adata_tf,base=10)
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
    hyp_deg_dict = {}
    for lineage in adata_tf.obs["lineage"].unique():
        print(lineage)
        figures_lin = f"{deg_dir}/{lineage}"
        os.makedirs(figures_lin, exist_ok=True)
        sc.settings.figdir = figures_lin
        lin_adata_tf = adata_tf[adata_tf.obs["lineage"] == lineage]
        lin_adata_tf_norm = lin_adata_tf[lin_adata_tf.obs["treatment"] == 'Normoxia']
        ## broad ct markers
        for treat in lin_adata_tf.obs['treatment'].cat.categories:
            treat_lin_adata = lin_adata_tf[lin_adata_tf.obs['treatment']==treat]
            subsampled = obs_key_wise_subsampling(treat_lin_adata, 'celltype',500)

            sc.tl.rank_genes_groups(
                subsampled,
                "celltype",
                method="wilcoxon",
                pts=True,
                key_added="rank_peaks_groups_celltype",
            )
            with pd.ExcelWriter(
                    f"{figures_lin}/{lineage}_{treat.lower()}_celltype_markers.xlsx", engine="xlsxwriter"
            ) as writer:
                for ct in subsampled.obs["celltype"].cat.categories:
                    df= sc.get.rank_genes_groups_df(
                        subsampled, key="rank_peaks_groups_celltype", group=ct
                    )
                    df.index = df['names']
                    df.to_excel(writer, sheet_name=f"{ct} v rest"[:31])
            ##  ct v every other ct
            figures_lin_treat_comp = f"{figures_lin}/cell_type_comparisons/{treat.lower()}"
            os.makedirs(figures_lin_treat_comp, exist_ok=True)
            for ct in treat_lin_adata.obs["celltype"].cat.categories:
                print(ct)
                with pd.ExcelWriter(
                        f"{figures_lin_treat_comp}/{ct}.xlsx", engine="xlsxwriter"
                ) as writer:
                    for ct2 in treat_lin_adata.obs["celltype"].cat.categories:
                        cts_adata_tf = treat_lin_adata[treat_lin_adata.obs["celltype"].isin([ct, ct2])]
                        try:
                            sc.tl.rank_genes_groups(
                                cts_adata_tf,
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
                    hyp_deg_dict[ct]=df
                except:
                    print(ct)
                    print('no hyperoxia comparison')
                    continue
    with pd.ExcelWriter(
            f"{deg_dir}/hyperoxia_datf_all.xlsx", engine="xlsxwriter"
    ) as writer:
        for ct in adata_tf.obs['celltype'].cat.categories:
            try:
                hyp_deg_dict[ct].to_excel(writer, sheet_name=f"{ct}")
            except:
                print(ct)
    print('DONE!')
