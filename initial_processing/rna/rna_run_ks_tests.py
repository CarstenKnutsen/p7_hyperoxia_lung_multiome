import numpy as np
import pandas as pd
import os
import scanpy as sc
import sys
import muon as mu

sys.path.insert(1, "/home/carsten/alvira_bioinformatics/lungsc_ck")
from ck_functions import anndataks_ck

figures = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/figures/tissue_embedding'
sc_file = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/single_cell_files'
deg_dir = f"{figures}/degs"
os.makedirs(deg_dir, exist_ok=True)
sc.set_figure_params(dpi=200, format="png")
sc.settings.figdir = figures

if __name__ == "__main__":
    mudata = mu.read(f'{sc_file}/multi_all_cells_raw.h5mu')
    adata = mudata.mod['rna']
    adata.X = adata.layers['soupx'].copy()
    adata.X = adata.X.todense()
    sc.pp.normalize_total(adata, target_sum=1e6)
    hyp_ks_stats = pd.DataFrame(
        index=adata.obs["celltype"].cat.categories,
        columns=["n_genes_up_Hyperoxia", "n_genes_dn_Hyperoxia"],
        data=np.nan,
    )
    hyp_deg_dict = {}
    dataset_dict = {}
    # for lineage in adata.obs["lineage"].unique():
    #     print(lineage)
    #     lin_degs = f"{deg_dir}/{lineage}"
    #     os.makedirs(lin_degs, exist_ok=True)
    #     lin_adata = adata[adata.obs["lineage"] == lineage]
    #     lin_adata_norm = lin_adata[lin_adata.obs["treatment"] == 'Normoxia']
    #     with pd.ExcelWriter(
    #         f"{lin_degs}/{lineage.lower()}_cellsubtype_Hyperoxia_degs.xlsx",
    #         engine="xlsxwriter",
    #     ) as writer:
    #         for ct in sorted(lin_adata.obs["celltype"].unique()):
    #             ct_adata = lin_adata[lin_adata.obs["celltype"] == ct].copy()
    #             hyp = ct_adata[ct_adata.obs["treatment"] == "Hyperoxia"].copy()
    #             nor = ct_adata[ct_adata.obs["treatment"] == "Normoxia"].copy()
    #             if len(nor.obs_names) == 0 or len(hyp.obs_names) == 0:
    #                 print(ct)
    #                 continue
    #             comp = anndataks_ck(
    #                 nor, hyp, log1p=False, labels=(f"Normoxia", "Hyperoxia")
    #             )
    #             comp = comp.sort_values("statistic", ascending=False)
    #             comp.to_excel(writer, sheet_name=ct[:31])
    #             comp_stats = comp[comp["statistic"] > 0.3]
    #             comp_stats_up = comp_stats[comp_stats["log2_fold_change"] > 0]
    #             comp_stats_dn = comp_stats[comp_stats["log2_fold_change"] < 0]
    #             hyp_ks_stats.at[ct, "n_genes_up_Hyperoxia"] = len(comp_stats_up.index)
    #             hyp_ks_stats.at[ct, "n_genes_dn_Hyperoxia"] = len(comp_stats_dn.index)
    #             hyp_ks_stats.at[ct, "n_cells_nor"] = len(nor.obs_names)
    #             hyp_ks_stats.at[ct, "n_cells_hyp"] = len(hyp.obs_names)
    #             hyp_deg_dict[ct] = comp
    #     with pd.ExcelWriter(f'{lin_degs}/{lineage}_markers.xlsx', engine='xlsxwriter') as writer:
    #         for ct in sorted(lin_adata_norm.obs['celltype'].unique()):
    #             ct_adata = lin_adata_norm[lin_adata_norm.obs['celltype'] == ct].copy()
    #             nct_adata = lin_adata_norm[lin_adata_norm.obs['celltype'] != ct].copy()
    #             if len(lin_adata_norm.obs_names) == 0 or len(lin_adata_norm.obs_names) == 0:
    #                 continue
    #             comp = anndataks_ck(nct_adata,
    #                                 ct_adata,
    #                                 log1p=False,
    #                                 labels=(f'rest', f'{ct}'))
    #             comp = comp.sort_values('statistic', ascending=False)
    #             comp.to_excel(writer, sheet_name=ct[:31])
    #
    #
    #     for ct in lin_adata_norm.obs['celltype'].cat.categories:
    #         print(ct)
    #         ct_fp = f'{lin_degs}/cell_type_comparisons/{ct}'
    #         os.makedirs(ct_fp, exist_ok=True)
    #         with pd.ExcelWriter(f'{ct_fp}/{ct}_comparisons.xlsx', engine="xlsxwriter") as writer:
    #             ct_adata = lin_adata_norm[lin_adata_norm.obs['celltype'] == ct].copy()
    #             for ct2 in lin_adata_norm.obs['celltype'].cat.categories:
    #                 if ct == ct2:
    #                     continue
    #                 else:
    #                     ct2_adata = lin_adata_norm[lin_adata_norm.obs['celltype'] == ct2].copy()
    #                     comp = anndataks_ck(ct2_adata,
    #                                         ct_adata,
    #                                         labels=(f'{ct2}', f'{ct}')
    #                                         )
    #                     comp = comp.sort_values('statistic', ascending=False)
    #                     comp.to_excel(writer, sheet_name=ct2[:31])
    #             writer.close()

    with pd.ExcelWriter(
        f"{deg_dir}/celltype_hyperoxia_degs.xlsx", engine="xlsxwriter"
    ) as writer:
        for lineage in adata.obs['lineage'].cat.categories:
            deg_dir_lin_fn = f'{deg_dir}/{lineage}/{lineage}_celltype_hyperoxia_degs.xlsx'
            hyper_df_dt = pd.read_excel(deg_dir_lin_fn, sheet_name=None, index_col=0, header=0)
            for key in sorted(hyper_df_dt.keys()):
                hyper_df_dt[key].to_excel(writer, sheet_name=key[:31])
    hyp_ks_stats.to_csv(f"{deg_dir}/Hyperoxia_deg_counts.csv")

    hyper_deg_fn = f'{deg_dir}/celltype_hyperoxia_degs.xlsx'
    final_dict = {}
    for ct in adata.obs['celltype'].cat.categories:
        if ct == 'Imm/endo cells':
            ct = 'Imm_endo cells'
        try:
            hyper_df = pd.read_excel(hyper_deg_fn, sheet_name=ct, index_col=0, header=0)
        except:
            print(f'{ct} not compared')
            continue
        final_df = hyper_df.copy()
        final_df = final_df.loc[(final_df['statistic'] > 0.2)
                                # | (final_df['pvalue'] < 0.05)
        ]
        final_df = final_df.loc[(final_df['log2_fold_change'].abs() > 0.5)]
        final_df = final_df.sort_values('log2_fold_change', ascending=False)
        final_df = final_df.sort_values('statistic', ascending=False)
        final_dict[ct] = final_df
    with pd.ExcelWriter(f'{deg_dir}/cellsubtype_hyperoxia_degs_curated.xlsx', engine='xlsxwriter') as writer:
        for ct in final_dict.keys():
            final_dict[ct].to_excel(writer, sheet_name=ct[:31])
