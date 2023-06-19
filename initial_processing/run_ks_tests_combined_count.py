import numpy as np
import pandas as pd
import os
import scanpy as sc
import sys

sys.path.insert(1, "/home/carsten/alvira_bioinformatics/lungsc_ck")
from ck_functions import anndataks_ck

figures = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/figures/tissue_embedding'
data = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/single_cell_files'
deg_dir = f"{figures}/degs"
os.makedirs(deg_dir, exist_ok=True)
sc.set_figure_params(dpi=200, format="png")
sc.settings.figdir = figures

if __name__ == "__main__":
    adata = sc.read(f"{data}/multiome_gex_processed_cell_typed_raw.gz.h5ad")
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
    for lineage in ['immune','mesenchymal']:
        print(lineage)
        lin_degs = f"{deg_dir}/{lineage}"
        os.makedirs(lin_degs, exist_ok=True)
        lin_adata = adata[adata.obs["lineage"] == lineage]
        with pd.ExcelWriter(
            f"{lin_degs}/{lineage.lower()}_cellsubtype_Hyperoxia_degs.xlsx",
            engine="xlsxwriter",
        ) as writer:
            for ct in sorted(lin_adata.obs["celltype"].unique()):
                ct_adata = lin_adata[lin_adata.obs["celltype"] == ct].copy()
                hyp = ct_adata[ct_adata.obs["treatment"] == "Hyperoxia"].copy()
                nor = ct_adata[ct_adata.obs["treatment"] == "Normoxia"].copy()
                if len(nor.obs_names) == 0 or len(hyp.obs_names) == 0:
                    print(ct)
                    continue
                comp = anndataks_ck(
                    nor, hyp, log1p=False, labels=(f"Normoxia", "Hyperoxia")
                )
                comp = comp.sort_values("statistic", ascending=False)
                comp["gene_type"] = adata.var["gene_type"]
                comp.to_excel(writer, sheet_name=ct[:31])
                comp_stats = comp[comp["statistic"] > 0.3]
                comp_stats_up = comp_stats[comp_stats["log2_fold_change"] > 0]
                comp_stats_dn = comp_stats[comp_stats["log2_fold_change"] < 0]
                hyp_ks_stats.at[ct, "n_genes_up_Hyperoxia"] = len(comp_stats_up.index)
                hyp_ks_stats.at[ct, "n_genes_dn_Hyperoxia"] = len(comp_stats_dn.index)
                hyp_ks_stats.at[ct, "n_cells_nor"] = len(nor.obs_names)
                hyp_ks_stats.at[ct, "n_cells_hyp"] = len(hyp.obs_names)
                hyp_deg_dict[ct] = comp
        with pd.ExcelWriter(f'{lin_degs}/{lineage}_markers.xlsx', engine='xlsxwriter') as writer:
            for ct in sorted(lin_adata.obs['celltype'].unique()):
                ct_adata = lin_adata[lin_adata.obs['celltype'] == ct].copy()
                nct_adata = lin_adata[lin_adata.obs['celltype'] != ct].copy()
                if len(ct_adata.obs_names) == 0 or len(nct_adata.obs_names) == 0:
                    continue
                comp = anndataks_ck(nct_adata,
                                    ct_adata,
                                    log1p=False,
                                    labels=(f'rest', f'{ct}'))
                comp = comp.sort_values('statistic', ascending=False)
                comp['gene_type'] = adata.var['gene_type']
                comp.to_excel(writer, sheet_name=ct[:31])


        # for ct in lin_adata.obs['celltype'].cat.categories:
        #     print(ct)
        #     ct_fp = f'{lin_degs}/cell_type_comparisons/{ct}'
        #     os.makedirs(ct_fp, exist_ok=True)
        #     with pd.ExcelWriter(f'{ct_fp}/{ct}_comparisons.xlsx', engine="xlsxwriter") as writer:
        #         ct_adata = lin_adata[lin_adata.obs['celltype'] == ct].copy()
        #         for ct2 in lin_adata.obs['celltype'].cat.categories:
        #             if ct == ct2:
        #                 continue
        #             else:
        #                 ct2_adata = lin_adata[lin_adata.obs['celltype'] == ct2].copy()
        #                 comp = anndataks_ck(ct2_adata,
        #                                     ct_adata,
        #                                     labels=(f'{ct2}', f'{ct}')
        #                                     )
        #                 comp = comp.sort_values('statistic', ascending=False)
        #                 comp['gene_type'] = adata.var['gene_type']
        #                 comp.to_excel(writer, sheet_name=ct2[:31])
        #         writer.close()

    with pd.ExcelWriter(
        f"{deg_dir}/cellsubtype_Hyperoxia_degs.xlsx", engine="xlsxwriter"
    ) as writer:
        for key in sorted(hyp_deg_dict.keys()):
            hyp_deg_dict[key].to_excel(writer, sheet_name=key[:31])
    hyp_ks_stats.to_csv(f"{deg_dir}/Hyperoxia_deg_counts.csv")
