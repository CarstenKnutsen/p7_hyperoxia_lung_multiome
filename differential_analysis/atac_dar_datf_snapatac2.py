"""
Goal: Run differential accesibility between peaks in peaks file, run motif enrichment on differentially accesible peaks
Author:Carsten Knutsen
Date:231017
conda_env:snapatac
"""

import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd
import snapatac2 as snap
import os

da_output = "/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/figures/atac/snapatac2_no_nor3/dap_datf_tile"
sc_file = "/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/single_cell_files"
genome = "/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa"
peaks_fn = "/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/figures/atac/snapatac2_no_nor3/peaks_df.csv"
peak_md_fn = "/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/figures/atac/snapatac2_no_nor3/peak_homer_annotation.txt"
peak_md_cols = [
    "Chr",
    "Start",
    "End",
    "Annotation",
    "Distance to TSS",
    "Nearest PromoterID",
    "Gene Name",
    "Gene Type",
]

os.makedirs(da_output, exist_ok=True)


if __name__ == "__main__":
    adata_peak = sc.read(f"{sc_file}/snapatac2_tile_matrix_no_nor3.h5ad")
    print(adata_peak)
    peak_md = pd.read_csv(peak_md_fn, sep="\t", index_col=0)
    peaks = pd.read_csv(peaks_fn, header=0, index_col=0)
    for lineage in sorted(adata_peak.obs["lineage"].unique()):
        da_lin_output = f"{da_output}/{lineage}"
        os.makedirs(da_lin_output, exist_ok=True)
        lin_adata = adata_peak[adata_peak.obs["lineage"] == lineage]
        # for treatment in lin_adata.obs['treatment'].unique():
        #     lin_treat = lin_adata[lin_adata.obs['treatment']==treatment]
        ### Hyperoxia
        da_lin_hyp_output = f"{da_lin_output}/hyperoxia"
        os.makedirs(da_lin_hyp_output, exist_ok=True)
        for ct in sorted(lin_adata.obs["celltype"].unique()):
            print(ct)
            ct_adata = lin_adata[lin_adata.obs["celltype"] == ct]
            ct_norm = ct_adata[ct_adata.obs["treatment"] == "Normoxia"]
            ct_hyper = ct_adata[ct_adata.obs["treatment"] == "Hyperoxia"]
            ct_peaks = peaks[(peaks[f'Normoxia_{ct}'] == True)|
            (peaks[f'Hyperoxia_{ct}'] == True)]["Peaks"].values
            if len(ct_norm.obs_names) < 10 or len(ct_hyper.obs_names) < 10:
                print(ct)
                print("Too few cells")
                continue
            diff_df = snap.tl.diff_test(
                adata_peak, ct_hyper.obs_names, ct_norm.obs_names,
            )
            diff_df_pd = diff_df.to_pandas()
            diff_df_pd.index = diff_df_pd["feature name"]
            diff_df_pd = diff_df_pd[diff_df_pd.columns.tolist()[1:]]
            # diff_df_pd[peak_md_cols] = peak_md.loc[diff_df_pd.index][peak_md_cols]
            diff_df_pd.to_csv(f"{da_lin_hyp_output}/{ct}_hyperoxia_dap.csv")
            diff_df_pd_t = diff_df_pd.loc[diff_df_pd["adjusted p-value"] < 0.05]
            diff_df_pd_t = diff_df_pd_t.loc[
                ~diff_df_pd_t.index.str.startswith("chrY")
            ]  # Drop any male genes
            df_pd_up = diff_df_pd_t[diff_df_pd_t["log2(fold_change)"] > 0].index.values
            df_pd_dn = diff_df_pd_t[diff_df_pd_t["log2(fold_change)"] < 0].index.values
            dict = {
                "Hyperoxia": df_pd_up,
                "Normoxia": df_pd_dn,
            }

            motifs = snap.tl.motif_enrichment(
                motifs=snap.datasets.Meuleman_2020(),
                regions=dict,
                background=ct_peaks,
                genome_fasta=genome,
                method="hypergeometric",
            )
            with pd.ExcelWriter(
                f"{da_lin_hyp_output}/{ct}_hyperoxia_datf.xlsx", engine="xlsxwriter"
            ) as writer:
                for k in motifs.keys():
                    motif_df = motifs[k].to_pandas().to_excel(writer, sheet_name=k)
