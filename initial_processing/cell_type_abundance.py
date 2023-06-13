import numpy as np
import pandas as pd
import os
import scanpy as sc
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import ticker


figures = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/figures/cell_type_abundance'
data = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/single_cell_files'
os.makedirs(figures, exist_ok=True)
sc.set_figure_params(dpi=200, format="png")
sc.settings.figdir = figures

if __name__ == "__main__":
    adata = sc.read(f"{data}/multiome_gex_processed_cell_typed_raw.gz.h5ad")
    adata.X = adata.X.todense()
    df = adata.obs
    df2 = df.groupby('treatment')['celltype'].value_counts(normalize=True).mul(100).rename(
        '% cells sampled').reset_index()
    df2['celltype'] = df2['level_1']
    for lineage in adata.obs['lineage'].cat.categories:
        lin_adata = adata[adata.obs['lineage'] == lineage].copy()
        df3 = df2.loc[df2['celltype'].isin(lin_adata.obs['celltype'].cat.categories.tolist())].copy()
        df3['celltype'] = df3['celltype'].astype('string')
        order = lin_adata.obs['celltype'].cat.categories
        sns.catplot(data=df3,
                    x='celltype',
                    y='% cells sampled',
                    hue='treatment',
                    order=order,
                    hue_order=['Normoxia', 'Hyperoxia'],
                    kind='bar')
        plt.xticks(rotation=90)
        plt.yscale('log')
        ax = plt.gca()
        ax.yaxis.set_major_formatter(ticker.ScalarFormatter())  # set regular formatting
        plt.title(f'{lineage}')
        plt.savefig(f'{figures}/{lineage}_abunance_by_treatment.png', bbox_inches='tight')
    df2 = df.groupby('treatment')['lineage'].value_counts(normalize=True).mul(100).rename(
        '% cells sampled').reset_index()
    df2['lineage'] = df2['level_1']
    df3 = df2.loc[df2['lineage'].isin(adata.obs['lineage'].cat.categories.tolist())].copy()
    df3['lineage'] = df3['lineage'].astype('string')
    sns.catplot(data=df3,
                x='lineage',
                y='% cells sampled',
                hue='treatment',
                hue_order=['Normoxia', 'Hyperoxia'],
                kind='bar')
    plt.xticks(rotation=90)
    plt.savefig(f'{figures}/lineage_abunance_by_treatment.png', bbox_inches='tight')

