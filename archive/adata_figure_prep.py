'''
Goal: Set categories and obs metadata orders on multiomic object
Date: July 11th 2023
Author: Carsten Knutsen
Conda env: multiome
'''

import os
import pandas as pd
import scanpy as sc
import numpy as np
import seaborn as sns
import itertools
import muon as mu


figures = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/figures/figures_for_paper/tables'
os.makedirs(figures,exist_ok=True)
sc_file = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/single_cell_files'
abv_dict = {'Arterial EC':'Art',
            'Cap1':'Ca1',
            'Cap2':'Ca2',
            'Intermediate cap':'ICa',
            'Lymphatic EC':'Lym',
            'Proliferating EC':'PEC',
            'Venous EC':'Ven',
            'AT1':'AT1',
            'AT1_AT2':'A12',
            'AT2':'AT2',
            'Basal':'Bas',
            'Ciliated':'Cil',
            'Club':'Clu',
            'Goblet':'Gob',
            'Neuroendocrine':'Neu',
            'Proliferating AT2':'PA2',
            'Alveolar macrophage':'AlM',
            'B cell':'BCe',
            'DC':'DC',
            'Imm_endo cells':'ImE',
            'Intermediate monocyte':'IMo',
            'Interstitial macrophage':'IMa',
            'Monocyte':'Mon',
            'Proliferating macrophage':'PMo',
            'T cell':'Tce',
            'Abberant muscle cell':'AMC',
            'Acta1+ mese cell':'HAc',
            'Adventitial fibroblast':'AdF',
            'Airway smooth muscle':'ASM',
            'Alveolar fibroblast':'AlF',
            'EpiMT':'EMT',
            'Mesothelial':'Mes',
            'Myofibroblast':'Myo',
            'Pericyte':'Per',
            'Proliferating fibroblast':'PFi',
            'Proliferating myofibroblast':'PMy',
            'Proliferating pericyte':'PPe',
            'Vascular smooth muscle':'VSM'
            }
if __name__ == '__main__':
    mdata = mu.read(f'{sc_file}/multi_all_cells_processed.h5mu')
    for mod in ['rna','atac']:
        adata = mdata.mod[mod]
        ct_order = []
        abv_order = []
        ct_number = 0
        ct_number_dict = {}
        adata.obs['celltype_abv'] = adata.obs['celltype']
        adata.obs['celltype_abv'].replace(abv_dict, inplace=True)
        for lin in adata.obs['lineage'].cat.categories:
            for ct in adata[adata.obs['lineage'] == lin].obs['celltype'].cat.categories:
                ct_order.append(ct)
                ct_number += 1
                ct_number_dict[ct] = str(ct_number)
            for ct in adata[adata.obs['lineage'] == lin].obs['celltype'].cat.categories:
                abv = adata[adata.obs['celltype']==ct].obs['celltype_abv'].values[0]
                abv_order.append(abv)
        adata.obs['ct_number_name'] = [f'{ct_number_dict[x]}. {x} ' for x in adata.obs['celltype']]
        adata.obs['ct_number'] = [f'{ct_number_dict[x]}' for x in adata.obs['celltype']]
        adata.obs['celltype'] = pd.Categorical(adata.obs['celltype'], categories=ct_order)
        adata.obs['celltype_abv'] = pd.Categorical(adata.obs['celltype_abv'], categories=abv_order)
        adata.obs['treatment'].cat.reorder_categories(['Normoxia', 'Hyperoxia'], inplace=True)
        adata.uns['treatment_colors'] = np.array(
            (sns.color_palette("vlag").as_hex()[0], sns.color_palette("vlag").as_hex()[-1]))
        adata.uns['celltype_colors'] = np.array(sc.pl.palettes.godsnot_102, dtype='object')[
                                           :len(adata.obs['celltype'].cat.categories)]
        adata.uns['ct_number_colors'] = adata.uns['celltype_colors'].copy()
        adata.uns['ct_number_name_colors'] = adata.uns['celltype_colors'].copy()
        adata.obs['capital_lineage'] = adata.obs['lineage'].str.capitalize()
        adata.uns['lineage_colors'] = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
        adata.uns['capital_lineage_colors'] = adata.uns['lineage_colors'].copy()
    rna = mdata.mod['rna']
    with pd.ExcelWriter(f'{figures}/metadata_counts.xlsx', engine='xlsxwriter') as writer:
        obs_list = ['lineage', 'treatment', 'sex', 'mouse','celltype']
        num_obs = len(obs_list) + 1
        for ind in range(0, num_obs):
            for subset in itertools.combinations(obs_list, ind):
                if len(subset) != 0:
                    subset = list(subset)
                    if len(subset) == 1:
                        key = subset[0]
                        rna.obs[key].value_counts().to_excel(writer, sheet_name=key)
                    else:
                        key = "_".join(subset)
                        rna.obs.groupby(subset[:-1])[subset[-1]].value_counts().to_excel(writer, sheet_name=key[:31])

    abv = pd.Series(pd.Series(rna.obs['celltype_abv'].values, index=rna.obs['celltype'].values).to_dict(),
                    name='abbreviation')
    num = pd.Series(pd.Series(rna.obs['ct_number'].values, index=rna.obs['celltype'].values).to_dict(),
                    name='number')
    lin = pd.Series(pd.Series(rna.obs['lineage'].values, index=rna.obs['celltype'].values).to_dict(),
                    name='lineage')
    df = pd.DataFrame([abv, num,lin]).T
    df['number'] = df['number'].astype('int')
    df = df.sort_values('number', ascending=True).to_csv(f'{figures}/names_table.csv')
    print(adata)
    mdata.write(f'{sc_file}/multi_all_cells_processed.h5mu')