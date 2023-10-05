import os
import scanpy as sc
import numpy as np
import pandas as pd
import string
import anndata
from gtfparse import read_gtf
from anndata import AnnData
gtf_fn = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf'
data = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/single_cell_files/soupx'
sc_file = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/single_cell_files'
qc = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/figures/qc'

os.makedirs(sc_file, exist_ok=True)
os.makedirs(qc, exist_ok=True)


if __name__ == '__main__':
    var = pd.read_csv(
        '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/cellranger_output/230609_aggregate/outs/filtered_feature_bc_matrix/features.tsv.gz',
        header=None,sep = '\t', index_col=None)
    var.columns = ['gene_id','gene_name_feature','measure','seqname','num1','num2']
    var.set_index('gene_name_feature',drop=False,inplace=True)
    var.rename_axis('feature',inplace=True)
    gtf = read_gtf(gtf_fn)
    print(gtf.head(1))
    adata = sc.read_10x_mtx(data)
    adata.var = var[['gene_id','gene_name_feature',]].loc[var['measure']=='Gene Expression']
    ## Add .var columns from gtf
    for column in ['gene_type','gene_name','seqname']:
        temp_dict = pd.Series(gtf[column].values, index=gtf['gene_id']).to_dict()
        adata.var[column] = [temp_dict[x] for x in adata.var['gene_id'].values]
    ## Add .obs columns custom
    mouse = []
    sex= []
    for x in adata.obs_names:
        if x.split('-')[1] == '1':
            mouse.append('nor-1')
            sex.append('F')
        elif x.split('-')[1] == '2':
            mouse.append('nor-2')
            sex.append('F')
        elif x.split('-')[1] == '3':
            mouse.append('nor-3')
            sex.append('M')
        elif x.split('-')[1] == '4':
            mouse.append('hyp-1')
            sex.append('F')
        elif x.split('-')[1] == '5':
            mouse.append('hyp-2')
            sex.append('M')
    adata.obs['mouse'] = mouse
    adata.obs['sex'] = sex
    adata.obs['treatment'] = ['Hyperoxia' if x.split('-')[0].startswith('h') else 'Normoxia' for x in adata.obs['mouse']]
    adata.obs['timepoint'] = 'P7'
    adata.write(f'{sc_file}/multiome_gex_all_cells_raw.gz.h5ad', compression='gzip')