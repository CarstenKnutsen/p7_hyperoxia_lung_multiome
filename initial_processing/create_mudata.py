'''
Goal: Import output of 'cellranger-arc aggregate' and create anndata and mudata objects of the ATAC/RNA data

'''

import scanpy as sc
import muon as mu
import os

gtf_fn = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf'
rna_data = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/single_cell_files/soupx'
sc_file = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/single_cell_files'
h5_file = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/cellranger_output/230609_aggregate/outs/filtered_feature_bc_matrix.h5'
tf_matrix_dir = "/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/cellranger_output/230609_aggregate/outs/analysis/tf_analysis/filtered_tf_bc_matrix"

os.makedirs(sc_file, exist_ok=True)


if __name__ == '__main__':
    adata = sc.read_10x_mtx(rna_data)
    mudata = mu.read_10x_h5("/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/cellranger_output/230609_aggregate/outs/filtered_feature_bc_matrix.h5")
    mudata['rna'].layers['raw'] = mudata['rna'].X.copy()
    mudata['rna'].X = adata.X.copy()
    mudata['rna'].var_names = adata.var_names
    mudata['rna'].obs_names = adata.obs_names
    ## Add .obs columns custom
    mouse = []
    sex= []
    for x in mudata.obs_names:
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
    mudata.obs['mouse'] = mouse
    mudata.obs['sex'] = sex
    mudata.obs['treatment'] = ['Hyperoxia' if x.split('-')[0].startswith('h') else 'Normoxia' for x in mudata.obs['mouse']]
    mudata.obs['timepoint'] = 'P7'
    mudata.mod['atac'].uns['atac'] = dict(mudata.mod['atac'].uns['atac'])  # anndata does not write ordered dict
    mudata['rna'].var_names_make_unique()
    for modality in mudata.mod.keys():
        mudata.mod[modality].obs = mudata.obs
        mudata[modality].write(f'{sc_file}/{modality}_all_cells_raw.gz.h5ad', compression='gzip')
    print(mudata)
    mudata.write(f'{sc_file}/multi_all_cells_raw.h5mu')
