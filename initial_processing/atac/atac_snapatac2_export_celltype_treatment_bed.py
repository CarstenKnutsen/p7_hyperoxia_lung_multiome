
'''
Goal: Run SnapATAC2 pipeline on multiome snATAC data and create tile, and peak matrices
Author:Carsten Knutsen
Date:231017
conda_env:snapatac
'''

import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd
import snapatac2 as snap
import os

figures = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/figures/atac/snapatac2'
sc_file = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/single_cell_files'
fragment_file = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/cellranger_output/230609_aggregate/outs/atac_fragments.tsv.gz'
output_f = f'{sc_file}/multiome_snapatac_pilot.h5ad'
os.makedirs(figures, exist_ok=True)
sc.set_figure_params(dpi=300, format="png")
sc.settings.figdir = figures

if __name__ == "__main__":
    frag_mat = snap.read(f'{sc_file}/multiome_snapatac_pilot.h5ad')
    frag_mat.obs['treatment_celltype'] = frag_mat.obs['treatment'] + '_' + frag_mat.obs['celltype']

















