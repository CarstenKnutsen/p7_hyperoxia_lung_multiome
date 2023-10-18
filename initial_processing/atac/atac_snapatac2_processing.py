
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
gtf = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf'
chrom_sizes_fn = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/refdata-cellranger-arc-mm10-2020-A-2.0.0/sizes.genome'
genome = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa'
output_f = f'{sc_file}/multiome_snapatac_pilot.h5ad'
obs_list = ['lineage','celltype','treatment','mouse','sex']
os.makedirs(figures, exist_ok=True)
sc.set_figure_params(dpi=300, format="png")
sc.settings.figdir = figures

if __name__ == "__main__":
    adata_rna = sc.read(f'{sc_file}/p7_multiome_rna_processed.gz.h5ad')
    chrom_sizes = pd.read_csv(chrom_sizes_fn, sep='\t', header=None, index_col=0).to_dict()
    chrom_sizes = chrom_sizes[1]
    data = snap.pp.import_data(
        fragment_file,
        chrom_sizes,
        file=output_f,  # Optional
        sorted_by_barcode=False,
    )
    fig = snap.pl.frag_size_distr(data, show=False)
    fig.update_yaxes(type="log")
    fig.write_image(f'{output_fol}/fragment_dist.png')
    snap.metrics.tsse(data, gtf)
    fig = snap.pl.tsse(data, interactive=False)
    fig.write_image(f'{output_fol}/tsse.png')
    snap.pp.filter_cells(data, min_counts=1000, min_tsse=0, max_counts=100000)
    overlap_names = list(set(adata_rna.obs_names) & set(data.obs_names))
    data.subset(obs_indices=overlap_names)
    for obs in obs_list:
        data.obs[obs] = adata_rna.obs[obs].loc[overlap_names]
    snap.pp.add_tile_matrix(data)
    nap.pp.select_features(data, n_features=500000)
    snap.tl.spectral(data)
    snap.pp.knn(data)
    snap.tl.leiden(data)
    snap.tl.umap(data)
    pd.DataFrame(data=data.obsm['X_umap'],
                 index=data.obs_names).to_csv(f'{output_fol}/tile_matrix_umap_coords.csv')
    for obs in obs_list:
        fig = snap.pl.umap(data, color='leiden', interactive=False, height=500)
        fig.write_image(f'{output_fol}/umap_{obs}.png')
    snap.tl.macs3(data, groupby='celltype')
    peaks = snap.tl.merge_peaks(data.uns['macs3'], chrom_sizes)
    peaks.to_pandas().to_csv(f'{figures}/peaks_df.csv')
    peak_matrix = snap.pp.make_peak_matrix(data, use_rep=peaks['Peaks'])
    peak_matrix.write(f"{sc_file}/snapatac2_peak_matrix.gz.h5ad", compression='gzip')
    gene_matrix = snap.pp.make_gene_matrix(data, gtf)
    data.close()
    sc.pp.filter_genes(gene_matrix, min_cells=5)
    sc.pp.normalize_total(gene_matrix, key_added=None, target_sum=1e6)
    sc.pp.log1p(gene_matrix, base=10)
    sce.pp.magic(gene_matrix, solver="approximate")
    gene_matrix.write(f"{sc_file}/snapatac2_gene_matrix.gz.h5ad", compression='gzip')


















