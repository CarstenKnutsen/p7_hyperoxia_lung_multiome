#!/bin/bash
echo 'starting'
python initial_processing/create_mudata.py
echo 'starting RNA'
python initial_processing/rna/rna_qc_and_normalization.py
python initial_processing/rna/rna_all_lineage_cluster_embed.py
python initial_processing/rna/rna_cell_typing.py
echo 'starting ATAC'
python initial_processing/atac/atac_qc_clustering_embedding.py
echo 'Exporting'
python initial_processing/export_mudata_objects.py
echo 'starting TF!!!!!!!!!!!!!!!'
python differential_analysis/atac_da_tests_tf.py
echo 'starting Peaks!!!!!!!!!!!!!'
python differential_analysis/atac_da_tests.py
echo 'starting GENES!!!!!!!!!!!!!'
python differential_analysis/rna_deg_tests.py
echo 'starting GENES KS!!!!!!!!!!!!!'
#python differential_analysis/rna_ks_tests.py
echo 'starting Pathway!!!!!!!!!!!!!'
#python differential_analysis/rna_pathway_analysis_msbio.py
