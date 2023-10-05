#!/bin/bash
echo 'starting TF!!!!!!!!!!!!!!!'
python differential_analysis/atac_da_tests_tf.py
echo 'starting Peaks!!!!!!!!!!!!!'
python differential_analysis/atac_da_tests.py
echo 'starting GENES!!!!!!!!!!!!!'
python differential_analysis/rna_deg_tests.py
echo 'starting GENES KS!!!!!!!!!!!!!'
python differential_analysis/rna_ks_tests.py
echo 'starting Pathway!!!!!!!!!!!!!'
python differential_analysis/rna_pathway_analysis_msbio.py
