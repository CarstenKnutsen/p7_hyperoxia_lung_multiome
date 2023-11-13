#!/bin/bash
/home/carsten/homer/bin/annotatePeaks.pl /home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/figures/atac/snapatac2_no_nor3/all_peaks.bed \
 /home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa \
  -gtf /home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf > /home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/figures/atac/snapatac2_no_nor3/peak_homer_annotation.txt
