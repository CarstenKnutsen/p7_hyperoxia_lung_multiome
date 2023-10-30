#!/bin/bash
annotatePeaks.pl /home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/pilot/231016_snapatac2_pilot/all_peaks.bed \
 /home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa \
  -gtf /home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf > /home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/pilot/231016_snapatac2_pilot/peak_homer_annotation.txt
