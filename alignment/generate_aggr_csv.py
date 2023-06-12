'''Write libraries files for cellranger-arc for multiome sequencing
author: Carsten Knutsen
Date: May 8th 2023
env: Run on sherlock
'''

import os
import csv
root_dir = '/scratch/users/cknutsen/2023_multiomic_sequencing/data/'
filename = '/scratch/users/cknutsen/2023_multiomic_sequencing/analysis/aggregate.csv'
if __name__ =='__main__':
    fields = ['library_id', 'atac_fragments', 'per_barcode_metrics','gex_molecule_info']
    rows = []
    for sample in os.listdir(root_dir):
        row = [f'{sample}',f'{root_dir}{sample}/outs/atac_fragments.tsv.gz',f'{root_dir}{sample}/outs/per_barcode_metrics.csv',f'{root_dir}{sample}/outs/gex_molecule_info.h5']
        rows.append(row)
    with open(filename, 'w') as csvfile:
        # creating a csv writer object
        csvwriter = csv.writer(csvfile)

        # writing the fields
        csvwriter.writerow(fields)

        # writing the data rows
        csvwriter.writerows(rows)




