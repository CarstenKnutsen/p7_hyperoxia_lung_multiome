'''Write libraries files for cellranger-arc for multiome sequencing
author: Carsten Knutsen
Date: May 8th 2023
env: Run on sherlock
'''

import os
import csv
root_dir = '/scratch/users/cknutsen/2023_multiomic_sequencing/data/'
def write_csv(sample, path, atac_id, rna_id):

    # field names
    fields = ['fastqs', 'sample', 'library_type']

    # data rows of csv file
    rows = [[f'{path}{sample}/rna', f'{rna_id}', 'Gene Expression'],
            [f'{path}{sample}/atac', f'{atac_id}', 'Chromatin Accessibility'],
            ]

    # name of csv file
    filename = f"{path}{sample}/{sample}_libraries.csv"

    # writing to csv file
    with open(filename, 'w') as csvfile:
        # creating a csv writer object
        csvwriter = csv.writer(csvfile)

        # writing the fields
        csvwriter.writerow(fields)

        # writing the data rows
        csvwriter.writerows(rows)
    return
if __name__ =='__main__':
    for sample in os.listdir(root_dir):
        atac_id = '_'.join(os.listdir(f'{root_dir}{sample}/atac')[0].split(".")[0].split('_')[:-4])
        rna_id  = '_'.join(os.listdir(f'{root_dir}{sample}/rna')[0].split(".")[0].split('_')[:-4])
        write_csv(sample,root_dir,atac_id,rna_id)




