#!/usr/bin/bash
#SBATCH --job-name=aggregate
#SBATCH --output=aggregate.%j.out
#SBATCH --error=aggregate.%j.err
#SBATCH --time=2-00:00:00
#SBATCH -p normal
#SBATCH -n 64
#SBATCH --mem=128GB	

cellranger-arc aggr --id=aggregate_w_secondary_no_nor2 \
                  --csv=/scratch/users/cknutsen/2023_multiomic_sequencing/analysis/aggregate.csv \
                  --normalize=none \
                  --reference=/home/groups/calvira/software/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
                  

