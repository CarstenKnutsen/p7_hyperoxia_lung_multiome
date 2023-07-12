#!/usr/bin/bash
#SBATCH --job-name=alignment
#SBATCH --output=alignment.%j.out
#SBATCH --error=alignment.%j.err
#SBATCH --time=2-00:00:00
#SBATCH -p normal
#SBATCH -n 64
#SBATCH --mem=128GB	

for dir in /scratch/users/cknutsen/2023_multiomic_sequencing/data/*/
do
    dir2=${dir%*/}   
    dir2="${dir2##*/}"
    echo $dir2   
    cellranger-arc count --id="${dir2}" \
                       --reference=/home/groups/calvira/software/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
                       --libraries=$dir"$dir2"_libraries.csv   
done
cellranger-arc aggr --id=aggregate \
                  --csv=/scratch/users/cknutsen/2023_multiomic_sequencing/analysis/aggregate.csv \
                  --normalize=none \
                  --reference=/home/groups/calvira/software/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
                  --nosecondary

