#!/usr/bin/bash
#SBATCH --job-name=google_drive_backup
#SBATCH --output=google_drive_backup.%j.out
#SBATCH --error=google_drive_backup.%j.err
#SBATCH --time=2-00:00:00
#SBATCH -p normal
#SBATCH -n 64
#SBATCH --mem=128GB	

ml load system rclone

rclone copy /scratch/users/cknutsen/2023_multiomic_sequencing/analysis/aggregate_w_secondary_no_nor2 cknutsen_google_drive:multiome/230609_aggregate
