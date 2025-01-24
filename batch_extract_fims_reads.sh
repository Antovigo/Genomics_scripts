#!/bin/bash
#SBATCH -c 4                               # Request 8 cores
#SBATCH -t 1-00:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem-per-cpu=8G                   # Memory total in GiB
#SBATCH -o output/%x_%j.out                 # File to which STDOUT will be written
#SBATCH -e output/%x_%j.err                 # File to which STDERR will be writtenh

# Folders
reads_folder="/n/scratch/users/a/anv888/Genomics/Illumina_DNA_Reads_LTEE_ecoli"
output_folder="/n/scratch/users/a/anv888/Genomics/FimS_Reads"

module load miniconda3
source activate breseq

/home/anv888/Code/Genomics_scripts/extract_fims.py $reads_folder $output_folder
