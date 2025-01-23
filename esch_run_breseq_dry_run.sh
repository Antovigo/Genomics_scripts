#!/bin/bash
#SBATCH -c 1                               # Request 8 cores
#SBATCH -t 0-00:30                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem-per-cpu=50M                   # Memory total in GiB
#SBATCH -o output/%x_%j.out                 # File to which STDOUT will be written
#SBATCH -e output/%x_%j.err                 # File to which STDERR will be writtenh

# Usage
# To process the folders whose names contain XYZ
# sbatch -J XYZ ~/Code/Genomics_scripts/esch_run_breseq.sh XYZ

# References
ecoli_chr="/n/scratch/users/a/anv888/Genomics/Reference_genomes/Escherichia_coli_MG1655_U00096.3.gb"
phi80_chr="/n/scratch/users/a/anv888/Genomics/Reference_genomes/phage_phi80_JX871397.gb"

vnat_chr1="/n/scratch/users/a/anv888/Genomics/Reference_genomes/NZ_CP009977.1.gb"
vnat_chr2="/n/scratch/users/a/anv888/Genomics/Reference_genomes/NZ_CP009978.1.gb"

# Folders
reads_folder="/n/scratch/users/a/anv888/Genomics/Illumina_DNA_Reads_LTEE_ecoli"
breseq_folder="/n/scratch/users/a/anv888/Genomics/Breseq_output_LTEE_ecoli"

lineage=$1

#module load miniconda3
#source activate breseq

# Clear empty folders
#find $breseq_folder -type d -empty -delete

echo "Files to process"
for i in $reads_folder/$lineage
	echo $i
