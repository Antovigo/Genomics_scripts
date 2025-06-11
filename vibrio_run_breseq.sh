#!/bin/bash
#SBATCH -c 8                               # Request 8 cores
#SBATCH -t 1-12:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem-per-cpu=1G                   # Memory total in GiB
#SBATCH -o output/%x_%j.out                 # File to which STDOUT will be written
#SBATCH -e output/%x_%j.err                 # File to which STDERR will be writtenh

# Usage
# To process the folders whose names contain XYZ
# sbatch -J XYZ ~/Code/Genomics_scripts/esch_run_breseq.sh XYZ

# References
#ecoli_chr="/n/scratch/users/a/anv888/Genomics/Reference_genomes/Escherichia_coli_MG1655_U00096.3.gb"
#phi80_chr="/n/scratch/users/a/anv888/Genomics/Reference_genomes/phage_phi80_JX871397.gb"

vnat_chr1="/n/scratch/users/a/anv888/Genomics/Reference_genomes/NZ_CP009977.1.gb"
vnat_chr2="/n/scratch/users/a/anv888/Genomics/Reference_genomes/NZ_CP009978.1.gb"

# Folders
reads_folder="/n/scratch/users/a/anv888/Genomics/Illumina_DNA_Reads_06_25"
breseq_folder="/n/scratch/users/a/anv888/Genomics/Breseq_output_06_25"

lineage=$1

module load miniconda3
source activate breseq

# Clear empty folders
#find $breseq_folder -type d -empty -delete

echo "Files to process (wild match):"
for i in $reads_folder/$lineage*
do
	echo $i
done

echo "Processing files:"
for i in $reads_folder/$lineage*
do
        name=$(basename $i)
        input_folder=$reads_folder/$name
        output_folder=$breseq_folder/$name

        # Skip if the output folder already exists
        if [ -d "$output_folder" ]; then
                echo "Output folder $output_folder already exists, skipping sample $name"
                continue
        fi

        start_time=$(date)
        echo ""
        echo "Starting the processing of sample $i ($start_time)."
        mkdir -p $output_folder

        # Run breseq
	breseq -r "$vnat_chr1" -r "$vnat_chr2" -j 8 -p "$input_folder"/*.fastq -o "$output_folder" --junction-minimum-pos-hash-score 2 --junction-minimum-side-match 5 --polymorphism-frequency-cutoff 0.01 --polymorphism-minimum-variant-coverage-each-strand 1 --polymorphism-minimum-variant-coverage 1 --polymorphism-reject-indel-homopolymer-length 8 --polymorphism-reject-surrounding-homopolymer-length 8  --polymorphism-score-cutoff 1.5
	#breseq -r "$ecoli_chr" -r "$phi80_chr" -j 8 -p "$input_folder"/*.fastq -o "$output_folder" --junction-minimum-pos-hash-score 2 --junction-minimum-side-match 5 --polymorphism-frequency-cutoff 0.01 --polymorphism-minimum-variant-coverage-each-strand 1 --polymorphism-minimum-variant-coverage 1 --polymorphism-reject-indel-homopolymer-length 8 --polymorphism-reject-surrounding-homopolymer-length 8  --polymorphism-score-cutoff 1.5

        end_time=$(date)
        echo "Finished the processing of sample $i ($end_time)."
        echo ""

        # Prevent R from accumulating disk space
        #killall -9 R
        #killall python
        #killall python3

done

echo 'Everything is done!'
