#!/bin/bash
#SBATCH -c 8                               # Request 8 cores
#SBATCH -t 2-12:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem-per-cpu=1G                   # Memory total in GiB
#SBATCH -o output/%x_%j.out                 # File to which STDOUT will be written
#SBATCH -e output/%x_%j.err                 # File to which STDERR will be writtenh

# References
#ecoli_chr="/n/scratch/users/a/anv888/Genomics/Reference_genomes/rpoS_reference.gb"
ecoli_chr="/n/scratch/users/a/anv888/Genomics/MG1655_U00096_with_phi80_genome.gb"
#vnat_chr1="/n/scratch/users/a/anv888/Genomics/Reference_genomes/NZ_CP009977.1.gb"
#vnat_chr2="/n/scratch/users/a/anv888/Genomics/Reference_genomes/NZ_CP009978.1.gb"

# Folders
reads_folder="/n/scratch/users/a/anv888/Genomics/Illumina_DNA_Reads_MFepsilon"
breseq_folder="/n/scratch/users/a/anv888/Genomics/Breseq_output_MFepsilon"

lineage=$1

module load miniconda3
source activate breseq

# Clear empty folders
#find $breseq_folder -type d -empty -delete

echo "Files to process:"
for i in $reads_folder/$lineage
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
        # E. coli version, fine, no polymorphisms (no -p option)
        breseq -r "$ecoli_chr" -j 8 "$input_folder"/*.fastq -o "$output_folder"

        end_time=$(date)
        echo "Finished the processing of sample $i ($end_time)."
        echo ""

        # Prevent R from accumulating disk space
        #killall -9 R
        #killall python
        #killall python3

done

echo 'Everything is done!'
