# Index the reference sequence
bwa index reference.fasta

# Pair-end
bwa mem reference.fasta reads_1.fastq reads_2.fastq > aligned_reads.sam

# Convert to binary and sort
samtools view -bS aligned_reads.sam | samtools sort - > aligned_reads.sorted.bam

# Index
samtools index aligned_reads.sorted.bam

# Extract the matches
samtools fastq -F 4 aligned_reads.sorted.bam > mapped_reads.fastq

# Extract the partial matches
samtools view -F 4 aligned_reads.sorted.bam | cut -f1-11 | grep "S" > partially_mapped.sam 

