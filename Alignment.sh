# Install HISAT2 if not already installed
sudo apt-get install hisat2

# Build HISAT2 index (if not already built)
hisat2-build genome.fa genome_index

# Create a directory for alignment results
mkdir hisat2_results

# Run HISAT2 for alignment
hisat2 -x genome_index -1 trimmed_reads/output_R1_paired.fastq -2 trimmed_reads/output_R2_paired.fastq -S hisat2_results/alignment.sam
