# Install FastQC if not already installed
sudo apt-get install fastqc

# Create a directory for FastQC results
mkdir fastqc_results

# Run FastQC on all FASTQ files
fastqc *.fastq -o fastqc_results/
