# Install SAMtools if not already installed
sudo apt-get install samtools

# Create a directory for BAM results
mkdir bam_results

# Convert SAM to sorted BAM
samtools view -Sb hisat2_results/alignment.sam | samtools sort -o bam_results/alignment_sorted.bam
