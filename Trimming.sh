# Install Trimmomatic if not already installed
sudo apt-get install trimmomatic

# Create a directory for trimmed reads
mkdir trimmed_reads

# Run Trimmomatic for paired-end reads
java -jar trimmomatic-0.39.jar PE -phred33 \
    input_R1.fastq input_R2.fastq \
    trimmed_reads/output_R1_paired.fastq trimmed_reads/output_R1_unpaired.fastq \
    trimmed_reads/output_R2_paired.fastq trimmed_reads/output_R2_unpaired.fastq \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
