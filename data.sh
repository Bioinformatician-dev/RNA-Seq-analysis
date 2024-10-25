# Check the quality of the raw reads
fastqc sample_R1.fastq sample_R2.fastq -o fastqc_reports/

# Trim the reads to remove low-quality sequences and adapters
trimmomatic PE -phred33 sample_R1.fastq sample_R2.fastq \
trimmed_sample_R1_paired.fastq trimmed_sample_R1_unpaired.fastq \
trimmed_sample_R2_paired.fastq trimmed_sample_R2_unpaired.fastq \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36

# Build the STAR index
STAR --runMode genomeGenerate --genomeDir ./genome_index \
--genomeFastaFiles genome.fa --sjdbGTFfile annotations.gtf \
--sjdbOverhang 99

# Build the STAR index
STAR --runMode genomeGenerate --genomeDir ./genome_index \
--genomeFastaFiles genome.fa --sjdbGTFfile annotations.gtf \
--sjdbOverhang 99

# Align the trimmed reads
STAR --runThreadN 4 --genomeDir ./genome_index \
--readFilesIn trimmed_sample_R1_paired.fastq trimmed_sample_R2_paired.fastq \
--outFileNamePrefix sample_ --outSAMtype BAM SortedByCoordinate

# Count reads that overlap with features in the GTF
featureCounts -a annotations.gtf -o counts.txt sample_Aligned.sortedByCoord.out.bam
