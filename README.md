# RNA-Seq-analysis
Analyzing RNA-Seq data typically involves several steps, from preprocessing raw reads to differential expression analysis. Here's a basic pipeline using common tools in a Linux environment.

## Prerequisites
Make sure you have the following tools installed:

* FastQC (for quality control)
* Trimmomatic (for trimming)
* STAR (for alignment)
* featureCounts (for counting reads)
* DESeq2 or edgeR (for differential expression analysis in R)

  ## Installation

* Install FastQC
```bash
  sudo apt-get update
  sudo apt-get install fastqc
```
* Install Trimmomatic

```bash
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
sudo mv Trimmomatic-0.39 /usr/local/bin/Trimmomatic
```

* Install STAR

  ```bash
    conda install -c bioconda star
  ```

 * Install featureCounts
   ```bash
     conda create -n featurecounts -c bioconda subread
   ```

*  Install R and Bioconductor for DESeq2 or edgeR 
  
  ```bash
 
   if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")}
 
     BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db"))
     install.packages(c("ggplot2", "pheatmap"))
 
 ```
