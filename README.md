# lri-rotation

# Oxford Nanopore Data Analysis Pipeline

### Step 1: QC of raw data

When you get the data from minION, you want to run some precursory QC on the data. There are a number of methods to use here: FastQC, PycoQC and MinION_QC. Here I use PycoQC. In order to use any of these, you must have the sequencing_summary.txt file that was output with the sequencing data.

```
# install PycoQC, and verify that it installed properly.
$ pip install pycoQC --user
$ pycoQC -h
# if you run into warnings, it may require you to add the pycoQC directory to the path

# next use pycoQC on the sequence_summary file. This will output an HTML file with interactive QC info.
$ pycoQC -f sequencing_summary.txt -o pycoQC_output.html

```

### Step 2: Trimming and Filtering

The data may still contain adaptor sequences, which are not part of the genome itself. To remove these you will want to use a program called Porechop. The input to porechop can be zipped or unzipped fasta or fastq files, or a directory of files. Porechop searches for the adaptor region and trims it out of all sequences. The final output is a single fasta or fastq file with all of the sequences.

```
# install porechop
$ conda activate /software/apps/porechop/0.2.3_seqan2.1.1

# trim off adaptors
$ porechop -i input_directory > output.fastq
```

Next you may want to filter out the low-quality reads. The suggested program is NanoFilt. With NanoFilt you can filter out reads by length and by quality, and you can trim off known low quality regions (usually the first handful of bases of each read). This step is not always used, as Guppy already does some basic filtering.

```
# install NanoFilt
$ pip install NanoFilt

# trim reads (<500 bp long, remove first 50 bp)
$ NanoFilt -l 500 --headcrop 50 input.fastq > output.fastq
```

### Step 3: Mapping to Reference genome

Download the fasta file of your reference genome from the internet. The human genome options are [here](https://www.ncbi.nlm.nih.gov/genome/guide/human/). This comes with scaffold or alternative chromosomes which may be useful or may clutter up the mapping, I just removed those with the commands below.

```
# download human genome GRCh38
$ wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz

# make .fa of NC chromosome sequences.
$ grep ">NC" GRCh38_latest_genomic.fna |sed "s/^.//g" > target.lst
  # troubleshoot installing seqtk, I had to add seqtk directory to my path
$ git clone https://github.com/lh3/seqtk.git;
$ cd seqtk; make
$ seqtk subseq GRCh38_latest_genomic.fna target.lst > human_NC_genome.fna
```

The actual mapping is performed using minimap2. Minimap2 takes a fastq or fasta file and a reference genome fasta file and maps all of the short reads to the reference genome. The output file is a SAM file.

This is a more computationally intensive step, so I would put this in a script to run in the background.

```
$ conda install minimap2
$ minimap2 -ax map-ont ref_genome.fna input.fastq > output_aln.sam
```

### Step 4: SAM/BAM File Formatting

SAM/BAM files are used for interpreting NGS data. However, some post-processing of the files is necessary before you can use it software. Firs you need to covert the SAM file to a binary BAM file. Then you sort the file in order by left-most coordinate. You may want to make a BAM index file too, this is useful for programs like IGV. Finally you want to get coverage summary using bedtools. And explanation of the summary is found at [here](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html).

```
# convert SAM -> BAM
$ module load samtools-1.13-gcc-8.2.0-nj2zx6n
$ samtools view -bT human_NC_genome.fna input.sam | samtools sort -o output.bam

# Make BAM index
$ samtools index input.bam output.bai

# Make a coverage summary file
$ conda install bedtools
$ genomeCoverageBed -ibam input.bam > coverage.txt
```

### Step 5: Generalized Script
Here I made a script that has 2 parameters, 1) a folder that contains the sequencing_summary file and the fastq_pass folder, and 2) the reference genome.
```
source ~/.bashrc

cd $1
mkdir output

# pycoQC
pycoQC -f sequencing_summary* -o output/pycoQC_output.html

# porechop
conda activate /software/apps/porechop/0.2.3_seqan2.1.1
porechop -i fastq_pass > chopped.fastq

# minimap2 + samTools
conda activate myenv
conda install minimap2
conda install bedtools
module load samtools-1.13-gcc-8.2.0-nj2zx6n

minimap2 -ax map-ont $2 chopped.fastq | samtools view -bT $2 | samtools sort -o output/aln.bam
samtools index output/aln.bam output/aln.bai
genomeCoverageBed -ibam output/aln.bam > output/aln_coverage.txt
```
I wrote another script which will run through a directory of input folders and run the script on each folder. You input 1) the directory of input folders and 2) the reference genome.
```
source ~/.bashrc

cd $1

for DIR in */; do sbatch ~/truncated_samples/gen_scripts/mapping.sh $DIR $2; done
```

### Step 6: Analysis
IGV:
Use your BAM and BAI files to visualize coverege in the IGV software.

R:
I wrote an Rscript to make coverage graphs and report the average and % coverage across the entire genome (this only works for simplified human genome I created in step 3). I have not yet figured out how to make a generalized script for this. 

I use the coverage file created by genomeCoverageBed. The columns in the file are as follows:
  1. chromosome
  2. depth of coverage from features in input file
  3. number of bases on chromosome with depth equal to column 2.
  4. size of chromosome in base pairs
  5. fraction of bases on chromosome with depth equal to column 2.

  I calculate the genome coverage in 2 ways. First, the average number of reads that cover each base, and second in the percent of bases covered at all.
  
  $$ Average \ Coverage = \frac{\sum (Base*\#\ of\ Reads)}{Total\ Bases}$$
  $$ Percent\ Coverage = (1 - \frac{Uncovered\ Bases}{Total\ Bases})*100\% $$

  The code creates a bargraph of the values per 

```
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)

# read in coverage file
coverage_NC <- read.table("path/to/output/aln_coverage.txt", header = FALSE, sep = "\t", col.names = c("chromosome", "depth", "numbases", "size", "percent_depth"))

# filter the rows that contain the word "genome," this is the values for the entire genome, not broken up by chromosome.
genome <- coverage_NC %>% filter(grepl('genome', chromosome))

# create a list with the ID of all the chromosomes in the data frame.
chr_index <- c("NC_000001", "NC_000002", "NC_000003", "NC_000004", "NC_000005", "NC_000006", "NC_000007", "NC_000008", "NC_000009", "NC_000010", "NC_000011", "NC_000012", "NC_000013", "NC_000014", "NC_000015", "NC_000016", "NC_000017", "NC_000018", "NC_000019", "NC_000020", "NC_000021", "NC_000022", "NC_000023", "NC_000024", "NC_012920", "genome")

# create a list of the reader-friendly names we will assign to the chromosomes
id_index <- (c((seq(1:22)),"X", "Y","m", "total"))

# create a for loop, where we search for the rows for each chromosome in the chr_index and then perform calculations to determine  average depth/base and the percent coverage.
avg_chr_cov <- vector()

for (i in seq(1:length(chr_index))){
  chro <- coverage_NC %>% filter(grepl(chr_index[i], chromosome))
  avg_chr_cov[i] <- sum(chro[,2]*chro[,3])/chro[1,4]
}

perc_chr_cov <- vector()

for (i in seq(1:length(chr_index))){
  chro <- coverage_NC %>% filter(grepl(chr_index[i], chromosome))
  perc_chr_cov[i] <- (1-chro[1,3]/chro[1,4])*100
}

# plot the data as bar plots
summary <- data.frame(id_index, avg_chr_cov, perc_chr_cov)
summary$id_index <- factor(summary$id_index, levels = summary$id_index)

avg_cov_in <- ggplot(data=summary, aes(x=id_index, y=avg_chr_cov)) +
  geom_bar(stat="identity") + ggtitle("Average Coverage by Chromosome (w mtDNA)") +
  xlab("Chromosome ID") + ylab("Average Coverage")
per_cov_in <- ggplot(data=summary, aes(x=id_index, y=perc_chr_cov)) +
  geom_bar(stat="identity")+ggtitle("Percent Coverage by Chromosome (w mtDNA)") +
  xlab("Chromosome ID") + ylab("Percent Coverage")

avg_cov_ex <- ggplot(data=summary[c(1:24,26),], aes(x=id_index, y=avg_chr_cov)) +
  geom_bar(stat="identity")+ ggtitle("Average Coverage by Chromosome")+
  xlab("Chromosome ID") + ylab("Average Coverage")
per_cov_ex <- ggplot(data=summary[c(1:24,26),], aes(x=id_index, y=perc_chr_cov)) +
  geom_bar(stat="identity")+ggtitle("Percent Coverage by Chromosome") +
  xlab("Chromosome ID") + ylab("Percent Coverage")

summary_plots <- list(avg_cov_in, per_cov_in, avg_cov_ex, per_cov_ex)

avg_total <- summary$avg_chr_cov[26]
per_total <- summary$perc_chr_cov[26]

# create a CSV file with the average and percent values for the whole genome.

totals <- as.data.frame(cbind(avg_total, per_total))
colnames(totals) <-(c("Average", "Percent"))

# save the figures and the CSV to the output folder
write.csv(totals, "/path/to/output/totals.csv", row.names = FALSE)

ggsave("/path/to/output/coverage_summary.pdf", width = 11, height = 7.5, marrangeGrob(summary_plots, nrow=2, ncol=2))
```

