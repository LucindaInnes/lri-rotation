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
Here I made a script (ONT_pipeline.sh) that has 2 parameters, 1) a folder that contains the sequencing_summary file and the fastq_pass folder, and 2) the reference genome.
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
I wrote another script (multi_ONT_pipeline.sh) which will run through a directory of input folders and run the script on each folder. You input 1) the directory of input folders, 2) the first script and 3) the reference genome.

```
source ~/.bashrc

cd $1

for DIR in */; do sbatch $2 $DIR $3; done
```

### Step 6: Analysis
IGV:
Use your BAM and BAI files to visualize coverege in the IGV software.

R:
I wrote an Rscript (summary_calc.R) to make coverage graphs and report the average and % coverage across the entire genome (this only works for simplified human genome I created in step 3). I have not yet figured out how to make a generalized script for this. 




