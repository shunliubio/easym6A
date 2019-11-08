## easym6A: Process m6A/MeRIP-seq data in single or batch job mode

[![release](https://img.shields.io/badge/release-v1.0-orange.svg)](https://img.shields.io/badge/release-v1.0-orange.svg)
[![license](https://img.shields.io/badge/license-GPLv3-green.svg)](https://img.shields.io/badge/license-GPLv3-green.svg)

easym6A creates a bash script to process m6A/MeRIP-seq data from adapter trimming to peak calling. It can also serve as a pipeline for RNA-seq data.

### Dependencies:

- **Cutadapt**, tested with 1.15
- **Samtools**, tested with 1.7
- **HISAT2**, tested with 2.1.0
- **Picard**, tested with 2.17.10
- **bedtools**, tested with 2.27.1
- **bedGraphToBigWig**
- **StringTie**, tested with 1.3.4d
- **prepDE.py**
- **gffcompare**, tested with 0.10.4
- **MACS2**, tested with 2.1.1.20160309
- **HOMER**, tested with 4.9
- **R**, tested with 3.3.3. Install three libraries **exomePeak**, **MeTPeak** and **MeTDiff** in R.

### Installation:

Add executable directories of above dependencies, easym6A.pl and 3peakSuit.R to PATH.

### Usage:

#### 1. Prepare the index for the reference genome and repetitive elements.

easym6A uses HISAT2 to map reads to the genome. One genome sequence file in fasta format and one genome annotation file in gtf format are required to build the index for the genome. More specification can be referred in the [HISAT2 manual](https://ccb.jhu.edu/software/hisat2/manual.shtml#the-hisat2-build-indexer).

```
# Given that hg38.fa and gencode.v27.annotation.gtf are in the directory geonome and annotation, respectively.

# Generate index files into genome/hg38_tran_index
mkdir -p genome/hg38_tran_index
# Extract exons from the gtf file
extract_exons.py annotation/gencode.v27.annotation.gtf > genome/hg38_tran_index/gencode.v27.annotation.exon
# Extract splicing sites from the gtf file
extract_splice_sites.py annotation/gencode.v27.annotation.gtf > genome/hg38_tran_index/gencode.v27.annotation.ss
# Build the index for the genome
hisat2-build --ss genome/hg38_tran_index/gencode.v27.annotation.ss --exon genome/hg38_tran_index/gencode.v27.annotation.exon genome/hg38.fa genome/hg38_tran_index/hg38_tran
```

As easym6A can filter reads mapped to repetitive elements, one repetitive element sequence file in fasta format is required to build the index for repetitive elements, unless the filter step is ignored.

```
# Given that one would like to filter reads mapped to rRNAs. hg38_rRNA.fa is in the directory genome.

# Generate index files into genome/hg38_rRNA_index
mkdir -p genome/hg38_rRNA_index
# Build the index for rRNAs
hisat2-build genome/hg38_rRNA.fa genome/hg38_rRNA_index/hg38_rRNA
```

#### 2. Edit one table file including m6A/MeRIP-seq sample information

The table file is consisted of 19 fields, which records basic information of samples such as their file names, file paths, library types, etc. One row for one sample. Here is an example

| Species | Cell Line | Dataset | Sample | Group | File Name | Fastq File Dir | Fastq File | Bam File Path | 5’ Adapter | 3’ Adapter | 5’ Barcode | 3’ Barcode | Q33 | Strandness | Fragment Length | Read Length | Seq Layout | ID |
|---------|-----------|---------|----------------|-----------------|----------------------|--------------------------|--------------------------------|----------------------------------------------|------------|------------|------------|------------|-----|------------|-----------------|-------------|------------|----|
| human | HeLa | - | siControl_rep1 | control_input | input_siControl_rep1 | /your/raw/data/file/path | input_siControl_chr22.fastq.gz | /your/bam/file/path/input_siControl_rep1.bam | - | - | - | - | N | R | 150 | 50 | SINGLE | 1 |
| human | HeLa | - | siControl_rep1 | control_ip | ip_siControl_rep1 | /your/raw/data/file/path | ip_siControl_chr22.fastq.gz | /your/bam/file/path/ip_siControl_rep1.bam | - | - | - | - | N | R | 150 | 50 | SINGLE | 2 |
| human | HeLa | - | siMETTL3_rep1 | treatment_input | input_siMETTL3_rep1 | /your/raw/data/file/path | input_siMETTL3_chr22.fastq.gz | /your/bam/file/path/input_siMETTL3_rep1.bam | - | - | - | - | N | R | 150 | 50 | SINGLE | 3 |
| human | HeLa | - | siMETTL3_rep1 | treatment_ip | ip_siMETTL3_rep1 | /your/raw/data/file/path | ip_siMETTL3_chr22.fastq.gz | /your/bam/file/path/ip_siMETTL3_rep1.bam | - | - | - | - | N | R | 150 | 50 | SINGLE | 4 |

```
single job mode (call peaks):

$ m6Aseq_workflow.pl -s <sampleList.txt> -c <configure.txt> -n <1,2> -e <siControl> [options]

single job mode (call differential peaks):

$ m6Aseq_workflow.pl -s <sampleList.txt> -c <configure.txt> -n <1,2,3,4> -e <siMETTL3> [options]

batch job mode:

$ m6Aseq_workflow.pl -s <sampleList.txt> -b <batchList.txt> -x <1> -y <10> [options]

 Options:
   -h/--help                   brief help message
   --man                       full documentation
   -s/--samplelist <file>      a specific-formated file that records sample infomation. (Required)
   -c/--configure <file>       a specific-formated file that records ouput path. (Requried)
   -n/--sampleno <string>      a comma-seperated ID list of samples that you want to analyze. (Required)
   -e/--runname <string>       a user-defined job name. (Required)
   -t/--threads <num>          number of cpu threads to run the program. (Default: 1)
   -l/--parallel               reduce the analysis time in parallel mode. (Default: off)
   -m/--method <string>        the peak calling tool(s) you want to use. (Default: all)
   -b/--batch <file>           a specific-formated file that records batch job infomation.
   -x/--bstart <int>           batch start ID. Required when -b/--batch is set. Used together with -y/--bend.
   -y/--bend <int>             batch end ID. Required when -b/--batch is set. Used together with -x/--bstart.
   -a/--onlybam                run gene expression analysis only. like RNA-seq pipleline. (Default: off)
   -k/--onlypeak               run peak calling analysis only. (Default: off)
   -p/--rmrep                  remove repetitive elements. (Default: off)
   -d/--rmdup                  remove PCR duplicates. (Default: off)
   -f/--keeptmp                keep intermediated files. (Default: off)
   -r/--run                    run the bash script(s) generated by the program. (Default: off)
```

