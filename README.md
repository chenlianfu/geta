# GETA
GETA is an automatic genome-wide annotation tool (GWAT) with improved accuracy and gene integrity for eukaryotes written by Lianfu Chen, Congcong Liu, Lei Deng and so on.

## 1. introduction
GETA can integrate the various evidence from ab initio gene finding and experimental data, including homologous proteins and RNA-Seq data.

GETA is a multi-threaded, automatic, and memory-saving GWAT, which effectively improved the accuracy and completeness of predicted annotation. The current version is v1.0 and later version will support more types of data like the single-molecular sequencing data and add other analysis modules on the basis of continuous improvement in accuracy and completeness.

## 2. installation
1.unpack

2.install dependencies:
    ParaFly
    java 1.8.0_144
    HISAT2 2.1.0
    samtools 1.3.1
    hmmer 3.1b2
    Augustus

3.add these directories of the executables to the PATH environment variable

## 3. the usage of the main script geta.pl
Usage:

    perl ./geta.pl [options]

For example:

    perl ./geta.pl --out_prefix out -1 reads.1 fastq -2 reads.2.fastq --protein homolog.fasta --cpu 80 --hisat2 " --min-intronlen 20 --max-intronlen 20000 --rna-strandness RF" --strand_specific --sam2transfrag " --fraction --min_expressed_base_depth 2 --max_expressed_base_depth 50 --min_junction_depth 2 --max_junction_depth 50" --species oryza_sativa_20171120 --pfam_db /opt/biosoft/hmmer-3.1b2/Pfam-AB.hmm --gene_prefix OS01Gene --genome genome.fasta

Parameters:

    --RM_species <string>    default: None
    species identifier for RepeatMasker.

    --genome <string>
    genome file in fasta format.

    --out_prefix <string>    default: out
    the prefix of outputs.
    
    -1 <string> -2 <string>
    fastq format files contain of paired-end RNA-seq data.

    -S <string>
    fastq format file contains of single-end RNA-seq data.

    --protein <string>
    homologous protein sequences (derived from multiple species would be recommended) file in fasta format.

    --cpu <int>    default: 4
    the number of threads.

    --trimmomatic <string>    default: "TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 TOPHRED33"
    the parameter of trimmomatic. The value must be surrounded by a pair of quotes.

    --hisat2 <int>    default: " --min-intronlen 20 --max-intronlen 20000 --rna-strandness unstranded --dta --score-min L,0.0,-0.4"
    the parameter of hisat2. The value must be surrounded by a pair of quotes, and the first quote is followed by a white space.

    --strand_specific
    enable the ability of analysing the strand-specific information provided by the tag "XS" of alignments. If this parameter was set, the paramter "--rna-strandness" should be set to "RF" usually, so the value of "--hisat2" would be " --min-intronlen 20 --max-intronlen 20000 --rna-strandness RF --dta --score-min L,0.0,-04".

    --sam2transfrag <string>    default: " --fraction 0.05 --min_expressed_base_depth 2 --max_expressed_base_depth 50 --min_junction_depth 2 --max_junction_depth 50"
    the parameter of sam2transfrag. The value must be surrounded by a pair of quotes, and the first quote is followed by a white space.
    the defalut value means: the dynamic coverage threshold of mapping regions (come frome SAM file) was determinated by the maximum base depth of each region * 0.05, as well as this threshold should between 2 and 50.

    --ORF2bestGeneModels <string>    default: " --min_cds_num 3 --min_cds_length 900 --min_cds_exon_ratio 0.60 --intron_length_fractile 0.95 --cds_length_fractile 0.95"
    the parameter of ORF2bestGeneModels. The value must be surrounded by a pair of quotes, and the first quote is followed by a white space.
    the defalut value means: a good gene model shold be completed, contain at least 3 CDS sequences whose total length shold >= 900bp, the length ratio of CDS/exon should >= 60%, intron and cds length should not great than a auto-calculated threashold.

    --species <string>
    species identifier for Augustus.  

    --pfam_db <string>    default: "/opt/biosoft/hmmer-3.1b2/Pfam-AB.hmm"
    the absolute path of Pfam database which was used for filtering of false positive gene models.

    --gene_prefix <string>    default: gene
    the prefix of gene id shown in output file.

This script was tested on CentOS 6.8 with such softwares can be run directly in terminal:
1. ParaFly
2. java version "1.8.0_144"
3. HISAT2 version 2.1.0
4. samtools Version: 1.3.1
5. hmmer-3.1b2

Version: 1.0

