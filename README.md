GETA
============
GETA (Genome-wide Electronic Tool for Annotation) is a pipeline with improved accuracy in the annotation of eukaryotic genome. This software was written with Perl by Lianfu Chen, Hansheng Zhao, Congcong Liu and Lei Deng.

Introduction
============
GETA can integrate the various evidence from ab initio gene finding and experimental data, including homologous proteins and RNA-Seq data.

GETA is a multi-threaded, automatic, and memory-saving GWAT, which effectively improved the accuracy and completeness of predicted annotation. The current version is v1.0 and later version will support more types of data like the single-molecular sequencing data and add other analysis modules on the basis of continuous improvement in accuracy and completeness.

Installation
============
1. unpack

2. install dependencies:

    RepeatMasker and RepeatModeler
    Java (version: jre1.8.0_45)
    HISAT2 (version: 2.1.0)
    Samtools (version: 1.8)
    genewise (version: 2.4.1)
    NCBI-Blast+ (version: 2.6.0+)
    AUGUSTUS (version: 3.3.1) *Note, bam2hints should be compiled correctly*
    HMMER (version: 3.1b2)
    ParaFly

3. add these directories of the executables to the PATH environment variable

4. the detail command lines of installation of all the dependency softwares and GETA: INSTALL.

Usage of the main script geta.pl
=================

    GETA (Genome-wide Electronic Tool for Annotation) is a pipeline software for predicting gene models
     from whole genome sequences. With one command, you can quickly obtain an accurate genome annotation
     GFF3 result file by providing the genome sequence, RNA-Seq raw data, and whole genome homologous protein
     sequences of closely related species. This program has two outstanding features: (1) accurate prediction,
     GETA outputs gene models with the normal amount, high BUSCO integrity, and accurate exon boundaries;
     (2) simple run, the operation is straightforward with a fully automated command that produces the final
     result.

    Current Version: 2.7.1

    Usage:
        perl /opt/biosoft/geta/bin/geta.pl [options]

    For example:
        perl /opt/biosoft/geta/bin/geta.pl --genome genome.fasta --RM_species_Dfam Embryophyta --RM_species_RepBase
     Embryophyta --pe1 libA.1.fq.gz,libB.1.fq.gz --pe2 libA.2.fq.gz,libB.2.fq.gz --protein homolog.fasta
     --augustus_species GETA_genus_species --HMM_db /opt/biosoft/bioinfomatics_databases/Pfam/Pfam-A.hmm
     --config /opt/biosoft/geta/conf_for_big_genome.txt --out_prefix out --gene_prefix GS01Gene --cpu 120


    Parameters:
    [INPUT]
        --genome <string>    default: None, required
        Enter the FASTA file of the genome sequence that you want to annotate. If the input genome file
     has repeats masked, you can skip the repeat sequence masking step by removing the --RM_species_Dfam,
      --RM_species_RepBase, --RM_lib parameters and adding the --no_RepeatModeler parameter. In this instance,
     it is advised to hard mask the transposable sequences using base N and to soft mask simple or tandem
     repeats using lowercase characters.

        --RM_species_Dfam | --RM_species <string>    default: None
        Enter the name of a species or class for RepeatMakser to perform a repeat sequence analysis for
     the genome using the HMM data of corresponding taxonomic species in the Dfam database. The file /opt
    /biosoft/geta/RepeatMasker_lineage.txt has the values that can be provided for this parameter to represent
     the class of species. For example, Eukaryota is for eukaryotes, Viridiplantae is for plants, Metazoa
     is for animals, and Fungi is for fungi. Before attempting to enter this parameter, the RepeatMasker
     program needed to be installed and the Dfam database needed to be configured. Note that due to its massive
     size, the Dfam database has been split up into nine partitions. By default, RepeatMasker only contains
     the zero root partition of the Dfam database, which is suitable for species such as mammals and fungi.
     If necessary, consider downloading the proper Dfam database partition and configuring it to RepeatMasker.
     For example, the 5th partition of the Dfam database is designated for Viridiplantae, while the 7th
     partition is for Hymenoptera.

        --RM_species_RepBase <string>    default: None
        Enter the name of a species or class for RepeatMakser to perform a repeat sequence analysis for
     the genome using the nucleotide sequences of corresponding taxonomic species in the RepBase database.
     The file /opt/biosoft/geta/RepeatMasker_lineage.txt has the values that can be provided for this parameter
     to represent the class of species. For example, Eukaryota is for eukaryotes, Viridiplantae is for plants,
     Metazoa is for animals, and Fungi is for fungi. Before attempting to enter this parameter, the RepeatMasker
     program needed to be installed and the RepBase database needed to be configured. Note that RepBase is
     no longer providing free downloads and that the most recent version of the database, 20181026, is older
     and contains few repetitive sequence data.

        --RM_lib <string>    default: None
        Enter a FASTA file and use the repetitive sequence to conduct genome-wide repeat analysis. This
     file is usually the output of RepeatModeler software's analysis of the entire genome sequence, indicating
     the repeated sequences across the genome. By default, the GETA program calls RepeatModeler to look
     up the entire genome sequence and acquire the species' repetitive sequence database. RepeatMakser
     is then called to search the repeated sequences. After adding this argument, the time-consuming RepeatModler
     step is skipped, which may significantly reduce the running time of the program. Additionally, the software
     supports the simultaneous use of the --RM_species_Dfam, --RM_species_RepBase, and --RM_lib arguments,
     so that multiple methods can be used for repeat sequence analysis, and eventually multiple results
     can be combined and the result of any method can be recognized.

        --no_RepeatModeler    default: None
        When this parameter is added, the program will no longer run the RepeatModeler step, which is
     suitable for cases where the repeats have been masked in the input genome file.

        --pe1 <string> --pe2 <string>    default: None
        Enter one or more pairs of FASTQ format files from Paired-End next-generation sequencing technology.
     This parameter supports the input of multiple pairs of FASTQ files, using commas to separate the FASTQ
     file paths of different libraries. This parameter also accepts compressed files in .gz format.

        --se <string>    default: None
        Enter one or more FASTQ format files from Single-End next-generation sequencing technology. This
     parameter supports the input of multiple Single-End FASTQ files, using commas to separate the FASTQ
     file paths of different libraries. This parameter also accepts compressed files in .gz format.

        --sam <string>    default: None
        Enter one or more SAM format files from the output of alignment software such as HISTA2. This
     parameter supports the input of multiple SAM files, using commas to separate the SAM file paths. This
     parameter also accepts compressed files in .bam format. In addition, the program allows for the full
     or partial use of the three parameters --pe1/--pe2, --se, and --sam, then all of the input data are
     used for genome alignment to generate the transcript sequence for gene model prediction.

        --strand_specific    default: None
        When this parameter is added, all input next-generation sequencing data are treated as strand
    -specific, and the program will predict gene models only on the forward strand of the transcript. When
     two neighboring genes overlap in the genome, strand-specific sequencing data and this parameter can
     help accurately estimate gene borders.

        --protein <string>    default: None
        Enter a FASTA file containing whole genome protein sequences from neighboring species. It is
     recommended to use whole genome homologous protein sequences from 3 ~ 10 different species. It is also
     recommended to modify the name of the protein sequence by appending the Species information, which begins
     with the species character, to the end of its original name. For example, if the protein sequence is
     XP_002436309.2, it will be better renamed XP_002436309_2_SpeciesSorghumBicolor. In this way, it is beneficial
     to retain the homologous matching results of more species in a gene region, improving gene prediction
     accuracy. The more species employed, the more accurate gene models may be predicted, but the computational
     time required increases. Note that evidence-supported gene prediction requires at least one type of
     homologous protein or next-generation sequencing data.

        --augustus_species <string>    default: None
        When an AUGUSTUS species name is provided, the program starts from an existing species model
     or retrains a new species model when performing AUGUSTUS Training using gene models predicted by transcripts
     or homologous proteins. If the input AUGUSTUS species model exists, its parameters will be optimized.
     If not, a new AUGUSTUS species HMM model will be trained and then its parameters will be optimized.
     The AUGUSTUS Training step requires the installation of AUGUSTUS software and configuration of the
     $AUGUSTUS_CONFIG_PATH environment variable. A species configuration folder from AUGUSTUS Training with
     the name provided in this parameter is generated in the temporary folder following the program's successful
     execution. If the user executing the program has write access, the produced species configuration folder
     can be copied to the species folder specified in $AUGUSTUS_CONFIG_PATH. If you do not enter this parameter,
     the program will automatically set the value of this parameter to "GETA + prefix of genome FASTA file
     name + date + process ID".

        --HMM_db <string>    default: None
        Enter one or more HMM databases, for filtering gene models. This parameter supports the input
     of multiple HMM databases, separated by commas. The program filters those gene models that do not match
     in all databases when using multiple HMM databases.

        --BLASTP_db <string>    default: None
        Enter one or more diamond databases, for filtering gene models. This parameter supports the input
     of multiple diamond databases, separated by commas. The program filters those gene models that do not
     match in all databases when using multiple diamond databases. When this parameter is left unset, the
     homologous proteins provided by the --protein parameter will be used to build the diamond database
     for filtering gene models.

        --config <string>    default: None
        Enter a parameter profile path to set the detailed parameters of other commands called by this
     program. If this parameter is left unset, When the genome size exceeds 1GB, the software installation
     directory's conf_for_big_genome.txt configuration file is automatically used. conf_for_small_genome.txt
     for genome size < 50MB, conf_all_defaults.txt for genome size between 50MB and 1GB.  Additionally, the
     thresholds for filtering the gene models typically need to be adjusted when GETA predicts an abnormally
     high number of genes. Then, the GETA pipeline can be rerun by setting this parameter to a new configuration
     file that is made by modifying the contents of the conf_all_defaults.txt file in the software installation
     directory.

        --BUSCO_lineage_dataset <string>    default: None
        Enter one or more BUSCO databases, the program will additionally perform BUSCO analysis on the
     whole genome protein sequences obtained by gene prediction. This parameter supports the input of multiple
     BUSCO databases, separated by commas. The information contained in the /opt/biosoft/geta/BUSCO_lineages_list.2021
    -12-14.txt file can be used to choose the proper BUSCO databases. Finally, the BUSCO results are exported
     to the 7.outputResults subdirectory and to the gene_prediction.summary file.

    [OUTPUT]
        --out_prefix <string>    default: out
        Enter a perfix of output files or temporary directory.

        --gene_prefix <string>    default: gene
        Enter the gene name prefix in the output GFF3 files.

        --chinese_help    default: None
        display the chinese usage and exit.

        --help    default: None
        display this help and exit.

    [Settings]
        --cpu <int>    default: 4
        Enter the number of CPU threads used by GETA or the called programes to run.

        --genetic_code <int>    default: 1
        Enter the genetic code. The values for this parameter can be found on the NCBI Genetic Codes
     website at: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi. This parameter is mainly effective
     for the gene prediction steps through homologous proteins, as well as the situation where start and
     stop codon information is used for filling the end of incomplete gene models.

        --homolog_prediction_method <string>    default: all
        Enter a method for gene prediction using homologous proteins. The value can be set to exonerate,
     genewise, gth, or all. This parameter supports the input of multiple methods, separated by commas. If
     the value was set to all, it indicates all three methods were used. The more methods you use, the more
     computation time you consume, but the better the result will be. Of the three methods, exonerate and
     genewise produced similar accuracy results, but gth showed a significant decrease in sensitivity and
     a significant increase in specificity. The following table shows the accuracy of the prediction results
     for the Oryza sativa genome using three methods. We compared the annotation results of 28736 gene models
     on NCBI to assess four accuracy metrics: gene level sensitivity, gene level specificity, exon level
     specificity, and exon level specificity. It is obvious that using multiple methods for gene prediction,
     combining results, and then filtering can result in a closer number of gene models to the actual number
     of genes and more accurate results. In addition, this parameter has a higher priority and can override
     the homolog_prediction parameter value in the parameter configuration file specified by --config.
        Method     Gene_num    gene_sensitivity    gene_specificity    exon_sensitivity    exon_specificity

        exonerate  38537       47.05%              35.09%              58.90%              73.55%
        genewise   40455       47.32%              33.61%              62.08%              71.27%
        gth        8888        19.80%              64.02%              30.43%              90.54%
        all        40538       48.54%              34.41%              63.85%              71.86%
        filtered   28184       45.62%              46.51%              61.26%              79.59%

        --optimize_augustus_method <int>    default: 1
        Enter the method for AUGUSTUS parameters optimization. 1, indicates that only BGM2AT.optimize_augustus
     is called for AUGUSTUS optimization, which can use all CPU threads to parallel test all parameters and
     is fast. 2, means that the script optimize_augustus.pl provided by the AUGUSTUS software was called
     sequentially for AUGUSTUS optimization after BGM2AT.optimize_augustus had finished its run. This second
     step is much slower, but probably more effective. This parameter has a higher priority and can override
     the BGM2AT parameter value in the parameter configuration file specified by --config.

        --no_alternative_splicing_analysis    default: None
        When this parameter is added, the program does not perform alternative splicing analysis. Note
     that GETA defaults to perform alternative splicing analysis based on intron and base sequencing depth
     information when NGS reads were input.

        --delete_unimportant_intermediate_files    defaults: None
        When this parameter is added and the program runs successfully, the insignificant intermediate
     files are deleted, leaving only the minimal, small, and important intermediate result files.


    This software has been tested and successfully run on Rocky 9.2 system using the following dependent
     software versions:

    01. ParaFly
    02. RepeatMasker (version: 4.1.5)
    03. RepeatModeler (version: 2.0.4)
    04. makeblastdb/rmblastn/tblastn/blastp (Version: 2.14.0)
    05. java (version: 1.8.0_282)
    06. hisat2 (version: 2.1.0)
    07. samtools (version: 1.17)
    08. mmseqs (version 15-6f452)
    09. genewise (version: 2.4.1)
    10. gth (Vesion: 1.7.3)
    11. exonerate (Vesion: 2.2.0)
    12. augustus/etraining (version: 3.5.0)
    13. diamond (version 2.1.8)
    14. hmmscan (version: 3.3.2)
    15. busco (Version: 5.4.7)

    Version of GETA: 2.7.1

Test of GETA
============
we have compared the accuracy of GETA with other methods according to several species. the command line and output GFF3 files can be found at website: <a href="http://122.205.95.116/geta/" target="_noblank">122.205.95.116/geta/</a>
