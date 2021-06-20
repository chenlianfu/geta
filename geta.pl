#!/usr/bin/perl
use strict;
use Getopt::Long;
use File::Basename;

my $usage = <<USAGE;
Usage:
    perl $0 [options]

For example:
    perl $0 --RM_species Embryophyta --genome genome.fasta -1 liba.1.fq.gz,libb.1.fq.gz -2 liba.2.fq.gz,libb.2.fq.gz --protein homolog.fasta --augustus_species oryza_sativa_20171120 --out_prefix out --config conf.txt --cpu 80 --gene_prefix OS01Gene --pfam_db /opt/biosoft/hmmer-3.1b2/Pfam-AB.hmm

Parameters:
[required]
    --RM_species <string>
    species identifier for RepeatMasker.

    --genome <string>
    genome file in fasta format.

    -1 <string> -2 <string>
    fastq format files contain of paired-end RNA-seq data. if you have data come from multi librarys, input multi fastq files separated by comma. the compress file format .gz also can be accepted.

    -S <string>
    fastq format file contains of single-end RNA-seq data. if you have data come from multi librarys, input multi fastq files separated by comma. the compress file format .gz also can be accepted.

    --protein <string>
    homologous protein sequences (derived from multiple species would be recommended) file in fasta format.

    --augustus_species <string>
    species identifier for Augustus. the relative hmm files of augustus training will be created with this prefix. if the relative hmm files of augustus training exists, the program will delete the hmm files directory firstly, and then start the augustus training steps.

    --use_existed_augustus_species <string>
    species identifier for Augustus. This parameter is conflict with --augustus_species. When this parameter set, the --augustus_species parameter will be invalid, and the relative hmm files of augustus training should exists, and the augustus training step will be skipped (this will save lots of runing time).

[other]
    --out_prefix <string>    default: out
    the prefix of outputs.

    --config <string>    default: None
    Input a file containing the parameters of several main programs (such as trimmomatic, hisat2 and augustus) during the pipeline. If you do not input this file, the default parameters should be suitable for most situation.
    
    --RM_lib <string>    default: None
    A fasta file of repeat sequences. Generally to be the result of RepeatModeler. If not set, RepeatModeler will be used for producing this file automaticly, which shall time-consuming.

    --cpu <int>    default: 4
    the number of threads.

    --strand_specific    default: False
    enable the ability of analysing the strand-specific information provided by the tag "XS" from SAM format alignments. If this parameter was set, the paramter "--rna-strandness" of hisat2 should be set to "RF" usually.

    --pfam_db <string>    default: None
    the absolute path of Pfam database which was used for filtering of false positive gene models.

    --gene_prefix <string>    default: gene
    the prefix of gene id shown in output file.

    --no_augustus_training_iteration    default: False


This script was tested on CentOS 6.8 with such softwares can be run directly in terminal:
1. ParaFly
2. java (version: 1.8.0_172)
3. hisat2 (version: 2.1.0)
4. samtools (version: 1.8)
5. hmmscan (version: 3.1b2)
6. makeblastdb/tblastn/blastp (version: 2.6.0)
7. RepeatMasker (version: 4.0.7)
8. RepeatModeler (version: 1.0.11)
9. genewise (version: 2.4.1)
10. augustus/etraining (version: 3.3.1)

Version: 2.4.12

USAGE
if (@ARGV==0){die $usage}

my ($RM_species, $RM_lib, $genome, $out_prefix, $pe1, $pe2, $single_end, $protein, $cpu, $trimmomatic, $strand_specific, $sam2transfrag, $ORF2bestGeneModels, $augustus_species, $pfam_db, $gene_prefix, $cmdString, $no_augustus_training_iteration, $config, $use_existed_augustus_species);
GetOptions(
    "RM_species:s" => \$RM_species,
    "RM_lib:s" => \$RM_lib,
    "genome:s" => \$genome,
    "out_prefix:s" => \$out_prefix,
    "1:s" => \$pe1,
    "2:s" => \$pe2,
    "S:s" => \$single_end,
    "protein:s" => \$protein,
    "cpu:s" => \$cpu,
    "strand_specific!" => \$strand_specific,
    "augustus_species:s" => \$augustus_species,
    "use_existed_augustus_species:s" => \$use_existed_augustus_species,
    "pfam_db:s" => \$pfam_db,
    "gene_prefix:s" => \$gene_prefix,
    "no_augustus_training_iteration!" => \$no_augustus_training_iteration,
    "config:s" => \$config,
);

# 检测依赖的软件
print STDERR "\n============================================\n";
print STDERR "Detecting the dependency softwares:\n";
# 检测RepeatMasker
my $software_info = `RepeatMasker`;
if ($software_info =~ m/RepeatMasker/) {
    print STDERR "RepeatMasker:\tOK\n";
}
else {
    die "RepeatMasker:\tFailed\n\n";
}
# 检测RepeatModeler
$software_info = `RepeatModeler`;
if ($software_info =~ m/RepeatModeler/) {
    print STDERR "RepeatModeler:\tOK\n";
}
else {
    die "RepeatModeler:\tFailed\n\n";
}
# 检测ParaFly
$software_info = `ParaFly 2>&1`;
if ($software_info =~ m/Usage: ParaFly/) {
    print STDERR "ParaFly:\tOK\n";
}
else {
    die "ParaFly:\tFailed\n\n";
}
# 检测JAVA
my $software_info = `java -version 2>&1`;
if ($software_info =~ m/version \"(1.(\d).*?)\"/) {
    if ($2 == 8) {
        print STDERR "java:\tOK\n";
    }
    else {
        print STDERR "java:\tthis java version $1 may not work properly, version 1.8 is desired\n";
    }
}
else {
    die "java:\tFailed\n\n";
}
# 检测HISAT2
$software_info = `hisat2 --version`;
if ($software_info =~ m/version 2.(\d+)\.(\d+)/) {
    if ($1 >= 1) {
        print STDERR "HISAT2:\tOK\n";
    }
    else {
        print STDERR "HISAT2:\tthis HISAT2 version 2.$1.$2 may not work properly, version 2.1.0 is desired\n";
    }
}
else {
    die "HISAT2:\tFailed\n\n";
}
# 检测samtools
$software_info = `samtools --version`;
if ($software_info =~ m/samtools 1.(\d+)/) {
    if ($1 >= 3) {
        print STDERR "samtools:\tOK\n";
    }
    else {
        print STDERR "samtools:\tthis samtools version 1.$1.$2 may not work properly, version 1.3.1 is desired\n";
    }
}
else {
    die "samtools:\tFailed\n\n";
}
# 检测hmmer
$software_info = `hmmscan -h`;
if ($software_info =~ m/HMMER 3.(\d+)/) {
    if ($1 >= 1) {
        print STDERR "hmmer:\tOK\n";
    }
    else {
        print STDERR "hmmer:\tthis hmmer version 3.$1 may not work properly, version 3.1b2 is desired\n";
    }
}
else {
    die "hmmer:\tFailed\n\n";
}

# 参数设置
# 检查输入文件
die "--RM_species shoud be set!" unless $RM_species;
my $pwd = `pwd`;
chomp($pwd);
die "No genome fasta input\n" unless $genome;
$genome =~ s/^/$pwd\// unless $genome =~ m/^\//;
$protein  =~ s/^/$pwd\// unless $protein =~ m/^\//;
unless (($pe1 && $pe2) or $single_end or $protein) {
    die "No RNA-Seq short reads or homologous proteins as input\n";
}
die "No Augustus species provided\n" unless ($augustus_species or $use_existed_augustus_species);
if ($use_existed_augustus_species) {
    my $species_config_dir = `echo \$AUGUSTUS_CONFIG_PATH`;
    chomp($species_config_dir);
    $species_config_dir = "$species_config_dir/species/$use_existed_augustus_species";
    if (-e $species_config_dir) {
        $augustus_species = $use_existed_augustus_species;
    }
    else {
        die "The AUGUSUTS HMM files of $use_existed_augustus_species does not exists!\n";
    }
}
if ($RM_lib) {
    $RM_lib =~ s/^/$pwd\// unless $RM_lib =~ m/^\//;
}
if ($config) {
    $config =~ s/^/$pwd\// unless $config =~ m/^\//;
}
my (%pe_reads, %se_reads);
if ($pe1 && $pe2) {
    my @pe1 = split /,/, $pe1;
    my @pe2 = split /,/, $pe2;
    my $pe1_num = @pe1;
    my $pe2_num = @pe2;
    if ($pe1_num != $pe2_num) { die "the input file number of -1 was not equal to -2.\n" };
    foreach (@pe1) {
        s/^/$pwd\// unless m/^\//;
        my $pe_file = $_;
        $_ = shift @pe2;
        s/^/$pwd\// unless m/^\//;
        $pe_file .= "\t$_";
        $pe_reads{$pe_file} = 1;
    }
}
if ($single_end) {
    my @se = split /,/, $single_end;
    foreach (@se) {
        s/^/$pwd\// unless m/^\//;
        $se_reads{$_} = 1;
    }
}

$out_prefix ||= "out";
$cpu ||= 4;
# 各个主要命令的参数设置
my %config = (
    'RepeatMasker' => '-e ncbi -gff',
    'trimmomatic' => 'TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 TOPHRED33',
    'hisat2-build' => '-p 1',
    'hisat2' => '--min-intronlen 20 --max-intronlen 20000 --dta --score-min L,0.0,-0.4',
    'sam2transfrag' => '--fraction 0.05 --min_expressed_base_depth 2 --max_expressed_base_depth 50 --min_junction_depth 2 --max_junction_depth 50 --min_fragment_count_per_transfrags 10 --min_intron_length 20',
    'TransDecoder.LongOrfs' => '-m 100 -G universal',
    'TransDecoder.Predict' => '--retain_long_orfs_mode dynamic',
    'homolog_genewise' => '--coverage_ratio 0.4 --evalue 1e-9',
    'homolog_genewiseGFF2GFF3' => '--min_score 15 --gene_prefix genewise --filterMiddleStopCodon',
    'geneModels2AugusutsTrainingInput' => '--min_evalue 1e-9 --min_identity 0.8 --min_coverage_ratio 0.6 --min_cds_num 2 --min_cds_length 600 --min_cds_exon_ratio 0.60',
    'BGM2AT' => '--min_gene_number_for_augustus_training 500 --gene_number_for_accuracy_detection 200 --min_gene_number_of_optimize_augustus_chunk 50 --max_gene_number_of_optimize_augustus_chunk 200',
    'prepareAugusutusHints' => '--margin 20',
    'paraAugusutusWithHints' => '--gene_prefix augustus --min_intron_len 30 --alternatives_from_evidence',
    'paraCombineGeneModels' => '--overlap 30 --min_augustus_transcriptSupport_percentage 10.0 --min_augustus_intronSupport_number 1 --min_augustus_intronSupport_ratio 0.01',
    'PfamValidateABinitio' => '--CDS_length 750 --CDS_num 2 --evalue 1e-5 --coverage 0.25',
    'remove_genes_in_repeats' => '--ratio 0.8',
    'remove_short_genes' => '--cds_length 150',
);
if ($config) {
    open IN, $config or die "Can not open file $config, $!\n";
    my $tag;
    while (<IN>) {
        next if m/^#/;
        next if m/^\s*$/;
        s/^\s+//;
        s/\n/ /;
        if (/\[(.*)\]/) { $tag = $1; delete $config{$1}; }
        else { $config{$tag} .= $_; }
    }
    close IN;
}

my $dirname = dirname($0);
$dirname =~ s/bin$//;

# 生成临时文件夹
mkdir "$out_prefix.tmp" unless -e "$out_prefix.tmp";
chdir "$out_prefix.tmp";
$pwd = `pwd`; print STDERR "PWD: $pwd";
unless (-e "genome.fasta") {
    open OUT, ">", "genome.fasta" or die "Can not create file genome.fasta, $!\n";
    open IN, $genome or die "Can not open file $genome, $!\n";
    $_ = <IN>;
    print OUT;
    while (<IN>) {
        if (m/^>/) {
            print OUT "\n$_";
        }
        else {
            s/\s+?$//g;
            print OUT;
        }
    }
    close IN;
    close OUT;
}
$pwd = `pwd`; chomp($pwd);
$genome = "$pwd/genome.fasta";

# Step 0: RepeatMasker and RepeatModeler
print STDERR "\n============================================\n";
print STDERR "Step 0: RepeatMasker and RepeatModeler " . "(" . (localtime) . ")" . "\n";
mkdir "0.RepeatMasker" unless -e "0.RepeatMasker";
unless (-e "0.RepeatMasker.ok") {
    chdir "0.RepeatMasker";
    $pwd = `pwd`; print STDERR "PWD: $pwd";
    
    # 进行RepeatMasker分析
    mkdir "repeatMasker" unless -e "repeatMasker";
    my $cpu_RepeatMasker = int($cpu / 4);
    $cmdString = "RepeatMasker $config{'RepeatMasker'} -pa $cpu_RepeatMasker -species $RM_species -dir repeatMasker/ $genome &> repeatmasker.log";
    unless (-e "RepeatMasker.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "RepeatMasker.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    # 进行RepeatModeler分析
    mkdir "repeatModeler" unless -e "repeatModeler";
    chdir "repeatModeler";
    $pwd = `pwd`; print STDERR "PWD: $pwd";
    unless ($RM_lib) {
        $cmdString = "BuildDatabase -name species -engine ncbi $genome";
        unless (-e "BuildDatabase.ok") {
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            open OUT, ">", "BuildDatabase.ok" or die $!; close OUT;
        }
        else {
            print STDERR "CMD(Skipped): $cmdString\n";
        }
    my $cpu_RepeatModeler = int($cpu / 4);
        $cmdString = "RepeatModeler -pa $cpu_RepeatModeler -database species -LTRStruct &> RepeatModeler.log";
        unless (-e "RepeatModeler.ok") {
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            open OUT, ">", "RepeatModeler.ok" or die $!; close OUT;
        }
        else {
            print STDERR "CMD(Skipped): $cmdString\n";
        }
        $cmdString = "RepeatMasker $config{'RepeatMasker'} -pa $cpu -lib RM_\*/\*.classified -dir ./ $genome &> repeatmasker.log";
    }
    else {
        $cmdString = "RepeatMasker $config{'RepeatMasker'} -pa $cpu -lib $RM_lib -dir ./ $genome &> repeatmasker.log";
    }
    unless (-e "RepeatMasker.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "RepeatMasker.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }
    chdir "../";
    $pwd = `pwd`; print STDERR "PWD: $pwd";

    # 合并RepeatMasker和RepeatModeler的结果
    $cmdString = "$dirname/bin/merge_repeatMasker_out.pl $genome repeatMasker/genome.fasta.out repeatModeler/genome.fasta.out > genome.repeat.stats";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    $cmdString = "$dirnam