#!/usr/bin/perl
use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw/abs_path getcwd cwd/;
use File::Basename;

my $command_line_geta = join " ", @ARGV;
my $dirname = dirname($0);
$dirname =~ s/\/bin$//;

my $usage = <<USAGE;
Usage:
    perl $0 [options] genome.fasta

For example:
    perl $0 --RM_species Embryophyta -1 liba.1.fq.gz,libb.1.fq.gz -2 liba.2.fq.gz,libb.2.fq.gz --protein homolog.fasta --augustus_species oryza_sativa_20220608 --out_prefix out --config conf.txt --cpu 80 --gene_prefix OS01Gene --HMM_db /opt/biosoft/hmmer-3.3.1/EggNOG_Eukaryota.hmm genome.fasta

Parameters:
[General]
    --pe1 <string> -pe2 <string>    default: None, Not Required but Recommened.
    fastq format files contain of paired-end RNA-seq data. if you have data come from multi librarys, input multi fastq files separated by comma. the compress file format .gz also can be accepted.

    --se <string>    default: None, Not Required.
    fastq format file contains of single-end RNA-seq data. if you have data come from multi librarys, input multi fastq files separated by comma. the compress file format .gz also can be accepted.

	--sam <string>    default: None, Not Required. --pe1/--pe2, --se, and --sam can be set together, all data will be used for analysis.
    SAM format file contains alignments of RNA-seq data. if you have data come from multi librarys, input multi SAM files separated by comma. the compress file format .bam also can be accepted.

    --protein <string>    default: None, Not Required but Recommened.
    homologous protein sequences (derived from multiple species would be recommended) file in fasta format. more protein sequences, more better.

    --augustus_species <string>    Required when --use_existed_augustus_species were not provided
    species identifier for Augustus. the relative hmm files of augustus training will be created with this prefix. if the relative hmm files of augustus training exists, the program will delete the hmm files directory firstly, and then start the augustus training steps.

    --use_existed_augustus_species <string>    Required when --augustus_species were not provided
    species identifier for Augustus. This parameter is conflict with --augustus_species. When this parameter set, the --augustus_species parameter will be invalid, and the relative hmm files of augustus training should exists, and the augustus training step will be skipped (this will save lots of runing time).

[other]
    --out_prefix <string>    default: out
    the prefix of outputs.

    --config <string>    default: None
    Input a file containing the parameters of several main programs (such as trimmomatic, hisat2 and augustus) during the pipeline. If you do not input this file, the default parameters should be suitable for most situation.

	--genetic_code <int>    default: 1
    设置遗传密码。该参数对应的值请参考NCBI Genetic Codes: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi。用于设置对基因模型强制补齐时的起始密码子和终止密码子。
    
    --RM_species <string>    default: None
    species identifier for RepeatMasker. The acceptable value of this parameter can be found in file $dirname/RepeatMasker_species.txt. Such as, Eukaryota for eucaryon, Fungi for fungi, Viridiplantae for plants, Metazoa for animals. The repeats in genome sequences would be searched aganist the Repbase database when this parameter set. 

    --RM_lib <string>    default: None
    A fasta file of repeat sequences. Generally to be the result of RepeatModeler. If not set, RepeatModeler will be used to product this file automaticly, which shall time-consuming.

    --augustus_species_start_from <string>    default: None
    species identifier for Augustus. The optimization step of Augustus training will start from the parameter file of this species, so it may save much time when setting a close species.

    --cpu <int>    default: 4
    the number of threads.

    --strand_specific    default: False
    enable the ability of analysing the strand-specific information provided by the tag "XS" from SAM format alignments. If this parameter was set, the paramter "--rna-strandness" of hisat2 should be set to "RF" usually.

    --HMM_db <string>    default: None
    the absolute path of protein family HMM database which was used for filtering of false positive gene models. multiple databases can be input, and the prefix of database files should be seperated by comma.

    --BLASTP_db <string>    default: None
    the absolute path of protein family diamond database which was used for filtering of false positive gene models. 若该参数没有设置，程序会以homologous protein构建diamond数据库，进行基因模型过滤。multiple databases can be input, and the prefix of database files should be seperated by comma.

    --gene_prefix <string>    default: gene
    the prefix of gene id shown in output file.

    --enable_augustus_training_iteration    default: False
    开启augustus_training_iteration，运行在第一次Augustus training后，根据基因预测的结果，选择有证据支持的基因模型，再一次进行Augustus training（迭代）。此举会消耗较多计算时间，且可能对基因预测没有改进，或产生不好的影响。

    --no_alternative_splicing_analysis    default: None
    添加该参数后，程序不会进行可变剪接分析。

    --delete_unimportant_intermediate_files    defaults: None
    添加该参数后，若程序运行成功，会删除不重要的中间文件，仅保留最少的、较小的、重要的中间结果文件。删除这些中间文件后，若程序重新运行后，能继续运行GETA的第六个步骤，即合并三种预测结果并对基因模型进行过滤。

    --keep_all_files_of_6thStep_combineGeneModels    defaults: None
    若使用--delete_unimportant_intermediate_files参数后，默认会删除GETA流程第6步生成的文件夹，再次运行程序会完整运行第6步流程。若再额外添加本参数，则会保留第6步所有数据，再次运行程序，则会根据第6步结果数据直接生成最终结果文件，可以用于修改基因ID前缀。

This script was tested on CentOS 8.4 with such softwares can be run directly in terminal:
01. ParaFly
02. java (version: 1.8.0_282)
03. hisat2 (version: 2.1.0)
04. samtools (version: 1.10)
05. hmmscan (version: 3.3.1)
06. makeblastdb/tblastn/blastp (version: 2.6.0)
07. RepeatMasker (version: 4.1.2-p1)
08. RepeatModeler (version: 2.0.3)
09. genewise (version: 2.4.1)
10. augustus/etraining (version: 3.4.0)
11. diamond (version 2.0.2.140)

Version: 2.5.7

USAGE
if (@ARGV==0){die $usage}

my ($RM_species, $RM_lib, $genome, $out_prefix, $pe1, $pe2, $single_end, $sam, $protein, $cpu, $trimmomatic, $strand_specific, $sam2transfrag, $ORF2bestGeneModels, $augustus_species, $HMM_db, $BLASTP_db, $gene_prefix, $cmdString, $enable_augustus_training_iteration, $config, $use_existed_augustus_species, $augustus_species_start_from, $no_alternative_splicing_analysis, $delete_unimportant_intermediate_files, $keep_all_files_of_6thStep_combineGeneModels, $genetic_code);
GetOptions(
    "RM_species:s" => \$RM_species,
    "RM_lib:s" => \$RM_lib,
    "genome:s" => \$genome,
    "out_prefix:s" => \$out_prefix,
    "pe1:s" => \$pe1,
    "pe2:s" => \$pe2,
    "S:s" => \$single_end,
	"sam:s" => \$sam,
    "protein:s" => \$protein,
    "cpu:i" => \$cpu,
    "strand_specific!" => \$strand_specific,
    "augustus_species:s" => \$augustus_species,
    "use_existed_augustus_species:s" => \$use_existed_augustus_species,
    "HMM_db:s" => \$HMM_db,
    "BLASTP_db:s" => \$BLASTP_db,
    "gene_prefix:s" => \$gene_prefix,
    "enable_augustus_training_iteration!" => \$enable_augustus_training_iteration,
    "config:s" => \$config,
    "augustus_species_start_from:s" => \$augustus_species_start_from,
    "no_alternative_splicing_analysis!" => \$no_alternative_splicing_analysis,
    "delete_unimportant_intermediate_files!" => \$delete_unimportant_intermediate_files,
    "keep_all_files_of_6thStep_combineGeneModels!" => \$keep_all_files_of_6thStep_combineGeneModels,
	"genetic_code:i" => \$genetic_code,
);

# 检测依赖的软件是否满足。
&detecting_dependent_softwares();

# 输入参数设置
# 输入基因组序列、同源蛋白序列、重复序列数据库、配置文件、转录组数据、Augustus物种信息、HMM数据库和BLASTP数据库。
$genome = abs_path($ARGV[0]);
$protein = abs_path($protein) if $protein;
$RM_lib = abs_path($RM_lib) if $RM_lib;
$config = abs_path($config) if $config;
if ( $pe1 && $pe2 ) {
	my @files; foreach ( split /,/, $pe1 ) { push @files, abs_path($_); } $pe1 = join ",", @files;
	my @files; foreach ( split /,/, $pe2 ) { push @files, abs_path($_); } $pe2 = join ",", @files;
}
if ( $single_end ) {
	my @files; foreach ( split /,/, $single_end ) { push @files, abs_path($_); } $single_end = join ",", @files;
}
if ( $sam ) {
	my @files; foreach ( split /,/, $sam ) { push @files, abs_path($_); } $sam = join ",", @files;
}
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
my (%HMM_db, %BLASTP_db);
if ( $HMM_db ) {
    foreach ( split /,/, $HMM_db ) {
        $_ = abs_path($_);
        $HMM_db{$_} = basename($_);
    }
}
if ( $BLASTP_db ) {
    foreach ( split /,/, $BLASTP_db ) {
        $BLASTP_db{$_} = basename($_);
    }
}
die "No genome sequences was found in file $genome\n" unless -s $genome;
die "No RNA-Seq short reads or homologous proteins was input\n" unless (($pe1 && $pe2) or $single_end or $sam or $protein);
die "No Augustus species provided\n" unless ($augustus_species or $use_existed_augustus_species);

# 其它参数
$genetic_code ||= 1;
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
    'homolog_prediction' => '--identity 0.2 --evalue 1e-9 --subject_coverage 0.3 --max_hits_num_per_match_region 10 --max_hit_num_per_single_species 2 --method all --genetic_code 1 --threshod_ratio_of_intron_Supported_times 0.5',
    'homolog_predictionGFF2GFF3' => '--min_score 15 --gene_prefix genewise --filterMiddleStopCodon',
    'geneModels2AugusutsTrainingInput' => '--min_evalue 1e-9 --min_identity 0.8 --min_coverage_ratio 0.8 --min_cds_num 2 --min_cds_length 450 --min_cds_exon_ratio 0.60',
    'BGM2AT' => '--min_gene_number_for_augustus_training 500 --gene_number_for_accuracy_detection 200 --min_gene_number_of_optimize_augustus_chunk 50 --max_gene_number_of_optimize_augustus_chunk 200',
    'prepareAugusutusHints' => '--margin 20',
    'paraAugusutusWithHints' => '--gene_prefix augustus --min_intron_len 20',
    'paraCombineGeneModels' => '--overlap 30 --min_augustus_transcriptSupport_percentage 10.0 --min_augustus_intronSupport_number 1 --min_augustus_intronSupport_ratio 0.01',
    'pickout_better_geneModels_from_evidence' => '--overlap_ratio 0.2 --ratio1 2 --ratio2 1.5 --ratio3 0.85 --ratio4 0.85',
    'fillingEndsOfGeneModels' => '--start_codon ATG --stop_codon TAG,TGA,TAA',
    'alternative_splicing_analysis' => '--min_intron_depth 1 --min_base_depth_ratio_for_ref_specific_intron 0.3 --min_intron_depth_ratio_for_evidence_specific_intron 0.2 --min_base_depth_ratio_for_common_intron 0.2 --min_gene_depth 10 --min_transcript_confidence_for_output 0.05 --transcript_num_for_output_when_all_low_confidence 8 --added_mRNA_ID_prefix t',
    'GFF3_extract_TranscriptID_for_filtering' => '--min_CDS_ratio 0.3 --min_CDS_length 600 --max_repeat_overlap_ratio 0.3 --ignore_repeat_Name Simple_repeat,Low_complexity,Satellite,Unknown,Tandem_repeat',
    'para_hmmscan' => '--evalue1 1e-5 --evalue2 1e-3 --hmm_length 80 --coverage 0.25 --no_cut_ga --chunk 20 --hmmscan_cpu 2',
    'diamond' => '--sensitive --max-target-seqs 20 --evalue 1e-5 --id 10 --index-chunks 1 --block-size 5',
    'parsing_blast_result.pl' => '--evalue 1e-9 --identity 0.1 --CIP 0.4 --subject-coverage 0.4 --query-coverage 0.4',
    'get_valid_geneModels' => '',
    'get_valid_transcriptID' => '--hmm_evalue 1e-7 --hmm_coverage 0.4 --blast_evalue 1e-10 --blast_CIP 0.5 --blast_coverage 0.5 --blast_evalue_for_genesie 1e-10 --blast_CIP_for_genewise 0.5 --blast_coverage_for_genewise 0.8',
    'remove_genes_in_repeats1' => '--ratio 0.3 --ignore_Simple_repeat --ignore_Unknown',
    'remove_genes_in_repeats2' => '--ratio 0.8',
    'remove_short_genes' => '--cds_length 300',
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

# 检测剩余可用内存量
my $MemAvailable = 0;
$MemAvailable = &get_MemAvailable();

# 生成临时文件夹
mkdir "$out_prefix.tmp" unless -e "$out_prefix.tmp";
chdir "$out_prefix.tmp";
my $pwd = `pwd`; print STDERR "PWD: $pwd";

# 准备基因组序列：读取FASTA序列以>开始的头部时，去除第一个空及之后的字符；去除基因组序列中尾部的换行符，将所有小写字符变换为大写字符。
unless (-e "genome.fasta") {
    open OUT, ">", "genome.fasta" or die "Can not create file genome.fasta, $!\n";
    open IN, $genome or die "Can not open file $genome, $!\n";
    $_ = <IN>;
    s/\s.*\s?$/\n/;
    print OUT;
    while (<IN>) {
        if (m/^>/) {
            s/\s.*\s?$/\n/;
            print OUT "\n$_";
        }
        else {
            s/\s+?$//g;
            $_ = uc($_);
            print OUT;
        }
    }
    print "\n";
    close IN;
    close OUT;
}
$pwd = `pwd`; chomp($pwd);
$genome = "$pwd/genome.fasta";

# 准备直系同源基因蛋白序列：读取FASTA序列以>开始的头部时，去除第一个空及之后的字符，将所有怪异字符变为下划线字符；去除序列中尾部的换行符, 将所有小写字符变换为大写字符。
unless (-e "homolog.fasta") {
    open OUT, ">", "homolog.fasta" or die "Can not create file homolog.fasta, $!\n";
    open IN, $protein or die "Can not open file $protein, $!\n";
    $_ = <IN>;
    s/\s.*\s?$/\n/; s/[^\w\n>]/_/g;
    print OUT;
    while (<IN>) {
        if (m/^>/) {
            s/\s.*\s?$/\n/; s/[^\w\n>]/_/g;
            print OUT "\n$_";
        }
        else {
            s/\s+?$//g;
        $_ = uc($_);
            print OUT;
        }
    }
    close IN;
    close OUT;
}
$pwd = `pwd`; chomp($pwd);
$protein = "$pwd/homolog.fasta";

# Step 0: RepeatMasker and RepeatModeler
print STDERR "Step 0: RepeatMasker and RepeatModeler " . "(" . (localtime) . ")" . "\n";
mkdir "0.RepeatMasker" unless -e "0.RepeatMasker";
unless (-e "0.RepeatMasker.ok") {
    chdir "0.RepeatMasker";
    mkdir "repeatMasker" unless -e "repeatMasker";
    chdir "repeatMasker";
    $pwd = `pwd`; print STDERR "PWD: $pwd";
    
    # 进行RepeatMasker分析
    if ( $RM_species ) {
        $cmdString = "para_RepeatMasker --species $RM_species --cpu $cpu --tmp_dir para_RepeatMasker.tmp $genome &> para_RepeatMasker.log";
    }
    else {
        $cmdString = "touch RepeatMasker_out.out";
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
        my $cpu_RepeatModeler = $cpu;
        $cmdString = "RepeatModeler -threads $cpu_RepeatModeler -database species -LTRStruct &> RepeatModeler.log";
        # 若RepeatModeler版本为2.0.3或更低，则需要将-threads参数换为-pa参数。
        my $RepeatModeler_info = `RepeatModeler`;
        my $RepeatModeler_version = $1 if $RepeatModeler_info =~ m/RepeatModeler - ([\d\.]+)/;
        if ( $RepeatModeler_version =~ m/(\d+)\.(\d+)\.(\d+)/ && ($1 < 2 or ( $1 == 2 && $1 == 0 && $3 <= 3)) )  {
			$cpu_RepeatModeler = int($cpu / 4);
			$cpu_RepeatModeler = 1 if $cpu_RepeatModeler < 1;
            $cmdString = "RepeatModeler -pa $cpu_RepeatModeler -database species -LTRStruct &> RepeatModeler.log";
        }
        unless (-e "RepeatModeler.ok") {
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            open OUT, ">", "RepeatModeler.ok" or die $!; close OUT;
        }
        else {
            print STDERR "CMD(Skipped): $cmdString\n";
        }
        $cmdString = "para_RepeatMasker --out_prefix RepeatModeler_out --lib RM_\*/\*.classified --cpu $cpu --tmp_dir para_RepeatMasker.tmp $genome &> para_RepeatMasker.log";
    }
    else {
        $cmdString = "para_RepeatMasker --out_prefix RepeatModeler_out --lib $RM_lib --cpu $cpu --tmp_dir para_RepeatMasker.tmp $genome &> para_RepeatMasker.log";
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
    $cmdString = "$dirname/bin/merge_repeatMasker_out.pl $genome repeatMasker/RepeatMasker_out.out repeatModeler/RepeatModeler_out.out > genome.repeat.stats";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    $cmdString = "$dirname/bin/maskedByGff.pl genome.repeat.gff3 $genome > genome.masked.fasta";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    #$cmdString = "$dirname/bin/maskedByGff.pl --mask_type softmask genome.repeat.gff3 $genome > genome.softmask.fasta";
    #print STDERR (localtime) . ": CMD: $cmdString\n";
    #system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    chdir "../";
    open OUT, ">", "0.RepeatMasker.ok" or die $!; close OUT;
}
else {
    print STDERR "Skip Step 0 for the file 0.RepeatMasker.ok exists\n";
}


# Step 1: NGSReads_predcition
print STDERR "\n============================================\n";
print STDERR "Step 1: NGSReads_predcition " . "(" . (localtime) . ")" . "\n";
my $pwd = `pwd`; print STDERR "PWD: $pwd";
mkdir "1.NGSReads_prediction" unless -e "1.NGSReads_prediction";
my $cmdString_Step1_paramter;
if (($pe1 && $pe2) or $single_end or $sam) {
	my @input_paramter;
	push @input_paramter, "--pe1 $pe1 --pe2 $pe2" if ($pe1 && $pe2);
	push @input_paramter, "--se $single_end" if $single_end;
	push @input_paramter, "--sam $sam" if $sam;
	push @input_paramter, "--config $config" if defined $config;
	push @input_paramter, "--strand_specific" if defined $strand_specific;
	push @input_paramter, "--genetic_code $genetic_code" if defined $genetic_code;
	$cmdString_Step1_paramter = join " ", @input_paramter;
	$cmdString = "$dirname/bin/NGSReads_prediction $cmdString_Step1_paramter --cpu $cpu --tmp_dir 1.NGSReads_prediction $genome > 1.NGSReads_prediction/NGSReads_prediction.gff3 2> 1.NGSReads_prediction/NGSReads_prediction.A.log";
	unless ( -e "1.NGSReads_prediction.ok" ) {
		print STDERR (localtime) . ": CMD: $cmdString\n";
		system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
		open OUT, ">", "1.NGSReads_prediction.ok" or die $!; close OUT;
	}
	else {
		 print STDERR "CMD(Skipped): $cmdString\n";
	}
}
else {
    open OUT, ">", "1.NGSReads_prediction.ok" or die "Can not create file 1.NGSReads_prediction.ok, $!"; close OUT;
}


# Step 2: homolog_prediction
print STDERR "\n============================================\n";
print STDERR "Step 2: Homolog_prediction" . "(" . (localtime) . ")" . "\n";
mkdir "2.homolog_prediction" unless -e "2.homolog_prediction";
unless (-e "2.homolog_prediction.ok") {
    chdir "2.homolog_prediction";
    $pwd = `pwd`; print STDERR "PWD: $pwd";

    open IN, "../1.NGSReads_prediction/NGSReads_prediction.gff3" or die "Can not open the file ../1.NGSReads_prediction/NGSReads_prediction.gff3, $!\n";
    my @gene_length;
    while (<IN>) {
        if (m/\tgene\t(\d+)\t(\d+)\t/) {
            push @gene_length, $2 - $1 + 1;
        }
    }
    close IN;
    @gene_length = sort {$a <=> $b} @gene_length;
    my $max_gene_length = 20000;
    $max_gene_length = $gene_length[@gene_length * 0.995] if $gene_length[@gene_length * 0.995];
    my ($segmentSize, $overlapSize) = (1000000, 100000);
    if ($max_gene_length * 4 > $overlapSize) {
        $overlapSize = $max_gene_length * 4;
        my $overlapSize_length = length($overlapSize);
        $overlapSize_length --;
        $overlapSize_length --;
        $overlapSize = int(($overlapSize / (10 ** $overlapSize_length)) + 1) * (10 ** $overlapSize_length);
        $segmentSize = $overlapSize * 10;
    }

    $cmdString = "$dirname/bin/homolog_prediction --tmp_dir homolog_prediction.tmp --cpu $cpu --segmentSize $segmentSize --overlapSize $overlapSize $config{'homolog_prediction'} $protein ../0.RepeatMasker/genome.masked.fasta > homolog_prediction.gff3 2> homolog_prediction.log";
    unless (-e "homolog_prediction.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "homolog_prediction.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    chdir "../";
    open OUT, ">", "2.homolog_prediction.ok" or die $!; close OUT;
}
else {
    print STDERR "Skip Step 2 for the file 2.homolog_prediction.ok exists\n";
}


# Step 3: Combine NGSReads and Homolog Predictions
print STDERR "\n============================================\n";
print STDERR "Step 3: Combine NGSReads and Homolog Predictions" . "(" . (localtime) . ")" . "\n";
mkdir "3.evidence_gene_models" unless -e "3.evidence_gene_models";
unless (-e "3.evidence_gene_models.ok") {
	chdir "3.evidence_gene_models";
	unlink "../1.NGSReads_prediction/c.transcript/10.FillingGeneModelsByHomolog.ok";
	unlink "../1.NGSReads_prediction/c.transcript/11.fillingEndsOfGeneModels.ok";
	unlink "../1.NGSReads_prediction/c.transcript/12.classGeneModels.ok";
	unlink "../1.NGSReads_prediction/c.transcript/FillingGeneModelsByHomolog_tmp/command.combineGeneModels.list.completed";
	$cmdString = "$dirname/bin/NGSReads_prediction $cmdString_Step1_paramter --cpu $cpu --tmp_dir ../1.NGSReads_prediction --homolog_gene_models ../2.homolog_prediction/homolog_prediction.gff3 $genome > NGSReads_prediction_FilledByHomolog.gff3 2> ../1.NGSReads_prediction/NGSReads_prediction.B.log";
	print STDERR (localtime) . ": CMD: $cmdString\n";
	system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

	open IN, "NGSReads_prediction_FilledByHomolog.gff3" or die "Can not create file NGSReads_prediction_FilledByHomolog.gff3, $!";
	open OUT1, ">", "NGSReads_prediction_FilledByHomolog.ABC.gff3" or die "Can not create file NGSReads_prediction_FilledByHomolog.ABC.gff3, $!";
	open OUT2, ">", "NGSReads_prediction_FilledByHomolog.D.gff3" or die "Can not create file NGSReads_prediction_FilledByHomolog.D.gff3, $!";
	$/ = "\n\n";
	while (<IN>) {
		if (m/Type=excellent/ or m/Type=good/ or m/Type=fair/) {
			print OUT1;
		}
		else {
			print OUT2;
		}
	}
	close IN; close OUT1; close OUT2;
	$/ = "\n";

	open IN, "../2.homolog_prediction/homolog_prediction.gff3" or die "Can not create file ../2.homolog_prediction/homolog_prediction.gff3, $!";
	open OUT1, ">", "homolog_prediction.AB.gff3" or die "Can not create file, homolog_prediction.AB.gff3 $!";
	open OUT2, ">", "homolog_prediction.CD.gff3" or die "Can not create file, homolog_prediction.CD.gff3 $!";
	$/ = "\n\n";
	while (<IN>) {
		if (m/Type=excellent/ or m/Type=good/) {
			print OUT1;
		}
		else {
			print OUT2;
		}
	}
	close IN; close OUT1; close OUT2;
	$/ = "\n";

	print STDERR (localtime) . ": Create files: NGSReads_prediction_FilledByHomolog.ABC.gff3, NGSReads_prediction_FilledByHomolog.D.gff3, homolog_prediction.AB.gff3, homolog_prediction.CD.gff3.\n";

	$cmdString = "$dirname/bin/GFF3Clear --genome $genome NGSReads_prediction_FilledByHomolog.ABC.gff3 homolog_prediction.AB.gff3 NGSReads_prediction_FilledByHomolog.D.gff3 homolog_prediction.CD.gff3 > evidence_gene_models.gff3 2> GFF3Clear.log";
	print STDERR (localtime) . ": CMD: $cmdString\n";
	system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

	open IN, "evidence_gene_models.gff3" or die "Can not open file evidence_gene_models.gff3, $!";
	my ( $num0, $num1, $num2, $num3, $num4, $num5, $num6 ) = (0, 0, 0, 0, 0, 0, 0);
	while (<IN>) {
    if (m/\tgene\t/) {
        $num0 ++;
		if ( m/ID=transfrag/ ) { $num5 ++; }
		elsif ( m/ID=homolog/ ) { $num6 ++; }
        if ( m/Type=excellent_gene_models_predicted_by_NGSReads/) { $num1 ++; }
        elsif ( m/Type=good_gene_models_predicted_by_NGSReads/) { $num2 ++; }
        elsif ( m/Type=fair_gene_models_predicted_by_NGSReads/) { $num3 ++; }
        elsif ( m/Type=poor_gene_models_predicted_by_NGSReads/) { $num4 ++; }
    }
}
close IN;
open OUT, ">", "evidence_gene_models.log" or die "Can not create file evidence_gene_models.log, $!";
print OUT "Total $num0 gene models were divided into 4 classes : A, $num1; B, $num2; C, $num3; D, $num4.\ngene models predicted by NGSReads: $num5; predicted by Homolog: $num6.\n";
print STDERR "\nFinally, total $num0 gene models were divided into 4 classes : A, $num1; B, $num2; C, $num3; D, $num4.\ngene models predicted by NGSReads: $num5; predicted by Homolog: $num6.\n";


	chdir "../";
	open OUT, ">", "3.evidence_gene_models.ok" or die $!; close OUT;
}
else {
	print STDERR "Skip Step 3 for the file 3.evidence_gene_models.ok exists\n";
}

# Step 4: Augustus gene prediction
print STDERR "\n============================================\n";
print STDERR "Step 4: Augustus/HMM Trainning " . "(" . (localtime) . ")" . "\n";
mkdir "4.augustus" unless -e "4.augustus";
unless (-e "4.augustus.ok") {
    chdir "4.augustus";
    $pwd = `pwd`; print STDERR "PWD: $pwd";

    if ($use_existed_augustus_species) {
        mkdir "training";
        open OUT, ">", "training.ok" or die $!;
        print STDERR "Skip Augustus training for --use_existed_augustus_species paramter set\n";

        chdir "training";
        $pwd = `pwd`; print STDERR "PWD: $pwd";

        # 准备Augustus training的输入文件
		$cmdString = "/bin/cp -a ../../3.evidence_gene_models/evidence_gene_models.gff3 geneModels.gff3";
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

        chdir "../";
    }

    # 第一次 Augustus HMM Training
    unless (-e "training") {
        mkdir "training";
        #my $species_config_dir = `echo \$AUGUSTUS_CONFIG_PATH`;
        #chomp($species_config_dir);
        #$species_config_dir = "$species_config_dir/species/$augustus_species";
        #$cmdString = "rm -rf $species_config_dir";
        #print STDERR "CMD: $cmdString\n";
        #(system $cmdString) == 0 or die "Failed to execute: $cmdString\n";
        unlink "training.ok" if (-e "training.ok");
    }
    unless (-e "training.ok") {
        chdir "training";
        $pwd = `pwd`; print STDERR "PWD: $pwd";

        # 准备Augustus training的输入文件
        $cmdString = "/bin/cp -a ../../3.evidence_gene_models/evidence_gene_models.gff3 geneModels.gff3";
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

        $cmdString = "$dirname/bin/geneModels2AugusutsTrainingInput $config{'geneModels2AugusutsTrainingInput'} --out_prefix ati --cpu $cpu geneModels.gff3 $genome &> geneModels2AugusutsTrainingInput.log";
        unless (-e "geneModels2AugusutsTrainingInput.ok") {
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

            # 若用于Augustus training的基因数量少于1000个，则重新运行geneModels2AugusutsTrainingInput，降低阈值来增加基因数量。
            open IN, "geneModels2AugusutsTrainingInput.log" or die "Can not open file geneModels2AugusutsTrainingInput.log, $!";
            my $training_genes_number = 0;
            while (<IN>) {
                $training_genes_number = $1 if m/Best gene Models number:\s+(\d+)/;
            }
            if ( $training_genes_number < 1000 ) {
                $cmdString = "$dirname/bin/geneModels2AugusutsTrainingInput --min_evalue 1e-9 --min_identity 0.9 --min_coverage_ratio 0.9 --min_cds_num 1 --min_cds_length 450 --min_cds_exon_ratio 0.40 --keep_ratio_for_excluding_too_long_gene 0.99 --out_prefix ati --cpu $cpu geneModels.gff3 $genome &> geneModels2AugusutsTrainingInput.log.Loose_thresholds";
                print STDERR (localtime) . ": CMD: $cmdString\n";
                system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            }

            open OUT, ">", "geneModels2AugusutsTrainingInput.ok" or die $!; close OUT;
        }
        else {
            print STDERR "CMD(Skipped): $cmdString\n";
        }

        # 进行Augustus training
        my (%gene_info, @flanking_length, $flanking_length, @gene_length);
        open IN, "geneModels.gff3" or die $!;
        while (<IN>) {
            if (m/\tgene\t/) {
                @_ = split /\t/;
                $gene_info{$_[0]}{$_[6]}{"$_[3]\t$_[4]"} = 1;
            }
        }
        close IN;
        foreach my $chr (keys %gene_info) {
            foreach my $strand (keys %{$gene_info{$chr}}) {
                my @region = sort {$a <=> $b} keys %{$gene_info{$chr}{$strand}};
                my $first_region = shift @region;
                my ($start, $end) = split /\t/, $first_region;
                push @gene_length, $end - $start + 1;
                foreach (@region) {
                    my ($aa, $bb) = split /\t/;
                    push @gene_length, $bb - $aa + 1;
                    my $distance = $aa - $end - 1;
                    push @flanking_length, $distance if $distance >= 50;
                    $end = $bb if $end < $bb;
                }
            }
        }
        @gene_length = sort {$a <=> $b} @gene_length;
        @flanking_length = sort {$a <=> $b} @flanking_length;
        $flanking_length = int($flanking_length[@flanking_length/2] / 8);
        $flanking_length = $gene_length[@gene_length/2] if $flanking_length >= $gene_length[@gene_length/2];
        $cmdString = "$dirname/bin/BGM2AT $config{'BGM2AT'} --augustus_species_start_from $augustus_species_start_from --flanking_length $flanking_length --CPU $cpu --onlytrain_GFF3 ati.filter1.gff3 ati.filter2.gff3 $genome $augustus_species &> BGM2AT.log";
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

        chdir "../";
        open OUT, ">", "training.ok" or die $!;
    }
    else {
        print STDERR "Skip Augustus training for file training.ok exists\n" unless $use_existed_augustus_species;
    }

    # Augustus Hint Preparing
    ########## version 2.5.1 ##########
    #$cmdString = "bam2hints --source=W --intronsonly --in=../2.hisat2/hisat2.sorted.bam --out=bam2intronHints.gff";
    #unless (($pe1 && $pe2) or $single_end) {
    #    open OUT, ">", "bam2hints.ok" or die $!; close OUT;
    #    open OUT, ">", "bam2intronHints.gff" or die $!; close OUT;
    #}
    #unless (-e "bam2hints.ok") {
    #    print STDERR (localtime) . ": CMD: $cmdString\n";
    #    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    #    open OUT, ">", "bam2hints.ok" or die $!; close OUT;
    #}
    #else {
    #    print STDERR "CMD(Skipped): $cmdString\n";
    #}
    #$cmdString = "$dirname/bin/prepareAugusutusHints $config{'prepareAugusutusHints'} bam2intronHints.gff ../1.NGSReads_prediction/c.transcript/transfrag.genome.gff3 ../2.homolog/genewise.gff3 ../2.homolog/genewise.start_stop_hints.gff > hints.gff";
    ########## version 2.5.1 ##########
    $cmdString = "$dirname/bin/prepareAugusutusHints $config{'prepareAugusutusHints'} ../1.NGSReads_prediction/c.transcript/intron.txt ../1.NGSReads_prediction/c.transcript/transfrag.genome.gff3 ../2.homolog_prediction/homolog_prediction.gff3 > hints.gff";
    unless (-e "prepareAugusutusHints.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "prepareAugusutusHints.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    # Get the longest gene length
    open IN, "training/geneModels.gff3" or die "Can not open the file training/geneModels.gff3, $!\n";
    my @gene_length;
    while (<IN>) {
        if (m/\tgene\t(\d+)\t(\d+)\t/) {
            push @gene_length, $2 - $1 + 1;
        }
    }
    @gene_length = sort {$a <=> $b} @gene_length;

    my ($segmentSize, $overlapSize) = (5000000, 100000);
    if ($gene_length[-1] * 4 > $overlapSize) {
        $overlapSize = $gene_length[-1] * 4;
        my $overlapSize_length = length($overlapSize);
        $overlapSize_length --;
        $overlapSize_length --;
        $overlapSize = int(($overlapSize / (10 ** $overlapSize_length)) + 1) * (10 ** $overlapSize_length);
        $segmentSize = $overlapSize * 50;
    }
    # 第一次 Augustus gene prediction
    $cmdString = "$dirname/bin/paraAugusutusWithHints $config{'paraAugusutusWithHints'} --species $augustus_species --cpu $cpu --segmentSize $segmentSize --overlapSize $overlapSize --tmp_dir aug_para_with_hints.tmp1  ../0.RepeatMasker/genome.masked.fasta hints.gff > augustus.1.gff3";
    unless (-e "first_augustus.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        unlink "command.augustus.list.completed"  if -e "command.augustus.list.completed";
        open OUT, ">", "first_augustus.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    # augustus_training_iteration
    if ($enable_augustus_training_iteration && (! $use_existed_augustus_species)) {
        unless (-e "training_again") {
            mkdir "training_again";
            my $species_config_dir = `echo \$AUGUSTUS_CONFIG_PATH`;
            chomp($species_config_dir);
            $species_config_dir = "$species_config_dir/species/$augustus_species";
            $cmdString = "rm -rf $species_config_dir && cp -a ./training/hmm_files_bak/ $species_config_dir";
            print STDERR "CMD: $cmdString\n";
            (system $cmdString) == 0 or die "Failed to execute: $cmdString\n";
            unlink "training_again.ok" if (-e "training_again.ok");
        }
        unless (-e "training_again.ok") {
            chdir "training_again";
            $pwd = `pwd`; print STDERR "PWD: $pwd";

            # 准备Augustus training的输入文件，注意需要去除Augustus结果中的可变剪接。
            open OUT, ">", "geneModels.gff3" or die "Cannot create file geneModels.gff3, $!\n";
            open IN, "../augustus.1.gff3" or die "Can not open file ../augustus.1.gff3, $!\n";
            my ($keep, $keep_num, $total_num) = (0, 0, 0);
            my (%augustus_gene_model_gene, %augustus_gene_model_mRNA, %augustus_mRNA_score, $augustus_gene_id, $augustus_mRNA_id, @keep);
            while (<IN>) {
                next if m/^\s*$/;
                next if m/^#/;
                if (m/\tgene\t/) {
                    $total_num ++;
                    @_ = split /\s+/;
                    $augustus_gene_id = $1 if $_[8] =~ /ID=([^;\s]+)/;
                    $augustus_gene_model_gene{$augustus_gene_id} = $_;
                    my %attr = $_[8] =~ m/([^;=]+)=([^;=]+)/g;
                    if ($attr{"exonHintRatio"} >= 95) {
                        $keep = 1;
                        $keep_num ++;
                    }
                    elsif ($attr{"intronSupport"} =~ m/(\d+)\/(\d+)/ && $1 == $2 && $1 > 0) {
                        $keep = 1;
                        $keep_num ++;
                    }
                    else {
                        $keep = 0;
                    }
                    push @keep, $augustus_gene_id if $keep == 1;
                }
                elsif (m/\tmRNA\t/) {
                    @_ = split /\s+/;
                    $augustus_mRNA_id = $1 if $_[8] =~ /ID=([^;\s]+)/;
                    $augustus_gene_model_mRNA{$augustus_gene_id}{$augustus_mRNA_id} .= $_;
                    $augustus_mRNA_score{$augustus_gene_id}{$augustus_mRNA_id} = $_[5];
                }
                else {
                    $augustus_gene_model_mRNA{$augustus_gene_id}{$augustus_mRNA_id} .= $_;
                }
            }
            close IN;
            foreach my $gene_id (@keep) {
                print OUT $augustus_gene_model_gene{$gene_id};
                my @mRNA_id = sort {$augustus_mRNA_score{$gene_id}{$b} <=> $augustus_mRNA_score{$gene_id}{$a}} keys %{$augustus_gene_model_mRNA{$gene_id}};
                print OUT $augustus_gene_model_mRNA{$gene_id}{$mRNA_id[0]};
                print OUT "\n";
            }
            close OUT;
            print STDERR "Total genes predicted by Augustus: $total_num\nGood genes picked for next traning: $keep_num\n";

            $cmdString = "$dirname/bin/geneModels2AugusutsTrainingInput --out_prefix ati --cpu $cpu $config{'geneModels2AugusutsTrainingInput'} geneModels.gff3 $genome &> geneModels2AugusutsTrainingInput.log";
            unless (-e "geneModels2AugusutsTrainingInput.ok") {
                print STDERR (localtime) . ": CMD: $cmdString\n";
                system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

                # 若用于Augustus training的基因数量少于1000个，则重新运行geneModels2AugusutsTrainingInput，降低阈值来增加基因数量。
                open IN, "geneModels2AugusutsTrainingInput.log" or die "Can not open file geneModels2AugusutsTrainingInput.log, $!";
                my $training_genes_number = 0;
                while (<IN>) {
                    $training_genes_number = $1 if m/Best gene Models number:\s+(\d+)/;
                }
                if ( $training_genes_number < 1000 ) {
                    $cmdString = "$dirname/bin/geneModels2AugusutsTrainingInput --min_evalue 1e-9 --min_identity 0.9 --min_coverage_ratio 0.9 --min_cds_num 1 --min_cds_length 450 --min_cds_exon_ratio 0.40 --keep_ratio_for_excluding_too_long_gene 0.99 --out_prefix ati --cpu $cpu geneModels.gff3 $genome &> geneModels2AugusutsTrainingInput.log.Loose_thresholds";
                    print STDERR (localtime) . ": CMD: $cmdString\n";
                    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
                }

                open OUT, ">", "geneModels2AugusutsTrainingInput.ok" or die $!; close OUT;
            }
            else {
                print STDERR "CMD(Skipped): $cmdString\n";
            }

            # 进行Augustus training
            my (%gene_info, @flanking_length, $flanking_length, @gene_length);
            open IN, "geneModels.gff3" or die $!;
            while (<IN>) {
                if (m/\tgene\t/) {
                    @_ = split /\t/;
                    $gene_info{$_[0]}{$_[6]}{"$_[3]\t$_[4]"} = 1;
                }
            }
            close IN;
            foreach my $chr (keys %gene_info) {
                foreach my $strand (keys %{$gene_info{$chr}}) {
                    my @region = sort {$a <=> $b} keys %{$gene_info{$chr}{$strand}};
                    my $first_region = shift @region;
                    my ($start, $end) = split /\t/, $first_region;
                    push @gene_length, $end - $start + 1;
                    foreach (@region) {
                        my ($aa, $bb) = split /\t/;
                        push @gene_length, $bb - $aa + 1;
                        my $distance = $aa - $end - 1;
                        push @flanking_length, $distance if $distance >= 50;
                        ($start, $end) = ($aa, $bb);
                    }
                }
            }
            @gene_length = sort {$a <=> $b} @gene_length;
            @flanking_length = sort {$a <=> $b} @flanking_length;
            $flanking_length = int($flanking_length[@flanking_length/2] / 8);
            $flanking_length = $gene_length[@gene_length/2] if $flanking_length >= $gene_length[@gene_length/2];
            $cmdString = "$dirname/bin/BGM2AT $config{'BGM2AT'} --flanking_length $flanking_length --CPU $cpu --onlytrain_GFF3 ati.filter1.gff3 --stopAfterFirstEtraining ati.filter2.gff3 $genome $augustus_species &> BGM2AT.log";
            unless (-e "firsttest.out") {
                print STDERR (localtime) . ": CMD: $cmdString\n";
                system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            }
            else {
                print STDERR "CMD(Skipped): $cmdString\n";
            }

            my $accuracy_value = 0;
            open IN, "../training/secondtest.out" or die "Can not open file ../training/secondtest.out, $!\n";
            while (<IN>) {
                if (m/^nucleotide level/) {
                    @_ = split /[\s\|]+/;
                    $accuracy_value += ($_[-2] * 3 + $_[-1] * 2);
                }
                elsif (m/^exon level/) {
                    @_ = split /[\s\|]+/;
                    $accuracy_value += ($_[-2] * 4 + $_[-1] * 3);
                }
                elsif (m/^gene level/) {
                    @_ = split /[\s\|]+/;
                    $accuracy_value += ($_[-2] * 2 + $_[-1] * 1);
                }
            }
            close IN;
            my $first_accuracy = $accuracy_value / 15;

            my $accuracy_value = 0;
            open IN, "firsttest.out" or die "Can not open file firsttest.out, $!\n";
            while (<IN>) {
                if (m/^nucleotide level/) {
                    @_ = split /[\s\|]+/;
                    $accuracy_value += ($_[-2] * 3 + $_[-1] * 2);
                }
                elsif (m/^exon level/) {
                    @_ = split /[\s\|]+/;
                    $accuracy_value += ($_[-2] * 4 + $_[-1] * 3);
                }
                elsif (m/^gene level/) {
                    @_ = split /[\s\|]+/;
                    $accuracy_value += ($_[-2] * 2 + $_[-1] * 1);
                }
            }
            close IN;
            my $second_accuracy = $accuracy_value / 15;
            print STDERR "The accuracy value of augustus training is: $first_accuracy\n";
            print STDERR "The accuracy value of augustus training iteration is: $second_accuracy\n";

            if ($second_accuracy > $first_accuracy) {
                $cmdString = "$dirname/bin/BGM2AT $config{'BGM2AT'} --flanking_length $flanking_length --CPU $cpu --onlytrain_GFF3 ati.filter1.gff3 ati.filter2.gff3 $genome $augustus_species &>> BGM2AT.log";
                print STDERR (localtime) . ": CMD: $cmdString\n";
                system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

                chdir "../";
                $pwd = `pwd`; print STDERR "PWD: $pwd";
                open OUT, ">", "training_again.ok" or die $!;

                # Get the longest gene length 
                open IN, "augustus.1.gff3" or die "Can not open the file augustus.1.gff3, $!\n";
                my @gene_length;
                while (<IN>) {
                    if (m/\tgene\t(\d+)\t(\d+)\t/) { 
                        push @gene_length, $2 - $1 + 1;
                    }
                }
                @gene_length = sort {$a <=> $b} @gene_length;

                my ($segmentSize, $overlapSize) = (5000000, 100000);
                if ($gene_length[-1] * 4 > $overlapSize) {
                    $overlapSize = $gene_length[-1] * 4;
                    my $overlapSize_length = length($overlapSize);
                    $overlapSize_length --;
                    $overlapSize_length --;
                    $overlapSize = int(($overlapSize / (10 ** $overlapSize_length)) + 1) * (10 ** $overlapSize_length);
                    $segmentSize = $overlapSize * 50;
                }

                # Augustus gene prediction
                $cmdString = "$dirname/bin/paraAugusutusWithHints $config{'paraAugusutusWithHints'} --species $augustus_species --cpu $cpu --segmentSize $segmentSize --overlapSize $overlapSize --tmp_dir aug_para_with_hints.tmp2 ../0.RepeatMasker/genome.masked.fasta hints.gff > augustus.2.gff3";
                print STDERR (localtime) . ": CMD: $cmdString\n";
                system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

                $cmdString = "$dirname/bin/addHintRatioToAugustusResult training/geneModels.gff3 hints.gff augustus.2.gff3 > augustus.gff3";
                print STDERR (localtime) . ": CMD: $cmdString\n";
                system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            }
            else {
                print STDERR "The iteration step could not increase the accuracy of augustus training! and the hmm files will be rolled back!\n";
                chdir "../";
                $pwd = `pwd`; print STDERR "PWD: $pwd";
                my $species_config_dir = `echo \$AUGUSTUS_CONFIG_PATH`;
                chomp($species_config_dir);
                $species_config_dir = "$species_config_dir/species/$augustus_species";
                $cmdString = "rm -rf $species_config_dir && cp -a ./training/hmm_files_bak/ $species_config_dir";
                print STDERR "CMD: $cmdString\n";
                (system $cmdString) == 0 or die "Failed to execute: $cmdString\n";

                open OUT, ">", "training_again.ok" or die $!;

                $cmdString = "$dirname/bin/addHintRatioToAugustusResult training/geneModels.gff3 hints.gff augustus.1.gff3 > augustus.gff3";
                print STDERR (localtime) . ": CMD: $cmdString\n";
                system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            }
        }
        else {
            print STDERR "Skip Augustus training again for file training_again.ok exists\n";

            unless (-e "augustus.gff3") {
                # Get the longest gene length
                open IN, "augustus.1.gff3" or die "Can not open the file augustus.1.gff3, $!\n";
                my @gene_length;
                while (<IN>) {
                    if (m/\tgene\t(\d+)\t(\d+)\t/) { 
                        push @gene_length, $2 - $1 + 1;
                    }
                }
                @gene_length = sort {$a <=> $b} @gene_length;

                my ($segmentSize, $overlapSize) = (5000000, 100000);
                if ($gene_length[-1] * 4 > $overlapSize) {
                    $overlapSize = $gene_length[-1] * 4;
                    my $overlapSize_length = length($overlapSize);
                    $overlapSize_length --;
                    $overlapSize_length --;
                    $overlapSize = int(($overlapSize / (10 ** $overlapSize_length)) + 1) * (10 ** $overlapSize_length);
                    $segmentSize = $overlapSize * 50;
                }

                # Augustus gene prediction
                $cmdString = "$dirname/bin/paraAugusutusWithHints $config{'paraAugusutusWithHints'} --species $augustus_species --cpu $cpu --segmentSize $segmentSize --overlapSize $overlapSize --tmp_dir aug_para_with_hints.tmp2 ../0.RepeatMasker/genome.masked.fasta hints.gff > augustus.2.gff3";
                print STDERR (localtime) . ": CMD: $cmdString\n";
                system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

                $cmdString = "$dirname/bin/addHintRatioToAugustusResult training/geneModels.gff3 hints.gff augustus.2.gff3 > augustus.gff3";
                print STDERR (localtime) . ": CMD: $cmdString\n";
                system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            }
        }
    }
    else {
        $cmdString = "$dirname/bin/addHintRatioToAugustusResult training/geneModels.gff3 hints.gff augustus.1.gff3 > augustus.gff3";
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    }

    chdir "../";
    open OUT, ">", "4.augustus.ok" or die $!; close OUT;
}
else {
    print STDERR "Skip Step 4 for the file 4.augustus.ok exists\n";
}


# Step 6: CombineGeneModels
print STDERR "\n============================================\n";
print STDERR "Step 6: CombineGeneModels " . "(" . (localtime) . ")" . "\n";
mkdir "6.combineGeneModels" unless -e "6.combineGeneModels";
unless (-e "6.combineGeneModels.ok") {
    chdir "6.combineGeneModels";
    $pwd = `pwd`; print STDERR "PWD: $pwd";

    open OUT, ">", "geneModels.Readme" or die "Can not create file geneModels.Readme, $!";
    print OUT "geneModels.a.gff3\t第一轮整合后获得的以AUGUSTUS结果为主的有足够证据支持的基因模型。
geneModels.b.gff3\t第一轮整合后获得的由AUGUSTUS预测的没有足够证据支持的基因模型。
geneModels.c.gff3\t以Transcript和Homolog预测结果整合得到的完全由证据支持的基因模型。
geneModels.d.gff3\t第二轮整合后以Transcript和Homolog预测的更优结果为主，优化后的基因模型。
geneModels.e.gff3\t对geneModels.d.gff3中不完整基因模型成功进行强制补齐后获得的完整基因模型。
geneModels.f.gff3\t对geneModels.d.gff3中不完整基因模型未能强制补齐后获得的不完整基因模型。
geneModels.gb_AS.gff3\t对geneModels.b.gff3基因模型进行可变剪接分析的结果，增加的转录本没有CDS信息。
geneModels.ge_AS.gff3\t对geneModels.e.gff3基因模型进行可变剪接分析的结果，增加的转录本没有CDS信息。
geneModels.gf_AS.gff3\t对geneModels.f.gff3基因模型进行可变剪接分析的结果，增加的转录本没有CDS信息。
geneModels.gb.gff3\t对geneModels.b.gff3基因模型进行可变剪接分析的结果，增加的转录本增加了CDS信息。
geneModels.ge.gff3\t对geneModels.e.gff3基因模型进行可变剪接分析的结果，增加的转录本增加了CDS信息。
geneModels.gf.gff3\t对geneModels.f.gff3基因模型进行可变剪接分析的结果，增加的转录本增加了CDS信息。
geneModels.h.coding.gff3\t对基因模型进行过滤，获得的编码蛋白基因模型。
geneModels.h.lncRNA.gff3\t对基因模型进行过滤，获得的lnc_RNA基因模型。
geneModels.h.lowQuality.gff3\t对基因模型进行过滤，获得的低质量基因模型。
geneModels.i.coding.gff3\t对geneModels.h.coding.gff3中的基因模型进行了强制补齐\n";
    close OUT;

    # 6.1 第一轮基因预测结果整合：以AUGUSTUS结果为主，进行三种基因预测结果的整合
    # 对三种基因预测结果进行第一轮整合，以Augustus结果为准。得到 combine.1.gff3 为有Evidence支持的结果，combine.2.gff3为支持不足的结果。
    my $cmdString1 = "$dirname/bin/paraCombineGeneModels $config{'paraCombineGeneModels'} --cpu $cpu ../4.augustus/augustus.gff3 ../1.NGSReads_prediction/c.transcript/transfrag.genome.gff3 ../2.homolog_prediction/homolog_prediction.gff3 ../4.augustus/hints.gff &> /dev/null";
    my $cmdString2 = "$dirname/bin/GFF3Clear --genome $genome --no_attr_add --coverage 0.8 combine.1.gff3 > geneModels.a.gff3 2> /dev/null";
    my $cmdString3 = "$dirname/bin/GFF3Clear --genome $genome --no_attr_add --coverage 0.8 combine.2.gff3 > geneModels.b.gff3 2> /dev/null";
    my $cmdString4 = "perl -p -i -e 's/(=[^;]+)\.t1/\$1.t01/g' geneModels.a.gff3 geneModels.b.gff3";
    unless (-e "01.paraCombineGeneModels.ok") {
        print STDERR (localtime) . ": CMD: $cmdString1\n";
        system("$cmdString1") == 0 or die "failed to execute: $cmdString1\n";
        print STDERR (localtime) . ": CMD: $cmdString2\n";
        system("$cmdString2") == 0 or die "failed to execute: $cmdString2\n";
        print STDERR (localtime) . ": CMD: $cmdString3\n";
        system("$cmdString3") == 0 or die "failed to execute: $cmdString3\n";
        print STDERR (localtime) . ": CMD: $cmdString4\n";
        system("$cmdString4") == 0 or die "failed to execute: $cmdString4\n";
        open OUT, ">", "01.paraCombineGeneModels.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString1\n";
        print STDERR "CMD(Skipped): $cmdString2\n";
        print STDERR "CMD(Skipped): $cmdString3\n";
        print STDERR "CMD(Skipped): $cmdString4\n";
    }

    # 6.2 第二轮基因预测结果整合：以转录本和同源蛋白预测结果为准，对上一步的基因模型进行优化。
    my $cmdString1 = "perl -p -e 's/(=[^;]+)\.t1/\$1.t01/g;' ../4.augustus/training/geneModels.gff3 > geneModels.c.gff3;";
    my $cmdString2 = "$dirname/bin/pickout_better_geneModels_from_evidence $config{'pickout_better_geneModels_from_evidence'} geneModels.a.gff3 geneModels.c.gff3 > picked_evidence_geneModels.gff3 2> picked_evidence_geneModels.log";
    my $cmdString3 = "$dirname/bin/GFF3Clear --genome $genome --no_attr_add picked_evidence_geneModels.gff3 geneModels.a.gff3 > geneModels.d.gff3 2> /dev/null";
    my $cmdString4 = "perl -p -i -e 's/Integrity=[^;]+;?//g' geneModels.d.gff3";
    unless (-e "02.pickout_better_geneModels_from_evidence.ok") {
        # 先挑选出更优的有Evidence支持的基因模型
        print STDERR (localtime) . ": CMD: $cmdString1\n";
        system("$cmdString1") == 0 or die "failed to execute: $cmdString1\n";
        # 再用更优的有Evidence支持的基因模型替换掉相应的基因模型
        print STDERR (localtime) . ": CMD: $cmdString2\n";
        system("$cmdString2") == 0 or die "failed to execute: $cmdString2\n";
        print STDERR (localtime) . ": CMD: $cmdString3\n";
        system("$cmdString3") == 0 or die "failed to execute: $cmdString3\n";
        print STDERR (localtime) . ": CMD: $cmdString4\n";
        system("$cmdString4") == 0 or die "failed to execute: $cmdString4\n";
        open OUT, ">", "02.pickout_better_geneModels_from_evidence.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString1\n";
        print STDERR "CMD(Skipped): $cmdString2\n";
        print STDERR "CMD(Skipped): $cmdString3\n";
        print STDERR "CMD(Skipped): $cmdString4\n";
    }

    # 6.3 对不完整基因模型进行首尾补齐。
    $cmdString = "$dirname/bin/fillingEndsOfGeneModels $config{'fillingEndsOfGeneModels'} --filling_need_transcriptID filling_need_transcriptID.txt --nonCompletedGeneModels geneModels.f.gff3 $genome geneModels.d.gff3 > geneModels.e.gff3 2> fillingEndsOfGeneModels.1.log";
    unless ( -e "03.fillingEndsOfGeneModels.ok" ) {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "03.fillingEndsOfGeneModels.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    # 6.4 分别对对基因模型 geneModels.b.gff3, geneModels.e.gff3 and geneModels.f.gff3 进行可变剪接分析
    unless ( $no_alternative_splicing_analysis ) {
        $cmdString1 = "$dirname/bin/paraAlternative_splicing_analysis $config{'alternative_splicing_analysis'} --tmp_dir paraAlternative_splicing_analysis.gb.tmp --cpu $cpu geneModels.b.gff3 ../1.NGSReads_prediction/c.transcript/intron.txt ../1.NGSReads_prediction/c.transcript/base_depth.txt > geneModels.gb_AS.gff3 2> alternative_splicing.gb.stats && $dirname/bin/GFF3_add_CDS_for_transcript $genome geneModels.gb_AS.gff3 > geneModels.gb.gff3";
        $cmdString2 = "$dirname/bin/paraAlternative_splicing_analysis $config{'alternative_splicing_analysis'} --tmp_dir paraAlternative_splicing_analysis.ge.tmp --cpu $cpu geneModels.e.gff3 ../1.NGSReads_prediction/c.transcript/intron.txt ../1.NGSReads_prediction/c.transcript/base_depth.txt > geneModels.ge_AS.gff3 2> alternative_splicing.ge.stats && $dirname/bin/GFF3_add_CDS_for_transcript $genome geneModels.ge_AS.gff3 > geneModels.ge.gff3";
        $cmdString3 = "$dirname/bin/paraAlternative_splicing_analysis $config{'alternative_splicing_analysis'} --tmp_dir paraAlternative_splicing_analysis.gf.tmp --cpu $cpu geneModels.f.gff3 ../1.NGSReads_prediction/c.transcript/intron.txt ../1.NGSReads_prediction/c.transcript/base_depth.txt > geneModels.gf_AS.gff3 2> alternative_splicing.gf.stats && $dirname/bin/GFF3_add_CDS_for_transcript $genome geneModels.gf_AS.gff3 > geneModels.gf.gff3";
        unless ( -e "04.alternative_splicing_analysis.ok" ) {
            print STDERR (localtime) . ": CMD: $cmdString1\n";
            system("$cmdString1") == 0 or die "failed to execute: $cmdString1\n";
            print STDERR (localtime) . ": CMD: $cmdString2\n";
            system("$cmdString2") == 0 or die "failed to execute: $cmdString2\n";
            print STDERR (localtime) . ": CMD: $cmdString3\n";
            system("$cmdString3") == 0 or die "failed to execute: $cmdString3\n";
            open OUT, ">", "04.alternative_splicing_analysis.ok" or die $!; close OUT;
        }
        else {
            print STDERR "CMD(Skipped): $cmdString1\n";
            print STDERR "CMD(Skipped): $cmdString2\n";
            print STDERR "CMD(Skipped): $cmdString3\n";
        }
    }
    else {
        $cmdString1 = "$dirname/bin/GFF3_add_CDS_for_transcript $genome geneModels.b.gff3 > geneModels.gb.gff3";
        $cmdString2 = "$dirname/bin/GFF3_add_CDS_for_transcript $genome geneModels.e.gff3 > geneModels.ge.gff3";
        $cmdString3 = "$dirname/bin/GFF3_add_CDS_for_transcript $genome geneModels.f.gff3 > geneModels.gf.gff3";
        unless ( -e "04.alternative_splicing_analysis.ok" ) {
            print STDERR (localtime) . ": CMD: $cmdString1\n";
            system("$cmdString1") == 0 or die "failed to execute: $cmdString1\n";
            print STDERR (localtime) . ": CMD: $cmdString2\n";
            system("$cmdString2") == 0 or die "failed to execute: $cmdString2\n";
            print STDERR (localtime) . ": CMD: $cmdString3\n";
            system("$cmdString3") == 0 or die "failed to execute: $cmdString3\n";
            open OUT, ">", "04.alternative_splicing_analysis.ok" or die $!; close OUT;
        }
        else {
            print STDERR "CMD(Skipped): $cmdString1\n";
            print STDERR "CMD(Skipped): $cmdString2\n";
            print STDERR "CMD(Skipped): $cmdString3\n";
        }
    }

    # 6.5 提取待过滤的转录本ID，对 geneModels.b.gff3, geneModels.e.gff3 and geneModels.f.gff3 中的基因模型进行过滤。
    #########################################################################
    # 以下几种类型的转录本用于过滤：
    # 1. CDS占转录本长度比例较低（< 30%）的转录本；
    # 2. 所有转录本的CDS长度较短（< 600 bp）的基因，对应的所有转录本；
    # 3. 所有转录本的CDS和重复序列区域重叠比例都较高（默认 >= 30%）的基因，对应的所有转录本；
    # 4. 没有足够证据支持的基因，对应的所有转录本；
    # 5. 没法填补完整的基因，对应的所有转录本；
    # 6. 通过填补而完整的基因，对应的所有转录本。
    #########################################################################
    # 提取CDS占转录本长度比例较低、CDS长度较短和CDS和重复序列区域重叠比例较高的转录本ID
    my $cmdString1 = "$dirname/bin/GFF3_extract_TranscriptID_for_filtering $config{'GFF3_extract_TranscriptID_for_filtering'} ../0.RepeatMasker/genome.repeat.gff3 geneModels.gb.gff3 geneModels.ge.gff3 geneModels.gf.gff3 > transcriptID_for_filtering.txt";
    # 提取没有足够证据基因的所有转录本ID
    my $cmdString2 = "perl -ne 'print \"\$1\\tNotEnoughEvidence\\n\" if m/ID=([^;]*\\.t\\d+);/;' geneModels.gb.gff3 >> transcriptID_for_filtering.txt";
    # 提取没法填补完整基因的所有转录本ID
    my $cmdString3 = "perl -ne 'print \"\$1\\tFilling2Uncomplete\\n\" if m/ID=([^;]*\\.t\\d+);/;' geneModels.gf.gff3 >> transcriptID_for_filtering.txt";
    # 提取通过填补而完整的基因的所有转录本ID
    my $cmdString4 = "perl -e 'open IN, \"filling_need_transcriptID.txt\"; while (<IN>) { s/.t\\d+\\n//; \$gene{\$_} = 1; } while (<>) { print \"\$1\\tFilling2Complete\\n\" if m/ID=(([^;]+)\\.t\\d+);/ && exists \$gene{\$2} }' geneModels.ge_AS.gff3 >> transcriptID_for_filtering.txt";
    unless ( -e "05.extract_TranscriptID_for_filtering.ok" ) {
        print STDERR (localtime) . ": CMD: $cmdString1\n";
        system("$cmdString1") == 0 or die "failed to execute: $cmdString1\n";
        print STDERR (localtime) . ": CMD: $cmdString2\n";
        system("$cmdString2") == 0 or die "failed to execute: $cmdString2\n";
        print STDERR (localtime) . ": CMD: $cmdString3\n";
        system("$cmdString3") == 0 or die "failed to execute: $cmdString3\n";
        print STDERR (localtime) . ": CMD: $cmdString4\n";
        system("$cmdString4") == 0 or die "failed to execute: $cmdString4\n";
        open OUT, ">", "05.extract_TranscriptID_for_filtering.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString1\n";
        print STDERR "CMD(Skipped): $cmdString2\n";
        print STDERR "CMD(Skipped): $cmdString3\n";
        print STDERR "CMD(Skipped): $cmdString4\n";
    }
    
    # 6.6 提取待过滤转录本的蛋白序列。
    my $cmdString1 = "$dirname/bin/gff3_to_protein.pl $genome geneModels.gb.gff3 geneModels.gf.gff3 geneModels.ge.gff3 > proteins_all.fasta 2> gff3_to_protein.log";
    my $cmdString2 = "perl -p -i -e 's/\\*\$//' proteins_all.fasta";
    my $cmdString3 = "$dirname/bin/fasta_extract_subseqs_from_list.pl proteins_all.fasta transcriptID_for_filtering.txt > proteins_for_filtering.fasta 2> fasta_extract_subseqs_from_list.log";
    unless ( -e "06.extract_proteins_for_filtering.ok" ) {
        print STDERR (localtime) . ": CMD: $cmdString1\n";
        system("$cmdString1") == 0 or die "failed to execute: $cmdString1\n";
        print STDERR (localtime) . ": CMD: $cmdString2\n";
        system("$cmdString2") == 0 or die "failed to execute: $cmdString2\n";
        print STDERR (localtime) . ": CMD: $cmdString3\n";
        system("$cmdString3") == 0 or die "failed to execute: $cmdString3\n";
        open OUT, ">", "06.extract_proteins_for_filtering.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString1\n";
        print STDERR "CMD(Skipped): $cmdString2\n";
        print STDERR "CMD(Skipped): $cmdString3\n";
    }
    
    # 6.7 对蛋白序列进行HMM和BLASTP验证。
    my ($cmdString1, $cmdString2, $cmdString3);
    if ( $HMM_db ) {
        my $hmmscan_cpu = 0;
        $hmmscan_cpu = $1 if $config{'para_hmmscan'} =~ m/--hmmscan_cpu\s+(\d+)/;
        my $para_hmmscan_cpu = $cpu;
        $para_hmmscan_cpu = int($cpu / $hmmscan_cpu + 0.5) if $hmmscan_cpu;
        $para_hmmscan_cpu = 1 if $para_hmmscan_cpu < 1;
        foreach ( sort keys %HMM_db ) {
            $cmdString1 .= "$dirname/bin/para_hmmscan $config{'para_hmmscan'} --outformat --cpu $para_hmmscan_cpu --no_cut_ga --hmm_db $_ --tmp_prefix $HMM_db{$_} proteins_for_filtering.fasta >> validation_hmmscan.tab 2>> para_hmmscan.1.log; $dirname/bin/para_hmmscan $config{'para_hmmscan'} --chunk 1 --outformat --cpu $para_hmmscan_cpu --no_cut_ga --hmm_db $_ --tmp_prefix $HMM_db{$_} proteins_for_filtering.fasta >> validation_hmmscan.tab 2>> para_hmmscan.2.log; ";
        }
    }
    if ( $BLASTP_db ) {
        foreach ( sort keys %BLASTP_db ) {
            $cmdString2 .= "diamond blastp $config{'diamond'} --outfmt 5 --db $_ --query proteins_for_filtering.fasta --out validation_blastp_$BLASTP_db{$_}.xml --threads $cpu &>> diamond_blastp.log; ";
            $cmdString3 = "$dirname/bin/parsing_blast_result.pl $config{'parsing_blast_result.pl'} --out-hit-confidence validation_blastp_$BLASTP_db{$_}.xml >> validation_blastp.tab; ";
        }
    }
    else {
        $cmdString2 = "diamond makedb --db homolog --in ../homolog.fasta &> diamond_makedb.log; diamond blastp $config{'diamond'} --outfmt 5 --db homolog --query proteins_for_filtering.fasta --out validation_blastp.xml --threads $cpu &> diamond_blastp.log";
        $cmdString3 = "$dirname/bin/parsing_blast_result.pl $config{'parsing_blast_result.pl'} --out-hit-confidence validation_blastp.xml > validation_blastp.tab";
    }
    unless ( -e "07.validating.ok" ) {
        open OUT, ">", "validation_hmmscan.tab" or die "Can not create file validation_hmmscan.tab, $!", close OUT;
        open OUT, ">", "validation_blastp.tab" or die "Can not create file validation_blastp.tab, $!", close OUT;
        if ( $HMM_db ) {
            print STDERR (localtime) . ": CMD: $cmdString1\n";
            system("$cmdString1") == 0 or die "failed to execute: $cmdString1\n";
        }
        print STDERR (localtime) . ": CMD: $cmdString2\n";
        system("$cmdString2") == 0 or die "failed to execute: $cmdString2\n";
        print STDERR (localtime) . ": CMD: $cmdString3\n";
        system("$cmdString3") == 0 or die "failed to execute: $cmdString3\n";
        open OUT, ">", "07.validating.ok" or die $!; close OUT;
    }
    else {
        if ( $HMM_db ) {
            print STDERR "CMD(Skipped): $cmdString1\n";
        }
        if ( $BLASTP_db ) {
            print STDERR "CMD(Skipped): $cmdString2\n";
            print STDERR "CMD(Skipped): $cmdString3\n";
        }
    }

    # 6.8 / 6.9 根据HMM和BLASTP验证结果对基因模型进行过滤。
    # 获得验证通过的转录本ID
    my $cmdString = "$dirname/bin/get_valid_transcriptID $config{'get_valid_transcriptID'} validation_hmmscan.tab validation_blastp.tab > transcriptID_validating_passed.tab 2> get_valid_transcriptID.log";
    unless ( -e "08.get_valid_transcriptID.ok" ) {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "08.get_valid_transcriptID.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }
    my $cmdString = "$dirname/bin/get_valid_geneModels $config{'get_valid_geneModels'} --out_prefix geneModels.h transcriptID_for_filtering.txt transcriptID_validating_passed.tab geneModels.gb.gff3 geneModels.ge.gff3 geneModels.gf.gff3 2> get_valid_geneModels.log";
    unless ( -e "09.get_valid_geneModels.ok" ) {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "09.get_valid_geneModels.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    # 6.10 再次对基因模型进行首尾填补。
    $cmdString = "$dirname/bin/fillingEndsOfGeneModels $config{'fillingEndsOfGeneModels'} $genome geneModels.h.coding.gff3 > geneModels.i.coding.gff3 2> fillingEndsOfGeneModels.2.log";
    unless ( -e "10.fillingEndsOfGeneModels" ) {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "10.fillingEndsOfGeneModels" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    chdir "../";
    open OUT, ">", "6.combineGeneModels.ok" or die $!; close OUT;
}
else {
    print STDERR "Skip Step 6 for the file 6.combineGeneModels.ok exists\n";
}

# Step 7: OutPut
print STDERR "\n============================================\n";
print STDERR "Step 7: OutPut " . "(" . (localtime) . ")" . "\n";
chdir "../";
$pwd = `pwd`; print STDERR "PWD: $pwd";

# 7.1 输出GFF3格式文件基因结构注释信息
$cmdString = "$dirname/bin/GFF3Clear --GFF3_source GETA --gene_prefix $gene_prefix --gene_code_length 6 --genome $genome $out_prefix.tmp/6.combineGeneModels/geneModels.i.coding.gff3 > $out_prefix.geneModels.gff3 2> /dev/null";
print STDERR (localtime) . ": CMD: $cmdString\n";
system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

$cmdString = "$dirname/bin/GFF3_extract_bestGeneModels $out_prefix.geneModels.gff3 > $out_prefix.bestGeneModels.gff3 2> $out_prefix.AS_num_of_codingTranscripts.stats";
print STDERR (localtime) . ": CMD: $cmdString\n";
system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

$cmdString = "$dirname/bin/GFF3Clear --GFF3_source GETA --gene_prefix ${out_prefix}ncGene  --gene_code_length 6 --genome $genome --no_attr_add $out_prefix.tmp/6.combineGeneModels/geneModels.h.lncRNA.gff3 > $out_prefix.geneModels_lncRNA.gff3 2> /dev/null";
print STDERR (localtime) . ": CMD: $cmdString\n";
system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

$cmdString = "$dirname/bin/GFF3Clear --GFF3_source GETA --gene_prefix ${out_prefix}lqGene --gene_code_length 6 --genome $genome --no_attr_add --coverage 0.6 $out_prefix.tmp/6.combineGeneModels/geneModels.h.lowQuality.gff3 > $out_prefix.geneModels_lowQuality.gff3 2> /dev/null";
print STDERR (localtime) . ": CMD: $cmdString\n";
system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

# 7.2 输出GTF文件和基因的序列信息
$cmdString = "$dirname/bin/gff3ToGtf.pl $genome $out_prefix.geneModels.gff3 > $out_prefix.geneModels.gtf 2> /dev/null";
print STDERR (localtime) . ": CMD: $cmdString\n";
system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

$cmdString = "$dirname/bin/gff3ToGtf.pl $genome $out_prefix.bestGeneModels.gff3 > $out_prefix.bestGeneModels.gtf 2> /dev/null";
print STDERR (localtime) . ": CMD: $cmdString\n";
system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

#$cmdString = "$dirname/bin/eukaryotic_gene_model_statistics.pl $out_prefix.bestGeneModels.gtf $genome $out_prefix &> $out_prefix.geneModels.stats";
$cmdString = "$dirname/bin/gff3_to_sequences.pl --out_prefix $out_prefix --only_gene_sequences --only_coding_gene_sequences --only_first_isoform --genetic_code 1 $genome $out_prefix.geneModels.gff3 > $out_prefix.geneModels.stats";
print STDERR (localtime) . ": CMD: $cmdString\n";
system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

# 7.3 输出重复序列信息及其统计结果、RepeatModeler软件构建的重复序列数据库和masked genome sequence
open IN, "$out_prefix.tmp/0.RepeatMasker/genome.masked.fasta" or die "Can not open file $out_prefix.tmp/0.RepeatMasker/genome.masked.fasta, $!";
open OUT, ">", "$out_prefix.maskedGenome.fasta" or die "Can not create file $out_prefix.maskedGenome.fasta, $!";
print OUT <IN>;
close IN; close OUT;

open IN, "$out_prefix.tmp/0.RepeatMasker/genome.repeat.stats" or die "Can not open file $out_prefix.tmp/0.RepeatMasker/genome.repeat.stats, $!";
open OUT, ">", "$out_prefix.repeat.stats" or die "Can not create file $out_prefix.repeat.stats, $!";
print OUT <IN>;
close IN; close OUT;

open IN, "$out_prefix.tmp/0.RepeatMasker/genome.repeat.gff3" or die "Can not open file $out_prefix.tmp/0.RepeatMasker/genome.repeat.gff3, $!";
open OUT, ">", "$out_prefix.repeat.gff3" or die "Can not create file $out_prefix.repeat.gff3, $!";
print OUT <IN>;
close IN; close OUT;

if ( $RM_lib ) {
    $cmdString = "ln -sf $RM_lib $out_prefix.repeat.lib";
    unless ( -s "$out_prefix.repeat.lib" ) {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    }
    else {
        print STDERR "Skipped CMD: $cmdString\n";
    }
}
else {
    open IN, "$out_prefix.tmp/0.RepeatMasker/repeatModeler/species-families.fa" or die "Can not open file $out_prefix.tmp/0.RepeatMasker/repeatModeler/species-families.fa, $!";
    open OUT, ">", "$out_prefix.repeat.lib" or die "Can not create file $out_prefix.repeat.lib, $!";
    print OUT <IN>;
    close IN; close OUT;
}

# 7.4 输出转录本、同源蛋白和Augustus的基因预测结果
if (($pe1 && $pe2) or $single_end) {
    open IN, "$out_prefix.tmp/1.NGSReads_prediction/c.transcript/transfrag.alignment.gff3" or die "Can not open file $out_prefix.tmp/1.NGSReads_prediction/c.transcript/transfrag.alignment.gff3, $!";
    open OUT, ">", "$out_prefix.transfrag_alignment.gff3" or die "Can not create file $out_prefix.transfrag_alignment.gff3, $!";
    print OUT <IN>;
    close IN; close OUT;

    open IN, "$out_prefix.tmp/1.NGSReads_prediction/c.transcript/transfrag.genome.gff3" or die "Can not open file $out_prefix.tmp/1.NGSReads_prediction/c.transcript/transfrag.genome.gff3, $!";
    open OUT, ">", "$out_prefix.transfrag_prediction.gff3" or die "Can not create file $out_prefix.transfrag_prediction.gff3, $!";
    print OUT <IN>;
    close IN; close OUT;
}

if ($protein) {
    open IN, "$out_prefix.tmp/2.homolog_prediction/homolog_prediction.gff3" or die "Can not open file $out_prefix.tmp/2.homolog/genewise.gff3, $!";
    open OUT, ">", "$out_prefix.homolog_prediction.gff3" or die "Can not create file $out_prefix.homolog_prediction.gff3, $!";
    print OUT <IN>;
    close IN; close OUT;
}

$cmdString = "$dirname/bin/GFF3Clear --genome $genome --no_attr_add $out_prefix.tmp/4.augustus/augustus.gff3 > $out_prefix.augustus_prediction.gff3";
print STDERR (localtime) . ": CMD: $cmdString\n";
system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

# 7.5 输出GETA流程信息，用于追踪基因预测结果的可靠性
# (1) 获取基因组重复序列统计信息
open OUT, ">", "$out_prefix.gene_prediction.summary" or die "Can not create file $out_prefix.gene_prediction.summary, $!";
open IN, "$out_prefix.tmp/0.RepeatMasker/genome.repeat.stats" or die "Can not open file $out_prefix.tmp/0.RepeatMasker/genome.repeat.stats, $!";
while (<IN>) {
    print OUT $_ if m/^Genome Size/;
    print OUT "$_\n" if m/^Repeat Ratio/;
}
close IN;
# (2) 获取转录组二代测序数据和基因组的匹配率
if ( -e "$out_prefix.tmp/2.hisat2/hisat2.log" ) {
    open IN, "$out_prefix.tmp/2.hisat2/hisat2.log" or die $!;
    while (<IN>) {
        print OUT "The alignment rate of RNA-Seq reads is: $1\n\n" if m/^(\S+) overall alignment rate/;
    }
    close IN;
# (3) 获取转录组数据预测的基因数量的统计信息
    open IN, "$out_prefix.tmp/1.NGSReads_prediction/c.transcript/transfrag.genome.gff3" or die "Can not open file $out_prefix.tmp/1.NGSReads_prediction/c.transcript/transfrag.genome.gff3, $!";
    my ($num1, $num2, $num3, $num4, $num5);
    while (<IN>) {
        if (m/\tgene\t/) {
            $num1 ++;
            if (m/Integrity=complete/) {
                $num2 ++;
            }
            elsif (m/Integrity=5prime_partial/) {
                $num3 ++;
            }
            elsif (m/Integrity=3prime_partial/) {
                $num4 ++;
            }
            elsif (m/Integrity=internal/) {
                $num5 ++;
            }
        }
    }
    print OUT "$num1 genes were predicted by Transfrag, including $num2 complete, $num3 5prime_partial, $num4 3prime_partial, $num5 internal genes.\n\n";
}
# (4) 获取同源蛋白预测的基因数量的统计信息
if ( $protein ) {
    open IN, "$out_prefix.tmp/2.homolog_prediction/homolog_prediction.gff3" or die "Can not open file $out_prefix.tmp/2.homolog/genewise.gff3, $!";
    my ($num1, $num2, $num3, $num4, $num5);
    while (<IN>) {
        if (m/\tgene\t/) {
            $num1 ++;
            if (m/Integrity=complete/) {
                $num2 ++;
            }
            elsif (m/Integrity=5prime_partial/) {
                $num3 ++;
            }
            elsif (m/Integrity=3prime_partial/) {
                $num4 ++;
            }
            elsif (m/Integrity=internal/) {
                $num5 ++;
            }
        }
    }
    print OUT "$num1 genes were predicted by Homolog, including $num2 complete, $num3 5prime_partial, $num4 3prime_partial, $num5 internal genes.\n\n";
}
if (-e "$out_prefix.tmp/4.augustus/augustus.2.gff3") {
    open IN, "$out_prefix.tmp/4.augustus/training_again/secondtest.out" or die "Can not open file $out_prefix.tmp/4.augustus/training_again/secondtest.out, $!";
}
else {
    open IN, "$out_prefix.tmp/4.augustus/training/secondtest.out" or die "Can not open file $out_prefix.tmp/4.augustus/training/secondtest.out, $!";
}
# (5) 获取AUGUSTUS Training的准确率信息和预测基因数量统计
my ($accuary1, $accuary2, $accuary3, $accuary4, $accuary5, $accuary6, $out);
while (<IN>) {
    if (m/^nucleotide level/) {
        @_ = m/([\d\.]+)/g;
        $out .= "nucleotide_level\t$_[0]\t$_[1]\n";
        ($accuary1, $accuary2) = ($_[-2], $_[-1]);
    }
    elsif (m/^exon level/) {
        @_ = m/([\d\.]+)/g;
        $out .= "exon_level\t$_[-2]\t$_[-1]\n";
        ($accuary3, $accuary4) = ($_[-2], $_[-1]);
    }
    elsif (m/^gene level/) {
        @_ = m/([\d\.]+)/g;
        $out .= "gene_level\t$_[-2]\t$_[-1]\n";
        ($accuary5, $accuary6) = ($_[-2], $_[-1]);
    }
}
my $accuary = ($accuary1 * 3 + $accuary2 * 2 + $accuary3 * 4 + $accuary4 * 3 + $accuary5 * 2 + $accuary6 * 1) / 15;
$accuary = int($accuary * 10000) / 100;
$out = "The accuary of AUGUSTUS Training is $accuary\%.\nLevel\tSensitivity\tSpecificity\n$out";
print OUT $out;
my $num_of_gene_predicted_by_AUGUSTUS = 0;
open IN, "$out_prefix.tmp/4.augustus/augustus.gff3" or die "Can not open file $out_prefix.tmp/4.augustus/augustus.gff3, $!";
while (<IN>) {
    $num_of_gene_predicted_by_AUGUSTUS ++ if m/\tgene\t/;
}
print OUT "$num_of_gene_predicted_by_AUGUSTUS genes were predicted by AUGUSTUS.\n\n";

# (6) 获取基因预测整合过滤的统计信息
print OUT "Statistics of the combination of 3 gene prediction methods and filtration of gene models:\n";
open IN, "$out_prefix.tmp/6.combineGeneModels/geneModels.a.gff3" or die "Can not open file $out_prefix.tmp/6.combineGeneModels/geneModels.a.gff3, $!";
my ($num_of_gene_a, $num_of_gene_b, $num_of_gene_c, $num_of_gene_d) = (0, 0, 0, 0);
while (<IN>) {
    $num_of_gene_a ++ if (m/\tgene\t/ && m/augustus/);
    $num_of_gene_c ++ if (m/\tgene\t/ && m/transfrag/);
    $num_of_gene_d ++ if (m/\tgene\t/ && m/genewise/);
}
close IN;
open IN, "$out_prefix.tmp/6.combineGeneModels/geneModels.b.gff3" or die "Can not open file $out_prefix.tmp/6.combineGeneModels/geneModels.b.gff3, $!";
while (<IN>) {
    $num_of_gene_b ++ if m/\tgene\t/;
}
close IN;
print OUT "(1) After first round of combination in which the AUGUSTUS results were mainly used, $num_of_gene_a genes models were supported by enough evidences, $num_of_gene_c genes models were come from transcript, $num_of_gene_d genes models were come from homolog, and $num_of_gene_b genes did not supported by enough evidences.\n";

open IN, "$out_prefix.tmp/6.combineGeneModels/picked_evidence_geneModels.log" or die "Can not open file $out_prefix.tmp/6.combineGeneModels/picked_evidence_geneModels.log, $!";
my ($number1, $number2, $number3) = (0, 0, 0);
$_ = <IN>; $number1 = $1 if m/(\d+)/;
<IN>; <IN>; <IN>;
$_ = <IN>; $number2 = $1 if m/^.*?\d+.*?(\d+)/;
$_ = <IN>; $number3 = $1 if m/(\d+)/;
close IN;
print OUT "(2) In the second round of combination, $number1 evidence gene models were processed, and $number3 accurate gene models were picked out and replaced the genes predicted by AUGUSTUS, $number2 of which had the same CDS structures with the gene models predicted by AUGUSTUS.\n";

open IN, "$out_prefix.tmp/6.combineGeneModels/fillingEndsOfGeneModels.1.log" or die "Can not open file $out_prefix.tmp/6.combineGeneModels/fillingEndsOfGeneModels.1.log, $!";
my @line = <IN>;
close IN;
my @number1 = $line[-2] =~ m/(\d+)/g;
my @number2 = $line[-1] =~ m/(\d+)/g;
print OUT "(3) After the two steps of combination, $number1[0] gene models with enough evidence supported were predicted. $number1[1] gene models were complete; $number1[2] gene models were uncomplete. $number2[0] uncomplete gene models can be filled to complete, and $number2[1] can not.\n";

open IN, "$out_prefix.tmp/6.combineGeneModels/get_valid_geneModels.log" or die "Can not open file $out_prefix.tmp/6.combineGeneModels/get_valid_geneModels.log, $!";
$_ = <IN>; 
close IN;
my @number = m/(\d+)/g;
print OUT "(4) HMM and BLASTP validation were performed to $number[2] protein sequences of $number[1] genes, and $number[5] protein sequences of $number[4] genes had valid alignment results. There are $number[3] accurate gene models which did not need validation.\n";

open IN, "$out_prefix.geneModels.gff3" or die "Can not open file $out_prefix.geneModels.gff3, $!";
my (%gene_ID, %gene2transcript, $num_augustus, $num_transfrag, $num_genewise);
while (<IN>) {
    if (m/\tgene\t.*ID=([^;]+)/) {
        $gene_ID{$1} = 1;
        if (m/Source=([a-zA-Z]+)/) {
            $num_augustus ++ if $1 eq 'augustus';
            $num_transfrag ++ if $1 eq 'transfrag';
            $num_genewise ++ if $1 eq 'genewise';
        }
    }
    elsif ( m/ID=([^;]+).*Parent=([^;]+)/ && exists $gene_ID{$2} ) {
        $gene2transcript{$2}{$1} = 1;
    }
}
close IN;
my $as_gene_num = 0;
foreach ( keys %gene2transcript ) {
    my @as_gene = keys %{$gene2transcript{$_}};
    $as_gene_num ++ if @as_gene >= 2;
}
my $num_of_gene_ID = 0; $num_of_gene_ID = %gene_ID;
print OUT "(5) Finally, $num_of_gene_ID coding gene models were obtained, $as_gene_num of which had alternative splicing, $num_augustus gene models were come from AUGUSTUS prediction, $num_transfrag gene models were come from transcript, $num_genewise gene models were come from homolog.\n";

my ($num_of_gene1, $num_of_gene2, $num_of_gene3) = (0, 0, 0);
open IN, "$out_prefix.geneModels_lncRNA.gff3" or die "Can not open file $out_prefix.geneModels_lncRNA.gff3, $!";
while (<IN>) {
    if (m/\tgene\t.*ID=([^;]+)/) {
        $num_of_gene1 ++;
    }
}
close IN;
open IN, "$out_prefix.geneModels_lowQuality.gff3" or die "Can not open file $out_prefix.geneModels_lowQuality.gff3, $!";
while (<IN>) {
    if (m/\tgene\t.*ID=([^;]+)/) {
        $num_of_gene2 ++;
    }
}
close IN;
$num_of_gene3 = $num_of_gene1 + $num_of_gene2;
print OUT "(6) $num_of_gene3 gene models were filtered, and $num_of_gene1 of which had lncRNA transcripts.\n";

if ( $delete_unimportant_intermediate_files ) {
    print STDERR "Due to the --delete_unimportant_intermediate_files was set, so the unimportant and large intermediate files are being deleted...\n";
    # 删除 0.RepeatMasker 文件夹下的数据
    $cmdString = "rm -rf $out_prefix.tmp/0.RepeatMasker/repeatMasker $out_prefix.tmp/0.RepeatMasker/repeatModeler";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or warn "failed to execute: $cmdString\n";
    # 删除 1.trimmomatic 文件夹下的数据
    $cmdString = "rm -rf $out_prefix.tmp/1.trimmomatic/*.fastq $out_prefix.tmp/1.trimmomatic/command*";
    # 删除 2.hisat2 文件夹下的数据
       $cmdString = "rm -rf $out_prefix.tmp/2.hisat2/genome* $out_prefix.tmp/2.hisat2/*.sam $out_prefix.tmp/2.hisat2/*.bam $out_prefix.tmp/2.hisat2/*.ok $out_prefix.tmp/2.hisat2/hisat2-build.log";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or warn "failed to execute: $cmdString\n";
    # 删除 1.NGSReads_prediction/c.transcript 文件夹下的数据
    $cmdString = "rm -rf $out_prefix.tmp/1.NGSReads_prediction/c.transcript/command* $out_prefix.tmp/1.NGSReads_prediction/c.transcript/*.ok $out_prefix.tmp/1.NGSReads_prediction/c.transcript/pipeliner* $out_prefix.tmp/1.NGSReads_prediction/c.transcript/proteins.fasta $out_prefix.tmp/1.NGSReads_prediction/c.transcript/split* $out_prefix.tmp/1.NGSReads_prediction/c.transcript/transdecoder2ORF.gff3 $out_prefix.tmp/1.NGSReads_prediction/c.transcript/transfrag.gff3 $out_prefix.tmp/1.NGSReads_prediction/c.transcript/transfrag.gtf $out_prefix.tmp/1.NGSReads_prediction/c.transcript/transfrag.noStrand.* $out_prefix.tmp/1.NGSReads_prediction/c.transcript/transfrag.strand.* $out_prefix.tmp/1.NGSReads_prediction/c.transcript/transfrag.stats $out_prefix.tmp/1.NGSReads_prediction/c.transcript/transfrag.transdecoder.gff3";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or warn "failed to execute: $cmdString\n";
    # 删除 2.homolog 文件夹下的数据
    $cmdString = "rm -rf $out_prefix.tmp/2.homolog/*.ok $out_prefix.tmp/2.homolog/blast* $out_prefix.tmp/2.homolog/command* $out_prefix.tmp/2.homolog/geneRegion_genewise.tmp $out_prefix.tmp/2.homolog/genewise.gff $out_prefix.tmp/2.homolog/genome* $out_prefix.tmp/2.homolog/homolog_gene_region.tab $out_prefix.tmp/2.homolog/makeblastdb.log $out_prefix.tmp/2.homolog/out* $out_prefix.tmp/2.homolog/tmp_for_genome_splitting_and_joining $out_prefix.tmp/2.homolog/para_blast.log";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or warn "failed to execute: $cmdString\n";
    # 删除 4.augustus 文件夹下的数据
    $cmdString = "rm -rf $out_prefix.tmp/4.augustus/training/ati* $out_prefix.tmp/4.augustus/training/badgenes.lst $out_prefix.tmp/4.augustus/training/blank* $out_prefix.tmp/4.augustus/training/com* $out_prefix.tmp/4.augustus/training/etraining* $out_prefix.tmp/4.augustus/training/genes.* $out_prefix.tmp/4.augustus/training/gff* $out_prefix.tmp/4.augustus/training/GFF3Clear.log $out_prefix.tmp/4.augustus/training/hmm* $out_prefix.tmp/4.augustus/training/new* $out_prefix.tmp/4.augustus/training/tmp* $out_prefix.tmp/4.augustus/training/training.gb.onlytrain $out_prefix.tmp/4.augustus/training/*.ok $out_prefix.tmp/4.augustus/aug_* $out_prefix.tmp/4.augustus/augustus.?.gff3 $out_prefix.tmp/4.augustus/hints.gff $out_prefix.tmp/4.augustus/*.ok $out_prefix.tmp/4.augustus/extrinsic.cfg";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or warn "failed to execute: $cmdString\n";
    # 删除 6.combineGeneModels 文件夹下的数据
    if ( $keep_all_files_of_6thStep_combineGeneModels ) {
        print STDERR "Due to the --keep_all_files_of_6thStep_combineGeneModels was set, the $out_prefix.tmp/FEN00I06.tmp/6.combineGeneModels Directory was keeped.\n";
    }
    else {
        $cmdString = "rm -rf $out_prefix.tmp/FEN00I06.tmp/6.combineGeneModels*";
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or warn "failed to execute: $cmdString\n";
    }
}

print STDERR "\n============================================\n";
print STDERR "GETA complete successfully! " . "(" . (localtime) . ")" . "\n\n";


sub detecting_dependent_softwares {
    # 检测依赖的软件
    print STDERR "\n============================================\n";
    print STDERR "Detecting the dependent softwares:\n";
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
    if ($software_info =~ m/Runtime Environment/) {
        print STDERR "java:\tOK\n";
    }
    else {
        die "java:\tFailed\n\n";
    }
    # 检测HISAT2
    $software_info = `hisat2 --version`;
    if ($software_info =~ m/version 2.(\d+)\.(\d+)/) {
        print STDERR "HISAT2:\tOK\n";
    }
    else {
        die "HISAT2:\tFailed\n\n";
    }
    # 检测samtools
    $software_info = `samtools --version`;
    if ($software_info =~ m/samtools 1.(\d+)/) {
        print STDERR "samtools:\tOK\n";
    }
    else {
        die "samtools:\tFailed\n\n";
    }
    # 检测hmmer
    $software_info = `hmmscan -h`;
    if ($software_info =~ m/HMMER 3.(\d+)/) {
        print STDERR "hmmer:\tOK\n";
    }
    else {
        die "hmmer:\tFailed\n\n";
    }
    # 检测diamond
    $software_info = `diamond version`;
    if ($software_info =~ m/diamond version/) {
        print STDERR "diamond:\tOK\n";
    }
    else {
        die "diamond:\tFailed\n\n";
    }
    print STDERR "============================================\n\n";
    my $pwd = `pwd`; chomp($pwd);
    print STDERR "PWD: $pwd\n";
    print STDERR (localtime) . ": CMD: $0 $command_line_geta\n\n"; 
}


# 检测服务器剩余可用内存容量，结果单位是 kB。
sub get_MemAvailable {
    open IN, "/proc/meminfo" or die "Can not open file /proc/meminfo, $!";
    my $MemAvailable;
    while (<IN>) {
        if (m/MemAvailable:\s*(\d+)\s*kB/) {
            $MemAvailable = $1;
            next;
        }
    }
    close IN;
    return $MemAvailable;
}
