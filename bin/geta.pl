#!/usr/bin/perl
use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw/abs_path getcwd cwd/;
use File::Basename;

my $command_line_geta = join " ", @ARGV;
my $bin_path = dirname($0);
my $software_dir = $bin_path; $software_dir =~ s/\/bin$//;

my $usage_chinese = <<USAGE;

GETA (Genome-wide Electronic Tool for Annotation) 是一个对全基因组序列进行基因预测的流程软件。输入基因组序列、转录组二代测序原始数据和临近物种全基因组蛋白序列，即可一个命令快速得到准确的基因组注释GFF3结果文件。本流程的两大特色：（1）预测准确，程序输出结果中基因数量正常、BUSCO检测的完整性高、基因模型的exon边界准确；（2）运行简单，一个命令全自动化运行得到最终结果。

软件当前版本：2.7.1。

Usage:
    perl $0 [options]

For example:
    perl $0 --genome genome.fasta --RM_species_Dfam Embryophyta --RM_species_RepBase Embryophyta --pe1 libA.1.fq.gz,libB.1.fq.gz --pe2 libA.2.fq.gz,libB.2.fq.gz --protein homolog.fasta --augustus_species GETA_genus_species --HMM_db /opt/biosoft/bioinfomatics_databases/Pfam/PfamA,/opt/biosoft/bioinfomatics_databases/Pfam/PfamB --config /opt/biosoft/geta/conf_for_big_genome.txt --out_prefix out --gene_prefix GS01Gene --cpu 120

Parameters:
[INPUT]
    --genome <string>    default: None, required
    输入需要进行注释的基因组序列FASTA格式文件。当输入了屏蔽重复序列的基因组文件时，可以不使用--RM_species_Dfam、--RM_species_RepBase和--RM_lib参数，并添加--no_RepeatModeler参数，从而跳过重复序列分析与屏蔽步骤。此时，推荐再输入FASTA文件中仅对转座子序列使用碱基N进行硬屏蔽，对简单重复或串联重复使用小写字符进行软屏蔽。

    --RM_species_Dfam | --RM_species <string>    default: None
    输入一个物种或分类名称，以利于调用RepeatMasker程序运用Dfam数据库中相应物种分类的HMM数据信息进行重复序列分析。可以根据文件$software_dir/RepeatMasker_lineage.txt的内容查找适用于本参数的物种分类名称。例如，设置Eukaryota适合真核生物、Viridiplantae适合植物、Metazoa适合动物、Fungi适合真菌。当使用了该参数，需要先安装好RepeatMasker软件并配置好Dfam数据库。需要注意的是：Dfam数据库因为较大被分隔成了9份；而RepeatMasker默认仅带有Dfam数据库第0个root部分，适合哺乳动物或真菌等物种；若有需要，则考虑下载适合待分析物种的Dfam数据库并配置到RepeatMasker；例如，对植物物种进行分析则需要使用Dfam数据库第5个Viridiplantae部分，对膜翅目昆虫进行分析则需要使用Dfam数据库第7个Hymenoptera部分。

    --RM_species_RepBase <string>    default: None
    输入一个物种名称，运用调用RepeatMasker程序运用RepBase数据库中相应物种分类的重复核酸序列信息进行重复序列分析。可以根据文件$software_dir/RepeatMasker_lineage.txt的内容查找适用于本参数的物种分类名称。例如，设置Eukaryota适合真核生物、Viridiplantae适合植物、Metazoa适合动物、Fungi适合真菌。当输入了该参数，需要先安装好RepeatMasker软件并配置好RepBase数据库。需要注意的是：RepBase数据库最后一个版本20181026已经较老，且其包含的重复序列数据量较少；RepBase也不再提供免费下载。

    --RM_lib <string>    default: None
    输入一个FASTA格式文件，运用其中的重复序列数据信息进行全基因组重复序列分析。该文件往往是RepeatModeler软件对全基因组序列进行分析的结果，表示全基因组上的重复序列。GETA软件默认会调用RepeatModeler对全基因组序列进行分析，得到物种自身的重复序列数据库后，再调用RepeatMakser进行重复序列分析。当添加该参数后，则程序跳过费时的RepeatModler步骤，能极大减少程序运行时间。此外，程序支持--RM_species_Dfam、--RM_species_RepBase和--RM_lib参数同时使用，则依次使用多种方法进行重复序列分析，最终能合并多个结果并认可任一方法的结果。

    --no_RepeatModeler    default: None
    添加该参数后，程序不会运行RepeatModeler步骤，适合直接输入屏蔽了重复序列基因组文件的情形。

    --pe1 <string> --pe2 <string>    default: None
    输入二代双末端测序的两个FASTQ格式文件。参数支持输入多对数据文件，使用逗号对不同文库的FASTQ文件路径进行分隔即可。参数也支持输入.gz格式的压缩文件。

    --se <string>    default: None
    输入二代单端测序的FASTQ格式文件。参数支持输入多个数据文件，使用逗号对不同文库的FASTQ文件路径进行分隔即可。参数也支持输入.gz格式的压缩文件。

    --sam <string>    default: None
    输入二代测序SAM格式数据文件。参数支持输入多个数据文件，使用逗号对不同的SAM文件路径进行分隔即可。参数也支持输入.bam格式的压缩文件。此外，程序支持--pe1/--pe2、--se和--sam这三个参数全部或部分使用，则利用所有输入的数据进行基因组比对，再获取转录本序列进行基因预测。

    --strand_specific    default: None
    添加该参数后，则认为所有输入的二代数据是链特异性测序数据，程序则仅在转录本正义链上预测基因。使用链特异性测序数据并设置本参数后，当两相邻基因有重叠时，能准确预测其边界。

    --protein <string>    default: None
    输入临近物种的全基因组蛋白序列。推荐使用多个（3~10个）物种的全基因组同源蛋白序列。推荐修改蛋白序列的名称，在原名称尾部追加以Species字符开头的物种信息。例如蛋白序列名称为XP_002436309.2，则修改为XP_002436309_2_SpeciesSorghumBicolor。这样有利于在一个基因区域保留更多物种的同源匹配结果，有利于基因预测的准确性。可以考虑使用GETA中附带的fasta_remove_redundancy.pl整合多个物种的全基因组蛋白序列，去冗余时也能附加物种信息。此外，使用的物种数量越多，则能预测的更多准确的基因模型，但是越消耗计算时间。需要注意的是，同源蛋白和二代测序数据，本软件要求至少需要输入其中一种数据用于有证据支持的基因预测。

    --augustus_species <string>    default: None
    输入一个AUGUSTUS物种名称，则程序在利用转录本或同源蛋白预测的基因模型进行AUGUSTUS Training时，从已有的物种模开始训练或重头训练新的物种模型。若输入的AUGUSTUS物种模型存在，则在其基础上对其进行优化；若不存在，则生成新的AUGUSTUS物种模型，并进行参数优化。程序进行AUGUSTUS Training需要安装好AUGUSTUS软件，并设置好\$AUGUSTUS_CONFIG_PATH环境变量。程序运行完毕后，在临时文件夹中有生成名称为本参数值的AUGUSTUS Training物种配置文件夹；若对\$AUGUSTUS_CONFIG_PATH路径中指定的species文件夹有写入权限，则将生成的物种配置文件夹拷贝过去。若不输入本参数信息，则程序自动设置本参数的值为“GETA + 基因组FASTA文件名称前缀 + 日期 + 进程ID”。

    --HMM_db <string>    default: None
    输入HMM数据库路径，用于对基因模型进行过滤。参数支持输入多个数据库路径，使用逗号进行分隔。当使用多个HMM数据库时，程序过滤在所有数据库中都没有匹配的基因模型。

    --BLASTP_db <string>    default: None
    输入diamond数据库路径，用于对基因模型进行过滤。参数支持输入多个数据库路径，使用逗号进行分隔。当使用多个diamond数据库时，程序过滤在所有数据库中都没有匹配的基因模型。若不设置该参数，则以--protein参数输入的同源蛋白序列构建diamond数据库，进行基因模型过滤。

    --config <string>    default: None
    输入一个参数配置文件路径，用于设置本程序调用的其它命令的详细参数。若不设置该参数，当基因组>1GB时，自动使用软件安装目录中的conf_for_big_genome.txt配置文件；当基因组<50MB时，自动使用软件安装目录中的conf_for_small_genome.txt配置文件；当基因组在50MB~1GB之间时，使用默认参数配置。此外，当软件预测的基因数量异常时往往要修改基因模型的过滤阈值。此时，通过修改软件安装目录中的conf_all_defaults.txt文件内容生成新的配置文件，并输入给本参来再次运行GETA流程。

    --BUSCO_lineage_dataset <string>    default: None
    输入BUSCO数据库路径，则程序额外对基因预测得到的全基因组蛋白序列进行BUSCO分析。本参数支持输入多个BUSCO数据库路径，使用逗号进行分隔，则分别利用多个数据库进行分析。可以根据$software_dir/BUSCO_lineages_list.2021-12-14.txt文件内容选择合适的BUSCO数据库。BUSCO的结果输出到7.output_gene_models子目录下和gene_prediction.summary文件中。

[OUTPUT]
    --out_prefix <string>    default: out
    设置输出文件或临时文件夹前缀。

    --gene_prefix <string>    default: gene
    设置输出GFF3文件中的基因名称前缀。

    --chinese_help    default: None
    使用该参数后，程序给出中文用法并退出。

    --help    default: None
    display this help and exit.

[Settings]
    --cpu <int>    default: 4
    设置程序运行使用的CPU线程数。

    --put_massive_temporary_data_into_memory    default: None
    设置将海量的临时文件存放到内存中。这样能避免磁盘I/O不足而造成程序运行减缓，但需要消耗更多内存。本流程在很多步骤中对数据进行了分割，再通过并行化来加速计算，但这对磁盘形成了极大的I/O负荷。因此，当磁盘性能较差时会严重影响计算速度。若系统内存充足，推荐添加本参数，从而将海量的临时数据存放到代表内存的/dev/shm文件夹下，以加速程序运行。此外，程序在数据分割和并行化步骤运行完毕后，会自动删除/dev/shm中的临时数据以释放内存。

    --gene_predicted_by_unmasked_genome    default: None
    添加该参数后，程序在利用NGS read、homology和AUGUSTUS进行基因预测时，使用输入的基因组序列进行基因预测。而默认程序对输入的基因组序列进行重复序列屏蔽，再使用屏蔽了的基因组序列采用三种算法进行基因预测。

    --genetic_code <int>    default: 1
    设置遗传密码。该参数对应的值请参考NCBI Genetic Codes: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi。本参数主要生效于同源蛋白进行基因预测的步骤，或对基因模型首尾进行强制补齐时使用的起始密码子和终止密码子信息的情形。

    --homolog_prediction_method <string>    default: all
    设置使用同源蛋白进行基因预测的方法。其值可以设定为exonerate、genewise、gth或all。若想使用多种方法进行分析，则输入使用逗号分割的多个值；若使用所有三种方法进行分析，可以设置--method参数值为all。使用的方法越多，越消耗计算时间，但结果更好。三种方法中：exonerate和genewise的准确性结果比较一致（），但gth方法预测基因模型的sensitivity下降很多，specificity提高很多。以三种方法对Oryza sativa基因组的预测为例，其预测结果的准确性如下表所示。和NCBI上标准的共28736个基因模型的注释结果进行比较，评估四个准确性值：基因水平sensitivity、基因水平specificity、exon水平sensitivity、exon水平specificity。可以看出，使用多种方法进行基因预测、合并结果后再过滤，得到的基因模型数量能接近真实的基因数量，且结果较准确。此外，使用本参数的优先级更高，能覆盖--config指定参数配置文件中homolog_prediction的参数值。
    方法       基因数量    gene_sensitivity    gene_specificity    exon_sensitivity    exon_specificity
    exonerate  31310       46.17%              42.37%              62.15%              78.45%
    genewise   32407       46.47%              41.21%              64.66%              77.15%
    gth        9116        20.52%              64.70%              32.71%              90.87%
    all        32293       48.19%              42.89%              66.91%              78.00%
    filtered   23857       45.19%              54.44%              64.44%              83.72%

    --optimize_augustus_method <int>    default: 3
    设置AUGUSTUS Training时的参数优化方法。1，表示仅调用BGM2AT.optimize_augustus进行优化，能充分利用所有CPU线程对所有参数并行化测试，速度快；2，表示仅调用AUGUSTUS软件自带的optimize_augustus.pl程序进行优化，该方法的速度较慢，但结果更好；3，表示先使用BGM2AT.optimize_augustus优化完毕后，再使用AUGUSTUS软件自带的optimize_augustus.pl程序接着再进行优化，同时兼顾运算速度和效果。使用本参数的优先级更高，能覆盖--config指定参数配置文件中BGM2AT相同参数的值。
    
    --no_alternative_splicing_analysis    default: None
    添加该参数后，程序不会进行可变剪接分析。需要注意的时，当输入了NGS reads数据时，程序默认会根据intron和碱基测序深度信息进行基因的可变剪接分析。

    --delete_unimportant_intermediate_files    defaults: None
    添加该参数后，若程序运行成功，会删除不重要的中间文件，仅保留最少的、较小的、重要的中间结果文件。


在Rocky 9.2系统使用以下依赖的软件版本对本软件进行了测试并运行成功。

01. ParaFly (Version 0.1.0)
02. GNU parallel (Version 20230722)
03. RepeatMasker (version: 4.1.6)
04. RepeatModeler (version: 2.0.5)
05. makeblastdb/rmblastn/tblastn/blastp (Version: 2.14.0)
06. java (version: jdk-20.0.1)
07. hisat2 (version: 2.1.0)
08. samtools (version: 1.17)
09. mmseqs (version 15-6f452)
10. genewise (version: 2.4.1)
11. gth (Vesion: 1.7.3)
12. exonerate (Vesion: 2.2.0)
13. augustus/etraining (version: 3.5.0)
14. diamond (version 2.1.8)
15. hmmscan (version: 3.3.2)
16. busco (Version: 5.4.7)

Version of GETA: 2.7.1

USAGE

my $usage_english = &get_usage_english();
if (@ARGV==0){die $usage_english}

my ($genome, $RM_species, $RM_species_Dfam, $RM_species_RepBase, $RM_lib, $no_RepeatModeler, $pe1, $pe2, $single_end, $sam, $strand_specific, $protein, $augustus_species, $HMM_db, $BLASTP_db, $config, $BUSCO_lineage_dataset);
my ($out_prefix, $gene_prefix, $chinese_help, $help);
my ($cpu, $put_massive_temporary_data_into_memory, $genetic_code, $homolog_prediction_method, $optimize_augustus_method, $no_alternative_splicing_analysis, $delete_unimportant_intermediate_files);
my ($cmdString, $cmdString1, $cmdString2, $cmdString3, $cmdString4, $cmdString5, @cmdString);
my ($start_codon, $stop_codon);
GetOptions(
    "genome:s" => \$genome,
    "RM_species:s" => \$RM_species,
    "RM_species_Dfam:s" => \$RM_species_Dfam,
    "RM_species_RepBase:s" => \$RM_species_RepBase,
    "RM_lib:s" => \$RM_lib,
    "no_RepeatModeler!" => \$no_RepeatModeler,
    "pe1:s" => \$pe1,
    "pe2:s" => \$pe2,
    "se:s" => \$single_end,
    "sam:s" => \$sam,
    "strand_specific!" => \$strand_specific,
    "protein:s" => \$protein,
    "augustus_species:s" => \$augustus_species,
    "HMM_db:s" => \$HMM_db,
    "BLASTP_db:s" => \$BLASTP_db,
    "config:s" => \$config,
    "BUSCO_lineage_dataset:s" => \$BUSCO_lineage_dataset,
    "out_prefix:s" => \$out_prefix,
    "gene_prefix:s" => \$gene_prefix,
    "chinese_help!" => \$chinese_help,
    "help!" => \$help,
    "cpu:i" => \$cpu,
    "put_massive_temporary_data_into_memory!" => \$put_massive_temporary_data_into_memory,
    "genetic_code:i" => \$genetic_code,
    "homolog_prediction_method:s" => \$homolog_prediction_method,
    "optimize_augustus_method:i" => \$optimize_augustus_method,
    "no_alternative_splicing_analysis!" => \$no_alternative_splicing_analysis,
    "delete_unimportant_intermediate_files!" => \$delete_unimportant_intermediate_files,
);

if ( $chinese_help ) { die $usage_chinese }
elsif ( $help ) { die $usage_english }

# Step 0: 程序运行前的准备工作
# 0.1 对输入参数进行分析，使用绝对路径。若有参数设置不正确，则程序拒绝运行。
my (%HMM_db, %BLASTP_db, %config, %BUSCO_db, @BUSCO_db);
&parsing_input_parameters();

# 0.2 检测依赖的软件是否满足。
&detecting_dependent_softwares();

# 0.3 生成临时文件夹
my $tmp_dir = abs_path("$out_prefix.tmp");
mkdir $tmp_dir unless -e $tmp_dir;
chdir "$tmp_dir"; print STDERR "\nPWD: $tmp_dir\n";

# 0.4 准备基因组序列：
# 读取FASTA序列以>开始的头部时，去除第一个空及之后的字符，按从长到短排序。若序列名有重复，则修正序列名并保留其序列。
$cmdString = "$bin_path/genome_seq_clear.pl --no_rename --no_change_bp $genome > genome.fasta 2> genome.size";
unless (-s "genome.fasta") {
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
}
else {
    print STDERR "CMD(Skipped): $cmdString\n";
}
$genome = "$tmp_dir/genome.fasta";

# 获取基因组大小
my $genome_size = 0;
unless ( -s "genome.size" ) {
    open IN, $genome or die "Error: Can not open file $genome, $!";
    while ( <IN> ) {
        next if /^>/;
        $genome_size += (length($_) - 1);
    }
    close IN;

    open OUT, ">", "genome.size" or die "Error: Can not create file $tmp_dir/genome.size, $!";
    print OUT $genome_size;
    close OUT;
}
else {
    open IN, "$tmp_dir/genome.size" or die "Error: Can not open file $tmp_dir/genome.size, $!";
    $genome_size = <IN>;
    close IN;
}
# 根据基因组大小选择软件自带的配置文件，或使用 --config 参数指定的配置文件。
&choose_config_file($genome_size);

# 0.5 准备直系同源基因蛋白序列：
# 读取FASTA序列以>开始的头部时，去除第一个空及之后的字符，将所有怪异字符变为下划线字符；若遇到 > 符号不在句首，则表示fasta文件格式有误，删除该 > 及到下一个 > 之前的数据；去除序列中尾部的换行符, 将所有小写字符氨基酸变换为大写字符。程序中断后再次运行时，则重新读取新的同源蛋白序列文件，利用新的文件进行后续分析。
$cmdString = "$bin_path/fasta_format_revising.pl --seq_type protein --line_length 80 $protein > homolog.fasta 2> homolog.fasta.fasta_format_revising.log";
print STDERR (localtime) . ": CMD: $cmdString\n";
system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
$protein = "$tmp_dir/homolog.fasta";


# Step 1: RepeatMasker and RepeatModeler
print STDERR "\n============================================\n";
print STDERR "Step 1: RepeatMasker and RepeatModeler " . "(" . (localtime) . ")" . "\n";
mkdir "$tmp_dir/1.RepeatMasker" unless -e "$tmp_dir/1.RepeatMasker";

# 1.1 进行RepeatMasker Dfam分析
mkdir "$tmp_dir/1.RepeatMasker/repeatMasker_Dfam" unless -e "$tmp_dir/1.RepeatMasker/repeatMasker_Dfam";
chdir "$tmp_dir/1.RepeatMasker/repeatMasker_Dfam";
print STDERR "\nPWD: $tmp_dir/1.RepeatMasker/repeatMasker_Dfam\n";
if ( $RM_species_Dfam ) {
    $cmdString = "$bin_path/para_RepeatMasker $config{'para_RepeatMasker'} --species $RM_species_Dfam --engine hmmer --cpu $cpu --tmp_dir para_RepeatMasker.tmp $genome &> para_RepeatMasker.log";
}
else {
    $cmdString = "touch RepeatMasker_out.out";
}

&execute_cmds($cmdString, "$tmp_dir/1.RepeatMasker/repeatMasker_Dfam.ok");

# 1.2 进行RepeatMasker RepBase分析
mkdir "$tmp_dir/1.RepeatMasker/repeatMasker_RepBase" unless -e "$tmp_dir/1.RepeatMasker/repeatMasker_RepBase";
chdir "$tmp_dir/1.RepeatMasker/repeatMasker_RepBase";
print STDERR "\nPWD: $tmp_dir/1.RepeatMasker/repeatMasker_RepBase\n";
if ( $RM_species_RepBase ) {
    $cmdString = "$bin_path/para_RepeatMasker $config{'para_RepeatMasker'} --species $RM_species_RepBase --engine ncbi --cpu $cpu --tmp_dir para_RepeatMasker.tmp $genome &> para_RepeatMasker.log";
}
else {
    $cmdString = "touch RepeatMasker_out.out";
}

&execute_cmds($cmdString, "$tmp_dir/1.RepeatMasker/repeatMasker_RepBase.ok");

# 1.3 进行RepeatModeler分析
mkdir "$tmp_dir/1.RepeatMasker/repeatModeler" unless -e "$tmp_dir/1.RepeatMasker/repeatModeler";
chdir "$tmp_dir/1.RepeatMasker/repeatModeler";
print STDERR "\nPWD: $tmp_dir/1.RepeatMasker/repeatModeler\n";
if ( $RM_lib ) {
    $cmdString = "$bin_path/para_RepeatMasker $config{'para_RepeatMasker'} --out_prefix RepeatModeler_out --lib $RM_lib --cpu $cpu --tmp_dir para_RepeatMasker.tmp $genome &> para_RepeatMasker.log";
}
elsif ( $no_RepeatModeler ) {
    $cmdString = "touch RepeatModeler_out.out";
}
else {
    $cmdString = "BuildDatabase -name species -engine ncbi $genome";
    &execute_cmds($cmdString, "BuildDatabase.ok");

    my $cpu_RepeatModeler = $cpu;
    $cmdString = "RepeatModeler -threads $cpu_RepeatModeler -database species -LTRStruct &> RepeatModeler.log";
    # 若RepeatModeler版本为2.0.3或更低，则需要将-threads参数换为-pa参数。
    # 感谢Toney823在2023.07.27日提交的BUG反馈，https://github.com/chenlianfu/geta/issues/24，修改了下一行代码。
    my $RepeatModeler_version = `RepeatModeler --version`;
    if ( $RepeatModeler_version =~ m/(\d+)\.(\d+)\.(\d+)/ && ($1 < 2 or ( $1 == 2 && $2 == 0 && $3 <= 3)) )  {
        $cpu_RepeatModeler = int($cpu / 4);
        $cpu_RepeatModeler = 1 if $cpu_RepeatModeler < 1;
        $cmdString = "RepeatModeler -pa $cpu_RepeatModeler -database species -LTRStruct &> RepeatModeler.log";
    }

    &execute_cmds($cmdString, "RepeatModeler.ok");

    $cmdString = "$bin_path/para_RepeatMasker $config{'para_RepeatMasker'} --out_prefix RepeatModeler_out --lib RM_\*/\*.classified --cpu $cpu --tmp_dir para_RepeatMasker.tmp $genome &> para_RepeatMasker.log";
}

&execute_cmds($cmdString, "$tmp_dir/1.RepeatMasker/repeatModeler.ok");

# 1.4 合并RepeatMasker和RepeatModeler的结果
chdir "$tmp_dir/1.RepeatMasker/"; print STDERR "\nPWD: $tmp_dir/1.RepeatMasker\n";
if ( -s "repeatMasker_Dfam/RepeatMasker_out.out" or -s "repeatMasker_RepBase/RepeatMasker_out.out" or -s "repeatModeler/RepeatModeler_out.out" ) {
    $cmdString = "rm -rf genome.masked.fasta; $bin_path/merge_repeatMasker_out.pl $genome repeatMasker_Dfam/RepeatMasker_out.out repeatMasker_RepBase/RepeatMasker_out.out repeatModeler/RepeatModeler_out.out > genome.repeat.stats; $bin_path/maskedByGff.pl genome.repeat.gff3 $genome > genome.masked.fasta";
}
else{
    $cmdString = "ln -sf $genome genome.masked.fasta && touch genome.repeat.gff3 genome.repeat.stats";
}

&execute_cmds($cmdString, "$tmp_dir/1.RepeatMasker.ok");


# Step 2: homolog_prediction
print STDERR "\n============================================\n";
print STDERR "Step 2: Homolog_prediction" . "(" . (localtime) . ")" . "\n";
chdir $tmp_dir; print STDERR "\nPWD: $tmp_dir\n";
mkdir "$tmp_dir/2.homolog_prediction" unless -e "$tmp_dir/2.homolog_prediction";

if ( $protein ) {
    $cmdString = "$bin_path/homolog_prediction --tmp_dir $tmp_dir/2.homolog_prediction --cpu $cpu $config{'homolog_prediction'} --genetic_code $genetic_code --output_alignment_GFF3 $tmp_dir/2.homolog_prediction/homolog_alignment.gff3 --output_raw_GFF3 $tmp_dir/2.homolog_prediction/homolog_prediction.raw.gff3 $protein $genome > $tmp_dir/2.homolog_prediction/homolog_prediction.gff3 2> $tmp_dir/2.homolog_prediction/homolog_prediction.log";
}
else {
    $cmdString = "touch $tmp_dir/2.homolog_prediction/homolog_prediction.gff3";
}

&execute_cmds($cmdString, "$tmp_dir/2.homolog_prediction.ok");


# Step 3: NGSReads_predcition
print STDERR "\n============================================\n";
print STDERR "Step 3: NGSReads_predcition " . "(" . (localtime) . ")" . "\n";
chdir $tmp_dir; print STDERR "\nPWD: $tmp_dir\n";
mkdir "$tmp_dir/3.NGSReads_prediction" unless -e "$tmp_dir/3.NGSReads_prediction";

if (($pe1 && $pe2) or $single_end or $sam) {
    my @input_parameter;
    push @input_parameter, "--pe1 $pe1 --pe2 $pe2" if ($pe1 && $pe2);
    push @input_parameter, "--se $single_end" if $single_end;
    push @input_parameter, "--sam $sam" if $sam;
    push @input_parameter, "--strand_specific" if defined $strand_specific;
    push @input_parameter, "--genetic_code $genetic_code" if defined $genetic_code;
    push @input_parameter, "--put_massive_temporary_data_into_memory" if defined $put_massive_temporary_data_into_memory;
	push @input_parameter, "--homolog_gene_models $tmp_dir/2.homolog_prediction/homolog_prediction.raw.gff3" if defined $protein;
    my $input_parameter = join " ", @input_parameter;
    $cmdString = "$bin_path/NGSReads_prediction $input_parameter --config $tmp_dir/config.txt --cpu $cpu --tmp_dir $tmp_dir/3.NGSReads_prediction --output_alignment_GFF3 $tmp_dir/3.NGSReads_prediction/NGSReads_alignment.gff3 --output_raw_GFF3 $tmp_dir/3.NGSReads_prediction/NGSReads_prediction.raw.gff3 --intron_info_out $tmp_dir/3.NGSReads_prediction/intron.txt --base_depth_out $tmp_dir/3.NGSReads_prediction/base_depth.txt $genome > $tmp_dir/3.NGSReads_prediction/NGSReads_prediction.gff3 2> $tmp_dir/3.NGSReads_prediction/NGSReads_prediction.log";
}
else {
    $cmdString = "touch $tmp_dir/3.NGSReads_prediction/NGSReads_prediction.gff3 $tmp_dir/3.NGSReads_prediction/intron.txt; touch $tmp_dir/3.NGSReads_prediction/base_depth.txt; touch $tmp_dir/3.NGSReads_prediction/NGSReads_prediction.raw.gff3";
}

&execute_cmds($cmdString, "$tmp_dir/3.NGSReads_prediction.ok");


# Step 4: Augustus gene prediction
print STDERR "\n============================================\n";
print STDERR "Step 4: Augustus/HMM Trainning " . "(" . (localtime) . ")" . "\n";
mkdir "$tmp_dir/4.augustus" unless -e "$tmp_dir/4.augustus";
chdir "$tmp_dir/4.augustus";

# 4.1 Augustus HMM Training
mkdir "$tmp_dir/4.augustus/training" unless -e "$tmp_dir/4.augustus/training";
chdir "$tmp_dir/4.augustus/training"; print STDERR "\nPWD: $tmp_dir/4.augustus/training\n";

# 4.1.1 合并Transcript和Homolog预测的基因模型
$cmdString = "$bin_path/GFF3_merging_and_removing_redundancy $config{'GFF3_merging_and_removing_redundancy'} $genome $tmp_dir/3.NGSReads_prediction/NGSReads_prediction.raw.gff3 $tmp_dir/2.homolog_prediction/homolog_prediction.raw.gff3 > evidence_gene_models.gff3 2> GFF3_merging_and_removing_redundancy.log";

&execute_cmds($cmdString, "01.evidence_gene_models.ok");

# 4.1.2 选择完整且准确的基因模型
unless ( -s "excellent.gff3" ) {
	my $input = "$tmp_dir/4.augustus/training/evidence_gene_models.gff3";
	open IN, $input or die "Error: Can not open file $input, $!";
	my $output = "$tmp_dir/4.augustus/training/evidence_gene_models.excellent.gff3";
	open OUT, ">", $output or die "Error: Can not create file $output, $!";
	$/ = "\n\n";
	while (<IN>) {
		print OUT if m/excellent/;
	}
	close IN; close OUT;
	open OUT, ">", "02.evidence_gene_models.excellent.ok" or die $!; close OUT;
}

# 4.1.3 选择CDS数量较多、CDS长度较长、CDS/exon比例较大且去冗余的基因模型。
$cmdString = "$bin_path/geneModels2AugusutsTrainingInput $config{'geneModels2AugusutsTrainingInput'} --out_prefix ati --cpu $cpu evidence_gene_models.excellent.gff3 $genome &> geneModels2AugusutsTrainingInput.log";
unless ( -e "03.geneModels2AugusutsTrainingInput.ok" ) {
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

    # 若用于Augustus training的基因数量少于1000个，则重新运行geneModels2AugusutsTrainingInput，降低阈值来增加基因数量。
    my $input_file = "$tmp_dir/4.augustus/training/geneModels2AugusutsTrainingInput.log";
    open IN, $input_file or die "Error: Can not open file $input_file, $!";
    my $training_genes_number = 0;
    while (<IN>) {
        $training_genes_number = $1 if m/Best gene Models number:\s+(\d+)/;
    }
    close IN;
    if ( $training_genes_number < 1000 ) {
        $cmdString = "$bin_path/geneModels2AugusutsTrainingInput --min_evalue 1e-9 --min_identity 0.9 --min_coverage_ratio 0.9 --min_cds_num 1 --min_cds_length 450 --min_cds_exon_ratio 0.40 --keep_ratio_for_excluding_too_long_gene 0.99 --out_prefix ati --cpu $cpu $cpu $tmp_dir/4.evidence_gene_models/evidence_gene_models.gff3 $genome &> geneModels2AugusutsTrainingInput.log.Loose_thresholds";
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    }

    open OUT, ">", "03.geneModels2AugusutsTrainingInput.ok" or die $!; close OUT;
}
else {
    print STDERR "CMD(Skipped): $cmdString\n";
}

# 4.1.4 分析基因长度和基因间区长度信息，从而确定Trainig时选择基因模型两侧翼序列长度。
my $flanking_length;
unless ( -e "04.get_flanking_length.ok" && -s "$tmp_dir/4.augustus/training/flanking_length.txt" ) {
    my (%gene_info, @intergenic_length, @gene_length);
    # 读取基因模型信息
    my $input_file = "$tmp_dir/4.augustus/training/evidence_gene_models.gff3";
    open IN, $input_file or die "Error: Can not open file $input_file, $!";
    while (<IN>) {
        if (m/\tgene\t/) {
            @_ = split /\t/;
            $gene_info{$_[0]}{$_[6]}{"$_[3]\t$_[4]"} = 1;
        }
    }
    close IN;
    # 分析基因长度和基因间区长度
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
                push @intergenic_length, $distance if $distance >= 50;
                $end = $bb if $end < $bb;
            }
        }
    }
    # 设置flanking_length长度为(基因间区长度中位数的八分之一，基因长度中位数)的较小值。
    @gene_length = sort {$a <=> $b} @gene_length;
    @intergenic_length = sort {$a <=> $b} @intergenic_length;
    $flanking_length = int($intergenic_length[@intergenic_length/2] / 8);
    $flanking_length = $gene_length[@gene_length/2] if $flanking_length >= $gene_length[@gene_length/2];
    # 输出flanking_length到指定文件中。
    my $outpu_file = "$tmp_dir/4.augustus/training/flanking_length.txt";
    open OUT, ">", $outpu_file or die "Can not create file $outpu_file, $!";
    print OUT $flanking_length;
    close OUT;

    open OUT, ">", "04.get_flanking_length.ok" or die $!; close OUT;
    print STDERR (localtime) . ": The flanking length was set to $flanking_length for AUGUSTUS Training.\n";
}
else {
    my $input_file = "$tmp_dir/4.augustus/training/flanking_length.txt";
    open IN, $input_file or die "Error: Can not open file $input_file, $!";
    $flanking_length = <IN>;
    close IN;
    print STDERR (localtime) . ": The flanking length was set to $flanking_length for AUGUSTUS Training.\n";
}

# 4.1.5 进行Augustus training
# 若--augustus_species参数设置的文件夹已经存在，则对其进行参数优化；不存在，则完全重新进行Training和优化。
my $augustus_species_start_from = "";
if ( -e "$ENV{'AUGUSTUS_CONFIG_PATH'}/species/$augustus_species/${augustus_species}_parameters.cfg" ) {
    $augustus_species_start_from = "--augustus_species_start_from $augustus_species";
}
$cmdString1 = "$bin_path/BGM2AT $config{'BGM2AT'} --AUGUSTUS_CONFIG_PATH $tmp_dir/4.augustus/config $augustus_species_start_from --flanking_length $flanking_length --CPU $cpu --onlytrain_GFF3 ati.filter1.gff3 ati.filter2.gff3 $genome $augustus_species &> BGM2AT.log";
$cmdString2 = "cp accuary_of_AUGUSTUS_HMM_Training.txt ../";
$cmdString3 = "cp -a $tmp_dir/4.augustus/config/species/$augustus_species $ENV{'AUGUSTUS_CONFIG_PATH'}/species/ || echo Can not create file in directory $ENV{'AUGUSTUS_CONFIG_PATH'}/species/";

&execute_cmds($cmdString1, $cmdString2, $cmdString3, "$tmp_dir/4.augustus/training.ok");

# 4.2 准备Hints信息
chdir "$tmp_dir/4.augustus"; print STDERR "\nPWD: $tmp_dir/4.augustus\n";
$cmdString = "$bin_path/prepareAugusutusHints $config{'prepareAugusutusHints'} $tmp_dir/4.augustus/evidence_gene_models.gff3 > hints.gff";

&execute_cmds($cmdString, "prepareAugusutusHints.ok");

# 4.3 进行Augustus基因预测
# 4.3.1 先分析基因长度信息，从而计算基因组序列打断的长度，以利于并行化运行augustus命令。
my ($segmentSize, $overlapSize) = (1000000, 100000);
unless ( -e "get_segmentSize.ok" ) {
    # 获取最长的基因长度
    my $input_file = "$tmp_dir/4.augustus/training/evidence_gene_models.gff3";
    open IN, $input_file or die "Error: Can not open file $input_file, $!";
    my @gene_length;
    while (<IN>) {
        if (m/\tgene\t(\d+)\t(\d+)\t/) {
            push @gene_length, $2 - $1 + 1;
        }
    }
    @gene_length = sort {$a <=> $b} @gene_length;
    # 打断序列时，要求相邻序列重叠区域至少为最长基因长度的2倍，并取整（保留前两个数字四舍五入，后面的全为0）；片段长度为重叠序列长度10倍。
    if ($gene_length[-1] * 4 > $overlapSize) {
        $overlapSize = $gene_length[-1] * 2;
        my $overlapSize_length = length($overlapSize);
        $overlapSize_length --;
        $overlapSize_length --;
        $overlapSize = int(($overlapSize / (10 ** $overlapSize_length)) + 1) * (10 ** $overlapSize_length);
        $segmentSize = $overlapSize * 10;
    }
    # 将overlapSize和segmentSize输出到文件segmentSize.txt文件中
    my $outpu_file = "$tmp_dir/4.augustus/segmentSize.txt";
    open OUT, ">", $outpu_file or die "Can not create file $outpu_file, $!";
    print OUT "$segmentSize\t$overlapSize";
    close OUT;

    open OUT, ">", "get_segmentSize.ok" or die $!; close OUT;
    print STDERR (localtime) . ": When executing augustus command lines using ParaFly, the genome sequences were split into segments. The longest segment spanned $segmentSize bp, and two adjacent segments overlapped by $overlapSize bp.\n";
}
else {
    my $input_file = "$tmp_dir/4.augustus/segmentSize.txt";
    open IN, $input_file or die "Error: Can not open file $input_file, $!";
    ($segmentSize, $overlapSize) = split /\t/, <IN>;
    close IN;
    print STDERR (localtime) . ": When executing augustus command lines using ParaFly, the genome sequences were split into segments. The longest segment spanned $segmentSize bp, and two adjacent segments overlapped by $overlapSize bp.\n";
}

# 4.3.2 Augustus gene prediction
$cmdString1 = "$bin_path/paraAugusutusWithHints $config{'paraAugusutusWithHints'} --species $augustus_species --AUGUSTUS_CONFIG_PATH $tmp_dir/4.augustus/config --cpu $cpu --segmentSize $segmentSize --overlapSize $overlapSize --tmp_dir aug_para_with_hints $genome hints.gff > augustus.raw.gff3";
$cmdString2 = "$bin_path/addHintRatioToAugustusResult $tmp_dir/4.augustus/training/evidence_gene_models.gff3 hints.gff augustus.raw.gff3 > augustus.gff3";

&execute_cmds($cmdString1, $cmdString2, "$tmp_dir/4.augustus.ok");


# Step 5: CombineGeneModels
print STDERR "\n============================================\n";
print STDERR "Step 6: CombineGeneModels " . "(" . (localtime) . ")" . "\n";
mkdir "$tmp_dir/6.combine_gene_models" unless -e "$tmp_dir/6.combine_gene_models";
chdir "$tmp_dir/6.combine_gene_models"; print STDERR "\nPWD: $tmp_dir/6.combine_gene_models\n";

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

# 6.1 使用AUGUSTUS预测的基因模型对NGSReads和Homolog预测的不完整基因模型进行填补。
@cmdString = ();
push @cmdString, "$bin_path/GFF3_filling_gene_models_Parallel --cpu $cpu --tmp_dir a.FillingGeneModelsByAugustus_for_homolog_prediction --ouput_filling_detail_tab a.FillingGeneModelsByAugustus_for_homolog_prediction.tab --start_codon $start_codon --stop_codon $stop_codon --attribute_for_filling_complete Filled_by_AUGUSTUS=True $genome $tmp_dir/2.homolog_prediction/homolog_prediction.raw.gff3 $tmp_dir/4.augustus/augustus.gff3 > homolog_prediction.FillByAugustus.gff3 2> homolog_prediction.FillByAugustus.log";
push @cmdString, "$bin_path/GFF3_filling_gene_models_Parallel --cpu $cpu --tmp_dir a.FillingGeneModelsByAugustus_for_NGSReads_prediction --ouput_filling_detail_tab a.FillingGeneModelsByAugustus_for_NGSReads_prediction.tab --start_codon $start_codon --stop_codon $stop_codon --attribute_for_filling_complete Filled_by_AUGUSTUS=True $genome $tmp_dir/3.NGSReads_prediction/NGSReads_prediction.raw.gff3 $tmp_dir/4.augustus/augustus.gff3 > NGSReads_prediction.FillByAugustus.gff3 2> NGSReads_prediction.FillByAugustus.log";

&execute_cmds(@cmdString, "01.FillingGeneModelsByAugustus.ok");

# 6.2 合并NGSReads、Homolog和Augustus预测的基因模型
@cmdString = ();
push @cmdString, "$bin_path/GFF3_merging_and_removing_redundancy $config{'GFF3_merging_and_removing_redundancy'} --start_codon $start_codon --stop_codon $stop_codon $genome NGSReads_prediction.FillByAugustus.gff3 homolog_prediction.FillByAugustus.gff3 > evidence_prediction.gff3 2> GFF3_merging_and_removing_redundancy.1.log";
push @cmdString, "$bin_path/GFF3_merging_and_removing_redundancy $config{'GFF3_merging_and_removing_redundancy'} --start_codon $start_codon --stop_codon $stop_codon $genome evidence_prediction.gff3 $tmp_dir/4.augustus/augustus.gff3 > all_prediction.gff3 2> GFF3_merging_and_removing_redundancy.2.log";

&execute_cmds(@cmdString, "02.GFF3_merging_and_removing_redundancy.ok");

# 6.3 对不完整基因模型进行首尾强制补齐。
$cmdString = "$bin_path/fillingEndsOfGeneModels $config{'fillingEndsOfGeneModels'} --filling_need_transcriptID filling_need_transcriptID.txt --nonCompletedGeneModels geneModels.f.gff3 $genome all_prediction.gff3 > all_prediction.filled.gff3 2> fillingEndsOfGeneModels.1.log";

&execute_cmds($cmdString, "03.fillingEndsOfGeneModels.ok");

exit;
# 6.1 第一轮基因预测结果整合：以AUGUSTUS结果为主，进行三种基因预测结果的整合
# 对三种基因预测结果进行第一轮整合，以Augustus结果为准。得到 combine.1.gff3 为有Evidence支持的结果，combine.2.gff3为支持不足的结果。
my $cmdString1 = "$bin_path/paraCombineGeneModels $config{'paraCombineGeneModels'} --cpu $cpu $tmp_dir/4.augustus/augustus.gff3 $tmp_dir/3.NGSReads_prediction/NGSReads_prediction.raw.gff3 $tmp_dir/2.homolog_prediction/homolog_prediction.gff3 $tmp_dir/4.augustus/hints.gff &> /dev/null";
my $cmdString2 = "$bin_path/GFF3Clear --genome $genome --no_attr_add --coverage 0.8 combine.1.gff3 > geneModels.a.gff3 2> /dev/null";
my $cmdString3 = "$bin_path/GFF3Clear --genome $genome --no_attr_add --coverage 0.8 combine.2.gff3 > geneModels.b.gff3 2> /dev/null";
my $cmdString4 = "perl -p -i -e 's/(=[^;]+)\.t1/\$1.t01/g' geneModels.a.gff3 geneModels.b.gff3";

&execute_cmds($cmdString1, $cmdString2, $cmdString3, $cmdString4, "01.paraCombineGeneModels.ok");

# 6.2 第二轮基因预测结果整合：以转录本和同源蛋白预测结果为准，对上一步的基因模型进行优化。
my $cmdString1 = "perl -p -e 's/(=[^;]+)\.t1/\$1.t01/g;' $tmp_dir/4.evidence_gene_models/evidence_gene_models.gff3 > geneModels.c.gff3;";
my $cmdString2 = "$bin_path/pickout_better_geneModels_from_evidence $config{'pickout_better_geneModels_from_evidence'} geneModels.a.gff3 geneModels.c.gff3 > picked_evidence_geneModels.gff3 2> picked_evidence_geneModels.log";
my $cmdString3 = "$bin_path/GFF3Clear --genome $genome --no_attr_add picked_evidence_geneModels.gff3 geneModels.a.gff3 > geneModels.d.gff3 2> /dev/null";
my $cmdString4 = "perl -p -i -e 's/Integrity=[^;]+;?//g' geneModels.d.gff3";

&execute_cmds($cmdString1, $cmdString2, $cmdString3, $cmdString4, "02.pickout_better_geneModels_from_evidence.ok");

# 6.3 对不完整基因模型进行首尾补齐。
$cmdString = "$bin_path/fillingEndsOfGeneModels $config{'fillingEndsOfGeneModels'} --filling_need_transcriptID filling_need_transcriptID.txt --nonCompletedGeneModels geneModels.f.gff3 $genome geneModels.d.gff3 > geneModels.e.gff3 2> fillingEndsOfGeneModels.1.log";

&execute_cmds($cmdString, "03.fillingEndsOfGeneModels.ok");

# 6.4 分别对对基因模型 geneModels.b.gff3, geneModels.e.gff3 and geneModels.f.gff3 进行可变剪接分析
unless ( $no_alternative_splicing_analysis ) {
    $cmdString1 = "$bin_path/paraAlternative_splicing_analysis $config{'alternative_splicing_analysis'} --tmp_dir paraAlternative_splicing_analysis.gb.tmp --cpu $cpu geneModels.b.gff3 $tmp_dir/3.NGSReads_prediction/intron.txt $tmp_dir/3.NGSReads_prediction/base_depth.txt > geneModels.gb_AS.gff3 2> alternative_splicing.gb.stats && $bin_path/GFF3_add_CDS_for_transcript $genome geneModels.gb_AS.gff3 > geneModels.gb.gff3";
    $cmdString2 = "$bin_path/paraAlternative_splicing_analysis $config{'alternative_splicing_analysis'} --tmp_dir paraAlternative_splicing_analysis.ge.tmp --cpu $cpu geneModels.e.gff3 $tmp_dir/3.NGSReads_prediction/intron.txt $tmp_dir/3.NGSReads_prediction/base_depth.txt > geneModels.ge_AS.gff3 2> alternative_splicing.ge.stats && $bin_path/GFF3_add_CDS_for_transcript $genome geneModels.ge_AS.gff3 > geneModels.ge.gff3";
    $cmdString3 = "$bin_path/paraAlternative_splicing_analysis $config{'alternative_splicing_analysis'} --tmp_dir paraAlternative_splicing_analysis.gf.tmp --cpu $cpu geneModels.f.gff3 $tmp_dir/3.NGSReads_prediction/intron.txt $tmp_dir/3.NGSReads_prediction/base_depth.txt > geneModels.gf_AS.gff3 2> alternative_splicing.gf.stats && $bin_path/GFF3_add_CDS_for_transcript $genome geneModels.gf_AS.gff3 > geneModels.gf.gff3";
}
else {
    $cmdString1 = "$bin_path/GFF3_add_CDS_for_transcript $genome geneModels.b.gff3 > geneModels.gb.gff3";
    $cmdString2 = "$bin_path/GFF3_add_CDS_for_transcript $genome geneModels.e.gff3 > geneModels.ge.gff3";
    $cmdString3 = "$bin_path/GFF3_add_CDS_for_transcript $genome geneModels.f.gff3 > geneModels.gf.gff3";
}

&execute_cmds($cmdString1, $cmdString2, $cmdString3, "04.alternative_splicing_analysis.ok");

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
$cmdString1 = "$bin_path/GFF3_extract_TranscriptID_for_filtering $config{'GFF3_extract_TranscriptID_for_filtering'} $tmp_dir/1.RepeatMasker/genome.repeat.gff3 geneModels.gb.gff3 geneModels.ge.gff3 geneModels.gf.gff3 > transcriptID_for_filtering.txt";
# 提取没有足够证据基因的所有转录本ID
$cmdString2 = "perl -ne 'print \"\$1\\tNotEnoughEvidence\\n\" if m/ID=([^;]*\\.t\\d+);/;' geneModels.gb.gff3 >> transcriptID_for_filtering.txt";
# 提取没法填补完整基因的所有转录本ID
$cmdString3 = "perl -ne 'print \"\$1\\tFilling2Uncomplete\\n\" if m/ID=([^;]*\\.t\\d+);/;' geneModels.gf.gff3 >> transcriptID_for_filtering.txt";
# 提取通过填补而完整的基因的所有转录本ID
$cmdString4 = "perl -e 'open IN, \"filling_need_transcriptID.txt\"; while (<IN>) { s/.t\\d+\\n//; \$gene{\$_} = 1; } while (<>) { print \"\$1\\tFilling2Complete\\n\" if m/ID=(([^;]+)\\.t\\d+);/ && exists \$gene{\$2} }' geneModels.ge_AS.gff3 >> transcriptID_for_filtering.txt";

&execute_cmds($cmdString1, $cmdString2, $cmdString3, $cmdString4, "05.extract_TranscriptID_for_filtering.ok");
    
# 6.6 提取待过滤转录本的蛋白序列。
my $cmdString1 = "$bin_path/gff3_to_protein.pl $genome geneModels.gb.gff3 geneModels.gf.gff3 geneModels.ge.gff3 > proteins_all.fasta 2> gff3_to_protein.log";
my $cmdString2 = "perl -p -i -e 's/\\*\$//' proteins_all.fasta";
my $cmdString3 = "$bin_path/fasta_extract_subseqs_from_list.pl proteins_all.fasta transcriptID_for_filtering.txt > proteins_for_filtering.fasta 2> fasta_extract_subseqs_from_list.log";

&execute_cmds($cmdString1, $cmdString2, $cmdString3, "06.extract_proteins_for_filtering.ok" );

# 6.7 对蛋白序列进行HMM和BLASTP验证。
$cmdString1 = "";
if ( $HMM_db ) {
    my $hmmscan_cpu = 0;
    $hmmscan_cpu = $1 if $config{'para_hmmscan'} =~ m/--hmmscan_cpu\s+(\d+)/;
    my $para_hmmscan_cpu = $cpu;
    $para_hmmscan_cpu = int($cpu / $hmmscan_cpu + 0.5) if $hmmscan_cpu;
    $para_hmmscan_cpu = 1 if $para_hmmscan_cpu < 1;

    foreach ( sort keys %HMM_db ) {
        $cmdString1 .= "$bin_path/para_hmmscan $config{'para_hmmscan'} --outformat --cpu $para_hmmscan_cpu --no_cut_ga --hmm_db $_ --tmp_prefix $HMM_db{$_} proteins_for_filtering.fasta >> validation_hmmscan.tab 2>> para_hmmscan.1.log; $bin_path/para_hmmscan $config{'para_hmmscan'} --chunk 1 --outformat --cpu $para_hmmscan_cpu --no_cut_ga --hmm_db $_ --tmp_prefix $HMM_db{$_} proteins_for_filtering.fasta >> validation_hmmscan.tab 2>> para_hmmscan.2.log; ";
    }
}
$cmdString2 = "";
$cmdString3 = "";
if ( $BLASTP_db ) {
    foreach ( sort keys %BLASTP_db ) {
        $cmdString2 .= "diamond blastp $config{'diamond'} --outfmt 5 --db $_ --query proteins_for_filtering.fasta --out validation_blastp_$BLASTP_db{$_}.xml --threads $cpu &>> diamond_blastp.log; ";
        $cmdString3 = "$bin_path/parsing_blast_result.pl $config{'parsing_blast_result.pl'} --out-hit-confidence validation_blastp_$BLASTP_db{$_}.xml >> validation_blastp.tab; ";
    }
}
$cmdString2 .= "diamond makedb --db homolog --in $tmp_dir/homolog.fasta &> diamond_makedb.log; diamond blastp $config{'diamond'} --outfmt 5 --db homolog --query proteins_for_filtering.fasta --out validation_blastp.xml --threads $cpu &> diamond_blastp.log";
$cmdString3 .= "$bin_path/parsing_blast_result.pl $config{'parsing_blast_result.pl'} --out-hit-confidence validation_blastp.xml >> validation_blastp.tab";

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

# 6.8 根据HMM和BLASTP验证结果对基因模型进行过滤。
# 获得验证通过的转录本ID
$cmdString1 = "$bin_path/get_valid_transcriptID $config{'get_valid_transcriptID'} validation_hmmscan.tab validation_blastp.tab > transcriptID_validating_passed.tab 2> get_valid_transcriptID.log";
$cmdString2 = "$bin_path/get_valid_geneModels $config{'get_valid_geneModels'} --out_prefix geneModels.h transcriptID_for_filtering.txt transcriptID_validating_passed.tab geneModels.gb.gff3 geneModels.ge.gff3 geneModels.gf.gff3 2> get_valid_geneModels.log";

&execute_cmds($cmdString1, $cmdString2, "08.filtering_geneModels.ok");

# 6.9 再次对基因模型进行首尾填补。
$cmdString = "$bin_path/fillingEndsOfGeneModels $config{'fillingEndsOfGeneModels'} $genome geneModels.h.coding.gff3 > geneModels.i.coding.gff3 2> fillingEndsOfGeneModels.2.log";

&execute_cmds($cmdString, "$tmp_dir/6.combine_gene_models.ok");


# Step 7: OutPut
print STDERR "\n============================================\n";
print STDERR "Step 7: Output gene models " . "(" . (localtime) . ")" . "\n";
mkdir "$tmp_dir/7.output_gene_models" unless -e "$tmp_dir/7.output_gene_models";
chdir "$tmp_dir/7.output_gene_models"; print STDERR "\nPWD: $tmp_dir/7.output_gene_models\n";

# 7.1 输出GFF3格式文件基因结构注释信息
$cmdString1 = "$bin_path/GFF3Clear --GFF3_source GETA --gene_prefix $gene_prefix --gene_code_length 6 --genome $genome $tmp_dir/6.combine_gene_models/geneModels.i.coding.gff3 > $out_prefix.geneModels.gff3 2> /dev/null";
$cmdString2 = "$bin_path/GFF3_extract_bestGeneModels $out_prefix.geneModels.gff3 > $out_prefix.bestGeneModels.gff3 2> $out_prefix.AS_num_of_codingTranscripts.stats";
$cmdString3 = "$bin_path/GFF3Clear --GFF3_source GETA --gene_prefix ${out_prefix}ncGene  --gene_code_length 6 --genome $genome --no_attr_add $tmp_dir/6.combine_gene_models/geneModels.h.lncRNA.gff3 > $out_prefix.geneModels_lncRNA.gff3 2> /dev/null";
$cmdString4 = "$bin_path/GFF3Clear --GFF3_source GETA --gene_prefix ${out_prefix}lqGene --gene_code_length 6 --genome $genome --no_attr_add --coverage 0.6 $tmp_dir/6.combine_gene_models/geneModels.h.lowQuality.gff3 > $out_prefix.geneModels_lowQuality.gff3 2> /dev/null";

&execute_cmds($cmdString1, $cmdString2, $cmdString3, $cmdString4, "1.output_GFF3.ok");

# 7.2 输出GTF文件和基因的序列信息
$cmdString1 = "$bin_path/gff3ToGtf.pl $genome $out_prefix.geneModels.gff3 > $out_prefix.geneModels.gtf 2> /dev/null";
$cmdString2 = "$bin_path/gff3ToGtf.pl $genome $out_prefix.bestGeneModels.gff3 > $out_prefix.bestGeneModels.gtf 2> /dev/null";
#$cmdString3 = "$bin_path/eukaryotic_gene_model_statistics.pl $out_prefix.bestGeneModels.gtf $genome $out_prefix &> $out_prefix.geneModels.stats";
$cmdString3 = "$bin_path/gff3_to_sequences.pl --out_prefix $out_prefix --only_gene_sequences --only_coding_gene_sequences --only_first_isoform --genetic_code 1 $genome $out_prefix.geneModels.gff3 > $out_prefix.geneModels.stats 2> /dev/null";
&execute_cmds($cmdString1, $cmdString2, $cmdString3, "2.output_GTF.ok");

# 7.3 输出重复序列信息及其统计结果、RepeatModeler软件构建的重复序列数据库和masked genome sequence
if ( $RM_species or $RM_species_Dfam or $RM_species_RepBase or $RM_lib or (! $no_RepeatModeler) ) {
    $cmdString1 = "cp $tmp_dir/1.RepeatMasker/genome.masked.fasta $out_prefix.maskedGenome.fasta";
    $cmdString2 = "cp $tmp_dir/1.RepeatMasker/genome.repeat.stats $out_prefix.repeat.stats";
    $cmdString3 = "cp $tmp_dir/1.RepeatMasker/genome.repeat.gff3 $out_prefix.repeat.gff3";

    $cmdString4 = "echo NO RepeatModeler Result !";
    if ( $RM_lib ) {
        $cmdString4 = "ln -sf $RM_lib $out_prefix.repeat.lib";
    }
    elsif ( ! $no_RepeatModeler ) {
        $cmdString4 = "cp $tmp_dir/1.RepeatMasker/repeatModeler/RM_\*/\*.classified $out_prefix.repeat.lib";
    }

    &execute_cmds($cmdString1, $cmdString2, $cmdString3, $cmdString4, "3.output_repeat.ok");
}
else {
    unless ( -e "3.output_repeat.ok" ) {
        open OUT, ">", "3.output_repeat.ok" or die $!; close OUT;
    }
}

# 7.4 输出转录本、同源蛋白和Augustus的基因预测结果
my @cmdString;
if ( ($pe1 && $pe2) or $single_end or $sam ) {
    push @cmdString, "cp $tmp_dir/3.NGSReads_prediction/NGSReads_alignment.gff3 $out_prefix.NGSReads_alignment.gff3";
    push @cmdString, "cp $tmp_dir/3.NGSReads_prediction/NGSReads_prediction.gff3 $out_prefix.NGSReads_prediction.gff3";
}
if ( $protein ) {
    push @cmdString, "cp $tmp_dir/2.homolog_prediction/homolog_alignment.gff3 $out_prefix.homolog_alignment.gff3";
    push @cmdString, "cp $tmp_dir/2.homolog_prediction/homolog_prediction.gff3 $out_prefix.homolog_prediction.gff3";
}
if ( (($pe1 && $pe2) or $single_end or $sam) && $protein ) {
    push @cmdString, "cp $tmp_dir/4.evidence_gene_models/evidence_gene_models.gff3 $out_prefix.evidence_prediction.gff3";
}
push @cmdString, "$bin_path/GFF3Clear --genome $genome --no_attr_add $tmp_dir/4.augustus/augustus.gff3 > $out_prefix.augustus_prediction.gff3";
push @cmdString, "4.output_methods_GFF3.ok";

&execute_cmds(@cmdString);

# 7.5 进行BUSCO分析
my @cmdString;
if ( $BUSCO_lineage_dataset ) {
    foreach my $BUSCO_db_path ( @BUSCO_db ) {
        my $BUSCO_db_name = $BUSCO_db{$BUSCO_db_path};

        my $predictd_type = "GETA";
        push @cmdString, "rm -rf BUSCO_OUT.$predictd_type.$BUSCO_db_name*; busco -i $out_prefix.protein.fasta -o BUSCO_OUT.$predictd_type.$BUSCO_db_name -m protein -l $BUSCO_db_path -c $cpu --offline &> BUSCO_OUT.$predictd_type.$BUSCO_db_name.log";

        $predictd_type = "AUGUSTUS";
        push @cmdString, "rm -rf BUSCO_OUT.$predictd_type.$BUSCO_db_name*; $bin_path/gff3_to_protein.pl $genome $out_prefix.augustus_prediction.gff3 > protein.augustus_prediction.fasta; busco -i protein.augustus_prediction.fasta -o BUSCO_OUT.$predictd_type.$BUSCO_db_name -m protein -l $BUSCO_db_path -c $cpu --offline &> BUSCO_OUT.$predictd_type.$BUSCO_db_name.log";

        if ( $protein ) {
            $predictd_type = "homolog";
            push @cmdString, "rm -rf BUSCO_OUT.$predictd_type.$BUSCO_db_name*; $bin_path/gff3_to_protein.pl $genome $out_prefix.homolog_prediction.gff3 > protein.homolog_prediction.fasta; busco -i protein.homolog_prediction.fasta -o BUSCO_OUT.$predictd_type.$BUSCO_db_name -m protein -l $BUSCO_db_path -c $cpu --offline &> BUSCO_OUT.$predictd_type.$BUSCO_db_name.log";
        }

        if ( ($pe1 && $pe2) or $single_end or $sam ) {
            $predictd_type = "NGSReads";
            push @cmdString, "rm -rf BUSCO_OUT.$predictd_type.$BUSCO_db_name*; $bin_path/gff3_to_protein.pl $genome $out_prefix.NGSReads_prediction.gff3 > protein.NGSReads_prediction.fasta; busco -i protein.NGSReads_prediction.fasta -o BUSCO_OUT.$predictd_type.$BUSCO_db_name -m protein -l $BUSCO_db_path -c $cpu --offline &> BUSCO_OUT.$predictd_type.$BUSCO_db_name.log";
        }
    }
    push @cmdString, "5.BUSCO.ok";

    &execute_cmds(@cmdString);

    # 合并多个数据库的BUSCO结果
    unless ( -e "BUSCO_results.txt" ) {
        my $BUSCO_results_OUT;

        my @predictd_type = ("GETA", "AUGUSTUS");
        push @predictd_type, "homolog" if $protein;
        push @predictd_type, "NGSReads" if ( ($pe1 && $pe2) or $single_end or $sam );
        foreach my $predictd_type ( @predictd_type ) {
            $BUSCO_results_OUT .= "The BUSCO result of whole genome proteins predicted by $predictd_type:\n";
            foreach my $BUSCO_db_path ( @BUSCO_db ) {
                my $BUSCO_db_name = $BUSCO_db{$BUSCO_db_path};
                my $input_file = "BUSCO_OUT.$predictd_type.$BUSCO_db_name.log";
                open IN, $input_file or die "Error: Can not open file $input_file, $!";
                while ( <IN> ) {
                    $BUSCO_results_OUT .= sprintf("%-25s$1\n", $BUSCO_db_name) if m/(C:\S+)/;
                }
                close IN;
            }
            $BUSCO_results_OUT .= "\n";
        }

        open OUT, ">", "BUSCO_results.txt" or die "Error: Can not create file BUSCO_results.txt, $!";
        print OUT $BUSCO_results_OUT;
        close OUT;
    }
}

# 7.6 输出GETA基因基因预测的各项流程统计信息，用于追踪基因预测结果的可靠性
unless ( -e "$tmp_dir/7.output_gene_models.ok" ) {
    open OUT, ">", "$out_prefix.gene_prediction.summary" or die "Can not create file $out_prefix.gene_prediction.summary, $!";
    # (1) 获取基因组重复序列统计信息
    if ( -e "$out_prefix.repeat.stats" ) {
        my $input_file = "$out_prefix.repeat.stats";
        open IN, $input_file or die "Error: Can not open file $input_file, $!";
        while (<IN>) {
            print OUT $_ if m/^Genome Size/;
            print OUT "$_\n" if m/^Repeat Ratio/;
        }
        close IN;
    }
    # (2) 获取转录组二代测序数据和基因组的匹配率
    if ( -e "$tmp_dir/3.NGSReads_prediction/b.hisat2/hisat2.log" ) {
        my $input_file = "$tmp_dir/3.NGSReads_prediction/b.hisat2/hisat2.log";
        open IN, $input_file or die "Error: Can not open file $input_file, $!";
        while (<IN>) {
            print OUT "The alignment rate of RNA-Seq reads is: $1\n\n" if m/^(\S+) overall alignment rate/;
        }
        close IN;
    # (3) 获取转录组数据预测的基因数量的统计信息
        my $input_file = "$tmp_dir/3.NGSReads_prediction/NGSReads_prediction.A.log";
        print OUT "The statistics of gene models predicted by NGS Reads:\n";
        print OUT `tail -n 1 $input_file`;
        print OUT "\n";
    }
    # (4) 获取同源蛋白预测的基因数量的统计信息
    if ( $protein ) {
        my $input_file = "$tmp_dir/2.homolog_prediction/homolog_prediction.log";
        print OUT "The statistics of gene models predicted by homolog:\n";
        print OUT `tail -n 3 $input_file`;
    }
    # (5) 获取AUGUSTUS基因预测数量
    if (-e "$tmp_dir/4.augustus/augustus.gff3") {
        my $input_file = "$tmp_dir/4.augustus/augustus.gff3";
        print OUT "The statistics of gene models predicted by AUGUSTUS:\n";
        my $gene_num = `grep -P "\tgene" $input_file | wc -l`; chomp($gene_num);
        print OUT "$gene_num genes were predicted by augustus.\n\n";
    }
    # (6) 获取AUGUSTUS Training的准确率信息和预测基因数量统计
    my $input_file = "$tmp_dir/4.augustus/accuary_of_AUGUSTUS_HMM_Training.txt";
    open IN, $input_file or die "Error: Can not open file $input_file, $!";
    print OUT <IN>;
    close IN;
    print OUT "\n";
    
    # (7) 获取基因预测整合过滤的统计信息
    print OUT &statistics_combination();

    # (8) 获取BUSCO分析结果
    if ( -e "$tmp_dir/7.output_gene_models/BUSCO_results.txt" ) {
        my $input_file = "$tmp_dir/7.output_gene_models/BUSCO_results.txt";
        open IN, $input_file or die "Error: Can not open file $input_file, $!";
        print OUT <IN>;
        close IN;
        print OUT "\n";
    }
    close OUT;

    open OUT, ">" , "$tmp_dir/7.output_gene_models.ok" or die $!; close OUT;
}

# 7.6 删除中间文件
if ( $delete_unimportant_intermediate_files ) {
    my @cmdString;
    print STDERR "Due to the --delete_unimportant_intermediate_files was set, so the unimportant and large intermediate files are being deleted...\n";
    # 删除 1.RepeatMasker 文件夹下的数据
    push @cmdString, "rm -rf $tmp_dir/1.RepeatMasker/repeatMasker_Dfam $tmp_dir/1.RepeatMasker/repeatMasker_RepBase $tmp_dir/1.RepeatMasker/repeatModeler";
    # 删除 3.NGSReads_prediction 文件夹下的数据
    push @cmdString, "rm -rf $tmp_dir/3.NGSReads_prediction/a.trimmomatic $tmp_dir/3.NGSReads_prediction/b.hisat2 $tmp_dir/3.NGSReads_prediction/c.transcript";
    # 删除 2.homolog_prediction 文件夹下的数据
    push @cmdString, "rm -rf $tmp_dir/2.homolog_prediction/a.MMseqs2CalHits $tmp_dir/2.homolog_prediction/b.hitToGenePrediction $tmp_dir/2.homolog_prediction/c.getGeneModels";
    # 删除 4.augustus 文件夹下的数据
    push @cmdString, "rm -rf $tmp_dir/4.augustus/aug_para_with_hints";
    # 删除 6.combine_gene_models 文件夹下的数据
    push @cmdString, "rm -rf $tmp_dir/6.combine_gene_models/combineGeneModels_tmp $tmp_dir/6.combine_gene_models/*tmp";
    # 删除 7.output_gene_models 文件夹下的数据
    push @cmdString, "rm -rf $tmp_dir/7.output_gene_models/BUSCO_OUT* rm -rf $tmp_dir/7.output_gene_models/busco_downloads";

    push @cmdString, "6.rm_unimportant_intermediate_files.ok";
    &execute_cmds(@cmdString);
}

print "\n============================================\n";
my $input_file = "$out_prefix.gene_prediction.summary";
open IN, $input_file or die "Error: Can not open file $input_file, $!";
print <IN>;
close IN;

print STDERR "\n============================================\n";
print STDERR "GETA complete successfully! " . "(" . (localtime) . ")" . "\n\n";


sub statistics_combination {
    my $output;

    $output .= "In the process of integrating and filtering the gene models predicted by NGS reads, homolog, and AUGUSTUS, the corresponding statistical data were shown as follows:\n";

    # （1）分析NGS reads、homolog和AUGUSTUS三种方法预测的基因数量，合并后的基因数量。
    my ($augustus_gene_num, $NGSreads_gene_num, $homolog_gene_num) = (0, 0, 0);
    $augustus_gene_num = `grep -P "\tgene\t" $tmp_dir/4.augustus/augustus.gff3 | wc -l`;
    $NGSreads_gene_num = `grep -P "\tgene\t" $tmp_dir/3.NGSReads_prediction/NGSReads_prediction.raw.gff3 | wc -l`;
    $homolog_gene_num = `grep -P "\tgene\t" $tmp_dir/2.homolog_prediction/homolog_prediction.gff3 | wc -l`;
    chomp($augustus_gene_num); chomp($NGSreads_gene_num); chomp($homolog_gene_num);
    my $input_file = "$tmp_dir/6.combine_gene_models/geneModels.a.gff3";
    open IN, $input_file or die "Error: Can not open file $input_file, $!";
    my ($num_of_gene_a, $num_of_gene_b, $num_of_gene_c, $num_of_gene_d, $num_of_gene_all) = (0, 0, 0, 0, 0);
    while (<IN>) {
        $num_of_gene_a ++ if (m/\tgene\t/ && m/augustus/);
        $num_of_gene_c ++ if (m/\tgene\t/ && m/transfrag/);
        $num_of_gene_d ++ if (m/\tgene\t/ && m/genewise/);
    }
    close IN;

    $input_file = "$tmp_dir/6.combine_gene_models/geneModels.b.gff3";
    open IN, $input_file or die "Error: Can not open file $input_file, $!";
    while (<IN>) {
        $num_of_gene_b ++ if m/\tgene\t/;
    }
    close IN;

    $num_of_gene_all = $num_of_gene_a + $num_of_gene_b + $num_of_gene_c + $num_of_gene_d;
    $output .= "(1) We combined the $NGSreads_gene_num, $homolog_gene_num, and $augustus_gene_num gene models predicted by NGS reads, homolog, and AUGUSTUS, respectively. When gene models predicted by different methods overlapped, AUGUSTUS result was retained first, followed by NGS reads, and finally homolog. After the first round of combination, a total of $num_of_gene_all gene models were obtained, including $num_of_gene_a AUGUSTUS gene models with sufficient evidence, $num_of_gene_c gene models predicted by NGS reads, $num_of_gene_d gene models predicted by homolog, and $num_of_gene_b gene models predicted by AUGUSTUS but lacking evidence.\n";
    
    # （2）分析NGS reads和homolog合并后的基因数量，替换AUGUSTUS基因模型的数量，。
    $input_file = "$tmp_dir/4.evidence_gene_models/evidence_gene_models.log";
    open IN, $input_file or die "Error: Can not open file $input_file, $!";
    my ( $evidence_gene_models_num_all, $evidence_gene_models_num_NGSReads, $evidence_gene_models_num_homolog ) = (0, 0, 0);
    while (<IN>) {
        $evidence_gene_models_num_all = $1 if m/Total (\d+) gene models were divided/;
        ($evidence_gene_models_num_NGSReads, $evidence_gene_models_num_homolog) = ($1, $2) if m/ene models predicted by NGSReads: (\d+); predicted by Homolog: (\d+)/;
    }
    close IN;

    $input_file = "$tmp_dir/6.combine_gene_models/picked_evidence_geneModels.log";
    open IN, $input_file or die "Error: Can not open file $input_file, $!";
    my ($evidence_gene_models_num_picked, $evidence_gene_models_num_consistent) = (0, 0);
    <IN>; <IN>; <IN>; <IN>;
    $_ = <IN>; $evidence_gene_models_num_consistent = $1 if m/^.*?\d+.*?\d+.*?(\d+)/;
    $_ = <IN>; $evidence_gene_models_num_picked = $1 if m/(\d+)/;
    close IN;

    my $gene_models_supported_by_evidence_num = `grep -P "\tgene\t" $tmp_dir/6.combine_gene_models/geneModels.d.gff3 | wc -l`;
    chomp($gene_models_supported_by_evidence_num);

    $output .= "(2) After combining the gene models predicted by NGS reads and homologs, $evidence_gene_models_num_all evidence gene models were produced, of which $evidence_gene_models_num_NGSReads came from NGS reads and $evidence_gene_models_num_homolog from homologs. In the second round of combination, $evidence_gene_models_num_picked evidence gene models with long CDS lengths were used to replace the gene models of the previous step. There were $evidence_gene_models_num_consistent evidence gene models with the same CDS structures as predicted by AUGUSTUS at the time of replacement. Finally, $gene_models_supported_by_evidence_num gene models with sufficient evidence were obtained.\n";
    
    # （3）分析HMM和BLASTP进行验证的基因状况
    $input_file = "$tmp_dir/6.combine_gene_models/get_valid_geneModels.log";
    open IN, $input_file or die "Error: Can not open file $input_file, $!";
    $_ = join "", <IN>; 
    close IN;
    my @number = m/(\d+) 个/g;
    $output .= "(3) The HMM and BLASTP algorithms were used to filter transcripts from the predicted gene models. Low confidence mRNAs, such as CDS length less than 600bp or CDSs base proportion less than 30% of exons, had their protein sequences aligned to HMM or BLASTP databases and would be retained if a match was found in any database. An analysis was conducted on $gene_models_supported_by_evidence_num gene models that had evidence support and $num_of_gene_b that did not. Of these, $number[-8] low-confidence gene models were extracted for HMM and BLASTP test, resulting in $number[-2] genes identified as lncRNA, $number[-1] genes filtered as low-quality, and $number[-6] effective gene models remaining after filtering. \n";
    
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
    $output .= "(4) Finally, the gene models overlapping in the CDS region were de-redundant, achieving $num_of_gene_ID gene models, including $num_augustus predicted by AUGUSTUS, $num_transfrag from NGS reads, and $num_genewise from homologs. In addition,  alternative splicing was found in $as_gene_num gene models.\n\n";
    
    return $output;
}

sub detecting_dependent_softwares {
    # 检测依赖的软件
    print STDERR "\n============================================\n";
    print STDERR "Detecting the dependent softwares:\n";
    my $software_info;

    # 检测ParaFly
    $software_info = `ParaFly 2>&1`;
    if ($software_info =~ m/Usage: ParaFly/) {
        print STDERR "ParaFly:\tOK\n";
    }
    else {
        die "ParaFly:\tFailed\n\n";
    }

    # 检测RepeatMasker / RepeatModeler / rmblastn
    if ( $RM_species or $RM_species_Dfam or $RM_species_RepBase or $RM_lib or (! $no_RepeatModeler) ) {
        $software_info = `RepeatMasker`;
        if ($software_info =~ m/RepeatMasker/) {
            print STDERR "RepeatMasker:\tOK\n";
        }
        else {
            die "RepeatMasker:\tFailed\n\n";
        }

        $software_info = `RepeatModeler`;
        if ($software_info =~ m/RepeatModeler/) {
            print STDERR "RepeatModeler:\tOK\n";
        }
        else {
            die "RepeatModeler:\tFailed\n\n";
        }

        $software_info = `rmblastn -version`;
        if ($software_info =~ m/rmblastn/) {
            print STDERR "RepeatModeler:\tOK\n";
        }
        else {
            die "RepeatModeler:\tFailed\n\n";
        }
    }

    # 检测JAVA / HISAT2 / samtools
    if ( ($pe1 && $pe2) or $single_end or $sam ) {
        $software_info = `java -version 2>&1`;
        if ($software_info =~ m/Runtime Environment/) {
            print STDERR "java:\tOK\n";
        }
        else {
            die "java:\tFailed\n\n";
        }

        $software_info = `hisat2 --version`;
        if ($software_info =~ m/version 2.(\d+)\.(\d+)/) {
            print STDERR "HISAT2:\tOK\n";
        }
        else {
            die "HISAT2:\tFailed\n\n";
        }

        $software_info = `samtools --version`;
        if ($software_info =~ m/samtools 1.(\d+)/) {
            print STDERR "samtools:\tOK\n";
        }
        else {
            die "samtools:\tFailed\n\n";
        }
    }

    # 检测mmseqs / genewise / gth / exonerate
    if ( $protein ) {
        # 检测parallel
        $software_info = `parallel --version`;
        if ($software_info =~ m/GNU parallel/) {
            print STDERR "parallel:\tOK\n";
        }
        else {
            die "parallel:\tFailed\n\n";
        }

        $software_info = `mmseqs -h`;
        if ($software_info =~ m/MMseqs2/) {
            print STDERR "mmseqs:\tOK\n";
        }
        else {
            die "mmseqs:\tFailed\n\n";
        }

        my $homolog_prediction_method = $1 if $config{"homolog_prediction"} =~ m/--method\s+(\S+)/;
        $homolog_prediction_method = "exonerate,genewise,gth" if $homolog_prediction_method eq "all";
        foreach ( split /,/, $homolog_prediction_method ) {
            $cmdString = "$_ -version";
            $cmdString =~ s/-version/--version/ if $_ eq "exonerate";
            $software_info = `$cmdString`;
            if ($software_info =~ m/$_/) {
                print STDERR "$_:\tOK\n";
            }
            else {
                die "$_:\tFailed\n\n";
            }
        }
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
    print STDERR "\nPWD: $pwd\n";
    print STDERR (localtime) . ": CMD: $0 $command_line_geta\n\n"; 
}


# 将英文的使用方法放到尾部
sub get_usage_english {

my $usage_english = <<USAGE;

GETA (Genome-wide Electronic Tool for Annotation) is a pipeline software for predicting gene models from whole genome sequences. With one command, you can quickly obtain an accurate genome annotation GFF3 result file by providing the genome sequence, RNA-Seq raw data, and whole genome homologous protein sequences of closely related species. This program has two outstanding features: (1) accurate prediction, GETA outputs gene models with the normal amount, high BUSCO integrity, and accurate exon boundaries; (2) simple run, the operation is straightforward with a fully automated command that produces the final result.

Current Version: 2.7.1

Usage:
    perl $0 [options]

For example:
    perl $0 --genome genome.fasta --RM_species_Dfam Embryophyta --RM_species_RepBase Embryophyta --pe1 libA.1.fq.gz,libB.1.fq.gz --pe2 libA.2.fq.gz,libB.2.fq.gz --protein homolog.fasta --augustus_species GETA_genus_species --HMM_db /opt/biosoft/bioinfomatics_databases/Pfam/PfamA,/opt/biosoft/bioinfomatics_databases/Pfam/PfamB --config /opt/biosoft/geta/conf_for_big_genome.txt --out_prefix out --gene_prefix GS01Gene --cpu 120

Parameters:
[INPUT]
    --genome <string>    default: None, required
    Enter the FASTA file of the genome sequence that you want to annotate. If the input genome file has repeats masked, you can skip the repeat sequence masking step by removing the --RM_species_Dfam,  --RM_species_RepBase, --RM_lib parameters and adding the --no_RepeatModeler parameter. In this instance, it is advised to hard mask the transposable sequences using base N and to soft mask simple or tandem repeats using lowercase characters.

    --RM_species_Dfam | --RM_species <string>    default: None
    Enter the name of a species or class for RepeatMakser to perform a repeat sequence analysis for the genome using the HMM data of corresponding taxonomic species in the Dfam database. The file $software_dir/RepeatMasker_lineage.txt has the values that can be provided for this parameter to represent the class of species. For example, Eukaryota is for eukaryotes, Viridiplantae is for plants, Metazoa is for animals, and Fungi is for fungi. Before attempting to enter this parameter, the RepeatMasker program needed to be installed and the Dfam database needed to be configured. Note that due to its massive size, the Dfam database has been split up into nine partitions. By default, RepeatMasker only contains the zero root partition of the Dfam database, which is suitable for species such as mammals and fungi. If necessary, consider downloading the proper Dfam database partition and configuring it to RepeatMasker. For example, the 5th partition of the Dfam database is designated for Viridiplantae, while the 7th  partition is for Hymenoptera.

    --RM_species_RepBase <string>    default: None
    Enter the name of a species or class for RepeatMakser to perform a repeat sequence analysis for the genome using the nucleotide sequences of corresponding taxonomic species in the RepBase database. The file $software_dir/RepeatMasker_lineage.txt has the values that can be provided for this parameter to represent the class of species. For example, Eukaryota is for eukaryotes, Viridiplantae is for plants, Metazoa is for animals, and Fungi is for fungi. Before attempting to enter this parameter, the RepeatMasker program needed to be installed and the RepBase database needed to be configured. Note that RepBase is no longer providing free downloads and that the most recent version of the database, 20181026, is older and contains few repetitive sequence data.

    --RM_lib <string>    default: None
    Enter a FASTA file and use the repetitive sequence to conduct genome-wide repeat analysis. This file is usually the output of RepeatModeler software's analysis of the entire genome sequence, indicating the repeated sequences across the genome. By default, the GETA program calls RepeatModeler to look up the entire genome sequence and acquire the species' repetitive sequence database. RepeatMakser is then called to search the repeated sequences. After adding this argument, the time-consuming RepeatModler step is skipped, which may significantly reduce the running time of the program. Additionally, the software supports the simultaneous use of the --RM_species_Dfam, --RM_species_RepBase, and --RM_lib arguments, so that multiple methods can be used for repeat sequence analysis, and eventually multiple results can be combined and the result of any method can be recognized.

    --no_RepeatModeler    default: None
    When this parameter is added, the program will no longer run the RepeatModeler step, which is suitable for cases where the repeats have been masked in the input genome file.

    --pe1 <string> --pe2 <string>    default: None
    Enter one or more pairs of FASTQ format files from Paired-End next-generation sequencing technology. This parameter supports the input of multiple pairs of FASTQ files, using commas to separate the FASTQ file paths of different libraries. This parameter also accepts compressed files in .gz format.

    --se <string>    default: None
    Enter one or more FASTQ format files from Single-End next-generation sequencing technology. This parameter supports the input of multiple Single-End FASTQ files, using commas to separate the FASTQ file paths of different libraries. This parameter also accepts compressed files in .gz format.

    --sam <string>    default: None
    Enter one or more SAM format files from the output of alignment software such as HISTA2. This parameter supports the input of multiple SAM files, using commas to separate the SAM file paths. This parameter also accepts compressed files in .bam format. In addition, the program allows for the full or partial use of the three parameters --pe1/--pe2, --se, and --sam, then all of the input data are used for genome alignment to generate the transcript sequence for gene model prediction.

    --strand_specific    default: None
    When this parameter is added, all input next-generation sequencing data are treated as strand-specific, and the program will predict gene models only on the forward strand of the transcript. When two neighboring genes overlap in the genome, strand-specific sequencing data and this parameter can help accurately estimate gene borders.

    --protein <string>    default: None
    Enter a FASTA file containing whole genome protein sequences from neighboring species. It is recommended to use whole genome homologous protein sequences from 3 ~ 10 different species. It is also recommended to modify the name of the protein sequence by appending the Species information, which begins with the species character, to the end of its original name. For example, if the protein sequence is XP_002436309.2, it will be better renamed XP_002436309_2_SpeciesSorghumBicolor. In this way, it is beneficial to retain the homologous matching results of more species in a gene region, improving gene prediction accuracy. The fasta_remove_redundancy.pl script included in GETA can be used to integrate the whole genome protein sequences from multiple species, and species information can be appended while redundancy is removed. The more species employed, the more accurate gene models may be predicted, but the computational time required increases. Note that evidence-supported gene prediction requires at least one type of homologous protein or next-generation sequencing data. 

    --augustus_species <string>    default: None
    When an AUGUSTUS species name is provided, the program starts from an existing species model or retrains a new species model when performing AUGUSTUS Training using gene models predicted by transcripts or homologous proteins. If the input AUGUSTUS species model exists, its parameters will be optimized. If not, a new AUGUSTUS species HMM model will be trained and then its parameters will be optimized. The AUGUSTUS Training step requires the installation of AUGUSTUS software and configuration of the \$AUGUSTUS_CONFIG_PATH environment variable. A species configuration folder from AUGUSTUS Training with the name provided in this parameter is generated in the temporary folder following the program's successful execution. If the user executing the program has write access, the produced species configuration folder can be copied to the species folder specified in \$AUGUSTUS_CONFIG_PATH. If you do not enter this parameter, the program will automatically set the value of this parameter to "GETA + prefix of genome FASTA file name + date + process ID".

    --HMM_db <string>    default: None
    Enter one or more HMM databases, for filtering gene models. This parameter supports the input of multiple HMM databases, separated by commas. The program filters those gene models that do not match in all databases when using multiple HMM databases.

    --BLASTP_db <string>    default: None
    Enter one or more diamond databases, for filtering gene models. This parameter supports the input of multiple diamond databases, separated by commas. The program filters those gene models that do not match in all databases when using multiple diamond databases. When this parameter is left unset, the homologous proteins provided by the --protein parameter will be used to build the diamond database for filtering gene models.

    --config <string>    default: None
    Enter a parameter profile path to set the detailed parameters of other commands called by this program. If this parameter is left unset, When the genome size exceeds 1GB, the software installation directory's conf_for_big_genome.txt configuration file is automatically used. conf_for_small_genome.txt for genome size < 50MB, conf_all_defaults.txt for genome size between 50MB and 1GB.  Additionally, the thresholds for filtering the gene models typically need to be adjusted when GETA predicts an abnormally high number of genes. Then, the GETA pipeline can be rerun by setting this parameter to a new configuration file that is made by modifying the contents of the conf_all_defaults.txt file in the software installation directory.

    --BUSCO_lineage_dataset <string>    default: None
    Enter one or more BUSCO databases, the program will additionally perform BUSCO analysis on the whole genome protein sequences obtained by gene prediction. This parameter supports the input of multiple BUSCO databases, separated by commas. The information contained in the $software_dir/BUSCO_lineages_list.2021-12-14.txt file can be used to choose the proper BUSCO databases. Finally, the BUSCO results are exported to the 7.output_gene_models subdirectory and to the gene_prediction.summary file.

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

    --put_massive_temporary_data_into_memory    default: None
    Set up massive temporary files to be stored in memory. This prevents the program from running slowly due to insufficient disk I/O, but it requires more RAM. Many steps in this pipeline would split the input data into numerous pieces and then parallelize its command lines to speed up the computation, although this results in a significant I/O load on the disk. Therefore, low disk performance has a significant impact on computation speed. If your system memory is sufficient, you are advised to add this parameter so that massive temporary data can be stored in the /dev/shm folder, which represents the memory, to speed up program execution. In addition, the program automatically deletes temporary data in /dev/shm to free up memory after the data splitting and parallelization steps are completed.

    --genetic_code <int>    default: 1
    Enter the genetic code. The values for this parameter can be found on the NCBI Genetic Codes website at: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi. This parameter is mainly effective for the gene prediction steps through homologous proteins, as well as the situation where start and stop codon information is used for filling the end of incomplete gene models.

    --homolog_prediction_method <string>    default: all
    Enter a method for gene prediction using homologous proteins. The value can be set to exonerate, genewise, gth, or all. This parameter supports the input of multiple methods, separated by commas. If the value was set to all, it indicates all three methods were used. The more methods you use, the more computation time you consume, but the better the result will be. Of the three methods, exonerate and genewise produced similar accuracy results, but gth showed a significant decrease in sensitivity and a significant increase in specificity. The following table shows the accuracy of the prediction results for the Oryza sativa genome using three methods. We compared the annotation results of 28736 gene models on NCBI to assess four accuracy metrics: gene level sensitivity, gene level specificity, exon level specificity, and exon level specificity. It is obvious that using multiple methods for gene prediction, combining results, and then filtering can result in a closer number of gene models to the actual number of genes and more accurate results. In addition, this parameter has a higher priority and can override the homolog_prediction parameter value in the parameter configuration file specified by --config.
    Method     Gene_num    gene_sensitivity    gene_specificity    exon_sensitivity    exon_specificity
    exonerate  31310       46.17%              42.37%              62.15%              78.45%
    genewise   32407       46.47%              41.21%              64.66%              77.15%
    gth        9116        20.52%              64.70%              32.71%              90.87%
    all        32293       48.19%              42.89%              66.91%              78.00%
    filtered   23857       45.19%              54.44%              64.44%              83.72%

    --optimize_augustus_method <int>    default: 3
    Enter the method for AUGUSTUS parameters optimization. 1, indicates that only BGM2AT.optimize_augustus is called for optimization, which can fully utilize all CPU threads and run at a fast speed; 2, indicates that only the optimize_augustus.pl program provided by the AUGUSTUS software is used for optimization, which is slower than the first method, but the results are better; 3, indicates that BGM2AT.optimize_augustus is optimized first, then the optimize_augustus.pl program provided by the AUGUSTUS software is used, followed by further optimization. This method prioritizes both speed and effectiveness. It can cover the values of the same parameters specified in the --config parameter configuration file for BGM2AT.
    
    --no_alternative_splicing_analysis    default: None
    When this parameter is added, the program does not perform alternative splicing analysis. Note that GETA defaults to perform alternative splicing analysis based on intron and base sequencing depth information when NGS reads were input.

    --delete_unimportant_intermediate_files    defaults: None
    When this parameter is added and the program runs successfully, the insignificant intermediate files are deleted, leaving only the minimal, small, and important intermediate result files.


This software has been tested and successfully run on Rocky 9.2 system using the following dependent software versions:

01. ParaFly (Version 0.1.0)
02. GNU parallel (Version 20230722)
03. RepeatMasker (version: 4.1.6)
04. RepeatModeler (version: 2.0.5)
05. makeblastdb/rmblastn/tblastn/blastp (Version: 2.14.0)
06. java (version: jdk-20.0.1)
07. hisat2 (version: 2.1.0)
08. samtools (version: 1.17)
09. mmseqs (version 15-6f452)
10. genewise (version: 2.4.1)
11. gth (Vesion: 1.7.3)
12. exonerate (Vesion: 2.2.0)
13. augustus/etraining (version: 3.5.0)
14. diamond (version 2.1.8)
15. hmmscan (version: 3.3.2)
16. busco (Version: 5.4.7)

Version of GETA: 2.7.1

USAGE

return $usage_english;
}

# 子程序，对输入参数进行解析
sub parsing_input_parameters {
    $genome = abs_path($genome);

    $RM_species_Dfam = $RM_species if ( (! $RM_species_Dfam) && $RM_species );
    $RM_lib = abs_path($RM_lib) if $RM_lib;

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

    $protein = abs_path($protein) if $protein;

    die "No genome sequences was found in file $genome\n" unless -s $genome;
    die "No RNA-Seq short reads or homologous proteins was input\n" unless (($pe1 && $pe2) or $single_end or $sam or $protein);
    # 检测AUGUSTUS的环境变量\$AUGUSTUS_CONFIG_PATH
    die "The directory assigned by \$AUGUSTUS_CONFIG_PATH was not exists.\n" unless -e $ENV{"AUGUSTUS_CONFIG_PATH"};

    my $date = `date +%Y%m%d%H%M%S`; chomp($date);
    unless ( $augustus_species ) {
        $augustus_species = basename($genome);
        $augustus_species =~ s/\.fa.*//;
        $augustus_species = "GETA_${augustus_species}${date}_$$";
    }

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
    if ( $BUSCO_lineage_dataset ) {
        foreach ( split /,/, $BUSCO_lineage_dataset ) {
            $_ = abs_path($_);
            push @BUSCO_db, $_ unless exists $BUSCO_db{$_};
            $BUSCO_db{$_} = basename($_);
        }
    }

    $config = abs_path($config) if $config;

    $out_prefix ||= "out";
    $out_prefix = abs_path($out_prefix);

    $gene_prefix ||= "gene";

    $cpu ||= 4;

    $genetic_code ||= 1;
	# 根据遗传密码得到起始密码子和终止密码子
	@_ = &codon_table("$tmp_dir/codon.table");
	my (%start_codon, %stop_codon);
	%start_codon = %{$_[1]};
	%stop_codon = %{$_[2]};
	$start_codon = join ",", sort keys %start_codon;
	$stop_codon = join ",", sort keys %stop_codon;

    die $chinese_help if $chinese_help;
    die $help if $help;

    return 1;
}

sub codon_table {
    my %code = (
        "TTT" => "F",
        "TTC" => "F",
        "TTA" => "L",
        "TTG" => "L",
        "TCT" => "S",
        "TCC" => "S",
        "TCA" => "S",
        "TCG" => "S",
        "TAT" => "Y",
        "TAC" => "Y",
        "TAA" => "X",
        "TAG" => "X",
        "TGT" => "C",
        "TGC" => "C",
        "TGA" => "X",
        "TGG" => "W",
        "CTT" => "L",
        "CTC" => "L",
        "CTA" => "L",
        "CTG" => "L",
        "CCT" => "P",
        "CCC" => "P",
        "CCA" => "P",
        "CCG" => "P",
        "CAT" => "H",
        "CAC" => "H",
        "CAA" => "Q",
        "CAG" => "Q",
        "CGT" => "R",
        "CGC" => "R",
        "CGA" => "R",
        "CGG" => "R",
        "ATT" => "I",
        "ATC" => "I",
        "ATA" => "I",
        "ATG" => "M",
        "ACT" => "T",
        "ACC" => "T",
        "ACA" => "T",
        "ACG" => "T",
        "AAT" => "N",
        "AAC" => "N",
        "AAA" => "K",
        "AAG" => "K",
        "AGT" => "S",
        "AGC" => "S",
        "AGA" => "R",
        "AGG" => "R",
        "GTT" => "V",
        "GTC" => "V",
        "GTA" => "V",
        "GTG" => "V",
        "GCT" => "A",
        "GCC" => "A",
        "GCA" => "A",
        "GCG" => "A",
        "GAT" => "D",
        "GAC" => "D",
        "GAA" => "E",
        "GAG" => "E",
        "GGT" => "G",
        "GGC" => "G",
        "GGA" => "G",
        "GGG" => "G",
    );
    my %start_codon;
    $start_codon{"ATG"} = 1;
    if ( $genetic_code == 1 ) {
        # The Standard Code
        #$start_codon{"TTG"} = 1;
        #$start_codon{"CTG"} = 1;
    }
    elsif ( $genetic_code == 2 ) {
        # The Vertebrate Mitochondrial Code
        $code{"AGA"} = "X";
        $code{"AGG"} = "X";
        $code{"ATA"} = "M";
        $code{"TGA"} = "W";
        $start_codon{"ATA"} = 1;
        $start_codon{"ATT"} = 1;
        $start_codon{"ATC"} = 1;
        $start_codon{"GTG"} = 1;
    }
    elsif ( $genetic_code == 3 ) {
        # The Yeast Mitochondrial Code
        $code{"ATA"} = "M";
        $code{"CTT"} = "T";
        $code{"CTC"} = "T";
        $code{"CTA"} = "T";
        $code{"CTG"} = "T";
        $code{"TGA"} = "W";
        $start_codon{"ATA"} = 1;
        $start_codon{"GTG"} = 1;
    }
    elsif ( $genetic_code == 4 ) {
        # The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
        $code{"TGA"} = "W";
        $start_codon{"ATA"} = 1;
        $start_codon{"ATT"} = 1;
        $start_codon{"ATC"} = 1;
        $start_codon{"GTG"} = 1;
        $start_codon{"CTG"} = 1;
        $start_codon{"TTA"} = 1;
        $start_codon{"TTG"} = 1;
    }
    elsif ( $genetic_code == 5 ) {
        # The Invertebrate Mitochondrial Code
        $code{"AGA"} = "S";
        $code{"AGG"} = "S";
        $code{"ATA"} = "M";
        $code{"TGA"} = "W";
        $start_codon{"ATA"} = 1;
        $start_codon{"ATT"} = 1;
        $start_codon{"ATC"} = 1;
        $start_codon{"GTG"} = 1;
        $start_codon{"TTG"} = 1;
    }
    elsif ( $genetic_code == 6 ) {
        # The Ciliate, Dasycladacean and Hexamita Nuclear Code
        $code{"TAA"} = "Q";
        $code{"TAG"} = "Q";
    }
    elsif ( $genetic_code == 9 ) {
        # The Echinoderm and Flatworm Mitochondrial Code
        $code{"AAA"} = "N";
        $code{"AGA"} = "S";
        $code{"AGG"} = "S";
        $code{"TGA"} = "W";
        $start_codon{"GTG"} = 1;
    }
    elsif ( $genetic_code == 10 ) {
        # The Euplotid Nuclear Code
        $code{"TGA"} = "C";
    }
    elsif ( $genetic_code == 11 ) {
        # The Bacterial, Archaeal and Plant Plastid Code
        $start_codon{"ATA"} = 1;
        $start_codon{"ATT"} = 1;
        $start_codon{"ATC"} = 1;
        $start_codon{"GTG"} = 1;
        $start_codon{"CTG"} = 1;
        $start_codon{"TTG"} = 1;
    }
    elsif ( $genetic_code == 12 ) {
        # The Alternative Yeast Nuclear Code
        $code{"CTG"} = "S";
        $start_codon{"CTG"} = 1;
    }
    elsif ( $genetic_code == 13 ) {
        # The Ascidian Mitochondrial Code
        $code{"AGA"} = "G";
        $code{"AGG"} = "G";
        $code{"ATA"} = "M";
        $code{"TGA"} = "W";
        $start_codon{"ATA"} = 1;
        $start_codon{"GTG"} = 1;
        $start_codon{"TTG"} = 1;
    }
    elsif ( $genetic_code == 14 ) {
        # The Alternative Flatworm Mitochondrial Code
        $code{"AAA"} = "N";
        $code{"AGA"} = "S";
        $code{"AGG"} = "S";
        $code{"TAA"} = "Y";
        $code{"TGA"} = "W";
    }
    elsif ( $genetic_code == 16 ) {
        # Chlorophycean Mitochondrial Code
        $code{"TAG"} = "L";
    }
    elsif ( $genetic_code == 21 ) {
        # Trematode Mitochondrial Code
        $code{"TGA"} = "W";
        $code{"ATA"} = "M";
        $code{"AGA"} = "S";
        $code{"AGG"} = "S";
        $code{"AAA"} = "N";
        $start_codon{"GTG"} = 1;
    }
    elsif ( $genetic_code == 22 ) {
        # Scenedesmus obliquus Mitochondrial Code
        $code{"TCA"} = "X";
        $code{"TAG"} = "L";
    }
    elsif ( $genetic_code == 23 ) {
        # Thraustochytrium Mitochondrial Code
        $code{"TTA"} = "X";
        $start_codon{"ATT"} = 1;
        $start_codon{"GTG"} = 1;
    }
    elsif ( $genetic_code == 24 ) {
        # Rhabdopleuridae Mitochondrial Code
        $code{"AGA"} = "S";
        $code{"AGG"} = "K";
        $code{"TGA"} = "W";
        $start_codon{"GTG"} = 1;
        $start_codon{"CTG"} = 1;
        $start_codon{"TTG"} = 1;
    }
    elsif ( $genetic_code == 25 ) {
        # Candidate Division SR1 and Gracilibacteria Code
        $code{"TGA"} = "G";
        $start_codon{"GTG"} = 1;
        $start_codon{"TTG"} = 1;
    }
    elsif ( $genetic_code == 26 ) {
        # Pachysolen tannophilus Nuclear Code
        # warning: The descritpions of initiation codons by 2 methods are confict according to the NCBI web site.
        $code{"CTG"} = "A";
        $start_codon{"GTG"} = 1;
        $start_codon{"TTG"} = 1;
    }
    elsif ( $genetic_code == 27 ) {
        # Karyorelict Nuclear Code
        $code{"TAG"} = "Q";
        $code{"TAA"} = "Q";
    }
    elsif ( $genetic_code == 29 ) {
        # Mesodinium Nuclear Code
        $code{"TAA"} = "Y";
        $code{"TAG"} = "Y";
    }
    elsif ( $genetic_code == 30 ) {
        # Peritrich Nuclear Code
        $code{"TAA"} = "E";
        $code{"TAG"} = "E";
    }
    elsif ( $genetic_code == 31 ) {
        # Blastocrithidia Nuclear Code
        $code{"TGA"} = "W";
    }
    elsif ( $genetic_code == 33 ) {
        # Cephalodiscidae Mitochondrial UAA-Tyr Code
        $code{"TAA"} = "Y";
        $code{"TGA"} = "Y";
        $code{"AGA"} = "S";
        $code{"AGG"} = "K";
    }

    my %stop_codon;
    foreach ( keys %code ) {
        $stop_codon{$_} = 1 if $code{$_} eq "X";
    }

    return (\%code, \%start_codon, \%stop_codon);
}


# 根据基因组大小选择程序自带的配置文件。
sub choose_config_file {
    my $genomeSize = $_[0];

    # 程序默认使用 conf_all_defaults.txt 配置文件
    my $config_file = "$software_dir/conf_all_defaults.txt";
    # 读取默认配置文件
    open IN, $config_file or die "Error: Can not open file $config_file, $!\n";
    my $tag;
    while (<IN>) {
        s/#.*//; next if m/^\s*$/; s/^\s+//; s/\s*$/ /;
        if (/\[(.*)\]/) {
            $tag = $1; delete $config{$1};
        }
        else {
            $config{$tag} .= $_;
        }
    }
    close IN;

    # 当基因组较大或较小时，自动选择相应的配置文件
    if ( $genomeSize > 1000000000 ) {
        $config_file = "$software_dir/conf_for_big_genome.txt";
    }
    elsif ( $genomeSize < 50000000 ) {
        $config_file = "$software_dir/conf_for_small_genome.txt";
    }
    # 若手动指定了配置文件，则其优先度最高
    $config_file = $config if $config;

    # 读取配置文件，读取参数信息
    open IN, $config_file or die "Can not open file $config_file, $!\n";
    my $tag;
    while (<IN>) {
        s/#.*//; next if m/^\s*$/; s/^\s+//; s/\s*$/ /;
        if (/\[(.*)\]/) {
            $tag = $1; delete $config{$1};
        }
        else {
            $config{$tag} .= $_;
        }
    }
    close IN;

    # 覆盖%config数据
    $homolog_prediction_method ||= "all";
    $homolog_prediction_method = "exonerate,genewise,gth" if $homolog_prediction_method eq "all";
    foreach ( split /,/, $homolog_prediction_method ) {
        die "Error: The supported value of --homolog_prediction_method should be exonerate, genewise or gth, you input '$_' was wrong.\n" unless ( $_ eq "exonerate" or $_ eq "genewise" or $_ eq "gth" );
    }
    $config{"homolog_prediction"} =~ s/--method \S+/--method $homolog_prediction_method/;

    $optimize_augustus_method ||= 3;
    unless ( $optimize_augustus_method == 1 or $optimize_augustus_method == 2 or $optimize_augustus_method == 3 ) {
        die "Error: The value of --optimize_augustus_method shoud be 1, 2 or 3\n";
    }
    $config{"BGM2AT"} =~ s/--optimize_augustus_method \d+/--optimize_augustus_method $optimize_augustus_method/;

    if ( defined $put_massive_temporary_data_into_memory ) {
        $config{"homolog_prediction"} =~ s/\s*$/ --put_massive_temporary_data_into_memory/;
		$config{"BGM2AT"} =~ s/\s*$/ --put_massive_temporary_data_into_memory/;
    }

    # 生成本次程序运行的配置文件
    unless ( -e "$tmp_dir/config.txt" ) {
        open OUT, ">", "$tmp_dir/config.txt" or die "Error: Can not create file $tmp_dir/config.txt, $!";
        foreach ( sort keys %config ) {
            print OUT "[$_]\n$config{$_}\n\n";
        }
        close OUT;
    }

    return 1;
}

# 子程序，用于执行调用的Linux命令。程序运行完毕后，生成.ok文件。
sub execute_cmds {
    my $ok_file = pop @_;

    if ( -e $ok_file ) {
        foreach ( @_ ) {
            print STDERR "CMD(Skipped): $_\n";
        }
    }
    else {
        foreach ( @_ ) {
            print STDERR (localtime) . ": CMD: $_\n";
            system($_) == 0 or die "failed to execute: $_\n";
        }
        open OUT, ">", "$ok_file" or die $!; close OUT;
    }

    return 1;
}
