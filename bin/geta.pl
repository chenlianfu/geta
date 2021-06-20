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

# Step 1: Trimmomatic
print STDERR "\n============================================\n";
print STDERR "Step 1: Trimmomatic " . "(" . (localtime) . ")" . "\n";
mkdir "1.trimmomatic" unless -e "1.trimmomatic";
unless (-e "1.trimmomatic.ok") {
    chdir "1.trimmomatic";
    $pwd = `pwd`; print STDERR "PWD: $pwd";
    if ($pe1 && $pe2) {
        my @pe_reads = sort keys %pe_reads;
        my $pe_reads_num = @pe_reads;
        my $number = 0;
        foreach (@pe_reads) {
            $number ++;
            my $code = "0" x ( length($pe_reads_num) - length($number) ) . $number;
            @_ = split /\t/;
            $cmdString = "java -jar $dirname/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads $cpu $_[0] $_[1] reads$code.1.fastq reads$code.1.unpaired.fastq reads$code.2.fastq reads$code.2.unpaired.fastq ILLUMINACLIP:$dirname/Trimmomatic-0.38/adapters/$config{'trimmomatic'} &> trimmomatic.pe.log";
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        }
        if ($pe_reads_num == 1) {
            $cmdString = "ln -sf reads1.1.fastq reads.1.fastq && ln -sf reads1.2.fastq reads.2.fastq";
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        }
        else {
            $cmdString = "cat reads*.1.fastq > reads.1.fastq && cat reads*.2.fastq > reads.2.fastq";
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        }
    }
    if ($single_end) {
        my @se_reads = sort keys %se_reads;
        my $se_reads_num = @se_reads;
        my $number = 0;
        foreach (@se_reads) {
            $number ++;
            my $code = "0" x ( length($se_reads_num) - length($number) ) . $number;
            @_ = split /\t/;
            $cmdString = "java -jar $dirname/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads $cpu $single_end reads$code.fastq ILLUMINACLIP:$dirname/Trimmomatic-0.38/adapters/$config{'trimmomatic'} &> trimmomatic.single.log";
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        }
        if ($se_reads_num == 1) {
            $cmdString = "ln -sf reads1.fastq reads.fastq";
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        }
        else {
            $cmdString = "cat reads*.fastq > reads.fastq";
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        }
    }
    chdir "../";
    open OUT, ">", "1.trimmomatic.ok" or die $!; close OUT;
}
else {
    print STDERR "Skip Step 1 for the file 1.trimmomatic.ok exists\n";
}

# Step 2: HISAT2
print STDERR "\n============================================\n";
print STDERR "Step 2: HISAT2 " . "(" . (localtime) . ")" . "\n";
mkdir "2.hisat2" unless -e "2.hisat2";
unless (($pe1 && $pe2) or $single_end) {
    open OUT, ">", "2.hisat2.ok" or die $!; close OUT;
}
unless (-e "2.hisat2.ok") {
    chdir "2.hisat2";
    $pwd = `pwd`; print STDERR "PWD: $pwd";

    # 构建基因组hisat2索引数据库
    $cmdString = "hisat2-build $config{'hisat2-build'} ../0.RepeatMasker/genome.masked.fasta genome &> hisat2-build.log\n";
    unless (-e "hisat2-build.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "hisat2-build.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    # 将RNA-Seq reads比对到参考基因组
    # 能处理单端，双端，链特异性与链非特异性数据
    my $input;
    if ($pe1 && $pe2) {
        $input = "-1 ../1.trimmomatic/reads.1.fastq -2 ../1.trimmomatic/reads.2.fastq";
    }
    if ($single_end) {
        $input .= " -U ../1.trimmomatic/reads.fastq"
    }
    if ($strand_specific) {
        $input .= " --rna-strandness RF";
    }
    $cmdString = "hisat2 -x genome -p $cpu $input -S hisat2.sam $config{'hisat2'} 2> hisat2.log";
    unless (-e "hisat2.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        print STDERR "Warning: The HISAT2 parameter --rna-strandness may not set corretly !" if ( ($config{'hisat2'} =~ m/RF/) == 0 && $strand_specific);
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "hisat2.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    $cmdString = "samtools sort  -o hisat2.sorted.bam -O BAM hisat2.sam";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

    $cmdString = "samtools view -h hisat2.sorted.bam > hisat2.sorted.sam";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

    chdir "../";
    open OUT, ">", "2.hisat2.ok" or die $!; close OUT;
}
else {
    print STDERR "Skip Step 2 for the file 2.hisat2.ok exists\n";
}

# Step 3: Transcript
print STDERR "\n============================================\n";
print STDERR "Step 3: Transcript " . "(" . (localtime) . ")" . "\n";
mkdir "3.transcript" unless -e "3.transcript";
unless (($pe1 && $pe2) or $single_end) {
    open OUT, ">", "3.transcript.ok" or die $!; close OUT;
    chdir "3.transcript";
    open OUT, ">", "transfrag.genome.gff3" or die $!; close OUT;
    chdir "../";
}
unless (-e "3.transcript.ok") {
    chdir "3.transcript";
    $pwd = `pwd`; print STDERR "PWD: $pwd";

    # 将SAM比对将诶过分成一个个比对区域
    $cmdString = "$dirname/bin/split_sam_from_non_aligned_region ../2.hisat2/hisat2.sorted.sam splited_sam_out 10 > splited_sam_files.list";
    unless (-e "split_sam_from_non_aligned_region.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "split_sam_from_non_aligned_region.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }


    # 每个比对区域生成一个sam2transfrag命令
    open IN, "splited_sam_files.list" or die $!;
    open CMD, ">", "command.sam2transfrag.list" or die $!;
    my $no_strand_specific = "--no_strand_specific";
    $no_strand_specific = "" if $strand_specific;
    while (<IN>) {
        s/\.sam\n//;
        print CMD "$dirname/bin/sam2transfrag $config{'sam2transfrag'} $no_strand_specific --intron_info_out $_.intron $_.sam > $_.gtf\n";
    }
    close CMD;
    close IN;

    # 批量并行化进行transcripts计算
    $cmdString = "ParaFly -c command.sam2transfrag.list -CPU $cpu &> /dev/null";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

    # 整合所有并行化结果，得到GTF文件和可信的Intron信息
    open IN, "splited_sam_files.list" or die $!;
    open OUT1, ">", "transfrag.gtf" or die $!;
    open OUT2, ">", "intron.txt"  or die $!;
    while (<IN>) {
        s/\.sam\n//;
        open IN1, "$_.gtf" or die "Cannot open file $_.gtf, $!\n";
        print OUT1 <IN1>;
        close IN1;
        if (-e "$_.intron") {
            open IN1, "$_.intron" or die $!;
            print OUT2 <IN1>;
            close IN1;
        }
    }
    close OUT2; close OUT1; close IN;

    # 将GTF文件转换成GFF3文件和transcripts序列
    # 若是非链特异性测序，则得到的single exon转录本序列是没有方向的。
    unless (-e "transfragDump.ok") {
        $cmdString = "$dirname/bin/transfragDump transfrag.gtf $genome 2> transfrag.stats";
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "transfragDump.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    # 对transcripts序列使用Transdecoder进行ORF分析
    unless (-e "TransDecoder.ok") {
        $cmdString = "$dirname/TransDecoder-v5.5.0/TransDecoder.LongOrfs $config{'TransDecoder.LongOrfs'} -t transfrag.strand.fasta -S &> /dev/null";
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        $cmdString = "$dirname/TransDecoder-v5.5.0/TransDecoder.Predict $config{'TransDecoder.Predict'}  -t transfrag.strand.fasta &> /dev/null";
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        $cmdString = "cp transfrag.strand.fasta.transdecoder.gff3 transfrag.transdecoder.gff3";
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        unless ($strand_specific) {
            $cmdString = "$dirname/TransDecoder-v5.5.0/TransDecoder.LongOrfs $config{'TransDecoder.LongOrfs'} -t transfrag.noStrand.fasta &> /dev/null";
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            $cmdString = "$dirname/TransDecoder-v5.5.0/TransDecoder.Predict $config{'TransDecoder.Predict'} -t transfrag.noStrand.fasta --train transfrag.strand.fasta.transdecoder_dir/longest_orfs.cds.top_500_longest &> /dev/null";
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            $cmdString = "cat transfrag.noStrand.fasta.transdecoder.gff3 >> transfrag.transdecoder.gff3";
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        }
        open OUT, ">", "TransDecoder.ok" or die $!; close OUT;
    }
    else {
        $cmdString = "$dirname/TransDecoder-v5.5.0/TransDecoder.LongOrfs $config{'TransDecoder.LongOrfs'} -t transfrag.strand.fasta -S &> /dev/null";
        print STDERR "CMD(Skipped): $cmdString\n";
        $cmdString = "$dirname/TransDecoder-v5.5.0/TransDecoder.Predict $config{'TransDecoder.Predict'} -t transfrag.strand.fasta &> /dev/null";
        print STDERR "CMD(Skipped): $cmdString\n";
        $cmdString = "cp transfrag.strand.fasta.transdecoder.gff3 transfrag.transdecoder.gff3";
        print STDERR "CMD(Skipped): $cmdString\n";
        unless ($strand_specific) {
            $cmdString = "$dirname/TransDecoder-v5.5.0/TransDecoder.LongOrfs $config{'TransDecoder.LongOrfs'} -t transfrag.noStrand.fasta &> /dev/null";
            print STDERR "CMD(Skipped): $cmdString\n";
            $cmdString = "$dirname/TransDecoder-v5.5.0/TransDecoder.Predict $config{'TransDecoder.Predict'} -t transfrag.noStrand.fasta --train transfrag.strand.fasta.transdecoder_dir/longest_orfs.cds.top_500_longest &> /dev/null";
            print STDERR "CMD(Skipped): $cmdString\n";
            $cmdString = "cat transfrag.noStrand.fasta.transdecoder.gff3 >> transfrag.transdecoder.gff3";
            print STDERR "CMD(Skipped): $cmdString\n";
        }
    }

    # 将transcripts的ORF预测结果映射到基因组序列上，得到transcripts的基因预测结果： transfrag.genome.gff3
    $cmdString = "$dirname/bin/transdecoder2ORF --out_protein proteins.fasta transfrag.gtf transfrag.transdecoder.gff3 $genome > transdecoder2ORF.gff3";
    unless (-e "transdecoder2ORF.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "transdecoder2ORF.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    # 对 transdecoder2ORF.gff3 进行 GFF3Clear，得到transcripts的基因预测结果： transfrag.genome.gff3
    $cmdString = "$dirname/bin/GFF3Clear --GFF3_source GETA --genome $genome --gene_prefix transfrag --no_attr_add transdecoder2ORF.gff3 > transfrag.genome.gff3 2> /dev/null";
    unless (-e "GFF3Clear.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "GFF3Clear.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

=cut
    # 提取完整的 transcripts 基因预测结果
    open IN, "transfrag.genome.gff3" or die $!;    
    open OUT, ">", "transfrag.genome.complete.gff3" or die $!;
    my $keep = 0;
    while (<IN>) {
        if (m/\tgene\t/) {
            if (m/Form=one_transcript_get_1_gene_model.*Integrity=complete/) {
                $keep = 1;
            }
            else {
                $keep = 0;
            }
        }
        print OUT if $keep == 1;
    }
    close OUT; close IN;

    # 从完整的 transcripts 基因预测结果中提取最优基因模型
    $cmdString = "$dirname/bin/ORF2bestGeneModels $ORF2bestGeneModels transfrag.genome.complete.gff3 > best_candidates.gff3 2> ORF2bestGeneModels.log";
    unless (-e "ORF2bestGeneModels.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "ORF2bestGeneModels.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    # 过滤一些基因模型，以保证两两基因模型之间在protein序列上Identity要低于80%
    $cmdString = "$dirname/bin/bestGeneModels2lowIdentity best_candidates.gff3 proteins.fasta $cpu 0.8 > best_candidates.lowIdentity.gff3 2> bestGeneModels2lowIdentity.log";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
=cut

    chdir "../";
    open OUT, ">", "3.transcript.ok" or die $!; close OUT;
}
else {
    print STDERR "Skip Step 3 for the file 3.transcript.ok exists\n";
}

# Step 4: Homolog
print STDERR "\n============================================\n";
print STDERR "Step 4: Homolog " . "(" . (localtime) . ")" . "\n";
mkdir "4.homolog" unless -e "4.homolog";
unless (-e "4.homolog.ok") {
    chdir "4.homolog";
    $pwd = `pwd`; print STDERR "PWD: $pwd";

    my $max_gene_length = 20000;
    open IN, "../3.transcript/transfrag.genome.gff3" or die "Can not open the file ../3.transcript/transfrag.genome.gff3, $!\n";
    my @gene_length;
    while (<IN>) {
        if (m/\tgene\t(\d+)\t(\d+)\t/) {
            push @gene_length, $2 - $1 + 1;
        }
    }
    @gene_length = sort {$a <=> $b} @gene_length;
    $max_gene_length = $gene_length[@gene_length * 0.99] if $gene_length[@gene_length * 0.99] > $max_gene_length;
    my ($segmentSize, $overlapSize) = (1000000, 100000);
    if ($max_gene_length * 4 > $overlapSize) {
        $overlapSize = $max_gene_length * 4;
        my $overlapSize_length = length($overlapSize);
        $overlapSize_length --;
        $overlapSize_length --;
        $overlapSize = int(($overlapSize / (10 ** $overlapSize_length)) + 1) * (10 ** $overlapSize_length);
        $segmentSize = $overlapSize * 10;
    }
    $cmdString = "$dirname/bin/homolog_genewise --cpu $cpu --max_gene_length $max_gene_length --segmentSize $segmentSize --overlapSize $overlapSize $config{'homolog_genewise'} $protein ../0.RepeatMasker/genome.masked.fasta &> homolog_genewise.log";
    unless (-e "homolog_genewise.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "homolog_genewise.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    $cmdString = "$dirname/bin/homolog_genewiseGFF2GFF3 $config{'homolog_genewiseGFF2GFF3'} --input_genewise_start_info genewise.start_info.txt --output_start_and_stop_hints_of_augustus genewise.start_stop_hints.gff --genome $genome genewise.gff > genewise.gff3 2> genewise.gene_id_with_stop_codon.txt";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

    $cmdString = "$dirname/bin/GFF3Clear --genome $genome --gene_prefix genewise --no_attr_add genewise.gff3 > out.gff3 2> /dev/null";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

    $cmdString = "mv out.gff3 genewise.gff3";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

    chdir "../";
    open OUT, ">", "4.homolog.ok" or die $!; close OUT;
}
else {
    print STDERR "Skip Step 4 for the file 4.homolog.ok exists\n";
}

# Step 5: Augustus gene prediction
print STDERR "\n============================================\n";
print STDERR "Step 5: Augustus/HMM Trainning " . "(" . (localtime) . ")" . "\n";
mkdir "5.augustus" unless -e "5.augustus";
unless (-e "5.augustus.ok") {
    chdir "5.augustus";
    $pwd = `pwd`; print STDERR "PWD: $pwd";

    if ($use_existed_augustus_species) {
        mkdir "training";
        open OUT, ">", "training.ok" or die $!;
        print STDERR "Skip Augustus training for --use_existed_augustus_species paramter set\n";
    }

    # 第一次 Augustus HMM Training
    unless (-e "training") {
        mkdir "training";
        my $species_config_dir = `echo \$AUGUSTUS_CONFIG_PATH`;
        chomp($species_config_dir);
        $species_config_dir = "$species_config_dir/species/$augustus_species";
        $cmdString = "rm -rf $species_config_dir";
        print STDERR "CMD: $cmdString\n";
        (system $cmdString) == 0 or die "Failed to execute: $cmdString\n";
        unlink "training.ok" if (-e "training.ok");
    }
    unless (-e "training.ok") {
        chdir "training";
        $pwd = `pwd`; print STDERR "PWD: $pwd";

        # 准备Augustus training的输入文件
        open OUT, ">", "blank.augustus.gff3" or die "Can not create file blank.augustus.gff3, $!\n";
        close OUT;
        open OUT, ">", "blank.intron.gff" or die "Can not create file blank.intron.gff, $!\n";
        close OUT;
        $cmdString = "$dirname/bin/paraCombineGeneModels --cpu $cpu $config{'paraCombineGeneModels'} blank.augustus.gff3 ../../3.transcript/transfrag.genome.gff3 ../../4.homolog/genewise.gff3 blank.intron.gff";
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

        $cmdString = "$dirname/bin/GFF3Clear --genome $genome combine.1.gff3 > geneModels.gff3 2> GFF3Clear.log";
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

        $cmdString = "$dirname/bin/geneModels2AugusutsTrainingInput $config{'geneModels2AugusutsTrainingInput'} --out_prefix ati --cpu $cpu geneModels.gff3 $genome &> geneModels2AugusutsTrainingInput.log";
        unless (-e "geneModels2AugusutsTrainingInput.ok") {
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
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
        $cmdString = "$dirname/bin/BGM2AT $config{'BGM2AT'} --flanking_length $flanking_length --CPU $cpu --onlytrain_GFF3 ati.filter1.gff3 ati.filter2.gff3 $genome $augustus_species &> BGM2AT.log";
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

        chdir "../";
        open OUT, ">", "training.ok" or die $!;
    }
    else {
        print STDERR "Skip Augustus training for file training.ok exists\n" unless $use_existed_augustus_species;
    }

    # Augustus Hint Preparing
    $cmdString = "bam2hints --source=W --intronsonly --in=../2.hisat2/hisat2.sorted.bam --out=bam2intronHints.gff";
    unless (($pe1 && $pe2) or $single_end) {
        open OUT, ">", "bam2hints.ok" or die $!; close OUT;
        open OUT, ">", "bam2intronHints.gff" or die $!; close OUT;
    }
    unless (-e "bam2hints.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "bam2hints.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }
    $cmdString = "$dirname/bin/prepareAugusutusHints $config{'prepareAugusutusHints'} bam2intronHints.gff ../3.transcript/transfrag.genome.gff3 ../4.homolog/genewise.gff3 ../4.homolog/genewise.start_stop_hints.gff > hints.gff";
    unless (-e "prepareAugusutusHints.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "prepareAugusutusHints.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    # Get the longest gene length
    open IN, "../3.transcript/transfrag.genome.gff3" or die "Can not open the file ../3.transcript/transfrag.genome.gff3, $!\n";
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
    unless ($no_augustus_training_iteration or $use_existed_augustus_species) {
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
                    if ($attr{"hintRatio"} >= 95) {
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

                $cmdString = "ln -sf augustus.2.gff3 augustus.gff3";
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

                $cmdString = "ln -sf augustus.1.gff3 augustus.gff3";
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

                $cmdString = "ln -sf augustus.2.gff3 augustus.gff3";
                print STDERR (localtime) . ": CMD: $cmdString\n";
                system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            }
        }
    }
    else {
        $cmdString = "ln -sf augustus.1.gff3 augustus.gff3";
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    }

    chdir "../";
    open OUT, ">", "5.augustus.ok" or die $!; close OUT;
}
else {
    print STDERR "Skip Step 5 for the file 5.augustus.ok exists\n";
}



# Step 6: CombineGeneModels
print STDERR "\n============================================\n";
print STDERR "Step 6: CombineGeneModels " . "(" . (localtime) . ")" . "\n";
mkdir "6.combineGeneModels" unless -e "6.combineGeneModels";
unless (-e "6.combineGeneModels.ok") {
    chdir "6.combineGeneModels";
    $pwd = `pwd`; print STDERR "PWD: $pwd";

    $cmdString = "$dirname/bin/paraCombineGeneModels $config{'paraCombineGeneModels'} --cpu $cpu ../5.augustus/augustus.gff3 ../3.transcript/transfrag.genome.gff3 ../4.homolog/genewise.gff3 ../5.augustus/hints.gff &> /dev/null";
    unless (-e "paraCombineGeneModels.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "paraCombineGeneModels.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    if ($pfam_db) {
        #$cmdString = "rm -rf command.hmmscan.list* hmmscan.tmp for_pfam_search.fasta";
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

        $cmdString = "$dirname/bin/PfamValidateABinitio --out_prefix combine2 --cpu $cpu --pfam_db $pfam_db $config{'PfamValidateABinitio'} combine.2.gff3 $genome 2> PfamValidateABinitio.1.log";
        unless (-e "PfamValidateABinitio.1.ok") {
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            open OUT, ">", "PfamValidateABinitio.1.ok" or die $!; close OUT;
        }
        else {
            print STDERR "CMD(Skipped): $cmdString\n";
        }
    }
    else {
        $cmdString = "cp combine.2.gff3 combine2.filter_pass.gff3";
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    }

    $cmdString = "$dirname/bin/GFF3Clear --gene_prefix $gene_prefix --genome $genome combine.1.gff3 combine2.filter_pass.gff3 > genome.gff3 2> GFF3Clear.1.log";
    unless (-e "GFF3Clear.1.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "GFF3Clear.1.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    open OUT1, ">", "genome.completed.gff3" or die $!;
    open OUT2, ">", "genome.partial.gff3" or die $!;
    # 注意，有部分AUGUSTUS预测的基因，其预测结果中不是以ATG作为起始密码子。
    open IN, "genome.gff3" or die $!;
    my $complete_keep;
    $/ = "\n\n";
    while (<IN>) {
        #next if m/^\s*$/;
        if (m/Integrity=complete/) {
            $complete_keep = 1;
        }
        elsif (m/source=augustus/) {
            $complete_keep = 1;
        } else {
            $complete_keep = 0;
        }
        print OUT1 if $complete_keep == 1;
        print OUT2 if $complete_keep == 0;
    }
    $/ = "\n";

    $cmdString = "$dirname/bin/remove_genes_in_repeats $config{'remove_genes_in_repeats'} --filtered_gene_models genome.completed.genes_in_repeats.gff3 ../0.RepeatMasker/genome.repeat.gff3 genome.completed.gff3 > genome.completed.rm_genes_in_repeats.gff3 2> remove_genes_in_repeats.txt";
    unless (-e "remove_genes_in_repeats.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "remove_genes_in_repeats.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    if ($pfam_db) {
        $cmdString = "rm -rf command.hmmscan.list* hmmscan.tmp for_pfam_search.fasta";
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

        $cmdString = "$dirname/bin/remove_short_genes $config{'remove_short_genes'} genome.completed.rm_genes_in_repeats.gff3 > genome.completed.rm_genes_in_repeats.remove_short_genes.gff3 2> genome.completed.rm_genes_in_repeats.short_genes.gff3";
        unless (-e "remove_short_genes.ok") {
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            open OUT, ">", "remove_short_genes.ok" or die $!; close OUT;
        }
        else {
            print STDERR "CMD(Skipped): $cmdString\n";
        }

        $cmdString = "$dirname/bin/PfamValidateABinitio --out_prefix remove_short_genes --cpu $cpu --pfam_db $pfam_db $config{'PfamValidateABinitio'} genome.completed.rm_genes_in_repeats.short_genes.gff3 $genome 2> PfamValidateABinitio.2.log";
        unless (-e "PfamValidateABinitio.2.ok") {
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            open OUT, ">", "PfamValidateABinitio.2.ok" or die $!; close OUT;
        }
        else {
            print STDERR "CMD(Skipped): $cmdString\n";
        }

        $cmdString = "$dirname/bin/GFF3Clear --gene_prefix $gene_prefix --genome $genome --no_attr_add genome.completed.rm_genes_in_repeats.remove_short_genes.gff3 remove_short_genes.filter_pass.gff3 > genome.filter.gff3 2> GFF3Clear.2.log";
        unless (-e "GFF3Clear.2.ok") {
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            open OUT, ">", "GFF3Clear.2.ok" or die $!; close OUT;
        }
        else {
            print STDERR "CMD(Skipped): $cmdString\n";
        }
    }
    else {
        $cmdString = "cp genome.completed.rm_genes_in_repeats.gff3 genome.filter.gff3";
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    }

    chdir "../";
    open OUT, ">", "6.combineGeneModels.ok" or die $!; close OUT;
}
else {
    print STDERR "Skip Step 6 for the file 6.combineGeneModels.ok exists\n";
}

=cut
# Step 7: Add Alternative Splicing
print STDERR "\n============================================\n";
print STDERR "Step 7: Add Alternative Splicing\n";
mkdir "7.addAlternativeSplicing" unless -e "7.addAlternativeSplicing";
unless (-e "7.addAlternativeSplicing.ok") {
    chdir "7.addAlternativeSplicing";
    $pwd = `pwd`; print STDERR "PWD: $pwd";

    $cmdString = "stringtie ../2.hisat2/hisat2.sorted.bam -o stringtie.gtf -p $cpu";
    unless (-e "stringtie.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "stringtie.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    $cmdString = "$dirname/bin/GFF3Clear --gene_prefix $gene_prefix --genome $genome ../6.combineGeneModels/genome.completed.gff3 > genome.gff3";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

    $cmdString = "$dirname/bin/addAS_from_stringtie genome.gff3 stringtie.gtf $genome > genome.addAS.gff3 2> AS_source.txt";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

    chdir "../";
    open OUT, ">", "7.addAlternativeSplicing.ok" or die $!; close OUT;
}
else {
    print STDERR "Skip Step 7 for the file 7.addAlternativeSplicing.ok exists\n";
}
=cut

# Step 7: OutPut
print STDERR "\n============================================\n";
print STDERR "Step 7: OutPut " . "(" . (localtime) . ")" . "\n";
chdir "../";
$pwd = `pwd`; print STDERR "PWD: $pwd";

$cmdString = "$dirname/bin/GFF3Clear --GFF3_source GETA --gene_prefix $gene_prefix --genome $genome --no_attr_add $out_prefix.tmp/6.combineGeneModels/genome.filter.gff3 > $out_prefix.GeneModels.gff3 2> /dev/null";
#$cmdString = "cp $out_prefix.tmp/7.addAlternativeSplicing/genome.addAS.gff3 $out_prefix.gff3";
print STDERR (localtime) . ": CMD: $cmdString\n";
system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

$cmdString = "$dirname/bin/GFF3Clear --GFF3_source GETA --gene_prefix ${gene_prefix}Broken --genome $genome --no_attr_add $out_prefix.tmp/6.combineGeneModels/genome.partial.gff3 > $out_prefix.Incomplete_GeneModels.gff3 2> /dev/null";
print STDERR (localtime) . ": CMD: $cmdString\n";
system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

$cmdString = "$dirname/bin/GFF3Clear --GFF3_source GETA --gene_prefix ${gene_prefix}InRepeatRegion --genome $genome --no_attr_add $out_prefix.tmp/6.combineGeneModels/genome.completed.genes_in_repeats.gff3 > $out_prefix.InRepeatRegion_GeneModels.gff3 2> /dev/null";
print STDERR (localtime) . ": CMD: $cmdString\n";
system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

if ($pfam_db) {
    $cmdString = "$dirname/bin/GFF3Clear --GFF3_source GETA --gene_prefix ${gene_prefix}ShortGene --genome $genome --no_attr_add $out_prefix.tmp/6.combineGeneModels/remove_short_genes.filter_out.gff3 > $out_prefix.ShortCDS_GeneModels.gff3 2> /dev/null";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
}

$cmdString = "$dirname/bin/gff3ToGtf.pl $genome $out_prefix.GeneModels.gff3 > $out_prefix.GeneModels.gtf 2> /dev/null";
print STDERR (localtime) . ": CMD: $cmdString\n";
system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

$cmdString = "$dirname/bin/eukaryotic_gene_model_statistics.pl $out_prefix.GeneModels.gtf $genome $out_prefix &> $out_prefix.GeneModels.stats";
print STDERR (localtime) . ": CMD: $cmdString\n";
system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

open IN, "$out_prefix.tmp/0.RepeatMasker/genome.repeat.stats" or die "Can not open file $out_prefix.tmp/0.RepeatMasker/genome.repeat.stats, $!";
open OUT, ">", "$out_prefix.repeat.stats" or die "Can not create file $out_prefix.repeat.stats, $!";
print OUT <IN>;
close IN; close OUT;

open IN, "$out_prefix.tmp/0.RepeatMasker/genome.repeat.gff3" or die "Can not open file $out_prefix.tmp/0.RepeatMasker/genome.repeat.gff3, $!";
open OUT, ">", "$out_prefix.repeat.gff3" or die "Can not create file $out_prefix.repeat.gff3, $!";
print OUT <IN>;
close IN; close OUT;

if (($pe1 && $pe2) or $single_end) {
    open IN, "$out_prefix.tmp/3.transcript/transfrag.alignment.gff3" or die "Can not open file $out_prefix.tmp/3.transcript/transfrag.alignment.gff3, $!";
    open OUT, ">", "$out_prefix.transfrag_alignment.gff3" or die "Can not create file $out_prefix.transfrag_alignment.gff3, $!";
    print OUT <IN>;
    close IN; close OUT;

    open IN, "$out_prefix.tmp/3.transcript/transfrag.genome.gff3" or die "Can not open file $out_prefix.tmp/3.transcript/transfrag.genome.gff3, $!";
    open OUT, ">", "$out_prefix.transfrag_prediction.gff3" or die "Can not create file $out_prefix.transfrag_prediction.gff3, $!";
    print OUT <IN>;
    close IN; close OUT;
}

if ($protein) {
    open IN, "$out_prefix.tmp/4.homolog/genewise.gff3" or die "Can not open file $out_prefix.tmp/4.homolog/genewise.gff3, $!";
    open OUT, ">", "$out_prefix.homolog_prediction.gff3" or die "Can not create file $out_prefix.homolog_prediction.gff3, $!";
    print OUT <IN>;
    close IN; close OUT;
}

open OUT, ">", "$out_prefix.gene_prediction.summary" or die "Can not create file $out_prefix.gene_prediction.summary, $!";
open IN, "$out_prefix.tmp/0.RepeatMasker/genome.repeat.stats" or die "Can not open file $out_prefix.tmp/0.RepeatMasker/genome.repeat.stats, $!";
while (<IN>) {
    print OUT $_ if m/^Genome Size/;
    print OUT "$_\n" if m/^Repeat Ratio/;
}
close IN;
if ( -e "$out_prefix.tmp/2.hisat2/hisat2.log" ) {
    open IN, "$out_prefix.tmp/2.hisat2/hisat2.log" or die $!;
    while (<IN>) {
        print OUT "The alignment rate of RNA-Seq reads is: $1\n\n" if m/^(\S+) overall alignment rate/;
    }
    close IN;
    open IN, "$out_prefix.tmp/3.transcript/transfrag.genome.gff3" or die "Can not open file $out_prefix.tmp/3.transcript/transfrag.genome.gff3, $!";
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
if ( $protein ) {
    open IN, "$out_prefix.tmp/4.homolog/genewise.gff3" or die "Can not open file $out_prefix.tmp/4.homolog/genewise.gff3, $!";
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
if (-e "$out_prefix.tmp/5.augustus/augustus.2.gff3") {
    open IN, "$out_prefix.tmp/5.augustus/training_again/secondtest.out" or die "Can not open file $out_prefix.tmp/5.augustus/training_again/secondtest.out, $!";
}
else {
    open IN, "$out_prefix.tmp/5.augustus/training/secondtest.out" or die "Can not open file $out_prefix.tmp/5.augustus/training/secondtest.out, $!";
}
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
$out = "The accuary of AUGUSTUS Training is $accuary\%.\nLevel\tSensitivity\tSpecificity\n$out\n";
print OUT $out;
my ($gene_num1, $gene_num2, $gene_num3, $gene_num4) = (0, 0, 0, 0);
open IN, "$out_prefix.GeneModels.gff3" or die "Can not open file $out_prefix.GeneModels.gff3, $!";
while (<IN>) {
    $gene_num1 ++ if m/\tgene\t/;
}
close IN;
open IN, "$out_prefix.Incomplete_GeneModels.gff3" or die "Can not open file $out_prefix.Incomplete_GeneModels.gff3, $!";
while (<IN>) {
    $gene_num2 ++ if m/\tgene\t/;
}
close IN;
open IN, "$out_prefix.InRepeatRegion_GeneModels.gff3" or die "Can not open file $out_prefix.InRepeatRegion_GeneModels.gff3, $!";
while (<IN>) {
    $gene_num3 ++ if m/\tgene\t/;
}
close IN;
open IN, "$out_prefix.ShortCDS_GeneModels.gff3" or die "Can not open file $out_prefix.ShortCDS_GeneModels.gff3, $!";
while (<IN>) {
    $gene_num4 ++ if m/\tgene\t/;
}
close IN;
print OUT "The number of final complete gene models is: $gene_num1\n";
print OUT "The number of Incomplete gene models is: $gene_num2\n";
print OUT "The number of gene models whose CDS region overlap to Repeat Region > $1 is: $gene_num3\n" if $config{"remove_genes_in_repeats"} =~ m/([\d\.]+)/;
print OUT "The number of gene models whose CDS length < $1 and cannot find ortholog in Pfam validation is: $gene_num4\n" if $config{"remove_short_genes"} =~ m/([\d\.]+)/;



print STDERR "\n============================================\n";
print STDERR "GETA complete successfully! " . "(" . (localtime) . ")" . "\n";
