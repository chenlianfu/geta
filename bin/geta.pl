#!/usr/bin/perl
use strict;
use Getopt::Long;
use File::Basename;

my $usage = <<USAGE;
Usage:
    perl $0 [options]

For example:
    perl $0 --out_prefix out -1 reads.1 fastq -2 reads.2.fastq --protein homolog.fasta --cpu 80 --hisat2 " --min-intronlen 20 --max-intronlen 20000 --rna-strandness RF" --strand_specific --sam2transfrag " --fraction --min_expressed_base_depth 2 --max_expressed_base_depth 50 --min_junction_depth 2 --max_junction_depth 50" --species oryza_sativa_20171120 --pfam_db /opt/biosoft/hmmer-3.1b2/Pfam-AB.hmm --gene_prefix OS01Gene --genome genome.fasta

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

USAGE
if (@ARGV==0){die $usage}

my ($RM_species, $genome, $out_prefix, $pe1, $pe2, $single_end, $protein, $cpu, $trimmomatic, $hisat2, $strand_specific, $sam2transfrag, $ORF2bestGeneModels, $species, $pfam_db, $gene_prefix, $cmdString);
GetOptions(
    "RM_species:s" => \$RM_species,
    "genome:s" => \$genome,
    "out_prefix:s" => \$out_prefix,
    "1:s" => \$pe1,
    "2:s" => \$pe2,
    "S:s" => \$single_end,
    "protein:s" => \$protein,
    "cpu:s" => \$cpu,
    "trimmomatic:s" => \$trimmomatic,
    "hisat2:s" => \$hisat2,
    "strand_specific!" => \$strand_specific,
    "sam2transfrag:s" => \$sam2transfrag,
    "ORF2bestGeneModels:s" => \$ORF2bestGeneModels,
    "species:s" => \$species,
    "pfam_db:s" => \$pfam_db,
    "gene_prefix:s" => \$gene_prefix,
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
if ($software_info =~ m/java version \"(1.(\d).*?)\"/) {
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
die "--RM_species shoud be set!" unless $RM_species;
my $pwd = `pwd`;
chomp($pwd);
die "No genome fasta input\n" unless $genome;
$genome =~ s/^/$pwd\// unless $genome =~ m/^\//;
die "No homolog fasta input\n" unless $protein;
$protein  =~ s/^/$pwd\// unless $protein =~ m/^\//;
die "No Augustus species provided\n" unless $species;
unless (($pe1 && $pe2) or $single_end) {
    die "No RNA-Seq short reads as input\n";
}
if ($pe1 && $pe2) {
    $pe1 =~ s/^/$pwd\// unless $pe1 =~ m/^\//;
    $pe2 =~ s/^/$pwd\// unless $pe2 =~ m/^\//;
}
if ($single_end) {
    $single_end =~ s/^/$pwd\// unless $single_end =~ m/^\//;
}

$out_prefix ||= "out";
$trimmomatic ||= "TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 TOPHRED33";
$hisat2 ||= " --min-intronlen 20 --max-intronlen 20000 --rna-strandness unstranded --dta --score-min L,0.0,-0.4";
$sam2transfrag ||= " --fraction 0.05 --min_expressed_base_depth 2 --max_expressed_base_depth 50 --min_junction_depth 2 --max_junction_depth 50";
$ORF2bestGeneModels ||= " --min_cds_num 3 --min_cds_length 900 --min_cds_exon_ratio 0.60 --intron_length_fractile 0.95 --cds_length_fractile 0.95";
$pfam_db ||= "/opt/biosoft/hmmer-3.1b2/Pfam-AB.hmm";

my $dirname = dirname($0);
$dirname =~ s/bin$//;

# 生成临时文件夹
mkdir "$out_prefix.tmp" unless -e "$out_prefix.tmp";
chdir "$out_prefix.tmp";

# Step 0: RepeatMasker and RepeatModeler
print STDERR "\n============================================\n";
print STDERR "Step 0: RepeatMasker and RepeatModeler\n";
mkdir "0.RepeatMasker" unless -e "0.RepeatMasker";
unless (-e "0.RepeatMasker.ok") {
    chdir "0.RepeatMasker";
    
    # 进行RepeatMasker分析
    mkdir "repeatMasker" unless -e "repeatMasker";
    $cmdString = "RepeatMasker -e ncbi -pa $cpu -species $RM_species -dir repeatMasker/ -gff $genome &> repeatmasker.log";
    unless (-e "RepeatMasker.ok") {
        print STDERR "CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "RepeatMasker.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    # 进行RepeatModeler分析
    mkdir "repeatModeler" unless -e "repeatModeler";
    chdir "repeatModeler";
    $cmdString = "BuildDatabase -name species -engine ncbi $genome";
    unless (-e "BuildDatabase.ok") {
        print STDERR "CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "BuildDatabase.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }
    $cmdString = "RepeatModeler -engine ncbi -pa $cpu -database species";
    unless (-e "RepeatModeler.ok") {
        print STDERR "CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "RepeatModeler.ok" or die $!; close OUT;
    }
    $cmdString = "RepeatMasker -e ncbi -pa $cpu -lib RM_\*/\*.classified -dir ./ -gff $genome &> repeatmasker.log";
    unless (-e "RepeatMasker.ok") {
        print STDERR "CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "RepeatMasker.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }
    chdir "../";

    # 合并RepeatMasker和RepeatModeler的结果
    $cmdString = "$dirname/bin/merge_repeatMasker_out.pl repeatMasker/genome.fasta.out repeatModeler/genome.fasta.out > genome.repeat.stats";
    print STDERR "CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    $cmdString = "$dirname/bin/maskedByGff.pl genome.repeat.gff3 $genome > genome.masked.fasta";
    print STDERR "CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    chdir "../";
    open OUT, ">", "0.RepeatMasker.ok" or die $!; close OUT;
}
else {
    print STDERR "Skip Step 0 for the file 0.RepeatMasker.ok exists\n";
}

# Step 1: Trimmomatic
print STDERR "\n============================================\n";
print STDERR "Step 1: Trimmomatic\n";
mkdir "1.trimmomatic" unless -e "1.trimmomatic";
unless (-e "1.trimmomatic.ok") {
    chdir "1.trimmomatic";
    if ($pe1 && $pe2) {
        $cmdString = "java -jar $dirname/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads $cpu $pe1 $pe2 reads.1.fastq reads.1.unpaired.fastq reads.2.fastq reads.2.unpaired.fastq ILLUMINACLIP:$dirname/Trimmomatic-0.36/adapters/$trimmomatic &> trimmomatic.pe.log";
        print STDERR "CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    }
    if ($single_end) {
        $cmdString = "java -jar $dirname/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads $cpu $single_end reads.fastq ILLUMINACLIP:$dirname/Trimmomatic-0.36/adapters/$trimmomatic &> trimmomatic.single.log";
        print STDERR "CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    }
    chdir "../";
    open OUT, ">", "1.trimmomatic.ok" or die $!; close OUT;
}
else {
    print STDERR "Skip Step 1 for the file 1.trimmomatic.ok exists\n";
}

# Step 2: HISAT2
print STDERR "\n============================================\n";
print STDERR "Step 2: HISAT2\n";
mkdir "2.hisat2" unless -e "2.hisat2";
unless (-e "2.hisat2.ok") {
    chdir "2.hisat2";

    # 构建基因组hisat2索引数据库
    $cmdString = "hisat2-build $genome genome &> hisat2-build.log\n";
    unless (-e "hisat2-build.ok") {
        print STDERR "CMD: $cmdString\n";
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
    $cmdString = "hisat2 -x genome -p $cpu $hisat2 $input -S hisat2.sam 2> hisat2.log";
    unless (-e "hisat2.ok") {
        print STDERR "CMD: $cmdString\n";
        print STDERR "Warning: The HISAT2 parameter --rna-strandness may not set corretly !" if ($hisat2 =~ m/unstranded/ && $strand_specific);
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "hisat2.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    $cmdString = "samtools sort -@ $cpu -o hisat2.sorted.sam -O sam hisat2.sam";
    print STDERR "CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

    chdir "../";
    open OUT, ">", "2.hisat2.ok" or die $!; close OUT;
}
else {
    print STDERR "Skip Step 2 for the file 2.hisat2.ok exists\n";
}

# Step 3: Transcript
print STDERR "\n============================================\n";
print STDERR "Step 3: Transcript\n";
mkdir "3.transcript" unless -e "3.transcript";
unless (-e "3.transcript.ok") {
    chdir "3.transcript";

    # 将SAM比对将诶过分成一个个比对区域
    $cmdString = "$dirname/bin/split_sam_from_non_aligned_region ../2.hisat2/hisat2.sorted.sam splited_sam_out 10 > splited_sam_files.list";
    unless (-e "split_sam_from_non_aligned_region.ok") {
        print STDERR "CMD: $cmdString\n";
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
        print CMD "$dirname/bin/sam2transfrag $sam2transfrag $no_strand_specific --intron_info_out $_.intron $_.sam > $_.gtf\n";
    }
    close CMD;
    close IN;

    # 批量并行化进行transcripts计算
    $cmdString = "ParaFly -c command.sam2transfrag.list -CPU $cpu &> /dev/null";
    print STDERR "CMD: $cmdString\n";
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
        print STDERR "CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "transfragDump.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    # 对transcripts序列使用Transdecoder进行ORF分析
    unless (-e "TransDecoder.ok") {
        $cmdString = "$dirname/TransDecoder-2.0.1/TransDecoder.LongOrfs -t transfrag.strand.fasta -S &> /dev/null";
        print STDERR "CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        $cmdString = "$dirname/TransDecoder-2.0.1/TransDecoder.Predict -t transfrag.strand.fasta &> /dev/null";
        print STDERR "CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        $cmdString = "cp transfrag.strand.fasta.transdecoder.gff3 transfrag.transdecoder.gff3";
        print STDERR "CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        unless ($strand_specific) {
            $cmdString = "$dirname/TransDecoder-2.0.1/TransDecoder.LongOrfs -t transfrag.noStrand.fasta &> /dev/null";
            print STDERR "CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            $cmdString = "$dirname/TransDecoder-2.0.1/TransDecoder.Predict -t transfrag.noStrand.fasta --train transfrag.strand.fasta.transdecoder_dir/longest_orfs.cds.top_500_longest &> /dev/null";
            print STDERR "CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            $cmdString = "cat transfrag.noStrand.fasta.transdecoder.gff3 >> transfrag.transdecoder.gff3";
            print STDERR "CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        }
        open OUT, ">", "TransDecoder.ok" or die $!; close OUT;
    }
    else {
        $cmdString = "$dirname/TransDecoder-2.0.1/TransDecoder.LongOrfs -t transfrag.strand.fasta -S &> /dev/null";
        print STDERR "CMD(Skipped): $cmdString\n";
        $cmdString = "$dirname/TransDecoder-2.0.1/TransDecoder.Predict -t transfrag.strand.fasta &> /dev/null";
        print STDERR "CMD(Skipped): $cmdString\n";
        $cmdString = "cp transfrag.strand.fasta.transdecoder.gff3 transfrag.transdecoder.gff3";
        print STDERR "CMD(Skipped): $cmdString\n";
        unless ($strand_specific) {
            $cmdString = "$dirname/TransDecoder-2.0.1/TransDecoder.LongOrfs -t transfrag.noStrand.fasta &> /dev/null";
            print STDERR "CMD(Skipped): $cmdString\n";
            $cmdString = "$dirname/TransDecoder-2.0.1/TransDecoder.Predict -t transfrag.noStrand.fasta --train transfrag.strand.fasta.transdecoder_dir/longest_orfs.cds.top_500_longest &> /dev/null";
            print STDERR "CMD(Skipped): $cmdString\n";
            $cmdString = "cat transfrag.noStrand.fasta.transdecoder.gff3 >> transfrag.transdecoder.gff3";
            print STDERR "CMD(Skipped): $cmdString\n";
        }
    }

    # 将transcripts的ORF预测结果映射到基因组序列上，得到transcripts的基因预测结果： transfrag.genome.gff3
    $cmdString = "$dirname/bin/transdecoder2ORF --out_protein proteins.fasta transfrag.gtf transfrag.transdecoder.gff3 $genome > transfrag.genome.gff3";
    unless (-e "transdecoder2ORF.ok") {
        print STDERR "CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "transdecoder2ORF.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

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
        print STDERR "CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "ORF2bestGeneModels.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    # 过滤一些基因模型，以保证两两基因模型之间在protein序列上Identity要低于80%
    $cmdString = "$dirname/bin/bestGeneModels2lowIdentity best_candidates.gff3 proteins.fasta $cpu 0.8 > best_candidates.lowIdentity.gff3 2> bestGeneModels2lowIdentity.log";
    print STDERR "CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

    chdir "../";
    open OUT, ">", "3.transcript.ok" or die $!; close OUT;
}
else {
    print STDERR "Skip Step 3 for the file 3.transcript.ok exists\n";
}

# Step 4: Homolog
print STDERR "\n============================================\n";
print STDERR "Step 4: Homolog\n";
mkdir "4.homolog" unless -e "4.homolog";
unless (-e "4.homolog.ok") {
    chdir "4.homolog";

    $cmdString = "$dirname/bin/homolog_genewise.pl $protein ../0.RepeatMasker/genome.masked.fasta $cpu 0.4 1e-9 &> homolog_genewise.log";
    unless (-e "homolog_genewise.ok") {
        print STDERR "CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "homolog_genewise.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    $cmdString = "$dirname/bin/genewise_filter.pl --pfam_db $pfam_db genewise.gff $genome 25 300 2 1e-6 0.30 $cpu > genewise.filter.gff";
    unless (-e "genewise_filter.ok") {
        print STDERR "CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "genewise_filter.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    $cmdString = "$dirname/bin/genewise2GFF3.pl genewise.filter.gff > genewise.filter.gff3";
    print STDERR "CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

    chdir "../";
    open OUT, ">", "4.homolog.ok" or die $!; close OUT;
}
else {
    print STDERR "Skip Step 4 for the file 4.homolog.ok exists\n";
}

# Step 5: Augustus gene prediction
print STDERR "\n============================================\n";
print STDERR "Step 5: Augustus/HMM Trainning\n";
mkdir "5.augustus" unless -e "5.augustus";
unless (-e "5.augustus.ok") {
    chdir "5.augustus";

    # Augustus HMM Training
    mkdir "training" unless -e "training";
    unless (-e "training.ok") {
        chdir "training";

        #my $flanking_length = `tail -n 1 ../../3.transcript/transfrag.stats | cut -f 2`;
        my (%gene_info, @flanking_length, $flanking_length);
        open IN, "../../3.transcript/transfrag.genome.gff3" or die $!;
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
                foreach (@region) {
                    my ($aa, $bb) = split /\t/;
                    my $distance = $aa - $end - 1;
                    push @flanking_length, $distance if $distance >= 50;
                    ($start, $end) = ($aa, $bb);
                }
            }
        }
        @flanking_length = sort {$a <=> $b} @flanking_length;
        $flanking_length = int($flanking_length[@flanking_length/2] / 8);
        $flanking_length = 2500 if $flanking_length >= 2500;
        $cmdString = "$dirname/bin/BGM2AT --flanking_length $flanking_length ../../3.transcript/best_candidates.lowIdentity.gff3 $genome $species &> BGM2AT.log";
        print STDERR "CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

        chdir "../";
        open OUT, ">", "training.ok" or die $!;
    }
    else {
        print STDERR "Skip Augustus training for file training.ok exists\n";
    }

    # Augustus Hint Preparing
    $cmdString = "$dirname/bin/prepareAugusutusHints ../3.transcript/transfrag.genome.gff3 ../3.transcript/intron.txt ../4.homolog/genewise.filter.gff3 > hints.gff; sort -k1,1 -k4n,4 hints.gff | uniq > 11; mv 11 hints.gff";
    unless (-e "prepareAugusutusHints.ok") {
        print STDERR "CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "prepareAugusutusHints.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    # Augustus gene prediction
    $cmdString = "$dirname/bin/paraAugusutusWithHints --species $species --cpu $cpu ../0.RepeatMasker/genome.masked.fasta hints.gff > augustus.gff3";
    print STDERR "CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

    open IN, "augustus.gff3" or die $!;
    my @augustus_gff3_line;
    while (<IN>) {
        push @augustus_gff3_line, $_;
    }
    close IN;
    open OUT, ">", "augustus.gff3" or die $!;
    foreach (@augustus_gff3_line) {
        s/\ttranscript\t/\tmRNA\t/;
        print OUT;
    }
    close OUT;

    chdir "../";
    open OUT, ">", "5.augustus.ok" or die $!; close OUT;
}
else {
    print STDERR "Skip Step 5 for the file 5.augustus.ok exists\n";
}

# Step 6: CombineGeneModels
print STDERR "\n============================================\n";
print STDERR "Step 6: CombineGeneModels\n";
mkdir "6.combineGeneModels" unless -e "6.combineGeneModels";
unless (-e "6.combineGeneModels.ok") {
    chdir "6.combineGeneModels";

    $cmdString = "$dirname/bin/paraCombineGeneModels --cpu $cpu ../5.augustus/augustus.gff3 ../3.transcript/transfrag.genome.gff3 ../4.homolog/genewise.filter.gff3 ../5.augustus/hints.gff &> /dev/null";
    unless (-e "paraCombineGeneModels.ok") {
        print STDERR "CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "paraCombineGeneModels.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    $cmdString = "$dirname/bin/PfamValidateABinitio --cpu $cpu --pfam_db $pfam_db combine.2.gff3 $genome > combine.filter.gff3 2> PfamValidateABinitio.log";
    unless (-e "PfamValidateABinitio.ok") {
        print STDERR "CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "PfamValidateABinitio.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    $cmdString = "$dirname/bin/GFF3Clear --gene_prefix $gene_prefix --genome $genome combine.1.gff3 combine.filter.gff3 > genome.gff3";
    print STDERR "CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

    chdir "../";
    open OUT, ">", "6.combineGeneModels.ok" or die $!; close OUT;;
}
else {
    print STDERR "Skip Step 6 for the file 6.combineGeneModels.ok exists\n";
}

# Step 7: OutPut
print STDERR "\n============================================\n";
print STDERR "Step 7: OutPut\n";
chdir "../";

$cmdString = "cp $out_prefix.tmp/6.combineGeneModels/genome.gff3 $out_prefix.gff3";
print STDERR "CMD: $cmdString\n";
system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

$cmdString = "$dirname/bin/gff3ToGtf.pl $genome $out_prefix.gff3 > $out_prefix.gtf 2> /dev/null";
print STDERR "CMD: $cmdString\n";
system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

$cmdString = "$dirname/bin/eukaryotic_gene_model_statistics.pl $out_prefix.gtf $genome $out_prefix &> eukaryotic_gene_model_statistics.stats";
print STDERR "CMD: $cmdString\n";
system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

print STDERR "\n============================================\n";
print STDERR "GETA complete successfully!\n";
