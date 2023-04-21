#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Cwd qw/abs_path getcwd cwd/;
use File::Basename;

my $usage = <<USAGE;
Usage:
    $0 [options] homolog_proteins.fasta genome_seq.fasta > out.gff

    --cpu <int>    default: 8
    设置diamond程序使用线程数，genewise或gth命令运行的并行数。

    --coverage_ratio <float>    default: 0.4
    对diamond blastx的结果进行过滤，要求对同源蛋白的覆盖率不小于该值。

    --evalue <float>    default: 1e-9
    对diamond blastx的结果进行过滤，要求对同源蛋白的evalue值不大于该值。

    --method <string>    default: gth
    设置使用同源蛋白进行基因预测的方法，可设置值为gth或genewise。一般情况下，genewise进行基因预测的耗时是gth的五倍。

    --tmp_dir <string>    default: tmp_\$date\$pid
    程序运行时临时文件夹名称。

    --help    default: None
    display this help and exit.

USAGE
if (@ARGV==0){die $usage}

my ($cpu, $coverage_ratio, $evalue, $method, $tmp_dir, $help_flag);
GetOptions(
    "cpu:i" => \$cpu,
    "coverage_ratio:f" => \$coverage_ratio,
    "evalue:f" => \$evalue,
    "method:s" => \$method,
    "tmp_dir:s" => \$tmp_dir,
    "help" => \$help_flag,
);
$cpu ||= 8;
$coverage_ratio ||= 0.4;
$evalue ||= 1e-9;
$method ||= "gth";
unless ($method eq "gth" or $method eq "genewise") {
    die "--method was set wrong, it should be set to gth or genewise. \n";
}
my $date = `date +%Y%m%d%H%M%S`; chomp($date);
$tmp_dir ||= "tmp_$date$$";
$tmp_dir = abs_path($tmp_dir);
mkdir $tmp_dir unless -e $tmp_dir;

if ( $help_flag ) { die $usage }

my $input_protein = abs_path($ARGV[0]);
my $input_genome = abs_path($ARGV[1]);

my $pwd = `pwd`; print STDERR "##########\nPWD (Current Directory): $pwd";
print STDERR (localtime) . "CMD (Main Program): $0 " . join(" ", @ARGV) . "\n##########\n\n";


# 1. 运行diamond分析，寻找蛋白序列和参考基因组的匹配位点。
print STDERR "1. align genome sequences against homolog proteins through diamond software.\n";
# 1.1 运行 diamond makedb 命令，将蛋白序列做成数据库。
print STDERR "1.1 run diamond makedb, get a protein database.\n";
chdir $tmp_dir;
my $pwd = `pwd`; print STDERR "PWD: $pwd";
my $cmdString = "diamond makedb --threads $cpu --db homolog --in $input_protein 2> 1.1.diamond_makedb.log";
unless (-e "1.1.diamond_makedb.ok") {
	print STDERR (localtime) . ": CMD: $cmdString\n";
	system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
	open OUT, ">", "1.1.diamond_makedb.ok" or die $!; close OUT;
}
else {
	print STDERR "CMD(Skipped): $cmdString\n";
}

# 1.2 运行 diamond blastx 命令将基因组序列比对到蛋白序列数据库上。注意 --max-target-seqs 参数设置为蛋白序列的数量。
print STDERR "1.2 split genome sequences, and get diamond blastx commands.\n";
unless ( -e "1.2.split_genome_seqs_and_get_cmds.ok" ) {
	# 获取蛋白序列数量
	open IN, $input_protein or die "Can not open file $input_protein, $!";
	my $protein_number = 0;
	while (<IN>) {
		$protein_number ++ if m/^>/;
	}
	close IN;
	# 将基因组序列分割成单独的序列，得到diamond blastx命令
	open IN, $input_genome or die "Can not open file $input_genome, !";
	my (%seq, $seq_ID);
	while (<IN>) {
		if (m/^>(\S+)/) {
			$seq_ID = $1;
		}
		else {
			chomp;
			$seq{$seq_ID} .= $_;
		}
	}
	close IN;
	mkdir "split_genome_seqs" unless -e "split_genome_seqs";
	open CMD1, ">", "split_genome_seqs/command.diamond_blastx.list" or die "Can not create file split_genome_seq/command.diamond_blastx.list, $!";
	open CMD2, ">", "split_genome_seqs/command.parsing_blast_result.list" or die "Can not create file split_genome_seq/command.parsing_blast_result.list, $!";
	foreach ( keys %seq ) {
		open OUT, ">", "split_genome_seqs/$_.fa" or die "Can not create file split_genome_seq/$_.fa, $!";
		print OUT ">$_\n$seq{$_}\n";
		close OUT;
		print CMD1 "diamond blastx --db homolog --query split_genome_seqs/$_.fa --out split_genome_seqs/$_.xml --outfmt 5 --sensitive --max-target-seqs $protein_number --evalue 1e-3 --id 10 2> split_genome_seqs/$_.log\n";
		print CMD2 "parsing_blast_result.pl --type xml --no-header --max-hit-num $protein_number --evalue $evalue --subject-coverage $coverage_ratio --query-coverage 0 split_genome_seqs/$_.xml > split_genome_seqs/$_.tab\n";
	}
	close CMD1; close CMD2;
	open OUT, ">", "1.2.split_genome_seqs_and_get_cmds.ok" or die $!; close OUT;
}
else {
	print STDERR "This step was skipped, for the file 1.2.split_genome_seqs_and_get_cmds.ok exists.\n";
}

# 1.3 并行化运行 diamond blastx 命令
# 按内存余量确定 diamond 的并行化数量。按单个 dianmond 命令峰值消耗 60 Gb 内存。
my $MemAvailable = &get_MemAvailable();
my $paraFly_CPU = 1;
$paraFly_CPU = $cpu / 4 if $paraFly_CPU < ($cpu / 4);
$paraFly_CPU = $MemAvailable / 60000000 if $paraFly_CPU > ($MemAvailable / 60000000);
$paraFly_CPU = 1 if $paraFly_CPU < 1;
my $cmdString = "ParaFly -c split_genome_seqs/command.diamond_blastx.list -CPU $paraFly_CPU 2> 1.3.ParaFly_diamond_blastx.log";
unless ( -e "1.3.ParaFly_diamond_blastx.ok" ) {
	print STDERR (localtime) . ": CMD: $cmdString\n";
	if ( system("$cmdString") != 0 ) {
		$cmdString = "ParaFly -c split_genome_seqs/command.diamond_blastx.list -CPU 1 2>> 1.3.ParaFly_diamond_blastx.log";
		print STDERR (localtime) . ": CMD: $cmdString\n";
		system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
	}
	open OUT, ">", "1.3.ParaFly_diamond_blastx.ok" or die $!; close OUT;
}
else {
	print STDERR "CMD(Skipped): $cmdString\n";
}

# 1.4 解析 xml 文件，得到基因组区域及其匹配的氨基酸序列名称。
$cmdString = "ParaFly -c split_genome_seqs/command.parsing_blast_result.list -CPU $cpu 2> 1.4.ParaFly_parsing_blast_result.log";
unless ( -e "1.4.ParaFly_parsing_blast_result.ok" ) {
	print STDERR (localtime) . ": CMD: $cmdString\n";
	if ( system("$cmdString") != 0 ) {
		$cmdString = "ParaFly -c split_genome_seqs/command.parsing_blast_result.list -CPU 1 2> 1.4.ParaFly_parsing_blast_result.log";
		print STDERR (localtime) . ": CMD: $cmdString\n";
		system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
	}
	open OUT, ">", "1.4.ParaFly_parsing_blast_result.ok" or die $!; close OUT;
}
else {
	print STDERR "CMD(Skipped): $cmdString\n";
}

# 2. 解析diamond比对结果，
# 2.1 根据同源蛋白的比对结果，推测出物种的基因长度信息。
print STDERR "2

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
