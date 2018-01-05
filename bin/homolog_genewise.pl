#!/usr/bin/perl
use strict;
use File::Basename;

my $usage = <<USAGE;
Usage:
    perl $0 homolog_proteins.fasta genome_seq.fasta CPUs coverage_ratio e-value

For example:
    perl $0 homolog_proteins.fasta genome_seq.fasta 24 0.4 1e-9

    程序分4步运行：
    1. 将genome_seq.fasta文件使用makeblastdb构建数据库。成功运行后生成文件makeblastdb.ok。
    2. 运行程序同目录下的blast.pl程序并行化(24个线程)调用tblastn将homolog_proteins.fasta比对到数据库（阈值为evalue<=1e-9且覆盖率>=0.4）。成功运行后生成文件blast.ok和结果文件out.tab。
    3. 运行程序同目录下的blast2geneRegion.pl对out.tab文件进行分析，找到gene region及其最优homolog。成功运行后生成文件blast2geneRegion.ok和结果文件homolog_gene_region.tab。
    4. 运行程序同目录下的geneRegion_genewise.pl提取基因组的gene region序列，并行化(24个线程)调用genewise和homolog进行比对，生成最终结果文件genewise.gff。

    .ok 后缀的文件用于表明其步骤运行完毕。若重新运行时，检测到 .ok 文件存在，则跳过相应的步骤。

USAGE
if (@ARGV==0){die $usage}

my $homologProteinFastaFile = $ARGV[0];
my $genomeFastaFile = $ARGV[1];
my $cpus = $ARGV[2];
my $path = dirname($0);

my $cmdString = "makeblastdb -in $genomeFastaFile -dbtype nucl -title genome -parse_seqids -out genome &> makeblastdb.log";
print "CMD 1/4: $cmdString\n";
unless (-e "makeblastdb.ok") {
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    open OUT, ">", "makeblastdb.ok" or die $!; close OUT;
    print "1. Successfully Make Blast Database of $genomeFastaFile\n\n";
}
else {
    print "1. Skip CMD 1/4 for file makeblastdb.ok exists\n\n";
}

$cmdString = "$path/blast.pl tblastn genome $homologProteinFastaFile 1e-5 $cpus out 6";
print "CMD 2/4: $cmdString\n";
unless (-e "blast.ok") {
    system("$cmdString") == 0 or die "failed to execute: $!\n";
    if ( -e "FailedCommands" ) { die "failed to execute: $cmdString\n"; }
    open OUT, ">", "blast.ok" or die $!; close OUT;
    print "2. Successfully align homolog to genome by tblastn\n\n";
}
else {
    print "2. Skip CMD 2/4: for file blast.ok exists\n\n";
}

$cmdString = "$path/blast2geneRegion.pl out.tab $homologProteinFastaFile $ARGV[3] $ARGV[4] > homolog_gene_region.tab";
print "CMD 3/4: $cmdString\n";
unless (-e "blast2geneRegion.ok") {
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    open OUT, ">", "blast2geneRegion.ok" or die $!; close OUT;
    print "3. Successfully get homolog gene region information\n\n";
}
else {
    print "3. Skip CMD 3/4: for file blast2geneRegion.ok  exists\n\n";
}

open IN, "homolog_gene_region.tab";
my ($gene_total_length, $gene_number, $gene_avg_length);
while (<IN>) {
    @_ = split /\t/;
    $gene_total_length += ($_[2] - $_[1] + 1);
    $gene_number ++;
}
$gene_avg_length = int($gene_total_length / $gene_number);
print "    $gene_number gene regions, with average gene length $gene_avg_length\n";

$cmdString = "$path/geneRegion_genewise.pl $genomeFastaFile $homologProteinFastaFile homolog_gene_region.tab $gene_avg_length $cpus > genewise.gff";
print "CMD 4/4: $cmdString\n";
system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
print "4. Successfully run genewise parally\n\n";
print "The final result file : genewise.gff\n\n";
