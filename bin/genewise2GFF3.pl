#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    perl $0 [options]  genewise.gff 15 > genewise.gff3

	示例中，过滤掉得分低于15的genwise预测结果，同时将genwise结果转换成GFF3格式，同时修改Gene ID。

	--min_score <int>    default: 15
	设置最低得分阈值。

	--gene_prefix <string>    default: genewise
	设置基因ID的前缀。

USAGE
if (@ARGV==0){die $usage}

my ($min_score, $gene_prefix);
GetOptions(
	"min_score:i" => \$min_score,
	"gene_prefix:s" => \$gene_prefix,
);

$min_score ||= 15;
$gene_prefix ||= "genewise";

my (%info, %score, %sort1, %sort2, %sort3, %sort4);
open IN, $ARGV[0] or die $!;
while (<IN>) {
	if (m/\tmatch\t.*ID=([^;]+);Name=([^;]+)/) {
		my $id = $1;
		@_ = split /\t/;
		$score{$id} = $_[5];
		$sort1{$id} = $_[0];
		$sort2{$id} = $_[3];
		$sort3{$id} = $_[4];
		$sort4{$id} = $_[6];
		s/\tmatch\t/\tgene\t/;
		my $mRNA = $_;
		$mRNA =~ s/\tgene\t/\tmRNA\t/;
		$mRNA =~ s/ID=([^;]+);(.*)/ID=$1.mRNA;Parent=$1;$2/;
		s/(\d+)\t(\d+)\t(.*?)\t\-/$2\t$1\t$3\t\-/;
		$mRNA =~ s/(\d+)\t(\d+)\t(.*?)\t\-/$2\t$1\t$3\t\-/;
		$info{$id}{"gene"} = "$_$mRNA";
	}
	elsif (m/\tcds\t.*ID=([^;]+);Name=([^;]+)/) {
		my $id = $1;
		s/\tcds\t/\tCDS\t/;
		s/ID=([^;]+);.*/ID=$1.mRNA.CDS;Parent=$1.mRNA;/;
		s/(\d+)\t(\d+)\t(.*?)\t\-/$2\t$1\t$3\t\-/;
		$info{$id}{"cds"} .= $_;
		#s/CDS\t/exon\t/;
		#s/CDS;/exon;/;
		#s/\d\tID/\.\tID/;
		#$info{$id}{"cds"} .= $_;
	}
}
close IN;

my $total_number = keys %info;
$total_number = length $total_number;
my $number;
foreach my $gene_id (sort {$sort1{$a} cmp $sort1{$b} or $sort2{$a} <=> $sort2{$b} or $sort3{$a} <=> $sort3{$b} or $sort4{$a} cmp $sort4{$b}} keys %info) {
	if ($score{$gene_id} >= $min_score) {
		$number ++;
		my $out = "$info{$gene_id}{'gene'}$info{$gene_id}{'cds'}";
		my $id_out = $gene_prefix . 0 x ($total_number - length($number)) . $number;
		$out =~ s/$gene_id/$id_out/g;
		print "$out\n";
	}
}
