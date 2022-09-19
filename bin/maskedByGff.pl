#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    perl $0 file.gff genome.fasta > out.masked.fasta

    --mask_type <string>    default: hardmaskN
    设置重复序列屏蔽的类型，有三个值可以选择：softmask(将原本需要硬隐蔽的重复序列字符由大写字母换成小写字母)、hardmaskX(将原本需要硬隐蔽的复序列字符替换成X)和hardmaskN(将原本需要硬隐蔽的重复序列字符替换成N)。

    --mask_all
    加入该参数后，则对GFF文件中的所有重复序列进行屏蔽。默认设置：不对Low_complexity、Simple_repeat、Satellite和Tandem_repeat类型的重复序列进行屏蔽，因为这种类型的重复序列长度较短，出现在基因区域的可能性较高；对Unknown和Other类型的重复序列进行软屏蔽；对其它类型的转座子序列进行硬屏蔽。

    --min_coverge_ratio <float>    default: 0.6
    设置最小覆盖率阈值。当需要进行硬屏蔽时，若匹配区域对目标重复序列的覆盖率>=此值，才继续进行硬屏蔽，否则修改为软屏蔽。程序读取GFF3文件最后一列中Ratio标签的值作为覆盖率；若GFF3文件中没有Ratio标签，则认为其Ratio的值为1。本参数的值设定范围为0~1。本参数作用：仅对覆盖率较高而比较明确的转座子序列进行硬屏蔽。

USAGE
if (@ARGV==0){die $usage}

my ($mask_type, $mask_all, $min_coverge_ratio);
GetOptions(
    "mask_type:s" => \$mask_type,
    "mask_all" => \$mask_all,
    "min_coverge_ratio:f" => \$min_coverge_ratio,
);
$mask_type ||= "hardmaskN";
$min_coverge_ratio ||= 0.6;

open GFF, '<', $ARGV[0];
open FASTA, '<', $ARGV[1];

my ($fasta_head, @fasta_head, %fasta);
while (<FASTA>) {
    chomp;
    if (/^>(\S+)/) {
        $fasta_head = $1;
        push @fasta_head, $fasta_head;
    }else {
        $fasta{$fasta_head} .= $_;
    }
}

while (<GFF>) {
    unless ($mask_all) {
        next if m/Name=Low_complexity/;
        next if m/Name=Simple_repeat/;
        #next if m/Name=.*RNA/;
        next if m/Name=Satellite/;
        next if m/Name=Tandem_repeat/;
    }
    my $ratio = 1;
    $ratio = $1 if m/Ratio=([^;\s]+)/i;

    if (/(\S+)\t.+\t.+\t(\d+)\t(\d+)\t/) {
        my $id = $1;
        my $start = $2 - 1;
        my $length = $3 - $start;
        if (m/Name=Unknown/) {
            substr($fasta{$id},$start,$length) =~ tr/ATCGN/atcgn/;
        }
        elsif (m/Name=Unspecified/) {
            substr($fasta{$id},$start,$length) =~ tr/ATCGN/atcgn/;
        }
        elsif (m/Name=Other/) {
            substr($fasta{$id},$start,$length) =~ tr/ATCGN/atcgn/;
        }
        elsif ($mask_type eq "softmask") {
            substr($fasta{$id},$start,$length) =~ tr/ATCGN/atcgn/;
        }
        elsif ($mask_type eq "hardmaskX") {
            if ( $ratio >= $min_coverge_ratio ) {
                substr($fasta{$id},$start,$length) =~ tr/ATCGNatcgn/XXXXXxxxxx/;
            }
            else {
                substr($fasta{$id},$start,$length) =~ tr/ATCGN/atcgn/;
            }
        }
        elsif ($mask_type eq "hardmaskN") {
            if ( $ratio >= $min_coverge_ratio ) {
                substr($fasta{$id},$start,$length) =~ tr/ATCGXatcgx/NNNNNnnnnn/;
            }
            else {
                substr($fasta{$id},$start,$length) =~ tr/ATCGN/atcgn/;
            }
        }
        else {
            die "The right mask type should be: softmask or hardmaskX or hardmaskN!\n";
        }
    }
}

foreach (@fasta_head) {
    my $seq = $fasta{$_};
    $seq =~ s/(\w{60})/$1\n/g;
    chomp($seq);
    print ">$_\n$seq\n";
}
