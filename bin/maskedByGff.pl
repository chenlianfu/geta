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

    --out_repeat_percentage_of_each_sequence <string>    default: None
    输出每条基因组序列的所有重复序列的百分比。

USAGE
if (@ARGV==0){die $usage}

my ($mask_type, $mask_all, $min_coverge_ratio, $out_repeat_percentage_of_each_sequence);
GetOptions(
    "mask_type:s" => \$mask_type,
    "mask_all" => \$mask_all,
    "min_coverge_ratio:f" => \$min_coverge_ratio,
    "out_repeat_percentage_of_each_sequence:s" => \$out_repeat_percentage_of_each_sequence,
);
$mask_type ||= "hardmaskN";
$min_coverge_ratio ||= 0.6;

open GFF, '<', $ARGV[0];
open FASTA, '<', $ARGV[1];

my ($fasta_head, @fasta_head, %fasta, %seq_length);
while (<FASTA>) {
    chomp;
    if (/^>(\S+)/) {
        $fasta_head = $1;
        push @fasta_head, $fasta_head;
    }else {
        $fasta{$fasta_head} .= $_;
        $seq_length{$fasta_head} += length($_);
    }
}

my %repeat_region;
while (<GFF>) {
    next if m/^#/;
    next if m/^\s*$/;
    @_ = split /\t/;
    $repeat_region{$_[0]}{"$_[3]\t$_[4]"} = 1;

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

if ( $out_repeat_percentage_of_each_sequence ) {
    open OUT, ">", $out_repeat_percentage_of_each_sequence or die "Can not create file $out_repeat_percentage_of_each_sequence, $!";
    foreach ( sort keys %repeat_region ) {
        my @region = keys %{$repeat_region{$_}};
        my $seq_length = $seq_length{$_};
        my $region_length = 0;
        $region_length = &match_length(@region);
        my $ratio = 0;
        $ratio = $region_length * 100 / $seq_length if $seq_length;
        $ratio = sprintf("%.2f", $ratio);
        print OUT "$_\t$ratio%\n";
    }
}

sub match_length {
    my @inter_sorted_site;
    foreach (@_) {
        my @aaa = $_ =~ m/(\d+)/g;
        @aaa = sort { $a <=> $b } @aaa;
        push @inter_sorted_site, "$aaa[0]\t$aaa[1]";
    }
    @inter_sorted_site = sort { $a <=> $b } @inter_sorted_site;

    my $out_site_number;
    my $former_region = shift @inter_sorted_site;
    my @aaa = $former_region =~ m/(\d+)/g;
    $out_site_number += ($aaa[1] - $aaa[0] + 1);
    foreach (@inter_sorted_site) {
        my @former_region = $former_region =~ m/(\d+)/g;
        my @present_region = $_ =~ m/(\d+)/g;

        if ($present_region[0] > $former_region[1]) {
            $out_site_number += ($present_region[1] - $present_region[0] + 1);
            $former_region = $_;
        }
        elsif ($present_region[1] > $former_region[1]) {
            $out_site_number += ($present_region[1] - $former_region[1]);
            $former_region = $_;
        }
        else {
            next
        }
    }
    return $out_site_number;
}
