#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    perl $0 file.gff genome.fasta > out.masked.fasta

	--mask_type <string>    default: hardmaskN
	    设置重复序列屏蔽的类型，有三个值可以选择：softmask(将重复序列字符由大写字母换成小写字母)、hardmaskX(将重复序列字符替换成X)和hardmaskN(将重复序列字符替换成N)。
	--mask_all
	    加入该参数后，则对GFF文件中的所有重复序列进行屏蔽。默认设置是不对Low_complexity和Simple_repeat类型的重复序列进行屏蔽。因为这种类型的重复序列长度较短，出现在基因区域的可能性较高。

USAGE
if (@ARGV==0){die $usage}

my ($mask_type, $mask_all);
GetOptions(
	"mask_type:s" => \$mask_type,
	"mask_all" => \$mask_all,
);
$mask_type ||= "hardmaskN";

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
		next if m/Low_complexity/;
		next if m/Simple_repeat/;
	}
    if (/(\S+)\t.+\t.+\t(\d+)\t(\d+)\t/) {
        my $id = $1;
        my $start = $2 - 1;
        my $length = $3 - $start;
        if ($mask_type eq "softmask") {
            substr($fasta{$id},$start,$length) =~ tr/ATCGN/atcgn/;
        }
        elsif ($mask_type eq "hardmaskX") {
            substr($fasta{$id},$start,$length) =~ tr/ATCGNatcgn/XXXXXxxxxx/;
        }
        elsif ($mask_type eq "hardmaskN") {
            substr($fasta{$id},$start,$length) =~ tr/ATCGXatcgx/NNNNNnnnnn/;
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
