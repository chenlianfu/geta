#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    $0 seq.fasta id.list > target.fasta

    默认设置下，分析id.list文件第一列表示序列ID的信息，并从seq.fasta中提取这些序列，输出FASTA格式。

    --reverse
    添加该参数，则提取在id.list中不存在的剩余序列。

USAGE
if (@ARGV==0){die $usage}

my $reverse;
GetOptions(
    "reverse!" => \$reverse
);

open IN, $ARGV[0] or die $!;
my (%seq, $seq_id);
while (<IN>) {
    chomp;
    if (m/^>(\S+)/) { $seq_id = $1; }
    else { $seq{$seq_id} .= $_; }
}
close IN;

open IN, $ARGV[1] or die $!;
my (@list, %list);
while (<IN>) {
    if (m/^(\S+)/ && ! exists $list{$1}) {
        $list{$1} = 1;
        push @list, $1 if exists $seq{$1};
    }
}
close IN;

my @list_number = keys %list;
my $list_number = @list_number;
my $list_num = @list;
print STDERR "Their are $list_number sequences ID in file $ARGV[1], $list_num of which can be found in file $ARGV[0].\n";

my $number;
if ($reverse) {
    foreach (sort keys %seq) {
        unless (exists $list{$_}) {
            print ">$_\n$seq{$_}\n";
            $number ++;
        }
    }
}
else {
    foreach (@list) {
        $number ++;
        print ">$_\n$seq{$_}\n";
    }
}
print STDERR "$number sequences were output.\n";
