#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Cwd qw/abs_path getcwd cwd/;
use File::Basename;

my $usage = <<USAGE;
Usage:
    $0 [options] seq1.fasta seq2.fasta ... > out.fasta
    
    The program reads the sequence in one or more fasta files sequentially, eliminates any redundant sequences, and reports the result. That is, when two sequences are identical, only the sequence read first is saved. The standard output contains the FASTA result, while the standard error output provides the number of redundant sequences that have been eliminated.

    --label <string>    default: None
    Set the labels for the input FASTA files. The value of this parameter is a list of comma-separated strings that correspond to the FASTA files entered afterward. After this parameter is introduced, the output sequence name is changed, and additional label information appears after the original sequence name to aid in the differentiation of sequence information from various species. If the number of labels in this parameter differs from the number of following input FASTA files, the program will still function normally. If there are fewer labels, some FASTA files do not add labels to the sequence name; if there are more labels, the labels in the last position do not take effect.

    --help    default: None
    display this help and exit.

USAGE
if (@ARGV==0){die $usage}

my ($label, $help_flag);
GetOptions(
    "label:s" => \$label,
    "help" => \$help_flag,
);

if ( $help_flag ) { die $usage }

# 分析标签
my @label;
if ( $label ) {
    @label = split /,/, $label;
}

my (%seq, $seq, $seq_ID);
# 读取输入的每个FASTA文件
foreach my $file ( @ARGV ) {
    my $species_name;
    $species_name = shift @label if @label;

    open IN, $file or die "Error: Can not open file $file, $!";
    # FASTA文件第一行，提取序列名，尾部追加标签信息
    $_ = <IN>;
    $seq_ID = $1 if m/^>(\S+)/;
    $seq_ID .= $species_name;
    while ( <IN> ) {
        # 每当读取到新的一条序列名称时，输出上一条序列。
        if (m/^>(\S+)/) {
            if ( ! exists $seq{$seq} ) {
                print ">$seq_ID\n$seq\n";
            }
            $seq{$seq} ++;
            $seq = "";
            $seq_ID = $1;
            $seq_ID .= $species_name;
        }
        else {
            chomp;
            $seq .= $_;
        }
    }
    if ( ! exists $seq{$seq} ) {
        print ">$seq_ID\n$seq\n";
    }
    $seq{$seq} ++;
    $seq = "";
    close IN;
}

my $num_remove = 0;
foreach ( keys %seq ) {
    $num_remove += ($seq{$_} - 1);
}
print STDERR "There are $num_remove redundant sequences were deleted.\n";
