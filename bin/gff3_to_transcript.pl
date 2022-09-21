#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    perl $0 genome.fasta input1.gff3 [input2.gff3 ...] > transcripts.fasta

    本程序用于根据GFF3中RNA信息的exon Feature，使用基因组序列转换出所有转录本序列。

USAGE
if (@ARGV==0){die $usage}

my ($help_flag);
GetOptions(
    "help" => \$help_flag,
);

my $genome_file = shift @ARGV;
open IN, $genome_file or die "Can not open file $genome_file, $!";
my (%seq, $seq_id);
while (<IN>) {
    chomp;
    if (m/^>(\S+)/) { $seq_id = $1; }
    else { $seq{$seq_id} .= $_; }
}
close IN;

my (%exon, %locus, %strand);
foreach ( @ARGV ) {
    open IN, $_ or die "Can not open file $_, $!";
    while (<IN>) {
        if (/\texon\t/) {
            @_ = split /\t/;
            if ( $_[8] =~ m/Parent=([^;\s]+)/ ) {
                $exon{$1}{"$_[3]\t$_[4]\t$_[7]"} = $_[3];
                $locus{$1} = $_[0];
                $strand{$1} = $_[6];
            }
        }
    }
    close IN;
}

foreach my $id (sort keys %exon) {
    my @exon = sort {$exon{$id}{$a} <=> $exon{$id}{$b}} keys %{$exon{$id}};
    my $exon_seq;
    foreach (@exon) {
        @_ = split /\t/, $_;
        $exon_seq .= substr($seq{$locus{$id}}, $_[0] - 1, $_[1] - $_[0] + 1);
    }
    if ( $strand{$id} eq "-" ) {
        $exon_seq = reverse $exon_seq;
        $exon_seq =~ tr/ATCGatcgn/TAGCTAGCN/;
    }
    print ">$id\n$exon_seq\n";
}
