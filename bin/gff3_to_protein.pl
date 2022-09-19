#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    perl $0 genome.fasta input1.gff3 [input2.gff3 ...] > proteins.fasta

    本程序用于根据GFF3中mRNA信息的CDS Feature，使用基因组序列转换出所有转录本的Protein序列。

    --out_CDS    default: None
    添加该参数后，程序输出CDS序列，而不是默认的Protein序列。

USAGE
if (@ARGV==0){die $usage}

my ($help_flag, $out_CDS);
GetOptions(
    "help" => \$help_flag,
    "out_CDS!" => \$out_CDS,
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

my (%cds, %locus, %strand);
foreach ( @ARGV ) {
    open IN, $_ or die "Can not open file $_, $!";
    while (<IN>) {
        if (/\tCDS\t/) {
            @_ = split /\t/;
            $_[8] =~ m/Parent=([^;\s]+)/;
            $cds{$1}{"$_[3]\t$_[4]\t$_[7]"} = $_[3];
            $locus{$1} = $_[0];
            $strand{$1} = $_[6];
        }
    }
    close IN;
}

foreach my $id (sort keys %cds) {
    my @cds = sort {$cds{$id}{$a} <=> $cds{$id}{$b}} keys %{$cds{$id}};
    my $cds_seq;
    foreach (@cds) {
        @_ = split /\t/, $_;
        $cds_seq .= substr($seq{$locus{$id}}, $_[0] - 1, $_[1] - $_[0] + 1);
    }
    if ($strand{$id} eq "+") {
        my $frame = 0;
        $frame = $1 if $cds[0] =~ m/(\d+)$/;
        print STDERR "Warning: $id\tthe frame not equal 0\n" if $frame != 0;
        my $pep = &cds2pep($cds_seq, $frame);
        if ( $out_CDS ) {
            print ">$id\n$cds_seq\n";
        }
        else {
            print ">$id\n$pep\n";
        }
    }
    elsif ($strand{$id} eq "-") {
        $cds_seq = &rc($cds_seq);
        my $frame = 0;
        $frame = $1 if $cds[-1] =~ m/(\d+)$/;
        print STDERR "Warning: $id\tthe frame not equal 0\n" if $frame != 0;
        my $pep = &cds2pep($cds_seq, $frame);
        if ( $out_CDS ) {
            print ">$id\n$cds_seq\n";
        }
        else {
            print ">$id\n$pep\n";
        }
    }
}

sub rc {
    my $seq = shift @_;
    $seq = reverse $seq;
    $seq =~ tr/ATCGatcgn/TAGCTAGCN/;
    return $seq;
}

sub cds2pep {
    my %cds2pep = (
        "TTT" => "F",
        "TTC" => "F",
        "TTA" => "L",
        "TTG" => "L",
        "TCT" => "S",
        "TCC" => "S",
        "TCA" => "S",
        "TCG" => "S",
        "TAT" => "Y",
        "TAC" => "Y",
        "TAA" => "*",
        "TAG" => "*",
        "TGT" => "C",
        "TGC" => "C",
        "TGA" => "*",
        "TGG" => "W",
        "CTT" => "L",
        "CTC" => "L",
        "CTA" => "L",
        "CTG" => "L",
        "CCT" => "P",
        "CCC" => "P",
        "CCA" => "P",
        "CCG" => "P",
        "CAT" => "H",
        "CAC" => "H",
        "CAA" => "Q",
        "CAG" => "Q",
        "CGT" => "R",
        "CGC" => "R",
        "CGA" => "R",
        "CGG" => "R",
        "ATT" => "I",
        "ATC" => "I",
        "ATA" => "I",
        "ATG" => "M",
        "ACT" => "T",
        "ACC" => "T",
        "ACA" => "T",
        "ACG" => "T",
        "AAT" => "N",
        "AAC" => "N",
        "AAA" => "K",
        "AAG" => "K",
        "AGT" => "S",
        "AGC" => "S",
        "AGA" => "R",
        "AGG" => "R",
        "GTT" => "V",
        "GTC" => "V",
        "GTA" => "V",
        "GTG" => "V",
        "GCT" => "A",
        "GCC" => "A",
        "GCA" => "A",
        "GCG" => "A",
        "GAT" => "D",
        "GAC" => "D",
        "GAA" => "E",
        "GAG" => "E",
        "GGT" => "G",
        "GGC" => "G",
        "GGA" => "G",
        "GGG" => "G",
    );
    my $seq = shift @_;
    my $fram = shift @_;
    my $gene = shift @_;
    $seq =~ s/\w{$fram}//;
    my $pep;
    while ((length $seq) >= 3) {
        $seq =~ s/(\w{3})//;
        if (exists $cds2pep{$1}) {
            $pep .= $cds2pep{$1};
        }
        else {
            $pep .= 'X';
        }
    }
    return $pep;
}
