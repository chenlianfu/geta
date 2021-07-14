#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 genome.fasta genome.gff3 > proteins.fasta

USAGE
if (@ARGV==0){die $usage}

open IN, $ARGV[0] or die $!;
my (%seq, $seq_id);
while (<IN>) {
    chomp;
    if (m/^>(\S+)/) { $seq_id = $1; }
    else { $seq{$seq_id} .= $_; }
}
close IN;

open IN, $ARGV[1] or die $!;
my (%cds, %locus, %strand);
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
        print ">$id\n$pep\n";
    }
    elsif ($strand{$id} eq "-") {
        $cds_seq = &rc($cds_seq);
        my $frame = 0;
        $frame = $1 if $cds[-1] =~ m/(\d+)$/;
        print STDERR "Warning: $id\tthe frame not equal 0\n" if $frame != 0;
        my $pep = &cds2pep($cds_seq, $frame);
        print ">$id\n$pep\n";
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
