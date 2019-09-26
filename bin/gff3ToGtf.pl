#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 genome.fasta file.gff3 > file.gtf

USAGE
if (@ARGV == 0) { die $usage }

my (%fasta, $fasta_id);
open FASTA, '<', $ARGV[0] or die $!;
while (<FASTA>) {
    chomp;
    if (/>(\S+)/) { $fasta_id = $1; }
    else { $_ = uc($_); $fasta{$fasta_id} .= $_; }
}
close FASTA;

my (%gene, %gff);
open GFF, '<', $ARGV[1] or die $!;
while ( <GFF> ) {
    if (/\tmRNA\t.*ID=([^;\n]+)/) {
        my $mRNA_id = $1;
        my $gene_id;
        if (/Parent=([^;\n]+)/) {
            $gene_id = $1;
            $gene{$gene_id}{$mRNA_id} = 1;
        }
    }
    if (/Parent=([^;\n]+)/) {
        $gff{$1} .= $_;
    }
}
close GFF;

foreach my $gene_id (sort keys %gene) {
    foreach my $mRNA_id (sort keys %{$gene{$gene_id}}) {
        my @output;

        my $gff_content = $gff{$mRNA_id};
        my @gff_content = split /\n/, $gff_content;

        my @utr = grep /UTR\t/, @gff_content;
        map {s/five_prime_UTR\t/5UTR\t/; s/three_prime_UTR\t/3UTR\t/} @utr;
        push @output, @utr;

        my @exon = grep /\texon\t/, @gff_content;
        push @output, @exon;
        
        my @cds = grep /\tCDS\t/, @gff_content;
        push @output, @cds;

        my %sort_cds;
        foreach (@cds) { @_ = split; $sort_cds{$_} = $_[3] }
        @cds = sort { $sort_cds{$a} <=> $sort_cds{$b} } @cds;

        my ($start_codon, $stop_codon);
        my @first = split /\t/, $cds[0];
        my @second = split /\t/, $cds[1];
        my @three = split /\t/, $cds[-2];
        my @four = split /\t/, $cds[-1];
        my $strand = $_[6];
        #print STDERR "$strand\t$mRNA_id\n";
        if ($strand eq "+") {
            # start codon
            my $start_codon_bases;
            my $length = $first[4] - $first[3] + 1;
            if ($length >= 3) {
                my $start = $first[3];
                my $end = $start + 2;
                $start_codon = "$first[0]\t$first[1]\tstart_codon\t$start\t$end\t\.\t$strand\t0\tgene_id \"$gene_id\"; transcript_id \"$mRNA_id\";";
                $start_codon_bases = substr($fasta{$first[0]}, $start - 1, 3);
            }
            else {
                print STDERR "$mRNA_id\t+\tstart_codon\tspliced\n";
                $start_codon = "$first[0]\t$first[1]\tstart_codon\t$first[3]\t$first[4]\t\.\t$strand\t0\tgene_id \"$gene_id\"; transcript_id \"$mRNA_id\";";
                my $end = $second[3] + (3 - $length - 1);
                $start_codon .= "\n$first[0]\t$first[1]\tstart_codon\t$second[3]\t$end\t\.\t$strand\t$length\tgene_id \"$gene_id\"; transcript_id \"$mRNA_id\";";
                $start_codon_bases = substr($fasta{$first[0]}, $first[3] - 1, $length);
                $start_codon_bases .= substr($fasta{$first[0]}, $second[3] - 1, 3 - $length);
            }
            if ($start_codon_bases eq "ATG") { push @output, $start_codon }
            else { print STDERR "$mRNA_id\t+\tstart_codon\t$start_codon_bases\n"; }

            # stop codon
            my $stop_codon_bases;
            $length = $four[4] - $four[3] + 1;
            if ($length >= 3) {
                my $end = $four[4];
                my $start = $end - 2;
                $stop_codon = "$four[0]\t$four[1]\tstop_codon\t$start\t$end\t\.\t$strand\t0\tgene_id \"$gene_id\"; transcript_id \"$mRNA_id\";";
                $stop_codon_bases = substr($fasta{$four[0]}, $start - 1, 3);
            }
            else {
                print STDERR "$mRNA_id\t+\tstop_codon\tspliced\n";
                my $frame = 3 - $length;
                $stop_codon = "$four[0]\t$four[1]\tstop_codon\t$four[3]\t$four[4]\t\.\t$strand\t$frame\tgene_id \"$gene_id\"; transcript_id \"$mRNA_id\";";
                my $start = $three[4] - (3 - $length - 1);
                $stop_codon .= "\n$four[0]\t$four[1]\tstop_codon\t$start\t$three[4]\t\.\t$strand\t0\tgene_id \"$gene_id\"; transcript_id \"$mRNA_id\";";
                $stop_codon_bases = substr($fasta{$four[0]}, $start - 1, 3 - $length);
                $stop_codon_bases .= substr($fasta{$four[0]}, $four[3] - 1, $length);
            }
            if ($stop_codon_bases =~ /^(TAA|TAG|TGA)$/) { push @output, $stop_codon; }
            else { print STDERR "$mRNA_id\t+\tstop_codon\t$stop_codon_bases\n"; }
            #print "$start_codon_bases\t$stop_codon_bases\n";
        }
        elsif ($strand eq "-") {
            # stop codon
            my $stop_codon_bases;
            my $length = $first[4] - $first[3] + 1;
            if ($length >= 3) {
                my $start = $first[3];
                my $end = $start + 2;
                $stop_codon = "\n$first[0]\t$first[1]\tstop_codon\t$start\t$end\t\.\t$strand\t0\tgene_id \"$gene_id\"; transcript_id \"$mRNA_id\";";
                $stop_codon_bases = substr($fasta{$first[0]}, $first[3] - 1, 3);
            }
            else {
                print STDERR "$mRNA_id\t-\tstop_codon\tspliced\n";
                my $frame = 3 - $length;
                $stop_codon = "$first[0]\t$first[1]\tstop_codon\t$first[3]\t$first[4]\t\.\t$strand\t$frame\tgene_id \"$gene_id\"; transcript_id \"$mRNA_id\";";
                my $end = $second[3] + (3 - $length - 1);
                $stop_codon .= "\n$first[0]\t$first[1]\tstop_codon\t$second[3]\t$end\t.\t$strand\t0\tgene_id \"$gene_id\"; transcript_id \"$mRNA_id\";";
                $stop_codon_bases = substr($fasta{$first[0]}, $first[3] - 1, $length);
                $stop_codon_bases .= substr($fasta{$first[0]}, $second[3] - 1, 3 - $length);
            }
            $stop_codon_bases = reverse $stop_codon_bases;
            $stop_codon_bases =~ tr/ATCGatcg/TAGCtagc/;
            if ($stop_codon_bases =~ /^(TAA|TAG|TGA)$/) { push @output, $stop_codon; }
            else { print STDERR "$mRNA_id\t-\tstop_codon\t$stop_codon_bases\n"; }

            # start codon
            my $start_codon_bases;
            $length = $four[4] - $four[3] + 1;
            if ($length >= 3) {
                my $end = $four[4];
                my $start = $end - 2;
                $start_codon = "$four[0]\t$four[1]\tstart_codon\t$start\t$end\t\.\t$strand\t0\tgene_id \"$gene_id\"; transcript_id \"$mRNA_id\";";
                $start_codon_bases = substr($fasta{$four[0]}, $start - 1, 3);
            }
            else {
                print STDERR "$mRNA_id\t-\tstart_codon\tspliced\n";
                $start_codon = "$four[0]\t$four[1]\tstart_codon\t$four[3]\t$four[4]\t\.\t$strand\t0\tgene_id \"$gene_id\"; transcript_id \"$mRNA_id\";";
                my $start = $three[4] - (3 - $length - 1);
                $start_codon .="\n$four[0]\t$four[1]\tstart_codon\t$start\t$three[4]\t\.\t$strand\t$length\tgene_id \"$gene_id\"; transcript_id \"$mRNA_id\";";
                $start_codon_bases = substr($fasta{$four[0]}, $start - 1, 3 - $length);
                $start_codon_bases .= substr($fasta{$four[0]}, $four[3] - 1, $length);
            }
            $start_codon_bases = reverse $start_codon_bases;
            $start_codon_bases =~ tr/ATCGatcg/TAGCtagc/;
            if ($start_codon_bases eq "ATG") { push @output, $start_codon }
            else { print STDERR "$mRNA_id\t-\tstart_codon\t$start_codon_bases\n"; }
            #print "$start_codon_bases\t$stop_codon_bases\n";
        }

        @output = map {split /\n/} @output;
        my (%output_sort, %output_sort1);
        foreach (@output) {
            @_ = split;
            $output_sort{$_} = $_[3];
            if (/\t5UTR\t/) { $output_sort1{$_} = 1; }
            elsif (/\tstart_codon\t/) { $output_sort1{$_} = 2; }
            elsif (/\t(CDS|exon)\t/) { $output_sort1{$_} = 3; }
            elsif (/\tstop_codon\t/) { $output_sort1{$_} = 4; }
            elsif (/\t3UTR\t/) { $output_sort1{$_} = 5; }
        }
        if ($strand eq "+") {
            @output = sort { $output_sort1{$a} <=> $output_sort1{$b} or $output_sort{$a} <=> $output_sort{$b} } @output;
        }
        else {
            @output = sort { $output_sort1{$a} <=> $output_sort1{$b} or $output_sort{$b} <=> $output_sort{$a} } @output;
        }
        map {s/(.*)\t.*/$1\tgene_id \"$gene_id\"; transcript_id \"$mRNA_id\"; gene_name \"$gene_id\";\n/} @output;
        foreach (@output) {
            print;
        }
    }
    print "\n";
}
