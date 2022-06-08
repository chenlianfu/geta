#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    perl $0 [options] old.gff3 > new.gff3

    --prefix <string>    default: g

USAGE
if (@ARGV == 0) { die $usage }

my ($genePrefix);
GetOptions(
    "prefix:s" => \$genePrefix,
);

$genePrefix ||= 'g';

open IN, '<', $ARGV[0] or die $!;
my (%gff3_gene_info, %gff3_mRNA_info, %gff3_featrue_info, %gff3_parent1, %gff3_parent2, %lines);
while ( <IN> ) {
    next if /^#/;
    next if /^\s/;
    next if exists $lines{$_};
    #print;
    $lines{$_} = 1;
    @_ = split /\t/;
    if ($_[2] eq "gene") {
        $_[8] =~ /ID=([^;\s]+)/;
        $gff3_gene_info{$1} .= $_;
    }
    elsif ($_[2] eq "mRNA") {
        my $mRNA_id = $1 if $_[8] =~ /ID=([^;\s]+)/;
        my $parent_id = $1 if $_[8] =~ /Parent=([^;\s]+)/;
        $gff3_mRNA_info{$mRNA_id} = $_;
        $gff3_parent1{$parent_id}{$mRNA_id} = $_[3];
        $gff3_parent2{$parent_id}{$mRNA_id} = $_[4];
    }
    else {
        my $parent_id = $1 if $_[8] =~ /Parent=([^;\s]+)/;
        $gff3_featrue_info{$parent_id} .= $_;
    }
}
close IN;

my @geneSym = keys %gff3_gene_info;
my (%geneSymSort1, %geneSymSort2, %geneSymSort3, %geneSymSort4, %geneSymSort5);
foreach ( @geneSym ) {
    @_ = split /\t/, $gff3_gene_info{$_};
    $geneSymSort1{$_} = $_[0];
    $geneSymSort2{$_} = $_[3];
    $geneSymSort3{$_} = $_[4];
    $geneSymSort4{$_} = $_[6];
    $geneSymSort5{$_} = $_[5];
}
@geneSym = sort { $geneSymSort1{$a} cmp $geneSymSort1{$b} or $geneSymSort2{$a} <=> $geneSymSort2{$b} or $geneSymSort3{$a} <=> $geneSymSort3{$b} or $geneSymSort4{$a} cmp $geneSymSort4{$b} or $geneSymSort5{$b} <=> $geneSymSort5{$a} } @geneSym;
my @geneSym_new;
foreach (@geneSym) {
    my @mRNA_id = keys %{$gff3_parent1{$_}};
    my $mRNA_info;
    foreach (@mRNA_id) { $mRNA_info .= $gff3_featrue_info{$_} };
    next unless $mRNA_info =~ m/\tCDS\t/;
    push @geneSym_new, $_ if @mRNA_id >= 1;
}

my $geneNum;
foreach my $gene_id (@geneSym_new) {
    #print "$gene_id\n";
    # 对gene信息进行分析
    $geneNum ++;
    my $gene_name = $genePrefix . '0' x (length(@geneSym) - length($geneNum)) . $geneNum;

    my $gene_info = $gff3_gene_info{$gene_id};
    @_ = split /\t/, $gene_info;
    my $strand = $_[6];
    my $last_filed = pop @_;
    $_ = join "\t", @_;
    print "$_\t";

    my @out = ("ID=$gene_name", "Name=$gene_name");
    $last_filed =~ s/[;\s]*$//;
    @_ = split /[;\s]/, $last_filed;
    foreach (@_) {
        next if m/ID=/; next if m/Name=/; next if m/Parent=/; push @out, $_;
    }
    my $out = join ";", @out;
    print "$out\n";

    # 对隶属于该gene的mRNAs信息进行分析
    my @mRNA_id = sort {$gff3_parent1{$gene_id}{$a} <=> $gff3_parent1{$gene_id}{$b} or $gff3_parent2{$gene_id}{$a} <=> $gff3_parent2{$gene_id}{$b}} keys %{$gff3_parent1{$gene_id}};
    my $transcriptIdNum;
    foreach my $mRNA_id (@mRNA_id) {
        my $mRNA_info = $gff3_featrue_info{$mRNA_id};
        next unless $mRNA_info =~ m/\tCDS\t/;
        $transcriptIdNum ++;
        my $mRNA_info = $gff3_mRNA_info{$mRNA_id};
        @_ = split /\t/, $mRNA_info;
        my $last_filed = pop @_;
        $_ = join "\t", @_;
        print "$_\t";

        my @out = ("ID=$gene_name.t$transcriptIdNum", "Parent=$gene_name");
        $last_filed =~ s/[;\s]*$//;
        @_ = split /[;\s]/, $last_filed;
        foreach (@_) {
            next if m/ID=/; next if m/Name=/; next if m/Parent=/; push @out, $_;
        }
        my $out = join ";", @out;
        print "$out\n";

        my @feature_info = split /\n/, $gff3_featrue_info{$mRNA_id};
        my (%featrue_sort1, %featrue_sort2, %featrue_sort3);
        foreach (@feature_info) {
            @_ = split /\t/;
            $featrue_sort1{$_} = $_[3];
            $featrue_sort2{$_} = $_[4];
            $featrue_sort3{$_} = $_[2];
        }
        @feature_info = sort { $featrue_sort1{$a} <=> $featrue_sort1{$b} or $featrue_sort2{$a} <=> $featrue_sort2{$b} or $featrue_sort3{$b} cmp $featrue_sort3{$a} } @feature_info if $strand eq "+";
        @feature_info = sort { $featrue_sort2{$b} <=> $featrue_sort2{$a} or $featrue_sort1{$b} <=> $featrue_sort1{$a} or $featrue_sort3{$b} cmp $featrue_sort3{$a} } @feature_info if $strand eq "-";


        my ($fiveUTR_Num, $threeUTR_Num, $exonNum, $cdsNum);
        foreach (@feature_info) {
            @_ = split /\t/;

            if ($_[2] eq "five_prime_UTR") {
                $fiveUTR_Num ++;
                my $last_filed = pop @_;
                $_ = join "\t", @_;
                print "$_\t";
                my @out = ("ID=$gene_name.t$transcriptIdNum.utr5p$fiveUTR_Num", "Parent=$gene_name.t$transcriptIdNum");
                $last_filed =~ s/[;\s]*$//;
                @_ = split /[;\s]/, $last_filed;
                foreach (@_) {
                    next if m/ID=/; next if m/Name=/; next if m/Parent=/; push @out, $_;
                }
                my $out = join ";", @out;
                print "$out\n";
            }
            elsif ($_[2] eq "three_prime_UTR") {
                $threeUTR_Num ++;
                my $last_filed = pop @_;
                $_ = join "\t", @_;
                print "$_\t";
                my @out = ("ID=$gene_name.t$transcriptIdNum.utr3p$threeUTR_Num", "Parent=$gene_name.t$transcriptIdNum");
                $last_filed =~ s/[;\s]*$//;
                @_ = split /[;\s]/, $last_filed;
                foreach (@_) {
                    next if m/ID=/; next if m/Name=/; next if m/Parent=/; push @out, $_;
                }
                my $out = join ";", @out;
                print "$out\n";
            }
            elsif ($_[2] eq "exon") {
                $exonNum ++;
                my $last_filed = pop @_;
                $_ = join "\t", @_;
                print "$_\t";
                my @out = ("ID=$gene_name.t$transcriptIdNum.exon$exonNum", "Parent=$gene_name.t$transcriptIdNum"); 
                $last_filed =~ s/[;\s]*$//;
                @_ = split /[;\s]/, $last_filed;
                foreach (@_) {
                    next if m/ID=/; next if m/Name=/; next if m/Parent=/; push @out, $_;
                }
                my $out = join ";", @out;
                print "$out\n";
            }
            elsif ($_[2] eq "CDS") {
                $cdsNum ++;
                my $last_filed = pop @_;
                $_ = join "\t", @_;
                print "$_\t";
                my @out = ("ID=$gene_name.t$transcriptIdNum.CDS$cdsNum", "Parent=$gene_name.t$transcriptIdNum");
                $last_filed =~ s/[;\s]*$//;
                @_ = split /[;\s]/, $last_filed;
                foreach (@_) {
                    next if m/ID=/; next if m/Name=/; next if m/Parent=/; push @out, $_;
                }
                my $out = join ";", @out;
                print "$out\n";
            }
        }
    }
    print "\n";
}
