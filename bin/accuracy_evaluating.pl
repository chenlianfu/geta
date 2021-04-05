#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 strandard.gff3 query.gff3

USAGE
if (@ARGV==0){die $usage}

# 读取标准的GFF3文件，该文件可以包含可变剪接
open IN, $ARGV[0] or die $!;
my (%feature, %gene);
while (<IN>) {
    next if m/^#/;
    next if m/^\s*$/;
    @_ = split /\t/;
    if ($_[8] =~ m/Parent=([^;,\s]+)/) {
        $feature{$1} .= $_;
    }
    else {
        if ($_[8] =~ m/ID=([^;\s]+)/) {
            my $id = $1;
            $gene{$id} = $_ if m/\tgene\t/;
        }
    }
}
close IN;

# 读取query GFF3文件，该文件可以包含可变剪接
open IN, $ARGV[1] or die $!;
my (%feature_query, %gene_query);
while (<IN>) {
    next if m/^#/;
    next if m/^\s*$/;
    @_ = split /\t/;
    if ($_[8] =~ m/Parent=([^;,\s]+)/) {
        $feature_query{$1} .= $_;
    }
    else {
        if ($_[8] =~ m/ID=([^;\s]+)/) {
            my $id = $1;
            $gene_query{$id} = $_ if m/\tgene\t/;
        }
    }
}
close IN;

# 对标准GFF3文件进行解析，得到gene的CDS信息，并进行INDEX数据库
my (%gene_info, %CDS_info, %index, $total_gene_number, $total_CDS_number, $total_bp_number);
foreach my $gene_id (keys %gene) {
    my $gene_info = $feature{$gene_id};
    $total_gene_number ++ if $gene_info =~ m/\tmRNA\t/;
    #print "OK:\t$gene_id\n";
    chomp($gene_info);
    my @gene_info = split /\n/, $gene_info;
    foreach (@gene_info) {
        my @CDS;
        if (m/\tmRNA\t/) {
            if (m/ID=([^;\s]+)/) {
                my $mRNA_info = $feature{$1};
                #print $mRNA_info;
                chomp($mRNA_info);
                my @mRNA_info = split /\n/, $mRNA_info;
                foreach (@mRNA_info) {
                    if (m/\tCDS\t/) {
                        @_ = split /\t/;
                        push @CDS, "$_[3]\t$_[4]";
                        $CDS_info{$gene_id}{"$_[3]\t$_[4]"} = 1;
                    }
                }
            }
        }
        @CDS = sort {$a <=> $b} @CDS;
        my $CDS = join "\n", @CDS;
        #print "$_\n$CDS\n";
        $gene_info{$gene_id}{$CDS} = 1;
    }

    my $gene_desc = $gene{$gene_id};
    @_ = split /\t/, $gene_desc;
    my $int1 = int($_[3] / 1000);
    my $int2 = int($_[4] / 1000);
    foreach ($int1 .. $int2) {
        $index{$_[0]}{$_[6]}{$_}{$gene_id} = 1;
    }

    my @CDSs = sort {$a <=> $b} keys %{$CDS_info{$gene_id}};
    $total_CDS_number += @CDSs;
    my ($start, $end) = split /\t/, $CDSs[0];
    $total_bp_number += $end - $start + 1;
    shift @CDSs;
    foreach (@CDSs) {
        @_ = split /\t/, $_;
        if ($_[0] > $end) {
            $total_bp_number += $_[1] - $_[0] + 1;
            ($start, $end) = @_;
        }
        else {
            if ($_[1] > $end) {
                $total_bp_number += $_[1] - $end;
                $end = $_[1];
            }
        }
    }
}

#print "$total_gene_number\t$total_CDS_number\t$total_bp_number\n";

# 进行positive分析
my ($query_gene_number, $positive_gene_number, $query_CDS_number, $positive_CDS_number, $query_bp_number, $positive_bp_number, $positive_bp_number_ref, %positive_gene);
foreach my $gene_id (keys %gene_query) {
    my $gene_desc = $gene_query{$gene_id};
    @_ = split /\t/, $gene_desc;
    my $strand = $_[6];
    my $int1 = int($_[3] / 1000);
    my $int2 = int($_[4] / 1000);
    #print "OK2:\t$gene_id\t$int1\t$int2\n";
    my (%target_gene_id, %target_CDS);
    foreach my $int ($int1 .. $int2) {
        my @id = keys %{$index{$_[0]}{$_[6]}{$int}};
        foreach (@id) {
            $target_gene_id{$_} = 1;
            foreach (keys %{$CDS_info{$_}}) {
                $target_CDS{$_} = 1;
            }
        }
    }
    #foreach (keys %target_gene_id) { print "gene_DB:\t$_\n"; }
    #foreach (keys %target_CDS) { print "CDS_DB:\t$_\n"; }

    my $gene_info = $feature_query{$gene_id};
    if ($gene_info =~ m/\tmRNA\t/) {
        $query_gene_number ++;
    }
    else {
        next;
    }

    chomp($gene_info);
    my @gene_info = split /\n/, $gene_info;
    my (%CDS, %transcript_CDS);
    foreach (@gene_info) {
        if (m/ID=([^;\s]+)/) {
            my $mRNA_info = $feature_query{$1};
            chomp($mRNA_info);
            my @mRNA_info = split /\n/, $mRNA_info;
            my @CDS;
            foreach (@mRNA_info) {
                if (m/\tCDS\t/) {
                     @_ = split /\t/;
                     push @CDS, "$_[3]\t$_[4]";
                     $CDS{"$_[3]\t$_[4]"} = 1;
                 }
             }
             @CDS = sort {$a <=> $b} @CDS;
             my $CDS = join "\n", @CDS;
             $transcript_CDS{$CDS} = 1;
        }
    }

    my @CDS = sort {$a <=> $b} keys %CDS;
    $query_CDS_number += @CDS;
    my ($start, $end) = split /\t/, $CDS[0];
    $query_bp_number += $end - $start + 1;
    shift @CDS;
    foreach (@CDS) {
        @_ = split /\t/, $_;
        if ($_[0] > $end) {
            $query_bp_number += $_[1] - $_[0] + 1;
            ($start, $end) = @_;
        }
        else {
            if ($_[1] > $end) {
                $query_bp_number += $_[1] - $end;
                $end = $_[1];
            }
        }
    }

    my (%positive_CDS, %positive_CDS_ref);
    my $same = 0;
    foreach (keys %transcript_CDS) {
        my $transcript_CDS =  $_;
        if ($strand eq "+") {
            $_ =~ s/\d+//;
        }
        elsif ($strand eq "-") {
            $_ =~ s/\d+$//;
        }
        my $transcript_CDS_for_validation = $_;

        CDS_COMPARE: foreach my $target_gene_id (keys %target_gene_id) {
            my @CDS_target = keys %{$gene_info{$target_gene_id}};
            foreach my $CDS_target (@CDS_target) {
                my $CDS_target_orig = $CDS_target;
                if ($strand eq "+") {
                    $CDS_target =~ s/\d+//;
                }
                elsif ($strand eq "-") {
                    $CDS_target =~ s/\d+$//;
                }

                if ($transcript_CDS_for_validation eq $CDS_target) {
                    $same = 1;
                    my @transcript_CDS = split /\n/, $transcript_CDS;
                    foreach (@transcript_CDS) {
                        $positive_CDS{$_} = 1;
                    }
                    foreach (split /\n/, $CDS_target_orig) {
                        $positive_CDS_ref{$_} = 1;
                    }

                    print STDERR "$gene_id\t$target_gene_id\n";
                    last CDS_COMPARE;
                }
            }
        }
    }
    $positive_gene_number ++ if $same == 1;

    my %bp_cds;
    foreach my $cds (keys %CDS) {
        foreach (keys %target_CDS) {
            @_ = split /\t/;
            my @cds = split /\t/, $cds;
            if ($_[0] <= $cds[1] && $cds[0] <= $_[1]) {
                my @cds_region;
                push @cds_region, @_;
                push @cds_region, @cds;
                @cds_region = sort {$a <=> $b} @cds_region;
                $bp_cds{"$cds_region[1]\t$cds_region[2]"} = 1;
            }

            if ($_ eq $cds) {
                $positive_CDS{$cds} = 1;
                last;
            }
        }
    }
    $positive_CDS_number += keys %positive_CDS;

    my (%bp_cds_query, %bp_cds_ref);
    foreach (keys %bp_cds) { $bp_cds_query{$_} = 1; $bp_cds_ref{$_} = 1; }
    foreach (keys %positive_CDS) { $bp_cds_query{$_} = 1; }
    foreach (keys %positive_CDS_ref) { $bp_cds_ref{$_} = 1; }

    my @bp_cds = sort {$a <=> $b} keys %bp_cds_query;
    my ($start, $end) = split /\t/, $bp_cds[0];
    $positive_bp_number += $end - $start + 1;
    shift @bp_cds;
    foreach (@bp_cds) {
        @_ = split /\t/, $_;
        if ($_[0] > $end) {
            $positive_bp_number += $_[1] - $_[0] + 1;
            ($start, $end) = @_;
        }
        else { 
            if ($_[1] > $end) {
                $positive_bp_number += $_[1] - $end;
                $end = $_[1];
            }
        }
    }

    my @bp_cds = sort {$a <=> $b} keys %bp_cds_ref;
    my ($start, $end) = split /\t/, $bp_cds[0];
    $positive_bp_number_ref += $end - $start + 1;
    shift @bp_cds;
    foreach (@bp_cds) {
        @_ = split /\t/, $_;
        if ($_[0] > $end) {
            $positive_bp_number_ref += $_[1] - $_[0] + 1;
            ($start, $end) = @_;
        }
        else { 
            if ($_[1] > $end) {
                $positive_bp_number_ref += $_[1] - $end;
                $end = $_[1];
            }
        }
    }
}

#print "$total_gene_number\t$total_CDS_number\t$total_bp_number\n";
#print "$query_gene_number\t$positive_gene_number\t$query_CDS_number\t$positive_CDS_number\t$query_bp_number\t$positive_bp_number\n";
my $sensitivity_gene = $positive_gene_number  * 100 / $total_gene_number;
my $sensitivity_CDS = $positive_CDS_number * 100 / $total_CDS_number;
my $sensitivity_bp = $positive_bp_number_ref * 100 / $total_bp_number;
my $specificity_gene = $positive_gene_number * 100 / $query_gene_number;
my $specificity_CDS = $positive_CDS_number * 100 / $query_CDS_number;
my $specificity_bp = $positive_bp_number * 100 / $query_bp_number;
=cut
print "\tTotal_number\tPredicted_number\tTrue_positive\tSensivitity\tSpecificity\n";
printf "Gene_level\t$total_gene_number\t$query_gene_number\t$positive_gene_number\t%.2f\%\t%.2f\%\t\n", $sensitivity_gene, $specificity_gene;
printf "CDS_level\t$total_CDS_number\t$query_CDS_number\t$positive_CDS_number\t%.2f\%\t%.2f\%\t\n", $sensitivity_CDS, $specificity_CDS;
printf "Nucleotide\t$total_bp_number\t$query_bp_number\t$positive_bp_number\t%.2f\%\t%.2f\%\t\n", $sensitivity_bp, $specificity_bp;
=cut
print "gene_prediction_method\tgene_standard\tgene_predicted\tgene_TP\tgene_sensitivity\tgene_specificity\tCDS_standard\tCDS_predicted\tCDS_TP\tCDS_sensitivity\tCDS_specificity\tbp_standard\tbp_predicted\tbp_TP\tbp_sensitivity\tbp_specificity\n";
printf "@ARGV[1]\t$total_gene_number\t$query_gene_number\t$positive_gene_number\t%.2f\%\t%.2f\%\t$total_CDS_number\t$query_CDS_number\t$positive_CDS_number\t%.2f\%\t%.2f\%\t$total_bp_number\t$query_bp_number\t$positive_bp_number\t%.2f\%\t%.2f\%\n", $sensitivity_gene, $specificity_gene, $sensitivity_CDS, $specificity_CDS, $sensitivity_bp, $specificity_bp;

foreach (sort keys %positive_gene) {
    print STDERR "$_\n";
}
