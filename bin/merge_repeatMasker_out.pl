#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    perl $0 repeatMasker/genome.fasta.out repeatModeler/genome.fasta.out > genome.repeat.stats

    --min_coverge_ratio <float>    default: 0.25
    设置最小覆盖率阈值。若匹配区域对目标重复序列的覆盖率低于此阈值，则过滤掉此结果。

USAGE
if (@ARGV==0){die $usage}
my $min_coverge_ratio;
GetOptions(
    "min_coverge_ratio:f" => \$min_coverge_ratio,
);
$min_coverge_ratio ||= 0.25;

my (%locus, %scaffold, @scaffold);
open OUT, '>', "genome.repeat.out" or die $!;
open GFF3, '>', "genome.repeat.gff3" or die $!;
print GFF3 "##gff-version 3\n";

foreach (@ARGV) {
    open IN, $_ or die $!;
    $_ = <IN>; print OUT;
    $_ = <IN>; print OUT;
    $_ = <IN>; print OUT;
    while (<IN>) {
        s/^/ /;
        my @aaa = split /\s+/;

        my $ratio;
        if ($aaa[12] =~ m/\((\d+)\)/) {
            $ratio = ($aaa[13] - $aaa[14] + 1) / ($1 + $aaa[13]);
        }
        elsif ($aaa[14] =~ m/\((\d+)\)/) {
            $ratio = ($aaa[13] - $aaa[12] + 1) / ($aaa[13] + $1);
        }
        next if $ratio < $min_coverge_ratio;

        $locus{$aaa[5]}{"$aaa[6]\t$aaa[7]"}{"description"} = $_;
        $locus{$aaa[5]}{"$aaa[6]\t$aaa[7]"}{"start"} = $aaa[6];
        $locus{$aaa[5]}{"$aaa[6]\t$aaa[7]"}{"strand"} = $aaa[9];
        $locus{$aaa[5]}{"$aaa[6]\t$aaa[7]"}{"name"} = $aaa[10];
        $locus{$aaa[5]}{"$aaa[6]\t$aaa[7]"}{"repeatClass"}  = $aaa[11];
        $locus{$aaa[5]}{"$aaa[6]\t$aaa[7]"}{"target"} = "$aaa[12]\t$aaa[13]\t$aaa[14]";

        push @scaffold, $aaa[5] unless exists $scaffold{$aaa[5]};
        $scaffold{$aaa[5]} = 1;
    }
}
close IN;

my (%stats1, %stats2, %repeatRegion);
foreach my $scaffold (@scaffold) {
    my %masked_site;
    my @region;
    my @locus = sort { $locus{$scaffold}{$a}{"start"} <=> $locus{$scaffold}{$b}{"start"} } keys %{$locus{$scaffold}};
    foreach my $locus (@locus) {
        print OUT $locus{$scaffold}{$locus}{"description"};

        my @start_end = split /\t/, $locus;
        my $name = $locus{$scaffold}{$locus}{"repeatClass"};
        my $target = &target($locus{$scaffold}{$locus}{"name"}, $locus{$scaffold}{$locus}{"target"});
        my $strand = '+';
        $strand = '-' if $locus{$scaffold}{$locus}{"strand"} eq 'C';
        print GFF3 "$scaffold\tRepeatMasker\tsimilarity\t$start_end[0]\t$start_end[1]\t\.\t$strand\t\.\tName=$name;Target=$target;\n";

        $repeatRegion{"$scaffold\t$start_end[0]\t$start_end[1]"} = 1;

        if ( $name =~ m/(\S+)\/(\S+)/ ) {
            $stats1{$1}{$2}{"num"} ++;
            $stats1{$1}{$2}{"regions"}{"$scaffold\t$start_end[0]\t$start_end[1]"} = 1;
        }
        else {
            $stats2{$name}{"num"} ++;
            $stats2{$name}{"regions"}{"$scaffold\t$start_end[0]\t$start_end[1]"} = 1;
        }
    }
}

my @repeatRegion = keys %repeatRegion;
my $repeat_bp_number = &collect_length(@repeatRegion);
print STDERR "Result files:\n\tgenome.repeat.out, genome.repeat.gff3.\n";
print "total $repeat_bp_number bp sites were identified to be repeat\n\n";

foreach my $class (sort keys %stats1) {
    print "$class:\n";
    foreach my $subclass (sort keys %{$stats1{$class}}) {
        print "\t$subclass\t$stats1{$class}{$subclass}{\"num\"}\t";
        my @length = keys %{$stats1{$class}{$subclass}{"regions"}};
        my $length = &collect_length(@length);
        print "$length\n";
    }
}
print "\n\n";
foreach my $class (sort keys %stats2) {
    print "$class\t$stats2{$class}{\"num\"}\t";
    my @length = keys %{$stats2{$class}{"regions"}};
    my $length = &collect_length(@length);
    print "$length\n";
}

sub collect_length {
    my (@region, %sort1, %sort2);
    foreach (@_) {
        my @field = split /\t/;
        $sort1{$_} = $field[0];
        if ($field[1] <= $field[2]) {
            push @region, "$field[0]\t$field[1]\t$field[2]";
            $sort2{$_} = $field[1];
        }
        else {
            push @region, "$field[0]\t$field[2]\t$field[1]";
            $sort2{$_} = $field[2];
        }
    }
    @region = sort {$sort1{$a} cmp $sort1{$b} or $sort2{$a} <=> $sort2{$b}} @region;

    my $length;
    my $first = shift @region;
    my ($seq, $start, $end) = split /\t/, $first;
    $length += $end - $start + 1;
    foreach (@region) {
        @_ = split /\t/;
        if ($_[0] eq $seq) {
            if ($_[1] > $end) {
                $length +=  $_[2] - $_[1] + 1;
                $end = $_[2];
            }
            elsif ($_[2] > $end) {
                $length += $_[2] - $end;
                $end = $_[2];
            }
        }
        else {
            $length += $_[2] - $_[1] + 1;
            $end = $_[2];
            $seq = $_[0];
        }
    }
    return $length;
}

sub target {
    my $input = $_[1];
    $input =~ s/\(\d+\)//;
    $input =~ s/^\s*//;
    $input =~ s/\s*$//;
    $input =~ s/\s+/ /g;
    my $out = "$_[0] $input";
    return $out;
}
