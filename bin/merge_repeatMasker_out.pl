#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    perl $0 genome.fasta repeatMasker/genome.fasta.out repeatModeler/genome.fasta.out > genome.repeat.stats

    --min_coverge_ratio <float>    default: 0.25
    设置最小覆盖率阈值。若匹配区域对目标重复序列的覆盖率低于此阈值，则过滤掉此结果。

USAGE
if (@ARGV==0){die $usage}
my $min_coverge_ratio;
GetOptions(
    "min_coverge_ratio:f" => \$min_coverge_ratio,
);
$min_coverge_ratio ||= 0.25;

my $genome_file = shift @ARGV;
my $genome_size;
open IN, $genome_file or die "Can not open file $genome_file, $!";
while (<IN>) {
    next if m/^>/;
    $genome_size += (length($_) - 1);
}
close IN;

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

        my $ratio = 0;
        if ($aaa[12] =~ m/\((\d+)\)/) {
            $ratio = ($aaa[13] - $aaa[14] + 1) / ($1 + $aaa[13]) if ($1 + $aaa[13]) != 0;
        }
        elsif ($aaa[14] =~ m/\((\d+)\)/) {
            $ratio = ($aaa[13] - $aaa[12] + 1) / ($aaa[13] + $1) if ($aaa[13] + $1) != 0;
        }
        next if $ratio < $min_coverge_ratio;
        $ratio = sprintf("%.4f", $ratio);

        $locus{$aaa[5]}{"$aaa[6]\t$aaa[7]"}{"description"} = $_;
        $locus{$aaa[5]}{"$aaa[6]\t$aaa[7]"}{"start"} = $aaa[6];
        $locus{$aaa[5]}{"$aaa[6]\t$aaa[7]"}{"strand"} = $aaa[9];
        $locus{$aaa[5]}{"$aaa[6]\t$aaa[7]"}{"name"} = $aaa[10];
        $locus{$aaa[5]}{"$aaa[6]\t$aaa[7]"}{"repeatClass"}  = $aaa[11];
        $locus{$aaa[5]}{"$aaa[6]\t$aaa[7]"}{"target"} = "$aaa[12]\t$aaa[13]\t$aaa[14]";
        $locus{$aaa[5]}{"$aaa[6]\t$aaa[7]"}{"ratio"} = $ratio;

        push @scaffold, $aaa[5] unless exists $scaffold{$aaa[5]};
        $scaffold{$aaa[5]} = 1;
    }
}
close IN;

# %stats1 统计子类的重复序列； %stats2统计大类的重复序列。
my (%stats1, %stats2, %Retroelements, %DNAelements, %Interspersed, %nonInterspersed, %total_repeat);
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
        my $Ratio = $locus{$scaffold}{$locus}{"ratio"};
        print GFF3 "$scaffold\tRepeatMasker\tsimilarity\t$start_end[0]\t$start_end[1]\t\.\t$strand\t\.\tName=$name;Target=$target;Ratio=$Ratio\n";

        $total_repeat{"num"} ++;
        $total_repeat{"regions"}{"$scaffold\t$start_end[0]\t$start_end[1]"} = 1;

        if ($name =~ m/LTR|SINE|LINE|Retroposon/) {
            $Retroelements{"num"} ++;
            $Retroelements{"regions"}{"$scaffold\t$start_end[0]\t$start_end[1]"} = 1;
        }
        elsif ($name =~ m/DNA/) {
            $DNAelements{"num"} ++;
            $DNAelements{"regions"}{"$scaffold\t$start_end[0]\t$start_end[1]"} = 1;
        }

        if ($name =~ m/Low_complexity|Satellite|Simple_repeat|rRNA/) {
            $nonInterspersed{"num"} ++;
            $nonInterspersed{"regions"}{"$scaffold\t$start_end[0]\t$start_end[1]"} = 1;
        }
        else {
            $Interspersed{"num"} ++;
            $Interspersed{"regions"}{"$scaffold\t$start_end[0]\t$start_end[1]"} = 1;
        }

        if ( $name =~ m/(\S+)\/(\S+)/ ) {
            $stats1{$1}{$2}{"num"} ++;
            $stats1{$1}{$2}{"regions"}{"$scaffold\t$start_end[0]\t$start_end[1]"} = 1;
            $stats2{$1}{"num"} ++;
            $stats2{$1}{"regions"}{"$scaffold\t$start_end[0]\t$start_end[1]"} = 1;
        }
        else {
            $stats2{$name}{"num"} ++;
            $stats2{$name}{"regions"}{"$scaffold\t$start_end[0]\t$start_end[1]"} = 1;
        }
    }
}

my @repeatRegion = keys %{$total_repeat{"regions"}};
my $total_repeat_num = $total_repeat{"num"};
my $total_repeat_size = &collect_length(@repeatRegion);
print STDERR "Result files:\n\tgenome.repeat.out, genome.repeat.gff3.\n";
my $total_repeat_ratio = int($total_repeat_size / $genome_size * 10000) / 100;
print "=======================================================\n";
printf "Genome Size: $genome_size bp\nRepeat Size: $total_repeat_size bp\nRepeat Ratio: %.2f\%\n", $total_repeat_ratio;
print "=======================================================\n";
print "Class\tNumber\tLength\tPercentage\n";
print "-------------------------------------------------------\n";
printf "Total\t$total_repeat_num\t$total_repeat_size\t%.2f\%\n", $total_repeat_ratio;

# 输出非散在重复序列的统计信息
my $nonInterspersed_repeat_num = $nonInterspersed{"num"};
my @nonInterspersed_region = keys %{$nonInterspersed{"regions"}};
my $nonInterspersed_size = &collect_length(@nonInterspersed_region);
my $nonInterspersed_ratio = int($nonInterspersed_size / $genome_size * 10000) / 100;
printf "    Non-interspersed Repeats\t$nonInterspersed_repeat_num\t$nonInterspersed_size\t%.2f\%\n", $nonInterspersed_ratio;

# 输出非散在重复序列的子类统计信息
my @class_nonInterspersed = qw/Simple_repeat Low_complexity Satellite rRNA/;
my %class_nonInterspersed_sort;
foreach (@class_nonInterspersed) {
    my @region = keys %{$stats2{$_}{"regions"}};
    $class_nonInterspersed_sort{$_} = &collect_length(@region);
}
foreach my $class (sort {$class_nonInterspersed_sort{$b} <=> $class_nonInterspersed_sort{$a}} keys %class_nonInterspersed_sort) {
    print "        $class\t$stats2{$class}{\"num\"}\t";
    my @length = keys %{$stats2{$class}{"regions"}};
    my $length = &collect_length(@length);
    my $ratio = int($length / $genome_size * 10000) / 100;
    printf "$length\t%.2f\%\n", $ratio;

    foreach my $subclass (sort keys %{$stats1{$class}}) {
        print "            $subclass\t$stats1{$class}{$subclass}{\"num\"}\t";
        my @length = keys %{$stats1{$class}{$subclass}{"regions"}};
        my $length = &collect_length(@length);
        my $ratio = int($length / $genome_size * 10000) / 100;
        printf "$length\t%.2f\%\n", $ratio;
    }
}

# 输出散在重复序列的统计信息
my $Interspersed_repeat_num = $Interspersed{"num"};
my @Interspersed_region = keys %{$Interspersed{"regions"}};
my $Interspersed_size = &collect_length(@Interspersed_region);
my $Interspersed_ratio = int($Interspersed_size / $genome_size * 10000) / 100;
printf "    Interspersed Repeats\t$Interspersed_repeat_num\t$Interspersed_size\t%.2f\%\n", $Interspersed_ratio;

# 输出逆转录转座子的统计信息
my $Retroelements_num = $Retroelements{"num"};
my @Retroelements_region = keys %{$Retroelements{"regions"}};
my $Retroelements_size = &collect_length(@Retroelements_region);
my $Retroelements_ratio = int($Retroelements_size / $genome_size * 10000) / 100;
printf "        Retroelements\t$Retroelements_num\t$Retroelements_size\t%.2f\%\n", $Retroelements_ratio;

# 输出逆转录转座子的子类统计信息
my @class = sort keys %stats2;
my %for_sort;
foreach (@class) {
    next unless m/LTR|SINE|LINE|Retroposon/;
    if ($_ eq "LTR") {
        $for_sort{$_} = 1;
    }
    elsif ($_ eq "LINE") {
        $for_sort{$_} = 2;
    }
    elsif ($_ eq "SINE") {
        $for_sort{$_} = 3;
    }
    elsif ($_ eq "Retroposon") {
        $for_sort{$_} = 4;
    }
    else {
        $for_sort{$_} = 5;
    }
}
@class = sort {$for_sort{$a} <=> $for_sort{$b} or $a cmp $b} keys %for_sort;
foreach my $class (@class) {
    print "            $class\t$stats2{$class}{\"num\"}\t";
    my @length = keys %{$stats2{$class}{"regions"}};
    my $length = &collect_length(@length);
    my $ratio = int($length / $genome_size * 10000) / 100;
    printf "$length\t%.2f\%\n", $ratio;

    my @subclass = keys %{$stats1{$class}};
    my %subclass_sort;
    foreach (@subclass) {
        my @length = keys %{$stats1{$class}{$_}{"regions"}};
        $subclass_sort{$_} = &collect_length(@length);
    }
    @subclass = sort {$subclass_sort{$b} <=> $subclass_sort{$a}} @subclass;
    foreach my $subclass (@subclass) {
        print "                $subclass\t$stats1{$class}{$subclass}{\"num\"}\t";
        my @length = keys %{$stats1{$class}{$subclass}{"regions"}};
        my $length = &collect_length(@length);
        my $ratio = int($length / $genome_size * 10000) / 100;
        printf "$length\t%.2f\%\n", $ratio;
    }
}

# 输出DNA转座子的统计信息
my $DNAelements_num = $DNAelements{"num"};
my @DNAelements_region = keys %{$DNAelements{"regions"}};
my $DNAelements_size = &collect_length(@DNAelements_region);
my $DNAelements_ratio = int($DNAelements_size / $genome_size * 10000) / 100;
printf "        DNA transposons\t$DNAelements_num\t$DNAelements_size\t%.2f\%\n", $DNAelements_ratio;

# 输出DNA转座子的子类统计信息
my @class_all = sort keys %stats2;
my @class;
foreach (@class_all) {
    push @class, $_ if m/DNA/;
    push @class, $_ if m/ARTEFACT/;
}
foreach my $class (@class) {
    if (exists $stats1{$class}) {
        my @subclass = keys %{$stats1{$class}};
        my %subclass_sort;
        foreach (@subclass) {
            my @length = keys %{$stats1{$class}{$_}{"regions"}};
            $subclass_sort{$_} = &collect_length(@length);
        }
        @subclass = sort {$subclass_sort{$b} <=> $subclass_sort{$a}} @subclass;
        foreach my $subclass (@subclass) {
            print "            $subclass\t$stats1{$class}{$subclass}{\"num\"}\t";
            my @length = keys %{$stats1{$class}{$subclass}{"regions"}};
            my $length = &collect_length(@length);
            my $ratio = int($length / $genome_size * 10000) / 100;
            printf "$length\t%.2f\%\n", $ratio;
        }
    }
    else {
        print "            $class\t$stats2{$class}{\"num\"}\t";
        my @length = keys %{$stats2{$class}{"regions"}};
        my $length = &collect_length(@length);
        my $ratio = int($length / $genome_size * 10000) / 100;
        printf "$length\t%.2f\%\n", $ratio;
    }
}

# 输出其它散在重复序列的统计信息
my @class = keys %stats2;
my %for_sort;
foreach (@class) {
    next if m/Simple_repeat|Low_complexity|Satellite|rRNA|LTR|SINE|LINE|Retroposon|DNA|ARTEFACT/;
    my @length = keys %{$stats2{$_}{"regions"}};
    $for_sort{$_} = &collect_length(@length);
}
@class = sort {$for_sort{$b} <=> $for_sort{$a}} keys %for_sort;
foreach my $class (@class) {
    my $class_out = $class;
    $class_out = "Rolling-circles" if $class eq "RC";
    print "        $class_out\t$stats2{$class}{\"num\"}\t";
    my @length = keys %{$stats2{$class}{"regions"}};
    my $length = &collect_length(@length);
    my $ratio = int($length / $genome_size * 10000) / 100;
    printf "$length\t%.2f\%\n", $ratio;

    foreach my $subclass (sort keys %{$stats1{$class}}) {
        print "            $subclass\t$stats1{$class}{$subclass}{\"num\"}\t";
        my @length = keys %{$stats1{$class}{$subclass}{"regions"}};
        my $length = &collect_length(@length);
        my $ratio = int($length / $genome_size * 10000) / 100;
        printf "$length\t%.2f\%\n", $ratio;
    }
}
print "=======================================================\n";

print "# 注意事项：
# 1. 某个大类的总个数可能超过其所属小类的总数之和，这是因为该大类下包含一些未能归于其子类的重复序列。
# 2. 重复序列的总碱基数会小于各类型重复序列碱基数之和，这是因为不同类型的重复序列之间可能有重叠。\n";


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
