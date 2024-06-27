#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 strandard.gff3 query.gff3 > out.txt 2> Consistent_gene_or_CDS.tab

    This program compares the gene models in two GFF3 files, using the gene models from the first input GFF3 file as a reference to check the accuracy of the gene models in the second GFF3 file.

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
my (%gene_info, %CDS_info, %index, $total_gene_number, %CDS_ref);
foreach my $gene_id (keys %gene) {
    # 计算ref中基因的总数
    my $gene_info = $feature{$gene_id};
    $total_gene_number ++ if $gene_info =~ m/\tmRNA\t/;
    #print "OK:\t$gene_id\n";
    
    # 分析基因的CDS信息，给所有可变剪接转录本的CDS去冗余获得，存放到 %CDS_info 中。
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

    # 给基因位置进行索引
    my $gene_desc = $gene{$gene_id};
    @_ = split /\t/, $gene_desc;
    my ($chr, $strand) = ($_[0], $_[6]);
    my $int1 = int($_[3] / 1000);
    my $int2 = int($_[4] / 1000);
    foreach ($int1 .. $int2) {
        $index{$_[0]}{$_[6]}{$_}{$gene_id} = 1;
    }

    # 得到ref中所有CDS信息。
    foreach ( keys %{$CDS_info{$gene_id}} ) {
        $CDS_ref{$chr}{$strand}{$_} = 1;
    }
}

my ($total_CDS_number, $total_bp_number) = &cal_CDS_and_bp_num(\%CDS_ref);

#print "$total_gene_number\t$total_CDS_number\t$total_bp_number\n";

# 进行positive分析
my ($query_gene_number, %positive_gene, %CDS_query, %CDS_same, %CDS_overlap);
foreach my $gene_id (keys %gene_query) {
    # 分析当前基因的位置
    my $gene_desc = $gene_query{$gene_id};
    @_ = split /\t/, $gene_desc;
    my ($chr, $strand) = ($_[0], $_[6]);
    my $int1 = int($_[3] / 1000);
    my $int2 = int($_[4] / 1000);
    #print "OK2:\t$gene_id\t$int1\t$int2\n";
    # 寻找重叠的ref基因ID
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

    # 检测query中基因的数量
    my $gene_info = $feature_query{$gene_id};
    if ($gene_info =~ m/\tmRNA\t/) {
        $query_gene_number ++;
    }
    else {
        next;
    }

    # 分析当前基因中所有转录本的CDS信息
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
                     $CDS_query{$chr}{$strand}{"$_[3]\t$_[4]"} = 1;
                 }
             }
             @CDS = sort {$a <=> $b} @CDS;
             my $CDS = join "\n", @CDS;
             $transcript_CDS{$CDS} = 1;
        }
    }

    # 检测是否有转录本的CDS信息和ref中的转录本一致。
    my $equal_ref_gene_id;
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
                if ($strand eq "+") {
                    $CDS_target =~ s/\d+//;
                }
                elsif ($strand eq "-") {
                    $CDS_target =~ s/\d+$//;
                }

                if ($transcript_CDS_for_validation eq $CDS_target) {
                    $same = 1;
                    $equal_ref_gene_id = $target_gene_id;

                    print STDERR "Same gene: $gene_id\t$target_gene_id\n";
                    last CDS_COMPARE;
                }
            }
        }
    }
    $positive_gene{$equal_ref_gene_id} = 1 if $same == 1;

    # 寻找共同或重叠的CDS
    foreach my $cds (keys %CDS) {
        foreach (keys %target_CDS) {
            @_ = split /\t/;
            my @cds = split /\t/, $cds;
            if ($_[0] <= $cds[1] && $cds[0] <= $_[1]) {
                my @cds_region;
                push @cds_region, @_;
                push @cds_region, @cds;
                @cds_region = sort {$a <=> $b} @cds_region;
                $CDS_overlap{$chr}{$strand}{"$cds_region[1]\t$cds_region[2]"} = 1;
            }

            if ($_ eq $cds) {
                $CDS_same{$chr}{$strand}{$cds} = 1;
                print STDERR "Same CDS: $chr\t$strand\t$cds\n";
                last;
            }
        }
    }
}

my $positive_gene_number = %positive_gene;
my $positive_CDS_number = &cal_CDS_num(\%CDS_same);
my ($query_CDS_number, $query_bp_number) = &cal_CDS_and_bp_num(\%CDS_query);
my ($positive_CDS_overlap_num, $positive_bp_number) = &cal_CDS_and_bp_num(\%CDS_overlap);
my $sensitivity_gene = $positive_gene_number  * 100 / $total_gene_number;
my $sensitivity_CDS = $positive_CDS_number * 100 / $total_CDS_number;
my $sensitivity_bp = $positive_bp_number * 100 / $total_bp_number;
my $specificity_gene = $positive_gene_number * 100 / $query_gene_number;
my $specificity_CDS = $positive_CDS_number * 100 / $query_CDS_number;
my $specificity_bp = $positive_bp_number * 100 / $query_bp_number;

print "gene_prediction_method\tgene_standard\tgene_predicted\tgene_TP\tgene_sensitivity\tgene_specificity\tCDS_standard\tCDS_predicted\tCDS_TP\tCDS_sensitivity\tCDS_specificity\tbp_standard\tbp_predicted\tbp_TP\tbp_sensitivity\tbp_specificity\n";
printf "@ARGV[1]\t$total_gene_number\t$query_gene_number\t$positive_gene_number\t%.2f\%\t%.2f\%\t$total_CDS_number\t$query_CDS_number\t$positive_CDS_number\t%.2f\%\t%.2f\%\t$total_bp_number\t$query_bp_number\t$positive_bp_number\t%.2f\%\t%.2f\%\n", $sensitivity_gene, $specificity_gene, $sensitivity_CDS, $specificity_CDS, $sensitivity_bp, $specificity_bp;


sub cal_CDS_and_bp_num {
    # 输入是一个哈希：染色体 -> 正负链 -> 制表符分隔的CDS起始结束 -> 1。
    my %CDS = %{$_[0]};
    my ($CDS_num, $bp_num) = (0, 0);

    foreach my $chr ( keys %CDS ) {
        foreach my $strand ( keys %{$CDS{$chr}} ) {
            my @CDS = sort { $a <=> $b } keys %{$CDS{$chr}{$strand}};
            $CDS_num += @CDS;

            my ($start, $end) = split /\t/, $CDS[0];
            $bp_num += $end - $start + 1;
            shift @CDS;
            foreach ( @CDS ) {
                @_ = split /\t/, $_;
                if ($_[0] > $end) {
                    $bp_num += $_[1] - $_[0] + 1;
                    ($start, $end) = @_;
                }
                else {
                    if ($_[1] > $end) {
                        $bp_num += $_[1] - $end;
                        $end = $_[1];
                    }
                }
            }

        }
    }

    return ($CDS_num, $bp_num);
}

sub cal_CDS_num {
    my %CDS = %{$_[0]};
    my $CDS_num = 0;

    foreach my $chr ( keys %CDS ) {
        foreach my $strand ( keys %{$CDS{$chr}} ) {
            my @CDS = keys %{$CDS{$chr}{$strand}};
            $CDS_num += @CDS;
        }
    }

    return $CDS_num;
}
