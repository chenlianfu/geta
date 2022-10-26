#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Cwd qw/abs_path getcwd cwd/;
use File::Basename;

my $usage = <<USAGE;
Usage:
    $0 [options] input.gtf > output.gff3

    本程序用于将GTF格式转换为GFF3格式。程序会忽略不包含gene_id和transcript_id的行；程序默认忽略不包含CDS信息的mRNA或gene，若需要保留非编码RNA或gene，请注意添加 --keep_NonCDS 参数。

    --gene_prefix <string>    default: None
    若设置该参数，则程序会对基因ID进行重命名，该参数用于设置gene ID前缀。若不设置该参数，则程序不会对基因进行重命名。

    --gene_code_length <int>    default: none
    设置基因数字编号的长度。若不添加改参数，则程序根据基因的数量自动计算出基因数字编号的长度。例如，基因总数量为10000~99999时，基因数字编号长度为5，于是第一个基因编号为00001；若基因总数量在1000~9999时，基因数字编号长度为4，第一个基因编号为0001。设置该参数用于强行指定基因编号的长度，从而决定基因编号前0的数量。若某个基因的编号数值长度 >= 本参数设置的值，则该基因编号数值前不加0。

    --keep_NonCDS    default: None
    程序默认情况下会去除所有不包含CDS信息的mRNA或Gene，若需要保留非编码RNA或gene，请添加本参数。

    --help    default: None
    display this help and exit.

USAGE
if (@ARGV==0){die $usage}

my ($gene_prefix, $gene_code_length, $keep_NonCDS, $help_flag);
GetOptions(
    "gene_prefix:s" => \$gene_prefix,
    "gene_code_length:i" => \$gene_code_length,
    "keep_NonCDS" => \$keep_NonCDS,
    "help" => \$help_flag,
);

if ( $help_flag ) { die $usage }

my (%gtf_info, %lines, %geneSort1, %geneSort2, %geneSort3, %geneSort4, %geneSort5, %geneExonPos, %geneCDSPos, %source, %score, %gff3_attr, %mRNAExonPos);
open IN, $ARGV[0] or die "Can not open file $ARGV[0], $!";
while (<IN>) {
    next if /^#/;
    next if /^\s/;
    next if exists $lines{$_};
    $lines{$_} = 1;

    @_ = split /\t/, $_;
    my $attr = pop @_;

    # 获取该行的gene_id和transcript_id
    my $line_ok = 1;
    my ($gene_id, $transcript_id);
    if ( $attr =~ s/gene_id \"(.*?)\";?// ) {
        $gene_id = $1;
    }
    else {
        $line_ok = 0;
    }
    if ( $attr =~ s/transcript_id \"(.*?)\";?// ) {
        $transcript_id = $1;
    }
    else {
        $line_ok = 0;
    }
    next unless $line_ok;

    # 获取 gene_id 和 transcript_id 的其它属性
    if ( $_[2] eq "gene" && $attr =~ m/\S+/ ) {
        my $attr_other;
        while ( $attr =~ s/(\S+) \"(.*?)\";?// ) {
            my ($tag, $value) = ($1, $2);
            $value =~ s/\s+/\%20/g;
            $attr_other .= "$tag=$value;";
        }
        $gff3_attr{$gene_id} = $attr_other;
    }
    if ( ($_[2] eq "mRNA" or $_[2] eq "transcript") && $attr =~ m/\S+/ ) {
        my $attr_other;
        while ( $attr =~ s/(\S+) \"(.*?)\";?// ) {
            my ($tag, $value) = ($1, $2);
            $value =~ s/\s+/\%20/g;
            $attr_other .= "$tag=$value;";
        }
        $gff3_attr{$transcript_id} = $attr_other;
    }

    # 得到 gene_id / transcript_id 的信息
    $gtf_info{$gene_id}{$transcript_id} .= $_;

    # 得到对gene进行排序的数据
    unless ( exists $geneSort1{$gene_id} ) {
        $geneSort1{$gene_id} = $_[0];
        $geneSort4{$gene_id} = $_[6];
    }
    $geneExonPos{$gene_id}{$_[3]} = 1;
    $geneExonPos{$gene_id}{$_[4]} = 1;
    $mRNAExonPos{$transcript_id}{$_[3]} = 1;
    $mRNAExonPos{$transcript_id}{$_[4]} = 1;
    $geneCDSPos{$gene_id}{$_[3]} = 1 if $_[2] eq "CDS";

    # 得到基因的source和score信息
    $source{$gene_id} = $_[1] if $_[2] eq "gene";
    $score{$gene_id} = $_[5] if $_[2] eq "gene";
    $source{$transcript_id} = $_[1] if ($_[2] eq "mRNA" or $_[2] eq "transcript");
    $score{$transcript_id} = $_[5] if ($_[2] eq "mRNA" or $_[2] eq "transcript");
}

# 得到基因按位置进行排序的数据
foreach my $gene_id ( keys %geneExonPos ) {
    my @pos = sort {$a <=> $b} keys %{$geneExonPos{$gene_id}};
    $geneSort2{$gene_id} = $pos[0];
    $geneSort3{$gene_id} = $pos[-1];
}
foreach my $gene_id ( keys %geneCDSPos ) {
    my @pos = sort {$a <=> $b} keys %{$geneCDSPos{$gene_id}};
    $geneSort5{$gene_id} = $pos[0];
}

# 对基因按基因组序列名、exon首尾位置、正负链和CDS首部位置进行排序。
my @gene_id = sort { $geneSort1{$a} cmp $geneSort1{$b} or $geneSort2{$a} <=> $geneSort2{$b} or $geneSort3{$a} <=> $geneSort3{$b} or $geneSort4{$a} cmp $geneSort4{$b}  or $geneSort5{$a} <=> $geneSort5{$b} } keys %gtf_info;
$gene_code_length ||= length(@gene_id);

my $geneNum = 0;
foreach my $gene_id ( @gene_id ) {
    # 默认情况下，程序忽略不包含CDS信息的gene。
    unless ($keep_NonCDS) {
        next unless exists $geneCDSPos{$gene_id};
    }
    # 得到输出GFF3文件中的gene_id信息。
    my $gene_name = $gene_id;
    if ( $gene_prefix ) {
        $geneNum ++;
        $gene_name = $gene_prefix . '0' x ($gene_code_length - length($geneNum)) . $geneNum;
    }

    # 输出GFF3文件的 gene feature 信息。
    my ($chr, $strand) = ($geneSort1{$gene_id}, $geneSort4{$gene_id});
    $source{$gene_id} = '.' unless $source{$gene_id};
    $score{$gene_id} = '.' unless $score{$gene_id};
    print "$chr\t$source{$gene_id}\tgene\t$geneSort2{$gene_id}\t$geneSort3{$gene_id}\t$score{$gene_id}\t$strand\t\.\tID=$gene_name;$gff3_attr{$gene_id}\n";

    # 对 mRNA 进行解析
    my $mRNA_number = 0;
    foreach my $mRNA_ID ( sort keys %{$gtf_info{$gene_id}} ) {
        my $mRNA_info = $gtf_info{$gene_id}{$mRNA_ID};

        # 默认情况下，程序忽略不包含CDS信息的mRNA。
        unless ($keep_NonCDS) {
            next unless $mRNA_info =~ m/\tCDS\t/;
        }

        # 输出GFF3文件的 mRNA feature 信息。
        $mRNA_number ++;
        my $mRNAID = $mRNA_ID;
        $mRNAID = "$gene_name.t$mRNA_number" if $gene_prefix;
        $source{$mRNA_ID} = '.' unless $source{$mRNA_ID};
        $score{$mRNA_ID} = '.' unless $score{$mRNA_ID};
        my @mRNAExonPos = sort {$a <=> $b} keys %{$mRNAExonPos{$mRNA_ID}};
        my $source = $source{$mRNA_ID};
        print "$chr\t$source\tmRNA\t$mRNAExonPos[0]\t$mRNAExonPos[-1]\t$score{$mRNA_ID}\t$strand\t\.\tID=$mRNAID;Parent=$gene_name;$gff3_attr{$mRNA_ID}\n";

        # 获取 CDS、exon、intron 和 UTR 信息。
        my (@CDS, @exon, @intron, @UTR);
        foreach ( split /\n/, $mRNA_info ) {
            @_ = split /\t/, $_;
            push @CDS, "$_[3]\t$_[4]\t$_[5]\t$_[6]\t$_[7]" if $_[2] eq "CDS";
            push @exon, "$_[3]\t$_[4]" if $_[2] eq "exon";
            push @intron, "$_[3]\t$_[4]" if $_[2] eq "intron";
            push @UTR, "five_prime_UTR\t$_[3]\t$_[4]" if $_[2] eq "5UTR";
            push @UTR, "three_prime_UTR\t$_[3]\t$_[4]" if $_[2] eq "3UTR";
        }
        @CDS = sort {$a <=> $b} @CDS;
        @exon = sort {$a <=> $b} @exon;
        # 若没有exon信息，则将CDS信息给予exon信息。
        unless ( @exon ) {
            foreach (@CDS) {
                @_ = split /\t/, $_;
                push @exon, "$_[0]\t$_[1]";
            }
        }
        # 若没有intron信息，则计算得到intron信息。
        @intron = &get_intron(\@exon, $mRNA_ID, 1) unless @intron;
        # 若没有UTR信息，则计算得到UTR信息。
        @UTR = &get_UTR(\@CDS, \@exon, $strand);

        # 输出转录本数据
        my (%sort, %sort_UTR, $CDS_num, $exon_num, $intron_num, $UTR3_num, $UTR5_num);
        if ($strand eq "+") {
            foreach (sort {$a <=> $b} @CDS) {
                $CDS_num ++;
                my $out = "$chr\t$source\tCDS\t$_\tID=$mRNAID.CDS$CDS_num;Parent=$mRNAID;\n";
                @_ = split /\t/;
                $sort{$out} = $_[0];
                $sort_UTR{$out} = 3;
            }
            foreach (sort {$a <=> $b} @exon) {
                $exon_num ++;
                my $out = "$chr\t$source\texon\t$_\t.\t$strand\t\.\tID=$mRNAID.exon$exon_num;Parent=$mRNAID;\n";
                @_ = split /\t/;
                $sort{$out} = $_[0];
                $sort_UTR{$out} = 2;
            }
            foreach (sort {$a <=> $b} @intron) {
                $intron_num ++;
                my $out = "$chr\t$source\tintron\t$_\t\.\t$strand\t\.\tID=$mRNAID.intron$intron_num;Parent=$mRNAID;\n";
                @_ = split /\t/;
                $sort{$out} = $_[0];
            }
            my (%UTR_sort, @UTR5, @UTR3);
            foreach (@UTR) {
                @_ = split /\t/;
                $UTR_sort{$_} = $_[1];
                push @UTR5, $_ if $_[0] eq "five_prime_UTR";
                push @UTR3, $_ if $_[0] eq "three_prime_UTR";
            }
            foreach (sort {$UTR_sort{$a} <=> $UTR_sort{$b}} @UTR5) {
                $UTR5_num ++;
                my $out = "$chr\t$source\t$_\t.\t$strand\t\.\tID=$mRNAID.utr5p$UTR5_num;Parent=$mRNAID;\n";
                @_ = split /\t/;
                $sort{$out} = $_[1];
                $sort_UTR{$out} = 1;
            }
            foreach (sort {$UTR_sort{$a} <=> $UTR_sort{$b}} @UTR3) {
                $UTR3_num ++;
                my $out = "$chr\t$source\t$_\t.\t$strand\t\.\tID=$mRNAID.utr3p$UTR3_num;Parent=$mRNAID;\n";
                @_ = split /\t/;
                $sort{$out} = $_[1];
                $sort_UTR{$out} = 4;
            }

            foreach (sort {$sort{$a} <=> $sort{$b} or $sort_UTR{$a} <=> $sort_UTR{$b}} keys %sort) {
                print;
            }
        }
        elsif ($strand eq "-") {
            foreach (sort {$b <=> $a} @CDS) {
                $CDS_num ++;
                my $out = "$chr\t$source\tCDS\t$_\tID=$mRNAID.CDS$CDS_num;Parent=$mRNAID;\n";
                @_ = split /\t/;
                $sort{$out} = $_[0];
                $sort_UTR{$out} = 3;
            }
            foreach (sort {$b <=> $a} @exon) {
                $exon_num ++;
                my $out = "$chr\t$source\texon\t$_\t.\t$strand\t\.\tID=$mRNAID.exon$exon_num;Parent=$mRNAID;\n";
                @_ = split /\t/;
                $sort{$out} = $_[0];
                $sort_UTR{$out} = 2;
            }
            foreach (sort {$b <=> $a} @intron) {
                $intron_num ++;
                my $out = "$chr\t$source\tintron\t$_\t\.\t$strand\t\.\tID=$mRNAID.intron$intron_num;Parent=$mRNAID;\n";
                @_ = split /\t/;
                $sort{$out} = $_[0];
            }
            my (%UTR_sort, @UTR5, @UTR3);
            foreach (@UTR) {
                @_ = split /\t/;
                $UTR_sort{$_} = $_[1];
                push @UTR5, $_ if $_[0] eq "five_prime_UTR";
                push @UTR3, $_ if $_[0] eq "three_prime_UTR";
            }
            foreach (sort {$UTR_sort{$b} <=> $UTR_sort{$a}} @UTR5) {
                $UTR5_num ++;
                my $out = "$chr\t$source\t$_\t.\t$strand\t\.\tID=$mRNAID.utr5p$UTR5_num;Parent=$mRNAID;\n";
                @_ = split /\t/;
                $sort{$out} = $_[1];
                $sort_UTR{$out} = 1;
            }
            foreach (sort {$UTR_sort{$a} <=> $UTR_sort{$b}} @UTR3) {
                $UTR3_num ++;
                my $out = "$chr\t$source\t$_\t.\t$strand\t\.\tID=$mRNAID.utr3p$UTR3_num;Parent=$mRNAID;\n";
                @_ = split /\t/;
                $sort{$out} = $_[1];
                $sort_UTR{$out} = 4;
            }

            foreach (sort {$sort{$b} <=> $sort{$a} or $sort_UTR{$a} <=> $sort_UTR{$b}} keys %sort) {
                print;
            }
        }
        print "\n";

        #print "$gtf_info{$gene_id}{$mRNA_ID}\n";
    }
}

sub get_UTR {
    my @cds = @{$_[0]};
    my @exon = @{$_[1]};
    my $strand = $_[2];

    my (@utr, %cds_pos);
    foreach (@cds) {
        @_ = split /\t/;
        $cds_pos{$_[0]} = 1;
        $cds_pos{$_[1]} = 1;
    }

    foreach (@exon) {
        my ($start, $end) = split /\t/;
        my $utr_keep = 1;
        foreach (@cds) {
            @_ = split /\t/;
            if ($_[0] <= $end && $_[1] >= $start) {
                $utr_keep = 0;
                if ($start < $_[0] && $end == $_[1]) {
                    my $utr_start = $start;
                    my $utr_end = $_[0] - 1;
                    push @utr, "$utr_start\t$utr_end";
                }
                elsif ($start == $_[0] && $end > $_[1]) {
                    my $utr_start = $_[1] + 1;
                    my $utr_end = $end;
                    push @utr, "$utr_start\t$utr_end";
                }
            }
        }
        push @utr, $_ if $utr_keep == 1;
    }

    my @out;
    my @cds_pos = sort {$a <=> $b} keys %cds_pos;
    if ($strand eq "+") {
        @utr = sort {$a <=> $b} @utr;
        foreach (@utr) {
            @_ = split /\t/;
            if ($_[1] <= $cds_pos[0]) {
                push @out, "five_prime_UTR\t$_";
            }
            elsif ($_[0] >= $cds_pos[1]) {
                push @out, "three_prime_UTR\t$_";
            }
        }
    }
    elsif ($strand eq "-") {
        @utr = sort {$b <=> $a} @utr;
        foreach (@utr) {
            @_ = split /\t/;
            if ($_[0] >= $cds_pos[1]) {
                push @out, "five_prime_UTR\t$_";
            }
            elsif ($_[1] <= $cds_pos[0]) {
                push @out, "three_prime_UTR\t$_";
            }
        }
    }

    return @out;
}

sub get_intron {
    my @exon = @{$_[0]};
    my $mRNA_ID = $_[1];
    my $intron_len = $_[2];
    @exon = sort {$a <=> $b} @exon;

    my @intron;
    my $first_exon = shift @exon;
    my ($last_start, $last_end) = split /\t/, $first_exon;
    foreach ( @exon ) {
        my ($start, $end) = split /\t/, $_;
        if ($start > $last_end + $intron_len) {
            my $intron_start = $last_end + 1;
            my $intron_stop = $start - 1;
            push @intron, "$intron_start\t$intron_stop";
        }
        else {
            my $value = $start - $last_end - 1;
            print STDERR "Warning: a intron length (value is $value) of mRNA $mRNA_ID < $intron_len was detected:\n\tThe former CDS/Exon: $last_start - $last_end\n\tThe latter CDS/Exon: $start - $end\n";
        }
        ($last_start, $last_end) = ($start, $end);
    }

    return @intron;
}
