#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    $0 [options] input.gtf > output.gff3

    本程序用于将GTF格式转换为GFF3格式。程序会忽略不包含gene_id的行，以及不属于gene要素但不包含transcript_id的行；程序默认忽略不包含CDS信息的mRNA或gene，若需要保留非编码RNA或gene，请注意添加 --keep_NonCDS 参数。

    --gene_prefix <string>    default: None
    若设置该参数，则程序会对基因ID进行重命名，该参数用于设置gene ID前缀。若不设置该参数，则程序不会对基因进行重命名。

    --gene_code_length <int>    default: none
    设置基因数字编号的长度。若不添加改参数，则程序根据基因的数量自动计算出基因数字编号的长度。例如，基因总数量为10000~99999时，基因数字编号长度为5，于是第一个基因编号为00001；若基因总数量在1000~9999时，基因数字编号长度为4，第一个基因编号为0001。设置该参数用于强行指定基因编号的长度，从而决定基因编号前0的数量。若某个基因的编号数值长度 >= 本参数设置的值，则该基因编号数值前不加0。

    --keep_NonCDS    default: None
    程序默认情况下会去除所有不包含CDS信息的mRNA或Gene，若需要保留非编码RNA或gene，请添加本参数。

    --help    default: None
    display this help and exit.

USAGE
if (@ARGV == 0) { die $usage }

my ($gene_prefix, $gene_code_length, $keep_NonCDS, $help_flag);
GetOptions(
    "gene_prefix:s"       => \$gene_prefix,
    "gene_code_length:i"  => \$gene_code_length,
    "keep_NonCDS"         => \$keep_NonCDS,
    "help"                => \$help_flag,
);

if ($help_flag) { print $usage; exit 0; }
die $usage unless defined $ARGV[0];

my (%gtf_info, %lines, %geneSort1, %geneSort2, %geneSort3, %geneSort4, %geneSort5,
    %geneExonPos, %geneCDSPos, %source, %score, %gff3_attr, %mRNAExonPos);

open my $IN, '<', $ARGV[0] or die "Can not open file $ARGV[0]: $!\n";
while (<$IN>) {
    next if /^#/;
    next if /^\s/;
    next if exists $lines{$_};
    $lines{$_} = 1;

    my @f = split /\t/, $_;
    next if @f < 9;
    my $attr = pop @f;
    my $type = $f[2];

    # gene_id is mandatory on every kept line
    my $gene_id;
    if ( $attr =~ s/gene_id \"(.*?)\";?// ) {
        $gene_id = $1;
    }
    else {
        next;
    }

    my $transcript_id;
    my $has_tx = ( $attr =~ s/transcript_id \"(.*?)\";?// );
    $transcript_id = $1 if $has_tx;
    next if $type ne 'gene' && !$has_tx;

    # 获取 gene_id / transcript_id 之外的其它属性
    if ( ($type eq 'gene' || $type eq 'mRNA' || $type eq 'transcript') && $attr =~ /\S/ ) {
        my $attr_other = '';
        while ( $attr =~ s/(\S+) \"(.*?)\";?// ) {
            my ($tag, $value) = ($1, $2);
            $value =~ s/\s+/\%20/g;
            $attr_other .= "$tag=$value;";
        }
        if ($type eq 'gene') { $gff3_attr{$gene_id} = $attr_other; }
        else                 { $gff3_attr{$transcript_id} = $attr_other; }
    }

    # 得到对gene进行排序的数据 (chr / strand)，取自第一条遇到的该基因的行
    unless ( exists $geneSort1{$gene_id} ) {
        $geneSort1{$gene_id} = $f[0];
        $geneSort4{$gene_id} = $f[6];
    }
    $source{$gene_id} = $f[1] if $type eq 'gene';
    $score{$gene_id}  = $f[5] if $type eq 'gene';

    next if $type eq 'gene';   # gene 行到此为止，不进入按转录本聚合的数据结构

    # 得到 gene_id / transcript_id 的信息
    $gtf_info{$gene_id}{$transcript_id} .= $_;

    $geneExonPos{$gene_id}{$f[3]} = 1;
    $geneExonPos{$gene_id}{$f[4]} = 1;
    $mRNAExonPos{$transcript_id}{$f[3]} = 1;
    $mRNAExonPos{$transcript_id}{$f[4]} = 1;
    $geneCDSPos{$gene_id}{$f[3]} = 1 if $type eq 'CDS';

    $source{$transcript_id} = $f[1] if ($type eq 'mRNA' or $type eq 'transcript');
    $score{$transcript_id}  = $f[5] if ($type eq 'mRNA' or $type eq 'transcript');
}
close $IN;

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
my @gene_id = sort {
       $geneSort1{$a} cmp $geneSort1{$b}
    or $geneSort2{$a} <=> $geneSort2{$b}
    or $geneSort3{$a} <=> $geneSort3{$b}
    or $geneSort4{$a} cmp $geneSort4{$b}
    or $geneSort5{$a} <=> $geneSort5{$b}
} keys %gtf_info;
$gene_code_length ||= length(scalar @gene_id);

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
        my $pad = $gene_code_length - length($geneNum);
        $pad = 0 if $pad < 0;
        $gene_name = $gene_prefix . ('0' x $pad) . $geneNum;
    }

    # 输出GFF3文件的 gene feature 信息。
    my ($chr, $strand) = ($geneSort1{$gene_id}, $geneSort4{$gene_id});
    $source{$gene_id} = '.' unless defined $source{$gene_id};
    $score{$gene_id}  = '.' unless defined $score{$gene_id};
    $gff3_attr{$gene_id} = '' unless defined $gff3_attr{$gene_id};
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
        $source{$mRNA_ID} = '.' unless defined $source{$mRNA_ID};
        $score{$mRNA_ID}  = '.' unless defined $score{$mRNA_ID};
        $gff3_attr{$mRNA_ID} = '' unless defined $gff3_attr{$mRNA_ID};
        my @mRNAExonPos = sort {$a <=> $b} keys %{$mRNAExonPos{$mRNA_ID}};
        my $source = $source{$mRNA_ID};
        print "$chr\t$source\tmRNA\t$mRNAExonPos[0]\t$mRNAExonPos[-1]\t$score{$mRNA_ID}\t$strand\t\.\tID=$mRNAID;Parent=$gene_name;$gff3_attr{$mRNA_ID}\n";

        # 获取 CDS、exon、intron 和 UTR 信息。
        my (@CDS, @exon, @intron, @UTR);
        foreach ( split /\n/, $mRNA_info ) {
            my @g = split /\t/, $_;
            push @CDS,    "$g[3]\t$g[4]\t$g[5]\t$g[6]\t$g[7]" if $g[2] eq "CDS";
            push @exon,   "$g[3]\t$g[4]" if $g[2] eq "exon";
            push @intron, "$g[3]\t$g[4]" if $g[2] eq "intron";
            push @UTR, "five_prime_UTR\t$g[3]\t$g[4]"  if $g[2] eq "5UTR";
            push @UTR, "three_prime_UTR\t$g[3]\t$g[4]" if $g[2] eq "3UTR";
        }
        @CDS  = sort { (split /\t/, $a)[0] <=> (split /\t/, $b)[0] } @CDS;
        @exon = sort { (split /\t/, $a)[0] <=> (split /\t/, $b)[0] } @exon;

        # 若没有exon信息，则用 CDS + 已知UTR 的并集来还原exon范围
        unless ( @exon ) {
            foreach (@CDS) {
                my @g = split /\t/;
                push @exon, "$g[0]\t$g[1]";
            }
            foreach (@UTR) {
                my @g = split /\t/;   # type, start, end
                push @exon, "$g[1]\t$g[2]";
            }
            @exon = merge_intervals(@exon);
        }

        # 若没有intron信息，则计算得到intron信息。
        @intron = &get_intron(\@exon, $mRNA_ID, 1) unless @intron;
        # 若没有UTR信息，则计算得到UTR信息。
        @UTR = &get_UTR(\@CDS, \@exon, $strand) unless @UTR;

        # 输出转录本数据
        my (%sort, %sort_UTR, $CDS_num, $exon_num, $intron_num, $UTR3_num, $UTR5_num);
        if ($strand eq "+") {
            foreach (sort { (split /\t/, $a)[0] <=> (split /\t/, $b)[0] } @CDS) {
                $CDS_num ++;
                my $out = "$chr\t$source\tCDS\t$_\tID=$mRNAID.CDS$CDS_num;Parent=$mRNAID;\n";
                my @g = split /\t/;
                $sort{$out} = $g[0];
                $sort_UTR{$out} = 3;
            }
            foreach (sort { (split /\t/, $a)[0] <=> (split /\t/, $b)[0] } @exon) {
                $exon_num ++;
                my $out = "$chr\t$source\texon\t$_\t.\t$strand\t\.\tID=$mRNAID.exon$exon_num;Parent=$mRNAID;\n";
                my @g = split /\t/;
                $sort{$out} = $g[0];
                $sort_UTR{$out} = 2;
            }
            foreach (sort { (split /\t/, $a)[0] <=> (split /\t/, $b)[0] } @intron) {
                $intron_num ++;
                my $out = "$chr\t$source\tintron\t$_\t\.\t$strand\t\.\tID=$mRNAID.intron$intron_num;Parent=$mRNAID;\n";
                my @g = split /\t/;
                $sort{$out} = $g[0];
                $sort_UTR{$out} = 2;
            }
            my (%UTR_sort, @UTR5, @UTR3);
            foreach (@UTR) {
                my @g = split /\t/;
                $UTR_sort{$_} = $g[1];
                push @UTR5, $_ if $g[0] eq "five_prime_UTR";
                push @UTR3, $_ if $g[0] eq "three_prime_UTR";
            }
            foreach (sort {$UTR_sort{$a} <=> $UTR_sort{$b}} @UTR5) {
                $UTR5_num ++;
                my $out = "$chr\t$source\t$_\t.\t$strand\t\.\tID=$mRNAID.utr5p$UTR5_num;Parent=$mRNAID;\n";
                my @g = split /\t/;
                $sort{$out} = $g[1];
                $sort_UTR{$out} = 1;
            }
            foreach (sort {$UTR_sort{$a} <=> $UTR_sort{$b}} @UTR3) {
                $UTR3_num ++;
                my $out = "$chr\t$source\t$_\t.\t$strand\t\.\tID=$mRNAID.utr3p$UTR3_num;Parent=$mRNAID;\n";
                my @g = split /\t/;
                $sort{$out} = $g[1];
                $sort_UTR{$out} = 4;
            }

            foreach (sort {$sort{$a} <=> $sort{$b} or $sort_UTR{$a} <=> $sort_UTR{$b}} keys %sort) {
                print;
            }
        }
        elsif ($strand eq "-") {
            foreach (sort { (split /\t/, $b)[0] <=> (split /\t/, $a)[0] } @CDS) {
                $CDS_num ++;
                my $out = "$chr\t$source\tCDS\t$_\tID=$mRNAID.CDS$CDS_num;Parent=$mRNAID;\n";
                my @g = split /\t/;
                $sort{$out} = $g[0];
                $sort_UTR{$out} = 3;
            }
            foreach (sort { (split /\t/, $b)[0] <=> (split /\t/, $a)[0] } @exon) {
                $exon_num ++;
                my $out = "$chr\t$source\texon\t$_\t.\t$strand\t\.\tID=$mRNAID.exon$exon_num;Parent=$mRNAID;\n";
                my @g = split /\t/;
                $sort{$out} = $g[0];
                $sort_UTR{$out} = 2;
            }
            foreach (sort { (split /\t/, $b)[0] <=> (split /\t/, $a)[0] } @intron) {
                $intron_num ++;
                my $out = "$chr\t$source\tintron\t$_\t\.\t$strand\t\.\tID=$mRNAID.intron$intron_num;Parent=$mRNAID;\n";
                my @g = split /\t/;
                $sort{$out} = $g[0];
                $sort_UTR{$out} = 2;
            }
            my (%UTR_sort, @UTR5, @UTR3);
            foreach (@UTR) {
                my @g = split /\t/;
                $UTR_sort{$_} = $g[1];
                push @UTR5, $_ if $g[0] eq "five_prime_UTR";
                push @UTR3, $_ if $g[0] eq "three_prime_UTR";
            }
            foreach (sort {$UTR_sort{$b} <=> $UTR_sort{$a}} @UTR5) {
                $UTR5_num ++;
                my $out = "$chr\t$source\t$_\t.\t$strand\t\.\tID=$mRNAID.utr5p$UTR5_num;Parent=$mRNAID;\n";
                my @g = split /\t/;
                $sort{$out} = $g[1];
                $sort_UTR{$out} = 1;
            }
            foreach (sort {$UTR_sort{$a} <=> $UTR_sort{$b}} @UTR3) {
                $UTR3_num ++;
                my $out = "$chr\t$source\t$_\t.\t$strand\t\.\tID=$mRNAID.utr3p$UTR3_num;Parent=$mRNAID;\n";
                my @g = split /\t/;
                $sort{$out} = $g[1];
                $sort_UTR{$out} = 4;
            }

            foreach (sort {$sort{$b} <=> $sort{$a} or $sort_UTR{$a} <=> $sort_UTR{$b}} keys %sort) {
                print;
            }
        }
        elsif ($strand eq '.') {
            foreach (sort { (split /\t/, $b)[0] <=> (split /\t/, $a)[0] } @exon) {
                $exon_num ++;
                my $out = "$chr\t$source\texon\t$_\t.\t$strand\t\.\tID=$mRNAID.exon$exon_num;Parent=$mRNAID;\n";
                my @g = split /\t/;
                $sort{$out} = $g[0];
                $sort_UTR{$out} = 2;
            }
            foreach (sort {$sort{$b} <=> $sort{$a} or $sort_UTR{$a} <=> $sort_UTR{$b}} keys %sort) {
                print;
            }
        }
        print "\n";
    }
}

# ------------------------------------------------------------------------
# 根据 CDS 与 exon 的坐标计算 UTR。
# 对每个 exon 片段，检查它与所有 CDS 片段的重叠情况：
#   - 完全不重叠 -> 整个exon都是UTR
#   - 部分重叠   -> 把exon里CDS前面和/或后面剩下的部分记为UTR
# ------------------------------------------------------------------------
sub get_UTR {
    my @cds    = @{ $_[0] };
    my @exon   = @{ $_[1] };
    my $strand = $_[2];

    my (@utr, %cds_pos);
    foreach (@cds) {
        my @g = split /\t/;
        $cds_pos{$g[0]} = 1;
        $cds_pos{$g[1]} = 1;
    }
    my @cds_pos = sort { $a <=> $b } keys %cds_pos;
    return () unless @cds_pos;
    my ($cds_min, $cds_max) = ($cds_pos[0], $cds_pos[-1]);

    foreach (@exon) {
        my ($start, $end) = split /\t/;
        my $overlap = 0;
        foreach my $c (@cds) {
            my ($cs, $ce) = split /\t/, $c;
            next unless $cs <= $end && $ce >= $start;
            $overlap = 1;
            push @utr, "$start\t" . ($cs - 1) if $start < $cs;   # 前端UTR
            push @utr, ($ce + 1) . "\t$end"   if $end   > $ce;   # 后端UTR
        }
        push @utr, "$start\t$end" unless $overlap;
    }

    my @out;
    if ($strand eq "+") {
        @utr = sort { (split /\t/, $a)[0] <=> (split /\t/, $b)[0] } @utr;
        foreach (@utr) {
            my @g = split /\t/;
            if    ( $g[1] <= $cds_min ) { push @out, "five_prime_UTR\t$_"; }
            elsif ( $g[0] >= $cds_max ) { push @out, "three_prime_UTR\t$_"; }
        }
    }
    elsif ($strand eq "-") {
        @utr = sort { (split /\t/, $b)[0] <=> (split /\t/, $a)[0] } @utr;
        foreach (@utr) {
            my @g = split /\t/;
            if    ( $g[0] >= $cds_max ) { push @out, "five_prime_UTR\t$_"; }
            elsif ( $g[1] <= $cds_min ) { push @out, "three_prime_UTR\t$_"; }
        }
    }

    return @out;
}

# 计算 intron (exon 之间的间隔)；$intron_len 是判定为intron所需的最小间隔长度。
sub get_intron {
    my @exon       = @{ $_[0] };
    my $mRNA_ID    = $_[1];
    my $intron_len = $_[2];
    @exon = sort { (split /\t/, $a)[0] <=> (split /\t/, $b)[0] } @exon;

    my @intron;
    my $first_exon = shift @exon;
    my ($last_start, $last_end) = split /\t/, $first_exon;
    foreach (@exon) {
        my ($start, $end) = split /\t/, $_;
        if ($start > $last_end + $intron_len) {
            my $intron_start = $last_end + 1;
            my $intron_stop  = $start - 1;
            push @intron, "$intron_start\t$intron_stop";
        }
        else {
            my $value = $start - $last_end - 1;
            print STDERR "Warning: an intron length (value is $value) of mRNA $mRNA_ID < $intron_len was detected:\n\tThe former CDS/Exon: $last_start - $last_end\n\tThe latter CDS/Exon: $start - $end\n";
        }
        ($last_start, $last_end) = ($start, $end);
    }

    return @intron;
}

# 合并有重叠/相邻(间隔<=0)的区间，输入输出都是 "start\tend" 字符串列表。
sub merge_intervals {
    my @intervals = sort { $a->[0] <=> $b->[0] }
                    map  { [ split /\t/, $_ ] } @_;
    my @merged;
    for my $iv (@intervals) {
        if ( @merged && $iv->[0] <= $merged[-1][1] + 1 ) {
            $merged[-1][1] = $iv->[1] if $iv->[1] > $merged[-1][1];
        }
        else {
            push @merged, [ @$iv ];
        }
    }
    return map { "$_->[0]\t$_->[1]" } @merged;
}
