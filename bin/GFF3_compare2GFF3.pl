#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Cwd qw/abs_path/;
use File::Temp qw/tempfile tempdir/;
use Storable qw/nstore retrieve/;
use POSIX qw/:sys_wait_h/;

my $usage = <<USAGE;
Usage:
    $0 [options] genome.fasta in1.gff3 in2.gff3 > compare_result.txt

    本程序用于比较两个GFF3文件中的基因模型，得到两者的重叠、包含关系统计信息，不做跨文件的去冗余（即不会删除任何一方的基因）。

    程序运行原理：
    (1) 对两个输入的GFF3文件，分别独立进行文件内部的去冗余处理（算法与GFF3_merging_and_removing_redundancy一致：按CDS长度、intron数量、基因完整性打分，同一文件内相互重叠比例超过阈值的基因模型只保留得分最高的一个），得到每个文件各自的非冗余基因集合。
    (2) 将两个文件各自的非冗余基因集合进行两两比较（只比较跨文件的基因对，不再比较同一文件内部），对CDS有重叠（重叠碱基数 / 较小基因模型CDS长度 > --overlap_coverage）的基因对，记录其匹配碱基数、各自的覆盖比例，并根据覆盖比例判定关系类型：
        - overlap（部分重叠）：两个基因模型的CDS互相重叠，但任意一方的CDS都没有被对方基本完全覆盖；
        - B_contained_in_A（文件2基因被文件1基因包含）：文件2基因的CDS基本完全落在文件1基因的CDS范围内（比例 >= --containment_ratio），而文件1基因未被完全覆盖；
        - A_contained_in_B（文件1基因被文件2基因包含）：与上一条方向相反；
        - identical（近似完全一致）：两个基因模型的CDS互相之间的覆盖比例都 >= --containment_ratio，即两者CDS区域几乎完全重合。
    (3) 最终输出：两个文件各自的原始基因数、去除内部冗余后的基因数；各自与对方完全不重叠的基因数；跨文件重叠的基因对总数及上述四类关系各自的数量；每一对重叠基因的详细信息（基因ID、染色体、正负链、CDS长度、匹配碱基数、覆盖比例、关系类型）；以及两个文件中各自独有（与对方无重叠）的基因列表。

    程序使用须知：
    (1) 程序需要输入基因组序列，用于文件内部去冗余打分时判断基因模型的完整性。
    (2) 本程序接受输入带可变剪接的GFF3文件，对一个基因的多个转录本打分后，选其中得分最高转录本的CDS信息代表该基因，再用于去冗余和跨文件比较。
    (3) 本程序只接受恰好两个GFF3文件作为输入。两个文件内部的基因ID不要求跨文件唯一（因为不做跨文件的合并去冗余），但同一个文件内部的基因ID必须唯一。

    --intron_score <float>    default: 0.3
    设置intron让基因模型得分提升比例，用于文件内部去冗余打分。基因模型的第一个intron让得分增加比例 = 该参数值，第二个及其后的intron让得分增加比例按 (1 - 该参数值) 的等比数列递减，所有intron让得分增加比例之和不超过100%。

    --complete5p_score <float>    default: 0.5
    设置基因模型在5'端完整时得分增加比例，用于文件内部去冗余打分。

    --complete3p_score <float>    default: 0.5
    设置基因模型在3'端完整时得分增加比例，用于文件内部去冗余打分。

    --overlap_coverage <float>    default: 0.30
    设置判定"两个基因模型存在重叠关系"的覆盖度阈值。当两个基因模型CDS重叠区碱基数 / 较小基因模型CDS长度 > 该阈值时，判定为重叠（该阈值同时用于文件内部去冗余，以及跨文件的重叠关系判定）。

    --containment_ratio <float>    default: 0.95
    设置判定"包含"关系的覆盖度阈值。当一个基因模型CDS的重叠碱基数 / 其自身CDS长度 >= 该阈值时，判定该基因模型基本完全被对方包含。该值应不小于--overlap_coverage。

    --start_codon <string>    default: ATG
    设置起始密码子。若有多个起始密码子，则使用逗号分割。

    --stop_codon <string>    default: TAA,TAG,TGA
    设置终止密码子。若有多个终止密码子，则使用逗号分割。

    --cpu <int>    default: 8
    设置文件内部去冗余步骤的并行worker数量。设为1则完全串行运行，不fork子进程。（跨文件比较步骤本身不fork，始终串行执行）

    --tmp_dir <string>    default: 系统临时目录下自动创建，程序结束后自动清理
    并行worker之间交换计算结果所用的临时目录（仅在--cpu > 1且分区任务数 > 1时使用）。

    --help    default: None
    display this help and exit.

    --chinese_help    default: None
    使用该参数后，程序给出中文用法并退出。

USAGE
my $usage_english = &get_usage_english();
if (@ARGV==0){die $usage_english}

my ($intron_score, $complete5p_score, $complete3p_score, $overlap_coverage, $containment_ratio, $start_codon, $stop_codon, $cpu, $tmp_dir, $help, $chinese_help);
GetOptions(
    "intron_score:f" => \$intron_score,
    "complete5p_score:f" => \$complete5p_score,
    "complete3p_score:f" => \$complete3p_score,
    "overlap_coverage:f" => \$overlap_coverage,
    "containment_ratio:f" => \$containment_ratio,
    "start_codon:s" => \$start_codon,
    "stop_codon:s" => \$stop_codon,
    "cpu:i" => \$cpu,
    "tmp_dir:s" => \$tmp_dir,
    "help" => \$help,
    "chinese_help" => \$chinese_help,
);
if ( $chinese_help ) { die $usage }
if ( $help ) { die $usage_english }

########### 解析参数 #################
my $input_genome = abs_path(shift @ARGV);
my @input_GFF3_display = @ARGV;
die "Error: 本程序只接受恰好两个GFF3文件作为输入（实际输入了" . scalar(@ARGV) . "个）\n" unless @ARGV == 2;
my @input_GFF3 = map { abs_path($_) } @ARGV;

$intron_score = 0.3 unless defined $intron_score;
$complete5p_score = 0.5 unless defined $complete5p_score;
$complete3p_score = 0.5 unless defined $complete3p_score;
$overlap_coverage = 0.3 unless defined $overlap_coverage;
$containment_ratio = 0.95 unless defined $containment_ratio;
$start_codon ||= "ATG";
$stop_codon ||= "TAA,TAG,TGA";
$cpu = 8 unless defined $cpu;
$cpu = 1 if $cpu < 1;
my (%start_codon, %stop_codon);
foreach (split /,/, $start_codon) { $start_codon{$_} = 1; }
foreach (split /,/, $stop_codon) { $stop_codon{$_} = 1; }

if ( $containment_ratio < $overlap_coverage ) {
    warn "Warning: --containment_ratio ($containment_ratio) is smaller than --overlap_coverage ($overlap_coverage), this is unusual and may make the 'overlap' category disappear.\n";
}

###############################

# 读取基因组序列（用于文件内部去冗余打分时的完整性判断）
my %seq;
my %missing_chr_warned; # 记录已经警告过的、genome.fasta中找不到的序列ID，避免同一个ID反复刷屏警告
{
    my $seq_id;
    open IN, $input_genome or die "Error: Can not open file $input_genome, $!";
    while (<IN>) {
        chomp;
        if ( m/^>(\S+)/ ) { $seq_id = $1; }
        else { tr/atcgn/ATCGN/; $seq{$seq_id} .= $_; }
    }
    close IN;
}

# 分别独立解析两个GFF3文件，得到各自的基因信息（互相独立的命名空间，不要求跨文件基因ID唯一）
my ($gene_info1, $gene_CDS1, $gene_score1, $gene_CDS_length1, $chr1, $strand1, $raw_count1) = &process_GFF3($input_GFF3[0]);
my ($gene_info2, $gene_CDS2, $gene_score2, $gene_CDS_length2, $chr2, $strand2, $raw_count2) = &process_GFF3($input_GFF3[1]);

# 基因组序列此后不再需要，及时释放内存（尤其在即将fork子进程之前）
%seq = ();

# 因缺少有效CDS而被整体跳过（不参与后续任何比较）的基因数量
my $skipped_no_cds1 = $raw_count1 - scalar(keys %$gene_CDS1);
my $skipped_no_cds2 = $raw_count2 - scalar(keys %$gene_CDS2);
if ( $skipped_no_cds1 > 0 ) {
    print STDERR "Note: $skipped_no_cds1 gene(s) in $input_GFF3_display[0] have no mRNA with a valid CDS and were skipped entirely (see warnings above for details).\n";
}
if ( $skipped_no_cds2 > 0 ) {
    print STDERR "Note: $skipped_no_cds2 gene(s) in $input_GFF3_display[1] have no mRNA with a valid CDS and were skipped entirely (see warnings above for details).\n";
}

##########################################################################################
# 第一步：两个文件分别独立进行文件内部去冗余（只在各自文件内部检测重叠冲突）
##########################################################################################
my (%partition1, %partition2);
foreach my $gene_ID ( keys %$gene_CDS1 ) {
    my $key = "$chr1->{$gene_ID}\t$strand1->{$gene_ID}";
    push @{$partition1{$key}}, $gene_ID;
}
foreach my $gene_ID ( keys %$gene_CDS2 ) {
    my $key = "$chr2->{$gene_ID}\t$strand2->{$gene_ID}";
    push @{$partition2{$key}}, $gene_ID;
}

my @jobs;
foreach my $key ( keys %partition1 ) {
    push @jobs, ["file1\t$key", $partition1{$key}, $gene_CDS1, $gene_score1, $gene_CDS_length1];
}
foreach my $key ( keys %partition2 ) {
    push @jobs, ["file2\t$key", $partition2{$key}, $gene_CDS2, $gene_score2, $gene_CDS_length2];
}
# 优先处理基因数量较多的分区，有利于负载均衡
@jobs = sort { scalar(@{$b->[1]}) <=> scalar(@{$a->[1]}) } @jobs;

my %results = &run_partitions(@jobs);

my (%deleted1, %deleted2);
foreach my $key ( sort keys %results ) {
    my $deleted = $results{$key}{"deleted"} || [];
    my $log = $results{$key}{"log"} || [];
    if ( $key =~ m/^file1\t/ ) {
        foreach ( @$deleted ) { $deleted1{$_} = 1; }
    }
    else {
        foreach ( @$deleted ) { $deleted2{$_} = 1; }
    }
    foreach my $line ( @$log ) {
        print STDERR "[$key] $line\n";
    }
}

my @survivors1 = sort grep { ! $deleted1{$_} } keys %$gene_CDS1;
my @survivors2 = sort grep { ! $deleted2{$_} } keys %$gene_CDS2;

print STDERR "# Statistics of redundancy removal within each individual input GFF3 file:\n";
print STDERR "File $input_GFF3_display[0]: $raw_count1 genes in total, " . scalar(@survivors1) . " genes remaining after removing redundancy within this file\n";
print STDERR "File $input_GFF3_display[1]: $raw_count2 genes in total, " . scalar(@survivors2) . " genes remaining after removing redundancy within this file\n\n";

##########################################################################################
# 第二步：跨文件比较（只比较file1的非冗余基因 vs file2的非冗余基因，不做删除，只统计关系）
##########################################################################################
my %indexB;
foreach my $gene_ID ( @survivors2 ) {
    my $key = "$chr2->{$gene_ID}\t$strand2->{$gene_ID}";
    foreach my $CDS ( split /\n/, $gene_CDS2->{$gene_ID} ) {
        my ($start, $end) = (split /\t/, $CDS)[0,1];
        my $index1 = int($start / 1000);
        my $index2 = int($end / 1000);
        foreach my $index ( $index1 .. $index2 ) {
            $indexB{$key}{$index}{"$start\t$end"}{$gene_ID} = 1;
        }
    }
}

my (@pairs, %has_overlap1, %has_overlap2);
foreach my $gene_ID ( @survivors1 ) {
    my $key = "$chr1->{$gene_ID}\t$strand1->{$gene_ID}";
    next unless exists $indexB{$key};

    my %candidate;
    foreach my $CDS ( split /\n/, $gene_CDS1->{$gene_ID} ) {
        my ($start, $end) = (split /\t/, $CDS)[0,1];
        my $index1 = int($start / 1000);
        my $index2 = int($end / 1000);
        foreach my $index ( $index1 .. $index2 ) {
            next unless exists $indexB{$key}{$index};
            foreach my $region ( keys %{$indexB{$key}{$index}} ) {
                my ($r_start, $r_end) = split /\t/, $region;
                next unless ( $r_end >= $start && $r_start <= $end );
                foreach my $other_gene ( keys %{$indexB{$key}{$index}{$region}} ) {
                    $candidate{$other_gene} = 1;
                }
            }
        }
    }

    foreach my $other_gene ( sort keys %candidate ) {
        my @CDS_A = split /\n/, $gene_CDS1->{$gene_ID};
        my @CDS_B = split /\n/, $gene_CDS2->{$other_gene};
        my $match_length = &get_match_length(\@CDS_A, \@CDS_B);
        next unless $match_length > 0;

        my $len_A = $gene_CDS_length1->{$gene_ID};
        my $len_B = $gene_CDS_length2->{$other_gene};
        my $ratio_A = $match_length / $len_A;
        my $ratio_B = $match_length / $len_B;
        my $ratio = $ratio_A > $ratio_B ? $ratio_A : $ratio_B;

        next unless $ratio > $overlap_coverage;

        my $relation;
        if ( $ratio_A >= $containment_ratio && $ratio_B >= $containment_ratio ) {
            $relation = "identical";          # 两者CDS区域几乎完全重合（互相包含）
        }
        elsif ( $ratio_A >= $containment_ratio ) {
            $relation = "A_contained_in_B";   # file1基因被file2基因包含
        }
        elsif ( $ratio_B >= $containment_ratio ) {
            $relation = "B_contained_in_A";   # file2基因被file1基因包含
        }
        else {
            $relation = "overlap";            # 部分重叠，互不包含
        }

        $has_overlap1{$gene_ID} = 1;
        $has_overlap2{$other_gene} = 1;
        push @pairs, {
            geneA => $gene_ID, lenA => $len_A,
            geneB => $other_gene, lenB => $len_B,
            match_length => $match_length,
            ratioA => $ratio_A, ratioB => $ratio_B,
            relation => $relation,
        };
    }
}

my @unique1 = grep { ! $has_overlap1{$_} } @survivors1;
my @unique2 = grep { ! $has_overlap2{$_} } @survivors2;

my $n_overlap   = scalar grep { $_->{relation} eq "overlap" } @pairs;
my $n_identical = scalar grep { $_->{relation} eq "identical" } @pairs;
my $n_A_in_B    = scalar grep { $_->{relation} eq "A_contained_in_B" } @pairs;
my $n_B_in_A    = scalar grep { $_->{relation} eq "B_contained_in_A" } @pairs;

##########################################################################################
# 输出比较结果
##########################################################################################
print "# GFF3 Comparison Report\n";
print "# File1 (A): $input_GFF3_display[0]\n";
print "# File2 (B): $input_GFF3_display[1]\n";
print "# overlap_coverage threshold: $overlap_coverage\n";
print "# containment_ratio threshold: $containment_ratio\n";
print "#\n";
print "# File1 raw gene count: $raw_count1\tafter internal redundancy removal: " . scalar(@survivors1) . "\n";
print "# File2 raw gene count: $raw_count2\tafter internal redundancy removal: " . scalar(@survivors2) . "\n";
print "#\n";
print "# File1 genes with NO overlap in File2: " . scalar(@unique1) . &format_pct(scalar(@unique1), scalar(@survivors1)) . "\n";
print "# File2 genes with NO overlap in File1: " . scalar(@unique2) . &format_pct(scalar(@unique2), scalar(@survivors2)) . "\n";
print "#\n";
print "# Total cross-file overlapping gene-model pairs (ratio > $overlap_coverage): " . scalar(@pairs) . "\n";
print "#   - overlap (partial overlap, neither side >= containment_ratio): $n_overlap" . &format_pct($n_overlap, scalar(@pairs)) . "\n";
print "#   - A_contained_in_B (File1 gene contained in File2 gene, i.e. File2 包含 File1): $n_A_in_B" . &format_pct($n_A_in_B, scalar(@pairs)) . "\n";
print "#   - B_contained_in_A (File2 gene contained in File1 gene, i.e. File1 包含 File2): $n_B_in_A" . &format_pct($n_B_in_A, scalar(@pairs)) . "\n";
print "#   - identical (CDS regions nearly identical / mutually contained): $n_identical" . &format_pct($n_identical, scalar(@pairs)) . "\n";
print "#\n";
print "# ---- Pair details (sorted by match length, descending) ----\n";
print join("\t", qw/GeneA(File1) ChrA StrandA CDS_len_A GeneB(File2) ChrB StrandB CDS_len_B Match_length RatioA RatioB Relation/), "\n";
foreach my $pair ( sort { $b->{match_length} <=> $a->{match_length} } @pairs ) {
    printf "%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%.4f\t%.4f\t%s\n",
        $pair->{geneA}, $chr1->{$pair->{geneA}}, $strand1->{$pair->{geneA}}, $pair->{lenA},
        $pair->{geneB}, $chr2->{$pair->{geneB}}, $strand2->{$pair->{geneB}}, $pair->{lenB},
        $pair->{match_length}, $pair->{ratioA}, $pair->{ratioB}, $pair->{relation};
}
print "#\n";
print "# ---- File1-only genes (no overlap with File2) ----\n";
print join("\t", qw/GeneID Chr Strand CDS_length/), "\n";
foreach my $gene_ID ( @unique1 ) {
    print "$gene_ID\t$chr1->{$gene_ID}\t$strand1->{$gene_ID}\t$gene_CDS_length1->{$gene_ID}\n";
}
print "#\n";
print "# ---- File2-only genes (no overlap with File1) ----\n";
print join("\t", qw/GeneID Chr Strand CDS_length/), "\n";
foreach my $gene_ID ( @unique2 ) {
    print "$gene_ID\t$chr2->{$gene_ID}\t$strand2->{$gene_ID}\t$gene_CDS_length2->{$gene_ID}\n";
}


# 生成形如 " (12.34%)" 的百分比字符串，分母为0时返回空字符串（避免除以0报错，同时避免输出无意义的百分比）
sub format_pct {
    my ($numerator, $denominator) = @_;
    return "" unless $denominator;
    return sprintf(" (%.2f%%)", 100 * $numerator / $denominator);
}


# 解析单个GFF3文件，对每个基因的所有转录本打分，选出得分最高转录本的CDS信息代表该基因。
# 返回：gene_info, gene_CDS, gene_score, gene_CDS_length, chr, strand （均为哈希引用，键均为该文件自己的基因ID）
sub process_GFF3 {
    my ($file) = @_;
    my %GFF3_info = &get_geneModels_from_GFF3($file);
    my (%gene_CDS, %gene_score, %gene_CDS_length, %chr, %strand, %gene_info);

    foreach my $gene_ID ( sort keys %GFF3_info ) {
        my @gene_header = split /\t/, $GFF3_info{$gene_ID}{"header"};
        my $this_gene_chr = $gene_header[0];
        my $this_gene_strand = $gene_header[6];

        my @mRNA_ID = @{$GFF3_info{$gene_ID}{"mRNA_ID"} || []};
        my (%mRNA_score, %mRNA_CDS, %mRNA_CDS_length);
        foreach my $mRNA_ID ( sort @mRNA_ID ) {
            my $mRNA_info = $GFF3_info{$gene_ID}{"mRNA_info"}{$mRNA_ID};
            my $mRNA_header = $GFF3_info{$gene_ID}{"mRNA_header"}{$mRNA_ID};
            my @mRNA_header = split /\t/, $mRNA_header;
            my ($this_chr, $this_strand) = ($mRNA_header[0], $mRNA_header[6]);

            my (@CDS, $CDS_length);
            $CDS_length = 0;
            foreach ( split /\n/, $mRNA_info || "" ) {
                my @field = split /\t/;
                if ( $field[2] eq "CDS" ) {
                    push @CDS, "$field[3]\t$field[4]\t$field[6]";
                    $CDS_length += (abs($field[4] - $field[3]) + 1);
                }
            }

            # 有些mRNA记录没有任何CDS子特征（如非编码转录本、或GFF3记录不完整），无法参与打分和CDS比较，跳过并警告，而不是让$CDS_length/$score以undef参与后续数值运算
            if ( @CDS == 0 ) {
                warn "Warning: mRNA $mRNA_ID (gene $gene_ID) in $file has no CDS feature, skipped.\n";
                next;
            }

            $mRNA_CDS{$mRNA_ID} = join "\n", @CDS;
            $mRNA_CDS_length{$mRNA_ID} = $CDS_length;

            my $integrity = &analysis_geneModels_integrity(\@CDS, $this_chr, $this_strand);

            my $score = $CDS_length;
            if ( @CDS == 2 ) {
                $score += $CDS_length * $intron_score;
            }
            elsif ( @CDS > 2 ) {
                my $add_score = $CDS_length * $intron_score;
                $score += $add_score;
                foreach ( 1 .. (@CDS - 2) ) {
                    $add_score = $add_score * ( 1 - $intron_score );
                    $score += $add_score;
                }
            }
            if ( $integrity eq "5prime_partial" ) {
                $score += $CDS_length * $complete3p_score;
            }
            elsif ( $integrity eq "3prime_partial" ) {
                $score += $CDS_length * $complete5p_score;
            }
            elsif ( $integrity eq "complete" ) {
                $score += $CDS_length * ($complete3p_score + $complete5p_score);
            }

            $mRNA_score{$mRNA_ID} = $score;
        }

        # 若该基因所有mRNA都没有CDS（即上面全部被跳过），此基因无法参与CDS比较，整体跳过并警告
        unless ( %mRNA_score ) {
            warn "Warning: gene $gene_ID in $file has no mRNA with a valid CDS, this gene is skipped entirely.\n";
            next;
        }

        my @scored_mRNA_ID = sort {$mRNA_score{$b} <=> $mRNA_score{$a} or $a cmp $b} keys %mRNA_score;
        $chr{$gene_ID} = $this_gene_chr;
        $strand{$gene_ID} = $this_gene_strand;
        $gene_info{$gene_ID} = $GFF3_info{$gene_ID};
        $gene_CDS{$gene_ID} = $mRNA_CDS{$scored_mRNA_ID[0]};
        $gene_score{$gene_ID} = sprintf("%d", $mRNA_score{$scored_mRNA_ID[0]});
        $gene_CDS_length{$gene_ID} = $mRNA_CDS_length{$scored_mRNA_ID[0]};
    }

    my $total_raw_gene_count = scalar keys %GFF3_info;
    return (\%gene_info, \%gene_CDS, \%gene_score, \%gene_CDS_length, \%chr, \%strand, $total_raw_gene_count);
}


# 对一批分区任务(jobs)分发执行文件内部去冗余，每个job为 [job_key, \@gene_ids, gene_CDS_ref, gene_score_ref, gene_CDS_length_ref]。
# cpu<=1或任务数<=1时直接串行执行；否则按$cpu并行fork子进程，通过Storable临时文件回传结果，逻辑与原合并程序的并行调度一致。
# 返回一个哈希：job_key => { deleted => \@deleted_gene_ids, log => \@log_lines }
sub run_partitions {
    my @jobs = @_;
    my %job_results;

    if ( $cpu <= 1 || @jobs <= 1 ) {
        foreach my $job ( @jobs ) {
            my ($key, $ids_ref, $cds_ref, $score_ref, $len_ref) = @$job;
            my ($deleted, $log) = &find_redundant_in_partition($ids_ref, $cds_ref, $score_ref, $len_ref, $overlap_coverage);
            $job_results{$key} = { deleted => $deleted, log => $log };
        }
        return %job_results;
    }

    my $tmp_dir_created = 0;
    if ( $tmp_dir ) {
        $tmp_dir = abs_path($tmp_dir);
        mkdir $tmp_dir unless -e $tmp_dir;
    }
    else {
        $tmp_dir = tempdir( "GFF3_compare_XXXXXX", TMPDIR => 1, CLEANUP => 1 );
        $tmp_dir_created = 1;
    }

    my @queue = @jobs;
    my %running; # pid => [job_key, tmpfile]

    while ( @queue || %running ) {
        while ( @queue && scalar(keys %running) < $cpu ) {
            my $job = shift @queue;
            my ($key, $ids_ref, $cds_ref, $score_ref, $len_ref) = @$job;
            my ($fh, $tmpfile) = tempfile( DIR => $tmp_dir, SUFFIX => ".dat", UNLINK => 0 );
            close $fh;

            my $pid = fork();
            die "Error: fork failed: $!" unless defined $pid;

            if ( $pid == 0 ) {
                my ($deleted, $log) = &find_redundant_in_partition($ids_ref, $cds_ref, $score_ref, $len_ref, $overlap_coverage);
                eval { nstore( { deleted => $deleted, log => $log }, $tmpfile ); };
                if ($@) {
                    warn "Warning: worker for job [$key] failed to write result to $tmpfile: $@";
                    exit 1;
                }
                exit 0;
            }
            else {
                $running{$pid} = [$key, $tmpfile];
            }
        }

        last unless %running;
        my $finished_pid = waitpid(-1, 0);
        next if $finished_pid <= 0;

        my $exit_code = $? >> 8;
        my $info = delete $running{$finished_pid};
        next unless $info;
        my ($key, $tmpfile) = @$info;

        if ( $exit_code != 0 ) {
            warn "Warning: worker process for job [$key] exited with code $exit_code; its results may be incomplete\n";
        }

        if ( -e $tmpfile ) {
            my $data = eval { retrieve($tmpfile) };
            $job_results{$key} = $data if $data;
            unlink $tmpfile;
        }
    }

    rmdir $tmp_dir if $tmp_dir_created && -d $tmp_dir;
    return %job_results;
}


# 对单个分区（同一文件内、同一条染色体、同一条链上的所有基因）做重叠检测与文件内部去冗余。
# 纯函数：只通过传入的哈希引用读取数据，不修改任何数据、不直接打印，把"应删除的基因ID列表"和"日志文本"作为返回值交给调用者处理。
# 与原GFF3_merging_and_removing_redundancy中的同名函数算法完全一致，仅将原来的全局变量改为显式传参，以支持同时处理两个互相独立的文件。
sub find_redundant_in_partition {
    my ($gene_ids_ref, $gene_CDS_ref, $gene_score_ref, $gene_CDS_length_ref, $this_overlap_coverage) = @_;
    my @gene_ids = @$gene_ids_ref;
    my (%cds2gene, %overlap);

    foreach my $gene_ID ( @gene_ids ) {
        foreach my $CDS ( split /\n/, $gene_CDS_ref->{$gene_ID} ) {
            my ($start, $end) = (split /\t/, $CDS)[0,1];
            my $index1 = int($start / 1000);
            my $index2 = int($end / 1000);
            foreach my $index ( $index1 .. $index2 ) {
                $cds2gene{$index}{"$start\t$end"}{$gene_ID} = 1;
            }
        }
    }

    foreach my $gene_ID ( @gene_ids ) {
        foreach my $CDS ( split /\n/, $gene_CDS_ref->{$gene_ID} ) {
            my ($start, $end) = (split /\t/, $CDS)[0,1];
            my $index1 = int($start / 1000);
            my $index2 = int($end / 1000);
            foreach my $index ( $index1 .. $index2 ) {
                foreach my $region ( keys %{$cds2gene{$index}} ) {
                    my ($r_start, $r_end) = split /\t/, $region;
                    next unless ( $r_end >= $start && $r_start <= $end );
                    foreach my $other_gene ( keys %{$cds2gene{$index}{$region}} ) {
                        next if $other_gene eq $gene_ID;
                        $overlap{$gene_ID}{$other_gene} = 1;
                        $overlap{$other_gene}{$gene_ID} = 1;
                    }
                }
            }
        }
    }

    my %cluster;
    while (%overlap) {
        my %cluster_one;
        my @stack = ( (keys %overlap)[0] );
        while (@stack) {
            my $one = shift @stack;
            next if $cluster_one{$one};
            $cluster_one{$one} = 1;
            if ( exists $overlap{$one} ) {
                push @stack, keys %{$overlap{$one}};
                delete $overlap{$one};
            }
        }
        $cluster{ join("\t", sort keys %cluster_one) } = 1;
    }

    my (@deleted, @log);
    foreach my $cluster ( sort keys %cluster ) {
        my @genes = split /\t/, $cluster;
        @genes = sort { $gene_score_ref->{$b} <=> $gene_score_ref->{$a} or $a cmp $b } @genes;
        my %genes = map { $_ => 1 } @genes;

        while (@genes) {
            my $gene = shift @genes;
            next unless $genes{$gene};
            my @gene_CDS_arr = split /\n/, $gene_CDS_ref->{$gene};
            my $gene_len = $gene_CDS_length_ref->{$gene};
            delete $genes{$gene};

            foreach my $target_gene ( sort { $gene_score_ref->{$b} <=> $gene_score_ref->{$a} } keys %genes ) {
                my @target_CDS_arr = split /\n/, $gene_CDS_ref->{$target_gene};
                my $target_len = $gene_CDS_length_ref->{$target_gene};

                my $match_length = &get_match_length(\@gene_CDS_arr, \@target_CDS_arr);
                my $ratio1 = $match_length / $gene_len;
                my $ratio2 = $match_length / $target_len;
                my $ratio = $ratio1 > $ratio2 ? $ratio1 : $ratio2;

                if ( $ratio > $this_overlap_coverage ) {
                    push @log, "Delete gene $target_gene (Score: $gene_score_ref->{$target_gene}), for its CDS coverage ratio with gene $gene (Score: $gene_score_ref->{$gene}) is: $ratio > $this_overlap_coverage";
                    push @deleted, $target_gene;
                    delete $genes{$target_gene};
                }
            }
            @genes = sort { $gene_score_ref->{$b} <=> $gene_score_ref->{$a} or $a cmp $b } keys %genes;
        }
    }

    return (\@deleted, \@log);
}


sub get_match_length {
    my @region1 = @{$_[0]};
    my @region2 = @{$_[1]};

    my $out_length;
    my @region_match;
    foreach my $region1 (@region1) {
        my ($start1, $end1) = split /\t/, $region1;
        foreach my $region2 (@region2) {
            my ($start2, $end2) = split /\t/, $region2;
            if ($start1 < $end2 && $start2 < $end1) {
                my ($start, $end) = ($start1, $end1);
                $start = $start2 if $start2 > $start1;
                $end = $end2 if $end2 < $end1;
                push @region_match, "$start\t$end";
            }
        }
    }

    $out_length = &get_gene_CDS_length(\@region_match);
    return $out_length;
}

sub get_gene_CDS_length {
    my @region = @{$_[0]};
    return 0 unless @region;
    @region = sort { (split /\t/, $a)[0] <=> (split /\t/, $b)[0] } @region;

    my $out_length;
    my $last_region = shift @region;
    my @last = split /\t/, $last_region;
    $out_length += ($last[1] - $last[0] + 1);
    foreach ( @region ) {
        my @last_region = split /\t/, $last_region;
        my @region = split /\t/;

        if ($region[0] > $last_region[1]) {
            $out_length += ($region[1] - $region[0] + 1);
            $last_region = $_;
        }
        elsif ($region[1] > $last_region[1]) {
            $out_length += ($region[1] - $last_region[1]);
            $last_region = $_;
        }
        else {
            next;
        }
    }

    return $out_length;
}


sub analysis_geneModels_integrity {
    my @CDS = @{$_[0]};
    my ($chr, $strand) = ($_[1], $_[2]);
    my $out;

    # genome.fasta中找不到该序列ID时，无法截取序列判断起止密码子，按"internal"处理（不给完整性加分），并只警告一次
    unless ( exists $seq{$chr} ) {
        unless ( $missing_chr_warned{$chr} ) {
            warn "Warning: sequence ID '$chr' was not found in the input genome fasta; gene models on it will get no completeness (5'/3' complete) score bonus.\n";
            $missing_chr_warned{$chr} = 1;
        }
        return "internal";
    }

    @CDS = sort { (split /\t/, $a)[0] <=> (split /\t/, $b)[0] } @CDS;
    my $seq = "";
    foreach ( @CDS ) {
        my @field = split /\t/;
        $seq .= substr($seq{$chr}, $field[0] - 1, $field[1] - $field[0] + 1);
    }
    if ( $strand eq "-" ) {
        $seq = reverse $seq;
        $seq =~ tr/ATCGatcgn/TAGCTAGCN/;
    }

    if ( $seq =~ m/^(\w{3})/ && exists $start_codon{$1} ) {
        if ( $seq =~ m/(\w{3})$/ && exists $stop_codon{$1} ) {
            $out = "complete";
        }
        else {
            $out = "3prime_partial";
        }
    }
    elsif ( $seq =~ m/(\w{3})$/ && exists $stop_codon{$1} ) {
        $out = "5prime_partial";
    }
    else {
        $out = "internal";
    }

    return $out;
}


sub get_geneModels_from_GFF3 {
    my %gene_info;
    my $input_file = $_[0];
    # 第一轮，找gene信息
    open IN, $input_file or die "Error: Can not open file $input_file, $!";
    while (<IN>) {
        if ( m/\tgene\t.*ID=([^;\s]+)/ ) {
            $gene_info{$1}{"header"} = $_;
        }
    }
    close IN;
    # 第二轮，找Parent值是geneID的信息，包含但不限于 mRNA 信息
    my %mRNA_ID2gene_ID;
    open IN, $input_file or die "Error: Can not open file $input_file, $!";
    while (<IN>) {
        if ( m/Parent=([^;\s]+)/ ) {
            my $parent = $1;
            if ( exists $gene_info{$parent} ) {
                if ( m/ID=([^;\s]+)/ ) {
                    push @{$gene_info{$parent}{"mRNA_ID"}}, $1;
                    $gene_info{$parent}{"mRNA_header"}{$1} = $_;
                    $mRNA_ID2gene_ID{$1} = $parent;
                }
            }
        }
    }
    close IN;
    # 第三轮，找Parent值不是geneID的信息
    open IN, $input_file or die "Error: Can not open file $input_file, $!";
    while (<IN>) {
        if ( m/Parent=([^;\s]+)/ && exists $mRNA_ID2gene_ID{$1} ) {
            my $parent = $1;
            $gene_info{$mRNA_ID2gene_ID{$1}}{"mRNA_info"}{$parent} .= $_;
        }
    }
    close IN;

    return %gene_info;
}


sub get_usage_english {

my $usage = <<USAGE;
Usage:
    $0 [options] genome.fasta in1.gff3 in2.gff3 > compare_result.txt

    This program compares the gene models in two GFF3 files and reports their overlap / containment relationships. It does NOT perform any cross-file redundancy removal (no gene from either file is deleted).

    How it works:
    (1) The two input GFF3 files are each independently deduplicated internally, using the same scoring/redundancy algorithm as GFF3_merging_and_removing_redundancy (CDS length, intron count, and completeness are used to score gene models; within a single file, when two gene models' CDS overlap ratio exceeds the threshold, only the higher-scoring one is kept).
    (2) The two files' resulting non-redundant gene sets are then compared pairwise, across files only (not within a file). For any pair whose CDS overlap ratio (overlap length / the shorter gene model's CDS length) exceeds --overlap_coverage, the program records the matched base pairs and each gene's own coverage ratio, and classifies the relationship as:
        - overlap: the two gene models' CDS overlap, but neither one is almost fully covered by the other;
        - B_contained_in_A: File2's gene model is almost entirely contained within File1's gene model's CDS (ratio >= --containment_ratio), while File1's gene is not fully covered;
        - A_contained_in_B: the reverse direction of the above;
        - identical: both genes' coverage ratios are >= --containment_ratio, i.e. their CDS regions are nearly identical.
    (3) Final output includes: each file's raw gene count and count after internal redundancy removal; the number of genes in each file that have no overlap with the other file; the total number of cross-file overlapping gene pairs and the count of each relationship type above; detailed per-pair information (gene IDs, chromosome, strand, CDS lengths, matched length, coverage ratios, relationship type); and the list of genes unique to each file (no overlap with the other file).

    Usage instructions:
    (1) The program needs a genome sequence as input, used for gene model completeness checks during internal-redundancy scoring.
    (2) This program accepts GFF3 files with alternative splicing; for each gene, all its transcripts are scored and the CDS information of the highest-scoring transcript represents that gene, for both internal deduplication and cross-file comparison.
    (3) This program only accepts exactly two GFF3 files as input. Gene IDs are not required to be unique across the two files (since no cross-file merging/deduplication is performed), but gene IDs must be unique within each individual file.

    --intron_score <float>    default: 0.3
    Sets the proportion by which introns increase a gene model's score, used for internal-redundancy scoring. See GFF3_merging_and_removing_redundancy for the detailed formula.

    --complete5p_score <float>    default: 0.5
    Sets the score increase proportion for a complete 5' end, used for internal-redundancy scoring.

    --complete3p_score <float>    default: 0.5
    Sets the score increase proportion for a complete 3' end, used for internal-redundancy scoring.

    --overlap_coverage <float>    default: 0.30
    Sets the coverage threshold for judging that two gene models "overlap". When (overlapping CDS base pairs) / (the shorter gene model's CDS length) > this threshold, the two are judged to overlap. This threshold is used both for internal (within-file) redundancy removal and for cross-file overlap judgement.

    --containment_ratio <float>    default: 0.95
    Sets the coverage threshold for judging a "containment" relationship. When a gene model's own coverage ratio (matched length / its own CDS length) >= this threshold, that gene model is judged to be almost entirely contained within the other. Should not be smaller than --overlap_coverage.

    --start_codon <string>    default: ATG
    Set the start codon. If there are multiple start codons, separate them with commas.

    --stop_codon <string>    default: TAA,TAG,TGA
    Set the stop codon. If there are multiple stop codons, separate them with commas.

    --cpu <int>    default: 8
    Number of parallel worker processes used for the internal (within-file) redundancy removal step. Set to 1 to run fully serially with no forking. (The cross-file comparison step itself is always run serially, without forking.)

    --tmp_dir <string>    default: auto-created under the system temp directory, and cleaned up automatically on exit
    Scratch directory used by worker processes to exchange results (only used when --cpu > 1 and there is more than one job).

    --help    default: None
    display this help and exit.

    --chinese_help    default: None
    display the chinese usage and exit.

USAGE

return $usage;
}
