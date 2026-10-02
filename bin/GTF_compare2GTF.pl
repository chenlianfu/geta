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
    $0 [options] in1.gtf in2.gtf > compare_result.txt

    本程序用于比较两个GTF文件中的基因模型，基于exon信息统计两者的重叠、包含关系。不做跨文件去冗余（不会删除任何一方的基因）。不需要基因组序列。

    程序运行原理：
    (1) 读取两个GTF文件中feature为exon的行，按gene_id / transcript_id组织基因模型。若一个基因有多个转录本（可变剪接），只选取exon总长度（重叠区域不重复计算）最长的转录本代表该基因，长度相同时选transcript_id字典序靠前者。
    (2) 两个文件各自独立进行文件内部去冗余：基因模型得分 = exon总长度 + intron加分（见--intron_score）；同一文件内、同一染色体同一条链上，exon重叠比例（重叠碱基数 / 较短基因模型的exon总长度）超过--overlap_coverage的基因模型互为冗余，只保留得分最高者。
    (3) 将两个文件各自去冗余后的基因集合做跨文件两两比较，对exon重叠比例超过--overlap_coverage的基因对，按覆盖比例分类：
        - overlap：部分重叠，双方都没有被对方基本完全覆盖；
        - A_contained_in_B：文件1的基因被文件2的基因包含（文件1基因的覆盖比例 >= --containment_ratio，文件2基因的覆盖比例 < --containment_ratio）；
        - B_contained_in_A：文件2的基因被文件1的基因包含（与上一条方向相反）；
        - identical：双方覆盖比例都 >= --containment_ratio，即exon区域几乎完全重合。
    (4) 输出：两个文件各自的原始基因数、去冗余后基因数；各自与对方无重叠的基因数（及百分比）；跨文件重叠基因对总数及四种关系的数量（及百分比）；每一对重叠基因的详细信息；两个文件各自独有基因的列表。文件内部去冗余删除基因的信息输出到标准错误。

    注意事项：
    (1) 只接受恰好两个GTF文件。同一文件内的gene_id必须唯一（同一gene_id被视为同一个基因）；两个文件之间的gene_id不要求相同或唯一。
    (2) exon行必须包含gene_id；缺少transcript_id时，以gene_id代替。

    --intron_score <float>    default: 0.3
    intron让基因模型得分提升的比例。第一个intron让得分增加 (exon总长度 * 该值)，之后每多一个intron，增加量按 (1 - 该值) 的等比数列递减。

    --overlap_coverage <float>    default: 0.30
    判定两个基因模型存在重叠的阈值：重叠碱基数 / 较短基因模型的exon总长度 > 该值。文件内部去冗余和跨文件比较均使用该值。

    --containment_ratio <float>    default: 0.95
    判定包含关系的阈值：某基因模型的重叠碱基数 / 其自身exon总长度 >= 该值，则认为它基本被对方包含。应不小于--overlap_coverage。

    --cpu <int>    default: 8
    文件内部去冗余步骤的并行进程数。设为1则串行运行。

    --tmp_dir <string>    default: 系统临时目录下自动创建并在结束后清理
    并行进程交换结果用的临时目录（仅在--cpu > 1且分区任务数 > 1时使用）。

    --help    显示英文用法并退出。
    --chinese_help    显示中文用法并退出。

USAGE
my $usage_english = &get_usage_english();
if (@ARGV == 0) { die $usage_english }

my ($intron_score, $overlap_coverage, $containment_ratio, $cpu, $tmp_dir, $help, $chinese_help);
GetOptions(
    "intron_score:f"      => \$intron_score,
    "overlap_coverage:f"  => \$overlap_coverage,
    "containment_ratio:f" => \$containment_ratio,
    "cpu:i"               => \$cpu,
    "tmp_dir:s"           => \$tmp_dir,
    "help"                => \$help,
    "chinese_help"        => \$chinese_help,
) or die $usage_english;
if ($chinese_help) { die $usage }
if ($help) { die $usage_english }

########### 解析参数 #################
die "Error: 本程序只接受恰好两个GTF文件作为输入（实际输入了" . scalar(@ARGV) . "个）\n" unless @ARGV == 2;
my @input_GTF_display = @ARGV;
my @input_GTF = map { abs_path($_) or die "Error: Can not find file $_\n" } @ARGV;

$intron_score      = 0.3  unless defined $intron_score;
$overlap_coverage  = 0.3  unless defined $overlap_coverage;
$containment_ratio = 0.95 unless defined $containment_ratio;
$cpu = 8 unless defined $cpu;
$cpu = 1 if $cpu < 1;

if ($containment_ratio < $overlap_coverage) {
    warn "Warning: --containment_ratio ($containment_ratio) is smaller than --overlap_coverage ($overlap_coverage), this is unusual and may make the 'overlap' category disappear.\n";
}
###############################

# 分别独立解析两个GTF文件（互相独立的命名空间）
my ($gene_exon1, $gene_score1, $gene_len1, $chr1, $strand1, $tx1, $raw_count1) = &process_GTF($input_GTF[0]);
my ($gene_exon2, $gene_score2, $gene_len2, $chr2, $strand2, $tx2, $raw_count2) = &process_GTF($input_GTF[1]);

##########################################################################################
# 第一步：两个文件分别独立进行文件内部去冗余
##########################################################################################
my (%partition1, %partition2);
foreach my $gene_ID (keys %$gene_exon1) {
    push @{$partition1{"$chr1->{$gene_ID}\t$strand1->{$gene_ID}"}}, $gene_ID;
}
foreach my $gene_ID (keys %$gene_exon2) {
    push @{$partition2{"$chr2->{$gene_ID}\t$strand2->{$gene_ID}"}}, $gene_ID;
}

my @jobs;
foreach my $key (keys %partition1) {
    push @jobs, ["file1\t$key", $partition1{$key}, $gene_exon1, $gene_score1, $gene_len1];
}
foreach my $key (keys %partition2) {
    push @jobs, ["file2\t$key", $partition2{$key}, $gene_exon2, $gene_score2, $gene_len2];
}
# 优先处理基因数量较多的分区，有利于负载均衡
@jobs = sort { scalar(@{$b->[1]}) <=> scalar(@{$a->[1]}) or $a->[0] cmp $b->[0] } @jobs;

my %results = &run_partitions(@jobs);

my (%deleted1, %deleted2);
foreach my $key (sort keys %results) {
    my $deleted = $results{$key}{"deleted"} || [];
    my $log     = $results{$key}{"log"} || [];
    if ($key =~ m/^file1\t/) { $deleted1{$_} = 1 foreach @$deleted; }
    else                     { $deleted2{$_} = 1 foreach @$deleted; }
    #print STDERR "[$key] $_\n" foreach @$log;
}

my @survivors1 = sort grep { !$deleted1{$_} } keys %$gene_exon1;
my @survivors2 = sort grep { !$deleted2{$_} } keys %$gene_exon2;

print STDERR "# Statistics of redundancy removal within each individual input GTF file:\n";
print STDERR "File $input_GTF_display[0]: $raw_count1 genes in total, " . scalar(@survivors1) . " genes remaining after removing redundancy within this file\n";
print STDERR "File $input_GTF_display[1]: $raw_count2 genes in total, " . scalar(@survivors2) . " genes remaining after removing redundancy within this file\n\n";

##########################################################################################
# 第二步：跨文件比较（file1的非冗余基因 vs file2的非冗余基因，只统计关系，不删除）
##########################################################################################
my %indexB;
foreach my $gene_ID (@survivors2) {
    my $key = "$chr2->{$gene_ID}\t$strand2->{$gene_ID}";
    foreach my $exon (split /\n/, $gene_exon2->{$gene_ID}) {
        my ($start, $end) = split /\t/, $exon;
        foreach my $index (int($start / 1000) .. int($end / 1000)) {
            $indexB{$key}{$index}{"$start\t$end"}{$gene_ID} = 1;
        }
    }
}

my (@pairs, %has_overlap1, %has_overlap2);
foreach my $gene_ID (@survivors1) {
    my $key = "$chr1->{$gene_ID}\t$strand1->{$gene_ID}";
    next unless exists $indexB{$key};

    my %candidate;
    foreach my $exon (split /\n/, $gene_exon1->{$gene_ID}) {
        my ($start, $end) = split /\t/, $exon;
        foreach my $index (int($start / 1000) .. int($end / 1000)) {
            next unless exists $indexB{$key}{$index};
            foreach my $region (keys %{$indexB{$key}{$index}}) {
                my ($r_start, $r_end) = split /\t/, $region;
                next unless ($r_end >= $start && $r_start <= $end);
                $candidate{$_} = 1 foreach keys %{$indexB{$key}{$index}{$region}};
            }
        }
    }

    my @exon_A = split /\n/, $gene_exon1->{$gene_ID};
    foreach my $other_gene (sort keys %candidate) {
        my @exon_B = split /\n/, $gene_exon2->{$other_gene};
        my $match_length = &get_match_length(\@exon_A, \@exon_B);
        next unless $match_length > 0;

        my $len_A = $gene_len1->{$gene_ID};
        my $len_B = $gene_len2->{$other_gene};
        my $ratio_A = $match_length / $len_A;
        my $ratio_B = $match_length / $len_B;
        my $ratio = $ratio_A > $ratio_B ? $ratio_A : $ratio_B;
        next unless $ratio > $overlap_coverage;

        my $relation;
        if    ($ratio_A >= $containment_ratio && $ratio_B >= $containment_ratio) { $relation = "identical"; }
        elsif ($ratio_A >= $containment_ratio) { $relation = "A_contained_in_B"; }
        elsif ($ratio_B >= $containment_ratio) { $relation = "B_contained_in_A"; }
        else                                   { $relation = "overlap"; }

        $has_overlap1{$gene_ID} = 1;
        $has_overlap2{$other_gene} = 1;
        push @pairs, {
            geneA => $gene_ID,    lenA => $len_A,
            geneB => $other_gene, lenB => $len_B,
            match_length => $match_length,
            ratioA => $ratio_A, ratioB => $ratio_B,
            relation => $relation,
        };
    }
}

my @unique1 = grep { !$has_overlap1{$_} } @survivors1;
my @unique2 = grep { !$has_overlap2{$_} } @survivors2;

my $n_total     = scalar @pairs;
my $n_overlap   = scalar grep { $_->{relation} eq "overlap" } @pairs;
my $n_identical = scalar grep { $_->{relation} eq "identical" } @pairs;
my $n_A_in_B    = scalar grep { $_->{relation} eq "A_contained_in_B" } @pairs;
my $n_B_in_A    = scalar grep { $_->{relation} eq "B_contained_in_A" } @pairs;

##########################################################################################
# 输出比较结果
##########################################################################################
print "# GTF Comparison Report\n";
print "# File1 (A): $input_GTF_display[0]\n";
print "# File2 (B): $input_GTF_display[1]\n";
print "# overlap_coverage threshold: $overlap_coverage\n";
print "# containment_ratio threshold: $containment_ratio\n";
print "#\n";
print "# File1 raw gene count: $raw_count1\tafter internal redundancy removal: " . scalar(@survivors1) . "\n";
print "# File2 raw gene count: $raw_count2\tafter internal redundancy removal: " . scalar(@survivors2) . "\n";
print "#\n";
print "# File1 genes with NO overlap in File2: " . scalar(@unique1) . &format_pct(scalar(@unique1), scalar(@survivors1)) . "\n";
print "# File2 genes with NO overlap in File1: " . scalar(@unique2) . &format_pct(scalar(@unique2), scalar(@survivors2)) . "\n";
print "#\n";
print "# Total cross-file overlapping gene-model pairs (ratio > $overlap_coverage): $n_total\n";
print "#   - overlap (partial overlap, neither side >= containment_ratio): $n_overlap" . &format_pct($n_overlap, $n_total) . "\n";
print "#   - A_contained_in_B (File1 gene contained in File2 gene, i.e. File2 包含 File1): $n_A_in_B" . &format_pct($n_A_in_B, $n_total) . "\n";
print "#   - B_contained_in_A (File2 gene contained in File1 gene, i.e. File1 包含 File2): $n_B_in_A" . &format_pct($n_B_in_A, $n_total) . "\n";
print "#   - identical (exon regions nearly identical / mutually contained): $n_identical" . &format_pct($n_identical, $n_total) . "\n";
print "#\n";
print "# ---- Pair details (sorted by match length, descending) ----\n";
print join("\t", qw/GeneA(File1) TranscriptA ChrA StrandA ExonLen_A GeneB(File2) TranscriptB ChrB StrandB ExonLen_B Match_length RatioA RatioB Relation/), "\n";
foreach my $pair (sort { $b->{match_length} <=> $a->{match_length} or $a->{geneA} cmp $b->{geneA} or $a->{geneB} cmp $b->{geneB} } @pairs) {
    my ($gA, $gB) = ($pair->{geneA}, $pair->{geneB});
    printf "%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%.4f\t%.4f\t%s\n",
        $gA, $tx1->{$gA}, $chr1->{$gA}, $strand1->{$gA}, $pair->{lenA},
        $gB, $tx2->{$gB}, $chr2->{$gB}, $strand2->{$gB}, $pair->{lenB},
        $pair->{match_length}, $pair->{ratioA}, $pair->{ratioB}, $pair->{relation};
}
print "#\n";
print "# ---- File1-only genes (no overlap with File2) ----\n";
print join("\t", qw/GeneID Transcript Chr Strand ExonLen/), "\n";
foreach my $g (@unique1) {
    print "$g\t$tx1->{$g}\t$chr1->{$g}\t$strand1->{$g}\t$gene_len1->{$g}\n";
}
print "#\n";
print "# ---- File2-only genes (no overlap with File1) ----\n";
print join("\t", qw/GeneID Transcript Chr Strand ExonLen/), "\n";
foreach my $g (@unique2) {
    print "$g\t$tx2->{$g}\t$chr2->{$g}\t$strand2->{$g}\t$gene_len2->{$g}\n";
}


sub format_pct {
    my ($num, $den) = @_;
    return "" unless $den;
    return sprintf(" (%.2f%%)", 100 * $num / $den);
}


# 解析单个GTF文件，仅使用exon行。每个基因选exon总长度最长的转录本作为代表。
# 返回：gene_exon, gene_score, gene_len, chr, strand, tx(代表转录本ID) 的哈希引用，以及基因总数
sub process_GTF {
    my ($file) = @_;
    my (%tx, $n_exon_lines, $n_no_gene, $n_no_tx);
    ($n_exon_lines, $n_no_gene, $n_no_tx) = (0, 0, 0);

    open my $in, '<', $file or die "Error: Can not open file $file, $!\n";
    while (<$in>) {
        chomp;
        s/\r$//;
        next if m/^\s*#/ || m/^\s*$/;
        my @f = split /\t/;
        next unless @f >= 9 && $f[2] eq "exon";
        $n_exon_lines++;
        my ($chr, $start, $end, $strand, $attr) = @f[0, 3, 4, 6, 8];
        ($start, $end) = ($end, $start) if $start > $end;

        my ($gid) = $attr =~ m/(?:^|[;\s])gene_id\s+"?([^";]+)"?/;
        my ($tid) = $attr =~ m/(?:^|[;\s])transcript_id\s+"?([^";]+)"?/;
        if (!defined $gid) { $n_no_gene++; next; }
        $gid =~ s/\s+$//;
        if (defined $tid) { $tid =~ s/\s+$//; }
        else { $tid = $gid; $n_no_tx++; }

        push @{$tx{$gid}{$tid}{"exons"}}, "$start\t$end";
        $tx{$gid}{$tid}{"chr"}    = $chr    unless defined $tx{$gid}{$tid}{"chr"};
        $tx{$gid}{$tid}{"strand"} = $strand unless defined $tx{$gid}{$tid}{"strand"};
    }
    close $in;

    warn "Warning: $n_no_gene exon line(s) in $file have no gene_id and were ignored.\n" if $n_no_gene;
    warn "Warning: $n_no_tx exon line(s) in $file have no transcript_id; gene_id was used instead.\n" if $n_no_tx;
    warn "Warning: no exon line found in $file.\n" unless $n_exon_lines;

    my (%gene_exon, %gene_score, %gene_len, %chr, %strand, %best_tx);
    foreach my $gid (sort keys %tx) {
        my ($best_tid, $best_len, @best_exons);
        foreach my $tid (sort keys %{$tx{$gid}}) {
            my @merged = &merge_regions($tx{$gid}{$tid}{"exons"});
            my $len = 0;
            foreach (@merged) { my ($s, $e) = split /\t/; $len += $e - $s + 1; }
            # 严格大于：长度相同时保留transcript_id字典序靠前者
            if (!defined $best_len || $len > $best_len) {
                ($best_tid, $best_len, @best_exons) = ($tid, $len, @merged);
            }
        }

        my $n_exon = scalar @best_exons;
        my $score = $best_len;
        if ($n_exon > 1) {
            my $add = $best_len * $intron_score;
            foreach (1 .. ($n_exon - 1)) {
                $score += $add;
                $add *= (1 - $intron_score);
            }
        }

        $gene_exon{$gid}  = join "\n", @best_exons;
        $gene_score{$gid} = $score;
        $gene_len{$gid}   = $best_len;
        $chr{$gid}        = $tx{$gid}{$best_tid}{"chr"};
        $strand{$gid}     = $tx{$gid}{$best_tid}{"strand"};
        $best_tx{$gid}    = $best_tid;
    }

    return (\%gene_exon, \%gene_score, \%gene_len, \%chr, \%strand, \%best_tx, scalar(keys %gene_exon));
}


# 对区间列表（"start\tend"）排序并合并相互重叠的区间（仅合并有碱基重叠的，不合并相邻的）
sub merge_regions {
    my @r = map { [split /\t/] } @{$_[0]};
    @r = sort { $a->[0] <=> $b->[0] or $a->[1] <=> $b->[1] } @r;
    my @out;
    foreach my $r (@r) {
        if (@out && $r->[0] <= $out[-1][1]) {
            $out[-1][1] = $r->[1] if $r->[1] > $out[-1][1];
        }
        else {
            push @out, [@$r];
        }
    }
    return map { "$_->[0]\t$_->[1]" } @out;
}


# 对一批分区任务分发执行文件内部去冗余，每个job为 [job_key, \@gene_ids, gene_exon_ref, gene_score_ref, gene_len_ref]
# 返回哈希：job_key => { deleted => \@deleted_gene_ids, log => \@log_lines }
sub run_partitions {
    my @jobs = @_;
    my %job_results;

    if ($cpu <= 1 || @jobs <= 1) {
        foreach my $job (@jobs) {
            my ($key, $ids_ref, $exon_ref, $score_ref, $len_ref) = @$job;
            my ($deleted, $log) = &find_redundant_in_partition($ids_ref, $exon_ref, $score_ref, $len_ref, $overlap_coverage);
            $job_results{$key} = { deleted => $deleted, log => $log };
        }
        return %job_results;
    }

    my $tmp_dir_created = 0;
    if ($tmp_dir) {
        $tmp_dir = abs_path($tmp_dir) if -e $tmp_dir;
        mkdir $tmp_dir unless -e $tmp_dir;
        $tmp_dir = abs_path($tmp_dir);
    }
    else {
        $tmp_dir = tempdir("GTF_compare_XXXXXX", TMPDIR => 1, CLEANUP => 1);
        $tmp_dir_created = 1;
    }

    my @queue = @jobs;
    my %running; # pid => [job_key, tmpfile]

    while (@queue || %running) {
        while (@queue && scalar(keys %running) < $cpu) {
            my $job = shift @queue;
            my ($key, $ids_ref, $exon_ref, $score_ref, $len_ref) = @$job;
            my ($fh, $tmpfile) = tempfile(DIR => $tmp_dir, SUFFIX => ".dat", UNLINK => 0);
            close $fh;

            my $pid = fork();
            die "Error: fork failed: $!" unless defined $pid;

            if ($pid == 0) {
                my ($deleted, $log) = &find_redundant_in_partition($ids_ref, $exon_ref, $score_ref, $len_ref, $overlap_coverage);
                eval { nstore({ deleted => $deleted, log => $log }, $tmpfile); };
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

        if ($exit_code != 0) {
            warn "Warning: worker process for job [$key] exited with code $exit_code; its results may be incomplete\n";
        }
        if (-e $tmpfile) {
            my $data = eval { retrieve($tmpfile) };
            $job_results{$key} = $data if $data;
            unlink $tmpfile;
        }
    }

    rmdir $tmp_dir if $tmp_dir_created && -d $tmp_dir;
    return %job_results;
}


# 对单个分区（同一文件、同一染色体、同一条链）做重叠检测与文件内部去冗余。
# 返回 (\@deleted, \@log)
sub find_redundant_in_partition {
    my ($gene_ids_ref, $gene_exon_ref, $gene_score_ref, $gene_len_ref, $this_overlap_coverage) = @_;
    my @gene_ids = @$gene_ids_ref;
    my (%exon2gene, %overlap);

    foreach my $gene_ID (@gene_ids) {
        foreach my $exon (split /\n/, $gene_exon_ref->{$gene_ID}) {
            my ($start, $end) = split /\t/, $exon;
            foreach my $index (int($start / 1000) .. int($end / 1000)) {
                $exon2gene{$index}{"$start\t$end"}{$gene_ID} = 1;
            }
        }
    }

    foreach my $gene_ID (@gene_ids) {
        foreach my $exon (split /\n/, $gene_exon_ref->{$gene_ID}) {
            my ($start, $end) = split /\t/, $exon;
            foreach my $index (int($start / 1000) .. int($end / 1000)) {
                foreach my $region (keys %{$exon2gene{$index}}) {
                    my ($r_start, $r_end) = split /\t/, $region;
                    next unless ($r_end >= $start && $r_start <= $end);
                    foreach my $other_gene (keys %{$exon2gene{$index}{$region}}) {
                        next if $other_gene eq $gene_ID;
                        $overlap{$gene_ID}{$other_gene} = 1;
                        $overlap{$other_gene}{$gene_ID} = 1;
                    }
                }
            }
        }
    }

    # 将相互重叠的基因归并成簇
    my %cluster;
    while (%overlap) {
        my %cluster_one;
        my @stack = ((sort keys %overlap)[0]);
        while (@stack) {
            my $one = shift @stack;
            next if $cluster_one{$one};
            $cluster_one{$one} = 1;
            if (exists $overlap{$one}) {
                push @stack, keys %{$overlap{$one}};
                delete $overlap{$one};
            }
        }
        $cluster{join("\t", sort keys %cluster_one)} = 1;
    }

    my (@deleted, @log);
    foreach my $cluster (sort keys %cluster) {
        my @genes = split /\t/, $cluster;
        @genes = sort { $gene_score_ref->{$b} <=> $gene_score_ref->{$a} or $a cmp $b } @genes;
        my %genes = map { $_ => 1 } @genes;

        while (@genes) {
            my $gene = shift @genes;
            next unless $genes{$gene};
            my @gene_exon_arr = split /\n/, $gene_exon_ref->{$gene};
            my $gene_len = $gene_len_ref->{$gene};
            delete $genes{$gene};

            foreach my $target_gene (sort { $gene_score_ref->{$b} <=> $gene_score_ref->{$a} or $a cmp $b } keys %genes) {
                my @target_exon_arr = split /\n/, $gene_exon_ref->{$target_gene};
                my $target_len = $gene_len_ref->{$target_gene};

                my $match_length = &get_match_length(\@gene_exon_arr, \@target_exon_arr);
                my $ratio1 = $match_length / $gene_len;
                my $ratio2 = $match_length / $target_len;
                my $ratio = $ratio1 > $ratio2 ? $ratio1 : $ratio2;

                if ($ratio > $this_overlap_coverage) {
                    push @log, "Delete gene $target_gene (Score: $gene_score_ref->{$target_gene}), for its exon coverage ratio with gene $gene (Score: $gene_score_ref->{$gene}) is: $ratio > $this_overlap_coverage";
                    push @deleted, $target_gene;
                    delete $genes{$target_gene};
                }
            }
            @genes = sort { $gene_score_ref->{$b} <=> $gene_score_ref->{$a} or $a cmp $b } keys %genes;
        }
    }

    return (\@deleted, \@log);
}


# 两组区间（各自内部互不重叠）的交集总长度
sub get_match_length {
    my @region1 = @{$_[0]};
    my @region2 = @{$_[1]};

    my @region_match;
    foreach my $region1 (@region1) {
        my ($start1, $end1) = split /\t/, $region1;
        foreach my $region2 (@region2) {
            my ($start2, $end2) = split /\t/, $region2;
            if ($start1 <= $end2 && $start2 <= $end1) {
                my ($start, $end) = ($start1, $end1);
                $start = $start2 if $start2 > $start1;
                $end = $end2 if $end2 < $end1;
                push @region_match, "$start\t$end";
            }
        }
    }
    return &get_region_length(\@region_match);
}


# 区间并集的总长度
sub get_region_length {
    my @merged = &merge_regions($_[0]);
    my $length = 0;
    foreach (@merged) {
        my ($s, $e) = split /\t/;
        $length += $e - $s + 1;
    }
    return $length;
}


sub get_usage_english {

my $usage = <<USAGE;
Usage:
    $0 [options] in1.gtf in2.gtf > compare_result.txt

    This program compares the gene models in two GTF files, based on exon information only, and reports their overlap / containment relationships. No genome sequence is needed. It does NOT perform cross-file redundancy removal (no gene from either file is deleted).

    How it works:
    (1) Only 'exon' lines are read, and grouped by gene_id / transcript_id. If a gene has several transcripts (alternative splicing), only the transcript with the longest total exon length (overlapping bases counted once) represents the gene; ties are broken by the alphabetically first transcript_id.
    (2) Each file is deduplicated independently. A gene model's score = total exon length + intron bonus (see --intron_score). Within one file, on the same chromosome and strand, gene models whose exon overlap ratio (overlapping bases / total exon length of the shorter model) exceeds --overlap_coverage are redundant, and only the highest-scoring one is kept.
    (3) The two non-redundant gene sets are compared across files. For each pair whose overlap ratio exceeds --overlap_coverage, the relationship is classified as:
        - overlap: partial overlap, neither gene is almost fully covered by the other;
        - A_contained_in_B: File1 gene is contained in File2 gene (File1 gene's coverage ratio >= --containment_ratio, File2 gene's < --containment_ratio);
        - B_contained_in_A: the reverse of the above;
        - identical: both coverage ratios >= --containment_ratio, i.e. exon regions are nearly identical.
    (4) Output: raw and post-deduplication gene counts of each file; the number (and percentage) of genes in each file with no overlap in the other file; the total number of overlapping pairs and the count (and percentage) of each relationship type; details of every pair; and the list of genes unique to each file. Genes deleted during the internal deduplication are reported to STDERR.

    Notes:
    (1) Exactly two GTF files are required. gene_id must identify one gene within a file; gene_id values need not match or be unique across the two files.
    (2) Exon lines must contain gene_id. If transcript_id is missing, gene_id is used instead.

    --intron_score <float>    default: 0.3
    Proportion by which introns raise a gene model's score. The first intron adds (total exon length * this value); each further intron adds an amount that shrinks by a factor of (1 - this value).

    --overlap_coverage <float>    default: 0.30
    Two gene models overlap when (overlapping bases) / (total exon length of the shorter model) > this value. Used both for within-file deduplication and for the cross-file comparison.

    --containment_ratio <float>    default: 0.95
    A gene model is considered contained in the other when (overlapping bases) / (its own total exon length) >= this value. Should not be smaller than --overlap_coverage.

    --cpu <int>    default: 8
    Number of parallel worker processes for the within-file deduplication. Set to 1 to run serially.

    --tmp_dir <string>    default: auto-created under the system temp directory and removed on exit
    Scratch directory for exchanging results between workers (only used when --cpu > 1 and there is more than one job).

    --help    display this help and exit.
    --chinese_help    display the Chinese usage and exit.

USAGE

return $usage;
}
