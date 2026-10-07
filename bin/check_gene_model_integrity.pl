#!/usr/bin/env perl
#===============================================================================
# check_gene_model_integrity.pl
# 功能：检测 GFF3 中基因模型（CDS）的完整性，并在标准输出中给出统计结果
#       检测类型：5'端缺失、3'端缺失、双端缺失、内部终止密码子、移码
#       每个问题转录本的详细报告通过 --out 参数输出到文件
#===============================================================================
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(max min);

########################################################################
# 中文用法说明（放在代码最前面；英文用法说明放在文件最末尾的__DATA__部分）
########################################################################
my $usage_cn = <<'USAGE';
Usage:
    perl check_gene_model_integrity.pl [options] genome.fa genes.gff3

    本程序用于检测GFF3文件中基因模型（CDS）的完整性，在标准输出中给出统计结果，并可通过--out参数将每个存在问题的转录本的详细报告输出到文件。检测的问题类型包括：5'端缺失、3'端缺失、双端缺失、内部终止密码子、移码。

    判定规则：

    (1) 程序按Parent将CDS归入各个转录本，按转录方向拼接CDS序列，并按第一个CDS的phase去掉序列开头多余的碱基。

    (2) 5'端缺失：第一个CDS的phase不为0，或者首密码子不是起始密码子。

    (3) 3'端缺失：末密码子不是终止密码子（且CDS之外紧邻的3 bp也不是终止密码子），或者CDS长度（去除phase后）不是3的倍数。若终止密码子位于CDS之外紧邻的3 bp，仍判定为3'端完整，并在统计中单独给出此类转录本的数量。

    (4) 双端缺失：同时满足5'端缺失和3'端缺失。

    (5) 内部终止密码子：除最后一个完整密码子之外，其它位置出现终止密码子。

    (6) 移码：相邻CDS的phase与按长度推算的值不一致，或者CDS长度（去除phase后）不是3的倍数。

    (7) 一个转录本可以同时属于多种类型；基因层面按"任一转录本有该问题"计数。

    使用须知：

    (1) 程序需要输入基因组序列和GFF3文件，两者的顺序若写反，程序会根据文件后缀自动交换；输入文件支持.gz压缩格式。

    (2) 起始密码子和终止密码子由--genetic_code参数指定的遗传密码自动确定。

    (3) CDS所在的序列不在基因组FASTA中，或者同一转录本的CDS分布于不同序列或不同链的转录本，将被跳过，并在标准错误中给出警告。

    (4) 标准输出为统计结果（包含转录本层面和基因层面）；若设置了--out参数，则每个问题转录本的详细报告输出到该文件，文件为制表符分隔的文本，第一行为以#开头的表头；若不设置--out参数，则不输出详细报告。

    --genetic_code <int>    default: 1
    设置遗传密码，程序据此自动确定起始密码子和终止密码子，用于判断基因模型的完整性以及检测内部终止密码子。该参数对应的值请参考NCBI Genetic Codes: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi 。支持的遗传密码编号为：1、2、3、4、5、6、9、10、11、12、13、14、15、16、21、22、23、24、25、26、27、28、29、30、31、32、33。例如：1为标准遗传密码（起始密码子为ATG，终止密码子为TAA、TAG、TGA），11为细菌/古菌/植物质体遗传密码，2为脊椎动物线粒体遗传密码，4为支原体/螺旋体等遗传密码（TGA编码色氨酸）。

    --out <string>    default: 无（不设置时不输出详细报告）
    设置一个文件路径，程序会将每个存在问题的转录本的详细报告输出到该文件中。报告各列依次为：基因ID、转录本ID、序列ID、链方向、CDS起始坐标、CDS终止坐标、CDS数量、去除phase后的CDS长度、末端状态（5'缺失/3'缺失/双端缺失/两端完整）、首密码子、末密码子、内部终止密码子数量、内部终止密码子位置（氨基酸序号）、是否移码、异常原因说明。

    --help    default: None
    Display the English usage and exit.

    --chinese_help    default: None
    显示中文用法说明并退出。

USAGE

use constant PROG => 'check_gene_model_integrity.pl';

########################################################################
# 读取英文用法说明（位于文件最末尾的__DATA__部分）
########################################################################
my $usage_en = do { local $/; <DATA> };

if ( @ARGV == 0 ) { print STDERR $usage_en; exit 1; }

my ($genetic_code, $out_file, $start_codon_opt, $stop_codon_opt, $help, $chinese_help);
GetOptions(
    "genetic_code:i" => \$genetic_code,
    "out:s"          => \$out_file,
    # 为与GFF3_merging_and_removing_redundancy保持一致而保留，不在用法说明中给出：若设置，则优先于--genetic_code所确定的密码子。
    "start_codon:s"  => \$start_codon_opt,
    "stop_codon:s"   => \$stop_codon_opt,
    "help"           => \$help,
    "chinese_help"   => \$chinese_help,
) or die $usage_en;
if ( $chinese_help ) { print $usage_cn; exit 0; }
if ( $help )         { print $usage_en; exit 0; }
die "Error: genome.fa and genes.gff3 are required.\n\n$usage_en" unless @ARGV == 2;

my ($fasta_file, $gff_file) = @ARGV;
# 若顺序写反，自动交换
if ($fasta_file =~ /\.gff3?(\.gz)?$/i && $gff_file !~ /\.gff3?(\.gz)?$/i) {
    ($fasta_file, $gff_file) = ($gff_file, $fasta_file);
}

#---------------------------- 遗传密码表 / 起始、终止密码子 -------------------
$genetic_code = 1 unless ( defined $genetic_code && $genetic_code > 0 );
my $table = $genetic_code;
my ($code_ref, $start_codon_ref, $stop_codon_ref) = codon_table($genetic_code);
my (%start_codon, %stop_codon);
if ( defined $start_codon_opt && $start_codon_opt ne "" ) {
    foreach ( split /,/, $start_codon_opt ) {
        s/\s+//g; $_ = uc($_); tr/U/T/;
        $start_codon{$_} = 1 if $_ ne "";
    }
}
else {
    %start_codon = %$start_codon_ref;
}
if ( defined $stop_codon_opt && $stop_codon_opt ne "" ) {
    foreach ( split /,/, $stop_codon_opt ) {
        s/\s+//g; $_ = uc($_); tr/U/T/;
        $stop_codon{$_} = 1 if $_ ne "";
    }
}
else {
    %stop_codon = %$stop_codon_ref;
}
die "Error: no start codon is available.\n" unless %start_codon;
die "Error: no stop codon is available.\n"  unless %stop_codon;

# 提前打开 --out 文件，路径有误时尽早报错
my $out_fh;
if ( defined $out_file && $out_file ne "" ) {
    open $out_fh, '>', $out_file or die "Error: Can not create file $out_file, $!\n";
}

#---------------------------- 读取基因组 --------------------------------------
my %genome;
{
    my $fh = open_in($fasta_file);
    my $id;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        if ($line =~ /^>(\S+)/) {
            $id = $1;
            $genome{$id} = '';
        } elsif (defined $id) {
            $genome{$id} .= uc $line;
        }
    }
    close $fh;
}

#---------------------------- 解析 GFF3 ---------------------------------------
my (%type_of, %parents_of, %cds_of);
{
    my $fh = open_in($gff_file);
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        last if $line =~ /^##FASTA/;
        next if $line =~ /^#/ || $line =~ /^\s*$/;
        my @f = split /\t/, $line;
        next if @f < 9;
        my ($seqid, undef, $type, $start, $end, undef, $strand, $phase, $attr) = @f;
        my %a = parse_attributes($attr);

        if (defined $a{ID}) {
            $type_of{$a{ID}} = $type;
            $parents_of{$a{ID}} = [split /,/, $a{Parent}] if defined $a{Parent};
        }
        if ($type eq 'CDS' && defined $a{Parent}) {
            for my $p (split /,/, $a{Parent}) {
                push @{ $cds_of{$p} }, {
                    seqid  => $seqid,
                    start  => $start,
                    end    => $end,
                    strand => $strand,
                    phase  => ($phase =~ /^[012]$/ ? $phase : undef),
                };
            }
        }
    }
    close $fh;
}

#---------------------------- 逐转录本检测 ------------------------------------
my ($total_tx, $skipped_noseq, $skipped_incons) = (0, 0, 0);
my %tx_cnt;          # 转录本层面各类计数
my %gene_cat;        # gene => { 类别 => 1 }
my %all_genes;       # 含 CDS 的基因
my %gene_has_perfect;
my @detail;
my $stop_outside = 0; # 终止密码子位于 CDS 之外的转录本数

for my $tid (sort keys %cds_of) {
    my @segs   = @{ $cds_of{$tid} };
    my $seqid  = $segs[0]{seqid};
    my $strand = $segs[0]{strand} eq '-' ? '-' : '+';

    if (grep { $_->{seqid} ne $seqid
               || ($_->{strand} eq '-' ? '-' : '+') ne $strand } @segs) {
        warn "警告: 转录本 $tid 的 CDS 分布于不同序列或链，已跳过\n";
        $skipped_incons++;
        next;
    }
    unless (exists $genome{$seqid}) {
        warn "警告: 转录本 $tid 所在序列 $seqid 不在 FASTA 中，已跳过\n";
        $skipped_noseq++;
        next;
    }

    $total_tx++;
    my $gene = find_gene($tid);
    $all_genes{$gene} = 1;

    @segs = sort { $a->{start} <=> $b->{start} } @segs;
    my @tx = $strand eq '-' ? reverse @segs : @segs;   # 转录方向上的 CDS 顺序
    my $cds_min = $segs[0]{start};
    my $cds_max = max(map { $_->{end} } @segs);

    # 拼接 CDS
    my $seq = '';
    for my $s (@segs) {
        $seq .= substr($genome{$seqid}, $s->{start} - 1, $s->{end} - $s->{start} + 1);
    }
    $seq = revcomp($seq) if $strand eq '-';

    # phase 处理
    my $has_phase = !grep { !defined $_->{phase} } @segs;
    my $p0 = $has_phase ? $tx[0]{phase} : 0;

    # 相邻 CDS phase 一致性（移码线索）
    my $phase_bad = 0;
    if ($has_phase) {
        my $L = 0;
        for my $k (0 .. $#tx) {
            my $exp = ($p0 - $L) % 3;
            $phase_bad++ if $k > 0 && $tx[$k]{phase} != $exp;
            $L += $tx[$k]{end} - $tx[$k]{start} + 1;
        }
    }

    $seq = (length($seq) > $p0) ? substr($seq, $p0) : '' if $p0 > 0;

    my $ncodon = int(length($seq) / 3);
    my $remain = length($seq) % 3;

    # 内部终止（由遗传密码确定的终止密码子）
    my @internal;
    for my $i (0 .. $ncodon - 1) {
        next unless exists $stop_codon{ substr($seq, $i * 3, 3) };
        push @internal, $i + 1 if $i < $ncodon - 1 || $remain > 0;
    }

    # 首、末密码子
    my $first_codon = $ncodon > 0 ? substr($seq, 0, 3) : '-';
    my $last_codon  = ($remain == 0 && $ncodon > 0) ? substr($seq, ($ncodon - 1) * 3, 3) : '-';

    # 5' 端
    my $five = 0;
    if ($p0 > 0 || $ncodon < 1 || !exists $start_codon{ substr($seq, 0, 3) }) {
        $five = 1;
    }

    # 3' 端
    my $three = 1;
    my $note3 = '';
    if ($remain == 0 && $ncodon > 0) {
        if (exists $stop_codon{$last_codon}) {
            $three = 0;
        } else {
            # 终止密码子可能未包含在 CDS 内，检查 CDS 之外紧邻的 3 bp
            my $glen = length $genome{$seqid};
            my $d = '';
            if ($strand eq '+') {
                my $e = $cds_max;
                $d = substr($genome{$seqid}, $e, 3) if $e + 3 <= $glen;
            } else {
                my $s = $cds_min;
                $d = revcomp(substr($genome{$seqid}, $s - 4, 3)) if $s - 1 >= 3;
            }
            if (length($d) == 3 && exists $stop_codon{$d}) {
                $three = 0;
                $stop_outside++;
                $note3 = "终止密码子($d)位于CDS之外";
            }
        }
    }

    # 移码
    my $fs_len = ($remain != 0) ? 1 : 0;
    my $fs     = ($phase_bad > 0 || $fs_len) ? 1 : 0;

    # 分类
    my $terminal = ($five && $three) ? "双端缺失"
                 : $five             ? "5'缺失"
                 : $three            ? "3'缺失"
                 :                     "两端完整";
    my $perfect = (!$five && !$three && !@internal && !$fs) ? 1 : 0;

    $tx_cnt{five}++       if $five;
    $tx_cnt{three}++      if $three;
    $tx_cnt{both}++       if $five && $three;
    $tx_cnt{only5}++      if $five && !$three;
    $tx_cnt{only3}++      if $three && !$five;
    $tx_cnt{termok}++     if !$five && !$three;
    $tx_cnt{internal}++   if @internal;
    $tx_cnt{fs}++         if $fs;
    $tx_cnt{fs_phase}++   if $phase_bad > 0;
    $tx_cnt{fs_len}++     if $fs_len;
    $tx_cnt{perfect}++    if $perfect;
    $tx_cnt{bad}++        if !$perfect;

    my $gc = $gene_cat{$gene} ||= {};
    $gc->{five}     = 1 if $five;
    $gc->{three}    = 1 if $three;
    $gc->{both}     = 1 if $five && $three;
    $gc->{only5}    = 1 if $five && !$three;
    $gc->{only3}    = 1 if $three && !$five;
    $gc->{internal} = 1 if @internal;
    $gc->{fs}       = 1 if $fs;
    $gene_has_perfect{$gene} = 1 if $perfect;

    unless ($perfect) {
        my @notes;
        push @notes, "ncodon<1,CDS过短"                          if $ncodon < 1;
        push @notes, "首CDS phase=$p0"                           if $p0 > 0;
        push @notes, "首密码子($first_codon)不是起始密码子"       if $p0 == 0 && $ncodon >= 1 && !exists $start_codon{$first_codon};
        push @notes, "末密码子($last_codon)不是终止密码子"        if $three && $last_codon ne '-' && !$note3;
        push @notes, $note3                                      if $note3;
        push @notes, "phase不连续($phase_bad处)"                 if $phase_bad;
        push @notes, "长度非3倍数(余$remain)"                    if $fs_len;
        push @notes, "含" . scalar(@internal) . "个内部终止密码子" if @internal;
        push @detail, join("\t", $gene, $tid, $seqid, $strand,
            $cds_min, $cds_max, scalar(@segs),
            length($seq), $terminal,
            $first_codon, $last_codon,
            scalar(@internal), (@internal ? join(',', @internal) : '-'),
            ($fs ? '是' : '否'), (@notes ? join(';', @notes) : '-'));
    }
}

#---------------------------- 输出详细报告（--out） ---------------------------
if ($out_fh) {
    print $out_fh join("\t", '#Gene', qw(Transcript Seqid Strand CDS_start CDS_end CDS_num
                                          CDS_len_after_phase Terminal_status First_codon Last_codon
                                          N_internal_stops Stop_positions_aa Frameshift Note)), "\n";
    print $out_fh "$_\n" for @detail;
    close $out_fh;
}

#---------------------------- 输出统计结果 ------------------------------------
my $n_all = scalar keys %all_genes;
my $n_perfect_gene = scalar grep { $gene_has_perfect{$_} } keys %all_genes;
my %gc_n;
for my $g (keys %gene_cat) {
    $gc_n{$_}++ for keys %{ $gene_cat{$g} };
}
my $T = $total_tx;

print "================ 基因模型完整性检测结果 ================\n";
printf "基因组文件                     : %s\n", $fasta_file;
printf "GFF3 文件                      : %s\n", $gff_file;
printf "遗传密码表                     : %d\n", $table;
printf "起始密码子                     : %s\n", join(',', sort keys %start_codon);
printf "终止密码子                     : %s\n", join(',', sort keys %stop_codon);
printf "异常转录本明细文件             : %s\n", $out_fh ? $out_file : '未设置(--out)';
printf "含 CDS 的基因总数              : %d\n", $n_all;
printf "检测的转录本总数               : %d\n", $T;
printf "跳过的转录本(缺少序列)         : %d\n", $skipped_noseq;
printf "跳过的转录本(CDS链/序列不一致) : %d\n", $skipped_incons;
print  "\n---------- 转录本层面（类型可重叠）----------\n";
printf "完整且无异常的转录本           : %d (%.2f%%)\n", $tx_cnt{perfect} // 0, pct($tx_cnt{perfect}, $T);
printf "存在任一异常的转录本           : %d (%.2f%%)\n", $tx_cnt{bad} // 0,     pct($tx_cnt{bad}, $T);
printf "  5'端缺失(含双端)             : %d (%.2f%%)\n", $tx_cnt{five} // 0,     pct($tx_cnt{five}, $T);
printf "  3'端缺失(含双端)             : %d (%.2f%%)\n", $tx_cnt{three} // 0,    pct($tx_cnt{three}, $T);
printf "  仅5'端缺失                   : %d (%.2f%%)\n", $tx_cnt{only5} // 0,    pct($tx_cnt{only5}, $T);
printf "  仅3'端缺失                   : %d (%.2f%%)\n", $tx_cnt{only3} // 0,    pct($tx_cnt{only3}, $T);
printf "  双端缺失                     : %d (%.2f%%)\n", $tx_cnt{both} // 0,     pct($tx_cnt{both}, $T);
printf "  两端完整                     : %d (%.2f%%)\n", $tx_cnt{termok} // 0,   pct($tx_cnt{termok}, $T);
printf "  含内部终止密码子             : %d (%.2f%%)\n", $tx_cnt{internal} // 0, pct($tx_cnt{internal}, $T);
printf "  移码(合计)                   : %d (%.2f%%)\n", $tx_cnt{fs} // 0,       pct($tx_cnt{fs}, $T);
printf "    其中 phase 不连续          : %d\n", $tx_cnt{fs_phase} // 0;
printf "    其中 长度非3的倍数         : %d\n", $tx_cnt{fs_len} // 0;
printf "  (终止密码子位于CDS之外而判为3'完整: %d)\n", $stop_outside;
print  "\n---------- 基因层面（任一转录本有该问题即计入）----------\n";
printf "有完整且无异常转录本的基因     : %d (%.2f%%)\n", $n_perfect_gene, pct($n_perfect_gene, $n_all);
printf "无任何完整转录本的基因         : %d (%.2f%%)\n", $n_all - $n_perfect_gene, pct($n_all - $n_perfect_gene, $n_all);
printf "  5'端缺失(含双端)             : %d (%.2f%%)\n", $gc_n{five} // 0,     pct($gc_n{five}, $n_all);
printf "  3'端缺失(含双端)             : %d (%.2f%%)\n", $gc_n{three} // 0,    pct($gc_n{three}, $n_all);
printf "  仅5'端缺失                   : %d (%.2f%%)\n", $gc_n{only5} // 0,    pct($gc_n{only5}, $n_all);
printf "  仅3'端缺失                   : %d (%.2f%%)\n", $gc_n{only3} // 0,    pct($gc_n{only3}, $n_all);
printf "  双端缺失                     : %d (%.2f%%)\n", $gc_n{both} // 0,     pct($gc_n{both}, $n_all);
printf "  含内部终止密码子             : %d (%.2f%%)\n", $gc_n{internal} // 0, pct($gc_n{internal}, $n_all);
printf "  移码                         : %d (%.2f%%)\n", $gc_n{fs} // 0,       pct($gc_n{fs}, $n_all);
print  "========================================================\n";

#============================== 子程序 ========================================
sub open_in {
    my ($f) = @_;
    my $fh;
    if ($f =~ /\.gz$/i) {
        open $fh, '-|', 'gzip', '-dc', $f or die "无法打开 $f: $!\n";
    } else {
        open $fh, '<', $f or die "无法打开 $f: $!\n";
    }
    return $fh;
}

sub pct {
    my ($a, $b) = @_;
    $a //= 0;
    return $b ? $a / $b * 100 : 0;
}

sub parse_attributes {
    my ($str) = @_;
    my %h;
    for my $kv (split /;/, $str) {
        next unless $kv =~ /=/;
        my ($k, $v) = split /=/, $kv, 2;
        $v =~ s/%([0-9A-Fa-f]{2})/chr(hex($1))/ge;
        $h{$k} = $v;
    }
    return %h;
}

# 沿 Parent 链向上查找 gene；找不到则以转录本自身作为基因
sub find_gene {
    my ($id) = @_;
    my %seen;
    my $cur = $id;
    while (defined $cur && !$seen{$cur}++) {
        return $cur if defined $type_of{$cur} && $type_of{$cur} eq 'gene';
        my $p = $parents_of{$cur};
        last unless $p && @$p;
        $cur = $p->[0];
    }
    return $id;
}

sub revcomp {
    my ($s) = @_;
    $s = reverse $s;
    $s =~ tr/ACGTacgtRYKMBDHVrykmbdhv/TGCAtgcaYRMKVHDByrmkvhdb/;
    return $s;
}

# 子程序，根据遗传密码编号返回（密码子表、起始密码子表、终止密码子表）的哈希引用。
# 遗传密码的设置参考NCBI Genetic Codes，代码来自 GFF3_merging_and_removing_redundancy。
# 密码子表中，终止密码子对应的值为 "X"。
sub codon_table {
    my ($genetic_code) = @_;
    my %code = (
        "TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
        "TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
        "TAT" => "Y", "TAC" => "Y", "TAA" => "X", "TAG" => "X",
        "TGT" => "C", "TGC" => "C", "TGA" => "X", "TGG" => "W",
        "CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
        "CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
        "CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
        "CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
        "ATT" => "I", "ATC" => "I", "ATA" => "I", "ATG" => "M",
        "ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
        "AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
        "AGT" => "S", "AGC" => "S", "AGA" => "R", "AGG" => "R",
        "GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
        "GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
        "GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
        "GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G",
    );
    my %start_table;
    $start_table{"ATG"} = 1;
    if ( $genetic_code == 1 ) {
        # The Standard Code
    }
    elsif ( $genetic_code == 2 ) {
        # The Vertebrate Mitochondrial Code
        $code{"AGA"} = "X"; $code{"AGG"} = "X"; $code{"ATA"} = "M"; $code{"TGA"} = "W";
        $start_table{$_} = 1 foreach qw/ATA ATT ATC GTG/;
    }
    elsif ( $genetic_code == 3 ) {
        # The Yeast Mitochondrial Code
        $code{"ATA"} = "M"; $code{"CTT"} = "T"; $code{"CTC"} = "T"; $code{"CTA"} = "T"; $code{"CTG"} = "T"; $code{"TGA"} = "W";
        $start_table{$_} = 1 foreach qw/ATA GTG/;
    }
    elsif ( $genetic_code == 4 ) {
        # The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
        $code{"TGA"} = "W";
        $start_table{$_} = 1 foreach qw/ATA ATT ATC GTG CTG TTA TTG/;
    }
    elsif ( $genetic_code == 5 ) {
        # The Invertebrate Mitochondrial Code
        $code{"AGA"} = "S"; $code{"AGG"} = "S"; $code{"ATA"} = "M"; $code{"TGA"} = "W";
        $start_table{$_} = 1 foreach qw/ATA ATT ATC GTG TTG/;
    }
    elsif ( $genetic_code == 6 ) {
        # The Ciliate, Dasycladacean and Hexamita Nuclear Code
        $code{"TAA"} = "Q"; $code{"TAG"} = "Q";
    }
    elsif ( $genetic_code == 9 ) {
        # The Echinoderm and Flatworm Mitochondrial Code
        $code{"AAA"} = "N"; $code{"AGA"} = "S"; $code{"AGG"} = "S"; $code{"TGA"} = "W";
        $start_table{"GTG"} = 1;
    }
    elsif ( $genetic_code == 10 ) {
        # The Euplotid Nuclear Code
        $code{"TGA"} = "C";
    }
    elsif ( $genetic_code == 11 ) {
        # The Bacterial, Archaeal and Plant Plastid Code
        $start_table{$_} = 1 foreach qw/ATA ATT ATC GTG CTG TTG/;
    }
    elsif ( $genetic_code == 12 ) {
        # The Alternative Yeast Nuclear Code
        $code{"CTG"} = "S";
        $start_table{"CTG"} = 1;
    }
    elsif ( $genetic_code == 13 ) {
        # The Ascidian Mitochondrial Code
        $code{"AGA"} = "G"; $code{"AGG"} = "G"; $code{"ATA"} = "M"; $code{"TGA"} = "W";
        $start_table{$_} = 1 foreach qw/ATA GTG TTG/;
    }
    elsif ( $genetic_code == 14 ) {
        # The Alternative Flatworm Mitochondrial Code
        $code{"AAA"} = "N"; $code{"AGA"} = "S"; $code{"AGG"} = "S"; $code{"TAA"} = "Y"; $code{"TGA"} = "W";
    }
    elsif ( $genetic_code == 15 ) {
        # Blepharisma Nuclear Code
        $code{"TAG"} = "Q";
    }
    elsif ( $genetic_code == 16 ) {
        # Chlorophycean Mitochondrial Code
        $code{"TAG"} = "L";
    }
    elsif ( $genetic_code == 21 ) {
        # Trematode Mitochondrial Code
        $code{"TGA"} = "W"; $code{"ATA"} = "M"; $code{"AGA"} = "S"; $code{"AGG"} = "S"; $code{"AAA"} = "N";
        $start_table{"GTG"} = 1;
    }
    elsif ( $genetic_code == 22 ) {
        # Scenedesmus obliquus Mitochondrial Code
        $code{"TCA"} = "X"; $code{"TAG"} = "L";
    }
    elsif ( $genetic_code == 23 ) {
        # Thraustochytrium Mitochondrial Code
        $code{"TTA"} = "X";
        $start_table{$_} = 1 foreach qw/ATT GTG/;
    }
    elsif ( $genetic_code == 24 ) {
        # Rhabdopleuridae Mitochondrial Code
        $code{"AGA"} = "S"; $code{"AGG"} = "K"; $code{"TGA"} = "W";
        $start_table{$_} = 1 foreach qw/GTG CTG TTG/;
    }
    elsif ( $genetic_code == 25 ) {
        # Candidate Division SR1 and Gracilibacteria Code
        $code{"TGA"} = "G";
        $start_table{$_} = 1 foreach qw/GTG TTG/;
    }
    elsif ( $genetic_code == 26 ) {
        # Pachysolen tannophilus Nuclear Code
        $code{"CTG"} = "A";
        $start_table{"CTG"} = 1;
    }
    elsif ( $genetic_code == 27 ) {
        # Karyorelict Nuclear Code
        # TAA和TAG编码Gln；TGA既可编码Trp也可作为终止密码子，这里仍将TGA视为终止密码子。
        $code{"TAG"} = "Q"; $code{"TAA"} = "Q";
    }
    elsif ( $genetic_code == 28 ) {
        # Condylostoma Nuclear Code
        # TAA、TAG编码Gln，TGA编码Trp，但三者均可根据上下文作为终止密码子，因此这里仍将TAA、TAG、TGA都视为终止密码子。
    }
    elsif ( $genetic_code == 29 ) {
        # Mesodinium Nuclear Code
        $code{"TAA"} = "Y"; $code{"TAG"} = "Y";
    }
    elsif ( $genetic_code == 30 ) {
        # Peritrich Nuclear Code
        $code{"TAA"} = "E"; $code{"TAG"} = "E";
    }
    elsif ( $genetic_code == 31 ) {
        # Blastocrithidia Nuclear Code
    }
    elsif ( $genetic_code == 32 ) {
        # Balanophoraceae Plastid Code
        $code{"TAG"} = "W";
        $start_table{$_} = 1 foreach qw/ATA ATT ATC GTG CTG TTG/;
    }
    elsif ( $genetic_code == 33 ) {
        # Cephalodiscidae Mitochondrial UAA-Tyr Code
        $code{"TAA"} = "Y"; $code{"TGA"} = "W"; $code{"AGA"} = "S"; $code{"AGG"} = "K";
        $start_table{$_} = 1 foreach qw/GTG CTG TTG/;
    }
    else {
        print STDERR "Warning: 不支持的遗传密码 $genetic_code ，起始和终止密码子按标准遗传密码设置。\n";
        $table = 1;
    }

    my %stop_table;
    foreach ( keys %code ) {
        $stop_table{$_} = 1 if $code{$_} eq "X";
    }

    return (\%code, \%start_table, \%stop_table);
}

########################################################################
# 英文用法说明（放在文件最末尾，运行时通过DATA文件句柄读取）
########################################################################
__DATA__
Usage:
    perl check_gene_model_integrity.pl [options] genome.fa genes.gff3

    This program checks the integrity of the gene models (CDS) in a GFF3 file, prints summary statistics to STDOUT, and can write a detailed report of every problematic transcript to a file with the --out option. The problem types checked are: 5' end partial, 3' end partial, both ends partial, internal stop codons, and frameshifts.

    Rules:

    (1) The program assigns CDS features to transcripts by Parent, concatenates the CDS sequence in the transcription direction, and removes the extra bases at the beginning according to the phase of the first CDS.

    (2) 5' partial: the phase of the first CDS is not 0, or the first codon is not a start codon.

    (3) 3' partial: the last codon is not a stop codon (and the 3 bp immediately downstream of the CDS is not a stop codon either), or the CDS length (after removing the phase) is not a multiple of 3. If the stop codon lies in the 3 bp immediately outside the CDS, the 3' end is still regarded as complete, and the number of such transcripts is reported separately in the statistics.

    (4) Both ends partial: both the 5' partial and the 3' partial conditions are met.

    (5) Internal stop codon: a stop codon appears at any position other than the last complete codon.

    (6) Frameshift: the phases of adjacent CDS features are inconsistent with the values deduced from their lengths, or the CDS length (after removing the phase) is not a multiple of 3.

    (7) A transcript can belong to several types at the same time; at the gene level a gene is counted if any of its transcripts has the problem.

    Usage instructions:

    (1) The genome sequence and the GFF3 file are both required; if their order is reversed the program swaps them automatically according to the file extension; gzip-compressed (.gz) input files are supported.

    (2) The start codons and stop codons are determined automatically from the genetic code given by --genetic_code.

    (3) Transcripts whose sequence is absent from the genome FASTA, or whose CDS features lie on different sequences or strands, are skipped, and a warning is printed to STDERR.

    (4) STDOUT contains the summary statistics (transcript level and gene level); if --out is set, the detailed report of every problematic transcript is written to that file as tab-separated text whose first line is a header starting with #; if --out is not set, no detailed report is written.

    --genetic_code <int>    default: 1
    Set the genetic code, from which the start codons and stop codons are determined automatically; they are used to judge gene integrity and to detect internal stop codons. For the corresponding values, please refer to NCBI Genetic Codes: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi . The supported genetic code numbers are: 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33. For example: 1 is the standard code (start codon ATG; stop codons TAA, TAG, TGA), 11 is the bacterial/archaeal/plant plastid code, 2 is the vertebrate mitochondrial code, and 4 is the Mycoplasma/Spiroplasma code (TGA codes for Trp).

    --out <string>    default: none (no detailed report is written)
    Set a file path; the program writes a detailed report of every problematic transcript to this file. The columns are, in order: gene ID, transcript ID, sequence ID, strand, CDS start, CDS end, number of CDS features, CDS length after removing the phase, terminal status (5' partial / 3' partial / both partial / both complete), first codon, last codon, number of internal stop codons, positions of internal stop codons (amino acid index), whether a frameshift is present, and a note describing the reasons.

    --help    default: None
    Display the English usage and exit.

    --chinese_help    default: None
    显示中文用法说明并退出。
