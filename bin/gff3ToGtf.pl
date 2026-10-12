#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use constant { C_S => 0, C_E => 1, C_SC => 2, C_PH => 3, C_SRC => 4, C_SEQ => 5, C_STR => 6, C_TY => 7 };

########################################################################
# 中文用法说明（放在代码最前面；英文用法说明放在文件最末尾的__DATA__部分）
########################################################################
my $usage_cn = <<'USAGE';
程序名称：gff3ToGtf.pl

程序功能：将GFF3基因模型转换为GTF格式，并根据基因组序列推断和校验start_codon、stop_codon，同时在标准错误输出中给出警告信息和统计信息。

用法：
    perl gff3ToGtf.pl -g file.gff3 -f genome.fasta -o file.gtf [其它参数]

    perl gff3ToGtf.pl genome.fasta file.gff3 > file.gtf

【输入输出参数】
    -g, --gff 文件
    输入的GFF3基因模型文件，支持.gz压缩格式，必须提供。

    -f, --fasta 文件
    基因组FASTA文件，支持.gz压缩格式，用于推断和校验起始密码子与终止密码子，除非使用--no_codon，否则必须提供。

    -o, --out 文件
    输出的GTF文件，默认输出到标准输出。

【密码子参数】
    --start 列表
    合法的起始密码子，多个用逗号分隔，默认为ATG，例如ATG,GTG,TTG。

    --stop 列表
    合法的终止密码子，多个用逗号分隔，默认为TAA,TAG,TGA；线粒体或特殊生物请按相应遗传密码表修改。

    --no_codon
    不读取基因组序列，不生成start_codon和stop_codon，只做格式转换。

【基因模型处理参数】
    --cds_stop 模式
    控制CDS是否包含终止密码子，可选keep（默认，保持输入不变）、include（CDS包含终止密码子，GFF3习惯）、exclude（CDS不包含终止密码子，GTF2.2和Ensembl习惯）。

    --infer_utr
    输入中没有UTR时，由exon与CDS的差集推断5UTR和3UTR。

    --add_gene_transcript
    额外输出gene行和transcript行（StringTie等工具需要）。

【运行与日志参数】
    --max_warn 整数
    每一类警告最多逐条显示的条数，默认为20，超出部分只计入末尾统计。

    --quiet
    不输出逐条警告，只输出末尾统计信息。

【帮助参数】
    -h, --help
    Show the usage in English and exit.

    --chinese_help
    显示中文用法说明并退出。
USAGE

# 英文用法说明在文件末尾的__DATA__部分，运行时通过DATA文件句柄读取
my $usage_en = do { local $/; <DATA> };

my $T0 = time;

# 直接运行（不带任何参数）时，给出英文用法说明
if (!@ARGV) { print STDERR $usage_en; exit 1; }

my ($opt_fasta, $opt_gff, $opt_out, $opt_help, $opt_chelp);
my $opt_start = 'ATG';
my $opt_stop  = 'TAA,TAG,TGA';
my ($opt_gt, $opt_utr, $opt_nocodon, $opt_quiet) = (0, 0, 0, 0);
my $opt_cdsstop = 'keep';
my $opt_maxwarn = 20;

GetOptions(
    'g|gff=s' => \$opt_gff, 'f|fasta=s' => \$opt_fasta, 'o|out=s' => \$opt_out,
    'start=s' => \$opt_start, 'stop=s' => \$opt_stop, 'no_codon' => \$opt_nocodon,
    'cds_stop=s' => \$opt_cdsstop, 'infer_utr' => \$opt_utr, 'add_gene_transcript' => \$opt_gt,
    'max_warn=i' => \$opt_maxwarn, 'quiet' => \$opt_quiet,
    'h|help' => \$opt_help, 'chinese_help' => \$opt_chelp,
) or do { print STDERR "\nRun with --help (English) or --chinese_help (中文) for usage.\n"; exit 1; };

if ($opt_help)  { print STDOUT $usage_en; exit 0; }
if ($opt_chelp) { print STDOUT $usage_cn; exit 0; }

# 兼容旧用法：perl gff3ToGtf.pl genome.fasta file.gff3
if (@ARGV == 2 && !defined $opt_fasta && !defined $opt_gff) { ($opt_fasta, $opt_gff) = @ARGV; }
elsif (@ARGV == 1 && $opt_nocodon && !defined $opt_gff)     { $opt_gff = $ARGV[0]; }
elsif (@ARGV) { print STDERR "ERROR: unexpected arguments: @ARGV\nRun with --help (English) or --chinese_help (中文) for usage.\n"; exit 1; }
die "ERROR: --gff is required (run with --help for usage)\n" unless defined $opt_gff;
die "ERROR: --fasta is required unless --no_codon is set (run with --help for usage)\n" unless $opt_nocodon || defined $opt_fasta;
die "ERROR: --cds_stop must be keep|include|exclude\n" unless $opt_cdsstop =~ /^(keep|include|exclude)$/;

my %START = map { uc($_) => 1 } grep { length } split /,/, $opt_start;
my %STOP  = map { uc($_) => 1 } grep { length } split /,/, $opt_stop;

# GFF3 子特征类型（大小写不敏感）-> 标准名
my %TYPE_NORM = (
    exon => 'exon', cds => 'CDS',
    five_prime_utr => 'five_prime_UTR', '5utr' => 'five_prime_UTR',
    three_prime_utr => 'three_prime_UTR', '3utr' => 'three_prime_UTR',
    utr => 'UTR', start_codon => 'start_codon', stop_codon => 'stop_codon',
);

# ----------------------------------------------------------------------
# 通用函数
# ----------------------------------------------------------------------
my (%stat, %wcount, %in_type);

sub warnx {    # 分类限流警告
    my ($cat, $msg) = @_;
    my $n = ++$wcount{$cat};
    return if $opt_quiet;
    if    ($n <= $opt_maxwarn)     { print STDERR "[WARN][$cat] $msg\n"; }
    elsif ($n == $opt_maxwarn + 1) { print STDERR "[WARN][$cat] 同类警告过多，后续不再逐条显示（见末尾统计）\n"; }
}

sub revcom {
    my $s = reverse $_[0];
    $s =~ tr/ACGTNacgtn/TGCANtgcan/;
    return $s;
}

sub decode { my $s = $_[0]; $s =~ s/%([0-9A-Fa-f]{2})/chr(hex($1))/ge; return $s; }
sub qv     { my $s = $_[0]; $s =~ s/"/'/g; return $s; }    # GTF 属性值内不允许双引号

sub parse_attr {
    my %h;
    for my $kv (split /;/, $_[0]) {
        next unless $kv =~ /^\s*([^=\s]+)\s*=\s*(.*?)\s*$/;
        $h{$1} = $2 unless exists $h{$1};
    }
    return \%h;
}

sub open_in {
    my ($f) = @_;
    my $fh;
    if ($f =~ /\.gz$/) { open($fh, '-|', 'gzip', '-dc', $f) or die "ERROR: cannot open $f: $!\n"; }
    else               { open($fh, '<', $f) or die "ERROR: cannot open $f: $!\n"; }
    return $fh;
}

# 以转录方向给出的 CDS 片段 [s,e] 列表，取前/后 n 个碱基对应的基因组区间（按转录方向排列）
sub head_pieces {
    my ($segs, $minus, $n) = @_;
    my @p;
    for my $g (@$segs) {
        last if $n <= 0;
        my $len = $g->[1] - $g->[0] + 1;
        my $t = $len < $n ? $len : $n;
        push @p, $minus ? [$g->[1] - $t + 1, $g->[1]] : [$g->[0], $g->[0] + $t - 1];
        $n -= $t;
    }
    return $n > 0 ? () : @p;
}
sub tail_pieces {
    my ($segs, $minus, $n) = @_;
    my @p;
    for my $g (reverse @$segs) {
        last if $n <= 0;
        my $len = $g->[1] - $g->[0] + 1;
        my $t = $len < $n ? $len : $n;
        unshift @p, $minus ? [$g->[0], $g->[0] + $t - 1] : [$g->[1] - $t + 1, $g->[1]];
        $n -= $t;
    }
    return $n > 0 ? () : @p;
}
# 取序列（转录方向）；越界返回 undef
sub pieces_seq {
    my ($sref, $pieces, $minus) = @_;
    my $out = '';
    my $L = length($$sref);
    for my $p (@$pieces) {
        return undef if $p->[0] < 1 || $p->[1] > $L;
        my $s = substr($$sref, $p->[0] - 1, $p->[1] - $p->[0] + 1);
        $s = revcom($s) if $minus;
        $out .= $s;
    }
    return $out;
}

# ----------------------------------------------------------------------
# 1. 读取 GFF3（一次扫描；子特征只用轻量正则取 Parent，不做完整属性解析）
# ----------------------------------------------------------------------
my (%feat, @feat_ids, %children, %is_parent, %need_seq, %gff_seq_order);
my $gff_n = 0;
{
    my $fh = open_in($opt_gff);
    my $ln = 0;
    while (my $line = <$fh>) {
        $ln++;
        last if $line =~ /^##FASTA/;          # GFF3 末尾内嵌的序列部分
        next if $line =~ /^#/;
        $line =~ s/[\r\n]+$//;
        next if $line =~ /^\s*$/;
        $stat{lines}++;

        my @c = split /\t/, $line, 9;
        if (@c < 9) {
            warnx('bad_line', "line $ln: 少于9列(需TAB分隔)，已跳过");
            $stat{bad_line}++;
            next;
        }
        my ($seq, $src, $type, $s, $e, $score, $strand, $phase, $attr) = @c;
        if ($s !~ /^\d+$/ || $e !~ /^\d+$/ || $s < 1) {
            warnx('bad_coord', "line $ln: 坐标非法($s,$e)，已跳过");
            $stat{bad_line}++;
            next;
        }
        if ($s > $e) { warnx('swap_coord', "line $ln: start>end，已交换"); ($s, $e) = ($e, $s); $stat{swapped}++; }
        $strand = '.' unless $strand eq '+' || $strand eq '-';
        $need_seq{$seq} = 1;
        $gff_seq_order{$seq} //= $gff_n++;
        $in_type{$type}++;

        my $norm = $TYPE_NORM{lc $type};
        if ($norm) {
            my ($pstr) = $attr =~ /(?:^|;)\s*Parent=([^;]*)/;
            if (!defined $pstr || $pstr eq '') {
                warnx('no_parent', "line $ln: $type 缺少 Parent，已跳过");
                $stat{child_no_parent}++;
                next;
            }
            my $rec = [$s, $e, $score, $phase, $src, $seq, $strand, $norm];
            push @{ $children{ decode($_) } }, $rec for split /,/, $pstr;
        }
        else {
            my $a  = parse_attr($attr);
            my $id = $a->{ID};
            if (!defined $id) {
                $stat{feature_no_id}++;
                warnx('no_id', "line $ln: $type 没有 ID，无法作为 gene/transcript，已忽略") if $type =~ /gene|RNA|transcript/i;
                next;
            }
            $id = decode($id);
            if (my $o = $feat{$id}) {          # 重复 ID：合并范围并提示
                warnx('dup_id', "line $ln: 重复 ID '$id'，已合并范围");
                $o->{start} = $s if $s < $o->{start};
                $o->{end}   = $e if $e > $o->{end};
                next;
            }
            my @par = defined $a->{Parent} ? map { decode($_) } split /,/, $a->{Parent} : ();
            $feat{$id} = {
                type => $type, seq => $seq, src => $src, start => $s, end => $e, strand => $strand,
                name => (defined $a->{Name} ? decode($a->{Name}) : undef), parents => \@par,
            };
            push @feat_ids, $id;
            $is_parent{$_} = 1 for @par;
        }
    }
    close $fh;
}

# ----------------------------------------------------------------------
# 2. 读取基因组（只载入 GFF3 用到的序列，节省内存；顺序用于输出排序）
# ----------------------------------------------------------------------
my (%fasta, %fasta_order);
unless ($opt_nocodon) {
    my $fh = open_in($opt_fasta);
    my ($id, @buf, $keep);
    my $n = 0;
    my $flush = sub { $fasta{$id} = uc join('', @buf) if defined $id && $keep; @buf = (); };
    while (<$fh>) {
        if (/^>(\S+)/) {
            $flush->();
            $id = $1;
            warnx('dup_fasta', "FASTA 中序列名重复: $id") if exists $fasta_order{$id};
            $fasta_order{$id} = $n++;
            $keep = $need_seq{$id};
            next;
        }
        next unless $keep;
        tr/ \t\r\n//d;
        push @buf, $_;
    }
    $flush->();
    close $fh;
    $stat{fasta_seqs_total}  = $n;
    $stat{fasta_seqs_loaded} = scalar keys %fasta;
}
sub seq_rank { my $s = shift; return exists $fasta_order{$s} ? $fasta_order{$s} : 1e9 + ($gff_seq_order{$s} // 0); }

# ----------------------------------------------------------------------
# 3. 单个转录本 -> GTF 行
# ----------------------------------------------------------------------
sub process_transcript {
    my ($tid, $gid, $gname) = @_;
    my $f    = $feat{$tid};
    my $kids = $children{$tid};

    my (@exon, @cds, @utr);
    for my $k (@$kids) {
        my $t = $k->[C_TY];
        if    ($t eq 'exon') { push @exon, $k; }
        elsif ($t eq 'CDS')  { push @cds, [@$k]; }                # 复制，避免多 Parent 共享时被改写
        elsif ($t eq 'start_codon' || $t eq 'stop_codon') { $stat{input_codon_ignored}++; }
        else { push @utr, $k; }
    }

    my $seqid  = $f->{seq};
    my $strand = $f->{strand};
    if ($strand eq '.') { for my $k (@$kids) { if ($k->[C_STR] ne '.') { $strand = $k->[C_STR]; last; } } }
    my $minus = $strand eq '-';
    warnx('no_strand', "$tid: 链方向未知('.')，不生成 start/stop codon") if $strand eq '.';
    my $mm = grep { $_->[C_STR] ne '.' && $_->[C_STR] ne $strand } @$kids;
    warnx('strand_conflict', "$tid: 有 $mm 个子特征的链方向与转录本不一致，统一使用转录本链方向") if $mm;

    @exon = sort { $a->[C_S] <=> $b->[C_S] } @exon;
    @cds  = sort { $a->[C_S] <=> $b->[C_S] } @cds;

    # GTF 必须有 exon：若 GFF3 只给了 CDS/UTR，则由它们合并推断
    if (!@exon) {
        my @iv = sort { $a->[0] <=> $b->[0] } map { [$_->[C_S], $_->[C_E]] } (@cds, @utr);
        if (!@iv) { warnx('no_feature', "$tid: 没有 exon/CDS/UTR，已跳过"); $stat{tx_skipped_empty}++; return undef; }
        my @m = ([@{ $iv[0] }]);
        for my $i (1 .. $#iv) {
            if ($iv[$i][0] <= $m[-1][1] + 1) { $m[-1][1] = $iv[$i][1] if $iv[$i][1] > $m[-1][1]; }
            else { push @m, [@{ $iv[$i] }]; }
        }
        @exon = map { [$_->[0], $_->[1], '.', '.', $f->{src}, $seqid, $strand, 'exon'] } @m;
        $stat{tx_exon_inferred}++;
        warnx('exon_inferred', "$tid: 无 exon，已由 CDS/UTR 合并推断");
    }

    # ---- CDS 校验、frame 计算（按转录方向累加；GTF frame = 到下一密码子起点需跳过的碱基数）----
    my @cds_t = $minus ? reverse @cds : @cds;
    my ($p0, $consumed) = (0, 0);
    if (@cds) {
        for my $i (1 .. $#cds) {
            warnx('cds_overlap', "$tid: CDS 片段重叠") if $cds[$i][C_S] <= $cds[$i - 1][C_E];
        }
        for my $c (@cds) {
            my $in = 0;
            for my $x (@exon) { if ($x->[C_S] <= $c->[C_S] && $x->[C_E] >= $c->[C_E]) { $in = 1; last; } }
            if (!$in) { warnx('cds_outside_exon', "$tid: CDS $c->[C_S]-$c->[C_E] 不在任何 exon 内"); $stat{cds_outside_exon}++; }
        }
        $p0 = $cds_t[0][C_PH] =~ /^[012]$/ ? $cds_t[0][C_PH] : 0;
        $consumed = (3 - $p0) % 3;
        for my $c (@cds_t) {
            my $ph = (3 - $consumed % 3) % 3;
            if ($c->[C_PH] =~ /^[012]$/ && $c->[C_PH] != $ph) {
                $stat{phase_corrected}++;
                warnx('phase_fix', "$tid: CDS $c->[C_S]-$c->[C_E] phase=$c->[C_PH] 与长度推算($ph)不符，已改为 $ph");
            }
            $c->[C_PH] = $ph;
            $consumed += $c->[C_E] - $c->[C_S] + 1;
        }
    }

    # ---- UTR：重命名/按位置判定 5'/3'；可选推断 ----
    my ($cmin, $cmax);
    if (@cds) { $cmin = $cds[0][C_S]; for (@cds) { $cmax = $_->[C_E] if !defined $cmax || $_->[C_E] > $cmax; } }
    my @utr_out;    # [type, s, e, score, src]
    for my $u (@utr) {
        my $ty = $u->[C_TY];
        my $o;
        if    ($ty eq 'five_prime_UTR')  { $o = '5UTR'; }
        elsif ($ty eq 'three_prime_UTR') { $o = '3UTR'; }
        elsif (@cds) { $o = ($minus ? $u->[C_S] > $cmax : $u->[C_E] < $cmin) ? '5UTR' : '3UTR'; }
        else { $stat{utr_dropped_noncoding}++; next; }
        push @utr_out, [$o, $u->[C_S], $u->[C_E], $u->[C_SC], $u->[C_SRC]];
    }
    if ($opt_utr && !@utr_out && @cds) {
        for my $x (@exon) {
            push @utr_out, [$minus ? '3UTR' : '5UTR', $x->[C_S], ($x->[C_E] < $cmin ? $x->[C_E] : $cmin - 1), '.', $f->{src}] if $x->[C_S] < $cmin;
            push @utr_out, [$minus ? '5UTR' : '3UTR', ($x->[C_S] > $cmax ? $x->[C_S] : $cmax + 1), $x->[C_E], '.', $f->{src}] if $x->[C_E] > $cmax;
        }
        $stat{tx_utr_inferred}++ if @utr_out;
    }

    my $attr = sprintf('gene_id "%s"; transcript_id "%s"; gene_name "%s";', qv($gid), qv($tid), qv($gname));
    my @rows;    # [rank, start, text, subrank]
    my $add = sub {
        my ($rank, $sub, $type, $s, $e, $score, $src, $frame) = @_;
        push @rows, [$rank, $s, join("\t", $seqid, $src, $type, $s, $e, $score, $strand, $frame, $attr), $sub];
    };

    # ---- start / stop codon ----
    if (@cds && !$opt_nocodon) {
        my $sref = exists $fasta{$seqid} ? \$fasta{$seqid} : undef;
        if ($strand ne '+' && $strand ne '-') { $stat{codon_skipped_nostrand}++; }
        elsif (!$sref) {
            warnx('missing_seq', "序列 '$seqid' 不在基因组 FASTA 中，相关转录本不生成 start/stop codon");
            $stat{codon_skipped_noseq}++;
        }
        else {
            my @segs = map { [$_->[C_S], $_->[C_E]] } @cds_t;
            my $emit = sub {
                my ($type, $rank, $pieces) = @_;
                my $used = 0;
                for my $p (@$pieces) {
                    $add->($rank, 0, $type, $p->[0], $p->[1], '.', $cds_t[0][C_SRC], (3 - $used) % 3);
                    $used += $p->[1] - $p->[0] + 1;
                }
            };

            # start codon：只有 5' 完整（首段 phase=0）的 CDS 才应有 start codon
            if ($p0 != 0) { $stat{partial5}++; }
            else {
                my @pc = head_pieces(\@segs, $minus, 3);
                my $sq = @pc ? pieces_seq($sref, \@pc, $minus) : undef;
                if (!defined $sq) { $stat{codon_out_of_range}++; warnx('codon_range', "$tid: start codon 区间越界/CDS<3bp"); }
                elsif ($START{$sq}) { $emit->('start_codon', 2, \@pc); $stat{start_codon}++; $stat{start_split}++ if @pc > 1; }
                else { $stat{start_not_valid}++; warnx('start_invalid', "$tid ($strand): 起始三联体 '$sq' 不在 {$opt_start}，未输出 start_codon"); }
            }

            # stop codon：先看 CDS 末 3bp（GFF3 习惯：CDS 含终止子），再看 CDS 之后紧邻的 3bp（GTF 习惯：CDS 不含终止子）
            if ($consumed % 3 != 0) {
                $stat{cds_not_mult3}++;
                warnx('cds_len', "$tid: CDS 总长(含 phase 偏移)不是 3 的倍数，视为 3' 不完整，未输出 stop_codon");
            }
            else {
                my (@sp, $in_cds);
                my @tp = tail_pieces(\@segs, $minus, 3);
                my $tq = @tp ? pieces_seq($sref, \@tp, $minus) : undef;
                if (defined $tq && $STOP{$tq}) { @sp = @tp; $in_cds = 1; }
                else {
                    my $g = $segs[-1];
                    my ($da, $db) = $minus ? ($g->[0] - 3, $g->[0] - 1) : ($g->[1] + 1, $g->[1] + 3);
                    my $ok = 0;
                    for my $x (@exon) {
                        if ($x->[C_S] <= $g->[0] && $x->[C_E] >= $g->[1]) { $ok = $minus ? $x->[C_S] <= $da : $x->[C_E] >= $db; last; }
                    }
                    if ($ok) {
                        my $dq = pieces_seq($sref, [[$da, $db]], $minus);
                        @sp = ([$da, $db]) if defined $dq && $STOP{$dq};
                    }
                }
                if (@sp) {
                    $emit->('stop_codon', 4, \@sp);
                    $stat{stop_codon}++;
                    $stat{stop_split}++ if @sp > 1;
                    $in_cds ? $stat{stop_in_cds}++ : $stat{stop_outside_cds}++;
                    if ($opt_cdsstop eq 'exclude' && $in_cds) {          # 从 CDS 末端削去 3bp
                        my $n = 3;
                        while ($n > 0 && @cds_t) {
                            my $l = $cds_t[-1];
                            my $len = $l->[C_E] - $l->[C_S] + 1;
                            if ($len <= $n) { pop @cds_t; $n -= $len; }
                            else { if ($minus) { $l->[C_S] += $n; } else { $l->[C_E] -= $n; } $n = 0; }
                        }
                        $stat{cds_trimmed}++;
                    }
                    elsif ($opt_cdsstop eq 'include' && !$in_cds) {      # 把终止子并入 CDS
                        my $l = $cds_t[-1];
                        if ($minus) { $l->[C_S] -= 3; } else { $l->[C_E] += 3; }
                        $stat{cds_extended}++;
                    }
                }
                else { $stat{no_stop}++; warnx('no_stop', "$tid ($strand): 未找到终止密码子(末3bp='" . (defined $tq ? $tq : 'NA') . "')，可能 3' 不完整或终止子被拆在剪接位点后"); }
            }
        }
    }

    # ---- 其它行 ----
    for my $x (@exon)  { $add->(3, 0, 'exon', $x->[C_S], $x->[C_E], $x->[C_SC], $x->[C_SRC], '.'); $stat{out_exon}++; }
    for my $c (@cds_t) { $add->(3, 1, 'CDS',  $c->[C_S], $c->[C_E], $c->[C_SC], $c->[C_SRC], $c->[C_PH]); $stat{out_CDS}++; }
    for my $u (@utr_out) {
        $add->($u->[0] eq '5UTR' ? 1 : 5, 0, $u->[0], $u->[1], $u->[2], $u->[3], $u->[4], '.');
        $stat{"out_$u->[0]"}++;
    }

    # 排序：5UTR -> start_codon -> exon/CDS -> stop_codon -> 3UTR；+链升序，-链降序
    @rows = $minus
        ? sort { $a->[0] <=> $b->[0] || $b->[1] <=> $a->[1] || $a->[3] <=> $b->[3] } @rows
        : sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] || $a->[3] <=> $b->[3] } @rows;

    my ($tmin, $tmax) = ($exon[0][C_S], $exon[0][C_E]);
    for (@exon) { $tmin = $_->[C_S] if $_->[C_S] < $tmin; $tmax = $_->[C_E] if $_->[C_E] > $tmax; }
    my $tline = join("\t", $seqid, $f->{src}, 'transcript', $tmin, $tmax, '.', $strand, '.', $attr);

    @cds ? $stat{coding_tx}++ : $stat{noncoding_tx}++;
    @exon > 1 ? $stat{multi_exon_tx}++ : $stat{single_exon_tx}++;
    $stat{tx_written}++;

    return { id => $tid, seq => $seqid, strand => $strand, start => $tmin, end => $tmax,
             text => join("\n", map { $_->[2] } @rows), tline => $tline };
}

# ----------------------------------------------------------------------
# 4. 组织 gene -> transcript 并输出（按基因组位置排序）
# ----------------------------------------------------------------------
for my $p (keys %children) {
    next if exists $feat{$p};
    warnx('missing_parent', "Parent '$p' 在 GFF3 中没有对应的 ID 行，其 " . scalar(@{ $children{$p} }) . " 个子特征被丢弃");
    $stat{orphan_children} += scalar @{ $children{$p} };
}

my %genes;
for my $tid (@feat_ids) {
    my $f = $feat{$tid};
    if (!exists $children{$tid}) {
        next if $is_parent{$tid};                       # 基因：由其转录本带出
        if ($f->{type} =~ /gene/i) {
            $stat{gene_without_transcript}++;
            warnx('gene_no_tx', "gene '$tid' 没有任何转录本/子特征，已忽略");
        }
        elsif ($f->{type} =~ /RNA|transcript/i) { warnx('tx_no_child', "$f->{type} '$tid' 没有 exon/CDS/UTR 子特征，已忽略"); $stat{tx_skipped_empty}++; }
        else { $stat{ignored_type}{ $f->{type} }++; }
        next;
    }
    my @par = @{ $f->{parents} };
    my $gid;
    if (@par) {
        $gid = $par[0];
        warnx('multi_parent', "$tid 有多个 Parent(" . join(',', @par) . ")，GTF 一个 transcript 只能属于一个基因，使用 '$gid'") if @par > 1;
    }
    else { $gid = $tid; $stat{tx_without_gene}++; }
    my $gname = (exists $feat{$gid} && defined $feat{$gid}{name}) ? $feat{$gid}{name} : $gid;

    my $r = process_transcript($tid, $gid, $gname) or next;
    my $g = ($genes{$gid} ||= { seq => $r->{seq}, strand => $r->{strand}, start => $r->{start}, end => $r->{end},
                                name => $gname, tx => [] });
    $g->{start} = $r->{start} if $r->{start} < $g->{start};
    $g->{end}   = $r->{end}   if $r->{end}   > $g->{end};
    push @{ $g->{tx} }, $r;
}

my $out_fh;
if (defined $opt_out) { open($out_fh, '>', $opt_out) or die "ERROR: cannot write $opt_out: $!\n"; }
else { $out_fh = \*STDOUT; }

my @gids = sort {
    seq_rank($genes{$a}{seq}) <=> seq_rank($genes{$b}{seq})
        || $genes{$a}{start} <=> $genes{$b}{start} || $a cmp $b
} keys %genes;

for my $gid (@gids) {
    my $g = $genes{$gid};
    my @tx = sort { $a->{start} <=> $b->{start} || $a->{id} cmp $b->{id} } @{ $g->{tx} };
    my $buf = '';
    if ($opt_gt) {
        my $gf = $feat{$gid};
        my ($gs, $ge, $src) = ($g->{start}, $g->{end}, $tx[0]{tline} =~ /^[^\t]+\t([^\t]*)/ ? $1 : '.');
        if ($gf) { ($gs, $ge, $src) = ($gf->{start}, $gf->{end}, $gf->{src}); }
        $buf .= join("\t", $g->{seq}, $src, 'gene', $gs, $ge, '.', $g->{strand}, '.',
                     sprintf('gene_id "%s"; gene_name "%s";', qv($gid), qv($g->{name}))) . "\n";
        $stat{out_gene_line}++;
    }
    for my $t (@tx) {
        $buf .= $t->{tline} . "\n" if $opt_gt;
        $buf .= $t->{text} . "\n";
    }
    print $out_fh $buf;       # 每个基因一次写出，减少 I/O 次数
    $stat{genes_written}++;
}
close $out_fh if defined $opt_out;

# ----------------------------------------------------------------------
# 5. 统计信息（STDERR）
# ----------------------------------------------------------------------
sub st { return $stat{ $_[0] } // 0; }
my $elapsed = time - $T0;
my $E = \*STDERR;
print $E "\n================ gff3ToGtf 转换统计 ================\n";
printf $E "%-38s %s\n", '输入 GFF3', $opt_gff;
printf $E "%-38s %d\n", '有效数据行', st('lines');
printf $E "%-38s %d\n", '格式错误/非法行(已跳过)', st('bad_line');
printf $E "%-38s %d\n", '起止坐标颠倒(已交换)', st('swapped');
if (%in_type) {
    print $E "输入特征类型计数:\n";
    printf $E "    %-34s %d\n", $_, $in_type{$_} for sort keys %in_type;
}
unless ($opt_nocodon) {
    printf $E "%-38s %d / %d\n", 'FASTA 序列(载入/总数)', st('fasta_seqs_loaded'), st('fasta_seqs_total');
}
print $E "--- 输出 ---\n";
printf $E "%-38s %d\n", '基因数', st('genes_written');
printf $E "%-38s %d\n", '转录本数', st('tx_written');
printf $E "%-38s %d / %d\n", '  编码 / 非编码', st('coding_tx'), st('noncoding_tx');
printf $E "%-38s %d / %d\n", '  单外显子 / 多外显子', st('single_exon_tx'), st('multi_exon_tx');
printf $E "%-38s %.2f\n", '平均每基因转录本数', st('genes_written') ? st('tx_written') / st('genes_written') : 0;
printf $E "%-38s %d\n", 'exon 行', st('out_exon');
printf $E "%-38s %d\n", 'CDS 行', st('out_CDS');
printf $E "%-38s %d / %d\n", '5UTR / 3UTR 行', st('out_5UTR'), st('out_3UTR');
printf $E "%-38s %d (跨内含子 %d)\n", 'start_codon 转录本数', st('start_codon'), st('start_split');
printf $E "%-38s %d (跨内含子 %d)\n", 'stop_codon 转录本数', st('stop_codon'), st('stop_split');
printf $E "%-38s %d / %d\n", '  终止子在CDS内 / 在CDS外', st('stop_in_cds'), st('stop_outside_cds');
printf $E "%-38s %d\n", 'gene/transcript 行', st('out_gene_line') + ($opt_gt ? st('tx_written') : 0) if $opt_gt;
print $E "--- 修复/推断 ---\n";
printf $E "%-38s %d\n", '无 exon，由 CDS/UTR 推断 exon', st('tx_exon_inferred');
printf $E "%-38s %d\n", '推断 UTR 的转录本', st('tx_utr_inferred');
printf $E "%-38s %d\n", 'CDS phase 被重算修正', st('phase_corrected');
printf $E "%-38s %s\n", "CDS 被削去终止子 / 被延伸含终止子", st('cds_trimmed') . ' / ' . st('cds_extended') if $opt_cdsstop ne 'keep';
printf $E "%-38s %d\n", '输入已有 start/stop_codon(忽略并重算)', st('input_codon_ignored');
print $E "--- 问题 ---\n";
printf $E "%-38s %d\n", "5' 不完整(首段phase≠0)", st('partial5');
printf $E "%-38s %d\n", '起始三联体不合法', st('start_not_valid');
printf $E "%-38s %d\n", '未找到终止密码子', st('no_stop');
printf $E "%-38s %d\n", 'CDS 长度非3倍数', st('cds_not_mult3');
printf $E "%-38s %d\n", 'CDS 不在 exon 内', st('cds_outside_exon');
printf $E "%-38s %d\n", '密码子区间越界', st('codon_out_of_range');
printf $E "%-38s %d\n", '因序列缺失跳过密码子', st('codon_skipped_noseq');
printf $E "%-38s %d\n", '因链未知跳过密码子', st('codon_skipped_nostrand');
printf $E "%-38s %d\n", '无子特征被跳过的转录本', st('tx_skipped_empty');
printf $E "%-38s %d\n", '无 Parent 的子特征', st('child_no_parent');
printf $E "%-38s %d\n", 'Parent 缺失而丢弃的子特征', st('orphan_children');
printf $E "%-38s %d\n", '无转录本的 gene', st('gene_without_transcript');
printf $E "%-38s %d\n", '无 gene 的独立转录本', st('tx_without_gene');
if ($stat{ignored_type}) {
    print $E "被忽略的无关特征类型:\n";
    printf $E "    %-34s %d\n", $_, $stat{ignored_type}{$_} for sort keys %{ $stat{ignored_type} };
}
if (%wcount) {
    print $E "警告分类计数:\n";
    printf $E "    %-34s %d\n", $_, $wcount{$_} for sort keys %wcount;
}
else { print $E "警告: 无\n"; }
printf $E "运行耗时: %d 秒\n", $elapsed;
print $E "====================================================\n";

########################################################################
# 英文用法说明（放在文件最末尾，运行时通过DATA文件句柄读取）
########################################################################
__DATA__
Program: gff3ToGtf.pl

Description: Convert gene models from GFF3 to GTF, infer and validate start_codon and stop_codon from the genome sequence, and report warnings and summary statistics to standard error.

Usage:
    perl gff3ToGtf.pl -g file.gff3 -f genome.fasta -o file.gtf [other options]

    perl gff3ToGtf.pl genome.fasta file.gff3 > file.gtf

[Input and Output Options]
    -g, --gff FILE
    Input GFF3 gene model file; .gz compressed files are supported; required.

    -f, --fasta FILE
    Genome FASTA file; .gz compressed files are supported; used to infer and validate start and stop codons; required unless --no_codon is set.

    -o, --out FILE
    Output GTF file; default is standard output.

[Codon Options]
    --start LIST
    Valid start codons, comma separated; default is ATG, e.g. ATG,GTG,TTG.

    --stop LIST
    Valid stop codons, comma separated; default is TAA,TAG,TGA; change it according to the genetic code table for mitochondria or special organisms.

    --no_codon
    Do not read the genome sequence and do not generate start_codon or stop_codon; only convert the format.

[Gene Model Processing Options]
    --cds_stop MODE
    Control whether CDS includes the stop codon: keep (default, keep the input unchanged), include (CDS includes the stop codon, GFF3 convention), exclude (CDS excludes the stop codon, GTF2.2 and Ensembl convention).

    --infer_utr
    Infer 5UTR and 3UTR from the difference between exon and CDS when the input has no UTR.

    --add_gene_transcript
    Additionally output gene and transcript lines (required by tools such as StringTie).

[Run and Log Options]
    --max_warn INT
    Maximum number of warnings displayed per warning category; default is 20; the rest are only counted in the final summary.

    --quiet
    Do not print individual warnings; only print the final summary.

[Help Options]
    -h, --help
    Show the usage in English and exit.

    --chinese_help
    显示中文用法说明并退出。

