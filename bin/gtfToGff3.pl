#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

########################################################################
# 中文用法说明（放在代码最前面；英文用法说明放在文件最末尾的__DATA__部分）
########################################################################
my $usage_cn = <<'USAGE';
程序名称：gtfToGff3.pl

程序功能：将GTF基因模型转换为GFF3格式（gff3ToGtf.pl的逆过程），自动补全gene、transcript（mRNA等）层级关系、ID和Parent，由exon和CDS推断phase，并在标准错误输出中给出警告信息和统计信息。

用法：
    perl gtfToGff3.pl -i file.gtf -o file.gff3 [其它参数]

    perl gtfToGff3.pl file.gtf > file.gff3

【输入输出参数】
    -i, --gtf 文件
    输入的GTF基因模型文件，支持.gz压缩格式，必须提供，也可以直接作为位置参数给出。

    -o, --out 文件
    输出的GFF3文件，默认输出到标准输出。

【基因模型转换参数】
    --cds_stop 模式
    控制CDS是否包含终止密码子，可选include（默认，GFF3习惯，CDS包含终止密码子，输入CDS不含时由stop_codon行补入）、keep（保持输入不变）、exclude（CDS不包含终止密码子，输入CDS含有时削去）。

    --keep_codon
    保留start_codon和stop_codon，作为转录本的子特征输出，默认丢弃（GFF3通常不单独记录起始和终止密码子）。

    --infer_utr
    输入中没有UTR时，由exon与CDS的差集推断five_prime_UTR和three_prime_UTR。

    --coding_type 字符串
    编码转录本（含CDS）的类型，默认为mRNA。

    --noncoding_type 字符串
    无法由GTF中的类型或biotype判断的非编码转录本的类型，默认为ncRNA。

    --no_extra_attr
    不保留GTF中除gene_id、transcript_id、gene_name、transcript_name之外的属性（如gene_biotype、tag、protein_id等），只输出ID、Parent和Name。

    --add_separator
    在每个基因之后输出GFF3的###分隔行，便于流式处理。

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

my ($opt_gtf, $opt_out, $opt_help, $opt_chelp);
my $opt_cdsstop = 'include';
my ($opt_codon, $opt_utr, $opt_noextra, $opt_sep, $opt_quiet) = (0, 0, 0, 0, 0);
my $opt_coding    = 'mRNA';
my $opt_noncoding = 'ncRNA';
my $opt_maxwarn   = 20;

GetOptions(
    'i|gtf=s' => \$opt_gtf, 'o|out=s' => \$opt_out,
    'cds_stop=s' => \$opt_cdsstop, 'keep_codon' => \$opt_codon, 'infer_utr' => \$opt_utr,
    'coding_type=s' => \$opt_coding, 'noncoding_type=s' => \$opt_noncoding,
    'no_extra_attr' => \$opt_noextra, 'add_separator' => \$opt_sep,
    'max_warn=i' => \$opt_maxwarn, 'quiet' => \$opt_quiet,
    'h|help' => \$opt_help, 'chinese_help' => \$opt_chelp,
) or do { print STDERR "\nRun with --help (English) or --chinese_help (中文) for usage.\n"; exit 1; };

if ($opt_help)  { print STDOUT $usage_en; exit 0; }
if ($opt_chelp) { print STDOUT $usage_cn; exit 0; }

if (@ARGV == 1 && !defined $opt_gtf) { $opt_gtf = $ARGV[0]; }
elsif (@ARGV) { print STDERR "ERROR: unexpected arguments: @ARGV\nRun with --help (English) or --chinese_help (中文) for usage.\n"; exit 1; }
die "ERROR: --gtf is required (run with --help for usage)\n" unless defined $opt_gtf;
die "ERROR: --cds_stop must be include|keep|exclude\n" unless $opt_cdsstop =~ /^(keep|include|exclude)$/;

# ----------------------------------------------------------------------
# 通用函数
# ----------------------------------------------------------------------
my (%stat, %wcount, %in_type, %ignored, %txtype_n);

sub warnx {    # 分类限流警告
    my ($cat, $msg) = @_;
    my $n = ++$wcount{$cat};
    return if $opt_quiet;
    if    ($n <= $opt_maxwarn)     { print STDERR "[WARN][$cat] $msg\n"; }
    elsif ($n == $opt_maxwarn + 1) { print STDERR "[WARN][$cat] 同类警告过多，后续不再逐条显示（见末尾统计）\n"; }
}

sub open_in {
    my ($f) = @_;
    my $fh;
    if ($f =~ /\.gz$/) { open($fh, '-|', 'gzip', '-dc', $f) or die "ERROR: cannot open $f: $!\n"; }
    else               { open($fh, '<', $f) or die "ERROR: cannot open $f: $!\n"; }
    return $fh;
}

# GFF3 转义：; = & , % 以及控制字符
sub enc { my $s = $_[0]; $s =~ s/([;=&,%\x00-\x1f\x7f])/sprintf('%%%02X', ord($1))/ge; return $s; }

# 解析 GTF 第9列：key "value"; key value; ... ；同名键允许多个值
sub parse_attr {
    my %h;
    my $s = $_[0];
    while ($s =~ /([^\s";]+)\s+(?:"([^"]*)"|([^;]*?))\s*(?:;|$)/g) {
        push @{ $h{$1} }, (defined $2 ? $2 : $3);
    }
    return \%h;
}

sub attr_str {
    my ($id, $parent, $name, $extras) = @_;
    my @a = ('ID=' . enc($id));
    push @a, 'Parent=' . enc($parent) if defined $parent;
    push @a, 'Name=' . enc($name) if defined $name && length $name;
    if (!$opt_noextra && $extras) {
        for my $k (sort keys %$extras) {
            next if $k =~ /^(?:ID|Parent|Name)$/;
            push @a, enc($k) . '=' . join(',', map { enc($_) } @{ $extras->{$k} });
        }
    }
    return join(';', @a);
}

# GTF 类型（小写）-> 内部标准类型
my %NORM = (
    gene => 'gene',
    transcript => 'transcript', mrna => 'transcript', ncrna => 'transcript', lncrna => 'transcript',
    lnc_rna => 'transcript', rrna => 'transcript', trna => 'transcript', mirna => 'transcript',
    snrna => 'transcript', snorna => 'transcript', primary_transcript => 'transcript',
    pseudogenic_transcript => 'transcript',
    exon => 'exon', cds => 'CDS',
    '5utr' => 'five_prime_UTR', five_prime_utr => 'five_prime_UTR',
    '3utr' => 'three_prime_UTR', three_prime_utr => 'three_prime_UTR',
    utr => 'UTR', start_codon => 'start_codon', stop_codon => 'stop_codon',
);
# 非编码 biotype -> GFF3 转录本类型
my %BIO = (
    lncrna => 'lnc_RNA', lincrna => 'lnc_RNA', lnc_rna => 'lnc_RNA', mirna => 'miRNA', rrna => 'rRNA',
    trna => 'tRNA', snrna => 'snRNA', snorna => 'snoRNA', scrna => 'scRNA',
);
my %NCFAMILY = map { $_ => 1 } qw(ncRNA lnc_RNA miRNA rRNA tRNA snRNA snoRNA scRNA);
my $SKIP_RE = qr/^(?:gene_id|transcript_id|exon_number|exon_id|exon_version|protein_id)$/;

# ----------------------------------------------------------------------
# 1. 读取 GTF（一次扫描；子特征行只用正则取 gene_id/transcript_id，每个转录本只完整解析一次属性）
# ----------------------------------------------------------------------
my (%gene, %tx, @gene_order, @tx_order, %seq_order);
my $seq_n = 0;

sub new_gene { my ($seq, $strand, $src) = @_; return { seq => $seq, strand => $strand, src => $src, s => undef, e => undef, has_line => 0, ltype => undef, name => undef, extras => {}, tx => [] }; }
sub new_tx   { my ($gid, $seq, $strand, $src) = @_; return { gid => $gid, seq => $seq, strand => $strand, src => $src, s => undef, e => undef, has_line => 0, ltype => undef, name => undef, pid => undef, done => 0, extras => {}, exon => [], cds => [], utr => [], cstart => [], cstop => [] }; }

sub absorb_gene {
    my ($g, $h) = @_;
    for my $k (keys %$h) {
        next if $k eq 'gene_id';
        if ($k eq 'gene_name') { $g->{name} = $h->{$k}[0] unless defined $g->{name}; }
        else { $g->{extras}{$k} = [ @{ $h->{$k} } ] unless exists $g->{extras}{$k}; }
    }
}
sub absorb_tx {
    my ($t, $g, $h, $force) = @_;
    for my $k (keys %$h) {
        next if $k =~ $SKIP_RE;
        if    ($k eq 'gene_name') { $g->{name} = $h->{$k}[0] unless defined $g->{name}; }
        elsif ($k =~ /^gene_/)    { $g->{extras}{$k} = [ @{ $h->{$k} } ] unless exists $g->{extras}{$k}; }
        elsif ($k eq 'transcript_name') { $t->{name} = $h->{$k}[0] if $force || !defined $t->{name}; }
        else { $t->{extras}{$k} = [ @{ $h->{$k} } ] if $force || !exists $t->{extras}{$k}; }
    }
    $t->{done} = 1;
}

{
    my $fh = open_in($opt_gtf);
    my $ln = 0;
    while (my $line = <$fh>) {
        $ln++;
        next if $line =~ /^#/;
        $line =~ s/[\r\n]+$//;
        next if $line =~ /^\s*$/;
        $stat{lines}++;

        my @c = split /\t/, $line, 9;
        if (@c < 9) { warnx('bad_line', "line $ln: 少于9列(需TAB分隔)，已跳过"); $stat{bad_line}++; next; }
        my ($seq, $src, $type, $s, $e, $score, $strand, $frame, $attr) = @c;
        if ($s !~ /^\d+$/ || $e !~ /^\d+$/ || $s < 1) { warnx('bad_coord', "line $ln: 坐标非法($s,$e)，已跳过"); $stat{bad_line}++; next; }
        if ($s > $e) { warnx('swap_coord', "line $ln: start>end，已交换"); ($s, $e) = ($e, $s); $stat{swapped}++; }
        $strand = '.' unless $strand eq '+' || $strand eq '-';
        $in_type{$type}++;

        my $norm = $NORM{ lc $type };
        if (!$norm) { $ignored{$type}++; next; }
        $seq_order{$seq} //= $seq_n++;

        my ($g1, $g2) = $attr =~ /(?:^|;)\s*gene_id\s+(?:"([^"]*)"|([^;\s]+))/;
        my $gid = defined $g1 ? $g1 : $g2;

        if ($norm eq 'gene') {
            if (!defined $gid || $gid eq '') { warnx('gene_no_id', "line $ln: gene 行缺少 gene_id，已跳过"); $stat{line_no_id}++; next; }
            my $g = $gene{$gid};
            if (!$g) { $g = $gene{$gid} = new_gene($seq, $strand, $src); push @gene_order, $gid; }
            warnx('dup_gene', "line $ln: gene '$gid' 有多条 gene 行，取最后一条坐标") if $g->{has_line};
            @$g{qw(has_line s e src seq strand)} = (1, $s, $e, $src, $seq, $strand);
            $g->{ltype} = $type if lc($type) ne 'gene';
            absorb_gene($g, parse_attr($attr));
            next;
        }

        my ($t1, $t2) = $attr =~ /(?:^|;)\s*transcript_id\s+(?:"([^"]*)"|([^;\s]+))/;
        my $tid = defined $t1 ? $t1 : $t2;
        if (!defined $tid || $tid eq '') { warnx('no_tid', "line $ln: $type 行缺少 transcript_id，已跳过"); $stat{line_no_id}++; next; }
        if (!defined $gid || $gid eq '') {
            $gid = $tid;
            $stat{gene_id_missing}++;
            warnx('no_gid', "line $ln: transcript '$tid' 缺少 gene_id，使用 transcript_id 作为 gene_id");
        }

        my $t = $tx{$tid};
        if (!$t) {
            $t = $tx{$tid} = new_tx($gid, $seq, $strand, $src);
            push @tx_order, $tid;
            my $g = $gene{$gid};
            if (!$g) { $g = $gene{$gid} = new_gene($seq, $strand, $src); push @gene_order, $gid; }
            push @{ $g->{tx} }, $tid;
        }
        else {
            warnx('tx_multi_gene', "line $ln: transcript '$tid' 出现在不同 gene_id ($t->{gid} / $gid) 下，使用 '$t->{gid}'") if $t->{gid} ne $gid;
            warnx('tx_multi_seq', "line $ln: transcript '$tid' 跨不同序列 ($t->{seq} / $seq)，使用 '$t->{seq}'") if $t->{seq} ne $seq;
            if ($strand ne '.' && $t->{strand} eq '.') { $t->{strand} = $strand; }
            elsif ($strand ne '.' && $strand ne $t->{strand}) { warnx('strand_conflict', "line $ln: transcript '$tid' 链方向不一致，使用 '$t->{strand}'"); }
        }
        my $g = $gene{ $t->{gid} };

        if ($norm eq 'transcript') {
            warnx('dup_tx', "line $ln: transcript '$tid' 有多条 transcript 行，取最后一条坐标") if $t->{has_line};
            @$t{qw(has_line s e src)} = (1, $s, $e, $src);
            $t->{ltype} = $type if lc($type) ne 'transcript';
            absorb_tx($t, $g, parse_attr($attr), 1);
            next;
        }

        my $rec = [$s, $e, $score, $frame, $src, $norm];
        if    ($norm eq 'exon') { push @{ $t->{exon} }, $rec; }
        elsif ($norm eq 'CDS') {
            push @{ $t->{cds} }, $rec;
            if (!defined $t->{pid}) { my ($p1, $p2) = $attr =~ /(?:^|;)\s*protein_id\s+(?:"([^"]*)"|([^;\s]+))/; $t->{pid} = defined $p1 ? $p1 : $p2; }
        }
        elsif ($norm eq 'start_codon') { push @{ $t->{cstart} }, $rec; }
        elsif ($norm eq 'stop_codon')  { push @{ $t->{cstop} },  $rec; }
        else { push @{ $t->{utr} }, $rec; }
        absorb_tx($t, $g, parse_attr($attr), 0) unless $t->{done};
    }
    close $fh;
}

# ----------------------------------------------------------------------
# 2. 单个转录本 -> GFF3 行
# ----------------------------------------------------------------------
sub by_start { $a->[0] <=> $b->[0] }

# 取 biotype（不使用哈希元素别名，避免自动创建空键）
sub biotype {
    for my $h (@_) {
        for my $k (qw(transcript_biotype transcript_type gene_biotype gene_type)) {
            return lc $h->{$k}[0] if exists $h->{$k};
        }
    }
    return '';
}

sub build_tx {
    my ($tid, $t) = @_;
    my $oid = $tid;
    if (exists $gene{$tid}) {       # GFF3 的 ID 必须唯一
        $oid = "$tid.t1";
        $stat{id_collision}++;
        warnx('id_collision', "transcript_id '$tid' 与 gene_id 相同，转录本 ID 改为 '$oid'");
    }
    my ($seq, $strand, $src) = @$t{qw(seq strand src)};
    warnx('no_strand', "$tid: 链方向未知('.')，按 + 链顺序处理") if $strand eq '.';
    my $minus = $strand eq '-';

    my @exon   = sort by_start @{ $t->{exon} };
    my @cds    = sort by_start map { [@$_] } @{ $t->{cds} };
    my @startc = sort by_start @{ $t->{cstart} };
    my @stopc  = sort by_start @{ $t->{cstop} };

    # GFF3 中 exon 虽非强制，但 GTF 的 exon 缺失时由 CDS/UTR/密码子合并推断
    if (!@exon) {
        my @iv = sort by_start map { [$_->[0], $_->[1]] } (@cds, @{ $t->{utr} }, @startc, @stopc);
        if (!@iv) { warnx('no_feature', "$tid: 没有 exon/CDS/UTR，已跳过"); $stat{tx_skipped_empty}++; return undef; }
        my @m = ([@{ $iv[0] }]);
        for my $i (1 .. $#iv) {
            if ($iv[$i][0] <= $m[-1][1] + 1) { $m[-1][1] = $iv[$i][1] if $iv[$i][1] > $m[-1][1]; }
            else { push @m, [@{ $iv[$i] }]; }
        }
        @exon = map { [$_->[0], $_->[1], '.', '.', $src, 'exon'] } @m;
        $stat{tx_exon_inferred}++;
        warnx('exon_inferred', "$tid: 无 exon，已由 CDS/UTR 合并推断");
    }

    # 首个 CDS（转录方向）的 frame 作为起始偏移（5'端不完整时为1或2）
    my $p0 = 0;
    if (@cds) { my $f0 = $minus ? $cds[-1][3] : $cds[0][3]; $p0 = $f0 =~ /^[012]$/ ? $f0 : 0; }

    # ---- 终止密码子与 CDS 的关系：GTF2.2 的 CDS 不含终止子，GFF3 习惯含终止子 ----
    if (@cds && @stopc) {
        my $all_in = 1;
        for my $p (@stopc) {
            my $in = 0;
            for my $c (@cds) { if ($c->[0] <= $p->[0] && $c->[1] >= $p->[1]) { $in = 1; last; } }
            $all_in = 0 unless $in;
        }
        if ($opt_cdsstop eq 'include' && !$all_in) {
            for my $p (@stopc) {
                my $done = 0;
                for my $c (@cds) {
                    if    ($c->[0] <= $p->[0] && $c->[1] >= $p->[1]) { $done = 1; last; }
                    elsif ($p->[0] == $c->[1] + 1) { $c->[1] = $p->[1]; $done = 1; last; }
                    elsif ($p->[1] == $c->[0] - 1) { $c->[0] = $p->[0]; $done = 1; last; }
                }
                push @cds, [$p->[0], $p->[1], '.', '.', $p->[4], 'CDS'] unless $done;
            }
            @cds = sort by_start @cds;
            $stat{cds_extended}++;
        }
        elsif ($opt_cdsstop eq 'exclude' && $all_in) {
            for my $p (reverse @stopc) {
                for my $c (@cds) {
                    next if $c->[6];
                    if ($c->[0] <= $p->[0] && $c->[1] >= $p->[1]) {
                        if    ($c->[0] == $p->[0] && $c->[1] == $p->[1]) { $c->[6] = 1; }
                        elsif ($c->[1] == $p->[1]) { $c->[1] = $p->[0] - 1; }
                        elsif ($c->[0] == $p->[0]) { $c->[0] = $p->[1] + 1; }
                        else { warnx('stop_inside', "$tid: 终止密码子位于CDS片段中部，无法削去"); }
                        last;
                    }
                }
            }
            @cds = grep { !$_->[6] } @cds;
            $stat{cds_trimmed}++;
        }
    }

    # ---- phase：按转录方向由 CDS 长度重新计算 ----
    my @cds_t = $minus ? reverse @cds : @cds;
    my $consumed = (3 - $p0) % 3;
    for my $c (@cds_t) {
        my $ph = (3 - $consumed % 3) % 3;
        if ($c->[3] =~ /^[012]$/ && $c->[3] != $ph) {
            $stat{phase_corrected}++;
            warnx('phase_fix', "$tid: CDS $c->[0]-$c->[1] frame=$c->[3] 与长度推算($ph)不符，已改为 $ph");
        }
        $c->[3] = $ph;
        $consumed += $c->[1] - $c->[0] + 1;
    }
    if (@cds && $consumed % 3) {
        $stat{cds_not_mult3}++;
        warnx('cds_len', "$tid: CDS 总长(含 frame 偏移)不是 3 的倍数，可能 3' 端不完整或终止子缺失");
    }
    my $coding = @cds ? 1 : 0;

    # ---- UTR：通用 UTR 按位置判定；可选由 exon 减 CDS 推断 ----
    my ($cmin, $cmax);
    if ($coding) { $cmin = $cds[0][0]; for (@cds) { $cmax = $_->[1] if !defined $cmax || $_->[1] > $cmax; } }
    my @utr_out;    # [type, s, e, score, src]
    for my $u (@{ $t->{utr} }) {
        my $o = $u->[5];
        if ($o eq 'UTR') {
            if ($coding) { $o = ($minus ? $u->[0] > $cmax : $u->[1] < $cmin) ? 'five_prime_UTR' : 'three_prime_UTR'; }
            else { $stat{utr_dropped_noncoding}++; next; }
        }
        push @utr_out, [$o, $u->[0], $u->[1], $u->[2], $u->[4]];
    }
    if ($opt_utr && !@utr_out && $coding) {
        for my $x (@exon) {
            push @utr_out, [$minus ? 'three_prime_UTR' : 'five_prime_UTR', $x->[0], ($x->[1] < $cmin ? $x->[1] : $cmin - 1), '.', $src] if $x->[0] < $cmin;
            push @utr_out, [$minus ? 'five_prime_UTR' : 'three_prime_UTR', ($x->[0] > $cmax ? $x->[0] : $cmax + 1), $x->[1], '.', $src] if $x->[1] > $cmax;
        }
        $stat{tx_utr_inferred}++ if @utr_out;
    }

    # ---- 转录本范围：取 transcript 行与子特征的并集 ----
    my ($ts, $te) = ($exon[0][0], $exon[0][1]);
    for (@exon, @cds) { $ts = $_->[0] if $_->[0] < $ts; $te = $_->[1] if $_->[1] > $te; }
    for (@utr_out)    { $ts = $_->[1] if $_->[1] < $ts; $te = $_->[2] if $_->[2] > $te; }
    if ($t->{has_line}) {
        if ($t->{s} != $ts || $t->{e} != $te) {
            warnx('tx_range', "$tid: transcript 行坐标 $t->{s}-$t->{e} 与子特征范围 $ts-$te 不一致，取并集");
            $stat{tx_range_mismatch}++;
        }
        $ts = $t->{s} if $t->{s} < $ts;
        $te = $t->{e} if $t->{e} > $te;
    }
    else { $stat{tx_line_inferred}++; }

    # ---- 转录本类型 ----
    my $g = $gene{ $t->{gid} };
    my $type;
    if (defined $t->{ltype}) { $type = $t->{ltype}; }
    elsif ($coding)          { $type = $opt_coding; }
    else {
        my $bt = biotype($t->{extras}, $g->{extras});
        $type = $bt =~ /pseudogene/ ? 'pseudogenic_transcript' : ($BIO{$bt} // $opt_noncoding);
    }
    $txtype_n{$type}++;

    # ---- 生成各行 ----
    my @rows;    # [start, rank, text]
    my $add = sub {
        my ($rank, $type2, $s, $e, $score, $src2, $phase, $attr) = @_;
        push @rows, [$s, $rank, join("\t", $seq, $src2, $type2, $s, $e, $score, $strand, $phase, $attr)];
    };
    my $n = @exon;
    for my $i (0 .. $#exon) {
        my $num = $minus ? $n - $i : $i + 1;
        $add->(0, 'exon', $exon[$i][0], $exon[$i][1], $exon[$i][2], $exon[$i][4], '.', attr_str("$oid.exon$num", $oid));
        $stat{out_exon}++;
    }
    my $cds_attr = attr_str("cds-$oid", $oid, undef, defined $t->{pid} ? { protein_id => [$t->{pid}] } : undef);
    for my $c (@cds) { $add->(2, 'CDS', $c->[0], $c->[1], $c->[2], $c->[4], $c->[3], $cds_attr); $stat{out_CDS}++; }
    for my $ty (qw(five_prime_UTR three_prime_UTR)) {
        my @u = sort { $a->[1] <=> $b->[1] } grep { $_->[0] eq $ty } @utr_out;
        my $tag = $ty eq 'five_prime_UTR' ? 'utr5p' : 'utr3p';
        for my $i (0 .. $#u) {
            my $num = $minus ? @u - $i : $i + 1;
            $add->($ty eq 'five_prime_UTR' ? 1 : 3, $ty, $u[$i][1], $u[$i][2], $u[$i][3], $u[$i][4], '.', attr_str("$oid.$tag$num", $oid));
            $stat{"out_$ty"}++;
        }
    }
    if ($opt_codon) {
        for my $p (@startc) { $add->(4, 'start_codon', $p->[0], $p->[1], $p->[2], $p->[4], ($p->[3] =~ /^[012]$/ ? $p->[3] : '.'), attr_str("$oid.start_codon", $oid)); $stat{out_start_codon}++; }
        for my $p (@stopc)  { $add->(5, 'stop_codon',  $p->[0], $p->[1], $p->[2], $p->[4], ($p->[3] =~ /^[012]$/ ? $p->[3] : '.'), attr_str("$oid.stop_codon",  $oid)); $stat{out_stop_codon}++; }
    }
    else { $stat{codon_dropped} += @startc + @stopc; }

    @rows = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @rows;
    my $tline = join("\t", $seq, $src, $type, $ts, $te, '.', $strand, '.', attr_str($oid, $t->{gid}, $t->{name}, $t->{extras}));

    $coding ? $stat{coding_tx}++ : $stat{noncoding_tx}++;
    @exon > 1 ? $stat{multi_exon_tx}++ : $stat{single_exon_tx}++;
    $stat{tx_written}++;
    return { id => $oid, type => $type, coding => $coding, strand => $strand, start => $ts, end => $te,
             text => $tline . "\n" . join("\n", map { $_->[2] } @rows) . "\n" };
}

# ----------------------------------------------------------------------
# 3. 按基因组位置组织并输出
# ----------------------------------------------------------------------
my %out;    # gid => { seq, start, text }
for my $gid (@gene_order) {
    my $g = $gene{$gid};
    my @res;
    for my $tid (@{ $g->{tx} }) { my $r = build_tx($tid, $tx{$tid}); push @res, $r if $r; }
    @res = sort { $a->{start} <=> $b->{start} || $a->{id} cmp $b->{id} } @res;

    my ($gs, $ge);
    if (@res) {
        ($gs, $ge) = ($res[0]{start}, $res[0]{end});
        for (@res) { $gs = $_->{start} if $_->{start} < $gs; $ge = $_->{end} if $_->{end} > $ge; }
    }
    if ($g->{has_line}) {
        if (@res && ($g->{s} != $gs || $g->{e} != $ge)) {
            warnx('gene_range', "gene '$gid': gene 行坐标 $g->{s}-$g->{e} 与转录本范围 $gs-$ge 不一致，取并集");
            $stat{gene_range_mismatch}++;
        }
        $gs = $g->{s} if !defined $gs || $g->{s} < $gs;
        $ge = $g->{e} if !defined $ge || $g->{e} > $ge;
    }
    else { $stat{gene_line_inferred}++; }
    if (!defined $gs) { warnx('gene_empty', "gene '$gid' 没有可用的转录本和坐标，已跳过"); next; }
    $stat{gene_without_tx}++ unless @res;

    my $strand = $g->{has_line} ? $g->{strand} : (@res ? $res[0]{strand} : '.');
    my $gtype = $g->{ltype} // 'gene';
    if (!defined $g->{ltype} && @res) {
        my $bt = biotype($g->{extras});
        if ($bt =~ /pseudogene/ || !grep { $_->{type} ne 'pseudogenic_transcript' } @res) { $gtype = 'pseudogene'; }
        elsif (!grep { $_->{coding} || !$NCFAMILY{ $_->{type} } } @res) { $gtype = 'ncRNA_gene'; }
    }
    my $src = $g->{has_line} ? $g->{src} : ($tx{ $g->{tx}[0] }{src} // '.');
    my $seq = $g->{has_line} ? $g->{seq} : $tx{ $g->{tx}[0] }{seq};
    my $text = join("\t", $seq, $src, $gtype, $gs, $ge, '.', $strand, '.', attr_str($gid, undef, $g->{name}, $g->{extras})) . "\n";
    $text .= $_->{text} for @res;
    $text .= "###\n" if $opt_sep;
    $out{$gid} = { seq => $seq, start => $gs, text => $text };
    $stat{genes_written}++;
    $stat{"gene_type_$gtype"}++;
}

my $out_fh;
if (defined $opt_out) { open($out_fh, '>', $opt_out) or die "ERROR: cannot write $opt_out: $!\n"; }
else { $out_fh = \*STDOUT; }
print $out_fh "##gff-version 3\n";
for my $gid (sort { ($seq_order{ $out{$a}{seq} } // 1e9) <=> ($seq_order{ $out{$b}{seq} } // 1e9) || $out{$a}{start} <=> $out{$b}{start} || $a cmp $b } keys %out) {
    print $out_fh $out{$gid}{text};
}
close $out_fh if defined $opt_out;

# ----------------------------------------------------------------------
# 4. 统计信息（STDERR）
# ----------------------------------------------------------------------
sub st { return $stat{ $_[0] } // 0; }
my $E = \*STDERR;
print $E "\n================ gtfToGff3 转换统计 ================\n";
printf $E "%-38s %s\n", '输入 GTF', $opt_gtf;
printf $E "%-38s %d\n", '有效数据行', st('lines');
printf $E "%-38s %d\n", '格式错误/非法行(已跳过)', st('bad_line');
printf $E "%-38s %d\n", '缺少 gene_id/transcript_id 的行', st('line_no_id');
printf $E "%-38s %d\n", '起止坐标颠倒(已交换)', st('swapped');
if (%in_type) {
    print $E "输入特征类型计数:\n";
    printf $E "    %-34s %d\n", $_, $in_type{$_} for sort keys %in_type;
}
print $E "--- 输出 ---\n";
printf $E "%-38s %d\n", '基因数', st('genes_written');
printf $E "%-38s %d\n", '转录本数', st('tx_written');
printf $E "%-38s %d / %d\n", '  编码 / 非编码', st('coding_tx'), st('noncoding_tx');
printf $E "%-38s %d / %d\n", '  单外显子 / 多外显子', st('single_exon_tx'), st('multi_exon_tx');
printf $E "%-38s %.2f\n", '平均每基因转录本数', st('genes_written') ? st('tx_written') / st('genes_written') : 0;
if (%txtype_n) {
    print $E "转录本类型计数:\n";
    printf $E "    %-34s %d\n", $_, $txtype_n{$_} for sort keys %txtype_n;
}
for my $k (sort grep { /^gene_type_/ } keys %stat) { (my $n = $k) =~ s/^gene_type_//; printf $E "    %-34s %d\n", "基因类型 $n", $stat{$k}; }
printf $E "%-38s %d\n", 'exon 行', st('out_exon');
printf $E "%-38s %d\n", 'CDS 行', st('out_CDS');
printf $E "%-38s %d / %d\n", 'five_prime_UTR / three_prime_UTR 行', st('out_five_prime_UTR'), st('out_three_prime_UTR');
if ($opt_codon) { printf $E "%-38s %d / %d\n", 'start_codon / stop_codon 行', st('out_start_codon'), st('out_stop_codon'); }
else { printf $E "%-38s %d\n", '丢弃的 start/stop_codon 行', st('codon_dropped'); }
print $E "--- 修复/推断 ---\n";
printf $E "%-38s %d\n", '推断 gene 行', st('gene_line_inferred');
printf $E "%-38s %d\n", '推断 transcript 行', st('tx_line_inferred');
printf $E "%-38s %d\n", '无 exon，由 CDS/UTR 推断 exon', st('tx_exon_inferred');
printf $E "%-38s %d\n", '推断 UTR 的转录本', st('tx_utr_inferred');
printf $E "%-38s %d\n", 'CDS phase 被重算修正', st('phase_corrected');
printf $E "%-38s %d / %d\n", 'CDS 被延伸含终止子 / 被削去终止子', st('cds_extended'), st('cds_trimmed');
printf $E "%-38s %d\n", '转录本ID与基因ID冲突(已改名)', st('id_collision');
printf $E "%-38s %d\n", '缺少 gene_id(用 transcript_id 代替)', st('gene_id_missing');
print $E "--- 问题 ---\n";
printf $E "%-38s %d\n", '无子特征被跳过的转录本', st('tx_skipped_empty');
printf $E "%-38s %d\n", 'CDS 长度非3倍数', st('cds_not_mult3');
printf $E "%-38s %d\n", 'transcript 行范围与子特征不一致', st('tx_range_mismatch');
printf $E "%-38s %d\n", 'gene 行范围与转录本不一致', st('gene_range_mismatch');
printf $E "%-38s %d\n", '无转录本的 gene', st('gene_without_tx');
printf $E "%-38s %d\n", '非编码转录本中被丢弃的通用 UTR', st('utr_dropped_noncoding');
if (%ignored) {
    print $E "被忽略的特征类型(GFF3中无对应层级或可由其它特征推出):\n";
    printf $E "    %-34s %d\n", $_, $ignored{$_} for sort keys %ignored;
}
if (%wcount) {
    print $E "警告分类计数:\n";
    printf $E "    %-34s %d\n", $_, $wcount{$_} for sort keys %wcount;
}
else { print $E "警告: 无\n"; }
printf $E "运行耗时: %d 秒\n", time - $T0;
print $E "====================================================\n";

########################################################################
# 英文用法说明（放在文件最末尾，运行时通过DATA文件句柄读取）
########################################################################
__DATA__
Program: gtfToGff3.pl

Description: Convert gene models from GTF to GFF3 (the reverse of gff3ToGtf.pl), automatically rebuild the gene and transcript (mRNA etc.) hierarchy, ID and Parent, derive phase from exon and CDS, and report warnings and summary statistics to standard error.

Usage:
    perl gtfToGff3.pl -i file.gtf -o file.gff3 [other options]

    perl gtfToGff3.pl file.gtf > file.gff3

[Input and Output Options]
    -i, --gtf FILE
    Input GTF gene model file; .gz compressed files are supported; required; it can also be given as a positional argument.

    -o, --out FILE
    Output GFF3 file; default is standard output.

[Gene Model Conversion Options]
    --cds_stop MODE
    Control whether CDS includes the stop codon: include (default, GFF3 convention, CDS includes the stop codon and it is added from stop_codon lines when the input CDS lacks it), keep (keep the input unchanged), exclude (CDS excludes the stop codon and it is trimmed when the input CDS contains it).

    --keep_codon
    Keep start_codon and stop_codon and output them as children of the transcript; by default they are dropped (GFF3 usually does not record start and stop codons separately).

    --infer_utr
    Infer five_prime_UTR and three_prime_UTR from the difference between exon and CDS when the input has no UTR.

    --coding_type STR
    Type of coding transcripts (with CDS); default is mRNA.

    --noncoding_type STR
    Type of non-coding transcripts whose type cannot be determined from the GTF feature type or biotype; default is ncRNA.

    --no_extra_attr
    Do not keep attributes other than gene_id, transcript_id, gene_name and transcript_name (such as gene_biotype, tag and protein_id); only output ID, Parent and Name.

    --add_separator
    Output a GFF3 ### separator line after each gene to facilitate streaming processing.

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

