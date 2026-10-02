#!/usr/bin/env perl
#===============================================================================
# check_internal_stop.pl
# 功能：检测 GFF3 中基因模型（CDS）内部出现提前终止密码子的基因数量
#
# 用法：
#   perl check_internal_stop.pl -gff genes.gff3 -fasta genome.fa [-table 1] [-out report.tsv]
#
# 参数：
#   -gff    GFF3 注释文件（必需）
#   -fasta  基因组 FASTA 文件（必需）
#   -table  遗传密码表：1=标准(默认), 2=脊椎动物线粒体, 4=支原体/螺旋体(TGA=W), 11=细菌/植物质体
#   -out    输出每个问题转录本的详细报告（TSV，可选）
#
# 判定规则：
#   1. 按 Parent 将 CDS 归组到转录本（mRNA 等），按链方向拼接 CDS，
#      并根据第一个 CDS 的 phase 去掉开头多余碱基
#   2. 翻译后，除最后一个密码子外，其余位置出现终止密码子即视为“内部终止”
#   3. 只要基因下任意一个转录本含内部终止，该基因即被计为异常基因
#===============================================================================
use strict;
use warnings;
use Getopt::Long;

my ($gff_file, $fasta_file, $out_file);
my $table = 1;
GetOptions(
    'gff=s'   => \$gff_file,
    'fasta=s' => \$fasta_file,
    'table=i' => \$table,
    'out=s'   => \$out_file,
) or die_usage();
die_usage() unless defined $gff_file && defined $fasta_file;

#---------------------------- 遗传密码表 --------------------------------------
my %codon_table = build_codon_table($table);

#---------------------------- 读取基因组 --------------------------------------
my %genome;
{
    open my $fh, '<', $fasta_file or die "无法打开 FASTA 文件 $fasta_file: $!\n";
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
    open my $fh, '<', $gff_file or die "无法打开 GFF3 文件 $gff_file: $!\n";
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
                    phase  => ($phase =~ /^[012]$/ ? $phase : 0),
                };
            }
        }
    }
    close $fh;
}

#---------------------------- 逐转录本检测 ------------------------------------
my $total_tx      = 0;
my $skipped_tx    = 0;
my $bad_tx        = 0;
my %all_genes;     # 所有含 CDS 的基因
my %bad_genes;     # 含内部终止的基因 => [转录本...]
my @report;

for my $tid (sort keys %cds_of) {
    my @segs = @{ $cds_of{$tid} };
    my $seqid  = $segs[0]{seqid};
    my $strand = $segs[0]{strand};

    unless (exists $genome{$seqid}) {
        warn "警告: 转录本 $tid 所在序列 $seqid 不在 FASTA 中，已跳过\n";
        $skipped_tx++;
        next;
    }

    $total_tx++;
    my $gene = find_gene($tid);
    $all_genes{$gene} = 1;

    @segs = sort { $a->{start} <=> $b->{start} } @segs;
    my $seq = '';
    for my $s (@segs) {
        $seq .= substr($genome{$seqid}, $s->{start} - 1, $s->{end} - $s->{start} + 1);
    }

    my $phase;
    if ($strand eq '-') {
        $seq   = revcomp($seq);
        $phase = $segs[-1]{phase};    # 转录方向上的第一个 CDS
    } else {
        $phase = $segs[0]{phase};
    }
    $seq = substr($seq, $phase) if $phase > 0 && length($seq) > $phase;

    my $ncodon = int(length($seq) / 3);
    my $remain = length($seq) % 3;
    my @internal;    # 内部终止密码子的位置（氨基酸序号，从 1 开始）
    for my $i (0 .. $ncodon - 1) {
        my $codon = substr($seq, $i * 3, 3);
        my $aa = $codon_table{$codon} // 'X';
        next unless $aa eq '*';
        push @internal, $i + 1 if $i < $ncodon - 1 || $remain > 0;
    }

    if (@internal) {
        $bad_tx++;
        push @{ $bad_genes{$gene} }, $tid;
        push @report, join("\t", $gene, $tid, $seqid, $strand, length($seq),
                           scalar(@internal), join(',', @internal));
    }
}

#---------------------------- 输出结果 ----------------------------------------
my $n_all = scalar keys %all_genes;
my $n_bad = scalar keys %bad_genes;

print "================ 内部终止密码子检测结果 ================\n";
printf "检测的遗传密码表           : %d\n", $table;
printf "含 CDS 的基因总数          : %d\n", $n_all;
printf "检测的转录本总数           : %d\n", $total_tx;
printf "跳过的转录本(缺少序列)     : %d\n", $skipped_tx;
printf "含内部终止密码子的转录本数 : %d\n", $bad_tx;
printf "含内部终止密码子的基因数   : %d (%.2f%%)\n",
       $n_bad, $n_all ? $n_bad / $n_all * 100 : 0;
print "========================================================\n";

if (defined $out_file) {
    open my $oh, '>', $out_file or die "无法写入 $out_file: $!\n";
    print $oh join("\t", qw(Gene Transcript Seqid Strand CDS_len_after_phase
                           N_internal_stops Stop_positions_aa)), "\n";
    print $oh "$_\n" for @report;
    close $oh;
    print "详细报告已写入: $out_file\n";
}

#============================== 子程序 ========================================
sub die_usage {
    die "用法: perl $0 -gff genes.gff3 -fasta genome.fa [-table 1] [-out report.tsv]\n";
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

sub build_codon_table {
    my ($t) = @_;
    my @b = qw(T C A G);
    my $aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
    my %tab;
    my $i = 0;
    for my $x (@b) { for my $y (@b) { for my $z (@b) {
        $tab{"$x$y$z"} = substr($aas, $i++, 1);
    } } }
    if    ($t == 1 || $t == 11) { }
    elsif ($t == 2) { $tab{TGA} = 'W'; $tab{ATA} = 'M'; $tab{AGA} = '*'; $tab{AGG} = '*'; }
    elsif ($t == 4) { $tab{TGA} = 'W'; }
    else  { die "暂不支持的遗传密码表: $t (支持 1, 2, 4, 11)\n"; }
    return %tab;
}
