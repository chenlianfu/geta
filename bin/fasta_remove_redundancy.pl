#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Digest::MD5 qw(md5);
use File::Basename qw(basename);

########################################################################
# 中文用法说明（放在代码最前面；英文用法说明放在文件最末尾的__DATA__部分）
########################################################################
my $usage_cn = <<'USAGE';
用法:
    {PROG} [options] seq1.fasta seq2.fasta ... > out.fasta

    程序依次读取一个或多个FASTA文件中的序列，剔除冗余序列并输出结果。即当两条序列完全相同时，只保留最先读取到的那一条。输出序列名只取原序列名中“>”后第一个空白符之前的部分。标准输出为FASTA格式的结果，标准错误输出为各文件及合计的统计信息（读取、保留、删除冗余、空序列的数量）。输入文件可以是普通文本文件或.gz压缩文件，文件名为“-”时从标准输入读取。

    --label <string>    默认值: 无
    设置输入FASTA文件的标签。参数值是以逗号分隔的字符串列表，依次对应其后输入的各个FASTA文件。设置该参数后，标签信息会被直接追加到原序列名的尾部（不会自动添加分隔符，例如标签写成“_SpeciesSorghumBicolor”），用于区分不同物种的序列。若标签数量少于FASTA文件数量，后面的文件不添加标签；若标签数量多于FASTA文件数量，多出的标签不生效，并给出警告。

    --ignore_case    默认值: 无（区分大小写）
    设置该参数后，比较序列是否相同时忽略大小写，例如ACGT与acgt视为同一条序列。

    --line_width <int>    默认值: 0
    设置输出序列每行的碱基或氨基酸个数。默认值为0，表示每条序列输出为一行，不换行。

    --help    默认值: 无
    Display the English usage and exit.

    --chinese_help    默认值: 无
    显示中文用法说明并退出。

USAGE

########################################################################
# 全局变量
########################################################################
my ($label, $ignore_case, $line_width, $help_flag, $chinese_help_flag);
$ignore_case = 0;
$line_width  = 0;

my %seen;          # 已输出序列的MD5摘要
my %names;         # 已输出的序列名
my $dup_names = 0; # 输出序列名重复的数量

########################################################################
# 函数
########################################################################
sub usage_en {
    local $/ = undef;
    my $text = <DATA>;
    $text = "" unless defined $text;
    my $prog = basename($0);
    $text =~ s/\{PROG\}/$prog/g;
    return $text;
}

sub usage_cn {
    my $text = $usage_cn;
    my $prog = basename($0);
    $text =~ s/\{PROG\}/$prog/g;
    return $text;
}

sub open_input {
    my ($file) = @_;
    my $fh;
    if ( $file eq '-' ) {
        return \*STDIN;
    }
    die "Error: Can not find file $file\n" unless -e $file;
    die "Error: $file is a directory, not a file\n" if -d $file;
    if ( $file =~ m/\.gz$/i ) {
        open( $fh, '-|', 'gzip', '-dc', $file ) or die "Error: Can not open gzip file $file, $!\n";
    }
    else {
        open( $fh, '<', $file ) or die "Error: Can not open file $file, $!\n";
    }
    return $fh;
}

sub format_seq {
    my ($seq) = @_;
    return $seq if $line_width <= 0;
    my @lines = $seq =~ m/(.{1,$line_width})/g;
    return join( "\n", @lines );
}

# 输出一条序列（若不冗余）。
sub flush_record {
    my ( $id, $seq, $st ) = @_;
    return unless defined $id;
    $st->{total}++;
    if ( $seq eq "" ) {
        $st->{empty}++;
        return;
    }
    my $key = md5( $ignore_case ? uc($seq) : $seq );
    if ( exists $seen{$key} ) {
        $st->{removed}++;
        return;
    }
    $seen{$key} = 1;
    $st->{kept}++;
    $dup_names++ if $names{$id}++;
    print ">$id\n", format_seq($seq), "\n";
}

########################################################################
# 主程序
########################################################################
# 直接运行（无参数）时，给出英文用法说明
if ( @ARGV == 0 ) {
    print STDERR usage_en();
    exit 1;
}

GetOptions(
    "label:s"      => \$label,
    "ignore_case!" => \$ignore_case,
    "line_width=i" => \$line_width,
    "help"         => \$help_flag,
    "chinese_help" => \$chinese_help_flag,
) or do {
    print STDERR "Error: invalid option(s).\n\n", usage_en();
    exit 1;
};

if ($help_flag)         { print usage_en(); exit 0; }
if ($chinese_help_flag) { print usage_cn(); exit 0; }

if ( @ARGV == 0 ) {
    print STDERR "Error: no input FASTA file.\n\n", usage_en();
    exit 1;
}
die "Error: --line_width must be a non-negative integer.\n" if $line_width < 0;

# 分析标签
my @labels;
if ( defined $label && $label ne "" ) {
    @labels = map { s/^\s+|\s+$//gr } split( /,/, $label, -1 );
    print STDERR "Warning: ", scalar(@labels), " labels were given but only ", scalar(@ARGV),
        " input files; the extra labels are ignored.\n" if @labels > @ARGV;
}

my %total = ( total => 0, kept => 0, removed => 0, empty => 0 );

foreach my $file (@ARGV) {
    my $tag = @labels ? shift @labels : "";
    my $fh  = open_input($file);
    my %st  = ( total => 0, kept => 0, removed => 0, empty => 0 );
    my ( $id, $seq ) = ( undef, "" );
    my ( $unnamed, $warned_before_header ) = ( 0, 0 );

    while ( my $line = <$fh> ) {
        $line =~ s/[\r\n]+\z//;
        if ( $line =~ m/^>/ ) {
            flush_record( $id, $seq, \%st );
            ($id) = $line =~ m/^>\s*(\S*)/;
            if ( $id eq "" ) {
                $id = "unnamed_" . ( ++$unnamed );
                print STDERR "Warning: a sequence without a name was found in $file, renamed to $id.\n";
            }
            $id .= $tag;
            $seq = "";
        }
        elsif ( $line =~ m/\S/ ) {
            if ( !defined $id ) {
                print STDERR "Warning: sequence lines before the first header in $file are ignored.\n"
                    unless $warned_before_header++;
                next;
            }
            $line =~ s/\s+//g;
            $seq .= $line;
        }
    }
    flush_record( $id, $seq, \%st );
    close $fh unless $file eq '-';

    print STDERR "Warning: no FASTA sequence was found in $file.\n" if $st{total} == 0;
    print STDERR "File $file: read $st{total}, kept $st{kept}, redundant $st{removed}, empty $st{empty}.\n";
    $total{$_} += $st{$_} for keys %st;
}

print STDERR "There are $total{removed} redundant sequences were deleted.\n";
print STDERR "Total: read $total{total}, kept $total{kept}, redundant $total{removed}, empty (skipped) $total{empty}.\n";
print STDERR "Warning: $dup_names output sequences share a name with an earlier output sequence; consider using --label.\n"
    if $dup_names;

########################################################################
# 英文用法说明（放在文件最末尾，运行时通过DATA文件句柄读取）
########################################################################
__DATA__
Usage:
    {PROG} [options] seq1.fasta seq2.fasta ... > out.fasta

    The program reads the sequences in one or more FASTA files sequentially, eliminates any redundant sequences, and reports the result. That is, when two sequences are identical, only the sequence read first is saved. The output sequence name keeps only the part of the original name between ">" and the first whitespace. The standard output contains the FASTA result, while the standard error output provides statistics for each file and in total (numbers of sequences read, kept, removed as redundant, and empty). Input files can be plain text or .gz compressed files; a file name of "-" means reading from standard input.

    --label <string>    default: None
    Set the labels for the input FASTA files. The value of this parameter is a list of comma-separated strings that correspond to the FASTA files entered afterward. After this parameter is introduced, the label is appended directly to the end of the original sequence name (no separator is added automatically, e.g. use "_SpeciesSorghumBicolor") to help distinguish sequences from different species. If there are fewer labels than FASTA files, the remaining files get no label; if there are more labels than files, the extra labels do not take effect and a warning is printed.

    --ignore_case    default: off (case-sensitive)
    When set, case is ignored when comparing sequences, e.g. ACGT and acgt are treated as the same sequence.

    --line_width <int>    default: 0
    Set the number of residues per line in the output. The default 0 means each sequence is written on a single line without wrapping.

    --help    default: None
    Display this help (English) and exit.

    --chinese_help    default: None
    显示中文用法说明并退出。

