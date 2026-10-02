#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename qw(basename);

# ======================================================================
# 中文用法说明（English usage is at the end of this file / 英文用法说明见代码最尾部）
# ======================================================================
my $usage_cn = <<'USAGE_CN';
Usage:
    perl __SCRIPT__ [options] input.fasta > output.fasta

程序功能:
    程序用于对FASTA文件的内容和格式进行修正：将序列名中的特殊字符进行替换，将序列中不正确的字符进行替换，过滤掉长度过短、不明确字符比例过大或ID重复的序列。修正内容如下：
    （1）序列名称必须以 > 符号开头。若 > 符号未出现在行首，则将 > 符号及其后的内容挪到下一行作为新的序列名称。
    （2）程序默认修改序列名称：以 > 符号后的第一段非空字符作为序列名称，再将名称中的特殊字符（除大小写字母、数字和下划线以外的其它字符）全部修改为下划线，序列名之后的其它字符不保留。可以添加--no_change_header参数不对序列名称进行修改。若修改后序列名称为空，则将其命名为unnamed_序号。
    （3）若某条序列的字符长度低于--min_length参数设定的阈值（默认值为1），则舍弃该序列。
    （4）删除FASTA文件中的空行，并去除序列中的空格、制表符和Windows换行符（\r）。
    （5）若序列是核酸序列，则将除ATCGNatcgn以外的其它字符全部修改为碱基N。
    （6）若序列是蛋白序列，则先去除尾部的终止密码子符号*，再将除20种氨基酸以外的其它字符全部修改为X。
    （7）程序能自动检测FASTA文件中序列的类型：检测FASTA文件前100条序列（或累计检测到1000万个字符时停止）中ATCGNatcgn字符的比例是否达到90%，达到则认为是核酸序列，否则认为是蛋白序列。可以通过--seq_type参数强行设置所有序列的类型，此时不进行自动检测。
    （8）程序默认将所有序列字符修改为大写，可以添加--no_change_to_UC参数不进行大写转换。
    （9）程序默认在一行上输出整条序列，可以添加--line_length参数设置在多行上输出序列和每行的字符长度。
    （10）若某个序列ID在FASTA文件中出现多次，仅保留第一次出现的序列。
    （11）若单条序列中不明确字符（核酸序列中除ATCGatcg以外的字符，包括N；蛋白序列中除20种氨基酸以外的字符，包括X）所占比例超过--max_unknown_character_ratio参数设定的阈值，则舍弃该序列。

参数说明:
    --no_change_header    default: None
    添加该参数后，则不对序列名称进行修改（仅去除首尾空白字符）。

    --min_length <int>    default: 1
    若某条序列的字符长度低于--min_length参数设定的阈值（默认值为1），则舍弃该序列。

    --max_unknown_character_ratio <float>    default: 0.05
    设置单条序列中最大允许的不明确字符的比例，取值范围为0到1。例如不明确氨基酸X占整条序列的比例超过该值时，则不输出该序列。

    --seq_type <string>    default: None
    设置FASTA文件中所有序列的类型，其参数的值为DNA或protein。该参数的值不区分大小写。若不设置该参数，则程序会自动检测FASTA文件的序列类型。

    --no_change_to_UC    default: None
    程序默认将所有序列字符修改为大写，可以添加--no_change_to_UC参数不进行大写转换。

    --line_length <int>    default: None
    程序默认在一行上输出整条序列，可以添加--line_length参数设置在多行上输出序列和每行的字符长度。

    --quiet    default: None
    添加该参数后，不在标准错误输出中打印针对每条序列的Warning信息，仅输出最后的统计信息。

    --help    default: None
    Display the English usage and exit.

    --chinese_help    default: None
    显示中文用法说明并退出。

USAGE_CN

# ======================================================================
# 参数解析
# ======================================================================
my $script = basename($0);
$usage_cn =~ s/__SCRIPT__/$script/g;

my ($help_flag, $chinese_help_flag, $no_change_header, $no_change_to_UC, $quiet, $seq_type, $line_length);
my $min_length                  = 1;
my $max_unknown_character_ratio = 0.05;

GetOptions(
    "no_change_header"              => \$no_change_header,
    "min_length=i"                  => \$min_length,
    "max_unknown_character_ratio=f" => \$max_unknown_character_ratio,
    "seq_type=s"                    => \$seq_type,
    "no_change_to_UC"               => \$no_change_to_UC,
    "line_length=i"                 => \$line_length,
    "quiet"                         => \$quiet,
    "help"                          => \$help_flag,
    "chinese_help"                  => \$chinese_help_flag,
) or do { print STDERR english_usage(); exit 1; };

if ( $chinese_help_flag ) { print $usage_cn;         exit 0; }
if ( $help_flag )         { print english_usage();   exit 0; }
if ( @ARGV != 1 )         { print STDERR english_usage(); exit 1; }

my $input_file = $ARGV[0];
die "Error: 文件 $input_file 不存在或不可读。\n" unless -r $input_file;

$min_length = 1 if $min_length < 1;
die "Error: --max_unknown_character_ratio 的取值范围应为0到1。\n"
    if $max_unknown_character_ratio < 0 or $max_unknown_character_ratio > 1;
$line_length = 0 unless defined $line_length && $line_length >= 1;

# 确定序列类型：优先使用--seq_type参数，否则自动检测。
if ( defined $seq_type && length $seq_type ) {
    my $t = uc($seq_type);
    if    ( $t eq "DNA" )     { $seq_type = "DNA"; }
    elsif ( $t eq "PROTEIN" ) { $seq_type = "protein"; }
    else { die "Error: --seq_type 参数的值只能是 DNA 或 protein，而不是 $seq_type。\n"; }
}
else {
    $seq_type = detect_seq_type($input_file);
}

# ======================================================================
# 主程序：单次扫描，逐行读取，遇到新的序列名时输出上一条序列
# ======================================================================
my %seen;
my ($seq_num, $num_output, $num_length_too_short, $num_unknown_character_ratio_too_high, $num_redundancy_ID) = (0, 0, 0, 0, 0);
my ($header, $seq, $have_header) = ("", "", 0);

open my $IN, "<", $input_file or die "Can not open file $input_file, $!\n";
while ( my $line = <$IN> ) {
    $line =~ s/[\r\n]+\z//;
    my $pos = index($line, ">");
    if ( $pos >= 0 ) {
        # > 符号不在行首时，> 符号之前的内容属于上一条序列
        if ( $pos > 0 ) {
            my $before = substr($line, 0, $pos);
            $before =~ s/\s+//g;
            $seq .= $before;
        }
        if ( $have_header ) {
            process_record($header, \$seq);
        }
        elsif ( length $seq ) {
            warn_msg("Warning: 第一个 > 符号之前存在不属于任何序列的内容，已舍弃。\n");
        }
        $seq = "";
        $header = substr($line, $pos + 1);
        $have_header = 1;
        $seq_num ++;
    }
    else {
        $line =~ s/\s+//g;     # 去除空格、制表符等；空行因此被删除
        $seq .= $line if length $line;
    }
}
close $IN;
process_record($header, \$seq) if $have_header;

print STDERR "读取了 $seq_num 条序列，输出了 $num_output 条序列。\n未能输出的序列中，有 $num_length_too_short 条序列长度低于 $min_length；有 $num_unknown_character_ratio_too_high 条序列包含的不明确字符比例超过 $max_unknown_character_ratio；有 $num_redundancy_ID 条序列由于重复ID被过滤。\n";

# ======================================================================
# 子程序
# ======================================================================
sub warn_msg {
    print STDERR @_ unless $quiet;
}

# 自动检测FASTA文件序列的类型：检测前100条序列（或累计1000万个字符）中ATCGNatcgn字符的比例
sub detect_seq_type {
    my ($file) = @_;
    open my $in, "<", $file or die "Can not open file $file, $!\n";
    my ($n_rec, $n_atcgn, $n_total) = (0, 0, 0);
    while ( my $l = <$in> ) {
        if ( substr($l, 0, 1) eq ">" ) {
            last if ++ $n_rec > 100;
            next;
        }
        $l =~ tr/ \t\r\n//d;
        next unless length $l;
        $n_total += length $l;
        $n_atcgn += ($l =~ tr/ATCGNatcgn//);
        last if $n_total >= 10_000_000;
    }
    close $in;

    if ( $n_total == 0 ) {
        print STDERR "Warning: 未能在FASTA文件中检测到任何序列字符，默认认定序列类型为DNA。\n";
        return "DNA";
    }
    my $bp_ratio = int($n_atcgn * 10000 / $n_total + 0.5) / 100;
    if ( $bp_ratio < 90 ) {
        print STDERR "未通过--seq_type参数设置FASTA文件中序列的类型为DNA或protein。程序通过检测FASTA文件前100条序列，检测到ATCGNatcgn字符的比例为$bp_ratio%，低于90%，因此认定FASTA文件序列类型为protein。\n";
        return "protein";
    }
    print STDERR "未通过--seq_type参数设置FASTA文件中序列的类型为DNA或protein。程序通过检测FASTA文件前100条序列，检测到ATCGNatcgn字符的比例为$bp_ratio%，不低于90%，因此认定FASTA文件序列类型为DNA。\n";
    return "DNA";
}

# 处理并输出一条序列。$seq_ref 为序列字符串的引用，避免大序列的复制。
sub process_record {
    my ($header, $seq_ref) = @_;

    # 修改序列名称
    my $header_name;
    if ( $no_change_header ) {
        $header =~ s/^\s+//;
        $header =~ s/\s+\z//;
        $header_name = $header;
    }
    else {
        $header =~ s/^\s+//;
        $header =~ s/\s.*//s;
        $header_name = $header;
        $header =~ tr/A-Za-z0-9_/_/c;
    }
    if ( $header eq "" ) {
        $header = $header_name = "unnamed_$seq_num";
        warn_msg("Warning: 第 $seq_num 条序列的名称为空，已命名为 $header。\n");
    }

    # 重复ID只保留第一次出现的序列
    if ( $seen{$header} ++ ) {
        warn_msg("Warning: 检测到序列 $header_name 在FASTA文件中出现的第 $seen{$header} 次，不输出该序列。\n");
        $num_redundancy_ID ++;
        return;
    }

    # 将序列字符变为大写
    $$seq_ref =~ tr/a-z/A-Z/ unless $no_change_to_UC;

    # 去除尾部终止密码子符号 *
    $$seq_ref =~ s/\*+\z//;

    # 长度过滤（先过滤可避免对被舍弃序列做无用的字符处理）
    my $length = length($$seq_ref);
    if ( $length < $min_length ) {
        warn_msg("Warning: 序列 $header_name 的长度为 $length，低于 $min_length，不输出该序列。\n");
        $num_length_too_short ++;
        return;
    }

    my ($unknown_num, $illegal_num, $illegal_regex);
    if ( $seq_type eq "DNA" ) {
        # 不明确字符：除ATCGatcg以外的字符（包括N）
        $unknown_num = ($$seq_ref =~ tr/ATCGatcg//c);
        # 非法字符：除ATCGNatcgn以外的字符，将被替换为N
        $illegal_num = ($$seq_ref =~ tr/ATCGNatcgn//c);
        if ( $illegal_num ) {
            report_illegal($header_name, $seq_ref, qr/[^ATCGNatcgn]/, $illegal_num, $length) unless $quiet;
            $$seq_ref =~ tr/ATCGNatcgn/N/c;
        }
    }
    else {
        # 20种氨基酸不包含B、J、O、U、X和Z。其中：X/Xaa/Unk表示任意氨基酸；B/Asx表示Asp或Asn；J/Xle表示Leu或Ile；Z/Glx表示Glu或Gln。
        # 有些FASTA蛋白文件中包含BJXZ字符，推荐统一换为X，否则有些程序（例如exonerate等）会运行出错。
        $unknown_num = ($$seq_ref =~ tr/ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy//c);
        $illegal_num = ($$seq_ref =~ tr/ACDEFGHIKLMNPQRSTVWXYacdefghiklmnpqrstvwxy//c);
        if ( $illegal_num ) {
            report_illegal($header_name, $seq_ref, qr/[^ACDEFGHIKLMNPQRSTVWXYacdefghiklmnpqrstvwxy]/, $illegal_num, $length) unless $quiet;
            $$seq_ref =~ tr/ACDEFGHIKLMNPQRSTVWXYacdefghiklmnpqrstvwxy/X/c;
        }
    }

    # 不明确字符比例过滤
    if ( $unknown_num > $max_unknown_character_ratio * $length ) {
        my $ratio = sprintf("%.2f", $unknown_num * 100 / $length);
        warn_msg("Warning: 序列 $header_name 的不明确字符比例为$ratio%，超过 " . ($max_unknown_character_ratio * 100) . "%，不输出该序列。\n");
        $num_unknown_character_ratio_too_high ++;
        return;
    }

    # 输出FASTA信息
    if ( $line_length ) {
        print ">$header\n", join("\n", unpack("(a$line_length)*", $$seq_ref)), "\n";
    }
    else {
        print ">$header\n$$seq_ref\n";
    }
    $num_output ++;
}

# 输出被替换的非法字符的统计信息
sub report_illegal {
    my ($name, $seq_ref, $regex, $illegal_num, $length) = @_;
    my %cha_num;
    $cha_num{$_} ++ foreach ( $$seq_ref =~ /($regex)/g );
    my $detail = join "、", map { "$cha_num{$_}个$_" } sort { $cha_num{$b} <=> $cha_num{$a} || $a cmp $b } keys %cha_num;
    my $ratio = sprintf("%.2f", $illegal_num * 100 / $length);
    print STDERR "Warning: 序列 $name 内部含有$detail，占整条序列比例为$ratio%，已被替换。\n";
}

# 读取位于代码最尾部的英文用法说明
sub english_usage {
    my $text = do { local $/; <DATA> };
    $text =~ s/__SCRIPT__/$script/g;
    return $text;
}

# ======================================================================
# English usage (kept at the very end of the code)
# ======================================================================
__DATA__
Usage:
    perl __SCRIPT__ [options] input.fasta > output.fasta

Description:
    This program corrects the content and format of a FASTA file: it replaces special characters in sequence names, replaces incorrect characters in sequences, and filters out sequences that are too short, contain too high a proportion of ambiguous characters, or have a duplicated ID. The corrections are as follows:
    (1) A sequence name must start with the ">" symbol. If ">" does not appear at the beginning of a line, the ">" and the content after it are moved to the next line as a new sequence name.
    (2) By default the sequence name is modified: the first non-blank string after ">" is taken as the sequence name, all special characters in it (anything other than letters, digits and underscore) are replaced with underscores, and everything after the name is discarded. Use --no_change_header to leave sequence names unmodified. If the name becomes empty, it is named unnamed_<index>.
    (3) A sequence whose length is below the threshold set by --min_length (default 1) is discarded.
    (4) Blank lines are removed, and spaces, tabs and Windows line endings (\r) inside sequences are removed.
    (5) For nucleotide sequences, all characters other than ATCGNatcgn are replaced with the base N.
    (6) For protein sequences, the trailing stop-codon symbol * is removed first, and then all characters other than the 20 amino acids are replaced with X.
    (7) The sequence type is detected automatically: the program checks whether ATCGNatcgn characters make up at least 90% of the first 100 sequences (or of the first 10 million characters, whichever comes first); if so the file is treated as nucleotide, otherwise as protein. Use --seq_type to force the type of all sequences, in which case no detection is done.
    (8) By default all sequence characters are converted to upper case. Use --no_change_to_UC to disable the conversion.
    (9) By default each sequence is written on a single line. Use --line_length to write sequences on multiple lines with the given number of characters per line.
    (10) If a sequence ID appears more than once in the FASTA file, only the first occurrence is kept.
    (11) A sequence is discarded if the proportion of ambiguous characters (in nucleotide sequences, anything other than ATCGatcg, including N; in protein sequences, anything other than the 20 amino acids, including X) exceeds the threshold set by --max_unknown_character_ratio.

Options:
    --no_change_header    default: None
    Do not modify sequence names (only leading and trailing whitespace is removed).

    --min_length <int>    default: 1
    Discard sequences whose length is below this threshold (default 1).

    --max_unknown_character_ratio <float>    default: 0.05
    Maximum proportion of ambiguous characters allowed in a single sequence, in the range 0 to 1. For example, a sequence whose proportion of the ambiguous amino acid X exceeds this value is not output.

    --seq_type <string>    default: None
    Set the type of all sequences in the FASTA file, either DNA or protein (case-insensitive). If not set, the type is detected automatically.

    --no_change_to_UC    default: None
    Do not convert sequence characters to upper case (they are converted by default).

    --line_length <int>    default: None
    Write each sequence on multiple lines with this many characters per line (default: the whole sequence on one line).

    --quiet    default: None
    Do not print per-sequence Warning messages to STDERR; only the final summary is printed.

    --help    default: None
    Display the English usage and exit.

    --chinese_help    default: None
    显示中文用法说明并退出。
