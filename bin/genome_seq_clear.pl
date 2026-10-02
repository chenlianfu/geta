#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw/max/;

my $usage = <<USAGE;
Usage:
    perl $0 [options] genome.fasta > new_genome.fasta

    程序能对基因组序列进行整理：（1）对基因组序列按从长到短排序；（2）对序列重命名，或取消序列名紧跟空格及其后字符信息；（3）对基因组序列全部修正为大写字符；（4）去除长度较短序列；（5）修正多次出现的序列名保留其序列，去除空行，每行输出指定字符长度等。
    程序在标准输出中给出结果 FASTA 信息，在标准错误输出中给出基因组大小。

    --sort_type <INT>    default: 1
    设置对序列的排序方式。1：对输入的序列按从长到短进行排序；2：按输入序列名称的ASCII 码排序，推荐结合 --no_rename 参数运行程序；3：和输入的 FASTA 文件序列顺序一致。

    --no_rename
    设置不对序列进行重命名。若添加该参数，表示不会对序列名称进行重命名，参数--seq_prefix、--old_and_new、--old_name_in_header 均不会生效。

    --no_change_bp
    设置不对碱基进行修改。默认设置下，程序会将序列中除 ATCGN 以外的其它字符变为碱基 N，并将小写字符变为大写。添加该参数则不对序列进行修改。

    --seq_prefix <String>  default: "scaffold"
    设置 --no_rename 参数后，该参数失效。重命名后的序列名称以该指定的参数为前缀，后接逐一递增的数字编号，编号前用数字 0 补齐以使所有序列数字编号的字符数一致。

    --old_and_new <string>    default: None
    输出旧名和新名对照表。输出文件包含两列内容，第一列为旧名，第二列为新名。当设置了 --no_rename 时该参数不生效。

    --old_name_in_header
    添加该参数后，在输出 FASTA 序列头部信息序列名后增加一个空格再增加 [formerName=xxx] 信息用记录其曾用名。当设置了 --no_rename 时该参数不生效。

    --chromosome_mode    default: None
    添加该参数后，程序能识别序列名前缀为 Chr（不区分大小写）的染色体名称，并对其进行以下特殊处理：（1）在没有设置 --no_rename 参数时，将这些序列重命名为 Chr 开头，其它序列重命名成 --seq_prefix 设置的前缀；（2）当 --sort_type 设置为 3时，优先输出染色体级别的序列；（3）染色体序列不会被 --min_length 参数设置的阈值过滤。不设置该参数时，名称以 chr 开头的序列不会得到任何特殊处理。

    --length_in_header
    添加该参数后，在输出 FASTA 序列头部信息后增加一个空格，再增加 [length=xxx] 信息用于记录序列长度。

    --added_info_to_header <string>    default: None
    使用该参数在 FASTA 序列头部信息后增加一个空格，再增加指定的信息。

    --min_length <INT>  default: 1000
    设置最短的序列长度。丢弃长度低于此阈值的序列（设置了 --chromosome_mode 时，染色体级别的序列不受此参数限制）。

    --line_length <INT>  default: 80
    设置输出的 fasta 文件中，序列在每行的最大字符长度。若该值 <= 0，则表示不对序列进行换行处理。

    --help
    显示本帮助信息并退出。

USAGE
if ( @ARGV == 0 ) { die $usage }

my ( $sort_type, $no_rename, $no_change_bp, $seq_prefix, $old_and_new, $old_name_in_header,
     $length_in_header, $chromosome_mode, $added_info_to_header, $min_length, $line_length, $help_flag );

GetOptions(
    "sort_type:i"            => \$sort_type,
    "no_rename"              => \$no_rename,
    "no_change_bp"           => \$no_change_bp,
    "seq_prefix:s"           => \$seq_prefix,
    "old_and_new:s"          => \$old_and_new,
    "old_name_in_header"     => \$old_name_in_header,
    "length_in_header"       => \$length_in_header,
    "added_info_to_header:s" => \$added_info_to_header,
    "min_length:i"           => \$min_length,
    "line_length:i"          => \$line_length,
    "chromosome_mode"        => \$chromosome_mode,
    "help"                   => \$help_flag,
) or die $usage;
die $usage if $help_flag;

$sort_type   = 1          unless defined $sort_type;
$seq_prefix  = "scaffold" unless defined $seq_prefix;
$min_length  = 1000       unless defined $min_length;
$line_length = 80         unless defined $line_length;

if ( $sort_type != 1 && $sort_type != 2 && $sort_type != 3 ) {
    warn "Warning: --sort_type 的值应为 1、2 或 3，收到的是 \"$sort_type\"，将按默认值 1（从长到短排序）处理。\n";
    $sort_type = 1;
}

my $genome_file = $ARGV[0];
die "Error: 未指定基因组 fasta 文件。\n$usage" unless defined $genome_file;

# ---------------------------------------------------------------------------
# 读取基因组序列
# ---------------------------------------------------------------------------
open my $IN, "<", $genome_file or die "Error: Can not open file $genome_file: $!\n";

my ( %seq, $seq_id, @seq_id, %seq_id_count, %seq_length, $seq_num );
while (<$IN>) {
    chomp;
    if (m/^>(\S+)/) {
        $seq_id = $1;

        if ( exists $seq_id_count{$seq_id} ) {
            # 保证重命名后的 ID 一定唯一，避免与文件中另一条真实同名序列冲突
            my $orig_id = $seq_id;
            my $suffix  = $seq_id_count{$orig_id};
            my $new_id  = "${orig_id}_$suffix";
            while ( exists $seq_id_count{$new_id} ) {
                $suffix++;
                $new_id = "${orig_id}_$suffix";
            }
            $seq_id_count{$orig_id}++;
            print STDERR "Warning: $orig_id appears $seq_id_count{$orig_id} times! forcibly rename this sequence id to $new_id\n";
            $seq_id = $new_id;
        }
        $seq_id_count{$seq_id}++;
        push @seq_id, $seq_id;
        $seq_num++;
        print STDERR "正读取第 $seq_num 条序列，$seq_id .\r";
    }
    elsif ( defined $seq_id && length($_) ) {
        my $seq = $_;
        unless ($no_change_bp) {
            $seq = uc($seq);
            $seq =~ s/[^ATCGN]/N/g;
        }
        $seq{$seq_id} .= $seq;
        $seq_length{$seq_id} += length($seq);
    }
    # 若文件首行不是以 > 开头（非法 FASTA），直接忽略该行，不中止程序
}
close $IN;
print STDERR "\n对基因组所有序列共 $seq_num 条，读取完毕。\n";

# ---------------------------------------------------------------------------
# 排序
# ---------------------------------------------------------------------------
if ( $sort_type == 1 ) {
    @seq_id = sort { $seq_length{$b} <=> $seq_length{$a} } @seq_id;
}
elsif ( $sort_type == 2 ) {
    @seq_id = sort { $a cmp $b } @seq_id;
}
# sort_type == 3：保持读入顺序不变

# 仅在设置了 --chromosome_mode 时，才识别 Chr 前缀并进行相关特殊处理
my @chr_id = $chromosome_mode ? ( grep { /^chr/i } @seq_id ) : ();
if ( $sort_type == 3 && $chromosome_mode && @chr_id ) {
    my %is_chr = map { $_ => 1 } @chr_id;
    @seq_id = ( @chr_id, grep { !$is_chr{$_} } @seq_id );
}

# ---------------------------------------------------------------------------
# 预先统计"最终会被编号输出"的染色体数 / scaffold 数，用于确定补零宽度
# （必须和下面主循环里判断是否保留、是否按染色体编号的逻辑完全一致）
# ---------------------------------------------------------------------------
my %is_chr_seq;
my ( $chr_kept_count, $scaffold_kept_count ) = ( 0, 0 );
foreach my $id (@seq_id) {
    my $chr_flag = ( $chromosome_mode && $id =~ m/^chr/i ) ? 1 : 0;
    my $kept     = $chr_flag || ( $seq_length{$id} >= $min_length );
    next unless $kept;
    if ($chr_flag) { $is_chr_seq{$id} = 1; $chr_kept_count++; }
    else            { $scaffold_kept_count++; }
}
my $chr_digits      = length( $chr_kept_count      || 1 );
my $scaffold_digits = length( $scaffold_kept_count || 1 );

# ---------------------------------------------------------------------------
# 仅在确实会生效时才创建/清空 --old_and_new 对照表文件，且全程只打开一次
# ---------------------------------------------------------------------------
my $OLD_AND_NEW_FH;
if ( $old_and_new && !$no_rename ) {
    open $OLD_AND_NEW_FH, ">", $old_and_new or die "Error: Can not create file $old_and_new: $!\n";
}

# ---------------------------------------------------------------------------
# 主输出循环
# ---------------------------------------------------------------------------
my ( $number1, $number2, $genome_size ) = ( 0, 0, 0 );
foreach my $id (@seq_id) {
    my $chr_flag = $is_chr_seq{$id} ? 1 : 0;
    next unless $chr_flag || $seq_length{$id} >= $min_length;

    $genome_size += $seq_length{$id};
    my $seq_name = $id;

    unless ($no_rename) {
        if ($chr_flag) {
            $number1++;
            $seq_name = "Chr" . ( "0" x max( 0, $chr_digits - length($number1) ) ) . $number1;
        }
        else {
            $number2++;
            $seq_name = $seq_prefix . ( "0" x max( 0, $scaffold_digits - length($number2) ) ) . $number2;
        }
        print $OLD_AND_NEW_FH "$id\t$seq_name\n" if $OLD_AND_NEW_FH;
        $seq_name = "$seq_name [formerName=$id]" if $old_name_in_header;
    }

    $seq_name = "$seq_name [length=$seq_length{$id}]" if $length_in_header;
    $seq_name = "$seq_name $added_info_to_header"      if $added_info_to_header;

    my $seq = $seq{$id};
    $seq = '' unless defined $seq;
    if ( $line_length > 0 ) {
        $seq =~ s/(.{$line_length})/$1\n/g;
        $seq =~ s/\n$//;
    }
    print ">$seq_name\n$seq\n";
}
close $OLD_AND_NEW_FH if $OLD_AND_NEW_FH;

print STDERR "$genome_size\n";
