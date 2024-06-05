#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    Perl $0 genome.fasta > new_genome.fasta

    程序能对基因组序列进行整理：（1）对基因组序列按从长到短排序；（2）对序列重命名，或取消序列名紧跟空及其后后的字符信息；（3）对基因组序列全部修正为大写字符；（4）去除长度较短序列；（4）修正多次出现的序列名保留其序列，去除空行，每行输出指定字符长度等。
    程序在标准输出中给出结果FASTA信息，在标准错误输出中给出基因组大小。

    --sort_type <INT>    default: 1
    设置对序列的排序方式。默认值1，对输入的序列按从长到短进行排序；2，按输入序列名称的ASCII码排序，推荐结合--no_rename参数运行程序；3，和输入的FASTA文件序列一致。

    --no_rename
    设置不对序列进行重命名。若添加该参数，表示不会对序列名称进行重命名，参数--seq_prefix不会生效。

    --no_change_bp
    设置不对碱基进行修改。默认设置下，程序会将序列中除ATCG以外的其它字符变为碱基N，并将小写字符变为大写。添加该参数则不对序列进行修改。

    --seq_prefix <String>  default: "scaffold"
    设置--no_rename参数后，该参数失效。重命名后的序列名称以该指定的参数为前缀，后接逐一递增的数字编号，编号前用数字0补齐以使所有序列数字编号的字符数一致。

    --old_and_new <string>    default: None
    输出旧名和新名对照表。输出文件包含两列内容，第一列为旧名，第二列为新名。当设置了--no_rename时该参数不生效。

    --old_name_in_header
    添加该参数后，在输出FASTA序列头部信息序列名后增加一个空格再增加[formerName=xxx]信息用记录其曾用名。当设置了--no_rename时该参数不生效。

    --chromosome_mode    default: None
    添加该参数后，程序能识别序列名前缀为Chr（不区分大小写）的染色体名称。在没有设置--no_rename参数时，程序对这些序列名重命名为Chr开头，对其它序列重命名成--seq_prefix设置的前缀。当设置--sort_type参数值为3时，程序优先输出染色体级别的序列。染色体序列不会被--min_length参数设置的阈值过滤。

    --length_in_header
    添加该参数后，在输出FASTA序列头部信息后增加一个空格，再增加[length=xxx]信息用于记录序列长度。

    --added_info_to_header <string>    default: None
    使用该参数在FASTA序列头部信息后增加一个空格，再增加指定的信息。

    --min_length <INT>  default: 1000
    设置最短的序列长度。丢弃长度低于此阈值的序列。

    --line_length <INT>  default: 80
    设置输出的fasta文件中，序列在每行的最大字符长度。若该值 < 0，则表示不对序列进行换行处理。
    
USAGE
if (@ARGV==0){die $usage}

my ($sort_type, $no_rename, $no_change_bp, $seq_prefix, $old_and_new, $old_name_in_header, $length_in_header, $chromosome_mode, $added_info_to_header, $min_length, $line_length);
GetOptions(
    "sort_type:i" => \$sort_type,
    "no_rename" => \$no_rename,
    "no_change_bp" => \$no_change_bp,
    "seq_prefix:s" => \$seq_prefix,
    "old_and_new:s" => \$old_and_new,
    "old_name_in_header" => \$old_name_in_header,
    "length_in_header" => \$length_in_header,
    "added_info_to_header:s" => \$added_info_to_header,
    "min_length:i" => \$min_length,
    "line_length:i" => \$line_length,
    "chromosome_mode" => \$chromosome_mode,
);
$sort_type ||= 1;
$seq_prefix ||= "scaffold";
$min_length ||= 1000;
$line_length ||= 80;

open IN, @ARGV[0] or die $!;
my (%seq, $seq_id, @seq_id, %seq_id, %seq_length, $seq_num);
while (<IN>) {
    chomp;
    if (m/^>(\S+)/) {
        $seq_id = $1;
        if (exists $seq_id{$seq_id}) {
            my $num = $seq_id{$seq_id};
            $seq_id{$seq_id} ++;
            print STDERR "Warning: $seq_id appears $seq_id{$seq_id} times! forcibly rename this sequence id to $seq_id\_$num\n";
            $seq_id = "$seq_id\_$num";
        }
        $seq_id{$seq_id} ++;
        push @seq_id, $seq_id;
        $seq_num ++;
        print STDERR "正读取第 $seq_num 条序列，$seq_id .\r";
    }
    else {
        my $seq = $_;
        unless ($no_change_bp) {
            $_ = uc($_);
            @_ = split //, $_;
            $seq = "";
            foreach (@_) {
                if ($_ =~ m/[ATCGN]/) {
                    $seq .= $_;
                }
                else {
                    $seq .= "N";
                }
            }
        }
        $seq{$seq_id} .= $seq;
        $seq_length{$seq_id} += length($seq);
    }
}
close IN;
print STDERR "对基因组所有序列共 $seq_num 条，读取完毕。\n";

# 对基因组序列按从长到短进行排序
if ( $sort_type == 1 ) {
    @seq_id = sort {$seq_length{$b} <=> $seq_length{$a}} @seq_id;
}
elsif ( $sort_type == 2 ) {
    @seq_id = sort {$a cmp $b} @seq_id;
}
# 将染色体序列排到最前面
my @chr_id = grep {/^chr/i} @seq_id;
if ( $sort_type == 3 ) {
    my @new_seq_id = @chr_id;
    foreach (@seq_id) {
        push @new_seq_id, $_ unless m/^chr/i;
    }
    @seq_id = @new_seq_id;
}

# 清空重命名对应文件内容
if ($old_and_new) {
    open OUT, ">", $old_and_new or die "Can not create file $old_and_new, $!";
    close OUT;
}

# 分别记录染色体编号和其它序列编号
my ($number1, $number2, $genome_size) = (0, 0, 0);
foreach my $id (@seq_id) {
    if ($seq_length{$id} >= $min_length or $id =~ m/^chr/i) {
        $genome_size += $seq_length{$id};
        if ($id =~ m/^chr/i && $chromosome_mode) {
            $number1 ++;
        }
        else {
            $number2 ++;
        }

        my $seq_name = $id;
        if ($length_in_header) {
            $seq_name = "$seq_name [length=$seq_length{$id}]";
        }
        unless ($no_rename) {
            if ($id =~ m/^chr/i && $chromosome_mode) {
                $seq_name = "Chr" . "0" x ((length @chr_id) - (length $number1)) . $number1;
            }
            else {
                $seq_name = $seq_prefix . "0" x ((length @seq_id) - (length $number2)) . $number2;
            }
            if ($old_and_new) {
                open OUT, ">>", $old_and_new or die "Can not create file $old_and_new, $!";
                print OUT "$id\t$seq_name\n";
                close OUT;
            }
            if ($old_name_in_header) {
                $seq_name = "$seq_name [formerName=$id]";
            }
            if ($length_in_header) {
                $seq_name = "$seq_name [length=$seq_length{$id}]";
            }
        }
        my $seq = $seq{$id};
        if ($line_length > 0) {
            $seq =~ s/(\w{$line_length})/$1\n/g;
            $seq =~ s/\n$//;
        }
        if ($added_info_to_header) {
            $seq_name = "$seq_name $added_info_to_header";
        }
        print ">$seq_name\n$seq\n";
    }
    else {
        next;
    }
}
print STDERR $genome_size;
