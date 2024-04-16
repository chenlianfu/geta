#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    $0 [options] input.fasta > output.fasta

    程序用于对FASTA文件内容和格式进行修正。将序列名中的特殊字符进行替换，将内容中的不正确字符进行替换，过滤掉长度过短、不明确字符比例过大或再次出现相同ID的序列。

    本程序用于对FASTA文件进行格式修正。修正内容如下：
    （1）序列名称必须以 > 符号开头。若 > 符号未出现行首，则将 > 符号挪到下一行。
    （2）程序默认修改序列名称，以 > 符号后的非空字符为序列名称，再将名称中的特殊字符（除大小写字母、数字和下划线以外的其它字符）全部修改为下划线。程序也不保留序列名之后的其它字符。可以添加--no_change_header参数不对序列头部进行修改。
    （3）若某条序列字符长度低于--min_length参数设定的阈值（默认值为1），则舍弃该序列。
    （4）删除FASTA文件中的空行。
    （5）若序列是核酸序列，则将除ATCGNatcgn以外的其它序列全部修改为碱基N。程序能自动检测FASTA文件中序列的类型。通过检测FASTA文件前100条序列中ATCGNatcgn字符的比例是否超过90%，确定是否属于核酸序列。可以添加参数--seq_type DNA强行设置所有序列为核酸序列，则不进行自动检测。
    （6）若序列是蛋白序列，则先去除尾部的终止密码子符号*，再将除20种氨基酸以外的其它字符全部修改为X。可以添加参数--seq_type protein强行设置所有序列为蛋白序列，则不进行自动检测。
    （7）程序默认将所有序列字符修改为大写，可以添加--no_change_to_UC参数不进行大写转换。
    （8）程序默认在一行上输出整条序列，可以添加--line_length参数设置在多行上输出序列和每行的字符长度。
    （9）若某个序列ID在FASTA文件中出现多次时，仅保留第一次出现的序列。

    --no_change_header    default: None
    添加该参数后，则不对序列头部进行修改。

    --min_length <int>    default: 1
    若某条序列字符长度低于--min_length参数设定的阈值（默认值为1），则舍弃该序列。

    --max_unknown_character_ratio <float>    default：0.05
    设置单条序列中最大允许的不明确字符的比例。例如不明确氨基酸X占整条序列比例过高时，则不输出其序列。

    --seq_type <string>    default: None
    设置FASTA文件中所有序列的类型，其参数的值为DNA或protein。该参数的值不区分大小写，程序会将其参数值的字符统一转换为大写后进行识别。若不设置该参数的值，则程序会自动检测FASTA文件前前100条序列中析ATCGNatcgn字符的比例是否超过90%，从而判断FASTA文件的类型。

    --no_change_to_UC    default: None
    程序默认将所有序列字符修改为大写，可以添加--no_change_to_UC参数不进行大写转换。

    --line_length    default: None
    程序默认在一行上输出整条序列，可以添加--line_length参数设置在多行上输出序列和每行的字符长度。

    --help    default: None
    display this help and exit.

USAGE
if (@ARGV==0){die $usage}

my ($help_flag, $no_change_header, $min_length, $max_unknown_character_ratio, $seq_type, $no_change_to_UC, $line_length);
GetOptions(
    "no_change_header!" => \$no_change_header,
    "min_length:i" => \$min_length,
    "max_unknown_character_ratio:f" => \$max_unknown_character_ratio,
    "seq_type:s" => \$seq_type,
    "no_change_to_UC!" => \$no_change_to_UC,
    "line_length:i" => \$line_length,
    "help" => \$help_flag,
);
$min_length ||= 1;
$max_unknown_character_ratio ||= 0.05;

if ( $help_flag ) { die $usage }

# 自动检测FASTA文件序列的类型。
$seq_type = uc($seq_type);
if ( $seq_type =~ m/DNA/ ) {
    $seq_type = "DNA";
}
elsif ( $seq_type =~ m/PROTEIN/ ) {
    $seq_type = "protein";
}
else {
    open IN, $ARGV[0] or die "Can not open file $ARGV[0], $!";
    my $seq_num = 0;
    my ($seq100_numerator, $seq100_denominator);
    while ( <IN> ) {
        if ( m/^>/ ) {
        $seq_num ++;
            last if $seq_num > 100;
        }
        else {
            chomp;
            @_ = $_ =~ m/([ATCGNatcgn])/g;
            $seq100_numerator += @_;
            $seq100_denominator += length($_);
        }
    }
    close IN;
    #print "$seq100_numerator\t$seq100_denominator\n";
    my $bp_ratio = int($seq100_numerator * 10000 / $seq100_denominator + 0.5) / 100;
    $seq_type = "DNA";
    $seq_type = "protein" if $bp_ratio < 90;
    if ( $bp_ratio < 90 ) {
        print STDERR "未通过--seq_type参数设置FASTA文件中序列的类型为DNA或protein。程序通过检测FASTA文件前100条序列，检测到ATCGNatcgn字符的比例为$bp_ratio%，低于90%，因此认证FASTA文件序列类型为protein。\n";
    }
    else {
        print STDERR "未通过--seq_type参数设置FASTA文件中序列的类型为DNA或protein。程序通过检测FASTA文件前100条序列，检测到ATCGNatcgn字符的比例为$bp_ratio%，高于90%，因此认证FASTA文件序列类型为DNA。\n";
    }
}

open IN, $ARGV[0] or die "Can not open file $ARGV[0], $!";
my ($seq, $header, %header);
my ($seq_num, $num_output, $num_length_too_short, $num_unknown_character_ratio_too_high, $num_redundancy_ID) = (0, 0, 0, 0, 0);
# 读取文件数据时，当读取到一条新的序列时，输出上一条序列的数据。
while ( <IN> ) {
    if ( m/^>(.*)/ ) {
        $seq_num ++;
        &printOUT($header, $seq) if $seq; $seq = "";
        $header = $1;
    }
    elsif ( m/(.*?)>(.*)/ ) {
        $seq_num ++;
        $seq .= $1;
        &printOUT($header, $seq) if $seq; $seq = "";
        $header = $2;
    }
    else {
        chomp;
        $seq .= $_;
    }
}
close IN;
&printOUT($header, $seq) if $seq; $seq = "";
print STDERR "读取了 $seq_num 条序列，输出了 $num_output 条序列。\n未能输出的序列中，有 $num_length_too_short 条序列长度低于 $min_length; 有 $num_unknown_character_ratio_too_high 条序列包含的不明确字符比例超过 $max_unknown_character_ratio；有 $num_redundancy_ID 条序列由于重复ID被过滤。";

sub printOUT {
    my ($header, $seq) = @_;

    # 修改序列头部
    my $header_name;
    unless ( $no_change_header ) {
        $header =~ s/\s+.*//;
        $header_name = $header;
        $header =~ s/\W/_/g;
    }
    if ( exists $header{$header} ) {
        $header{$header} ++;
        print STDERR "Warning: 检测到序列 $header 在FASTA文件中出现的第 $header{$header} 次，不输出该序列。\n";
        $num_redundancy_ID ++;
        return;
    }
    $header{$header} ++;

    # 将序列字符变为大写
    unless ( $no_change_to_UC ) {
        $seq = uc($seq);
    }

    # 去除尾部字符 * 
    $seq =~ s/\*$//;

    # 将DNA序列中除ATCGNatcgn以外的其它序列全部修改为碱基N
    my $length = length($seq);
    my ($unkown_cha_ratio, @cha);
    if ( $seq_type eq "DNA" ) {
        if ( @cha = $seq =~ m/([^ATCGatcg])/g ) {
            my (%cha_num, @cha_num);
            foreach ( @cha ) {
                $cha_num{$_} ++;
            }
            foreach ( sort {$cha_num{$b} <=> $cha_num{$a}} keys %cha_num ) {
                push @cha_num, "$cha_num{$_}个$_";
            }
            my $cha_num = join "、", @cha_num;
            $unkown_cha_ratio = int(@cha / $length * 10000 + 0.5) / 100;
            print STDERR "Warning: 序列 $header_name 内部含有$cha_num，占整条序列比例为$unkown_cha_ratio%。\n";
            $seq =~ s/[^ATCGNatcgn]/N/g;
        }
    }

    # 将蛋白序列中除20种氨基酸以外的其它字符全部修改为X
    elsif ( $seq_type eq "protein" ) {
        # 20种氨基酸不包含B、J、O、U、X和Z。其中：X/Xaa/Unk表示任意氨基酸；B/Asx表示天冬氨酸（Asp、D）或天冬酰胺（Asn、N）；J/Xle可代表亮氨酸（Leu、L）或异亮氨酸（Ile、I）；Z/Glx可代表谷氨酸（Glu、E）或谷氨酰胺（Gln、Q）。
        # 有些FASTA蛋白文件中包含BJXZ字符，推荐进行修改，统一换为X，否则有些程序，例如exonerate等会运行出错。
        if ( @cha = $seq =~ m/([^ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy])/g ) {
            my (%cha_num, @cha_num);
            foreach ( @cha ) {
                $cha_num{$_} ++;
            }
            foreach ( sort {$cha_num{$b} <=> $cha_num{$a}} keys %cha_num ) {
                push @cha_num, "$cha_num{$_}个$_";
            }
            my $cha_num = join "、", @cha_num;
            $unkown_cha_ratio = int(@cha / $length * 10000 + 0.5) / 100;
            print STDERR "Warning: 序列 $header_name 内部含有$cha_num，占整条序列比例为$unkown_cha_ratio%。\n";
            $seq =~ s/[^ACDEFGHIKLMNPQRSTVWXYacdefghiklmnpqrstvwxy]/X/g;
        }
    }

    # 序列输出到多行；
    if ( $line_length >= 1 ) {
        $line_length = int($line_length);
        $seq =~ s/(\w{$line_length})/$1\n/g;
    }

    # 输出FASTA信息
    if ($length >= $min_length) {
        if ($unkown_cha_ratio <= $max_unknown_character_ratio * 100) {
            print ">$header\n$seq\n";
            $num_output ++;
        }
        else {
            $num_unknown_character_ratio_too_high ++;
        }
    }
    else {
        $num_length_too_short ++;
    }
}
