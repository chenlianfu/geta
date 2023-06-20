#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    perl $0 --out_prefix out genome.fasta input1.gff3 [input2.gff3 ...] > statitics.txt

    本程序用于根据GFF3的内容信息和基因组序列，转换出对应的序列信息，并对各种类型的Features进行统计。程序能对各种类型的Feature都进行序列转换，每种类型的Feature各得到一个Fasta文件，gene类型的Feature能得到CDS、cDNA和Protein三个Fasta文件。

    程序支持输入多个GFF3文件，并根据其中的Feature ID输出序列。所以输入的GFF3文件第9列一定得要有ID信息。若一个文件中同一个ID出现多次，则仅使用其ID后出现的数据信息；若多个文件中出现相同的ID，则使用输入文件靠最前的GFF3文件中的数据信息。程序最终输出序列信息时，按照输入GFF3文件先后顺序和GFF3文件内容中出现的ID先后顺序输出序列。

    若GFF3文件中包含编码基因信息，则对这些编码基因的CDS长度、基因的cDNA长度、基因的intron长度、基因的gene长度、基因的CDS个数、基因的exon个数、基因的intron个数、单个CDS长度、单个exon长度、单个intron长度和基因间区长度进行了统计，并将结果输入到out.codingGeneModels.stats文件中。

    --out_prefix <string>    default: out
    设置程序输出的序列文件前缀。程序根据GFF3文件的Feature Name生成对应的Fasta文件，out.FeatureName.fasta。若Feature Name为gene，则额外生成out.CDS.fasta，out.cDNA.fasta和out.protein.fasta文件。若有编码基因信息，即gene中有CDS feature，则额外生成out.codingGeneModels.stats文件。

    --only_gene_sequences    default: None
    添加该参数后，仅输出GFF3文件中属于gene类型的序列。

    --only_coding_gene_sequences    default: None
    添加该参数后，仅输出GFF3文件中编码基因的序列信息。

    --only_first_isoform    default: None
    添加该参数后，若一个基因有多个可变剪接，则仅选择在GFF3文件中出现的第一个isoform进行统计和序列输出。和--only_longest_isoform参数只能有一个生效。当两个参数同时设置时，本参数是有效参数。来自GETA软件预测的基因模型，一般第一个isoform的表达量占比最大。

    --only_longest_isoform    default: None
    添加该参数后，若一个基因有多个可变剪接，则仅选择其CDS最长或exon最长的isoform进行统计和序列输出。

    --sort_isoforms    default: None
    添加该参数后，对基因模型的多个可变剪接的转录本进行排序后再输出序列。优先按CDS长度从长到短，然后按cDNA长度从长到短，最后按ID的ASCII编码从小到大进行排序。程序默认输出所有的可变剪接序列，按其在GFF3文件中出现的顺序输出序列。
    
    --genetic_code <int>    default: 1
    设置遗传密码。该参数对应的值请参考NCBI Genetic Codes: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi。该参数在将CDS转换为Protein序列时生效。

USAGE
if (@ARGV==0){die $usage}

my ($help_flag, $out_prefix, $only_gene_sequences, $only_coding_gene_sequences, $only_longest_isoform, $only_first_isoform, $sort_isoform, $genetic_code);
GetOptions(
    "help" => \$help_flag,
    "out_prefix:s" => \$out_prefix,
    "only_gene_sequences" => \$only_gene_sequences,
    "only_coding_gene_sequences" => \$only_coding_gene_sequences,
    "only_longest_isoform" => \$only_longest_isoform,
    "only_first_isoform" => \$only_first_isoform,
    "sort_isoforms" => \$sort_isoform,
    "genetic_code:i" => \$genetic_code,
);
$out_prefix ||= "out";

$genetic_code ||= 1;
my (%genetic_code, %start_codon, %stop_codon);
@_ = &codon_table($genetic_code);
%genetic_code = %{$_[0]};
%start_codon = %{$_[1]};
%stop_codon = %{$_[2]};

# 读取基因组序列
my $genome_file = shift @ARGV;
my @GFF3_file = @ARGV;
open IN, $genome_file or die "Can not open file $genome_file, $!";
my (%genome_seq, $genome_seq_id);
while (<IN>) {
    chomp;
    if (m/^>(\S+)/) { $genome_seq_id = $1; }
    else { $genome_seq{$genome_seq_id} .= $_; }
}
close IN;

# 解析GFF3文件内容
my %gff3_info;
foreach ( reverse @ARGV ) {
    my %get_info_from_GFF3 = &get_info_from_GFF3($_);
    foreach ( keys %get_info_from_GFF3 ) {
        $gff3_info{$_} = $get_info_from_GFF3{$_};
    }
}

# 获得GFF3文件中ID的顺序
my (@geneID, %geneID, %feature_name);
foreach my $input_file ( @GFF3_file ) {
    open IN, $input_file or die "Can not open file $input_file, $!";
    while (<IN>) {
        next if /^\s*$/;
        next if /^#/;
        if (/ID=([^;]+)/ && exists $gff3_info{$1}) {
            push @geneID, $1 unless exists $geneID{$1};
            $geneID{$1} = 1;
            @_ = split /\t/;
            $feature_name{$_[2]} = 1;
        }
    }
    close IN;
}

# 删除可能存在的输出文件
foreach ( keys %feature_name ) {
    unlink "$out_prefix.$_.fasta" if -e "$out_prefix.$_.fasta";
    if ( $_ eq "gene" ) {
        unlink "$out_prefix.CDS.fasta"; unlink "$out_prefix.cDNA.fasta"; unlink "$out_prefix.protein.fasta";
    }
}

# 对每个基因进行分析
my (%stats, %collected_output_file_name, %feature_num);
foreach my $gene_ID ( @geneID ) {
    my $header = $gff3_info{$gene_ID}{"header"};
    my @header_field = split /\t/, $header;
    $feature_num{$header_field[2]} ++;

    # 若添加 --only_gene_sequences 或 --only_coding_gene_sequences 参数时，略过feature name不是gene的信息。
    if ( $header_field[2] ne "gene" ) {
        next if $only_gene_sequences;
        next if $only_coding_gene_sequences;
    }

    # 当添加 --only_coding_gene_sequences 参数时，还需要检测 gene 模型是否包含CDS信息。
    my $if_contain_CDS = 0;

    # 若GFF3的Feature名称为gene，则需要分析第二层和第三层信息，输出cDNA(exon)、CDS和Protein序列。
    if ( $header_field[2] eq "gene" ) {
        my %mRNA_info;
        my @mRNA_ID = @{$gff3_info{$gene_ID}{"mRNA_ID"}};

        foreach my $mRNA_ID ( @mRNA_ID ) {
            my $mRNA_header = $gff3_info{$gene_ID}{"mRNA_header"}{$mRNA_ID};
            my @mRNA_header = split /\t/, $mRNA_header;
            my $mRNA_info = $gff3_info{$gene_ID}{"mRNA_info"}{$mRNA_ID};
            $if_contain_CDS = 1 if $mRNA_info =~ m/\tCDS\t/;

            # 计算转录本的CDS长度
            my @CDS = &get_feature($mRNA_info, "CDS");
            my $CDS_length = &cal_length(\@CDS);
            $mRNA_info{$mRNA_ID}{"CDS_length"} = $CDS_length;
            foreach (@CDS) { @_ = split /\t/; my $length = abs($_[1] - $_[0]) + 1; push @{$mRNA_info{$mRNA_ID}{"single_CDS_length"}}, $length; }

            # 计算转录本的exon长度
            my @exon = &get_feature($mRNA_info, "exon");
            @exon = @CDS unless @exon;
            my $exon_length = &cal_length(\@exon);
            $mRNA_info{$mRNA_ID}{"exon_length"} = $exon_length;
            foreach (@exon) { @_ = split /\t/; my $length = abs($_[1] - $_[0]) + 1; push @{$mRNA_info{$mRNA_ID}{"single_exon_length"}}, $length; }

            # 计算转录本的intron长度
            my @intron = &exons2intron(@exon);
            my $intron_length = 0;
            $intron_length = &cal_length(\@intron);
            $mRNA_info{$mRNA_ID}{"intron_length"} = $intron_length;
            foreach (@intron) { @_ = split /\t/; my $length = abs($_[1] - $_[0]) + 1; push @{$mRNA_info{$mRNA_ID}{"single_intron_length"}}, $length; }
        }

        # 若是编码基因，则统计具有可变剪接的基因数量，基因的可变剪接数量，基因的位置信息（用于计算基因间区）。
        if ( $if_contain_CDS == 1 ) {
            $feature_num{"coding_gene"} ++;
            $stats{"coding_gene_num"} ++;
            push @{$stats{"gene_length"}}, abs($header_field[4] - $header_field[3]) + 1;
            my $isoform_num = @mRNA_ID;
            push @{$stats{"isoform_num"}}, $isoform_num;
            $stats{"AS_gene_num"} ++ if $isoform_num >= 2;
            push @{$stats{"gene_pos"}}, "$header_field[0]\t$header_field[3]\t$header_field[4]";
        }

        # 若添加 --only_coding_gene_sequences 参数时，略过没有CDS信息的基因模型。
        if ( $only_coding_gene_sequences ) {
            next if $if_contain_CDS == 0;
        }

        # 对转录本进行排序，优先按CDS长度，然后按exon长度，最后按ID的ASCII编码。
        my @sort_mRNA_ID = sort { $mRNA_info{$b}{"CDS_length"} <=> $mRNA_info{$a}{"CDS_length"}
            or $mRNA_info{$b}{"exon_length"} <=> $mRNA_info{$a}{"exon_length"}
            or $a cmp $b } keys %mRNA_info;

        # 若添加 --only_first_isoform 参数，则仅对第一个可变剪接模型进行分析。
        if ( $only_first_isoform ) {
            @sort_mRNA_ID = ($mRNA_ID[0]);
        }
        # 若添加 --only_longest_isoform 参数，则仅对最优模型进行分析。
        elsif ( $only_longest_isoform ) {
            @sort_mRNA_ID = shift @sort_mRNA_ID;
        }
        elsif ( $sort_isoform ) {
        }
        else {
            @sort_mRNA_ID = @mRNA_ID;
        }

        foreach my $mRNA_ID ( @sort_mRNA_ID ) {
            my $mRNA_header = $gff3_info{$gene_ID}{"mRNA_header"}{$mRNA_ID};
            my @mRNA_header = split /\t/, $mRNA_header;
            my $mRNA_info = $gff3_info{$gene_ID}{"mRNA_info"}{$mRNA_ID};

            my $strand = $mRNA_header[6];
            my $genome_seq = $genome_seq{$mRNA_header[0]};

            # 计算 CDS, exon 和 Protein 序列
            my ($seq_CDS, $seq_cDNA, $seq_protein) = &get_seqs($mRNA_info, $genome_seq, $strand);

            # 输出 CDS, exon 和 Protein 序列
            my $header_output = &get_fasta_header($mRNA_ID, $mRNA_header[8]);
            if ( $seq_CDS ) {
                open OUT, ">>", "$out_prefix.CDS.fasta" or die "Can not create file $out_prefix.CDS.fasta, $!";
                print OUT "$header_output$seq_CDS\n";
                close OUT;
                $collected_output_file_name{"$out_prefix.CDS.fasta"} = 1;
            }
            if ( $seq_cDNA ) {
                open OUT, ">>", "$out_prefix.cDNA.fasta" or die "Can not create file $out_prefix.cDNA.fasta, $!";
                print OUT "$header_output$seq_cDNA\n";
                close OUT;
                $collected_output_file_name{"$out_prefix.cDNA.fasta"} = 1;
            }
            if ( $seq_protein ) {
                open OUT, ">>", "$out_prefix.protein.fasta" or die "Can not create file $out_prefix.protein.fasta, $!";
                print OUT "$header_output$seq_protein\n";
                close OUT;
                $collected_output_file_name{"$out_prefix.protein.fasta"} = 1;
            }

            if ( $if_contain_CDS == 1 ) {
                # 统计：single_CDS, single_exon, single_intron
                foreach ( @{$mRNA_info{$mRNA_ID}{"single_CDS_length"}} ) { push @{$stats{"single_CDS_length"}}, $_; }
                foreach ( @{$mRNA_info{$mRNA_ID}{"single_exon_length"}} ) { push @{$stats{"single_exon_length"}}, $_; }
                foreach ( @{$mRNA_info{$mRNA_ID}{"single_intron_length"}} ) { push @{$stats{"single_intron_length"}}, $_; }

                # 统计：CDS长度、exon长度、intron长度
                push @{$stats{"CDS_length"}}, $mRNA_info{$mRNA_ID}{"CDS_length"} if exists $mRNA_info{$mRNA_ID}{"CDS_length"};
                push @{$stats{"exon_length"}}, $mRNA_info{$mRNA_ID}{"exon_length"} if exists $mRNA_info{$mRNA_ID}{"exon_length"};
                push @{$stats{"intron_length"}}, $mRNA_info{$mRNA_ID}{"intron_length"} if exists $mRNA_info{$mRNA_ID}{"intron_length"};

                # 统计：CDS个数, exon个数, intron个数
                my @CDS_num = @{$mRNA_info{$mRNA_ID}{"single_CDS_length"}}; my $CDS_num = 0; $CDS_num = @CDS_num if @CDS_num;
                my @exon_num = @{$mRNA_info{$mRNA_ID}{"single_exon_length"}}; my $exon_num = 0; $exon_num = @exon_num if @exon_num;
                my @intron_num = @{$mRNA_info{$mRNA_ID}{"single_intron_length"}}; my $intron_num = 0; $intron_num = @intron_num if @intron_num;
                push @{$stats{"CDS_num"}}, $CDS_num;
                push @{$stats{"exon_num"}}, $exon_num;
                push @{$stats{"intron_num"}}, $intron_num;
            }
        }
    }

    # 输出GFF3第一层的序列信息。
    my $start_site = $header_field[3] - 1;
    my $seq_length = abs($header_field[4] - $header_field[3]) + 1;
    push @{$stats{$header_field[2]}}, $seq_length;
    my $sequence_output = substr($genome_seq{$header_field[0]}, $start_site, $seq_length);
    open OUT, ">>", "$out_prefix.$header_field[2].fasta" or die "Can not create file $out_prefix.$header_field[2].fasta, $!";
    $collected_output_file_name{"$out_prefix.$header_field[2].fasta"} = 1 if $header_field[2];
    $sequence_output = &rc($sequence_output) if $header_field[6] eq "-";
    my $header_output = &get_fasta_header($gene_ID, $header_field[8]);
    print OUT "$header_output$sequence_output\n";
    close OUT;
}

if ( -e "$out_prefix.CDS.fasta" ) {
    # 计算基因间区长度
    @_ = &cal_intergenic_length(@{$stats{"gene_pos"}});
    my @intergenic_length1 = @{$_[0]};
    my @intergenic_length2 = @{$_[1]};
    foreach ( @intergenic_length1 ) { push @{$stats{"intergenic_length >= 0"}}, $_; }
    foreach ( @intergenic_length2 ) { push @{$stats{"intergenic_length < 0"}}, $_; }
    my $intergenic_length1_num = @intergenic_length1;
    my $intergenic_length2_num = @intergenic_length2;

    # 输出编码基因的统计结果。
    open OUT, ">", "$out_prefix.codingGeneModels.stats" or die "Can not create file $out_prefix.codingGeneModels.stats, $!";
    $collected_output_file_name{"$out_prefix.codingGeneModels.stats"} = 1;
    my ($coding_gene_num, $AS_gene_num) = (0, 0);
    $coding_gene_num = $stats{"coding_gene_num"} if exists $stats{"coding_gene_num"};
    $AS_gene_num = $stats{"AS_gene_num"} if exists $stats{"AS_gene_num"};
    printf OUT "%30s \t$coding_gene_num\n", "coding_gene number:";
    printf OUT "%30s \t$AS_gene_num\n", "AS_gene number:";
    printf OUT "%30s \t$intergenic_length1_num\n", "intergenic_length >= 0 number:"; 
    printf OUT "%30s \t$intergenic_length2_num\n\n", "intergenic_length < 0 number:"; 
    print OUT " " x 23 . "\tMedian  \tMean\n";

    my @item = ("isoform_num", "gene_length", "exon_length", "CDS_length", "intron_length", "CDS_num", "exon_num", "intron_num", "single_CDS_length", "single_exon_length", "single_intron_length", "intergenic_length >= 0", "intergenic_length < 0");
    foreach ( @item ) {
        my @intput_data;
        foreach ( @{$stats{$_}} ) { push @intput_data, $_; }
        @intput_data = sort {$a <=> $b} @intput_data;
        my $median = 0;
        $median = $intput_data[@intput_data/2] if @intput_data;
        my $total = 0;
        foreach ( @intput_data ) { $total += $_; }
        my $mean = 0;
        $mean = int($total / @intput_data * 100 + 0.5) / 100 if @intput_data > 0;
        printf OUT "%22s \t%-7s \t$mean\n", $_, $median;
    }
    close OUT;
}

# 输出各个 Feature 的数量信息。
if ( %feature_num ) {
    foreach ( sort keys %feature_num ) {
        print STDERR "The GFF3 files contain $feature_num{$_} $_\n";
    }
    print STDERR "\n";
}

# 输出结果文件名
if ( %collected_output_file_name ) {
    print STDERR "The output files are :\n";
    foreach ( sort keys %collected_output_file_name ) {
        print STDERR "\t$_\n";
    }
}
else {
    print STDERR "Warning: none files were output.\n";
}

sub cal_intergenic_length {
    my @input = @_;
    my (%gene_in_chr1, %gene_in_chr2);
    foreach ( @input ) {
        @_ = split /\t/;
        $gene_in_chr1{$_[0]}{"$_[1]\t$_[2]"} = $_[1];
        $gene_in_chr2{$_[0]}{"$_[1]\t$_[2]"} = $_[2];
    }

    my (@intergenic_length1, @intergenic_length2);
    foreach my $chr ( keys %gene_in_chr1 ) {
        my @region = sort {$gene_in_chr1{$chr}{$a} <=> $gene_in_chr1{$chr}{$b} or $gene_in_chr2{$chr}{$a} <=> $gene_in_chr2{$chr}{$b}} keys %{$gene_in_chr1{$chr}};
        my $first_region = shift @region;
        my ($last_start, $last_end) = split /\t/, $first_region;
        foreach ( @region ) {
            @_ = split /\t/, $_;
            if ( $_[0] > $last_end ) {
                push @intergenic_length1, $_[0] - $last_end - 1;
                $last_end = $_[1];
            }
            else {
                push @intergenic_length2, $last_end - $_[0] + 1;
                $last_end = $_[1] if $_[1] > $last_end;
            }
        }
    }

    return (\@intergenic_length1, \@intergenic_length2);
}

sub get_seqs {
    my ($info, $genome_seq, $strand) = @_;

    # 提取 CDS / exon 信息
    my (@CDS, @exon);
    foreach ( split /\n/, $info ) {
        @_ = split /\t/;
        ($_[4], $_[3]) = ($_[3], $_[4]) if $_[3] > $_[4];

        if ( $_[2] eq "CDS") {
            push @CDS, "$_[3]\t$_[4]\t$_[7]";
        }
        elsif ( $_[2] eq "exon" ) {
            push @exon, "$_[3]\t$_[4]";
        }
    }

    # 给 CDS / exon 排序
    @CDS = sort {$a <=> $b} @CDS; @exon = sort {$a <=> $b} @exon;

    # 得到 CDS 和 exon 序列
    my ($seq_CDS, $seq_exon);
    foreach ( @CDS ) { @_ = split /\t/; $seq_CDS .= substr($genome_seq, $_[0] - 1, $_[1] - $_[0] + 1); }
    foreach ( @exon ) { @_ = split /\t/; $seq_exon .= substr($genome_seq, $_[0] - 1, $_[1] - $_[0] + 1); }
    $seq_CDS = &rc($seq_CDS) if $strand eq "-";
    $seq_exon = &rc($seq_exon) if $strand eq "-";

    # 获得起始CDS的frame值
    @CDS = reverse @CDS if $strand eq "-";
    my @frame = split /\t/, $CDS[0];
    my $frame = 0;
    $frame = $frame[2] if $frame[2];

    # 得到 Protein 序列
    my $seq_protein;
    my $seq = $seq_CDS;
    # 处理第一个CDS第一个密码子。若第一个CDS的移码框不是0，表示基因不完整，使用目标遗传密码翻译；若移码框是0，则看第一个密码子是否是起始密码子，若是起始密码子，则翻译为M，否则按目标遗传密码翻译。
    if ( $frame > 0 ) {
        $seq =~ s/\w{$frame}//;
        if ( length ($seq) >= 3 ) {
            $seq =~ s/(\w{3})//;
            if (exists $genetic_code{$1}) {
                $seq_protein .= $genetic_code{$1};
            }
            else {
                $seq_protein .= "X";
            }
        }
    }
    else {
        if ( length ($seq) >= 3 ) {
            $seq =~ s/(\w{3})//;
            if (exists $start_codon{$1}) {
                $seq_protein .= "M";
            }
            elsif ( exists $genetic_code{$1} ) {
                $seq_protein .= $genetic_code{$1};
            }
            else {
                $seq_protein .= "X";
            }
        }
    }
    # 对除第一个和最后一个的其余密码子，按目标遗传密码翻译。
    while ((length $seq) > 3) {
        $seq =~ s/(\w{3})//;
        if (exists $genetic_code{$1}) {
            $seq_protein .= $genetic_code{$1};
        }
        else {
            $seq_protein .= "X";
        }
    }
    # 对最后一个密码子，看是否是终止密码子。
    if ( (length $seq) == 3 ) {
        if (exists $stop_codon{$seq}) {
            $seq_protein .= "*";
        }
        elsif ( exists $genetic_code{$seq} ) {
            $seq_protein .= $genetic_code{$seq};
        }
        else {
            $seq_protein .= "X";
        }
    }

    return ($seq_CDS, $seq_exon, $seq_protein);
}

sub get_feature {
    my ($info, $feature_name) = @_;
    my @out;
    foreach ( split /\n/, $info ) {
        @_ = split /\t/;
        if ( $_[2] eq $feature_name ) {
            push @out, "$_[3]\t$_[4]";
        }
    }
    return @out;
}

sub cal_length {
    my @region = @{$_[0]};
    my $out = 0;
    foreach (@region) {
        @_ = split /\t/;
        $out += abs($_[1] - $_[0]) + 1;
    }
    return $out;
}

sub exons2intron {
    my %intron;
    my @exon = sort {$a <=> $b} @_;
    my $first_exon = shift @exon;
    my ($last_start, $last_end) = split /\t/, $first_exon;
    foreach ( @exon ) {
        my ($start, $end) = split /\t/, $_;
        if ($start > $last_end) {
            my $intron_start = $last_end + 1;
            my $intron_stop = $start - 1;
            $intron{"$intron_start\t$intron_stop"} = 1 if $intron_stop >= $intron_start;
        }
        ($last_start, $last_end) = ($start, $end);
    }
    my @intron = sort {$a <=> $b} keys %intron;
    return @intron;
}

sub get_fasta_header {
    my ($id, $input) = ($_[0], $_[1]); chomp($input);
    my $out = ">$id";
    my @out;
    foreach ( split /;/, $input ) {
        push @out, "[$_]" unless m/ID=/;
    }
    my $out2 = join " ", @out;
    $out .= " $out2\n" if $out2;
    return $out;
}

sub rc {
    my $seq = shift @_;
    $seq = reverse $seq;
    $seq =~ tr/ATCGatcgn/TAGCTAGCN/;
    return $seq;
}

# 子程序，返回基因的GFF3哈希信息：
# gene_ID => "header" => gene_header
# gene_ID => "mRNA_ID" => 数组
# gene_ID => "mRNA_header" => mRNA_ID => mRNA_header
# gene_ID => "mRNA_info" => mRNA_ID => mRNA_Info
sub get_info_from_GFF3 {
    my @input = @_;
    my %gene_info;
    # 第一轮，找第一层信息。即没有Parent标签的信息。
    open IN, $input[0] or die "Can not open file $input[0], $!";
    while (<IN>) {
        if ( ! m/Parent=/ && m/ID=([^;\s]+)/ ) {
            $gene_info{$1}{"header"} = $_;
        }
    }
    close IN;
    # 第二轮，找第二层信息。Parent值是geneID的信息，包含但不限于 mRNA 信息
    my %mRNA_ID2gene_ID;
    open IN, $input[0] or die "Can not open file $input[0], $!";
    while (<IN>) {
        if ( m/Parent=([^;\s]+)/ ) {
            my $parent = $1;
            if ( exists $gene_info{$parent} ) {
                if ( m/ID=([^;\s]+)/ ) {
                    push @{$gene_info{$parent}{"mRNA_ID"}}, $1;
                    $gene_info{$parent}{"mRNA_header"}{$1} = $_;
                    $mRNA_ID2gene_ID{$1} = $parent;
                }
            }
        }
    }
    close IN;
    # 第三轮，找第三层信息。找Parent值不是geneID的信息
    open IN, $input[0] or die "Can not open file $input[0], $!";
    while (<IN>) {
        if ( m/Parent=([^;\s]+)/ && exists $mRNA_ID2gene_ID{$1} ) {
            my $parent = $1;
            $gene_info{$mRNA_ID2gene_ID{$1}}{"mRNA_info"}{$parent} .= $_;
        }
    }
    close IN;

    return %gene_info;
}

sub codon_table {
    my $genetic_code = $_[0];
    my %code = (
        "TTT" => "F",
        "TTC" => "F",
        "TTA" => "L",
        "TTG" => "L",
        "TCT" => "S",
        "TCC" => "S",
        "TCA" => "S",
        "TCG" => "S",
        "TAT" => "Y",
        "TAC" => "Y",
        "TAA" => "X",
        "TAG" => "X",
        "TGT" => "C",
        "TGC" => "C",
        "TGA" => "X",
        "TGG" => "W",
        "CTT" => "L",
        "CTC" => "L",
        "CTA" => "L",
        "CTG" => "L",
        "CCT" => "P",
        "CCC" => "P",
        "CCA" => "P",
        "CCG" => "P",
        "CAT" => "H",
        "CAC" => "H",
        "CAA" => "Q",
        "CAG" => "Q",
        "CGT" => "R",
        "CGC" => "R",
        "CGA" => "R",
        "CGG" => "R",
        "ATT" => "I",
        "ATC" => "I",
        "ATA" => "I",
        "ATG" => "M",
        "ACT" => "T",
        "ACC" => "T",
        "ACA" => "T",
        "ACG" => "T",
        "AAT" => "N",
        "AAC" => "N",
        "AAA" => "K",
        "AAG" => "K",
        "AGT" => "S",
        "AGC" => "S",
        "AGA" => "R",
        "AGG" => "R",
        "GTT" => "V",
        "GTC" => "V",
        "GTA" => "V",
        "GTG" => "V",
        "GCT" => "A",
        "GCC" => "A",
        "GCA" => "A",
        "GCG" => "A",
        "GAT" => "D",
        "GAC" => "D",
        "GAA" => "E",
        "GAG" => "E",
        "GGT" => "G",
        "GGC" => "G",
        "GGA" => "G",
        "GGG" => "G",
    );
    my %start_codon;
    $start_codon{"ATG"} = 1;
    if ( $genetic_code == 1 ) {
        # The Standard Code
        $start_codon{"TTG"} = 1;
        $start_codon{"CTG"} = 1;
    }
    elsif ( $genetic_code == 2 ) {
        # The Vertebrate Mitochondrial Code
        $code{"AGA"} = "X";
        $code{"AGG"} = "X";
        $code{"ATA"} = "M";
        $code{"TGA"} = "W";
        $start_codon{"ATA"} = 1;
        $start_codon{"ATT"} = 1;
        $start_codon{"ATC"} = 1;
        $start_codon{"GTG"} = 1;
    }
    elsif ( $genetic_code == 3 ) {
        # The Yeast Mitochondrial Code
        $code{"ATA"} = "M";
        $code{"CTT"} = "T";
        $code{"CTC"} = "T";
        $code{"CTA"} = "T";
        $code{"CTG"} = "T";
        $code{"TGA"} = "W";
        $start_codon{"ATA"} = 1;
        $start_codon{"GTG"} = 1;
    }
    elsif ( $genetic_code == 4 ) {
        # The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
        $code{"TGA"} = "W";
        $start_codon{"ATA"} = 1;
        $start_codon{"ATT"} = 1;
        $start_codon{"ATC"} = 1;
        $start_codon{"GTG"} = 1;
        $start_codon{"CTG"} = 1;
        $start_codon{"TTA"} = 1;
        $start_codon{"TTG"} = 1;
    }
    elsif ( $genetic_code == 5 ) {
        # The Invertebrate Mitochondrial Code
        $code{"AGA"} = "S";
        $code{"AGG"} = "S";
        $code{"ATA"} = "M";
        $code{"TGA"} = "W";
        $start_codon{"ATA"} = 1;
        $start_codon{"ATT"} = 1;
        $start_codon{"ATC"} = 1;
        $start_codon{"GTG"} = 1;
        $start_codon{"TTG"} = 1;
    }
    elsif ( $genetic_code == 6 ) {
        # The Ciliate, Dasycladacean and Hexamita Nuclear Code
        $code{"TAA"} = "Q";
        $code{"TAG"} = "Q";
    }
    elsif ( $genetic_code == 9 ) {
        # The Echinoderm and Flatworm Mitochondrial Code
        $code{"AAA"} = "N";
        $code{"AGA"} = "S";
        $code{"AGG"} = "S";
        $code{"TGA"} = "W";
        $start_codon{"GTG"} = 1;
    }
    elsif ( $genetic_code == 10 ) {
        # The Euplotid Nuclear Code
        $code{"TGA"} = "C";
    }
    elsif ( $genetic_code == 11 ) {
        # The Bacterial, Archaeal and Plant Plastid Code
        $start_codon{"ATA"} = 1;
        $start_codon{"ATT"} = 1;
        $start_codon{"ATC"} = 1;
        $start_codon{"GTG"} = 1;
        $start_codon{"CTG"} = 1;
        $start_codon{"TTG"} = 1;
    }
    elsif ( $genetic_code == 12 ) {
        # The Alternative Yeast Nuclear Code
        $code{"CTG"} = "S";
        $start_codon{"CTG"} = 1;
    }
    elsif ( $genetic_code == 13 ) {
        # The Ascidian Mitochondrial Code
        $code{"AGA"} = "G";
        $code{"AGG"} = "G";
        $code{"ATA"} = "M";
        $code{"TGA"} = "W";
        $start_codon{"ATA"} = 1;
        $start_codon{"GTG"} = 1;
        $start_codon{"TTG"} = 1;
    }
    elsif ( $genetic_code == 14 ) {
        # The Alternative Flatworm Mitochondrial Code
        $code{"AAA"} = "N";
        $code{"AGA"} = "S";
        $code{"AGG"} = "S";
        $code{"TAA"} = "Y";
        $code{"TGA"} = "W";
    }
    elsif ( $genetic_code == 16 ) {
        # Chlorophycean Mitochondrial Code
        $code{"TAG"} = "L";
    }
    elsif ( $genetic_code == 21 ) {
        # Trematode Mitochondrial Code
        $code{"TGA"} = "W";
        $code{"ATA"} = "M";
        $code{"AGA"} = "S";
        $code{"AGG"} = "S";
        $code{"AAA"} = "N";
        $start_codon{"GTG"} = 1;
    }
    elsif ( $genetic_code == 22 ) {
        # Scenedesmus obliquus Mitochondrial Code
        $code{"TCA"} = "X";
        $code{"TAG"} = "L";
    }
    elsif ( $genetic_code == 23 ) {
        # Thraustochytrium Mitochondrial Code
        $code{"TTA"} = "X";
        $start_codon{"ATT"} = 1;
        $start_codon{"GTG"} = 1;
    }
    elsif ( $genetic_code == 24 ) {
        # Rhabdopleuridae Mitochondrial Code
        $code{"AGA"} = "S";
        $code{"AGG"} = "K";
        $code{"TGA"} = "W";
        $start_codon{"GTG"} = 1;
        $start_codon{"CTG"} = 1;
        $start_codon{"TTG"} = 1;
    }
    elsif ( $genetic_code == 25 ) {
        # Candidate Division SR1 and Gracilibacteria Code
        $code{"TGA"} = "G";
        $start_codon{"GTG"} = 1;
        $start_codon{"TTG"} = 1;
    }
    elsif ( $genetic_code == 26 ) {
        # Pachysolen tannophilus Nuclear Code
        # warning: The descritpions of initiation codons by 2 methods are confict according to the NCBI web site.
        $code{"CTG"} = "A";
        $start_codon{"GTG"} = 1;
        $start_codon{"TTG"} = 1;
    }
    elsif ( $genetic_code == 27 ) {
        # Karyorelict Nuclear Code
        $code{"TAG"} = "Q";
        $code{"TAA"} = "Q";
    }
    elsif ( $genetic_code == 29 ) {
        # Mesodinium Nuclear Code
        $code{"TAA"} = "Y";
        $code{"TAG"} = "Y";
    }
    elsif ( $genetic_code == 30 ) {
        # Peritrich Nuclear Code
        $code{"TAA"} = "E";
        $code{"TAG"} = "E";
    }
    elsif ( $genetic_code == 31 ) {
        # Blastocrithidia Nuclear Code
        $code{"TGA"} = "W";
    }
    elsif ( $genetic_code == 33 ) {
        # Cephalodiscidae Mitochondrial UAA-Tyr Code
        $code{"TAA"} = "Y";
        $code{"TGA"} = "Y";
        $code{"AGA"} = "S";
        $code{"AGG"} = "K";
    }

    my %stop_codon;
    foreach ( keys %code ) {
        $stop_codon{$_} = $code{$_} if $code{$_} eq "X";
    }

    my %genetic_codon;
    foreach ( keys %code ) {
        if ( $code{$_} eq "X" ) {
            $genetic_codon{$_} = "*";
        }
        else {
            $genetic_codon{$_} = $code{$_};
        }
    }

    return (\%genetic_codon, \%start_codon, \%stop_codon);
}
