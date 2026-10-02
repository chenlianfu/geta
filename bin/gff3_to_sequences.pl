#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $usage_chinese = <<USAGE;
Usage:
    perl $0 --out_prefix out genome.fasta input1.gff3 [input2.gff3 ...] > statitics.txt

    本程序用于根据GFF3的内容信息和基因组序列，转换出对应的序列信息，并对各种类型的Features进行统计。程序能对各种类型的Feature都进行序列转换，每种类型的Feature各得到一个Fasta文件，gene类型的Feature能得到CDS、cDNA和Protein三个Fasta文件。

    关于CDS和Protein序列的说明：程序根据GFF3中每一个CDS行第8列的phase信息来拼接CDS并翻译蛋白。沿转录方向，起始CDS的phase不为0时，去除开头对应的碱基，从第一个完整codon的第一个碱基开始输出CDS序列。对于后续的CDS，程序根据前面CDS累计长度推算其期望的phase：若与GFF3标注一致，则直接拼接(跨外显子的codon正常衔接)；若不一致(如存在移码frameshift)，则按GFF3标注的phase重新对齐读码框，被移码打断的不完整codon在Protein中记为X，在CDS中记为NNN，并在STDERR输出警告信息。

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
    设置遗传密码，程序据此自动确定密码子的翻译方式、起始密码子和终止密码子。该参数对应的值请参考NCBI Genetic Codes: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi 。该参数在将CDS转换为Protein序列时生效：终止密码子翻译为*；当基因5'端完整（起始CDS的phase为0）且第一个密码子是起始密码子时，该密码子翻译为M。支持的遗传密码编号为：1、2、3、4、5、6、9、10、11、12、13、14、15、16、21、22、23、24、25、26、27、28、29、30、31、32、33。例如：1为标准遗传密码（起始密码子为ATG，终止密码子为TAA、TAG、TGA），11为细菌/古菌/植物质体遗传密码，2为脊椎动物线粒体遗传密码。

    --help    default: None
    Display this English usage and exit.

    --chinese_help    default: None
    显示中文用法说明并退出。

USAGE
my $usage_english = &get_usage_english();
if (@ARGV==0){die $usage_english}

my ($help_flag, $chinese_help_flag, $out_prefix, $only_gene_sequences, $only_coding_gene_sequences, $only_longest_isoform, $only_first_isoform, $sort_isoform, $genetic_code);
GetOptions(
    "help" => \$help_flag,
    "chinese_help" => \$chinese_help_flag,
    "out_prefix:s" => \$out_prefix,
    "only_gene_sequences" => \$only_gene_sequences,
    "only_coding_gene_sequences" => \$only_coding_gene_sequences,
    "only_longest_isoform" => \$only_longest_isoform,
    "only_first_isoform" => \$only_first_isoform,
    "sort_isoforms" => \$sort_isoform,
    "genetic_code:i" => \$genetic_code,
) or die $usage_english;
if ( $help_flag ) { print STDERR $usage_english; exit 0; }
if ( $chinese_help_flag ) { print STDERR $usage_chinese; exit 0; }
$out_prefix ||= "out";

# 遗传密码：未设置或设置值不是正整数时，按标准遗传密码（1）处理
$genetic_code = 1 unless ( defined $genetic_code && $genetic_code > 0 );
my ( $genetic_code_ref, $start_codon_ref, $stop_codon_ref ) = codon_table($genetic_code);
my %genetic_code = %$genetic_code_ref;
my %start_codon  = %$start_codon_ref;
my %stop_codon   = %$stop_codon_ref;
die "Error: no start codon is available.\n" unless %start_codon;
die "Error: no stop codon is available.\n" unless %stop_codon;

# 读取基因组序列（整体读入后按记录切分处理，避免对每一行都做正则匹配+tr）
my $genome_file = shift @ARGV;
my @GFF3_file = @ARGV;
if ( @GFF3_file == 0 ) { die "Error: genome.fasta and at least one GFF3 file are required.\n\n$usage_english"; }
my %genome_seq = read_genome($genome_file);

# 解析每个GFF3文件：每个文件只读取一次磁盘，其余的"多趟扫描"都在内存中的数组上完成
my ( @file_gene_info, @file_order );
foreach my $file ( @GFF3_file ) {
    my ( $gene_info_ref, $order_ref ) = get_info_from_GFF3($file);
    push @file_gene_info, $gene_info_ref;
    push @file_order, $order_ref;
}

# 合并各文件的内容信息：若多个文件中出现相同ID，使用最靠前文件中的数据
# 做法：按文件逆序合并（靠后的文件先合并，靠前的文件后合并覆盖），与原逻辑等价
my %gff3_info;
for ( my $i = $#GFF3_file; $i >= 0; $i-- ) {
    my %gi = %{ $file_gene_info[$i] };
    foreach my $id ( keys %gi ) {
        $gff3_info{$id} = $gi{$id};
    }
}

# 确定最终的输出顺序：按文件先后顺序 + 文件内ID首次出现的顺序
my ( @geneID, %geneID_seen, %feature_name );
for ( my $i = 0; $i <= $#GFF3_file; $i++ ) {
    foreach my $id ( @{ $file_order[$i] } ) {
        next unless exists $gff3_info{$id};
        next if $geneID_seen{$id};
        push @geneID, $id;
        $geneID_seen{$id} = 1;
    }
}
foreach my $id ( @geneID ) {
    my $header = $gff3_info{$id}{"header"};
    my $type = (split /\t/, $header, 4)[2];
    $feature_name{$type} = 1 if defined $type;
}

# 生成空白的输出文件
my %output_fh;
foreach ( keys %feature_name ) {
    open $output_fh{$_}, ">", "$out_prefix.$_.fasta" or die "Can not create file $out_prefix.$_.fasta, $!";
    if ( $_ eq "gene" ) {
        open $output_fh{"CDS"}, ">", "$out_prefix.CDS.fasta"    or die "Can not create file $out_prefix.CDS.fasta, $!";
        open $output_fh{"cDNA"}, ">", "$out_prefix.cDNA.fasta"   or die "Can not create file $out_prefix.cDNA.fasta, $!";
        open $output_fh{"protein"}, ">", "$out_prefix.protein.fasta" or die "Can not create file $out_prefix.protein.fasta, $!";
    }
}

# 对每个基因进行分析
my ( %stats, %collected_output_file_name, %feature_num );
my $phase_inconsistent_num = 0;
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
        my %mRNA_stats;
        my %parsed_features;   # 缓存每个mRNA解析出的CDS/exon结构，避免重复parse
        next unless exists $gff3_info{$gene_ID}{"mRNA_ID"};
        my @mRNA_ID = @{ $gff3_info{$gene_ID}{"mRNA_ID"} };

        foreach my $mRNA_ID ( @mRNA_ID ) {
            my $mRNA_info = $gff3_info{$gene_ID}{"mRNA_info"}{$mRNA_ID};
            $mRNA_info = '' unless defined $mRNA_info;
            $if_contain_CDS = 1 if $mRNA_info =~ m/\tCDS\t/;

            # 一次性解析出CDS/exon（若无显式exon行，用CDS兜底），后面统计和输出序列都复用这份数据
            my ( $CDS_ref, $exon_ref ) = parse_mRNA_features($mRNA_info);
            $parsed_features{$mRNA_ID} = [ $CDS_ref, $exon_ref ];

            my $CDS_length = cal_length($CDS_ref);
            $mRNA_stats{$mRNA_ID}{"CDS_length"} = $CDS_length;
            $mRNA_stats{$mRNA_ID}{"single_CDS_length"} = [ map { $_->[1] - $_->[0] + 1 } @$CDS_ref ];

            my $exon_length = cal_length($exon_ref);
            $mRNA_stats{$mRNA_ID}{"exon_length"} = $exon_length;
            $mRNA_stats{$mRNA_ID}{"single_exon_length"} = [ map { $_->[1] - $_->[0] + 1 } @$exon_ref ];

            my $intron_ref = exons2intron($exon_ref);
            my $intron_length = cal_length($intron_ref);
            $mRNA_stats{$mRNA_ID}{"intron_length"} = $intron_length;
            $mRNA_stats{$mRNA_ID}{"single_intron_length"} = [ map { $_->[1] - $_->[0] + 1 } @$intron_ref ];
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
        my @sort_mRNA_ID = sort { $mRNA_stats{$b}{"CDS_length"} <=> $mRNA_stats{$a}{"CDS_length"}
            or $mRNA_stats{$b}{"exon_length"} <=> $mRNA_stats{$a}{"exon_length"}
            or $a cmp $b } keys %mRNA_stats;

        # 若添加 --only_first_isoform 参数，则仅对第一个可变剪接模型进行分析。
        if ( $only_first_isoform ) {
            @sort_mRNA_ID = ($mRNA_ID[0]);
        }
        # 若添加 --only_longest_isoform 参数，则仅对最优模型进行分析。
        elsif ( $only_longest_isoform ) {
            @sort_mRNA_ID = (shift @sort_mRNA_ID);
        }
        elsif ( $sort_isoform ) {
        }
        else {
            @sort_mRNA_ID = @mRNA_ID;
        }

        foreach my $mRNA_ID ( @sort_mRNA_ID ) {
            my $mRNA_header = $gff3_info{$gene_ID}{"mRNA_header"}{$mRNA_ID};
            my @mRNA_header = split /\t/, $mRNA_header;

            my $strand = $mRNA_header[6];
            my $genome_seq = $genome_seq{$mRNA_header[0]};

            # 复用前面已经解析好的CDS/exon结构，不再重新split
            my ( $CDS_ref, $exon_ref ) = @{ $parsed_features{$mRNA_ID} };
            my ( $seq_CDS, $seq_cDNA, $seq_protein, $phase_inconsistent ) = get_seqs($CDS_ref, $exon_ref, $genome_seq, $strand);

            if ( $phase_inconsistent ) {
                $phase_inconsistent_num ++;
                print STDERR "Warning: the phase of CDS in $mRNA_ID is inconsistent with the cumulative CDS length (possible frameshift); the reading frame is re-aligned according to the phase column of GFF3, and the broken codon is translated as X.\n";
            }

            # 输出 CDS, exon 和 Protein 序列
            my $header_output = get_fasta_header($mRNA_ID, $mRNA_header[8]);
            if ( $seq_CDS ) {
                print {$output_fh{"CDS"}} "$header_output$seq_CDS\n";
                $collected_output_file_name{"$out_prefix.CDS.fasta"} = 1;
            }
            if ( $seq_cDNA ) {
                print {$output_fh{"cDNA"}} "$header_output$seq_cDNA\n";
                $collected_output_file_name{"$out_prefix.cDNA.fasta"} = 1;
            }
            if ( $seq_protein ) {
                print {$output_fh{"protein"}} "$header_output$seq_protein\n";
                $collected_output_file_name{"$out_prefix.protein.fasta"} = 1;
            }

            if ( $if_contain_CDS == 1 ) {
                # 统计：single_CDS, single_exon, single_intron
                foreach ( @{$mRNA_stats{$mRNA_ID}{"single_CDS_length"}} )    { push @{$stats{"single_CDS_length"}}, $_; }
                foreach ( @{$mRNA_stats{$mRNA_ID}{"single_exon_length"}} )   { push @{$stats{"single_exon_length"}}, $_; }
                foreach ( @{$mRNA_stats{$mRNA_ID}{"single_intron_length"}} ) { push @{$stats{"single_intron_length"}}, $_; }

                # 统计：CDS长度、exon长度、intron长度
                push @{$stats{"CDS_length"}}, $mRNA_stats{$mRNA_ID}{"CDS_length"} if exists $mRNA_stats{$mRNA_ID}{"CDS_length"};
                push @{$stats{"exon_length"}}, $mRNA_stats{$mRNA_ID}{"exon_length"} if exists $mRNA_stats{$mRNA_ID}{"exon_length"};
                push @{$stats{"intron_length"}}, $mRNA_stats{$mRNA_ID}{"intron_length"} if exists $mRNA_stats{$mRNA_ID}{"intron_length"};

                # 统计：CDS个数, exon个数, intron个数
                my $CDS_num    = scalar @{$mRNA_stats{$mRNA_ID}{"single_CDS_length"}};
                my $exon_num   = scalar @{$mRNA_stats{$mRNA_ID}{"single_exon_length"}};
                my $intron_num = scalar @{$mRNA_stats{$mRNA_ID}{"single_intron_length"}};
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
    $collected_output_file_name{"$out_prefix.$header_field[2].fasta"} = 1 if $header_field[2];
    $sequence_output = rc($sequence_output) if $header_field[6] eq "-";
    my $header_output = get_fasta_header($gene_ID, $header_field[8]);
    print {$output_fh{$header_field[2]}} "$header_output$sequence_output\n";
}
foreach ( keys %output_fh ) {
    close $output_fh{$_};
}

if ( -e "$out_prefix.CDS.fasta" ) {
    # 计算基因间区长度
    my ( $intergenic_length1_ref, $intergenic_length2_ref ) = cal_intergenic_length(@{$stats{"gene_pos"}});
    my @intergenic_length1 = @$intergenic_length1_ref;
    my @intergenic_length2 = @$intergenic_length2_ref;
    foreach ( @intergenic_length1 ) { push @{$stats{"intergenic_length >= 0"}}, $_; }
    foreach ( @intergenic_length2 ) { push @{$stats{"intergenic_length < 0"}}, $_; }
    my $intergenic_length1_num = @intergenic_length1;
    my $intergenic_length2_num = @intergenic_length2;

    # 输出编码基因的统计结果。
    open my $STATS_OUT, ">", "$out_prefix.codingGeneModels.stats" or die "Can not create file $out_prefix.codingGeneModels.stats, $!";
    $collected_output_file_name{"$out_prefix.codingGeneModels.stats"} = 1;
    my ( $coding_gene_num, $AS_gene_num ) = (0, 0);
    $coding_gene_num = $stats{"coding_gene_num"} if exists $stats{"coding_gene_num"};
    $AS_gene_num = $stats{"AS_gene_num"} if exists $stats{"AS_gene_num"};
    printf $STATS_OUT "%30s \t$coding_gene_num\n", "coding_gene number:";
    printf $STATS_OUT "%30s \t$AS_gene_num\n", "AS_gene number:";
    printf $STATS_OUT "%30s \t$intergenic_length1_num\n", "intergenic_length >= 0 number:";
    printf $STATS_OUT "%30s \t$intergenic_length2_num\n\n", "intergenic_length < 0 number:";
    print $STATS_OUT " " x 23 . "\tMedian  \tMean\n";

    my @item = ("isoform_num", "gene_length", "exon_length", "CDS_length", "intron_length", "CDS_num", "exon_num", "intron_num", "single_CDS_length", "single_exon_length", "single_intron_length", "intergenic_length >= 0", "intergenic_length < 0");
    foreach ( @item ) {
        my @input_data = exists $stats{$_} ? @{$stats{$_}} : ();
        @input_data = sort { $a <=> $b } @input_data;
        my $median = 0;
        $median = $input_data[@input_data/2] if @input_data;
        my $total = 0;
        foreach ( @input_data ) { $total += $_; }
        my $mean = 0;
        $mean = int($total / @input_data * 100 + 0.5) / 100 if @input_data > 0;
        printf $STATS_OUT "%22s \t%-7s \t$mean\n", $_, $median;
    }
    close $STATS_OUT;
}

# 输出各个 Feature 的数量信息。
if ( %feature_num ) {
    foreach ( sort keys %feature_num ) {
        print STDERR "The GFF3 files contain $feature_num{$_} $_\n";
    }
    print STDERR "\n";
}
if ( $phase_inconsistent_num > 0 ) {
    print STDERR "Warning: $phase_inconsistent_num transcripts have CDS phases inconsistent with the cumulative CDS length (possible frameshift).\n\n";
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

# ------------------------------------------------------------------------------
# 子程序
# ------------------------------------------------------------------------------

# 读取基因组FASTA文件。
sub read_genome {
    my ($file) = @_;
    open my $fh, "<", $file or die "Can not open file $file, $!";
    local $/ = undef;
    my $content = <$fh>;
    close $fh;
    return () unless defined $content;
    $content =~ s/\r\n/\n/g;    # 兼容 Windows 换行符，避免混入序列中

    my %seq;
    my @records = split /^>/m, $content;
    shift @records if @records && $records[0] eq '';
    foreach my $rec ( @records ) {
        next unless length $rec;
        my ( $header, $rest ) = split /\n/, $rec, 2;
        next unless defined $header;
        my ($id) = split /\s+/, $header;
        next unless defined $id && length $id;
        my $s = defined $rest ? $rest : '';
        $s =~ s/\n//g;
        $s =~ tr/atcgn/ATCGN/;
        $seq{$id} = $s;
    }
    return %seq;
}

# 计算基因间区长度
sub cal_intergenic_length {
    my @input = @_;
    my ( %gene_in_chr1, %gene_in_chr2 );
    foreach my $line ( @input ) {
        my @f = split /\t/, $line;
        $gene_in_chr1{$f[0]}{"$f[1]\t$f[2]"} = $f[1];
        $gene_in_chr2{$f[0]}{"$f[1]\t$f[2]"} = $f[2];
    }

    my ( @intergenic_length1, @intergenic_length2 );
    foreach my $chr ( keys %gene_in_chr1 ) {
        my @region = sort { $gene_in_chr1{$chr}{$a} <=> $gene_in_chr1{$chr}{$b}
                          or $gene_in_chr2{$chr}{$a} <=> $gene_in_chr2{$chr}{$b} } keys %{$gene_in_chr1{$chr}};
        my $first_region = shift @region;
        my ( $last_start, $last_end ) = split /\t/, $first_region;
        foreach my $region ( @region ) {
            my ( $start, $end ) = split /\t/, $region;
            if ( $start > $last_end ) {
                push @intergenic_length1, $start - $last_end - 1;
                $last_end = $end;
            }
            else {
                push @intergenic_length2, $last_end - $start + 1;
                $last_end = $end if $end > $last_end;
            }
        }
    }

    return ( \@intergenic_length1, \@intergenic_length2 );
}

# 解析一个mRNA下属的CDS/exon信息（$info为多行GFF3文本，每行以\n结尾）。
# 返回 (\@CDS, \@exon)，元素分别为 [start, end, phase] 和 [start, end]，均已按start升序排序。
# 第8列phase若不是0/1/2（如"."），按0处理。
# 若该mRNA下没有显式的exon行，则用CDS的坐标作为exon坐标兜底。
sub parse_mRNA_features {
    my ($info) = @_;
    my ( @CDS, @exon );
    foreach my $line ( split /\n/, $info ) {
        next unless length $line;
        my @f = split /\t/, $line;
        next unless @f >= 8;
        my ( $start, $end ) = ( $f[3], $f[4] );
        ( $start, $end ) = ( $end, $start ) if $start > $end;
        if ( $f[2] eq "CDS" ) {
            my $phase = ( $f[7] =~ /^[012]$/ ) ? $f[7] : 0;
            push @CDS, [ $start, $end, $phase ];
        }
        elsif ( $f[2] eq "exon" ) {
            push @exon, [ $start, $end ];
        }
    }
    @CDS  = sort { $a->[0] <=> $b->[0] } @CDS;
    @exon = sort { $a->[0] <=> $b->[0] } @exon;
    @exon = map { [ $_->[0], $_->[1] ] } @CDS unless @exon;   # 无显式exon行时用CDS兜底
    return ( \@CDS, \@exon );
}

# 计算一组区间（arrayref of [start,end,...]）的总长度
sub cal_length {
    my ($regions) = @_;
    my $total = 0;
    foreach my $r ( @$regions ) {
        $total += abs($r->[1] - $r->[0]) + 1;
    }
    return $total;
}

# 由exon坐标（已排序的 arrayref of [start,end]）推算intron坐标，返回 arrayref of [start,end]
sub exons2intron {
    my ($exon_ref) = @_;
    return [] unless @$exon_ref;
    my @exon = @$exon_ref;
    my ( $last_start, $last_end ) = @{ shift @exon };
    my %intron;
    foreach my $e ( @exon ) {
        my ( $start, $end ) = @$e;
        if ( $start > $last_end ) {
            my $intron_start = $last_end + 1;
            my $intron_stop  = $start - 1;
            $intron{"$intron_start\t$intron_stop"} = 1 if $intron_stop >= $intron_start;
        }
        ( $last_start, $last_end ) = ( $start, $end );
    }
    my @intron = map { [ split /\t/, $_ ] } keys %intron;
    @intron = sort { $a->[0] <=> $b->[0] } @intron;
    return \@intron;
}

# 根据已解析好的CDS/exon结构生成CDS、cDNA(exon)、Protein序列。
#
# GFF3第8列phase的含义：该CDS片段开头有phase个碱基，属于上一个(跨外显子的)密码子的后半部分，
# 需要跳过这些碱基才能到达下一个完整密码子的第一个碱基。
#
# 处理方法（沿转录方向，正链按坐标升序，负链按坐标降序，每个片段单独取序列并在负链时单独反向互补）：
#   1. 第一个CDS：去除开头phase个碱基(它们属于上游不完整的密码子)。
#   2. 后续CDS：根据前面残留的不足3个碱基的数量，推算期望phase = (3 - 残留数) % 3。
#      - 若GFF3标注的phase与期望值一致：是正常的跨外显子衔接，不删除任何碱基，直接与残留碱基拼接。
#      - 若不一致(移码frameshift)：残留的不完整密码子记为X(Protein)和NNN(CDS)，
#        再按GFF3标注的phase去除本片段开头phase个碱基，重新对齐读码框。
# 返回 ($seq_CDS, $seq_exon, $seq_protein, $phase_inconsistent)
sub get_seqs {
    my ( $CDS_ref, $exon_ref, $genome_seq, $strand ) = @_;
    my @CDS  = @$CDS_ref;    # 已按起始坐标升序排序，元素为 [start, end, phase]
    my @exon = @$exon_ref;   # 已按起始坐标升序排序，元素为 [start, end]

    # cDNA(exon)序列
    my $seq_exon = join( '', map { substr( $genome_seq, $_->[0] - 1, $_->[1] - $_->[0] + 1 ) } @exon );
    $seq_exon = rc($seq_exon) if $strand eq "-";

    # 转为转录方向
    @CDS = reverse @CDS if $strand eq "-";

    my @codons;                  # 已确定的密码子；undef 表示被移码打断的不完整密码子
    my $carry = '';              # 当前残留的、不足3个碱基的序列
    my $phase_inconsistent = 0;
    my $first_phase = 0;
    my $is_first = 1;

    foreach my $c ( @CDS ) {
        my $s = substr( $genome_seq, $c->[0] - 1, $c->[1] - $c->[0] + 1 );
        $s = rc($s) if $strand eq "-";
        my $phase = $c->[2];
        my $pending;

        if ( $is_first ) {
            $first_phase = $phase;
            $is_first = 0;
            $pending = ( length($s) > $phase ) ? substr( $s, $phase ) : '';
        }
        else {
            my $expected = ( 3 - length($carry) ) % 3;
            if ( $phase == $expected ) {
                # 正常衔接：残留碱基 + 本片段完整序列
                $pending = $carry . $s;
            }
            else {
                # 移码：残留的不完整密码子无法补全，记为undef；按GFF3标注的phase重新对齐
                $phase_inconsistent = 1;
                push @codons, undef if length($carry) > 0;
                $pending = ( length($s) > $phase ) ? substr( $s, $phase ) : '';
            }
        }

        # 从pending中切出完整的密码子，不足3个的碱基作为新的残留
        my $n = int( length($pending) / 3 );
        for ( my $i = 0; $i < $n; $i++ ) {
            push @codons, substr( $pending, $i * 3, 3 );
        }
        $carry = substr( $pending, $n * 3 );
    }

    # 组装CDS（被移码打断的密码子用NNN占位；末尾不足3个的残余碱基保留在CDS末尾）
    my $seq_CDS = join( '', map { defined $_ ? $_ : "NNN" } @codons ) . $carry;

    # 翻译Protein
    my $seq_protein = '';
    for ( my $i = 0; $i < @codons; $i++ ) {
        my $codon = $codons[$i];
        if ( !defined $codon ) {
            $seq_protein .= "X";
        }
        elsif ( $i == 0 && $first_phase == 0 ) {
            # 第一个CDS的phase为0，说明5'端完整：若是起始密码子则翻译为M，否则按遗传密码翻译
            if ( exists $start_codon{$codon} )     { $seq_protein .= "M"; }
            elsif ( exists $genetic_code{$codon} ) { $seq_protein .= $genetic_code{$codon}; }
            else                                   { $seq_protein .= "X"; }
        }
        else {
            # 第一个CDS的phase>0表示基因5'端不完整，第一个密码子不当作起始密码子；
            # 其余密码子（包括末尾的终止密码子，翻译为"*"）按遗传密码翻译
            $seq_protein .= exists $genetic_code{$codon} ? $genetic_code{$codon} : "X";
        }
    }

    return ( $seq_CDS, $seq_exon, $seq_protein, $phase_inconsistent );
}

# 生成FASTA header行（保证末尾一定带换行符）
sub get_fasta_header {
    my ( $id, $input ) = @_;
    $input = '' unless defined $input;
    $input =~ s/\r?\n$//;
    my @extra;
    foreach ( split /;/, $input ) {
        push @extra, "[$_]" unless /ID=/;
    }
    my $out = ">$id";
    $out .= " " . join( " ", @extra ) if @extra;
    $out .= "\n";
    return $out;
}

# 取反向互补序列
sub rc {
    my ($seq) = @_;
    $seq = reverse $seq;
    $seq =~ tr/ATCGatcgn/TAGCTAGCN/;
    return $seq;
}

# 解析一个GFF3文件，返回 (\%gene_info, \@order)：
#   \%gene_info : gene_ID => "header" => gene_header
#                 gene_ID => "mRNA_ID" => 数组
#                 gene_ID => "mRNA_header" => mRNA_ID => mRNA_header
#                 gene_ID => "mRNA_info" => mRNA_ID => mRNA_Info（多行文本，每行以\n结尾）
#   \@order     : 该文件中第一层（无Parent）Feature ID首次出现的顺序
sub get_info_from_GFF3 {
    my ($file) = @_;
    open my $fh, "<", $file or die "Can not open file $file, $!";
    my @lines;
    while ( <$fh> ) {
        next if /^\s*$/;
        next if /^#/;
        s/\r?\n$//;
        push @lines, $_;
    }
    close $fh;

    # 预先解析每一行的 ID / Parent / Feature类型，供下面三趟处理复用
    my @parsed;
    foreach my $line ( @lines ) {
        my @f = split /\t/, $line, 9;
        my $type = $f[2];
        my $attr = defined $f[8] ? $f[8] : '';
        my $id     = ( $attr =~ /(?:^|;)ID=([^;\s]+)/ )     ? $1 : undef;
        my $parent = ( $attr =~ /(?:^|;)Parent=([^;\s]+)/ ) ? $1 : undef;
        push @parsed, [ $line, $id, $parent, $type ];
    }

    my ( %gene_info, @order, %seen_top );

    # 第一层：无 Parent 且有 ID
    foreach my $rec ( @parsed ) {
        my ( $line, $id, $parent, $type ) = @$rec;
        next unless defined $id;
        next if defined $parent;
        $gene_info{$id}{"header"} = $line;
        unless ( $seen_top{$id} ) {
            push @order, $id;
            $seen_top{$id} = 1;
        }
    }

    # 第二层：Parent 为第一层 ID（如 mRNA）
    my %mRNA_ID2gene_ID;
    foreach my $rec ( @parsed ) {
        my ( $line, $id, $parent, $type ) = @$rec;
        next unless defined $parent;
        next unless exists $gene_info{$parent};
        next unless defined $id;
        push @{ $gene_info{$parent}{"mRNA_ID"} }, $id
            unless exists $gene_info{$parent}{"mRNA_header"}{$id};
        $gene_info{$parent}{"mRNA_header"}{$id} = $line;
        $mRNA_ID2gene_ID{$id} = $parent;
    }

    # 第三层：Parent 为第二层 ID（如 CDS、exon）
    foreach my $rec ( @parsed ) {
        my ( $line, $id, $parent, $type ) = @$rec;
        next unless defined $parent;
        next unless exists $mRNA_ID2gene_ID{$parent};
        $gene_info{ $mRNA_ID2gene_ID{$parent} }{"mRNA_info"}{$parent} .= "$line\n";
    }

    return ( \%gene_info, \@order );
}

# 子程序，根据遗传密码编号返回（密码子翻译表、起始密码子表、终止密码子表）的哈希引用。
# 其中密码子翻译表中终止密码子翻译为"*"。遗传密码的设置参考NCBI Genetic Codes，代码来自 GFF3_merging_and_removing_redundancy。
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
        #$start_table{"TTG"} = 1;
        #$start_table{"CTG"} = 1;
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
    }

    my %stop_table;
    foreach ( keys %code ) {
        $stop_table{$_} = 1 if $code{$_} eq "X";
    }

    # 翻译表：终止密码子翻译为"*"
    my %translate_table;
    foreach ( keys %code ) {
        $translate_table{$_} = ( $code{$_} eq "X" ) ? "*" : $code{$_};
    }

    return (\%translate_table, \%start_table, \%stop_table);
}

# 英文用法说明
sub get_usage_english {

my $usage = <<USAGE;
Usage:
    perl $0 --out_prefix out genome.fasta input1.gff3 [input2.gff3 ...] > statitics.txt

    This program converts the content of GFF3 files and the genome sequence into the corresponding sequences, and makes statistics for the various types of Features. The program can convert the sequences of Features of every type, and one Fasta file is obtained for each type of Feature; for Features of the gene type, three Fasta files (CDS, cDNA and Protein) are obtained.

    Notes on the CDS and Protein sequences: the program joins the CDS and translates the protein according to the phase in column 8 of every CDS line in the GFF3. Along the transcription direction, if the phase of the first CDS is not 0, the corresponding bases at the beginning are removed, and the CDS sequence is output starting from the first base of the first complete codon. For each following CDS, the program derives the expected phase from the accumulated length of the preceding CDS: if it agrees with the GFF3 annotation, the sequences are joined directly (codons spanning exons are connected normally); if it does not agree (e.g. there is a frameshift), the reading frame is re-aligned according to the phase annotated in the GFF3, the incomplete codon broken by the frameshift is recorded as X in the Protein and as NNN in the CDS, and a warning is output to STDERR.

    The program accepts multiple GFF3 files and outputs sequences according to the Feature IDs in them, so column 9 of the input GFF3 files must contain ID information. If the same ID occurs several times in one file, only the data of its last occurrence is used; if the same ID occurs in several files, the data in the GFF3 file that is the earliest in the input order is used. When finally outputting the sequences, the program follows the order of the input GFF3 files and the order of appearance of the IDs in the GFF3 files.

    If the GFF3 files contain coding gene information, statistics are made for these coding genes on the CDS length, the cDNA length, the intron length, the gene length, the number of CDS, the number of exons, the number of introns, the single CDS length, the single exon length, the single intron length and the intergenic length, and the results are written to the file out.codingGeneModels.stats.

    --out_prefix <string>    default: out
    Set the prefix of the output sequence files. The program generates the corresponding Fasta files out.FeatureName.fasta according to the Feature Name of the GFF3 files. If the Feature Name is gene, the files out.CDS.fasta, out.cDNA.fasta and out.protein.fasta are additionally generated. If there is coding gene information, i.e. the genes contain CDS features, the file out.codingGeneModels.stats is additionally generated.

    --only_gene_sequences    default: None
    When this option is added, only the sequences of the gene type in the GFF3 files are output.

    --only_coding_gene_sequences    default: None
    When this option is added, only the sequence information of the coding genes in the GFF3 files is output.

    --only_first_isoform    default: None
    When this option is added, if a gene has several alternative splicing isoforms, only the first isoform that appears in the GFF3 file is selected for statistics and sequence output. Only one of this option and --only_longest_isoform can take effect; when both are set, this option is the effective one. For gene models predicted by the GETA software, the first isoform generally has the largest expression proportion.

    --only_longest_isoform    default: None
    When this option is added, if a gene has several alternative splicing isoforms, only the isoform with the longest CDS or the longest exons is selected for statistics and sequence output.

    --sort_isoforms    default: None
    When this option is added, the multiple alternative splicing transcripts of a gene model are sorted before their sequences are output. They are sorted first by CDS length from long to short, then by cDNA length from long to short, and finally by the ASCII code of the ID from small to large. By default the program outputs all alternative splicing sequences in the order in which they appear in the GFF3 file.

    --genetic_code <int>    default: 1
    Set the genetic code, from which the translation of codons, the start codons and the stop codons are determined automatically. For the corresponding values, please refer to NCBI Genetic Codes: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi . This option takes effect when the CDS is converted into the Protein sequence: stop codons are translated as *; when the 5' end of the gene is complete (the phase of the first CDS is 0) and the first codon is a start codon, that codon is translated as M. The supported genetic code numbers are: 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33. For example: 1 is the standard code (start codon ATG; stop codons TAA, TAG, TGA), 11 is the bacterial/archaeal/plant plastid code, and 2 is the vertebrate mitochondrial code.

    --help    default: None
    Display this English usage and exit.

    --chinese_help    default: None
    显示中文用法说明并退出。

USAGE

return $usage;
}
