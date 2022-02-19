#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 geneModels.gtf genome.fasta outPrefix > statistics.txt

    通过真核生物GTF文件对基因模型进行统计。

    程序统计了基因水平的信息：基因长度、起始位点、结束位点、正负链方向、转录本的个数、和前一个基因的距离。并得到基因正义链上的序列。统 计基因中拥有最长CDS转录本的信息：cDNA长度、exon个数、CDS长度、CDS个数、intron长度、intron个数、exon长度信息、cds长度信息、intron长度信 息。并得到该转录本的cDNA序列、CDS序列和氨基酸序列。

    得到几个文件：
    outPrefix.geneModels.statistics    对每个基因的统计信息，表格文件
    outPrefix.gene.fasta               基因正义链上的序列
    outPrefix.cDNA.fasta               cDNA序列
    outPrefix.cds.fasta                CDS序列 (若转录本的第一个CDS的Frame值不是0，这时程序会去除相应的1~2个碱基，以保证CDS序列的第一个碱基是密码子的第一个位点。这样能和蛋白质的氨基酸序列一致，有利于后续的一些运用。)
    outPrefix.protein.fasta                蛋白序列

    此外，程序对以上统计信息的中位数和算术平均数进行了分析，并将结果输出到标准输出。

USAGE
if (@ARGV==0){die $usage}

open IN, $ARGV[0] or die $!;
my (%scaffold, %exon_info, %cds_info, %exon_length, %cds_length, %boundary, %strand);
while (<IN>) {
    if (m/\texon\t.*gene_id\s+?\"(\S+?)\".*transcript_id\s+?\"(\S+?)\"/) {
        @_ = split /\t/;
        $scaffold{$_[0]}{$1}{$2} = 1;
        $exon_info{$1}{$2}{"$_[3]\t$_[4]\t$_[6]"} = 1;
        my $length = abs($_[4] - $_[3]) + 1;
        $exon_length{$2} += $length;
        $boundary{$1}{$_[3]} = 1;
        $boundary{$1}{$_[4]} = 1;
        $strand{$1} = $_[6];
    }
    elsif (m/\tCDS\t.*gene_id\s+?\"(\S+?)\".*transcript_id\s+?\"(\S+?)\"/) {
        @_ = split /\t/;
        $cds_info{$1}{$2}{"$_[3]\t$_[4]\t$_[7]"} = 1;
        my $length = abs($_[4] - $_[3]) + 1;
        $cds_length{$2} += $length;
    }
}
close IN;

open IN, $ARGV[1] or die $!;
my (%seq, $seq_id);
while (<IN>) {
    chomp;
    if (m/^>(\S+)/) { $seq_id = $1; }
    else { $_ = uc($_); $seq{$seq_id} .= $_; }
}
close IN;

open STATS, ">", "$ARGV[2].geneModels.statistics" or die $!;
open GENE, ">", "$ARGV[2].gene.fasta" or die $!;
open CDNA, ">", "$ARGV[2].cDNA.fasta" or die $!;
open CDS, ">", "$ARGV[2].CDS.fasta" or die $!;
open PEP, ">", "$ARGV[2].protein.fasta" or die $!;
print STATS "Gene_ID\tChromosome_ID\tStart_site\tEnd_site\tGene_length\tStrand\tIsoforms_number\tIntergenic_length\tcDNA_size\tExon_num\tCDS_size\tCDS_num\tIntron_size\tIntron_num\tExons_length\tCDS_length\tIntron_length\n";
foreach my $scaffold_id (sort keys %scaffold) {
    my @gene_id = sort keys %{$scaffold{$scaffold_id}};

    # 先对gene_id进行排序，依次按起始位点、结束位点、正负链方向进行排序。
    my (%gene_id_sort1, %gene_id_sort2, %gene_id_sort3);
    foreach my $gene_id (@gene_id) {
        my @boundary = sort {$a <=> $b} keys %{$boundary{$gene_id}};
        $gene_id_sort1{$gene_id} = $boundary[0];
        $gene_id_sort2{$gene_id} = $boundary[-1];
        $gene_id_sort3{$gene_id} = $strand{$gene_id};
    }
    @gene_id = sort { $gene_id_sort1{$a} <=> $gene_id_sort1{$b} or $gene_id_sort2{$a} <=> $gene_id_sort2{$b} or $gene_id_sort3{$a} cmp $gene_id_sort3{$b} } @gene_id;

    my $last_gene_end;
    foreach my $gene_id (@gene_id) {
        my $out_stats;
        # 进行基因水平的统计：基因长度、起始位点、结束位点、正负链方向、转录本的个数、和前一个基因的距离。
        my $gene_start = $gene_id_sort1{$gene_id};
        my $gene_end = $gene_id_sort2{$gene_id};
        my $gene_length = $gene_end - $gene_start + 1;
        my $gene_seq = substr($seq{$scaffold_id}, $gene_start - 1, $gene_length);
        my $strand = $strand{$gene_id};
        my @isoform = sort {$cds_length{$b} <=> $cds_length{$a} or $exon_length{$b} <=> $exon_length{$a} or $a cmp $b} keys %{$scaffold{$scaffold_id}{$gene_id}};
        my $isoform_number = @isoform;
        my $isoform = $isoform[0];
        my $intergenic_length = $gene_start - $last_gene_end - 1;
        if ($last_gene_end) {
            $out_stats = "$gene_id\t$scaffold_id\t$gene_start\t$gene_end\t$gene_length\t$strand\t$isoform_number\t$intergenic_length";
        }
        else {
            $out_stats = "$gene_id\t$scaffold_id\t$gene_start\t$gene_end\t$gene_length\t$strand\t$isoform_number\tNULL";
        }
        $last_gene_end = $gene_end;

        my @exon = &sort_exon_cds(keys %{$exon_info{$gene_id}{$isoform}});
        my @cds = &sort_exon_cds(keys %{$cds_info{$gene_id}{$isoform}});

        # 统计基因中拥有最长CDS转录本的信息：cDNA长度、exon个数、CDS长度、CDS个数、intron长度、intron个数、exon长度信息、cds长度信息、intron长度信息。
        my ($cDNA_size, $exon_num, $cds_size, $cds_num, $intron_size, $intron_num) = (0, 0, 0, 0, 0, 0);
        my ($exon_len, $cds_len, $intron_len, $intron_len, $strand, $last_exon_end);
        # 得到cDNA、cds和氨基酸序列。
        my ($cDNA_seq, $cds_seq);

        foreach (@exon) {
            @_ = split /\t/, $_;
            $strand = $_[2];
            my $exon_length = $_[1] - $_[0] + 1;
            $cDNA_size += $exon_length;
            $exon_num ++;
            $exon_len .= "$exon_length,";
            if ($last_exon_end) {
                my $intron_length = $_[0] - $last_exon_end - 1;
                if ($intron_length >= 0) {
                    $intron_size += $intron_length;
                    $intron_num ++;
                    $intron_len .= "$intron_length,";
                }
                else {
                    print STDERR "$gene_id: intron_length was calculated to a minus num\n";
                }
            }
            $last_exon_end = $_[1];

            $cDNA_seq .= substr($seq{$scaffold_id}, $_[0] - 1, $exon_length);
        }
        $exon_len =~ s/,$//;
        $intron_len =~ s/,$//;

        foreach (@cds) {
            @_ = split /\t/, $_;
            my $cds_length = $_[1] - $_[0] + 1;
            $cds_size += $cds_length;
            $cds_num ++;
            $cds_len .= "$cds_length,";

            $cds_seq .= substr($seq{$scaffold_id}, $_[0] - 1, $cds_length);
        }
        $cds_len =~ s/,$//;

        if ($intron_len) {
            $out_stats .= "\t$cDNA_size\t$exon_num\t$cds_size\t$cds_num\t$intron_size\t$intron_num\t$exon_len\t$cds_len\t$intron_len\n";
        }
        else {
            $out_stats .= "\t$cDNA_size\t$exon_num\t$cds_size\t$cds_num\t$intron_size\t$intron_num\t$exon_len\t$cds_len\tNULL\n";
        }

        print STATS $out_stats;

        if ($strand eq "+") {
            print GENE ">$gene_id\n$gene_seq\n";
            print CDNA ">$gene_id\n$cDNA_seq\n";
            my $frame = $1 if $cds[0] =~ m/(\d+)$/;
            my $pep = &cds2pep($cds_seq, $frame);
            print PEP ">$gene_id\n$pep\n";
            $cds_seq =~ s/\w{$frame}// if $frame > 0;
            print CDS ">$gene_id\n$cds_seq\n";
        }
        elsif ($strand eq "-") {
            $gene_seq = &rc($gene_seq);
            $cDNA_seq = &rc($cDNA_seq);
            $cds_seq = &rc($cds_seq);
            print GENE ">$gene_id\n$gene_seq\n";
            print CDNA ">$gene_id\n$cDNA_seq\n";
            my $frame = $1 if $cds[-1] =~ m/(\d+)$/;
            my $pep = &cds2pep($cds_seq, $frame, $gene_id);
            print PEP ">$gene_id\n$pep\n";
            $cds_seq =~ s/\w{$frame}// if $frame > 0;
            print CDS ">$gene_id\n$cds_seq\n";
        }
    }
}
close STATS;
close GENE;
close CDNA;
close CDS;
close PEP;

open IN, "$ARGV[2].geneModels.statistics" or die $!;
<IN>;
my (@gene_length, @isoform_number, @intergenic_length, @exon_size, @exon_num, @cds_size, @cds_num, @intron_size, @intron_num, @exon_len, @cds_len, @intron_len, $overlap_gene_number);
while (<IN>) {
    chomp;
    @_ = split /\t/;
    push @gene_length, $_[4];
    push @isoform_number, $_[6];
    push @intergenic_length, $_[7] if $_[7] ne "NULL" && $_[7] >= 0;
    $overlap_gene_number ++ if $_[7] ne "NULL" && $_[7] < 0;
    push @exon_size, $_[8];
    push @exon_num, $_[9];
    push @cds_size, $_[10];
    push @cds_num, $_[11];
    push @intron_size, $_[12];
    push @intron_num, $_[13];
    push @exon_len, split /,/, $_[14];
    push @cds_len, split /,/, $_[15];
    push @intron_len, split /,/, $_[16] if $_[16] ne "NULL";
}
close IN;

print "Item\tMedian\tMean\n";
my $median = &median(@gene_length);
my $mean = &mean(@gene_length);
print "Gene_length\t$median\t$mean\n";
$median = &median(@isoform_number);
$mean = &mean(@isoform_number);
print "Isoform_number\t$median\t$mean\n";
my $as_gene_num = 0;
foreach (@isoform_number) { $as_gene_num ++ if $_ > 1; }
print "AS_gene_number: $as_gene_num\n";
$median = &median(@intergenic_length);
$mean = &mean(@intergenic_length);
print "Intergenic_length > 0\t$median\t$mean\n";
print "Intergenic_length < 0 number: $overlap_gene_number\n";
$median = &median(@exon_size);
$mean = &mean(@exon_size);
print "cDNA_length\t$median\t$mean\n";
$median = &median(@exon_num);
$mean = &mean(@exon_num);
print "exon_number\t$median\t$mean\n";
$median = &median(@cds_size);
$mean = &mean(@cds_size);
print "CDS_length\t$median\t$mean\n";
$median = &median(@cds_num);
$mean = &mean(@cds_num);
print "CDS_num\t$median\t$mean\n";
$median = &median(@intron_size);
$mean = &mean(@intron_size);
print "Intron_length\t$median\t$mean\n";
$median = &median(@exon_len);
$mean = &mean(@exon_len);
print "Single_exon_length\t$median\t$mean\n";
$median = &median(@cds_len);
$mean = &mean(@cds_len);
print "Single_CDS_length\t$median\t$mean\n";
$median = &median(@intron_len);
$mean = &mean(@intron_len);
print "Single_intron_length\t$median\t$mean\n";

sub mean {
    my $number = @_;
    my $total = 0;
    foreach (@_) { $total += $_; }
    my $mean = $total / $number if  $number != 0;
    return $mean;
}

sub median {
    my @median = sort {$a <=> $b} @_;
    my $median = $median[@median/2];
    return $median;
}

sub sort_exon_cds {
    my @region = @_;
    my %sort;
    foreach (@region) {
        $sort{$_} = $1 if m/(\d+)/;
    }
    @region = sort {$sort{$a} <=> $sort{$b}} @region;
    return @region;
}

sub rc {
    my $seq = shift @_;
    $seq = reverse $seq;
    $seq =~ tr/ATCGatcgn/TAGCTAGCN/;
    return $seq;
}

sub cds2pep {
    my %cds2pep = (
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
        "TAA" => "*",
        "TAG" => "*",
        "TGT" => "C",
        "TGC" => "C",
        "TGA" => "*",
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
    my $seq = shift @_;
    $seq = uc($seq);
    my $fram = shift @_;
    my $gene = shift @_;
    $seq =~ s/\w{$fram}// if $fram > 0;
    my $pep;
    while ((length $seq) >= 3) {
        $seq =~ s/(\w{3})//;
        if (exists $cds2pep{$1}) {
            $pep .= $cds2pep{$1};
        }
        else {
            $pep .= 'X';
        }
    }
    return $pep;
}
