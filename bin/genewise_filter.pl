#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    perl $0 [options] genewise.gff genome_seq.fasta Min_score Min_cds_length Min_exon_number hmmscan_Evalue hmmscan_coverage_ratio CPUs > genewise.filter.gff

For example:
    perl $0 genewise.gff genome_seq.fasta 25 300 2 1e-6 0.30 24 > genewise.filter.gff

    --pfam_db <string>    default: "/opt/biosoft/hmmer-3.1b2/Pfam-AB.hmm"

    示例中是推荐设置。即表示程序按如下流程运行：
    genewise得分低于15的直接过滤掉;
    cds长度低于90的基因直接过滤掉;
    genewise得分[15,25)的基因对应的proteins比对到Pfam数据库进行筛选；
    cds长度[90,300)的基因对应的proteins比对到Pfam数据库进行筛选；
    single exon基因对应的proteins比对到Pfam数据库进行筛选；

    程序调用了hmmscan对pfam数据库进行比对。因此，需要能直接运行hmmscan命令。同时，使用的数据库默认是"/opt/biosoft/hmmer-3.1b2/Pfam-AB.hmm"。
    同时程序调用了ParaFly命令进行并行化运算hmmscan，需要能直接运行ParaFly命令。示例中的并行化数目是24.

USAGE
if (@ARGV==0){die $usage}

my $Pfam_db;
GetOptions(
    "pfam_db:s" => \$Pfam_db,
);

$Pfam_db ||= "/opt/biosoft/hmmer-3.1b2/Pfam-AB.hmm";

open IN, $ARGV[0] or die $!;
my (%gff, @id, %id);
while (<IN>) {
    if (m/ID=(\S+?);/) {
        push @id, $1 unless exists $id{$1};
        $gff{$1} .= $_;
        $id{$1} = 1;
    }
}
close IN;

my (%gene_ok, %gene_pfam);
my ($num_score_lt15, $num_score_low, $num_cds_lt90, $num_cds_short, $num_cds_lack, $num_validate, $num_total);
foreach (@id) {
    my $gff = $gff{$_};
    $num_total ++;

    $gff =~ m/match\t\d+?\t\d+?\t(\S+?)\t/;
    my $score = $1;
    my @gff = split /\n/, $gff;
    my ($cds_length, $cds_number);
    foreach (@gff) {
        if (m/cds\t(\d+)\t(\d+)/) {
            $cds_length += abs($2 - $1 + 1);
            $cds_number ++;
        }
    }

    if ($score < 15) {
        $num_score_lt15 ++;
        next;
    }
    elsif ($cds_length < 90) {
        $num_cds_lt90 ++;
        next;
    }
    elsif ($score < $ARGV[2]) {
        $gene_pfam{$_} = 1;
        $num_score_low ++;
    }
    elsif ($cds_length < $ARGV[3]) {
        $gene_pfam{$_} = 1;
        $num_cds_short ++;
    }
    elsif ($cds_number < $ARGV[4]) {
        $gene_pfam{$_} = 1;
        $num_cds_lack ++;
    }
    else {
        $gene_ok{$_} = 1;
        $num_validate ++;
    }
}

print STDERR "The gene number statistics:
Total:               $num_total;
Genewise score < 15: $num_score_lt15;
CDS length < 90:     $num_cds_lt90;
Genewise score < $ARGV[2]: $num_score_low;
CDS length < $ARGV[3]:    $num_cds_short;
Exon number < $ARGV[4]:     $num_cds_lack;
High quality:        $num_validate;\n";

open IN, $ARGV[1] or die $!;
my (%seq, $seq_id);
while (<IN>) {
    chomp;
    if (/^>(\S+)/) { $seq_id = $1; }
    else { $seq{$seq_id} .= $_; }
}
close IN;

open FASTA, ">", "for_pfam_search.fasta" or die $!;
open COMMAND, ">", "command.hmmscan.list" or die $!;
mkdir "hmmscan.tmp" unless -e "hmmscan.tmp";
foreach my $gene_id (keys %gene_pfam) {
    my @gff = split /\n/, $gff{$gene_id};
    my @cds;
    foreach (@gff) {
        @_ = split /\t/;
        push @cds, "$_[0]\t$_[3]\t$_[4]\t$_[6]\t$_[7]" if m/\tcds\t/i;
    }
    @_ = split /\t/, $cds[0];
    my ($seqID, $strand, $frame) = ($_[0], $_[3], $_[4]);

    my $cds_seq;
    if ($strand eq "+") {
        foreach (@cds) {
            @_ = split /\t/;
            my $start = $_[1] - 1;
            my $len = $_[2] - $start;
            $cds_seq .= substr($seq{$seqID}, $start, $len);
        }
    }
    elsif ($strand eq "-") {
        my %sort;
        foreach (@cds) { @_ = split /\t/; $sort{$_} = $_[1]; }
        @cds = sort {$sort{$a} <=> $sort{$b}} @cds;
        foreach (@cds) {
            @_ = split /\t/;
            my $start = $_[2] - 1;
            my $len = $_[1] - $start;
            $cds_seq .= substr($seq{$seqID}, $start, $len);
        }
        $cds_seq = reverse $cds_seq;
        $cds_seq =~ tr/ATCGatcg/TAGCtagc/;
    }

    $cds_seq =~ s/^\w{$frame}//;
    #print ">$gene_id\n$cds_seq\n";
    my $pep_seq = &cds2pep($cds_seq, $gene_id);
    print FASTA ">$gene_id\n$pep_seq\n";

    open OUT, ">", "hmmscan.tmp/$gene_id.fasta" or die $!;
    print OUT ">$gene_id\n$pep_seq\n";
    close OUT;
    print COMMAND "hmmscan -o hmmscan.tmp/$gene_id.txt --cpu 1 -E $ARGV[5] --domE $ARGV[5] --tblout hmmscan.tmp/$gene_id.tbl --domtblout hmmscan.tmp/$gene_id.domtbl $Pfam_db hmmscan.tmp/$gene_id.fasta\n";
}
close FASTA;
close COMMAND;

my $cmdString = "ParaFly -c command.hmmscan.list -CPU $ARGV[7] &> command.hmmscan.log";
system ($cmdString) == 0 or die "Failed to execute: $cmdString\n$!";

my (%num_pfam_ok, $num_pfam_ok);
foreach my $gene_id (keys %gene_pfam) {
    open IN, "hmmscan.tmp/$gene_id.domtbl" or die $!;
    while (<IN>) {
        next if m/^#/;
        @_ = split /\s+/;
        my $ratio1 = abs($_[15] - $_[16]) / $_[2];
        my $ratio2 = abs($_[17] - $_[18]) / $_[5];
        if ($ratio1 >= $ARGV[6] or $ratio2 >= $ARGV[6]) {
            $gene_ok{$gene_id} = 1;
            $num_pfam_ok{$gene_id} = 1;
        }
    }
    close IN;
}
$num_pfam_ok = keys %num_pfam_ok;
print STDERR "Pfam Search passed:  $num_pfam_ok.\n";

foreach (@id) {
    print "$gff{$_}\n" if exists $gene_ok{$_};
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
    my $gene = shift @_;
    my $pep;
    while ((length $seq) >= 3) {
        $seq =~ s/(\w{3})//;
        $pep .= $cds2pep{$1};
    }
    #print STDERR "Warning: CDS length of $gene is not multiple of 3\n" if (length $seq) > 0;
    #print STDERR "Warning: Stop Codon appear in the middle of $gene\n" if $pep =~ m/\*\w/;
    return $pep;
}
