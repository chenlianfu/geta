#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl homolog_paraBlast blastProgramType blastDB fastaFile evalue threads outPrefix outfmt
	7 parameters should be given, and the final result is outPrefix.xml or outPrefix.tab

For example:
    perl homolog_paraBlast blastp nr proteins.fa 1e-3 24 nr 5

USAGE
if (@ARGV != 7){die $usage}

my ($program, $db, $fasta, $evalue, $threads, $outPrefix, $outfmt) = @ARGV;
open IN, $fasta or die $!;
my ($seqID, @seqID, %seq);
while (<IN>) {
    chomp;
    if (/^>(\S+)/) { $seqID = $1; push @seqID, $seqID; }
    else           { $seq{$seqID} .= $_ }
}

mkdir "$outPrefix.tmp" unless -e "$outPrefix.tmp";
open COM, '>', "$outPrefix.command" or die $!;

foreach my $seq_id (@seqID) {
    $_ = $seq_id;
    s/\|/\\\|/;
    open FASTA, '>', "$outPrefix.tmp/$seq_id.fa" or die $!;
    print FASTA ">$seq_id\n$seq{$seq_id}\n";
    close FASTA;

    print COM "$program -query $outPrefix.tmp/$_.fa -db $db -evalue $evalue -num_threads 1 -outfmt $outfmt -out $outPrefix.tmp/$_.out\n";
}

system `ParaFly -c $outPrefix.command -CPU $threads`;

unless ( -e "FailedCommands" ) {
    if ($outfmt == 5) { open OUT, '>', "$outPrefix.xml" or die $! }
    else              { open OUT, '>', "$outPrefix.tab" or die $! }

    foreach (@seqID) {
        open IN, "$outPrefix.tmp/$_.out" or die $!;
        print OUT <IN>;
    }

    system `rm -rf $outPrefix.tmp`;
    unlink "$outPrefix.command";
    unlink "$outPrefix.command.completed";
}
