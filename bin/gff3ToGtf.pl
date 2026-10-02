#!/usr/bin/perl
use strict;
use warnings;

my $usage = <<USAGE;
Usage:
    perl $0 genome.fasta file.gff3 > file.gtf

Convert a GFF3 gene-model file into GTF, adding start_codon /
stop_codon features inferred from the CDS coordinates and the
genome sequence.
USAGE

die $usage unless @ARGV == 2;
my ($fasta_file, $gff_file) = @ARGV;

# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------

# reverse complement
sub revcom {
    my ($seq) = @_;
    $seq = reverse $seq;
    $seq =~ tr/ATCGatcg/TAGCtagc/;
    return $seq;
}

# pull a "key=value" out of a GFF3 column-9 attribute string
sub get_attr {
    my ($attr, $key) = @_;
    return $attr =~ /\b\Q$key\E=([^;\n]+)/ ? $1 : undef;
}

# ----------------------------------------------------------------------
# 1. load genome sequence
# ----------------------------------------------------------------------
my %fasta;
{
    open my $fh, '<', $fasta_file or die "ERROR: cannot open $fasta_file: $!\n";
    my $id;
    while (<$fh>) {
        chomp;
        if (/^>(\S+)/) { $id = $1; $fasta{$id} = ''; }
        elsif (defined $id) { $fasta{$id} .= uc($_); }
    }
    close $fh;
}

# ----------------------------------------------------------------------
# 2. load GFF3
#    - %gene       : gene_id -> { mRNA_id -> 1 }
#    - %gff        : mRNA_id -> concatenated child feature lines
#    - %strand_of  : mRNA_id -> strand, read directly off the mRNA line
#                    (the original script inferred this from whatever
#                    was left in @_ after a previous loop -- fragile and
#                    wrong whenever a transcript had no CDS at all)
# ----------------------------------------------------------------------
my (%gene, %gff, %strand_of);
{
    open my $fh, '<', $gff_file or die "ERROR: cannot open $gff_file: $!\n";
    while (<$fh>) {
        chomp;
        next if /^#/ || /^\s*$/;
        my @col = split /\t/;
        next if @col < 9;

        my ($type, $strand, $attr) = @col[2, 6, 8];

        if ($type eq 'mRNA') {
            my $mRNA_id = get_attr($attr, 'ID');
            my $parent  = get_attr($attr, 'Parent');
            if (!defined $mRNA_id || !defined $parent) {
                warn "WARN: mRNA line missing ID/Parent, skipped:\n$_\n";
                next;
            }
            $gene{$_}{$mRNA_id} = 1 for split /,/, $parent;
            $strand_of{$mRNA_id} = $strand;
            next; # mRNA line itself is not a child feature we need to keep
        }

        my $parent = get_attr($attr, 'Parent');
        next unless defined $parent;
        my $line = join("\t", @col);
        $gff{$_} .= "$line\n" for split /,/, $parent;
    }
    close $fh;
}

# ----------------------------------------------------------------------
# 3. compute start_codon / stop_codon lines for one mRNA
#    returns a list of GTF lines (no gene_id/transcript_id yet - that is
#    added by the caller together with all other feature lines)
# ----------------------------------------------------------------------
sub build_codons {
    my ($gene_id, $mRNA_id, $strand, $cds_ref, $fasta_ref) = @_;
    my @cds = @$cds_ref;

    my @first  = split /\t/, $cds[0];
    my @four   = split /\t/, $cds[-1];
    my @second = @cds > 1 ? split(/\t/, $cds[1])  : ();
    my @three  = @cds > 1 ? split(/\t/, $cds[-2]) : @first;

    my $tail = "gene_id \"$gene_id\"; transcript_id \"$mRNA_id\";";
    my @out;

    if ($strand eq '+') {
        # ---- start codon: lowest-coordinate CDS segment ----
        my $length = $first[4] - $first[3] + 1;
        my ($start_codon, $start_bases);
        if ($length >= 3) {
            my ($s, $e) = ($first[3], $first[3] + 2);
            $start_codon = join("\t", $first[0], $first[1], 'start_codon', $s, $e, '.', $strand, 0, $tail);
            $start_bases = substr($fasta_ref->{$first[0]}, $s - 1, 3);
        }
        elsif (@second) {
            my $e = $second[3] + (3 - $length - 1);
            $start_codon = join("\t", $first[0], $first[1], 'start_codon', $first[3], $first[4], '.', $strand, 0, $tail)
                . "\n" . join("\t", $first[0], $first[1], 'start_codon', $second[3], $e, '.', $strand, $length, $tail);
            $start_bases  = substr($fasta_ref->{$first[0]}, $first[3] - 1, $length);
            $start_bases .= substr($fasta_ref->{$first[0]}, $second[3] - 1, 3 - $length);
        }
        else {
            warn "WARN: $mRNA_id (+) start_codon spans a splice but only one CDS segment exists, skipped\n";
        }
        if (defined $start_bases) {
            if ($start_bases eq 'ATG') { push @out, $start_codon; }
            else { warn "WARN: $mRNA_id (+) start_codon bases '$start_bases' != ATG, skipped\n"; }
        }

        # ---- stop codon: highest-coordinate CDS segment ----
        $length = $four[4] - $four[3] + 1;
        my ($stop_codon, $stop_bases);
        if ($length >= 3) {
            my ($s, $e) = ($four[4] - 2, $four[4]);
            $stop_codon = join("\t", $four[0], $four[1], 'stop_codon', $s, $e, '.', $strand, 0, $tail);
            $stop_bases = substr($fasta_ref->{$four[0]}, $s - 1, 3);
        }
        elsif (@cds > 1) {
            my $frame = 3 - $length;
            my $s = $three[4] - (3 - $length - 1);
            $stop_codon = join("\t", $four[0], $four[1], 'stop_codon', $four[3], $four[4], '.', $strand, $frame, $tail)
                . "\n" . join("\t", $four[0], $four[1], 'stop_codon', $s, $three[4], '.', $strand, 0, $tail);
            $stop_bases  = substr($fasta_ref->{$four[0]}, $s - 1, 3 - $length);
            $stop_bases .= substr($fasta_ref->{$four[0]}, $four[3] - 1, $length);
        }
        else {
            warn "WARN: $mRNA_id (+) stop_codon spans a splice but only one CDS segment exists, skipped\n";
        }
        if (defined $stop_bases) {
            if ($stop_bases =~ /^(TAA|TAG|TGA)$/) { push @out, $stop_codon; }
            else { warn "WARN: $mRNA_id (+) stop_codon bases '$stop_bases' not a stop codon, skipped\n"; }
        }
    }
    elsif ($strand eq '-') {
        # ---- stop codon: lowest-coordinate CDS segment ----
        my $length = $first[4] - $first[3] + 1;
        my ($stop_codon, $stop_bases);
        if ($length >= 3) {
            my ($s, $e) = ($first[3], $first[3] + 2);
            # (original script had a stray leading "\n" on this line, which
            #  produced a spurious blank line in the GTF output -- removed)
            $stop_codon = join("\t", $first[0], $first[1], 'stop_codon', $s, $e, '.', $strand, 0, $tail);
            $stop_bases = substr($fasta_ref->{$first[0]}, $s - 1, 3);
        }
        elsif (@second) {
            my $frame = 3 - $length;
            my $e = $second[3] + (3 - $length - 1);
            $stop_codon = join("\t", $first[0], $first[1], 'stop_codon', $first[3], $first[4], '.', $strand, $frame, $tail)
                . "\n" . join("\t", $first[0], $first[1], 'stop_codon', $second[3], $e, '.', $strand, 0, $tail);
            $stop_bases  = substr($fasta_ref->{$first[0]}, $first[3] - 1, $length);
            $stop_bases .= substr($fasta_ref->{$first[0]}, $second[3] - 1, 3 - $length);
        }
        else {
            warn "WARN: $mRNA_id (-) stop_codon spans a splice but only one CDS segment exists, skipped\n";
        }
        if (defined $stop_bases) {
            my $b = revcom($stop_bases);
            if ($b =~ /^(TAA|TAG|TGA)$/) { push @out, $stop_codon; }
            else { warn "WARN: $mRNA_id (-) stop_codon bases '$b' not a stop codon, skipped\n"; }
        }

        # ---- start codon: highest-coordinate CDS segment ----
        $length = $four[4] - $four[3] + 1;
        my ($start_codon, $start_bases);
        if ($length >= 3) {
            my ($s, $e) = ($four[4] - 2, $four[4]);
            $start_codon = join("\t", $four[0], $four[1], 'start_codon', $s, $e, '.', $strand, 0, $tail);
            $start_bases = substr($fasta_ref->{$four[0]}, $s - 1, 3);
        }
        elsif (@cds > 1) {
            my $s = $three[4] - (3 - $length - 1);
            $start_codon = join("\t", $four[0], $four[1], 'start_codon', $four[3], $four[4], '.', $strand, 0, $tail)
                . "\n" . join("\t", $four[0], $four[1], 'start_codon', $s, $three[4], '.', $strand, $length, $tail);
            $start_bases  = substr($fasta_ref->{$four[0]}, $s - 1, 3 - $length);
            $start_bases .= substr($fasta_ref->{$four[0]}, $four[3] - 1, $length);
        }
        else {
            warn "WARN: $mRNA_id (-) start_codon spans a splice but only one CDS segment exists, skipped\n";
        }
        if (defined $start_bases) {
            my $b = revcom($start_bases);
            if ($b eq 'ATG') { push @out, $start_codon; }
            else { warn "WARN: $mRNA_id (-) start_codon bases '$b' != ATG, skipped\n"; }
        }
    }
    else {
        warn "WARN: $mRNA_id has unknown strand '$strand', start/stop codon skipped\n";
    }

    return @out;
}

# ----------------------------------------------------------------------
# 4. emit GTF, gene by gene / transcript by transcript
# ----------------------------------------------------------------------
my %rank = ('5UTR' => 1, start_codon => 2, exon => 3, CDS => 3, stop_codon => 4, '3UTR' => 5);

for my $gene_id (sort keys %gene) {
    for my $mRNA_id (sort keys %{ $gene{$gene_id} }) {

        my $content = $gff{$mRNA_id};
        if (!defined $content) {
            warn "WARN: $mRNA_id has no child features (exon/CDS/UTR), skipped\n";
            next;
        }
        my @lines = split /\n/, $content;

        my @utr = grep { /UTR\t/ } @lines;
        for (@utr) {
            s/\tfive_prime_UTR\t/\t5UTR\t/;
            s/\tthree_prime_UTR\t/\t3UTR\t/;
        }
        my @exon = grep { /\texon\t/ } @lines;
        my @cds  = grep { /\tCDS\t/ }  @lines;

        my $strand = $strand_of{$mRNA_id};
        if (!defined $strand) {
            warn "WARN: $mRNA_id has no recorded strand, skipped\n";
            next;
        }

        # sort CDS ascending by genomic start coordinate
        @cds = sort { (split /\t/, $a)[3] <=> (split /\t/, $b)[3] } @cds;

        my @output = (@utr, @exon, @cds);
        if (@cds) {
            push @output, build_codons($gene_id, $mRNA_id, $strand, \@cds, \%fasta);
        }
        else {
            warn "WARN: $mRNA_id has no CDS, no start/stop codon computed\n";
        }

        # tag each line with (sort_rank, genomic_start) instead of using the
        # line text itself as a hash key, then sort 5'UTR -> start_codon ->
        # exon/CDS -> stop_codon -> 3'UTR, ascending on + strand and
        # descending on - strand
        my @tagged;
        for my $line (@output) {
            my @f = split /\t/, $line;
            next unless @f >= 8;
            my $r = $rank{ $f[2] };
            next unless defined $r;
            push @tagged, [ $line, $r, $f[3] ];
        }
        @tagged = $strand eq '-'
            ? sort { $a->[1] <=> $b->[1] or $b->[2] <=> $a->[2] } @tagged
            : sort { $a->[1] <=> $b->[1] or $a->[2] <=> $b->[2] } @tagged;

        for my $t (@tagged) {
            (my $line = $t->[0]) =~ s/^(.*)\t[^\t]*$/$1\tgene_id "$gene_id"; transcript_id "$mRNA_id"; gene_name "$gene_id";/;
            print "$line\n";
        }
    }
    print "\n";
}
