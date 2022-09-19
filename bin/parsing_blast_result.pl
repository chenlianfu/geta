#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    $0 [options] blast.out > blast.tab
    
    对BLAST的xml或tab格式的结果进行解析和过滤，得到更准确的BLAST结果。结果为表格形式（BLAST outfmt6），结果按query序列的ID排序，每个query序列的比对结果按得分排序。

    --type <string>    default: xml
    设置输入BLAST结果文件的类型。可以设置为xml或tab两种类型。
    若是tab格式，则BLAST结果中没有query与subject的序列长度信息，默认设置下无法使用--subject-coverage和--query-coverage参数的覆盖率阈值对结果进行过滤。在设置--db-subject输入数据库FASTA文件后可以使用--subject-coverage参数进行过滤；在设置--db-query输入query序列FASTA文件后可以使用--query-coverage参数进行过滤。
    若是xml格式，结果文件中包含query和subject长度信息，从而不需要使用--db-subject和--db-query参数输入FASTA序列文件。

    --no-header
    添加该参数则不输出表头。

    --max-hit-num <int>    default: 20
    设置允许的最大hit数量。

    --evalue <float>    default: 1e-5
    设置HSP的evalue阈值。

    --identity <float>    default: 0.05
    设置HSP的identity阈值。

    --CIP <float>    default: 0.2
    设置cumulative identity percentage阈值（这里依然使用了比值，单位不是%，所以其值要设置不大于1，默认值0.2表示20%阈值），对Hit进行过滤。CIP = 所有HSPs的一致位点之和 / 所有HSPs的比对长度之和。

    --subject-coverage <float>    default: 0.2
    设置所有HSPs对subject序列总体的覆盖率阈值。该参数阈值在文献中也被称为CALP（cumulative alignment length percentage），即 sum of all HSPs / subject length。

    --db-subject <string>
    输入数据库的FASTA文件，以获取subject序列长度信息。

    --query-coverage <float>    default: 0.2
    设置所有HSPs对query序列总体的覆盖率阈值。该参数阈值在文献中也被称为CALP（cumulative alignment length percentage），即 sum of all HSPs / query length。

    --db-query <string>
    输入query序列的FASTA文件，以获取query序列长度信息。

    --percentage-of-top-bitscore <int>    default: 100
    使用bitscore得分对hit进行过滤，设置输出hits的bitscore得分和最高得分相差不超过最高得分的百分数。hit若有多个HSPs，则取最高的HSP得分作为hit的得分；若数据库非常大，则推荐将设置该参数值设置为10，则能极大减少比对结果，保留最准确的结果；若数据库比较小，则推荐设置该参数值为50，或使用默认值；使用该参数来减少比对结果，优于仅使用最优比对结果。

    --HSP-num <int>    default: max
    若一个hit有多个HSPs，该参数设置输出得分指定数目个最高的HSPs。默认输出所有的HSPs。

    --out-hit-confidence
    添加该参数，则在表格结果第13、14和15列分别输出Hit的CIP、CALP_query、CALP_subject值。

    --suject-annotation
    若--type参数的值是xml，添加该参数可以生效，则额外增加最后一列suject annotation注释结果。


USAGE
if (@ARGV==0) {die $usage}

my ($type, $no_header, $maxHitNum, $evalue, $identity, $CIP, $subjectCoverage, $queryCoverage, $dbSubject, $dbQuery, $percentageOfTopBitscore, $HSPNum, $out_hit_confidence, $subjectAnnotation);
GetOptions(
    "type:s" => \$type,
    "no-header!" => \$no_header,
    "max-hit-num:i" => \$maxHitNum,
    "evalue:f" => \$evalue,
    "identity:f" => \$identity,
    "CIP:f" => \$CIP,
    "subject-coverage:f" => \$subjectCoverage,
    "query-coverage:f" => \$queryCoverage,
    "db-subject:s" => \$dbSubject,
    "db-query:s" => \$dbQuery,
    "percentage-of-top-bitscore:i" => \$percentageOfTopBitscore,
    "HSP-num:i" => \$HSPNum,
    "out-hit-confidence!" => \$out_hit_confidence,
    "suject-annotation!" => \$subjectAnnotation,
);
$type ||= "xml";
$maxHitNum ||= 20;
$evalue ||= 1e-5;
$identity ||= 0.05;
$CIP ||= 0.2;
$subjectCoverage ||= 0.2;
$queryCoverage ||= 0.2;
$percentageOfTopBitscore ||= 100;

my $header = "#QueryID\tSubjectID\tIdentity\tMatchLength\tMismatchLength\tGaps\tQueryStart\tQueryEnd\tSubjectStart\tSubjectEnd\tEvalue\tBitScore";
$header .= "\tCIP\tCALP_query\tCALP_subject" if $out_hit_confidence;
$header .= "\tSubjectAnnotation" if ($type eq "xml" && $subjectAnnotation);
print "$header\n" unless $no_header;

my (%seqLength_subject, %seqLength_query, $seq_id);
if ($type eq "xml") {
    open IN, $ARGV[0] or die "Can not open file $ARGV[0], $!\n";
    my $one;
    while (<IN>) {
        $one .= $_;
        if (/<\/Iteration>/) {
            my $query_id = $1 if $one =~ m#<Iteration_query-def>(.+?)</Iteration_query-def>#;
            my $query_length = $1 if $one =~ m#<Iteration_query-len>(\d+?)</Iteration_query-len>#;

            my @hit = $one =~ m#<Hit>(.*?)</Hit>#gs;
            unless (@hit) {
                $one = "";
                #print STDERR "Filtered by Nohit: $query_id\n";
                next;
            }
            my $hitNum = 0;

            my (%hit_score, %hit_hsp);
            # 解析hit的得分和其HSP信息
            foreach (@hit) {
                my $subject_id = $1 if $_ =~ m#<Hit_id>(.+?)</Hit_id>#;
                my $subject_length = $1 if $_ =~ m#<Hit_len>(\d+?)</Hit_len>#;
                my $subject_annotation = $1 if $_ =~ m#<Hit_def>(.*?)</Hit_def>#;

                my (@hsp_query, @hsp_subject, %hsp, @score, $hit_total_identity, $hit_total_align_len);
                my @hsp = $_ =~ m#<Hsp>(.*?)</Hsp>#gs;
                next unless @hsp;
                foreach (@hsp) {
                    my $hsp_evalue = $1 if $_ =~ m#<Hsp_evalue>(\S+?)</Hsp_evalue>#;
                    next if $hsp_evalue > $evalue;

                    my $hsp_identity = $1 if $_ =~ m#<Hsp_identity>(\d+?)</Hsp_identity>#;
                    my $hsp_identity_num = $hsp_identity;
                    my $hsp_align_len = $1 if $_ =~ m#<Hsp_align-len>(\d+?)</Hsp_align-len>#;
                    my $hsp_gaps = $1 if $_ =~ m#<Hsp_gaps>(\d+?)</Hsp_gaps>#;
                    my $hsp_mismatch = $hsp_align_len - $hsp_identity - $hsp_gaps;
                    $hsp_identity = int(($hsp_identity / $hsp_align_len) * 10000 + 0.5) / 10000;
                    next if $hsp_identity < $identity;
                    $hsp_identity = $hsp_identity * 100;
                    $hit_total_identity += $hsp_identity_num;
                    $hit_total_align_len += $hsp_align_len;

                    my $query_start = $1 if $_ =~ m#<Hsp_query-from>(\d+?)</Hsp_query-from>#;
                    my $query_end = $1 if $_ =~ m#<Hsp_query-to>(\d+?)</Hsp_query-to>#;
                    my $subject_start = $1 if $_ =~ m#<Hsp_hit-from>(\d+?)</Hsp_hit-from>#;
                    my $subject_end = $1 if $_ =~ m#<Hsp_hit-to>(\d+?)</Hsp_hit-to>#;
                    my $hsp_score = $1 if $_ =~ m#<Hsp_bit-score>(\S+?)</Hsp_bit-score>#;
                    push @score, $hsp_score;

                    $hsp{"$query_id\t$subject_id\t$hsp_identity\t$hsp_align_len\t$hsp_mismatch\t$hsp_gaps\t$query_start\t$query_end\t$subject_start\t$subject_end\t$hsp_evalue\t$hsp_score"} = $hsp_score;

                    my $query_from = $1 if m#<Hsp_query-from>(\d+?)</Hsp_query-from>#;
                    my $query_to = $1 if m#<Hsp_query-to>(\d+?)</Hsp_query-to>#;
                    my $hit_from = $1 if m#<Hsp_hit-from>(\d+?)</Hsp_hit-from>#;
                    my $hit_to = $1 if m#<Hsp_hit-to>(\d+?)</Hsp_hit-to>#;
                    push @hsp_query, "$query_from\t$query_to";
                    push @hsp_subject, "$hit_from\t$hit_to";
                }
                #print STDERR "Filtered by Identity: $query_id\t$subject_id\n" unless @score;
                @score = sort {$b <=> $a} @score;
                my $match_length_query = &match_length(@hsp_query);
                my $match_length_subject = &match_length(@hsp_subject);
                my $match_coverage_query = $match_length_query / $query_length;
                my $match_coverage_subject = $match_length_subject / $subject_length;
                my $cip = 0;
                $cip = $hit_total_identity / $hit_total_align_len if $hit_total_identity;
                $match_coverage_query = int($match_coverage_query * 10000 + 0.5) / 100;
                $match_coverage_subject = int($match_coverage_subject * 10000 + 0.5) / 100;
                $cip = int($cip * 10000 + 0.5) / 100;
                #print STDERR "$query_id\t$subject_id\t$hit_total_identity\t$hit_total_align_len\t$cip\n";

                if ($match_coverage_query >= ($queryCoverage * 100) && $match_coverage_subject >= ($subjectCoverage * 100) && $cip >= ($CIP * 100)) {
                    $hit_score{$subject_id} = $score[0];
                    $hitNum ++;
                    last if $hitNum > $maxHitNum;
                    
                    foreach (keys %hsp) {
                        $_ .= "\t$cip\t$match_coverage_query\t$match_coverage_subject" if $out_hit_confidence;
                        $_ .= "\t$subject_annotation" if $subjectAnnotation;
                        $_ .= "\n";
                        $hit_hsp{$subject_id}{$_} = $hsp{$_};
                    }
                }
                else {
                    #print STDERR "Filtered by Coverage: $query_id\t$subject_id\t$match_coverage_query\t$match_coverage_subject\t$match_length_subject\t$subject_length\n";
                }
            }

            my @subject_id = sort {$hit_score{$b} <=> $hit_score{$a}} keys %hit_score;
            my $highest_score = $hit_score{$subject_id[0]};
            my $lowest_score = $highest_score * ( 1 - $percentageOfTopBitscore / 100 );
            #print STDERR "$query_id\t$subject_id[0]\t$highest_score\t$lowest_score\n";

            foreach my $id (@subject_id) {
                last if $hit_score{$id} < $lowest_score;
                my $hsp_num = 0;
                foreach (sort {$hit_hsp{$id}{$b} <=> $hit_hsp{$id}{$a}} keys %{$hit_hsp{$id}}) {
                    print;
                    $hsp_num ++;
                    last if ($HSPNum && ($hsp_num >= $HSPNum));
                }
            }

            $one = "";
        }
    }
}
elsif ($type eq "tab") {
    #my (%seqLength_subject, %seqLength_query, $seq_id);
    if ($dbSubject) {
        open IN, $dbSubject or die "Can not open file $dbSubject, $!\n";
        while (<IN>) {
            chomp;
            if (m/^>(\S+)/) { $seq_id = $1; }
            else { $seqLength_subject{$seq_id} += length($_); }
        }
        close IN;
    }
    if ($dbQuery) {
        open IN, $dbQuery or die "Can not open file $dbQuery, $!\n";
        while (<IN>) {
            chomp;
            if (m/^>(\S+)/) { $seq_id = $1; }
            else { $seqLength_query{$seq_id} += length($_); }
        }
    }

    open IN, $ARGV[0] or die "Can not open file $ARGV[0], $!\n";
    my ($one, $geneID);
    while (<IN>) {
        next if m/^#/;
        if (m/^(\S+?)\t/) {
            my $geneID_new = $1;

            if ($geneID_new eq $geneID) {
                $one .= $_;
            }
            else {
                #&parsing_gene_annot($one, \%seqLength_subject, \%seqLength_query) if $one;
                &parsing_gene_annot($one) if $one;
                $one = $_;
                $geneID = $geneID_new;
            }
        }
    }
    #&parsing_gene_annot($one, \%seqLength_subject, \%seqLength_query) if $one;
    &parsing_gene_annot($one) if $one;
}
else {
    die "The parameter --type should be set to xml or tab.\n";
}

sub parsing_gene_annot {
    #my %seqLength_subject = %{$_[1]};
    #my %seqLength_query = %{$_[2]};
    my @line = split /\n/, $_[0];
    my (%hit_hsp_raw, %hit_score_raw);
    foreach (@line) {
        next if m/^#/;
        @_ = split /\t/, $_;
        next if ($_[2] / 100) < $identity;
        next if $_[10] > $evalue;
        $hit_hsp_raw{$_[0]}{$_[1]}{$_} = $_[11];
        $hit_score_raw{$_[0]}{$_[1]} = $_[11] if $hit_score_raw{$_[0]}{$_[1]} < $_[11];

    }

    my (%hit_hsp, %hit_score);
    foreach my $query_id (keys %hit_score_raw) {
        my $hitNum = 0;
        foreach my $subject_id (sort {$hit_score_raw{$query_id}{$b} <=> $hit_score_raw{$query_id}{$a}} keys %{$hit_score_raw{$query_id}}) {
            my @hsp = keys %{$hit_hsp_raw{$query_id}{$subject_id}};
            my (@hsp_query, @hsp_subject, $hit_total_identity, $hit_total_align_len);
            foreach (@hsp) {
                @_ = split /\t/; 
                push @hsp_query, "$_[6]\t$_[7]";
                push @hsp_subject, "$_[8]\t$_[9]"; 
                $hit_total_identity += $_[2] / 100 * $_[3];
                $hit_total_align_len += $_[3];

            }

            my $match_length_query = &match_length(@hsp_query);
            my $match_length_subject = &match_length(@hsp_subject);
            my ($cip, $match_coverage_query, $match_coverage_subject) = (0, "-", "-"); 
            $cip = int($hit_total_identity / $hit_total_align_len * 10000 + 0.5) / 100 if $hit_total_identity;

            my $keep = 1;
            if (exists $seqLength_query{$query_id}) {
                $match_coverage_query = int($match_length_query / $seqLength_query{$query_id} * 10000 + 0.5) / 100;
                $keep = 0 if ($match_coverage_query < $queryCoverage * 100);
            }
            if (exists $seqLength_subject{$subject_id}) {
                $match_coverage_subject = int($match_length_subject / $seqLength_subject{$subject_id} * 10000 + 0.5) / 100;
                $keep = 0 if ($match_coverage_subject < $subjectCoverage * 100);
            }
            $keep = 0 if $cip < $CIP * 100;

            if ($keep == 1) {
                foreach (@hsp) {
                    my $out = $_;
                    $out .= "\t$cip\t$match_coverage_query\t$match_coverage_subject" if $out_hit_confidence;
                    $out .= "\n";
                    $hit_hsp{$query_id}{$subject_id}{$out} = $hit_hsp_raw{$query_id}{$subject_id}{$_};
                }
                $hit_score{$query_id}{$subject_id} = $hit_score_raw{$query_id}{$subject_id};
                $hitNum ++;
                last if $hitNum > $maxHitNum;
            }
        }
    }

    foreach my $query_id (sort keys %hit_score){
        my @subject_id = sort {$hit_score{$query_id}{$b} <=> $hit_score{$query_id}{$a}} keys %{$hit_score{$query_id}};
        my $highest_score = $hit_score{$query_id}{$subject_id[0]};
        my $lowest_score = $highest_score * ( 1 - $percentageOfTopBitscore / 100 );

        foreach my $subject_id (@subject_id) {
            last if $hit_score{$query_id}{$subject_id} < $lowest_score;
            my $hsp_num = 0;
            foreach (sort {$hit_hsp{$query_id}{$subject_id}{$b} <=> $hit_hsp{$query_id}{$subject_id}{$a}} keys %{$hit_hsp{$query_id}{$subject_id}}) {
                print;
                $hsp_num ++;
                last if ($HSPNum && ($hsp_num >= $HSPNum));
            }
        }
    }
}

sub match_length {
    my @inter_sorted_site;
    foreach (@_) {
        my @aaa = $_ =~ m/(\d+)/g;
        @aaa = sort { $a <=> $b } @aaa;
        push @inter_sorted_site, "$aaa[0]\t$aaa[1]";
    }
    @inter_sorted_site = sort { $a <=> $b } @inter_sorted_site;

    my $out_site_number;
    my $former_region = shift @inter_sorted_site;
    my @aaa = $former_region =~ m/(\d+)/g;
    $out_site_number += ($aaa[1] - $aaa[0] + 1);
    foreach (@inter_sorted_site) {
        my @former_region = $former_region =~ m/(\d+)/g;
        my @present_region = $_ =~ m/(\d+)/g;
        
        if ($present_region[0] > $former_region[1]) {
            $out_site_number += ($present_region[1] - $present_region[0] + 1);
            $former_region = $_;
        }
        elsif ($present_region[1] > $former_region[1]) {
            $out_site_number += ($present_region[1] - $former_region[1]);
            $former_region = $_;
        }
        else {
            next
        }
    }
    return $out_site_number;
}
