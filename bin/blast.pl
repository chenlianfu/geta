#!/usr/bin/perl
use strict;
use Getopt::Long;
use Cwd qw/abs_path getcwd cwd/;
use File::Basename;

my $usage = <<USAGE;
Usage:
    $0 [options] BLAST_DB file.fasta > out.txt

    --tmp-prefix <string>    default: blast
    设置临时文件或文件夹前缀。默认设置下，程序生成command.blast.list，blast.tmp/等临时文件或目录。

    --chunk <int>    default: 10
    设置每个数据块的序列条数。程序会将输入FASTA文件中的序列从前往后分割成多份，每10条相邻的序列分配到一个FASTA文件中；在blast.tmp/临时文件夹下生成次级文件夹，每个文件夹做多放置10个FASTA文件；每个fasta文件写出一条BLAST命令到command.blast.list文件中；然后程序调用ParaFly进行并行化计算。
    请注意：若数据块的数量超过100万个，默认设置下blast.tmp/文件夹中的目录数量太多（超过1万个），导致文件系统运行缓慢，ParaFly程序运行效率低下，无法充分利用服务器计算资源。此时推荐设置--chunk参数值为100。

    --blast-program <string>    default: blastp
    设置运行的BLAST命令，支持的命令有：blastn, blastp, blastx, tblastn, tblastx。

    --CPU <int>    default: 1
    设置并行运行的BLAST程序个数。

    --blast-threads <int>    default: 1
    设置BLAST命令的-num_threads参数值。该参数让每个BLAST命令可以多线程运行。
    请注意：--blast-threads参数值和--CPU参数值的乘积不要超过服务器的CPU总计算线程数。

    --evalue <float>    default: 1e-3
    设置BLAST命令的-evalue参数值。

    --outfmt <int>    default: 5
    设置BLAST命令的-outfmt参数值。输出方式。若为5，则输出xml格式结果，若为6或7，则输出表格结果。

    --max-target-seqs <int>    default: 20
    设置BLAST命令的-max_target_seqs参数值。该参数设置BLAST最多能匹配数据库中的序列数量。

    --completed_ratio <float>    default: 1
    设置程序要求的最低完成度。程序进行了数据分割，得到多个命令，要求对这些命令完成的比例要不小于该值，否则程序会运行失败。默认设置下要求所有的序列都完成比对。当本程序被集成到其它流程中时，可能个别命令执行失败对总体影响不大，考虑到整体的问题运行，考虑将该参数设置小于1，以让整个流程能顺利运行。

    --clean
    若添加该参数，则在运行程序成功后，会删除临时文件或文件夹。

USAGE
if (@ARGV==0){die $usage}

my ($tmpPrefix, $chunk, $blastProgram, $CPU, $blastThreads, $evalue, $outfmt, $maxTargetSeqs, $completed_ratio, $clean);
GetOptions(
    "tmp-prefix:s" => \$tmpPrefix,
    "chunk:i" => \$chunk,
    "blast-program:s" => \$blastProgram,
    "CPU:i" => \$CPU,
    "blast-threads:i" => \$blastThreads,
    "evalue:f" => \$evalue,
    "outfmt:i" => \$outfmt,
    "max-target-seqs:i" => \$maxTargetSeqs,
    "completed_ratio:f" => \$completed_ratio,
    "clean!" => \$clean,
);
$tmpPrefix ||= "blast";
$tmpPrefix = abs_path($tmpPrefix);
$chunk ||= 10;
$blastProgram ||= "blastp";
$CPU ||= 1;
$blastThreads ||= 1;
$evalue ||= 1e-3;
$outfmt ||= 5;
$maxTargetSeqs ||= 20;
$completed_ratio ||= 1;

my %blastProgram = ("blastn", 1, "blastp", 1, "blastx", 1, "tblastn", 1, "tblastx", 1);
if (! exists $blastProgram{$blastProgram}) {
    die "$blastProgram was not a supportted BLAST command (blastn, blastp, blastx, tblastn, tblastx)!\n";
}

mkdir "$tmpPrefix.tmp" unless -e "$tmpPrefix.tmp";
open IN, $ARGV[1] or die "Can not open file $ARGV[1], $!\n";
open CMD, ">", "command.$tmpPrefix.list" or die "Cannot create file command.$tmpPrefix.list, $!\n";

my ($fasta, @chunk);
my ($number_L1, $number_L2, $chunk_L1_number, $chunk_L2_number) = (1, 0, 1, 0);
my $chunk_L1_name = "chunk" . '0' x (5 - length($chunk_L1_number)) . $chunk_L1_number;
mkdir "$tmpPrefix.tmp/$chunk_L1_name" unless -e "$tmpPrefix.tmp/$chunk_L1_name";
while (<IN>) {
    if (m/^>/) {
        $number_L2 ++;

        if ($number_L2 > $chunk) {
            $chunk_L2_number ++;
            my $chunk_L2_name = '0' x (length($chunk) - length($chunk_L2_number)) . $chunk_L2_number;
            push @chunk, "$tmpPrefix.tmp/$chunk_L1_name/$chunk_L2_name";
            open OUT, ">", "$tmpPrefix.tmp/$chunk_L1_name/$chunk_L2_name.fasta" or die "Can not create file $tmpPrefix.tmp/$chunk_L1_name/$chunk_L2_name.fasta, $!\n";
            print OUT $fasta;
            print CMD "$blastProgram -query $tmpPrefix.tmp/$chunk_L1_name/$chunk_L2_name.fasta -db $ARGV[0] -num_threads $blastThreads -evalue $evalue -outfmt $outfmt -max_target_seqs $maxTargetSeqs -out $tmpPrefix.tmp/$chunk_L1_name/$chunk_L2_name.out\n";
            close OUT;
            $number_L2 = 1;
            $number_L1 ++;
            $fasta = "";
        }

        if ($number_L1 > $chunk) {
            $chunk_L1_number ++;
            $chunk_L1_name = "chunk" . '0' x (5 - length($chunk_L1_number)) . $chunk_L1_number;
            mkdir "$tmpPrefix.tmp/$chunk_L1_name" unless -e "$tmpPrefix.tmp/$chunk_L1_name";
            $number_L1 = 1;
            $chunk_L2_number = 0;
        }
    }
    $fasta .= $_;
}
if ($fasta) {
    $chunk_L2_number ++;
    my $chunk_L2_name = '0' x (length($chunk) - length($chunk_L2_number)) . $chunk_L2_number;
    push @chunk, "$tmpPrefix.tmp/$chunk_L1_name/$chunk_L2_name";
    open OUT, ">", "$tmpPrefix.tmp/$chunk_L1_name/$chunk_L2_name.fasta" or die "Can not create file $tmpPrefix.tmp/$chunk_L1_name/$chunk_L2_name.fasta, $!\n";
    print OUT $fasta;
    print CMD "$blastProgram -query $tmpPrefix.tmp/$chunk_L1_name/$chunk_L2_name.fasta -db $ARGV[0] -num_threads $blastThreads -evalue $evalue -outfmt $outfmt -max_target_seqs $maxTargetSeqs -out $tmpPrefix.tmp/$chunk_L1_name/$chunk_L2_name.out\n";
    close OUT;
}
close CMD;

my $cmdString = "ParaFly -c command.$tmpPrefix.list -CPU $CPU &> /dev/null";
print STDERR "CMD: $cmdString\n";
(system $cmdString) == 0 or warn "Warning: Failed to execute: $cmdString\n";

my ($cmd_num, $completed_num);
open IN, "command.$tmpPrefix.list" or die "Can not open file command.$tmpPrefix.list, $!";
while (<IN>) {
    $cmd_num ++ if m/\S+/;
}
close IN;
open IN, "command.$tmpPrefix.list.completed" or die "Can not open file command.$tmpPrefix.list.completed, $!";
while (<IN>) {
    $completed_num ++ if m/\S+/;
}
if ( $completed_num / $cmd_num < $completed_ratio ) {
    die "程序将$blastProgram任务分割成 $cmd_num 份，目前仅完成了 $completed_num 份，未能最低的完成比例需求 $completed_ratio。\n";
}

foreach (@chunk) {
    if ($outfmt == 5) {
        open IN, "$_.out" or die "Can not open file $_.out, $!\n";
        my $start .= <IN>;
        foreach (1 .. 19) {
            $start .= <IN>;
        }
        my @lines = <IN>;
        my $end = pop @lines;
        $end = (pop @lines) . $end;
        $end = (pop @lines) . $end;
        #print "$start$end";

        my (@out, $one_interation);
        foreach (@lines) {
            $one_interation .= $_;
            if (m/^<\/Iteration>/) {
                $one_interation =~ s#<Iteration_iter-num>\d+</Iteration_iter-num>#<Iteration_iter-num>1</Iteration_iter-num>#;
                $one_interation =~ s#<Iteration_query-ID>.*?</Iteration_query-ID>#<Iteration_query-ID>Query_1</Iteration_query-ID>#;
                my $query_def = $1 if $one_interation =~ m#<Iteration_query-def>(.*?)</Iteration_query-def>#;
                my $query_len = $1 if $one_interation =~ m#<Iteration_query-len>(\d+)</Iteration_query-len>#;
                my $start_new = $start;
                $start_new =~ s#<BlastOutput_query-def>.*?</BlastOutput_query-def>#<BlastOutput_query-def>$query_def</BlastOutput_query-def>#;
                $start_new =~ s#<BlastOutput_query-len>\d+</BlastOutput_query-len>#<BlastOutput_query-len>$query_len</BlastOutput_query-len>#;
                push @out, "$start_new$one_interation$end";
                $one_interation = "";
            }
        }

        foreach (@out) {
            print;
        }
        close IN;
    }
    else {
        open IN, "$_.out" or die "Can not open file $_.out, $!\n";
        print <IN>;
        close IN;
    }
}

if ($clean) {
    unlink "command.$tmpPrefix.list";
    unlink "command.$tmpPrefix.list.completed";
    system `rm -rf $tmpPrefix.tmp`;
}
