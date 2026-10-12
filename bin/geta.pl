#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename qw(basename dirname);
use File::Path qw(make_path remove_tree);
use File::Spec;
use Cwd qw(abs_path getcwd);
use POSIX qw(strftime);

$| = 1;

########################################################################
# 中文用法说明（放在代码最前面；英文用法说明放在文件最末尾的__DATA__部分）
########################################################################
my $usage_cn = <<'USAGE';
用法:
    perl __PROG__ [options]

说明:
    GETA (Genome-wide Electronic Tool for Annotation) 是一个对全基因组序列进行基因预测的流程软件。输入基因组序列、转录组二代测序原始数据和临近物种全基因组蛋白序列，即可一个命令快速得到准确的基因组注释GFF3结果文件。本流程的两大特色：（1）预测准确，程序输出结果中基因数量正常、BUSCO检测的完整性高、基因模型的exon边界准确；（2）运行简单，一个命令全自动化运行得到最终结果。
    软件当前版本：2.7.1。

示例:
    perl __PROG__ --genome genome.fasta --RM_species_Dfam Embryophyta --RM_species_RepBase Embryophyta --pe1 libA.1.fq.gz,libB.1.fq.gz --pe2 libA.2.fq.gz,libB.2.fq.gz --protein homolog.fasta --augustus_species GETA_genus_species --HMM_db /opt/biosoft/bioinfomatics_databases/Pfam/PfamA,/opt/biosoft/bioinfomatics_databases/Pfam/PfamB --config /opt/biosoft/geta/config.txt --out_prefix out --gene_prefix GS01Gene --cpu 120

运行流程:
    Step 0 准备工作：检查参数和依赖的软件及子程序，调用genome_seq_clear.pl整理基因组序列，调用fasta_format_revising.pl修正同源蛋白序列的格式，并读取配置文件（默认为软件安装目录中的config.txt，或--config指定的配置文件）。
    Step 1 重复序列分析：调用runRepeatMaskerModeler，使用RepeatMasker（Dfam数据库、RepBase数据库或--RM_lib）和RepeatModeler对基因组进行重复序列分析，得到重复序列注释和屏蔽了重复序列的基因组序列。仅设置了--no_RepeatModeler而没有设置其它重复序列分析参数时，跳过本步骤。
    Step 2 同源蛋白预测：调用homolog_prediction，利用临近物种的全基因组蛋白序列预测基因模型。
    Step 3 转录组预测：调用NGSReads_prediction，利用RNA-Seq数据预测基因模型，同时得到可信intron信息和碱基测序深度信息。
    Step 4 AUGUSTUS预测：调用runAugustus，以转录组和同源蛋白预测的基因模型为证据训练AUGUSTUS，并结合hints预测基因模型。
    Step 5 整合与过滤：调用curateGeneModels，整合三种方法预测的基因模型，去除位于转座子上的基因模型，使用HMM和/或蛋白数据库验证不可信的基因模型，添加可变剪接并进行ORF分析，去除不完整的转录本，得到最终的基因模型。
    Step 6 输出结果：调用outputResults，输出基因模型的GFF3、GTF和序列文件，以及低质量基因模型、重复序列注释、原始预测结果、BUSCO评估结果和各步骤的统计汇总。
    程序在${out_prefix}.tmp文件夹中运行，各步骤完成后生成.ok文件。程序中断后使用相同的命令重新运行，会跳过已完成的步骤；若想重新运行某个步骤，需要删除该步骤及其后续步骤对应的.ok文件。
    默认使用未屏蔽重复序列的全基因组序列进行同源蛋白预测、转录组预测和AUGUSTUS预测，后续整合过滤时会去除与转座子序列重叠的基因模型；整合、过滤和输出基因模型也使用未屏蔽的基因组序列。若希望程序运行更快，完全不在转座子区域进行基因预测，最终获得更少的基因模型，则添加--gene_predicted_by_masked_genome参数，使用屏蔽了重复序列的基因组序列进行上述三种预测。

输出文件:
    ${out_prefix}.geneModels.gff3、${out_prefix}.geneModels.gtf    最终的基因模型（含可变剪接转录本）；同时输出${out_prefix}.protein.fasta等基因序列文件。
    ${out_prefix}.bestGeneModels.gff3、${out_prefix}.bestGeneModels.gtf    每个基因仅保留一个转录本的基因模型。
    ${out_prefix}.geneModels_lowQuality.gff3    低质量的基因模型，包括未通过HMM或BLASTP验证的、位于转座子上的、5'端或3'端不完整以及没有CDS的基因模型。
    ${out_prefix}.repeat.gff3、${out_prefix}.repeat.stats、${out_prefix}.repeat.lib、${out_prefix}.maskedGenome.fasta    重复序列注释、统计结果、RepeatModeler构建的重复序列库和屏蔽了重复序列的基因组序列。
    ${out_prefix}.homolog_*.gff3、${out_prefix}.NGSReads_*.gff3、${out_prefix}.augustus_prediction.gff3、${out_prefix}.evidence_prediction.gff3、${out_prefix}.All_prediction.gff3    同源蛋白、转录组、AUGUSTUS以及整合过程中的原始基因预测结果。
    ${out_prefix}.gene_prediction.summary    各步骤的统计汇总，用于追踪基因预测结果的可靠性，同时输出到标准输出。
    ${out_prefix}.tmp    临时文件夹，保存各步骤的中间结果，其中6.output_gene_models子文件夹保存BUSCO分析结果。

依赖软件:
    在Rocky 9.2系统使用以下依赖的软件版本对本软件进行了测试并运行成功（按程序调用的先后顺序排列）。
    01. RepeatMasker (version: 4.1.6)
    02. HMMER (nhmmer/hmmscan, version: 3.3.2)
    03. makeblastdb/rmblastn (Version: 2.14.0)
    04. RepeatModeler (version: 2.0.5)
    05. mmseqs (version 15-6f452)
    06. gth (Version: 1.7.3)
    07. exonerate (Version: 2.2.0)
    08. genewise (version: 2.4.1)
    09. miniprot (Version: 0.18)
    10. spaln (Version: 3.0.9)
    11. java (version: jdk-20.0.1)
    12. STAR (Version: 2.7.11b)
    13. hisat2 (version: 2.1.0)
    14. samtools (version: 1.17)
    15. stringtie (Version: 2.1.3)
    16. TIdes (Version: 1.0.0)
    17. diamond (version 2.1.8)
    18. cd-hit (Version: 4.8.1-2019-0228)
    19. Barrnap (Version: 1.10.5)
    20. augustus/etraining (version: 3.5.0)
    21. busco (Version: 5.4.7)
    其中，HMMER中的nhmmer在使用Dfam数据库进行重复序列分析时需要，hmmscan在设置了--HMM_db时需要；RepeatModeler仅在没有设置--RM_lib且没有设置--no_RepeatModeler时需要；mmseqs仅在--homolog_prediction_method使用了gth、exonerate或genewise方法时需要；gth、exonerate、genewise、miniprot和spaln仅在--homolog_prediction_method使用了相应方法时需要，并确保能在终端中直接运行；使用二代测序数据时需要stringtie、samtools，使用FASTQ数据时还需要java（运行软件安装目录中的Trimmomatic-0.38）和比对软件STAR（默认）或hisat2/hisat2-build（--aligner hisat2）；有同源蛋白输入时使用TIdeS（tides.py或tides命令及其依赖的DIAMOND、cd-hit和Barrnap等）预测ORF，否则使用软件安装目录中的TransDecoder-v5.5.0；busco仅在设置了--BUSCO_lineage_dataset时需要。
    AUGUSTUS需要设置环境变量$AUGUSTUS_CONFIG_PATH，并确保augustus、etraining、optimize_augustus.pl、new_species.pl、gff2gbSmallDNA.pl、filterGenes.pl、randomSplit.pl能在终端中直接运行。

使用须知:
    (1) 程序需要调用与本程序位于同一目录下的子程序和脚本：runRepeatMaskerModeler、homolog_prediction、NGSReads_prediction、runAugustus、curateGeneModels、outputResults、GFF3Clear、genome_seq_clear.pl、fasta_format_revising.pl、gff3_to_sequences.pl、gff3ToGtf.pl和check_gene_model_integrity.pl。
    (2) 所有输入文件、数据库以及输出路径中不能包含空格或shell特殊字符（只允许字母、数字和 _-./+:@,= ）。
    (3) 同源蛋白（--protein）和二代测序数据（--pe1/--pe2、--se、--sam）至少需要输入其中一种；--HMM_db、--BLASTP_db和--protein至少需要设置其中一个，用于验证不可信的基因模型。
    (4) 没有输入二代测序数据时，不进行可变剪接分析；没有输入同源蛋白时，不进行同源蛋白预测。

[输入参数]
    --genome <string>    默认值: 无，必须设置
    输入需要进行注释的基因组序列FASTA格式文件。当输入了屏蔽重复序列的基因组文件时，可以不使用--RM_species_Dfam、--RM_species_RepBase和--RM_lib参数，并添加--no_RepeatModeler参数，从而跳过重复序列分析与屏蔽步骤。此时，推荐在FASTA文件中仅对转座子序列使用碱基N进行硬屏蔽，对简单重复或串联重复使用小写字符进行软屏蔽。需要注意的是，此时程序没有重复序列注释信息，因此不会去除位于转座子上的基因模型。

    --pe1 <string> --pe2 <string>    默认值: 无
    输入二代双末端测序的两个FASTQ格式文件。参数支持输入多对数据文件，使用逗号对不同文库的FASTQ文件路径进行分隔即可，--pe1与--pe2中的文件顺序须一一对应。参数也支持输入.gz格式的压缩文件。

    --se <string>    默认值: 无
    输入二代单端测序的FASTQ格式文件。参数支持输入多个数据文件，使用逗号对不同文库的FASTQ文件路径进行分隔即可。参数也支持输入.gz格式的压缩文件。

    --sam <string>    默认值: 无
    输入二代测序SAM格式数据文件。参数支持输入多个数据文件，使用逗号对不同的SAM文件路径进行分隔即可，文件必须以.sam、.bam或.cram结尾。此外，程序支持--pe1/--pe2、--se和--sam这三个参数全部或部分使用，则利用所有输入的数据进行基因组比对，再获取转录本序列进行基因预测。未设置--strand_specific时，通过本参数输入的数据应保留XS标签（如使用HISAT2的--dta或STAR的--outSAMstrandField intronMotif参数进行比对），且其中的参考序列名称需要与基因组FASTA文件一致。

    --strand_specific    默认值: 无
    添加该参数后，则认为所有输入的二代数据是链特异性测序数据，程序则仅在转录本正义链上预测基因。本参数仅支持常规的RF（dUTP）建库方案，即read1比对到转录本的负义链（如Illumina TruSeq Stranded）；所有输入数据必须具有相同的链特异性类型。使用链特异性测序数据并设置本参数后，当两相邻基因有重叠时，能准确预测其边界。

    --protein <string>    默认值: 无
    输入临近物种的全基因组蛋白序列。推荐使用多个（3~10个）物种的全基因组同源蛋白序列。推荐修改蛋白序列的名称，在原名称尾部追加以Species字符开头的物种信息。例如蛋白序列名称为XP_002436309.2，则修改为XP_002436309_2_SpeciesSorghumBicolor。这样有利于在一个基因区域保留更多物种的同源匹配结果，有利于基因预测的准确性。可以考虑使用GETA中附带的fasta_remove_redundancy.pl整合多个物种的全基因组蛋白序列，去冗余时也能附加物种信息。此外，使用的物种数量越多，则能预测的更多准确的基因模型，但是越消耗计算时间。需要注意的是，同源蛋白和二代测序数据，本软件要求至少需要输入其中一种数据用于有证据支持的基因预测。

[重复序列分析参数]
    --RM_species_Dfam <string>    默认值: 无
    输入一个物种或分类名称，以利于调用RepeatMasker程序运用Dfam数据库中相应物种分类的HMM数据信息进行重复序列分析。可以根据文件__SOFTWARE_DIR__/RepeatMasker_lineage.txt的内容查找适用于本参数的物种分类名称。例如，设置Eukaryota适合真核生物、Viridiplantae适合植物、Metazoa适合动物、Fungi适合真菌。当使用了该参数，需要先安装好RepeatMasker软件并配置好Dfam数据库。需要注意的是：Dfam数据库因为较大被分隔成了9份；而RepeatMasker默认仅带有Dfam数据库第0个root部分，适合哺乳动物或真菌等物种；若有需要，则考虑下载适合待分析物种的Dfam数据库并配置到RepeatMasker；例如，对植物物种进行分析则需要使用Dfam数据库第5个Viridiplantae部分，对膜翅目昆虫进行分析则需要使用Dfam数据库第7个Hymenoptera部分。

    --RM_species_RepBase | --RM_species <string>    默认值: 无
    输入一个物种名称，运用调用RepeatMasker程序运用RepBase数据库中相应物种分类的重复核酸序列信息进行重复序列分析。可以根据文件__SOFTWARE_DIR__/RepeatMasker_lineage.txt的内容查找适用于本参数的物种分类名称。当输入了该参数，需要先安装好RepeatMasker软件并配置好RepBase数据库。需要注意的是：RepBase数据库最后一个版本20181026已经较老，且其包含的重复序列数据量较少；RepBase也不再提供免费下载。本参数与--RM_species_Dfam相互独立，两个参数可以设置为不同的值，分别用于RepBase和Dfam数据库的分析。

    --RM_lib <string>    默认值: 无
    输入一个FASTA格式文件，运用其中的重复序列数据信息进行全基因组重复序列分析。该文件往往是RepeatModeler软件对全基因组序列进行分析的结果，表示全基因组上的重复序列。GETA软件默认会调用RepeatModeler对全基因组序列进行分析，得到物种自身的重复序列数据库后，再调用RepeatMasker进行重复序列分析。当添加该参数后，则程序跳过费时的RepeatModeler步骤，能极大减少程序运行时间。此外，程序支持--RM_species_Dfam、--RM_species_RepBase和--RM_lib参数同时使用，则依次使用多种方法进行重复序列分析，最终能合并多个结果并认可任一方法的结果。

    --no_RepeatModeler    默认值: 无
    添加该参数后，程序不会运行RepeatModeler步骤。若同时没有设置--RM_species_Dfam、--RM_species_RepBase和--RM_lib参数，则认为输入的基因组序列已经屏蔽了重复序列，程序跳过整个重复序列分析与屏蔽步骤。

    --gene_predicted_by_masked_genome    默认值: 无
    默认情况下，程序使用输入的全基因组序列（不屏蔽重复序列）进行同源蛋白、转录组和AUGUSTUS三种方法的基因预测，后续整合过滤时会去除与转座子序列重叠的基因模型。添加该参数后，程序使用屏蔽了重复序列的基因组序列进行这三种基因预测，完全不在转座子区域进行基因预测，程序运行更快，但最终获得的基因模型更少。当没有进行重复序列分析（输入的基因组已屏蔽重复序列）时，本参数不生效。

[基因预测参数]
    --homolog_prediction_method <string>    默认值: miniprot
    设置使用同源蛋白进行基因预测的方法，与homolog_prediction程序的--method参数保持一致。其值可以设定为exonerate、genewise、gth、miniprot、spaln或all。若想使用多种方法进行分析，则输入使用逗号分割的多个值，例如exonerate,miniprot,spaln；设置为all则使用全部五种方法。各方法的说明如下：（1）exonerate、genewise和gth：先调用MMseqs2将同源蛋白比对到基因组以确定基因区域，再调用相应软件对基因结构进行精细预测；这三种方法在流程内部已按基因区域整合为一个来源（mmseqs），使用它们需要安装MMseqs2；（2）miniprot：直接将同源蛋白比对到基因组并预测基因模型，作为一个来源；（3）spaln：直接将同源蛋白比对到基因组并预测基因模型，作为另一个来源。程序将各来源最原始的基因模型一起整合去冗余，再按对同源蛋白的覆盖度等指标分为A、B、C、D四类，然后根据gene和intron的Supported_times（被各来源、各同源蛋白支持的次数）过滤低可信的C和D类基因模型，并对支持度不足的基因模型标记Low_support。使用的方法越多，则越消耗计算时间，但结果往往更好。各方法对应的软件（gth、exonerate、genewise、miniprot、spaln、mmseqs）仅在选用了相应方法时才需要安装。

    --aligner <string>    默认值: NGSReads_prediction的默认值（star）
    设置RNA-Seq二代测序数据比对基因组所用的软件，可选star或hisat2，传递给NGSReads_prediction的--aligner参数。使用star需要安装STAR，使用hisat2则需要安装hisat2和hisat2-build。仅输入了--sam数据时，不进行比对，本参数不生效。

    --max_used_read_num <int>    默认值: 自动计算
    设置程序使用的NGSreads数据的最大数量。当程序输入了过量的NGSreads时，会自动根据基因组大小选择一定数据量的read。若添加该参数值，则指定使用的双末端测序的reads对数量或单端测序的read数量。若不设置该参数，程序自动计算使用的read对数量 = ( 2 ** (log10(genome_size / 1,000,000) - 1) )  * 50 M 。即10M的基因组最多使用50M个reads对，PE150测序数据量15G；100M基因组使用100M个reads对，PE150测序数据量30G；1G基因组使用200M个reads对，PE150测序数据量60G。

    --augustus_species <string>    默认值: 自动设置
    输入一个AUGUSTUS物种名称，只能含字母、数字、下划线、点和连字符。则程序在利用转录本或同源蛋白预测的基因模型进行AUGUSTUS Training时，从已有的物种模型开始训练或重头训练新的物种模型。若输入的AUGUSTUS物种模型存在，则在其基础上对其进行优化；若不存在，则生成新的AUGUSTUS物种模型，并进行参数优化。程序进行AUGUSTUS Training需要安装好AUGUSTUS软件，并设置好$AUGUSTUS_CONFIG_PATH环境变量。程序运行完毕后，在临时文件夹的4.augustus/config/species文件夹中有生成名称为本参数值的AUGUSTUS Training物种配置文件夹；若对$AUGUSTUS_CONFIG_PATH路径中指定的species文件夹有写入权限，则将生成的物种配置文件夹拷贝过去。若不输入本参数信息，则程序自动设置本参数的值为“GETA + 基因组FASTA文件名称前缀 + 日期 + 进程ID”，该名称保存在临时文件夹的augustus_species.txt文件中，程序中断后重新运行时会继续使用。

    --genetic_code <int>    默认值: 1
    设置遗传密码，传递给所有调用的子程序，用于预测基因、判断基因模型的完整性、补齐基因模型的首尾以及提取基因序列。该参数对应的值请参考NCBI Genetic Codes: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi 。支持的遗传密码表编号为：1~6、9~16、21~33。

[基因模型整合与过滤参数]
    --HMM_db <string>    默认值: 无
    输入HMM数据库路径（已用hmmpress处理），用于对基因模型进行验证过滤。参数支持输入多个数据库路径，使用逗号进行分隔。当使用多个HMM数据库时，程序过滤在所有数据库中都没有匹配的基因模型。--HMM_db、--BLASTP_db和--protein至少需要设置其中一个。

    --BLASTP_db <string>    默认值: 无
    输入diamond数据库路径，用于对基因模型进行验证过滤。参数支持输入多个数据库路径，使用逗号进行分隔。当使用多个diamond数据库时，程序过滤在所有数据库中都没有匹配的基因模型。若--HMM_db和本参数都没有设置，则使用--protein参数输入的同源蛋白序列作为diamond数据库，进行基因模型过滤。

    --no_alternative_splicing_analysis    默认值: 无
    添加该参数后，程序不会进行可变剪接分析。需要注意的是，当输入了NGS reads数据时，程序默认会根据intron和碱基测序深度信息进行基因的可变剪接分析。

    --config <string>    默认值: 软件安装目录中的config.txt
    输入一个参数配置文件路径，用于设置本程序调用的其它命令的详细参数。若不设置该参数，则使用软件安装目录中的config.txt文件；该文件列出了各子程序的参数及其默认值，并附有详细的注释说明。若需要修改参数，建议将config.txt复制一份并修改，再通过本参数输入；此时程序只读取本参数指定的文件，文件中没有出现的小节或参数，由各子程序使用其自身的默认值。配置文件以[程序名]作为小节标题，小节下面写附加传递给该程序的命令行参数，支持的小节有：[runRepeatMaskerModeler]、[homolog_prediction]、[NGSReads_prediction]、[runAugustus]、[curateGeneModels]和[outputResults]；以#开头的内容为注释。本程序命令行中设置的参数，优先于配置文件中的同名参数。程序不再根据基因组大小自动切换配置文件。当软件预测的基因数量异常时往往要修改基因模型的过滤阈值，此时修改配置文件并通过本参数输入，再次运行GETA流程（需要删除相应步骤的.ok文件）。同源蛋白预测程序homolog_prediction的更多参数（如--max_intron_size、--miniprot_options、--spaln_species、--spaln_options等）可以写在配置文件的[homolog_prediction]部分；curateGeneModels的--validate_options、--alt_splicing_options等参数可以写在[curateGeneModels]部分。

[输出参数]
    --out_prefix <string>    默认值: out
    设置输出文件或临时文件夹前缀。临时文件夹为${out_prefix}.tmp，各步骤的中间结果保存在其中。

    --gene_prefix <string>    默认值: gene
    设置输出GFF3文件中的基因名称前缀，基因ID由该前缀加6位数字编号组成。低质量基因的ID前缀为该前缀去掉"gene"后加"LqGene"。

    --BUSCO_lineage_dataset <string>    默认值: 无
    输入BUSCO数据库路径，则程序额外对最终基因模型以及同源蛋白、转录组、AUGUSTUS预测和整合过程中的各个原始预测结果的蛋白序列进行BUSCO分析，用于比较各种基因预测方法的结果。本参数支持输入多个BUSCO数据库路径，使用逗号进行分隔，则分别利用多个数据库进行分析。可以根据__SOFTWARE_DIR__/BUSCO_lineages_list.2021-12-14.txt文件内容选择合适的BUSCO数据库，BUSCO以离线模式运行，因此需要预先下载好数据库。BUSCO的结果输出到6.output_gene_models子目录下和gene_prediction.summary文件中。

    --delete_unimportant_intermediate_files    默认值: 无
    添加该参数后，若程序运行成功，会删除不重要的中间文件（如比对BAM文件、RepeatMasker/RepeatModeler、MMseqs2、miniprot、spaln、AUGUSTUS预测、HMM/BLASTP验证和可变剪接分析的临时文件夹），仅保留最少的、较小的、重要的中间结果文件。

[运行参数]
    --cpu <int>    默认值: 4
    设置程序运行使用的CPU线程数，传递给所有调用的子程序。

[帮助参数]
    --help    默认值: 无
    Display the English usage and exit.

    --chinese_help    默认值: 无
    显示中文用法说明并退出。
USAGE

########################################################################
# 全局变量
########################################################################
my $prog_name    = basename($0);
my $self_path    = abs_path($0) or die "Error: can not resolve the real path of $0.\n";
my $bin_path     = dirname($self_path);
my $software_dir = $bin_path; $software_dir =~ s/\/bin$//;
my $pwd          = getcwd();
my $command_line = "$0 " . join(" ", map { shq($_) } @ARGV);

# 英文用法说明位于文件最末尾的__DATA__部分
my $usage_en = do { local $/; <DATA> };
$usage_en = '' unless defined $usage_en;
foreach my $u ($usage_cn, $usage_en) {
    $u =~ s/__PROG__/$0/g;
    $u =~ s/__SOFTWARE_DIR__/$software_dir/g;
}

# 程序直接运行（无任何参数）时，给出英文用法说明
if (@ARGV == 0) { print STDERR $usage_en; exit 1; }

# homolog_prediction支持的同源蛋白基因预测方法；exonerate、genewise和gth属于同一个来源(mmseqs)
my @ALL_HOMOLOG_METHODS = qw/gth exonerate genewise miniprot spaln/;
my %MMSEQS_METHOD       = map { $_ => 1 } qw/gth exonerate genewise/;
# 各子程序支持的NCBI遗传密码表编号
my %SUPPORTED_GENETIC_CODE = map { $_ => 1 } (1 .. 6, 9 .. 16, 21 .. 33);

########################################################################
# 解析并检查参数
########################################################################
my %o;
Getopt::Long::Configure("no_auto_abbrev");
GetOptions(\%o,
    'genome=s', 'pe1=s', 'pe2=s', 'se=s', 'sam=s', 'strand_specific', 'protein=s',
    'RM_species_Dfam=s', 'RM_species_RepBase=s', 'RM_lib=s', 'no_RepeatModeler', 'gene_predicted_by_masked_genome',
    'homolog_prediction_method=s', 'aligner=s', 'max_used_read_num=i', 'augustus_species=s', 'genetic_code=i',
    'HMM_db=s', 'BLASTP_db=s', 'no_alternative_splicing_analysis', 'config=s',
    'out_prefix=s', 'gene_prefix=s', 'BUSCO_lineage_dataset=s', 'delete_unimportant_intermediate_files',
    'cpu=i',
    'help|h', 'chinese_help',
) or do { print STDERR "\nRun '$prog_name --help' for usage.\n"; exit 1; };

if ($o{chinese_help}) { print $usage_cn; exit 0; }
if ($o{help})         { print $usage_en; exit 0; }
die "Error: unrecognized argument(s): @ARGV\n\nRun '$prog_name --help' for usage.\n" if @ARGV;

$o{cpu}                       //= 4;
$o{out_prefix}                //= "out";
$o{gene_prefix}               //= "gene";
$o{genetic_code}              //= 1;
$o{homolog_prediction_method} //= "miniprot";
my $cpu          = $o{cpu};
my $genetic_code = $o{genetic_code};
my $gene_prefix  = $o{gene_prefix};
die "Error: --cpu must be a positive integer.\n" if $cpu < 1;
die "Error: --genetic_code $genetic_code is not supported, the supported genetic code tables are: 1-6, 9-16, 21-33.\n" unless $SUPPORTED_GENETIC_CODE{$genetic_code};
die "Error: --gene_prefix must not be empty or contain whitespace.\n" if $gene_prefix eq '' || $gene_prefix =~ /\s/;
die "Error: --max_used_read_num must be a positive integer.\n" if defined $o{max_used_read_num} && $o{max_used_read_num} < 1;

# ---- 基因组、二代测序数据和同源蛋白 ----
die "Error: --genome is required.\n\nRun '$prog_name --help' for usage.\n" unless defined $o{genome};
die "Error: genome file $o{genome} does not exist or is empty.\n" unless -s $o{genome};
my $genome_input = abs_path($o{genome});

die "Error: --pe1 and --pe2 must be set together.\n" if ( defined $o{pe1} xor defined $o{pe2} );
my (@pe1, @pe2, @se, @sam);
if (defined $o{pe1}) {
    @pe1 = abs_file_list($o{pe1}, '--pe1');
    @pe2 = abs_file_list($o{pe2}, '--pe2');
    die "Error: the number of files in --pe1 (" . scalar(@pe1) . ") and --pe2 (" . scalar(@pe2) . ") must be equal.\n" unless @pe1 == @pe2;
}
@se  = abs_file_list($o{se},  '--se')  if defined $o{se};
@sam = abs_file_list($o{sam}, '--sam') if defined $o{sam};
foreach my $f (@sam) {
    die "Error: --sam input $f should be a .sam, .bam or .cram file.\n" unless $f =~ m/\.(?:sam|bam|cram)$/i;
}
my $has_reads      = ( @pe1 || @se || @sam ) ? 1 : 0;
my $need_alignment = ( @pe1 || @se ) ? 1 : 0;

my $protein_input;
if (defined $o{protein}) {
    die "Error: homologous protein file $o{protein} does not exist or is empty.\n" unless -s $o{protein};
    $protein_input = abs_path($o{protein});
}
die "Error: no RNA-Seq reads or homologous proteins were input, please use --pe1/--pe2, --se, --sam or --protein.\n" unless ($has_reads || defined $protein_input);

if (defined $o{aligner}) {
    $o{aligner} = lc $o{aligner};
    die "Error: --aligner must be either 'star' or 'hisat2'.\n" unless $o{aligner} eq 'star' || $o{aligner} eq 'hisat2';
}

# ---- 重复序列分析 ----
# --RM_species_Dfam和--RM_species_RepBase相互独立，可以设置为不同的值
my $RM_dfam    = $o{RM_species_Dfam};
my $RM_repbase = $o{RM_species_RepBase};
$RM_dfam    = undef if defined $RM_dfam    && $RM_dfam    eq '';
$RM_repbase = undef if defined $RM_repbase && $RM_repbase eq '';
my $RM_lib;
if (defined $o{RM_lib}) {
    die "Error: repeat library file $o{RM_lib} does not exist or is empty.\n" unless -s $o{RM_lib};
    $RM_lib = abs_path($o{RM_lib});
}
my $no_RepeatModeler = $o{no_RepeatModeler} ? 1 : 0;
my $use_modeler      = ( !defined $RM_lib && !$no_RepeatModeler ) ? 1 : 0;
# 输入的基因组已屏蔽重复序列：不设置任何重复序列分析时，跳过整个重复序列分析步骤
my $skip_repeat      = ( $no_RepeatModeler && !defined $RM_dfam && !defined $RM_repbase && !defined $RM_lib ) ? 1 : 0;

# ---- 同源蛋白预测方法 ----
my %sel_method;
foreach my $m ( split /[\s,]+/, lc($o{homolog_prediction_method}) ) {
    next if $m eq "";
    if ( $m eq "all" ) { $sel_method{$_} = 1 foreach @ALL_HOMOLOG_METHODS; }
    elsif ( grep { $_ eq $m } @ALL_HOMOLOG_METHODS ) { $sel_method{$m} = 1; }
    else { die "Error: the supported value of --homolog_prediction_method should be exonerate, genewise, gth, miniprot, spaln or all, you input '$m' was wrong.\n"; }
}
my @homolog_methods    = grep { $sel_method{$_} } @ALL_HOMOLOG_METHODS;
die "Error: No method was specified by --homolog_prediction_method.\n" unless @homolog_methods;
my $homolog_method_str = join ",", @homolog_methods;

# ---- 基因模型验证数据库（转为绝对路径） ----
my ($HMM_db, $BLASTP_db);
if (defined $o{HMM_db}) {
    my @list;
    foreach my $p ( grep { length $_ } split /,/, $o{HMM_db} ) {
        my $a = File::Spec->rel2abs($p);
        die "Error: HMM database $p was not found.\n" unless ( -e $a || -e "$a.hmm" || -e "$a.h3m" );
        push @list, $a;
    }
    $HMM_db = join ",", @list if @list;
}
if (defined $o{BLASTP_db}) {
    my @list;
    foreach my $p ( grep { length $_ } split /,/, $o{BLASTP_db} ) {
        my $a = File::Spec->rel2abs($p);
        die "Error: diamond database $p was not found.\n" unless ( -e $a || -e "$a.dmnd" );
        push @list, $a;
    }
    $BLASTP_db = join ",", @list if @list;
}
die "Error: at least one of --HMM_db, --BLASTP_db and --protein is required to validate the gene models.\n" unless ( $HMM_db || $BLASTP_db || defined $protein_input );

# ---- 配置文件、BUSCO数据库、输出前缀 ----
# 默认使用软件安装目录中的config.txt；--config指定的文件则取代它（只读取该文件）
my $config_input;
if (defined $o{config}) {
    die "Error: configuration file $o{config} does not exist or is empty.\n" unless -s $o{config};
    $config_input = abs_path($o{config});
}
else {
    $config_input = "$software_dir/config.txt";
    die "Error: the default configuration file $config_input does not exist or is empty, please check the installation of GETA or use --config to input a configuration file.\n" unless -s $config_input;
}
my @BUSCO_db;
if (defined $o{BUSCO_lineage_dataset}) {
    foreach my $p ( grep { length $_ } split /,/, $o{BUSCO_lineage_dataset} ) {
        my $a = File::Spec->rel2abs($p);
        die "Error: BUSCO lineage dataset directory $p was not found.\n" unless -d $a;
        push @BUSCO_db, $a;
    }
}
my $out_prefix = File::Spec->rel2abs($o{out_prefix});
make_path( dirname($out_prefix) ) unless -d dirname($out_prefix);
my $tmp_dir = "$out_prefix.tmp";

# 各子程序拼接shell命令时不对路径加引号，因此路径中不能含空格或shell特殊字符
check_path_chars( $genome_input, @pe1, @pe2, @se, @sam, $bin_path, $out_prefix, $tmp_dir, $config_input );
check_path_chars( $protein_input )  if defined $protein_input;
check_path_chars( $RM_lib )         if defined $RM_lib;
check_path_chars( split /,/, $HMM_db )     if $HMM_db;
check_path_chars( split /,/, $BLASTP_db )  if $BLASTP_db;
check_path_chars( @BUSCO_db );

# ---- 检测依赖的软件和子程序 ----
&detecting_dependent_softwares();

########################################################################
# Step 0: 程序运行前的准备工作
########################################################################
make_path($tmp_dir) unless -d $tmp_dir;
die "Error: can not create the temporary directory $tmp_dir.\n" unless -d $tmp_dir;
$tmp_dir = abs_path($tmp_dir);
chdir $tmp_dir or die "Error: can not enter directory $tmp_dir, $!\n";
print STDERR "\n#################### $prog_name | Start ####################\n";
print STDERR "PWD: $pwd\n" . localtime() . ": CMD: $command_line\n";
print STDERR "\nPWD: $tmp_dir\n";

# AUGUSTUS物种名称：未设置时自动生成，并保存到临时文件夹，保证程序中断后重新运行时使用同一个名称
my $augustus_species;
{
    my $species_file = "$tmp_dir/augustus_species.txt";
    if (defined $o{augustus_species}) {
        $augustus_species = $o{augustus_species};
        die "Error: --augustus_species can only contain letters, digits, underscore, dot and hyphen.\n" unless $augustus_species =~ /^[A-Za-z0-9_.\-]+$/;
        write_file($species_file, "$augustus_species\n");
    }
    elsif ( -s $species_file ) {
        open my $in, "<", $species_file or die "Error: Can not open file $species_file, $!\n";
        $augustus_species = <$in>; close $in;
        chomp $augustus_species;
    }
    unless ( defined $augustus_species && length $augustus_species ) {
        my $name = basename($genome_input);
        $name =~ s/\.fa.*//;
        $name =~ s/[^A-Za-z0-9_.\-]/_/g;
        my $date = strftime("%Y%m%d%H%M%S", localtime);
        $augustus_species = "GETA_${name}${date}_$$";
        write_file($species_file, "$augustus_species\n");
    }
}

# 0.1 准备基因组序列：
# 读取FASTA序列以>开始的头部时，去除第一个空及之后的字符，按从长到短排序。若序列名有重复，则修正序列名并保留其序列。
my $genome = "$tmp_dir/genome.fasta";
&execute_cmds("$bin_path/genome_seq_clear.pl --no_rename --no_change_bp $genome_input > $genome 2> $tmp_dir/genome_seq_clear.log", "$tmp_dir/genome.fasta.ok");
die "Error: $genome is empty, please check $tmp_dir/genome_seq_clear.log .\n" unless -s $genome;
my $genome_size = get_genome_size($genome, "$tmp_dir/genome.size");

# 0.2 读取配置文件（默认为软件安装目录中的config.txt，或--config指定的文件），得到各子程序的附加参数
my %config = choose_config_file();

# 0.3 准备直系同源基因蛋白序列：
# 读取FASTA序列以>开始的头部时，去除第一个空及之后的字符，将所有怪异字符变为下划线字符；若遇到 > 符号不在句首，则表示fasta文件格式有误，删除该 > 及到下一个 > 之前的数据；去除序列中尾部的换行符, 将所有小写字符氨基酸变换为大写字符。
my $protein;
if (defined $protein_input) {
    $protein = "$tmp_dir/homolog.fasta";
    &execute_cmds("$bin_path/fasta_format_revising.pl --seq_type protein --line_length 80 $protein_input > $protein 2> $tmp_dir/homolog.fasta.fasta_format_revising.log", "$tmp_dir/homolog.fasta.ok");
    die "Error: $protein is empty, please check $tmp_dir/homolog.fasta.fasta_format_revising.log .\n" unless -s $protein;
}

# 0.4 设置最大使用的二代测序数据量 = ( 2 ** (log10(genome_size / 1,000,000) - 1) ) * 50 M reads pair
# 即10M的基因组最多使用50M个reads对，PE150测序数据量15G；100M基因组使用100M个reads对，PE150测序数据量30G；1G基因组使用200M个reads对，PE150测序数据量60G。
my $max_used_read_num = 50_000_000;
if (defined $o{max_used_read_num}) {
    $max_used_read_num = $o{max_used_read_num};
}
else {
    my $x = log($genome_size / 1_000_000) / log(10) - 1;
    $x = 0 if $x < 0;
    my $n = int( (2 ** $x) * 50_000_000 );
    $max_used_read_num = $n if $n > $max_used_read_num;
}

########################################################################
# Step 1: RepeatMasker and RepeatModeler
########################################################################
print STDERR "\n============================================\n";
print STDERR "Step 1: RepeatMasker and RepeatModeler (" . localtime() . ")\n";
my $repeat_dir = "$tmp_dir/1.RepeatMasker";
if ($skip_repeat) {
    print STDERR "Note: --no_RepeatModeler is set without --RM_species_Dfam, --RM_species_RepBase or --RM_lib, so the input genome is regarded as repeat-masked and this step is skipped.\n";
}
else {
    make_path($repeat_dir) unless -d $repeat_dir;
    my @p;
    push @p, "--RM_species_Dfam $RM_dfam"       if defined $RM_dfam;
    push @p, "--RM_species_RepBase $RM_repbase" if defined $RM_repbase;
    push @p, "--RM_lib $RM_lib" if defined $RM_lib;
    push @p, "--no_RepeatModeler" if $no_RepeatModeler;
    push @p, "--cpu $cpu", "--out_repeat_percentage", "--tmp_dir $repeat_dir";
    my $cmd = "$bin_path/runRepeatMaskerModeler " . extra_opts("runRepeatMaskerModeler") . join(" ", @p) . " $genome";
    &execute_cmds($cmd, "$tmp_dir/1.RepeatMasker.ok");
}

# 基因预测使用的基因组序列：默认使用未屏蔽重复序列的全基因组序列；
# 设置--gene_predicted_by_masked_genome后，使用屏蔽了重复序列的基因组序列
my $predict_genome = $genome;
if ($o{gene_predicted_by_masked_genome}) {
    my $masked_genome = "$repeat_dir/genome.masked.fasta";
    if ( -s $masked_genome ) {
        $predict_genome = $masked_genome;
    }
    else {
        print STDERR "Note: --gene_predicted_by_masked_genome is set but no repeat-masked genome was generated, so the input genome is used for gene prediction.\n";
    }
}
print STDERR "The genome used for gene prediction by homolog, NGS reads and AUGUSTUS: $predict_genome\n";

########################################################################
# Step 2: homolog_prediction
########################################################################
print STDERR "\n============================================\n";
print STDERR "Step 2: Homolog_prediction (" . localtime() . ")\n";
my $homolog_dir = "$tmp_dir/2.homolog_prediction";
my $homolog_raw = "$homolog_dir/homolog_prediction.raw.gff3";
make_path($homolog_dir) unless -d $homolog_dir;
{
    my $cmd;
    if (defined $protein) {
        $cmd = "$bin_path/homolog_prediction " . extra_opts("homolog_prediction") . "--ignore_signature --tmp_dir $homolog_dir --cpu $cpu --genetic_code $genetic_code --output_alignment_GFF3 $homolog_dir/homolog_alignment.gff3 --output_raw_GFF3 $homolog_raw $protein $predict_genome > $homolog_dir/homolog_prediction.gff3 2> $homolog_dir/homolog_prediction.log";
    }
    else {
        print STDERR "Note: no homologous protein was input, this step is skipped.\n";
        $cmd = "touch $homolog_raw";
    }
    &execute_cmds($cmd, "$tmp_dir/2.homolog_prediction.ok");
}

########################################################################
# Step 3: NGSReads_prediction
########################################################################
print STDERR "\n============================================\n";
print STDERR "Step 3: NGSReads_prediction (" . localtime() . ")\n";
my $ngs_dir = "$tmp_dir/3.NGSReads_prediction";
my $ngs_raw = "$ngs_dir/NGSReads_prediction.raw.gff3";
make_path($ngs_dir) unless -d $ngs_dir;
{
    my $cmd;
    if ($has_reads) {
        my @p;
        push @p, "--pe1 " . join(",", @pe1) . " --pe2 " . join(",", @pe2) if @pe1;
        push @p, "--se "  . join(",", @se)  if @se;
        push @p, "--sam " . join(",", @sam) if @sam;
        push @p, "--strand_specific" if $o{strand_specific};
        push @p, "--aligner $o{aligner}" if defined $o{aligner};
        push @p, "--genetic_code $genetic_code";
        push @p, "--homolog_protein $protein" if defined $protein;
        push @p, "--pe_used_pair_num $max_used_read_num --se_used_read_num $max_used_read_num";
        $cmd = "$bin_path/NGSReads_prediction " . extra_opts("NGSReads_prediction") . join(" ", @p) . " --cpu $cpu --tmp_dir $ngs_dir --output_alignment_GFF3 $ngs_dir/NGSReads_alignment.gff3 --output_raw_GFF3 $ngs_raw --intron_info_out $ngs_dir/intron.txt --base_depth_out $ngs_dir/base_depth.txt $predict_genome > $ngs_dir/NGSReads_prediction.gff3 2> $ngs_dir/NGSReads_prediction.log";
    }
    else {
        print STDERR "Note: no RNA-Seq data was input, this step is skipped.\n";
        $cmd = "touch $ngs_raw";
    }
    &execute_cmds($cmd, "$tmp_dir/3.NGSReads_prediction.ok");
}

########################################################################
# Step 4: AUGUSTUS Training and gene prediction
########################################################################
print STDERR "\n============================================\n";
print STDERR "Step 4: Augustus/HMM Training (" . localtime() . ")\n";
my $augustus_dir = "$tmp_dir/4.augustus";
make_path($augustus_dir) unless -d $augustus_dir;
my $augustus_gff3 = "$augustus_dir/augustus_prediction.gff3";
{
    # 作为证据的基因模型：NGSReads预测的优先，其次是同源蛋白预测的
    my @evidence;
    push @evidence, $ngs_raw    if $has_reads && -s $ngs_raw;
    push @evidence, $homolog_raw if defined $protein && -s $homolog_raw;
    die "Error: no gene model predicted by NGS reads or homologous proteins is available as the evidence of AUGUSTUS Training.\n" unless @evidence;
    my $intron_tab = ( $has_reads && -s "$ngs_dir/intron.txt" ) ? " --intron_tab $ngs_dir/intron.txt" : "";
    my $cmd = "$bin_path/runAugustus " . extra_opts("runAugustus") . "--cpu $cpu --genetic_code $genetic_code --input_evidence_gff3 " . join(",", @evidence) . " --augustus_species $augustus_species --tmp_dir $augustus_dir$intron_tab $predict_genome > $augustus_gff3 2> $augustus_dir/augustus_prediction.log";
    &execute_cmds($cmd, "$tmp_dir/4.augustus.ok");
}

########################################################################
# Step 5: curateGeneModels
########################################################################
print STDERR "\n============================================\n";
print STDERR "Step 5: CurateGeneModels (" . localtime() . ")\n";
my $curate_dir = "$tmp_dir/5.curate_gene_models";
make_path($curate_dir) unless -d $curate_dir;
{
    my @p;
    push @p, "--HMM_db $HMM_db" if $HMM_db;
    if    ($BLASTP_db)  { push @p, "--diamond_db $BLASTP_db"; }
    elsif (!$HMM_db)    { push @p, "--diamond_db $protein"; }
    push @p, "--genetic_code $genetic_code";
    my $repeat_gff3 = "$repeat_dir/genome.repeat.gff3";
    push @p, "--repeat $repeat_gff3" if -s $repeat_gff3;
    if ( $has_reads && !$o{no_alternative_splicing_analysis} && -s "$ngs_dir/intron.txt" && -s "$ngs_dir/base_depth.txt" ) {
        push @p, "--intron $ngs_dir/intron.txt --base_depth $ngs_dir/base_depth.txt";
    }
    my $cmd = "$bin_path/curateGeneModels " . extra_opts("curateGeneModels") . join(" ", @p) . " --cpu $cpu --tmp_dir $curate_dir $genome $ngs_raw $homolog_raw $augustus_gff3 > $curate_dir/curate_gene_models.gff3 2> $curate_dir/curate_gene_models.log";
    &execute_cmds($cmd, "$tmp_dir/5.curate_gene_models.ok");
}

########################################################################
# Step 6: Output gene models
########################################################################
print STDERR "\n============================================\n";
print STDERR "Step 6: Output gene models (" . localtime() . ")\n";
make_path("$tmp_dir/6.output_gene_models") unless -d "$tmp_dir/6.output_gene_models";
{
    my $busco = @BUSCO_db ? " --BUSCO_lineage_dataset " . join(",", @BUSCO_db) . " --busco_for_all_gff3" : "";
    my $cmd = "$bin_path/outputResults " . extra_opts("outputResults") . "--cpu $cpu --genome $genome --genetic_code $genetic_code --out_prefix $out_prefix --gene_prefix $gene_prefix$busco $tmp_dir 2> $tmp_dir/6.output_gene_models/outputResults.log";
    &execute_cmds($cmd, "$tmp_dir/6.output_gene_models.ok");
}

# 删除不重要的中间文件
&delete_unimportant_intermediate_files() if $o{delete_unimportant_intermediate_files};

print "\n============================================\n";
my $summary_file = "$out_prefix.gene_prediction.summary";
if (-s $summary_file) {
    open my $in, "<", $summary_file or die "Error: Can not open file $summary_file, $!\n";
    print while <$in>;
    close $in;
}

print STDERR "\n============================================\n";
print STDERR "GETA complete successfully! (" . localtime() . ")\n\n";
exit 0;

########################################################################
# 子程序
########################################################################

# 对shell命令中的参数加单引号（仅用于记录命令行）
sub shq {
    my ($s) = @_;
    return "''" if $s eq "";
    return $s if $s =~ m/^[\w\/.,:=+\-]+$/;
    $s =~ s/'/'\\''/g;
    return "'$s'";
}

sub write_file {
    my ($file, $text) = @_;
    open my $out, ">", $file or die "Error: Can not create file $file, $!\n";
    print $out $text;
    close $out;
    return 1;
}

# 检测命令是否可以在终端中直接运行（存在于PATH环境变量的目录中，且具有可执行权限）。
sub find_in_path {
    my ($cmd) = @_;
    foreach my $dir ( split /:/, ($ENV{PATH} // "") ) {
        next if $dir eq "";
        my $path = "$dir/$cmd";
        return $path if ( -f $path && -x $path );
    }
    return;
}

# 检查路径中是否含有空格或shell特殊字符
sub check_path_chars {
    foreach my $p (@_) {
        next unless defined $p;
        die "Error: the path '$p' contains spaces or special characters, please use a path only with letters, digits and _-./+:@,=\n" if $p =~ /[^\w\-\.\/\+:@,=]/;
    }
    return 1;
}

# 将逗号分隔的文件列表转换为绝对路径列表，并检查文件是否存在
sub abs_file_list {
    my ($list, $what) = @_;
    my @out;
    foreach my $f ( split /,/, $list ) {
        $f =~ s/^\s+|\s+$//g;
        next if $f eq '';
        die "Error: $what file $f does not exist or is empty.\n" unless -s $f;
        push @out, abs_path($f);
    }
    die "Error: no file was found in the $what parameter.\n" unless @out;
    return @out;
}

# 获取基因组大小（所有序列的总碱基数），结果缓存在size_file中。
sub get_genome_size {
    my ($fasta, $size_file) = @_;
    if ( -s $size_file && (stat($size_file))[9] >= (stat($fasta))[9] ) {
        open my $in, "<", $size_file or die "Error: Can not open file $size_file, $!\n";
        my $line = <$in>;
        close $in;
        return $1 if ( defined $line && $line =~ m/^(\d+)\s*$/ && $1 > 0 );
    }
    my $size = 0;
    open my $in, "<", $fasta or die "Error: Can not open file $fasta, $!\n";
    while (<$in>) {
        next if m/^>/;
        s/\s+//g;
        $size += length($_);
    }
    close $in;
    die "Error: No genome sequence was found in file $fasta\n" unless $size > 0;
    write_file($size_file, "$size\n");
    return $size;
}

# 子程序，读取配置文件。以[程序名]为小节标题，#开头（或行内#之后）为注释；
# 文件中出现的同名小节的内容会累加到一起。
sub read_config_file {
    my ($file, $config) = @_;
    open my $in, "<", $file or die "Error: Can not open configuration file $file, $!\n";
    my $tag;
    while (<$in>) {
        s/#.*//; next if m/^\s*$/; s/^\s+//; s/\s*$/ /;
        if ( m/^\[(.*?)\]/ ) {
            $tag = $1;
            $tag =~ s/^\s+|\s+$//g;
            next;
        }
        next unless defined $tag;
        $config->{$tag} .= $_;
    }
    close $in;
    return 1;
}

# 读取配置文件：默认为软件安装目录中的config.txt，--config指定了文件时只读取该文件。
# 配置文件中没有出现的小节或参数，由各子程序使用其自身的默认值。
sub choose_config_file {
    my %config;

    &read_config_file($config_input, \%config);
    print STDERR "Configuration file: $config_input\n";

    # --homolog_prediction_method参数的优先级高于配置文件中homolog_prediction的--method参数。
    # 若配置文件中没有--method参数，则追加之，避免homolog_prediction使用其自身的默认方法。
    my $method_option = "--method $homolog_method_str";
    if ( defined $config{homolog_prediction} && $config{homolog_prediction} =~ m/--method\s+\S+/ ) {
        $config{homolog_prediction} =~ s/--method\s+\S+/$method_option/;
    }
    else {
        $config{homolog_prediction} = ( $config{homolog_prediction} // "" ) . " $method_option ";
    }

    # 生成本次程序运行使用的配置文件，用于记录
    open my $out, ">", "$tmp_dir/config.txt" or die "Error: Can not create file $tmp_dir/config.txt, $!\n";
    foreach ( sort keys %config ) {
        print $out "[$_]\n$config{$_}\n\n";
    }
    close $out;

    return %config;
}

# 返回配置文件中某个程序的附加参数（以空格结尾，位于命令行的最前面，因此命令行中后面的同名参数优先）。
sub extra_opts {
    my ($name) = @_;
    my $x = $config{$name};
    return "" unless defined $x && $x =~ /\S/;
    $x =~ s/^\s+//; $x =~ s/\s+$//;
    return "$x ";
}

# 子程序，用于执行调用的Linux命令。程序运行完毕后，生成.ok文件。
sub execute_cmds {
    my $ok_file = pop @_;

    if ( -e $ok_file ) {
        foreach ( @_ ) {
            print STDERR "CMD(Skipped): $_\n";
        }
    }
    else {
        foreach ( @_ ) {
            print STDERR scalar(localtime) . ": CMD: $_\n";
            my $ret = system($_);
            if ( $ret != 0 ) {
                my $reason = $ret == -1 ? "failed to start: $!" : ( $ret & 127 ) ? "killed by signal " . ( $ret & 127 ) : "exit status " . ( $ret >> 8 );
                die "Error: failed to execute ($reason): $_\n";
            }
        }
        open my $ok, ">", $ok_file or die "Error: Can not create file $ok_file, $!\n"; close $ok;
    }

    return 1;
}

# 子程序，检测依赖的软件和子程序（软件按程序调用的先后顺序检测）。在正式运行前一次性列出所有缺失项。
sub detecting_dependent_softwares {
    print STDERR "\n============================================\n";
    print STDERR "Detecting the dependent softwares:\n";
    my @missing;
    my $need_cmd = sub {
        my ($cmd, $why) = @_;
        if ( find_in_path($cmd) ) { print STDERR "$cmd:\tOK\n"; }
        else { push @missing, "  - '$cmd' ($why) can not be executed directly in the terminal, please make sure it is installed and added to the PATH environment variable."; }
    };
    my $need_file = sub {
        my ($file, $why) = @_;
        if ( -e $file ) { print STDERR "$file:\tOK\n"; }
        else { push @missing, "  - '$file' ($why) was not found."; }
    };

    # 子程序和脚本
    my @scripts = qw/genome_seq_clear.pl runAugustus curateGeneModels outputResults GFF3Clear gff3_to_sequences.pl gff3ToGtf.pl check_gene_model_integrity.pl/;
    push @scripts, "runRepeatMaskerModeler" unless $skip_repeat;
    push @scripts, "fasta_format_revising.pl", "homolog_prediction" if defined $protein_input;
    push @scripts, "NGSReads_prediction" if $has_reads;
    foreach my $s (@scripts) {
        my $p = "$bin_path/$s";
        if ( -f $p && -x $p ) { print STDERR "$s:\tOK\n"; }
        else { push @missing, "  - '$s' was not found or is not executable in the directory of this program: $bin_path ."; }
    }

    # Step 1: RepeatMasker / HMMER(nhmmer) / rmblastn / RepeatModeler
    unless ($skip_repeat) {
        $need_cmd->("RepeatMasker", "repeat analysis");
        $need_cmd->("nhmmer", "HMMER, RepeatMasker with the Dfam database") if defined $RM_dfam;
        $need_cmd->("rmblastn", "RepeatMasker with RepBase, repeat library or RepeatModeler") if ( defined $RM_repbase || defined $RM_lib || $use_modeler );
        if ($use_modeler) {
            $need_cmd->("BuildDatabase", "RepeatModeler");
            $need_cmd->("RepeatModeler", "RepeatModeler");
        }
    }

    # Step 2: 同源蛋白预测（各软件仅在--homolog_prediction_method使用了相应方法时需要）
    if ( defined $protein_input ) {
        my @mmseqs_methods = grep { $MMSEQS_METHOD{$_} } @homolog_methods;
        $need_cmd->("mmseqs", "--homolog_prediction_method " . join(",", @mmseqs_methods)) if @mmseqs_methods;
        $need_cmd->($_, "--homolog_prediction_method $_") foreach @homolog_methods;
    }

    # Step 3: 转录组预测（java/Trimmomatic、比对软件、samtools、stringtie、ORF预测软件）
    if ($has_reads) {
        if ($need_alignment) {
            $need_cmd->("java", "Trimmomatic");
            $need_file->("$software_dir/Trimmomatic-0.38/trimmomatic-0.38.jar", "Trimmomatic");
            my $aligner = defined $o{aligner} ? $o{aligner} : "star";
            if ( $aligner eq 'hisat2' ) {
                $need_cmd->("hisat2", "--aligner hisat2");
                $need_cmd->("hisat2-build", "--aligner hisat2");
            }
            else {
                $need_cmd->("STAR", "--aligner star");
            }
        }
        $need_cmd->("samtools", "RNA-Seq alignment processing");
        $need_cmd->("stringtie", "transcript assembly");
        # ORF预测软件：有同源蛋白输入时使用TIdeS，否则使用TransDecoder
        if ( defined $protein_input ) {
            if ( find_in_path("tides.py") || find_in_path("tides") ) { print STDERR "tides:\tOK\n"; }
            else { push @missing, "  - 'tides.py' or 'tides' (TIdeS, ORF prediction) can not be executed directly in the terminal."; }
        }
        else {
            $need_file->("$software_dir/TransDecoder-v5.5.0/TransDecoder.LongOrfs", "TransDecoder");
            $need_file->("$software_dir/TransDecoder-v5.5.0/TransDecoder.Predict", "TransDecoder");
        }
    }

    # Step 4: AUGUSTUS（diamond同时用于Step 4中证据基因模型的挑选和Step 5中的BLASTP验证）
    unless ( defined $ENV{AUGUSTUS_CONFIG_PATH} && -d "$ENV{AUGUSTUS_CONFIG_PATH}/species" ) {
        push @missing, "  - the environment variable \$AUGUSTUS_CONFIG_PATH is not set correctly (no species directory was found).";
    }
    $need_cmd->($_, "AUGUSTUS Training and prediction") foreach qw/new_species.pl gff2gbSmallDNA.pl filterGenes.pl randomSplit.pl etraining optimize_augustus.pl augustus/;
    $need_cmd->("diamond", "gene model picking in Step 4 and BLASTP validation in Step 5");

    # Step 5: 基因模型验证（HMMER的hmmscan）
    $need_cmd->("hmmscan", "HMMER, --HMM_db") if $HMM_db;

    # Step 6: BUSCO
    $need_cmd->("busco", "--BUSCO_lineage_dataset") if @BUSCO_db;

    print STDERR "============================================\n";
    die "Error: the following required softwares or programs are missing:\n" . join("\n", @missing) . "\n" if @missing;
    return 1;
}

# 子程序，删除不重要的中间文件，仅保留最少的、较小的、重要的中间结果文件。
sub delete_unimportant_intermediate_files {
    print STDERR "\nDeleting the unimportant intermediate files ...\n";
    my @dirs = (
        "$tmp_dir/1.RepeatMasker/repeatMasker_Dfam",
        "$tmp_dir/1.RepeatMasker/repeatMasker_RepBase",
        glob("$tmp_dir/1.RepeatMasker/repeatModeler/RM_*/round-*"),
        "$tmp_dir/2.homolog_prediction/a.MMseqs2CalHits",
        "$tmp_dir/2.homolog_prediction/b.hitToGenePrediction",
        "$tmp_dir/2.homolog_prediction/d.miniprot",
        "$tmp_dir/2.homolog_prediction/e.spaln",
        "$tmp_dir/2.homolog_prediction/f.merge.workers",
        "$tmp_dir/3.NGSReads_prediction/a.RNASeq2Sam",
        "$tmp_dir/3.NGSReads_prediction/b.Sam2Transfrag",
        "$tmp_dir/3.NGSReads_prediction/c.stringtie",
        "$tmp_dir/3.NGSReads_prediction/d.Transfrag2ORF_Sam2Transfrag",
        "$tmp_dir/3.NGSReads_prediction/d.Transfrag2ORF_stringtie",
        "$tmp_dir/4.augustus/aug_para_with_hints",
        "$tmp_dir/5.curate_gene_models/04.ValidateGenes.tmp",
        "$tmp_dir/5.curate_gene_models/06.AlternativeSplicing.tmp",
    );
    my @files = (
        glob("$tmp_dir/1.RepeatMasker/repeatModeler/species.*"),
        "$tmp_dir/2.homolog_prediction/mmseqs.filteredHits.tab",
        "$tmp_dir/2.homolog_prediction/mmseqs.geneModels.tab",
        "$tmp_dir/3.NGSReads_prediction/results/alignment.bam",
        "$tmp_dir/3.NGSReads_prediction/results/alignment.bam.bai",
    );
    foreach my $d (@dirs)  { remove_tree($d) if -d $d; }
    foreach my $f (@files) { unlink $f if -e $f; }
    return 1;
}

########################################################################
# 英文用法说明（放在文件最末尾，运行时通过DATA文件句柄读取）
########################################################################
__DATA__
Usage:
    perl __PROG__ [options]

Description:
    GETA (Genome-wide Electronic Tool for Annotation) is a pipeline software for predicting gene models from whole genome sequences. With one command, you can quickly obtain an accurate genome annotation GFF3 result file by providing the genome sequence, RNA-Seq raw data, and whole genome homologous protein sequences of closely related species. This program has two outstanding features: (1) accurate prediction, GETA outputs gene models with the normal amount, high BUSCO integrity, and accurate exon boundaries; (2) simple run, the operation is straightforward with a fully automated command that produces the final result.
    Current Version: 2.7.1

Example:
    perl __PROG__ --genome genome.fasta --RM_species_Dfam Embryophyta --RM_species_RepBase Embryophyta --pe1 libA.1.fq.gz,libB.1.fq.gz --pe2 libA.2.fq.gz,libB.2.fq.gz --protein homolog.fasta --augustus_species GETA_genus_species --HMM_db /opt/biosoft/bioinfomatics_databases/Pfam/PfamA,/opt/biosoft/bioinfomatics_databases/Pfam/PfamB --config /opt/biosoft/geta/config.txt --out_prefix out --gene_prefix GS01Gene --cpu 120

Workflow:
    Step 0 Preparation: check the parameters and the dependent softwares and sub-programs, tidy the genome sequence with genome_seq_clear.pl, correct the format of the homologous protein sequences with fasta_format_revising.pl, and read the configuration file (config.txt in the software installation directory by default, or the file given by --config).
    Step 1 Repeat analysis: runRepeatMaskerModeler is called to analyse the repeats of the genome with RepeatMasker (the Dfam database, the RepBase database or --RM_lib) and RepeatModeler, giving the repeat annotation and the repeat-masked genome sequence. This step is skipped when only --no_RepeatModeler is set without any other repeat analysis parameter.
    Step 2 Homolog prediction: homolog_prediction is called to predict gene models using the whole genome protein sequences of closely related species.
    Step 3 NGS reads prediction: NGSReads_prediction is called to predict gene models using RNA-Seq data, and to obtain the reliable intron information and the base depth information.
    Step 4 AUGUSTUS prediction: runAugustus is called to train AUGUSTUS with the gene models predicted from the transcripts and homologous proteins as evidence, and to predict gene models with hints.
    Step 5 Integration and filtering: curateGeneModels is called to integrate the gene models predicted by the three methods, remove the gene models located in transposable elements, validate the unreliable gene models with HMM and/or protein databases, add alternative splicing transcripts, perform ORF analysis and remove the incomplete transcripts, giving the final gene models.
    Step 6 Output: outputResults is called to output the GFF3, GTF and sequence files of the gene models, the low-quality gene models, the repeat annotation, the raw prediction results, the BUSCO assessment results and the statistics of all the steps.
    The program runs in the ${out_prefix}.tmp directory and an .ok file is generated when a step is finished. When the program is interrupted, run the same command again to skip the finished steps; to run a step again, delete the .ok files of that step and the following steps.
    By default, the whole genome sequence without repeat masking is used for the homolog, NGS reads and AUGUSTUS predictions, and the gene models overlapping with transposable element sequences are removed in the later integration and filtering; the unmasked genome is also used for the integration, filtering and output of the gene models. To make the program run faster and not to predict genes in transposable element regions at all, which gives fewer gene models in the end, add the --gene_predicted_by_masked_genome parameter so that the repeat-masked genome is used for the three predictions.

Output files:
    ${out_prefix}.geneModels.gff3, ${out_prefix}.geneModels.gtf    The final gene models (including alternative splicing transcripts); gene sequence files such as ${out_prefix}.protein.fasta are also output.
    ${out_prefix}.bestGeneModels.gff3, ${out_prefix}.bestGeneModels.gtf    The gene models with only one transcript kept for each gene.
    ${out_prefix}.geneModels_lowQuality.gff3    The low-quality gene models, including those failing the HMM or BLASTP validation, located in transposable elements, incomplete at the 5' or 3' end, or without CDS.
    ${out_prefix}.repeat.gff3, ${out_prefix}.repeat.stats, ${out_prefix}.repeat.lib, ${out_prefix}.maskedGenome.fasta    The repeat annotation, its statistics, the repeat library built by RepeatModeler and the repeat-masked genome sequence.
    ${out_prefix}.homolog_*.gff3, ${out_prefix}.NGSReads_*.gff3, ${out_prefix}.augustus_prediction.gff3, ${out_prefix}.evidence_prediction.gff3, ${out_prefix}.All_prediction.gff3    The raw gene predictions based on homologous proteins, transcriptome reads and AUGUSTUS, and those produced during the integration.
    ${out_prefix}.gene_prediction.summary    The statistics of all the steps, used to trace the reliability of the gene prediction results; it is also written to standard output.
    ${out_prefix}.tmp    The temporary directory holding the intermediate results of all the steps; its 6.output_gene_models sub-directory holds the BUSCO results.

Dependencies:
    This software has been tested and successfully run on Rocky 9.2 system using the following dependent software versions (listed in the order in which they are called by the pipeline):
    01. RepeatMasker (version: 4.1.6)
    02. HMMER (nhmmer/hmmscan, version: 3.3.2)
    03. makeblastdb/rmblastn (Version: 2.14.0)
    04. RepeatModeler (version: 2.0.5)
    05. mmseqs (version 15-6f452)
    06. gth (Version: 1.7.3)
    07. exonerate (Version: 2.2.0)
    08. genewise (version: 2.4.1)
    09. miniprot (Version: 0.18)
    10. spaln (Version: 3.0.9)
    11. java (version: jdk-20.0.1)
    12. STAR (Version: 2.7.11b)
    13. hisat2 (version: 2.1.0)
    14. samtools (version: 1.17)
    15. stringtie (Version: 2.1.3)
    16. TIdes (Version: 1.0.0)
    17. diamond (version 2.1.8)
    18. cd-hit (Version: 4.8.1-2019-0228)
    19. Barrnap (Version: 1.10.5)
    20. augustus/etraining (version: 3.5.0)
    21. busco (Version: 5.4.7)
    Among them, nhmmer of HMMER is required when the Dfam database is used for the repeat analysis, and hmmscan of HMMER is required when --HMM_db is set; RepeatModeler is required only when neither --RM_lib nor --no_RepeatModeler is set; mmseqs is required only when the gth, exonerate or genewise method is used by --homolog_prediction_method; gth, exonerate, genewise, miniprot and spaln are each required only when the corresponding method is used by --homolog_prediction_method, and they must be executable directly in the terminal; when next-generation sequencing data are used, stringtie and samtools are required, and with FASTQ data java (to run Trimmomatic-0.38 in the software installation directory) and the aligner STAR (default) or hisat2/hisat2-build (--aligner hisat2) are also required; TIdeS (the tides.py or tides command and its dependencies DIAMOND, cd-hit and Barrnap, etc.) is used to predict ORFs when homologous proteins are input, otherwise TransDecoder-v5.5.0 in the software installation directory is used; busco is required only when --BUSCO_lineage_dataset is set.
    AUGUSTUS requires the environment variable $AUGUSTUS_CONFIG_PATH, and augustus, etraining, optimize_augustus.pl, new_species.pl, gff2gbSmallDNA.pl, filterGenes.pl and randomSplit.pl must be executable directly in the terminal.

Usage notes:
    (1) The following sub-programs and scripts must be in the same directory as this program: runRepeatMaskerModeler, homolog_prediction, NGSReads_prediction, runAugustus, curateGeneModels, outputResults, GFF3Clear, genome_seq_clear.pl, fasta_format_revising.pl, gff3_to_sequences.pl, gff3ToGtf.pl and check_gene_model_integrity.pl.
    (2) The paths of all input files, databases and outputs can not contain spaces or shell special characters (only letters, digits and _-./+:@,= are allowed).
    (3) At least one of homologous proteins (--protein) and next-generation sequencing data (--pe1/--pe2, --se, --sam) must be input; at least one of --HMM_db, --BLASTP_db and --protein must be set to validate the unreliable gene models.
    (4) When no next-generation sequencing data are input, the alternative splicing analysis is not performed; when no homologous protein is input, the homolog prediction is not performed.

[Input options]
    --genome <string>    default: None, required
    Enter the FASTA file of the genome sequence that you want to annotate. If the input genome file has repeats masked, you can skip the repeat sequence masking step by removing the --RM_species_Dfam, --RM_species_RepBase, --RM_lib parameters and adding the --no_RepeatModeler parameter. In this instance, it is advised to hard mask the transposable sequences using base N and to soft mask simple or tandem repeats using lowercase characters. Note that no repeat annotation is then available, so the gene models located in transposable elements are not removed.

    --pe1 <string> --pe2 <string>    default: None
    Enter one or more pairs of FASTQ format files from Paired-End next-generation sequencing technology. This parameter supports the input of multiple pairs of FASTQ files, using commas to separate the FASTQ file paths of different libraries; the order of the files in --pe1 and --pe2 must correspond. This parameter also accepts compressed files in .gz format.

    --se <string>    default: None
    Enter one or more FASTQ format files from Single-End next-generation sequencing technology. This parameter supports the input of multiple Single-End FASTQ files, using commas to separate the FASTQ file paths of different libraries. This parameter also accepts compressed files in .gz format.

    --sam <string>    default: None
    Enter one or more SAM format files from the output of alignment software such as HISAT2. This parameter supports the input of multiple files, using commas to separate the file paths; the files must end with .sam, .bam or .cram. In addition, the program allows for the full or partial use of the three parameters --pe1/--pe2, --se, and --sam, then all of the input data are used for genome alignment to generate the transcript sequence for gene model prediction. Without --strand_specific, the data input by this parameter should keep the XS tag (e.g. aligned by HISAT2 with --dta or STAR with --outSAMstrandField intronMotif), and the reference names must be the same as those in the genome FASTA file.

    --strand_specific    default: None
    When this parameter is added, all input next-generation sequencing data are treated as strand-specific, and the program will predict gene models only on the forward strand of the transcript. This parameter only supports the standard RF (dUTP) library type, i.e. read1 maps to the minus strand of the transcript (e.g. Illumina TruSeq Stranded); all input data must have the same strandedness. When two neighboring genes overlap in the genome, strand-specific sequencing data and this parameter can help accurately estimate gene borders.

    --protein <string>    default: None
    Enter a FASTA file containing whole genome protein sequences from neighboring species. It is recommended to use whole genome homologous protein sequences from 3 ~ 10 different species. It is also recommended to modify the name of the protein sequence by appending the Species information, which begins with the species character, to the end of its original name. For example, if the protein sequence is XP_002436309.2, it will be better renamed XP_002436309_2_SpeciesSorghumBicolor. In this way, it is beneficial to retain the homologous matching results of more species in a gene region, improving gene prediction accuracy. The fasta_remove_redundancy.pl script included in GETA can be used to integrate the whole genome protein sequences from multiple species, and species information can be appended while redundancy is removed. The more species employed, the more accurate gene models may be predicted, but the computational time required increases. Note that evidence-supported gene prediction requires at least one type of homologous protein or next-generation sequencing data.

[Repeat analysis options]
    --RM_species_Dfam <string>    default: None
    Enter the name of a species or class for RepeatMakser to perform a repeat sequence analysis for the genome using the HMM data of corresponding taxonomic species in the Dfam database. The file __SOFTWARE_DIR__/RepeatMasker_lineage.txt has the values that can be provided for this parameter to represent the class of species. For example, Eukaryota is for eukaryotes, Viridiplantae is for plants, Metazoa is for animals, and Fungi is for fungi. Before attempting to enter this parameter, the RepeatMasker program needed to be installed and the Dfam database needed to be configured. Note that due to its massive size, the Dfam database has been split up into nine partitions. By default, RepeatMasker only contains the zero root partition of the Dfam database, which is suitable for species such as mammals and fungi. If necessary, consider downloading the proper Dfam database partition and configuring it to RepeatMasker. For example, the 5th partition of the Dfam database is designated for Viridiplantae, while the 7th partition is for Hymenoptera.

    --RM_species_RepBase <string>    default: None
    Enter the name of a species or class for RepeatMakser to perform a repeat sequence analysis for the genome using the nucleotide sequences of corresponding taxonomic species in the RepBase database. The file __SOFTWARE_DIR__/RepeatMasker_lineage.txt has the values that can be provided for this parameter to represent the class of species. Before attempting to enter this parameter, the RepeatMasker program needed to be installed and the RepBase database needed to be configured. Note that RepBase is no longer providing free downloads and that the most recent version of the database, 20181026, is older and contains few repetitive sequence data. This parameter is independent of --RM_species_Dfam, and the two parameters can be set to different values, which are used for the RepBase and Dfam analyses respectively.

    --RM_lib <string>    default: None
    Enter a FASTA file and use the repetitive sequence to conduct genome-wide repeat analysis. This file is usually the output of RepeatModeler software's analysis of the entire genome sequence, indicating the repeated sequences across the genome. By default, the GETA program calls RepeatModeler to look up the entire genome sequence and acquire the species' repetitive sequence database. RepeatMakser is then called to search the repeated sequences. After adding this argument, the time-consuming RepeatModler step is skipped, which may significantly reduce the running time of the program. Additionally, the software supports the simultaneous use of the --RM_species_Dfam, --RM_species_RepBase, and --RM_lib arguments, so that multiple methods can be used for repeat sequence analysis, and eventually multiple results can be combined and the result of any method can be recognized.

    --no_RepeatModeler    default: None
    When this parameter is added, the program will no longer run the RepeatModeler step. If none of --RM_species_Dfam, --RM_species_RepBase and --RM_lib is set at the same time, the input genome is regarded as repeat-masked and the whole repeat analysis and masking step is skipped.

    --gene_predicted_by_masked_genome    default: None
    By default, the program uses the input whole genome sequence (without repeat masking) to predict gene models by homologous proteins, NGS reads and AUGUSTUS, and the gene models overlapping with transposable element sequences are removed in the later integration and filtering. When this parameter is added, the program uses the repeat-masked genome sequence for the three gene predictions, so that no gene is predicted in transposable element regions; the program runs faster, but fewer gene models are obtained in the end. This parameter has no effect when no repeat analysis is performed (the input genome is already repeat-masked).

[Gene prediction options]
    --homolog_prediction_method <string>    default: miniprot
    Enter the method(s) for gene prediction using homologous proteins. It is consistent with the --method parameter of the homolog_prediction program. The value can be set to exonerate, genewise, gth, miniprot, spaln, or all. This parameter supports the input of multiple methods, separated by commas, e.g. exonerate,miniprot,spaln; if the value was set to all, all five methods are used. The methods are: (1) exonerate, genewise and gth: MMseqs2 first aligns the homologous proteins to the genome to locate the gene regions, and then the corresponding software predicts the gene structures in detail. These three methods are already integrated region by region inside the pipeline, so they are regarded as one source (mmseqs). Using them requires the installation of MMseqs2; (2) miniprot: aligns the homologous proteins to the genome directly and predicts gene models; it is regarded as one source; (3) spaln: aligns the homologous proteins to the genome directly and predicts gene models; it is regarded as another source. The program merges the most original gene models of all sources and removes redundancy, divides them into four classes A, B, C and D according to indices such as the coverage of the homologous protein, filters the low-confidence class C and D gene models according to the Supported_times (the number of times supported by all sources and all homologous proteins) of the gene lines and intron lines, and marks the gene models with insufficient support with Low_support. The more methods you use, the more computation time you consume, but the result is often better. The software of each method (gth, exonerate, genewise, miniprot, spaln, mmseqs) is only required when the corresponding method is selected.

    --aligner <string>    default: the default of NGSReads_prediction (star)
    Enter the software used to align the RNA-Seq data to the genome, star or hisat2, passed to the --aligner parameter of NGSReads_prediction. STAR must be installed to use star, and hisat2 and hisat2-build must be installed to use hisat2. When only --sam data are input, no alignment is performed and this parameter has no effect.

    --max_used_read_num <int>    default: calculated automatically
    Set the maximum number of NGSreads used by the program. When the program is given too many NGSreads, it will automatically select a certain amount of data based on the size of the genome. If you add this parameter value, it specifies the number of paired reads or single-end reads used in Paired-end sequencing or single-end sequencing, respectively. If you do not set this parameter, the program automatically calculates the number of read pairs used, which is equal to ((2 ** (log10(genome_size / 1,000,000) - 1)) * 50 M). That is, 10M genome uses up to 50M read pairs, PE150 sequencing data volume of 15G; 100M genome uses up to 100M read pairs, PE150 sequencing data volume of 30G; 1G genome uses up to 200M read pairs, PE150 sequencing data volume of 60G.

    --augustus_species <string>    default: set automatically
    When an AUGUSTUS species name is provided (only letters, digits, underscore, dot and hyphen are allowed), the program starts from an existing species model or retrains a new species model when performing AUGUSTUS Training using gene models predicted by transcripts or homologous proteins. If the input AUGUSTUS species model exists, its parameters will be optimized. If not, a new AUGUSTUS species HMM model will be trained and then its parameters will be optimized. The AUGUSTUS Training step requires the installation of AUGUSTUS software and configuration of the $AUGUSTUS_CONFIG_PATH environment variable. A species configuration folder from AUGUSTUS Training with the name provided in this parameter is generated in the 4.augustus/config/species folder of the temporary directory following the program's successful execution. If the user executing the program has write access, the produced species configuration folder can be copied to the species folder specified in $AUGUSTUS_CONFIG_PATH. If you do not enter this parameter, the program will automatically set the value of this parameter to "GETA + prefix of genome FASTA file name + date + process ID"; this name is saved in the augustus_species.txt file of the temporary directory and is reused when the program is run again after an interruption.

    --genetic_code <int>    default: 1
    Enter the genetic code, passed to all the called sub-programs, used for gene prediction, judging the integrity of gene models, filling the ends of gene models and extracting gene sequences. The values for this parameter can be found on the NCBI Genetic Codes website at: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi. The supported genetic code table numbers are: 1-6, 9-16 and 21-33.

[Gene model integration and filtering options]
    --HMM_db <string>    default: None
    Enter one or more HMM databases (processed by hmmpress), for validating and filtering gene models. This parameter supports the input of multiple HMM databases, separated by commas. The program filters those gene models that do not match in all databases when using multiple HMM databases. At least one of --HMM_db, --BLASTP_db and --protein must be set.

    --BLASTP_db <string>    default: None
    Enter one or more diamond databases, for validating and filtering gene models. This parameter supports the input of multiple diamond databases, separated by commas. The program filters those gene models that do not match in all databases when using multiple diamond databases. When neither --HMM_db nor this parameter is set, the homologous proteins provided by the --protein parameter will be used as the diamond database for filtering gene models.

    --no_alternative_splicing_analysis    default: None
    When this parameter is added, the program does not perform alternative splicing analysis. Note that GETA defaults to perform alternative splicing analysis based on intron and base sequencing depth information when NGS reads were input.

    --config <string>    default: config.txt in the software installation directory
    Enter a parameter profile path to set the detailed parameters of other commands called by this program. If this parameter is not set, config.txt in the software installation directory is used; it lists the parameters of the sub-programs with their default values and detailed comments. To change parameters, it is recommended to copy config.txt, modify the copy and input it with this parameter; the program then reads only the file given by this parameter, and the sections or parameters absent from that file use the defaults of the sub-programs themselves. The profile uses [program name] as the section title, and the lines below it are the extra command line parameters passed to that program; the supported sections are [runRepeatMaskerModeler], [homolog_prediction], [NGSReads_prediction], [runAugustus], [curateGeneModels] and [outputResults]; the content after # is a comment. The parameters set on the command line of this program have a higher priority than those of the same name in the profile. The profile is no longer chosen automatically according to the genome size. When GETA predicts an abnormal number of genes, the thresholds for filtering the gene models typically need to be adjusted: modify the profile, input it with this parameter and rerun the GETA pipeline (the .ok files of the corresponding steps need to be deleted). More parameters of the homolog_prediction program (such as --max_intron_size, --miniprot_options, --spaln_species and --spaln_options) can be written in the [homolog_prediction] section of the profile, and the parameters of curateGeneModels such as --validate_options and --alt_splicing_options can be written in the [curateGeneModels] section.

[Output options]
    --out_prefix <string>    default: out
    Enter a prefix of output files or temporary directory. The temporary directory is ${out_prefix}.tmp, which keeps the intermediate results of all the steps.

    --gene_prefix <string>    default: gene
    Enter the gene name prefix in the output GFF3 files; a gene ID consists of this prefix followed by a 6-digit number. The prefix of the low-quality gene IDs is this prefix with "gene" removed and "LqGene" appended.

    --BUSCO_lineage_dataset <string>    default: None
    Enter one or more BUSCO databases, the program will additionally perform BUSCO analysis on the protein sequences of the final gene models and of the raw predictions based on homologous proteins, transcriptome reads and AUGUSTUS and those produced during the integration, so that the gene prediction methods can be compared. This parameter supports the input of multiple BUSCO databases, separated by commas. The information contained in the __SOFTWARE_DIR__/BUSCO_lineages_list.2021-12-14.txt file can be used to choose the proper BUSCO databases. BUSCO is run in offline mode, so the databases must be downloaded in advance. Finally, the BUSCO results are exported to the 6.output_gene_models subdirectory and to the gene_prediction.summary file.

    --delete_unimportant_intermediate_files    default: None
    When this parameter is added and the program runs successfully, the insignificant intermediate files (such as the alignment BAM file and the temporary directories of RepeatMasker/RepeatModeler, MMseqs2, miniprot, spaln, AUGUSTUS prediction, HMM/BLASTP validation and alternative splicing analysis) are deleted, leaving only the minimal, small, and important intermediate result files.

[Running options]
    --cpu <int>    default: 4
    Enter the number of CPU threads used by GETA or the called programes to run; it is passed to all the called sub-programs.

[Help options]
    --help    default: None
    Display the English usage and exit.

    --chinese_help    default: None
    显示中文用法说明并退出。

