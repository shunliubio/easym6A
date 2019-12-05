#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use POSIX;

my $sampleList;
my $sampleNo;
my $configFile;
my $runName;
my $threads=1;
my $parallel;
my $method="all";
my $batch;
my $bstart;
my $bend;
my $onlybam;
my $daq;#double assembly quant
my $onlypeak;
my $rmrep;
my $rmdup;
my $keeptmp;
my $run;
my $help;
my $man;

GetOptions(
	"samplelist|s=s"     => \$sampleList,
	"sampleno|n=s"       => \$sampleNo,
	"configure|c=s"      => \$configFile,
	"runname|e=s"        => \$runName,
	"threads|t=i"        => \$threads,
	"parallel|l!"        => \$parallel,
	"method|m=s"         => \$method,
	"batch|b=s"          => \$batch,
	"bstart|x=i"         => \$bstart,
	"bend|y=i"           => \$bend,
	"onlybam|a!"         => \$onlybam,
	"daq|u!"             => \$daq,
	"onlypeak|k!"        => \$onlypeak,
	"rmrep|p!"           => \$rmrep,
	"rmdup|d!"           => \$rmdup,
	"keeptmp|f!"         => \$keeptmp,
	"run|r!"             => \$run,
	"help|h!"            => \$help,
	"man!"               => \$man
) or pod2usage(2);
pod2usage(-verbose => 0) if $help;
pod2usage(-verbose => 2) if $man;
pod2usage("$0: -ob/--onlybam and -op/--onlypeak are mutually exclusive!") && pod2usage(1) if ($onlybam && $onlypeak);


my %fileDir;
my %fileName;
my %fqFile;
my %bamFile;
my %adapter5;
my %adapter3;
my %barcode5;
my %barcode3;
my %Q33;
my %strandness;
my %fragmentLength;
my %readLength;
my %seqLayout;
my %sampleNo;
my %stringtieStrandness=("R" => "--rf", "RF" => "--rf", "F" => "--fr", "FR" => "--fr", "U" => "");
my %fqFileCheck;
my $rmRep=$rmrep ? "rmRep." : "";
my $rmDup=$rmdup ? "rmDup." : "";

open IN,"<$sampleList" or die "Can't open $sampleList:$!";
<IN>;
while (my $line=<IN>) {
	chomp($line);
	my @line=split /\t/,$line;
	$fileName{$line[18]}=$line[5];
	$fileDir{$line[5]}=$line[6];
	$fqFile{$line[5]}=$line[7];
	$bamFile{$line[5]}=$line[8];
	$adapter5{$line[5]}=$line[9];
	$adapter3{$line[5]}=$line[10];
	$barcode5{$line[5]}=$line[11];
	$barcode3{$line[5]}=$line[12];
	$Q33{$line[5]}=$line[13];
	$strandness{$line[5]}=$line[14];
	$fragmentLength{$line[5]}=$line[15];
	$readLength{$line[5]}=$line[16];
	$seqLayout{$line[5]}=$line[17];
	$sampleNo{$line[5]}=$line[18];
}
close IN;

if ($batch) {
	pod2usage("$0: The option bstart and bend are required in batch job mod!") && pod2usage(1) unless ($bstart && $bend);
	my $batchLine=1;
	open IN,"<$batch" or die "Can't open $batch:$!";
	<IN>;
	while (my $line=<IN>) {
		chomp($line);
		my @line=split /\t/,$line;
		my $sample_no=$line[0];
		my $config_file=$line[1];
		my $job_name=$line[2];
		my $peak_method=$line[3];
		my $batch_id=$line[4];
		unless ($sampleList && $batch) {
			pod2usage("$0: In batch job mode, both a sample list file and a batch list file are required!") && pod2usage(1);
		} else {
			unless ($sample_no && $config_file && $job_name && $peak_method) {
				die "In the batch list line $batchLine: incorrect settings in comma-seperated lists of sample ID, configuration files path, user-defined run names and (or) peak calling tools!\n";
			}
		}
		&workFlowPrint($sample_no,$config_file,$job_name,$peak_method) if ($batch_id >= $bstart && $batch_id <= $bend);
		$batchLine++;
	}
	close IN;
} else {
	unless ($sampleList && $sampleNo && $configFile && $runName) {
		pod2usage("$0: In single job mode, a sample list file, a configuration file, a comma-seperated list of sample ID and a user-defined run name are required!") && pod2usage(1);
	}
	&workFlowPrint($sampleNo,$configFile,$runName,$method);
}


###################################################sub definition##############################################################


sub workFlowPrint {
	my $fileList=shift;
	my $configure_file=shift;
	my $run_name=shift;
	my $tools=shift;
	my @fileList=split /,/,$fileList;
	my $control_input_fq_file;
	my $control_ip_fq_file;
	my $treatment_input_fq_file;
	my $treatment_ip_fq_file;
	my $control_input_bam_bash;
	my $control_ip_bam_bash;
	my $treatment_input_bam_bash;
	my $treatment_ip_bam_bash;
	my $bash_bam_body="";
	my $bash_peak_body;
	my $bash_whole;
	#read configuration file
	my $bash_log_dir;
	my $cutadapt_out_dir;
	my $hisat2_index_repBase;
	my $hisat2_index_genome;
	my $hisat2_out_dir;
	my $stringtie_out_dir;
	my $chrom_size;
	my $tx_size;
	my $gtf_file;
	my $bed12_file;
	my $peak_out_dir;
	my $genome_fasta_file;
	open INPUT,"<$configure_file" or die "Can't open $configure_file:$!";
	<INPUT>;
	while (my $line=<INPUT>) {
		chomp($line);
		my @line=split /\t/,$line;
		if ($line[0] eq "bash_log_dir") {
			$bash_log_dir=$line[1];
		} elsif ($line[0] eq "cutadapt_out_dir") {
			$cutadapt_out_dir=$line[1];
		} elsif ($line[0] eq "hisat2_index_repBase") {
			$hisat2_index_repBase=$line[1];
		} elsif ($line[0] eq "hisat2_index_genome") {
			$hisat2_index_genome=$line[1];
		} elsif ($line[0] eq "hisat2_out_dir") {
			$hisat2_out_dir=$line[1];
		} elsif ($line[0] eq "stringtie_out_dir") {
			$stringtie_out_dir=$line[1];
		} elsif ($line[0] eq "chrom_size") {
			$chrom_size=$line[1];
		} elsif ($line[0] eq "transcriptome_size") {
			$tx_size=$line[1];
		} elsif ($line[0] eq "gtf_file") {
			$gtf_file=$line[1];
		} elsif ($line[0] eq "bed12_file") {
			$bed12_file=$line[1];
		} elsif ($line[0] eq "peak_out_dir") {
			$peak_out_dir=$line[1];
		} elsif ($line[0] eq "genome_fasta_file") {
			$genome_fasta_file=$line[1];
		}
	}
	close INPUT;
	my $hisat2_index_name_repBase=(split /\//,$hisat2_index_repBase)[-1];
	my $hisat2_index_name_genome=(split /\//,$hisat2_index_genome)[-1];
	#check control fq files
	&check_fq_file($fileName{$fileList[0]});
	&check_fq_file($fileName{$fileList[1]});
	#sample bam bash setting
	if (! $onlypeak) {
		$control_input_bam_bash=&mapping_assembly_quant($fileName{$fileList[0]},"control_input");
		$control_ip_bam_bash=&mapping_assembly_quant($fileName{$fileList[1]},"control_ip");
		$bash_bam_body.=$control_input_bam_bash.$control_ip_bam_bash;
	}
	#if 4 files
	if (@fileList == 2) {
		#peak_calling mode
		$bash_peak_body=&call_peak($fileName{$fileList[0]},$fileName{$fileList[1]},"control_input","control_ip",$run_name,$tools) if ! $onlybam;
	} elsif (@fileList == 4) {
		#diff_peak_calling mode
		#check treatment fq files
		&check_fq_file($fileName{$fileList[2]});
		&check_fq_file($fileName{$fileList[3]});
		if (! $onlypeak) {
			$treatment_input_bam_bash=&mapping_assembly_quant($fileName{$fileList[2]},"treatment_input");
			$treatment_ip_bam_bash=&mapping_assembly_quant($fileName{$fileList[3]},"treatment_ip");
			$bash_bam_body.=$treatment_input_bam_bash.$treatment_ip_bam_bash;
		}
		$bash_peak_body=&call_diff_peak($fileName{$fileList[0]},$fileName{$fileList[1]},$fileName{$fileList[2]},$fileName{$fileList[3]},"control_input","control_ip","treatment_input","treatment_ip",$run_name,$tools)  if ! $onlybam;
	} else {
		die "Only 2 or 4 samples should be used:$!";
	}
	#make bash log dir
	if (! -e $bash_log_dir) {
		mkdir($bash_log_dir,0700) or die "Can't create $bash_log_dir:$!";
	}
	#bash main setting
	my $bash_header="#!/bin/bash

#SBATCH --job-name=m6Aseq
#SBATCH --output=$bash_log_dir/$run_name.log
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --ntasks-per-node=24
#SBATCH --mem-per-cpu=2000


######################
# Begin work section #
######################


echo Starting Time is \`date \"+%Y-%m-%d %H:%M:%S\"\`
start=\$(date +%s)
";
	my $bash_tailer="
wait

echo Ending Time is \`date \"+%Y-%m-%d %H:%M:%S\"\`
end=\$(date +%s)
time=\$(( (\$end - \$start) / 60 ))
echo Used Time is \$time mins
";
	my $bash_main_setting="
#main setting
cutadapt_out_dir=$cutadapt_out_dir
hisat2_index_name_rep=$hisat2_index_name_repBase
hisat2_index_rep=$hisat2_index_repBase
hisat2_index_name=$hisat2_index_name_genome
hisat2_index=$hisat2_index_genome
hisat2_out_dir=$hisat2_out_dir
stringtie_out_dir=$stringtie_out_dir
gtf_file=$gtf_file
bed12_file=$bed12_file
chrom_size=$chrom_size
tx_size=$tx_size
peak_out_dir=$peak_out_dir
genome_fa=$genome_fasta_file
ncpus=$threads

";
	if ($onlybam) {
		if ($bash_bam_body) {
			$bash_whole=$bash_header.$bash_main_setting.$bash_bam_body.$bash_tailer;
		} else {
			die "can't build bam bash lines!\n";
		}
	} elsif ($onlypeak) {
		if ($bash_peak_body) {
			$bash_whole=$bash_header.$bash_main_setting.$bash_peak_body.$bash_tailer;
		} else {
			die "can't build peak bash lines!\n";
		}
	} else {
		if ($bash_bam_body && $bash_peak_body) {
			$bash_whole=$bash_header.$bash_main_setting.$bash_bam_body.$bash_peak_body.$bash_tailer;
		} else {
			die "can't build bam and peak bash lines!\n";
		}
	}
	open OUT,">$bash_log_dir/${run_name}.sh" or die "Can't write to $bash_log_dir/${run_name}.sh:$!";
	print OUT $bash_whole;
	close OUT;
	#if ($run) {
	#	system("nohup bash $bash_log_dir/${run_name}.sh > $bash_log_dir/${run_name}.log 2>&1 &");
	#	my $pid=`ps ux | grep ${run_name}.sh | grep -v grep | awk {'print \$2'}`;
	#	print "$bash_log_dir/${run_name}.sh $pid";
	#}
	system("sbatch $bash_log_dir/${run_name}.sh") if $run;
}


sub mapping_assembly_quant {
	my $fileName=shift;
	my $label=shift;
	my $samflag=$seqLayout{$fileName} eq "SINGLE" ? "-F 1548" : "-F 1548 -f 2";
	my $run_cutadapt;
	my $assembly_quant;
	my $keepTmp;
	my $bash_link;
	my $hisat2_strandness_option;
	#control seqLayout, Q33 and adaprter setting
	my ($Q33_cutadapt,$Q33_hisat2,$adapter3,$hisat2_rmRep_IO,$hisat2_mapping_IO)=&seqLayout_adapter_Q33($fileName,$label,$strandness{$fileName});
	my $bash_sample_setting="
#${label} -- mapping, assembly and quantification
${label}=$fileName
${label}_read_length=$readLength{$fileName}
${label}_hisat2_strandness=\"$strandness{$fileName}\"
${label}_stringtie_strandness=\"$stringtieStrandness{$strandness{$fileName}}\"
";
	if (exists($fqFileCheck{$fileName}[1])) {
		$bash_sample_setting.="
${label}_r1_fq_file=$fileDir{$fileName}/$fqFileCheck{$fileName}[0]
${label}_r2_fq_file=$fileDir{$fileName}/$fqFileCheck{$fileName}[1]

echo -e \"\\n\$$label\\n\"
";
		if ($strandness{$fileName} eq "RF" || $strandness{$fileName} eq "FR") {
			$hisat2_strandness_option="--rna-strandness \$${label}_hisat2_strandness";
		} elsif ($strandness{$fileName} eq "U") {
			$hisat2_strandness_option="";
		} else {
			die "sample $fileName: strandness must be FR, RF or U for PE data!\n";
		}
	} else {
		$bash_sample_setting.="
${label}_r1_fq_file=$fileDir{$fileName}/$fqFileCheck{$fileName}[0]

echo -e \"\\n\$$label\\n\"
";
		if ($strandness{$fileName} eq "R" || $strandness{$fileName} eq "F") {
			$hisat2_strandness_option="--rna-strandness \$${label}_hisat2_strandness";
		} elsif ($strandness{$fileName} eq "U") {
			$hisat2_strandness_option="";
		} else {
			die "sample $fileName: strandness must be F, R or U for SE data!\n";
		}
	}
	if ($adapter3) {
		$run_cutadapt="
echo -e \"\\nadapter trimming\\n\"
if [[ ! -d \$cutadapt_out_dir/\$$label ]];then mkdir -p \$cutadapt_out_dir/\$$label;fi
time cutadapt -e 0.1 -n 1 -O 1 -q 10 -m 16 $Q33_cutadapt $adapter3 > \$cutadapt_out_dir/\$$label/\$$label.cutadapt.log
";
	} else {
		$run_cutadapt="";
	}
	my $mapping="
if [[ -d \$hisat2_out_dir/\$$label ]];then rm -rf \$hisat2_out_dir/\$$label;fi
mkdir -p \$hisat2_out_dir/\$$label
";
	if ($rmrep) {
		$mapping.="
echo -e \"\\nhisat2 align -- remove repetitive elements\\n\"
time hisat2 --no-spliced-alignment --no-softclip --norc -p \$ncpus --time --reorder --no-unal $hisat2_strandness_option -x \$hisat2_index_rep $hisat2_rmRep_IO | \\
	samtools view $samflag -Shub - | samtools sort -T \$hisat2_out_dir/\$$label -@ \$ncpus -o \$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name_rep.align.sorted.bam -
";
	}
	$mapping.="
echo -e \"\\nhisat2 align -- genome mapping\\n\"
time hisat2 -p \$ncpus --time --reorder --dta --no-unal --pen-noncansplice 12 $hisat2_strandness_option -x \$hisat2_index $hisat2_mapping_IO | \\
	samtools view -@ \$ncpus $samflag -Shub - | samtools sort -T \$hisat2_out_dir/\$$label -@ \$ncpus -o \$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}align.sorted.bam -

echo -e \"\\nmark duplicates\\n\"
java -Xmx4G -jar \$(which picard.jar) MarkDuplicates I=\$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}align.sorted.bam \\
	O=\$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}align.sorted.dupmark.bam M=\$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}align.sorted.dupmark.dup.qc \\
	VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
";
	if ($seqLayout{$fileName} eq "SINGLE") {
		$mapping.="
echo -e \"\\nPBC file output\\n\"
bamToBed -i \$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}align.sorted.dupmark.bam | awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {print \$1,\$2,\$3,\$6}' | sort | uniq -c | \\
	awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} (\$1==1){m1=m1+1} (\$1==2){m2=m2+1} {m0=m0+1} {mt=mt+\$1} END{printf \"%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n\",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > \\
	\$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}align.sorted.dupmark.pbc.qc
";
	} else {
		$mapping.="
echo -e \"\\nPBC file output\\n\"
samtools sort -n -T \$hisat2_out_dir/\$$label -@ \$ncpus \$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}align.sorted.dupmark.bam | bamToBed -bedpe -i - | awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {print \$1,\$2,\$4,\$6,\$9,\$10}' | \\
	sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} (\$1==1){m1=m1+1} (\$1==2){m2=m2+1} {m0=m0+1} {mt=mt+\$1} END{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; printf \"%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n\",mt,m0,m1,m2,m0/mt,m1/m0,m1_m2}' > \\
	\$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}align.sorted.dupmark.pbc.qc
";
	}
	if ($rmdup) {
		$mapping.="
echo -e \"\\nRemove duplicates\\n\"
samtools view -@ \$ncpus $samflag -Shub \$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}align.sorted.dupmark.bam | \\
	samtools sort -T \$hisat2_out_dir/\$$label -@ \$ncpus -o \$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.bam -
";
	}
	$mapping.="rm \$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}align.sorted.dupmark.bam

echo -e \"\\nsamtools index\\n\"
time samtools index \$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.bam

libsize=\$(awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {print \$1}' \$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}align.sorted.dupmark.pbc.qc)
scalefactor=\$(awk -v p=\$libsize -F '\\t' 'BEGIN {OFS=\"\\t\";print 1000000/p}')
echo -e \"\\nThe lib size is \$libsize and the scale factor is \$scalefactor\\n\"

echo -e \"\\nbedtools genomeCoverage\\n\"
time genomeCoverageBed -scale \$scalefactor -bg -split -ibam \$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.bam > \\
	\$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.bedgraph

echo -e \"\\nsort bedgraph\\n\"
time LC_COLLATE=C sort -k1,1 -k2,2n \$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.bedgraph -o \$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.bedgraph

echo -e \"\\nbedGraphToBigWig\\n\"
time bedGraphToBigWig \$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.bedgraph \$chrom_size \$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.bw
";
	if ($label eq "control_input" || $label eq "treatment_input" || $daq) {#$label eq "control_ip" || $label eq "treatment_ip"
		$assembly_quant="
if [[ ! -d \$stringtie_out_dir ]];then mkdir -p \$stringtie_out_dir;fi
if [[ -d \$stringtie_out_dir/\$$label ]];then rm -rf \$stringtie_out_dir/\$$label;fi
mkdir -p \$stringtie_out_dir/\$$label

echo -e \"\\nstringtie de novo assembly for novel transcripts\\n\"
time stringtie \$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.bam -p \$ncpus -G \$gtf_file -o \$stringtie_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}denovo.gtf -l \$$label \$${label}_stringtie_strandness

echo -e \"\\nstringtie quantification for known transcripts\\n\"
time stringtie \$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.bam -p \$ncpus -G \$gtf_file \\
			-o \$stringtie_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}known.gtf -e -B \\
			-A \$stringtie_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}geneAbund.txt -C \$stringtie_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}readCov.gtf \$${label}_stringtie_strandness

echo -e \"\\ngene and transcript counts for DE\\n\"
echo -e \"\$$label\\t\$stringtie_out_dir/\$$label/\$$label.\${hisat2_index_name}.${rmRep}known.gtf\" > \$stringtie_out_dir/\$$label/prepDE_gtf.txt
time prepDE.py -i \$stringtie_out_dir/\$$label/prepDE_gtf.txt -g \$stringtie_out_dir/\$$label/gene_count.csv -t \$stringtie_out_dir/\$$label/transcript_count.csv -l \$${label}_read_length

echo -e \"\\nnovel transcript statistics\\n\"
cd \$stringtie_out_dir/\$$label
time gffcompare -r \$gtf_file -G \$$label.\$hisat2_index_name.${rmRep}denovo.gtf

echo -e \"\\n\$$label mapping, assembly and quantification done!\\n\"
";
	} else {
		$assembly_quant="echo -e \"\\n\$$label mapping, assembly and quantification done!\\n\"
";
	}
	if ($keeptmp) {
		$keepTmp="";
	} else {
		$keepTmp="
#rm \$hisat2_out_dir/\$$label/\$$label.\$hisat2_index_name_rep.align.sorted.bam
";
		if ($rmdup) {
			$keepTmp.="
rm \$hisat2_out_dir/\$$label/*.${rmRep}fastq.gz
";
		}
		if ($label eq "control_input" || $label eq "treatment_input" || $daq) {#$label eq "control_ip" || $label eq "treatment_ip"
			$keepTmp.="rm \$stringtie_out_dir/\$$label/prepDE_gtf.txt
rm \$stringtie_out_dir/\$$label/*.ctab
rm \$stringtie_out_dir/\$$label/gffcmp.annotated.gtf
rm \$stringtie_out_dir/\$$label/gffcmp.loci
rm \$stringtie_out_dir/\$$label/gffcmp.tracking
rm \$stringtie_out_dir/\$$label/*.refmap
rm \$stringtie_out_dir/\$$label/*.tmap
rm \$stringtie_out_dir/\$$label/\$$label.\$hisat2_index_name.${rmRep}readCov.gtf
";
		}
	}
	#link bash cmd
	$bash_link=$bash_sample_setting.$run_cutadapt.$mapping.$assembly_quant.$keepTmp;
	return $bash_link;
}

sub seqLayout_adapter_Q33 {
	my $fileName=shift;
	my $label=shift;
	my $strandness=shift;
	my $Q33_cutadapt;
	my $Q33_hisat2;
	my $cutadapt_adapter;
	my $hisat2_rmRep_IO;
	my $hisat2_mapping_IO;
	if ($seqLayout{$fileName} eq "SINGLE") {
		#adapter setting
		if ($adapter3{$fileName} =~ /[A-Z]/) {
			my @adapter3=split /,/,$adapter3{$fileName};
			if (@adapter3==1) {
				if ($adapter5{$fileName} =~ /[A-Z]/) {
					my @adapter5=split /,/,$adapter5{$fileName};
					if (@adapter5==1) {
						$cutadapt_adapter="-g $adapter5[0] -a $adapter3[0] -o \$cutadapt_out_dir/\$$label/\$$label.clean.fastq.gz \$${label}_r1_fq_file";
					} elsif (@adapter5==2) {
						die "number of 5' adapter don't match seqLayout! Should be 1\n";
					} else {
						die "number of 5' adapter is either 1 for single-end or 2 for paired-end!\n";
					}
				} elsif ($adapter5{$fileName} eq "-") {
					$cutadapt_adapter="-a $adapter3[0] -o \$cutadapt_out_dir/\$$label/\$$label.clean.fastq.gz \$${label}_r1_fq_file";
				} else {
					die "sample $fileName: 5' adapter must be [A-Z] or -!\n";
				}
				if ($rmrep) {
					$hisat2_rmRep_IO="-U \$cutadapt_out_dir/\$$label/\$$label.clean.fastq.gz --un-gz \$hisat2_out_dir/\$$label/\$$label.clean.${rmRep}fastq.gz";
					$hisat2_mapping_IO="-U \$hisat2_out_dir/\$$label/\$$label.clean.${rmRep}fastq.gz --un-gz \$hisat2_out_dir/\$$label/\$$label.clean.${rmRep}unmapped.fastq.gz";
				} else {
					$hisat2_rmRep_IO="";
					$hisat2_mapping_IO="-U \$cutadapt_out_dir/\$$label/\$$label.clean.fastq.gz --un-gz \$hisat2_out_dir/\$$label/\$$label.clean.${rmRep}unmapped.fastq.gz";
				}
			} elsif (@adapter3==2) {
				die "number of 3' adapter don't match seqLayout! Should be 1\n";
			}  else {
				die "number of 3' adapter is either 1 for single-end or 2 for paired-end!\n";
			}
		} elsif ($adapter3{$fileName} eq "-") {
			if ($adapter5{$fileName} =~ /[A-Z]/) {
				my @adapter5=split /,/,$adapter5{$fileName};
				if (@adapter5==1) {
					$cutadapt_adapter="-g $adapter5[0] -o \$cutadapt_out_dir/\$$label/\$$label.clean.fastq.gz \$${label}_r1_fq_file";
					if ($rmrep) {
						$hisat2_rmRep_IO="-U \$cutadapt_out_dir/\$$label/\$$label.clean.fastq.gz --un-gz \$hisat2_out_dir/\$$label/\$$label.clean.${rmRep}fastq.gz";
						$hisat2_mapping_IO="-U \$hisat2_out_dir/\$$label/\$$label.clean.${rmRep}fastq.gz --un-gz \$hisat2_out_dir/\$$label/\$$label.clean.${rmRep}unmapped.fastq.gz";
					} else {
						$hisat2_rmRep_IO="";
						$hisat2_mapping_IO="-U \$cutadapt_out_dir/\$$label/\$$label.clean.fastq.gz --un-gz \$hisat2_out_dir/\$$label/\$$label.clean.${rmRep}unmapped.fastq.gz";
					}
				} elsif (@adapter5==2) {
					die "number of 5' adapter don't match seqLayout! Should be 1\n";
				} else {
					die "number of 5' adapter is either 1 for single-end or 2 for paired-end!\n";
				}
			} elsif ($adapter5{$fileName} eq "-") {
				$cutadapt_adapter=0;
				if ($rmrep) {
					$hisat2_rmRep_IO="-U \$${label}_r1_fq_file --un-gz \$hisat2_out_dir/\$$label/\$$label.${rmRep}fastq.gz";
					$hisat2_mapping_IO="-U \$hisat2_out_dir/\$$label/\$$label.${rmRep}fastq.gz --un-gz \$hisat2_out_dir/\$$label/\$$label.${rmRep}unmapped.fastq.gz";
				} else {
					$hisat2_rmRep_IO="";
					$hisat2_mapping_IO="-U \$${label}_r1_fq_file --un-gz \$hisat2_out_dir/\$$label/\$$label.${rmRep}unmapped.fastq.gz";
				}
			} else {
				die "sample $fileName: 5' adapter must be [A-Z] or -!\n";
			}
		} else {
			die "sample $fileName: 3' adapter must be [A-Z] or -!\n";
		}
	} elsif ($seqLayout{$fileName} eq "PAIRED") {
		if ($adapter3{$fileName} =~ /[A-Z]/) {
			my @adapter3=split /,/,$adapter3{$fileName};
			if (@adapter3==1) {
				die "number of 3' adapter don't match seqLayout! Should be 2\n";
			} elsif (@adapter3==2) {
				if (exists($fqFileCheck{$fileName}[1])) {
					if ($adapter5{$fileName} =~ /[A-Z]/) {
						my @adapter5=split /,/,$adapter5{$fileName};
						if (@adapter5==1) {
							die "number of 5' adapter don't match seqLayout! Should be 2\n";
						} elsif (@adapter5==2) {
							$cutadapt_adapter="-g $adapter5[0] -G $adapter5[1] -a $adapter3[0] -A $adapter3[1] -o \$cutadapt_out_dir/\$$label/\$$label.1.clean.fastq.gz -p \$cutadapt_out_dir/\$$label/\$$label.2.clean.fastq.gz \$${label}_r1_fq_file \$${label}_r2_fq_file";
						} else {
							die "number of 5' adapter is either 1 for single-end or 2 for paired-end!\n";
						}
					} elsif ($adapter5{$fileName} eq "-") {
						$cutadapt_adapter="-a $adapter3[0] -A $adapter3[1] -o \$cutadapt_out_dir/\$$label/\$$label.1.clean.fastq.gz -p \$cutadapt_out_dir/\$$label/\$$label.2.clean.fastq.gz \$${label}_r1_fq_file \$${label}_r2_fq_file";
					} else {
						die "sample $fileName: 5' adapter must be [A-Z] or -!\n";
					}
					if ($rmrep) {
						$hisat2_rmRep_IO="-1 \$cutadapt_out_dir/\$$label/\$$label.1.clean.fastq.gz -2 \$cutadapt_out_dir/\$$label/\$$label.2.clean.fastq.gz --un-conc-gz \$hisat2_out_dir/\$$label/\$$label.%.clean.${rmRep}fastq.gz";
						$hisat2_mapping_IO="-1 \$hisat2_out_dir/\$$label/\$$label.1.clean.${rmRep}fastq.gz -2 \$hisat2_out_dir/\$$label/\$$label.2.clean.${rmRep}fastq.gz --un-conc-gz \$hisat2_out_dir/\$$label/\$$label.%.clean.${rmRep}unmapped.fastq.gz";
					} else {
						$hisat2_rmRep_IO="";
						$hisat2_mapping_IO="-1 \$cutadapt_out_dir/\$$label/\$$label.1.clean.fastq.gz -2 \$cutadapt_out_dir/\$$label/\$$label.2.clean.fastq.gz --un-conc-gz \$hisat2_out_dir/\$$label/\$$label.%.clean.${rmRep}unmapped.fastq.gz";
					}
				} else {
					die "number of fq file don't match number of 3' or 5' adapter";
				}
			} else {
				die "number of 3' adapter is either 1 for single-end or 2 for paired-end!\n";
			}
		} elsif ($adapter3{$fileName} eq "-") {
			if ($adapter5{$fileName} =~ /[A-Z]/) {
				my @adapter5=split /,/,$adapter5{$fileName};
				if (@adapter5==1) {
					die "number of 5' adapter don't match seqLayout! Should be 2\n";
				} elsif (@adapter5==2) {
					if (exists($fqFileCheck{$fileName}[1])) {
						$cutadapt_adapter="-g $adapter5[0] -G $adapter5[1] -o \$cutadapt_out_dir/\$$label/\$$label.1.clean.fastq.gz -p \$cutadapt_out_dir/\$$label/\$$label.2.clean.fastq.gz \$${label}_r1_fq_file \$${label}_r2_fq_file";
						if ($rmrep) {
							$hisat2_rmRep_IO="-1 \$cutadapt_out_dir/\$$label/\$$label.1.clean.fastq.gz -2 \$cutadapt_out_dir/\$$label/\$$label.2.clean.fastq.gz --un-conc-gz \$hisat2_out_dir/\$$label/\$$label.%.clean.${rmRep}fastq.gz";
							$hisat2_mapping_IO="-1 \$hisat2_out_dir/\$$label/\$$label.1.clean.${rmRep}fastq.gz -2 \$hisat2_out_dir/\$$label/\$$label.2.clean.${rmRep}fastq.gz --un-conc-gz \$hisat2_out_dir/\$$label/\$$label.%.clean.${rmRep}unmapped.fastq.gz";
						} else {
							$hisat2_rmRep_IO="";
							$hisat2_mapping_IO="-1 \$cutadapt_out_dir/\$$label/\$$label.1.clean.fastq.gz -2 \$cutadapt_out_dir/\$$label/\$$label.2.clean.fastq.gz --un-conc-gz \$hisat2_out_dir/\$$label/\$$label.%.clean.${rmRep}unmapped.fastq.gz";
						}
					} else {
						die "number of fq file don't match number of 3' or 5' adapter";
					}
				} else {
					die "number of 5' adapter is either 1 for single-end or 2 for paired-end!\n";
				}
			} elsif ($adapter5{$fileName} eq "-") {
				$cutadapt_adapter=0;
				if ($rmrep) {
					$hisat2_rmRep_IO="-1 \$${label}_r1_fq_file -2 \$${label}_r2_fq_file --un-conc-gz \$hisat2_out_dir/\$$label/\$$label.%.clean.${rmRep}fastq.gz";
					$hisat2_mapping_IO="-1 \$hisat2_out_dir/\$$label/\$$label.1.clean.${rmRep}fastq.gz -2 \$hisat2_out_dir/\$$label/\$$label.2.clean.${rmRep}fastq.gz --un-conc-gz \$hisat2_out_dir/\$$label/\$$label.%.clean.${rmRep}unmapped.fastq.gz";
				} else {
					$hisat2_rmRep_IO="";
					$hisat2_mapping_IO="-1 \$${label}_r1_fq_file -2 \$${label}_r2_fq_file --un-conc-gz \$hisat2_out_dir/\$$label/\$$label.%.clean.${rmRep}unmapped.fastq.gz";
				}
			} else {
				die "sample $fileName: 5' adapter must be [A-Z] or -!\n";
			}
		} else {
			die "sample $fileName: 3' adapter must be [A-Z] or -!\n";
		}
	} else {
		die "sample $fileName: seqLayout must be SINGLE or PAIRED!\n";
	}
	#Q33 setting
	if ($Q33{$fileName} eq "Y") {
		$Q33_cutadapt="";
		$Q33_hisat2="";
	} elsif ($Q33{$fileName} eq "N") {
		$Q33_cutadapt="--quality-base 64";
		$Q33_hisat2="--phred64";
	} else {
		die "sample $fileName: Q33 label must be Y or N!\n";
	}
	return $Q33_cutadapt,$Q33_hisat2,$cutadapt_adapter,$hisat2_rmRep_IO,$hisat2_mapping_IO;
}

sub check_fq_file {
	my $fileName=shift;
	my @fqFile=split /,/,$fqFile{$fileName};
	if (@fqFile==1) {
		if (-e $fileDir{$fileName}."/".$fqFile[0]) {
			$fqFileCheck{$fileName}=[$fqFile[0]];
		} else {
			die "sample $fileName: fq file not exist!\n";
		}
	} elsif (@fqFile==2) {
		if (-e $fileDir{$fileName}."/".$fqFile[0] && -e $fileDir{$fileName}."/".$fqFile[1]) {
			$fqFileCheck{$fileName}=[$fqFile[0],$fqFile[1]];
		} else {
			die "sample $fileName: both read1 or read2 fq file need exist!\n";
		}
	} else {
		die "number of fq file is either 1 for single-end or 2 for paired-end!\n";
	}
}

sub check_bam_file {
	my $fileName=shift;
	die "sample $fileName: bam file not exist!\n" unless (-e $bamFile{$fileName});
}

sub call_peak {
	my $inputFileName=shift;
	my $ipFileName=shift;
	my $inputLabel=shift;
	my $ipLabel=shift;
	my $run_name=shift;
	my $peak_tool=shift;
	my $peak_calling_3peakSuite_divbam;
	my $peak_calling_3peakSuite_divbam_samtools;
	my $peak_calling_3peakSuite_divbam_rm;
	my $bash_link="";
	my $bam_setting="
#peak calling setting
${inputLabel}=$inputFileName
${ipLabel}=$ipFileName
compare=$run_name
read_length=$readLength{$inputFileName}
fragment_length=$fragmentLength{$inputFileName}
strandness=\"$strandness{$inputFileName}\"
";
	if ($onlypeak) {
		&check_bam_file($inputFileName);
		&check_bam_file($ipFileName);
		$bam_setting.="
hisat2_out_bam_input=$bamFile{$inputFileName}
hisat2_out_bam_ip=$bamFile{$ipFileName}
";
		if ($strandness{$inputFileName} eq "R" || $strandness{$inputFileName} eq "F" || $strandness{$inputFileName} eq "RF" || $strandness{$inputFileName} eq "FR") {
			my @suffix=qw(.bam .sam);
			my ($input_bam_name,$input_bam_dir)=fileparse($bamFile{$inputFileName},@suffix);
			my ($ip_bam_name,$ip_bam_dir)=fileparse($bamFile{$ipFileName},@suffix);
			$peak_calling_3peakSuite_divbam="
#divide bam files into plus and minus ones
hisat2_out_bam_input_plus=\$peak_out_dir/\$compare/$input_bam_name.plus.bam
hisat2_out_bam_input_minus=\$peak_out_dir/\$compare/$input_bam_name.minus.bam
hisat2_out_bam_ip_plus=\$peak_out_dir/\$compare/$ip_bam_name.plus.bam
hisat2_out_bam_ip_minus=\$peak_out_dir/\$compare/$ip_bam_name.minus.bam
";
			if ($strandness{$inputFileName} eq "R" || $strandness{$inputFileName} eq "F") {
				$peak_calling_3peakSuite_divbam_samtools="
if [[ ! -d \$peak_out_dir/\$compare ]];then mkdir -p \$peak_out_dir/\$compare;fi
samtools view -hub -f 0 -F 16 \$hisat2_out_bam_input | samtools sort -T \$peak_out_dir/\$compare -@ \$ncpus -o \$hisat2_out_bam_input_plus -
samtools view -hub -f 16 \$hisat2_out_bam_input | samtools sort -T \$peak_out_dir/\$compare -@ \$ncpus -o \$hisat2_out_bam_input_minus -
samtools view -hub -f 0 -F 16 \$hisat2_out_bam_ip | samtools sort -T \$peak_out_dir/\$compare -@ \$ncpus -o \$hisat2_out_bam_ip_plus -
samtools view -hub -f 16 \$hisat2_out_bam_ip | samtools sort -T \$peak_out_dir/\$compare -@ \$ncpus -o \$hisat2_out_bam_ip_minus -";
			} elsif ($strandness{$inputFileName} eq "RF") {
#https://www.biostars.org/p/92935/
				$peak_calling_3peakSuite_divbam_samtools="
if [[ ! -d \$peak_out_dir/\$compare ]];then mkdir -p \$peak_out_dir/\$compare;fi
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_input_plus.tmp1 \$hisat2_out_bam_input
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_input_plus.tmp2 \$hisat2_out_bam_input
samtools merge -f \$hisat2_out_bam_input_plus \$hisat2_out_bam_input_plus.tmp1 \$hisat2_out_bam_input_plus.tmp2
rm \$hisat2_out_bam_input_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_input_minus.tmp1 \$hisat2_out_bam_input
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_input_minus.tmp2 \$hisat2_out_bam_input
samtools merge -f \$hisat2_out_bam_input_minus \$hisat2_out_bam_input_minus.tmp1 \$hisat2_out_bam_input_minus.tmp2
rm \$hisat2_out_bam_input_minus.tmp[12]
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_ip_plus.tmp1 \$hisat2_out_bam_ip
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_ip_plus.tmp2 \$hisat2_out_bam_ip
samtools merge -f \$hisat2_out_bam_ip_plus \$hisat2_out_bam_ip_plus.tmp1 \$hisat2_out_bam_ip_plus.tmp2
rm \$hisat2_out_bam_ip_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_ip_minus.tmp1 \$hisat2_out_bam_ip
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_ip_minus.tmp2 \$hisat2_out_bam_ip
samtools merge -f \$hisat2_out_bam_ip_minus \$hisat2_out_bam_ip_minus.tmp1 \$hisat2_out_bam_ip_minus.tmp2
rm \$hisat2_out_bam_ip_minus.tmp[12]";
			} else {
				$peak_calling_3peakSuite_divbam_samtools="
if [[ ! -d \$peak_out_dir/\$compare ]];then mkdir -p \$peak_out_dir/\$compare;fi
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_input_plus.tmp1 \$hisat2_out_bam_input
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_input_plus.tmp2 \$hisat2_out_bam_input
samtools merge -f \$hisat2_out_bam_input_plus \$hisat2_out_bam_input_plus.tmp1 \$hisat2_out_bam_input_plus.tmp2
rm \$hisat2_out_bam_input_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_input_minus.tmp1 \$hisat2_out_bam_input
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_input_minus.tmp2 \$hisat2_out_bam_input
samtools merge -f \$hisat2_out_bam_input_minus \$hisat2_out_bam_input_minus.tmp1 \$hisat2_out_bam_input_minus.tmp2
rm \$hisat2_out_bam_input_minus.tmp[12]
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_ip_plus.tmp1 \$hisat2_out_bam_ip
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_ip_plus.tmp2 \$hisat2_out_bam_ip
samtools merge -f \$hisat2_out_bam_ip_plus \$hisat2_out_bam_ip_plus.tmp1 \$hisat2_out_bam_ip_plus.tmp2
rm \$hisat2_out_bam_ip_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_ip_minus.tmp1 \$hisat2_out_bam_ip
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_ip_minus.tmp2 \$hisat2_out_bam_ip
samtools merge -f \$hisat2_out_bam_ip_minus \$hisat2_out_bam_ip_minus.tmp1 \$hisat2_out_bam_ip_minus.tmp2
rm \$hisat2_out_bam_ip_minus.tmp[12]";
			}
			$peak_calling_3peakSuite_divbam_samtools.="
wait
samtools index \$hisat2_out_bam_input_plus
samtools index \$hisat2_out_bam_input_minus
samtools index \$hisat2_out_bam_ip_plus
samtools index \$hisat2_out_bam_ip_minus
wait
";
			$peak_calling_3peakSuite_divbam_rm="
rm \$hisat2_out_bam_input_plus
rm \$hisat2_out_bam_input_minus
rm \$hisat2_out_bam_ip_plus
rm \$hisat2_out_bam_ip_minus
rm \$hisat2_out_bam_input_plus.bai
rm \$hisat2_out_bam_input_minus.bai
rm \$hisat2_out_bam_ip_plus.bai
rm \$hisat2_out_bam_ip_minus.bai
#mv \$peak_out_dir/\$compare/\$tool/plus/\$tool.Rdata \$peak_out_dir/\$compare/\$tool/\$tool.plus.Rdata
#mv \$peak_out_dir/\$compare/\$tool/minus/\$tool.Rdata \$peak_out_dir/\$compare/\$tool/\$tool.minus.Rdata
#rm -rf \$peak_out_dir/\$compare/\$tool/plus
#rm -rf \$peak_out_dir/\$compare/\$tool/minus
";
		} else {
			$peak_calling_3peakSuite_divbam="";
			$peak_calling_3peakSuite_divbam_samtools="";
			$peak_calling_3peakSuite_divbam_rm="";
		}
	} else {
		$bam_setting.="
hisat2_out_bam_input=\$hisat2_out_dir/\$$inputLabel/\$$inputLabel.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.bam
hisat2_out_bam_ip=\$hisat2_out_dir/\$$ipLabel/\$$ipLabel.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.bam
";
		if ($strandness{$inputFileName} eq "R" || $strandness{$inputFileName} eq "F" || $strandness{$inputFileName} eq "RF" || $strandness{$inputFileName} eq "FR") {
			$peak_calling_3peakSuite_divbam="
#divide bam files into plus and minus ones
hisat2_out_bam_input_plus=\$peak_out_dir/\$compare/\$$inputLabel.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.plus.bam
hisat2_out_bam_input_minus=\$peak_out_dir/\$compare/\$$inputLabel.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.minus.bam
hisat2_out_bam_ip_plus=\$peak_out_dir/\$compare/\$$ipLabel.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.plus.bam
hisat2_out_bam_ip_minus=\$peak_out_dir/\$compare/\$$ipLabel.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.minus.bam";
			if ($strandness{$inputFileName} eq "R" || $strandness{$inputFileName} eq "F") {
				$peak_calling_3peakSuite_divbam_samtools="
if [[ ! -d \$peak_out_dir/\$compare ]];then mkdir -p \$peak_out_dir/\$compare;fi
samtools view -hub -f 0 -F 16 \$hisat2_out_bam_input | samtools sort -T \$peak_out_dir/\$compare -@ \$ncpus -o \$hisat2_out_bam_input_plus -
samtools view -hub -f 16 \$hisat2_out_bam_input | samtools sort -T \$peak_out_dir/\$compare -@ \$ncpus -o \$hisat2_out_bam_input_minus -
samtools view -hub -f 0 -F 16 \$hisat2_out_bam_ip | samtools sort -T \$peak_out_dir/\$compare -@ \$ncpus -o \$hisat2_out_bam_ip_plus -
samtools view -hub -f 16 \$hisat2_out_bam_ip | samtools sort -T \$peak_out_dir/\$compare -@ \$ncpus -o \$hisat2_out_bam_ip_minus -";
			} elsif ($strandness{$inputFileName} eq "RF") {
				$peak_calling_3peakSuite_divbam_samtools="
if [[ ! -d \$peak_out_dir/\$compare ]];then mkdir -p \$peak_out_dir/\$compare;fi
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_input_plus.tmp1 \$hisat2_out_bam_input
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_input_plus.tmp2 \$hisat2_out_bam_input
samtools merge -f \$hisat2_out_bam_input_plus \$hisat2_out_bam_input_plus.tmp1 \$hisat2_out_bam_input_plus.tmp2
rm \$hisat2_out_bam_input_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_input_minus.tmp1 \$hisat2_out_bam_input
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_input_minus.tmp2 \$hisat2_out_bam_input
samtools merge -f \$hisat2_out_bam_input_minus \$hisat2_out_bam_input_minus.tmp1 \$hisat2_out_bam_input_minus.tmp2
rm \$hisat2_out_bam_input_minus.tmp[12]
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_ip_plus.tmp1 \$hisat2_out_bam_ip
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_ip_plus.tmp2 \$hisat2_out_bam_ip
samtools merge -f \$hisat2_out_bam_ip_plus \$hisat2_out_bam_ip_plus.tmp1 \$hisat2_out_bam_ip_plus.tmp2
rm \$hisat2_out_bam_ip_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_ip_minus.tmp1 \$hisat2_out_bam_ip
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_ip_minus.tmp2 \$hisat2_out_bam_ip
samtools merge -f \$hisat2_out_bam_ip_minus \$hisat2_out_bam_ip_minus.tmp1 \$hisat2_out_bam_ip_minus.tmp2
rm \$hisat2_out_bam_ip_minus.tmp[12]";
			} else {
				$peak_calling_3peakSuite_divbam_samtools="
if [[ ! -d \$peak_out_dir/\$compare ]];then mkdir -p \$peak_out_dir/\$compare;fi
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_input_plus.tmp1 \$hisat2_out_bam_input
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_input_plus.tmp2 \$hisat2_out_bam_input
samtools merge -f \$hisat2_out_bam_input_plus \$hisat2_out_bam_input_plus.tmp1 \$hisat2_out_bam_input_plus.tmp2
rm \$hisat2_out_bam_input_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_input_minus.tmp1 \$hisat2_out_bam_input
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_input_minus.tmp2 \$hisat2_out_bam_input
samtools merge -f \$hisat2_out_bam_input_minus \$hisat2_out_bam_input_minus.tmp1 \$hisat2_out_bam_input_minus.tmp2
rm \$hisat2_out_bam_input_minus.tmp[12]
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_ip_plus.tmp1 \$hisat2_out_bam_ip
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_ip_plus.tmp2 \$hisat2_out_bam_ip
samtools merge -f \$hisat2_out_bam_ip_plus \$hisat2_out_bam_ip_plus.tmp1 \$hisat2_out_bam_ip_plus.tmp2
rm \$hisat2_out_bam_ip_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_ip_minus.tmp1 \$hisat2_out_bam_ip
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_ip_minus.tmp2 \$hisat2_out_bam_ip
samtools merge -f \$hisat2_out_bam_ip_minus \$hisat2_out_bam_ip_minus.tmp1 \$hisat2_out_bam_ip_minus.tmp2
rm \$hisat2_out_bam_ip_minus.tmp[12]";
			}
			$peak_calling_3peakSuite_divbam_samtools.="
wait
samtools index \$hisat2_out_bam_input_plus
samtools index \$hisat2_out_bam_input_minus
samtools index \$hisat2_out_bam_ip_plus
samtools index \$hisat2_out_bam_ip_minus
wait
";
			$peak_calling_3peakSuite_divbam_rm="
rm \$hisat2_out_bam_input_plus
rm \$hisat2_out_bam_input_minus
rm \$hisat2_out_bam_ip_plus
rm \$hisat2_out_bam_ip_minus
rm \$hisat2_out_bam_input_plus.bai
rm \$hisat2_out_bam_input_minus.bai
rm \$hisat2_out_bam_ip_plus.bai
rm \$hisat2_out_bam_ip_minus.bai
#mv \$peak_out_dir/\$compare/\$tool/plus/\$tool.Rdata \$peak_out_dir/\$compare/\$tool/\$tool.plus.Rdata
#mv \$peak_out_dir/\$compare/\$tool/minus/\$tool.Rdata \$peak_out_dir/\$compare/\$tool/\$tool.minus.Rdata
#rm -rf \$peak_out_dir/\$compare/\$tool/plus
#rm -rf \$peak_out_dir/\$compare/\$tool/minus
";
		} else {
			$peak_calling_3peakSuite_divbam="";
			$peak_calling_3peakSuite_divbam_samtools="";
			$peak_calling_3peakSuite_divbam_rm="";
		}
	}
	my $peak_calling_exomePeak="";
	my $peak_calling_MeTPeak="";
	my $peak_calling_macs2="";
	if ($peak_tool eq "all") {
		$peak_calling_macs2=&macs2_peak($inputFileName,$ipFileName);
		$peak_calling_exomePeak=&three_peak_suite_peak("exomePeak",$inputFileName,$ipFileName);
		$peak_calling_MeTPeak=&three_peak_suite_peak("MeTPeak",$inputFileName,$ipFileName);
		$bash_link=$bam_setting.$peak_calling_3peakSuite_divbam.$peak_calling_3peakSuite_divbam_samtools.$peak_calling_macs2.$peak_calling_exomePeak.$peak_calling_MeTPeak.$peak_calling_3peakSuite_divbam_rm;
		return $bash_link;
	} else {
		my @tools=split /,/,$peak_tool;
		my %tool_count=();
		@tools=grep {++$tool_count{$_}<2} @tools;
		for (my $i = 0; $i < @tools; $i++) {
			if ($tools[$i] eq "exomePeak") {
				$peak_calling_exomePeak=&three_peak_suite_peak("exomePeak",$inputFileName,$ipFileName);
			} elsif ($tools[$i] eq "MeTPeak") {
				$peak_calling_MeTPeak=&three_peak_suite_peak("MeTPeak",$inputFileName,$ipFileName);
			} elsif ($tools[$i] eq "MACS2") {
				$peak_calling_macs2=&macs2_peak($inputFileName,$ipFileName);
			}
		}
		unless ($peak_calling_exomePeak eq "" && $peak_calling_MeTPeak eq "" && $peak_calling_macs2 eq "") {
			$bash_link=$bam_setting.$peak_calling_3peakSuite_divbam.$peak_calling_3peakSuite_divbam_samtools.$peak_calling_macs2.$peak_calling_exomePeak.$peak_calling_MeTPeak.$peak_calling_3peakSuite_divbam_rm;;
			return $bash_link;
		} else {
			die "unrecognized tools!\n";
		}
	}
}

sub macs2_peak {
	my $inputFileName=shift;
	my $ipFileName=shift;
	my $macs2_pipeline;
	if ($strandness{$inputFileName} eq "R" || $strandness{$inputFileName} eq "F" || $strandness{$inputFileName} eq "RF" || $strandness{$inputFileName} eq "FR") {
		$macs2_pipeline="
echo -e \"\\npeak calling -- MACS2\\n\"
echo -e \"\\n\$compare\\n\"
if [[ ! -d \$peak_out_dir/\$compare/macs2/plus ]];then mkdir -p \$peak_out_dir/\$compare/macs2/plus;fi
if [[ ! -d \$peak_out_dir/\$compare/macs2/minus ]];then mkdir -p \$peak_out_dir/\$compare/macs2/minus;fi
macs2 callpeak -t \$hisat2_out_bam_ip_plus -c \$hisat2_out_bam_input_plus -f BAM -g \$tx_size -n \$peak_out_dir/\$compare/macs2/plus/\$compare \\
	-B --SPMR --nomodel --tsize \$read_length --extsize \$fragment_length --keep-dup all > \$peak_out_dir/\$compare/macs2/run.plus.log 2>&1
macs2 callpeak -t \$hisat2_out_bam_ip_minus -c \$hisat2_out_bam_input_minus -f BAM -g \$tx_size -n \$peak_out_dir/\$compare/macs2/minus/\$compare \\
	-B --SPMR --nomodel --tsize \$read_length --extsize \$fragment_length --keep-dup all > \$peak_out_dir/\$compare/macs2/run.minus.log 2>&1";
		if ($strandness{$inputFileName} eq "F" || $strandness{$inputFileName} eq "FR" || $strandness{$inputFileName} eq "RF") {
			$macs2_pipeline.="
awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {if(ARGIND==1) {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_plus_\"a[length(a)],\$5,\"+\"} else {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_minus_\"a[length(a)],\$5,\"-\"}}' \\
	\$peak_out_dir/\$compare/macs2/plus/\${compare}_summits.bed \$peak_out_dir/\$compare/macs2/minus/\${compare}_summits.bed | sort -k 1,1V -k 2,2n > \$peak_out_dir/\$compare/macs2/\${compare}_summits.bed
awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {if(ARGIND==1) {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_plus_\"a[length(a)],\$5,\"+\",\$7,\$8,\$9,\$10} else {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_minus_\"a[length(a)],\$5,\"-\",\$7,\$8,\$9,\$10}}' \\
	\$peak_out_dir/\$compare/macs2/plus/\${compare}_peaks.narrowPeak \$peak_out_dir/\$compare/macs2/minus/\${compare}_peaks.narrowPeak | sort -k 1,1V -k 2,2n > \$peak_out_dir/\$compare/macs2/\${compare}_peaks.narrowPeak";
		} else {
			$macs2_pipeline.="
awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {if(ARGIND==1) {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_plus_\"a[length(a)],\$5,\"-\"} else {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_minus_\"a[length(a)],\$5,\"+\"}}' \\
	\$peak_out_dir/\$compare/macs2/plus/\${compare}_summits.bed \$peak_out_dir/\$compare/macs2/minus/\${compare}_summits.bed | sort -k 1,1V -k 2,2n > \$peak_out_dir/\$compare/macs2/\${compare}_summits.bed
awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {if(ARGIND==1) {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_plus_\"a[length(a)],\$5,\"-\",\$7,\$8,\$9,\$10} else {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_minus_\"a[length(a)],\$5,\"+\",\$7,\$8,\$9,\$10}}' \\
	\$peak_out_dir/\$compare/macs2/plus/\${compare}_peaks.narrowPeak \$peak_out_dir/\$compare/macs2/minus/\${compare}_peaks.narrowPeak | sort -k 1,1V -k 2,2n > \$peak_out_dir/\$compare/macs2/\${compare}_peaks.narrowPeak";
		}
		$macs2_pipeline.="
cat \$peak_out_dir/\$compare/macs2/run.plus.log \$peak_out_dir/\$compare/macs2/run.minus.log > \$peak_out_dir/\$compare/macs2/run.log";
	} elsif ($strandness{$inputFileName} eq "U") {
		$macs2_pipeline="
echo -e \"\\npeak calling -- MACS2\\n\"
echo -e \"\\n\$compare\\n\"
if [[ ! -d \$peak_out_dir/\$compare/macs2 ]];then mkdir -p \$peak_out_dir/\$compare/macs2;fi
macs2 callpeak -t \$hisat2_out_bam_ip -c \$hisat2_out_bam_input -f BAM -g \$tx_size -n \$peak_out_dir/\$compare/macs2/\$compare \\
	-B --SPMR --nomodel --tsize \$read_length --extsize \$fragment_length --keep-dup all > \$peak_out_dir/\$compare/macs2/run.log 2>&1
awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_\"a[length(a)],\$5,\".\"}' \$peak_out_dir/\$compare/macs2/\${compare}_summits.bed | sort -k 1,1V -k 2,2n > \$peak_out_dir/\$compare/macs2/\${compare}_summits_unannot.bed
intersectBed -wo -f 0.5 -a \$peak_out_dir/\$compare/macs2/\${compare}_summits_unannot.bed -b \$bed12_file | awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {print \$1,\$2,\$3,\$4,\$5,\$12}' | uniq | sort -k 1,1V -k 2,2n > \\
	\$peak_out_dir/\$compare/macs2/\${compare}_summits.bed
awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_\"a[length(a)],\$5,\$6,\$7,\$8,\$9,\$10}' \$peak_out_dir/\$compare/macs2/\${compare}_peaks.narrowPeak | sort -k 1,1V -k 2,2n > \$peak_out_dir/\$compare/macs2/\${compare}_peaks_unannot.narrowPeak";
	} else {
		die "exp $runName: strandness must be F, R, FR, RF or U!\n";
	}
	$macs2_pipeline.="
peak_macs2=`cat \$peak_out_dir/\$compare/macs2/\${compare}_summits.bed | wc -l`
echo -e \"\$peak_macs2 peaks were found\" >> \$peak_out_dir/\$compare/macs2/run.log
slopBed -b 25 -g \$chrom_size -i \$peak_out_dir/\$compare/macs2/\${compare}_summits.bed > \$peak_out_dir/\$compare/macs2/\${compare}_summits_mid50.bed
fastaFromBed -name+ -split -s -fi \$genome_fa -bed \$peak_out_dir/\$compare/macs2/\${compare}_summits_mid50.bed > \$peak_out_dir/\$compare/macs2/\${compare}_summits_mid50.fa
shuffleBed -incl \$bed12_file -seed 12345 -noOverlapping -i \$peak_out_dir/\$compare/macs2/\${compare}_summits_mid50.bed -g \$chrom_size > \$peak_out_dir/\$compare/macs2/\${compare}_summits_mid50_random.bed
fastaFromBed -name+ -split -s -fi \$genome_fa -bed \$peak_out_dir/\$compare/macs2/\${compare}_summits_mid50_random.bed > \$peak_out_dir/\$compare/macs2/\${compare}_summits_mid50_random.fa
findMotifs.pl \$peak_out_dir/\$compare/macs2/\${compare}_summits_mid50.fa fasta \$peak_out_dir/\$compare/macs2/homer -fasta \$peak_out_dir/\$compare/macs2/\${compare}_summits_mid50_random.fa \\
	-p \$ncpus -len 5,6,7,8 -S 10 -rna -dumpFasta > \$peak_out_dir/\$compare/macs2/homer_run.log 2>&1
echo -e \"\\nmacs2 done!\\n\"
";
	return $macs2_pipeline;
}

sub three_peak_suite_peak {
	my $tool=shift;
	my $inputFileName=shift;
	my $ipFileName=shift;
	my $peak_calling_3peakSuite;
	if ($strandness{$inputFileName} eq "R" || $strandness{$inputFileName} eq "F" || $strandness{$inputFileName} eq "RF" || $strandness{$inputFileName} eq "FR") {
		$peak_calling_3peakSuite="
echo -e \"\\npeak calling -- $tool\\n\"
echo -e \"\\n\$compare\\n\"
if [[ ! -d \$peak_out_dir/\$compare/$tool/plus ]];then mkdir -p \$peak_out_dir/\$compare/$tool/plus;fi
if [[ ! -d \$peak_out_dir/\$compare/$tool/minus ]];then mkdir -p \$peak_out_dir/\$compare/$tool/minus;fi
3peakSuite.R $tool \$gtf_file \$peak_out_dir/\$compare/$tool plus \$fragment_length \$read_length \$hisat2_out_bam_input_plus \$hisat2_out_bam_ip_plus > \\
	\$peak_out_dir/\$compare/$tool/run.plus.log 2>&1 &
3peakSuite.R $tool \$gtf_file \$peak_out_dir/\$compare/$tool minus \$fragment_length \$read_length \$hisat2_out_bam_input_minus \$hisat2_out_bam_ip_minus > \\
	\$peak_out_dir/\$compare/$tool/run.minus.log 2>&1
wait
";
		if ($tool eq "exomePeak") {
			$peak_calling_3peakSuite.="resfiles=(con_peak.bed con_peak.xls peak.bed peak.xls)"
		} else {
			$peak_calling_3peakSuite.="resfiles=(peak.bed peak.xls)"
		}
		if ($strandness{$inputFileName} eq "F" || $strandness{$inputFileName} eq "FR" || $strandness{$inputFileName} eq "RF") {
			$peak_calling_3peakSuite.="
for i in \${resfiles[@]}; do
	awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {if(ARGIND==1) {if(\$6 == \"+\") print \$0} else {if(\$6 == \"-\") print \$0}}' \$peak_out_dir/\$compare/$tool/plus/\$i \\
		\$peak_out_dir/\$compare/$tool/minus/\$i > \$peak_out_dir/\$compare/$tool/\$i
done";
		} else {
			$peak_calling_3peakSuite.="
for i in \${resfiles[@]}; do
	awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {if(ARGIND==1) {if(\$6 == \"-\") print \$0} else {if(\$6 == \"+\") print \$0}}' \$peak_out_dir/\$compare/$tool/plus/\$i \\
		\$peak_out_dir/\$compare/$tool/minus/\$i > \$peak_out_dir/\$compare/$tool/\$i
done";
		}
	} elsif ($strandness{$inputFileName} eq "U") {
		$peak_calling_3peakSuite="
echo -e \"\\npeak calling -- $tool\\n\"
echo -e \"\\n\$compare\\n\"
if [[ ! -d \$peak_out_dir/\$compare/$tool ]];then mkdir -p \$peak_out_dir/\$compare/$tool;fi
3peakSuite.R $tool \$gtf_file \$peak_out_dir/\$compare $tool \$fragment_length \$read_length \$hisat2_out_bam_input \$hisat2_out_bam_ip > \\
	\$peak_out_dir/\$compare/$tool/run.log 2>&1
";
	} else {
		die "exp $runName: strandness must be F, R, FR, RF or U!\n";
	}
	$peak_calling_3peakSuite.="
fastaFromBed -name+ -split -s -fi \$genome_fa -bed \$peak_out_dir/\$compare/$tool/peak.bed > \$peak_out_dir/\$compare/$tool/peak.fa
shuffleBed -incl \$bed12_file -seed 12345 -noOverlapping -i \$peak_out_dir/\$compare/$tool/peak.bed -g \$chrom_size > \$peak_out_dir/\$compare/$tool/random_peak.bed
fastaFromBed -name+ -split -s -fi \$genome_fa -bed \$peak_out_dir/\$compare/$tool/random_peak.bed > \$peak_out_dir/\$compare/$tool/random_peak.fa
findMotifs.pl \$peak_out_dir/\$compare/$tool/peak.fa fasta \$peak_out_dir/\$compare/$tool/homer -fasta \$peak_out_dir/\$compare/$tool/random_peak.fa \\
	-p \$ncpus -len 5,6,7,8 -S 10 -rna -dumpFasta > \$peak_out_dir/\$compare/$tool/homer_run.log 2>&1
echo -e \"\\n$tool done!\\n\"
";
	return $peak_calling_3peakSuite;
}

sub call_diff_peak {
	my $controlInputFileName=shift;
	my $controlIPFileName=shift;
	my $treatmentInputFileName=shift;
	my $treatmentIPFileName=shift;
	my $controlInputLabel=shift;
	my $controlIPLabel=shift;
	my $treatmentInputLabel=shift;
	my $treatmentIPLabel=shift;
	my $run_name=shift;
	my $peak_tool=shift;
	my $read_length;
	my $fragment_length;
	my $strandness;
	my $peak_calling_3peakSuite_divbam;
	my $peak_calling_3peakSuite_divbam_samtools;
	my $peak_calling_3peakSuite_divbam_rm;
	my $bash_link;
	#check consistency of read length, fragment length, strandness between control and treatment experiments
	if ($readLength{$controlInputFileName} == $readLength{$treatmentInputFileName}) {
		$read_length=$readLength{$controlInputFileName};
	} else {
		warn "read length of control and treatment experiments are not the same! Mean is used.\n";
		$read_length=($readLength{$controlInputFileName} + $readLength{$treatmentInputFileName}) / 2;
	}
	if ($fragmentLength{$controlInputFileName} == $fragmentLength{$treatmentInputFileName}) {
		$fragment_length=$fragmentLength{$controlInputFileName};
	} else {
		warn "fragment length of control and treatment experiments are not the same! Mean is used.\n";
		$fragment_length=($fragmentLength{$controlInputFileName} + $fragmentLength{$treatmentInputFileName}) / 2;
	}
	if ($strandness{$controlInputFileName} eq $strandness{$treatmentInputFileName}) {
		$strandness=$strandness{$controlInputFileName};
	} else {
		warn "strandness of control and treatment experiments are not the same! U (unstranded) is used.\n";
		$strandness="U";
	}
	my $bam_setting="
#peak calling setting
${controlInputLabel}=$controlInputFileName
${controlIPLabel}=$controlIPFileName
${treatmentInputLabel}=$treatmentInputFileName
${treatmentIPLabel}=$treatmentIPFileName
compare=$run_name
control_read_length=$readLength{$controlInputFileName}
control_fragment_length=$fragmentLength{$controlInputFileName}
treatment_read_length=$readLength{$treatmentInputFileName}
treatment_fragment_length=$fragmentLength{$treatmentInputFileName}
control_strandness=\"$strandness{$controlInputFileName}\"
treatment_strandness=\"$strandness{$treatmentInputFileName}\"
read_length=$read_length
fragment_length=$fragment_length
strandness=$strandness
";
	if ($onlypeak) {
		&check_bam_file($controlInputFileName);
		&check_bam_file($controlIPFileName);
		&check_bam_file($treatmentInputFileName);
		&check_bam_file($treatmentIPFileName);
		$bam_setting.="
hisat2_out_bam_control_input=$bamFile{$controlInputFileName}
hisat2_out_bam_control_ip=$bamFile{$controlIPFileName}
hisat2_out_bam_treatment_input=$bamFile{$treatmentInputFileName}
hisat2_out_bam_treatment_ip=$bamFile{$treatmentIPFileName}
";
		if ($strandness eq "R" || $strandness eq "F" || $strandness eq "RF" || $strandness eq "FR") {
			my @suffix=qw(.bam .sam);
			my ($control_input_bam_name,$control_input_bam_dir)=fileparse($bamFile{$controlInputFileName},@suffix);
			my ($control_ip_bam_name,$control_ip_bam_dir)=fileparse($bamFile{$controlIPFileName},@suffix);
			my ($treatment_input_bam_name,$treatment_input_bam_dir)=fileparse($bamFile{$treatmentInputFileName},@suffix);
			my ($treatment_ip_bam_name,$treatment_ip_bam_dir)=fileparse($bamFile{$treatmentIPFileName},@suffix);
			$peak_calling_3peakSuite_divbam="
#divide bam files into plus and minus ones
hisat2_out_bam_control_input_plus=\$peak_out_dir/\$compare/$control_input_bam_name.plus.bam
hisat2_out_bam_control_input_minus=\$peak_out_dir/\$compare/$control_input_bam_name.minus.bam
hisat2_out_bam_control_ip_plus=\$peak_out_dir/\$compare/$control_ip_bam_name.plus.bam
hisat2_out_bam_control_ip_minus=\$peak_out_dir/\$compare/$control_ip_bam_name.minus.bam
hisat2_out_bam_treatment_input_plus=\$peak_out_dir/\$compare/$treatment_input_bam_name.plus.bam
hisat2_out_bam_treatment_input_minus=\$peak_out_dir/\$compare/$treatment_input_bam_name.minus.bam
hisat2_out_bam_treatment_ip_plus=\$peak_out_dir/\$compare/$treatment_ip_bam_name.plus.bam
hisat2_out_bam_treatment_ip_minus=\$peak_out_dir/\$compare/$treatment_ip_bam_name.minus.bam
";
			if ($strandness eq "R" || $strandness eq "F") {
				$peak_calling_3peakSuite_divbam_samtools="
if [[ ! -d \$peak_out_dir/\$compare ]];then mkdir -p \$peak_out_dir/\$compare;fi
samtools view -hub -f 0 -F 16 \$hisat2_out_bam_control_input | samtools sort -T \$hisat2_out_dir/\$$controlInputLabel -@ \$ncpus -o \$hisat2_out_bam_control_input_plus -
samtools view -hub -f 16 \$hisat2_out_bam_control_input | samtools sort -T \$hisat2_out_dir/\$$controlInputLabel -@ \$ncpus -o \$hisat2_out_bam_control_input_minus -
samtools view -hub -f 0 -F 16 \$hisat2_out_bam_control_ip | samtools sort -T \$hisat2_out_dir/\$$controlIPLabel -@ \$ncpus -o \$hisat2_out_bam_control_ip_plus -
samtools view -hub -f 16 \$hisat2_out_bam_control_ip | samtools sort -T \$hisat2_out_dir/\$$controlIPLabel -@ \$ncpus -o \$hisat2_out_bam_control_ip_minus -
samtools view -hub -f 0 -F 16 \$hisat2_out_bam_treatment_input | samtools sort -T \$hisat2_out_dir/\$$treatmentInputLabel -@ \$ncpus -o \$hisat2_out_bam_treatment_input_plus -
samtools view -hub -f 16 \$hisat2_out_bam_treatment_input | samtools sort -T \$hisat2_out_dir/\$$treatmentInputLabel -@ \$ncpus -o \$hisat2_out_bam_treatment_input_minus -
samtools view -hub -f 0 -F 16 \$hisat2_out_bam_treatment_ip | samtools sort -T \$hisat2_out_dir/\$$treatmentIPLabel -@ \$ncpus -o \$hisat2_out_bam_treatment_ip_plus -
samtools view -hub -f 16 \$hisat2_out_bam_treatment_ip | samtools sort -T \$hisat2_out_dir/\$$treatmentIPLabel -@ \$ncpus -o \$hisat2_out_bam_treatment_ip_minus -";
			} elsif ($strandness eq "RF") {
				$peak_calling_3peakSuite_divbam_samtools="
if [[ ! -d \$peak_out_dir/\$compare ]];then mkdir -p \$peak_out_dir/\$compare;fi
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_control_input_plus.tmp1 \$hisat2_out_bam_control_input
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_control_input_plus.tmp2 \$hisat2_out_bam_control_input
samtools merge -f \$hisat2_out_bam_control_input_plus \$hisat2_out_bam_control_input_plus.tmp1 \$hisat2_out_bam_control_input_plus.tmp2
rm \$hisat2_out_bam_control_input_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_control_input_minus.tmp1 \$hisat2_out_bam_control_input
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_control_input_minus.tmp2 \$hisat2_out_bam_control_input
samtools merge -f \$hisat2_out_bam_control_input_minus \$hisat2_out_bam_control_input_minus.tmp1 \$hisat2_out_bam_control_input_minus.tmp2
rm \$hisat2_out_bam_control_input_minus.tmp[12]
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_control_ip_plus.tmp1 \$hisat2_out_bam_control_ip
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_control_ip_plus.tmp2 \$hisat2_out_bam_control_ip
samtools merge -f \$hisat2_out_bam_control_ip_plus \$hisat2_out_bam_control_ip_plus.tmp1 \$hisat2_out_bam_control_ip_plus.tmp2
rm \$hisat2_out_bam_control_ip_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_control_ip_minus.tmp1 \$hisat2_out_bam_control_ip
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_control_ip_minus.tmp2 \$hisat2_out_bam_control_ip
samtools merge -f \$hisat2_out_bam_control_ip_minus \$hisat2_out_bam_control_ip_minus.tmp1 \$hisat2_out_bam_control_ip_minus.tmp2
rm \$hisat2_out_bam_control_ip_minus.tmp[12]
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_treatment_input_plus.tmp1 \$hisat2_out_bam_treatment_input
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_treatment_input_plus.tmp2 \$hisat2_out_bam_treatment_input
samtools merge -f \$hisat2_out_bam_treatment_input_plus \$hisat2_out_bam_treatment_input_plus.tmp1 \$hisat2_out_bam_treatment_input_plus.tmp2
rm \$hisat2_out_bam_treatment_input_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_treatment_input_minus.tmp1 \$hisat2_out_bam_treatment_input
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_treatment_input_minus.tmp2 \$hisat2_out_bam_treatment_input
samtools merge -f \$hisat2_out_bam_treatment_input_minus \$hisat2_out_bam_treatment_input_minus.tmp1 \$hisat2_out_bam_treatment_input_minus.tmp2
rm \$hisat2_out_bam_treatment_input_minus.tmp[12]
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_treatment_ip_plus.tmp1 \$hisat2_out_bam_treatment_ip
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_treatment_ip_plus.tmp2 \$hisat2_out_bam_treatment_ip
samtools merge -f \$hisat2_out_bam_treatment_ip_plus \$hisat2_out_bam_treatment_ip_plus.tmp1 \$hisat2_out_bam_treatment_ip_plus.tmp2
rm \$hisat2_out_bam_treatment_ip_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_treatment_ip_minus.tmp1 \$hisat2_out_bam_treatment_ip
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_treatment_ip_minus.tmp2 \$hisat2_out_bam_treatment_ip
samtools merge -f \$hisat2_out_bam_treatment_ip_minus \$hisat2_out_bam_treatment_ip_minus.tmp1 \$hisat2_out_bam_treatment_ip_minus.tmp2
rm \$hisat2_out_bam_treatment_ip_minus.tmp[12]";
			} else {
				$peak_calling_3peakSuite_divbam_samtools="
if [[ ! -d \$peak_out_dir/\$compare ]];then mkdir -p \$peak_out_dir/\$compare;fi
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_control_input_plus.tmp1 \$hisat2_out_bam_control_input
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_control_input_plus.tmp2 \$hisat2_out_bam_control_input
samtools merge -f \$hisat2_out_bam_control_input_plus \$hisat2_out_bam_control_input_plus.tmp1 \$hisat2_out_bam_control_input_plus.tmp2
rm \$hisat2_out_bam_control_input_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_control_input_minus.tmp1 \$hisat2_out_bam_control_input
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_control_input_minus.tmp2 \$hisat2_out_bam_control_input
samtools merge -f \$hisat2_out_bam_control_input_minus \$hisat2_out_bam_control_input_minus.tmp1 \$hisat2_out_bam_control_input_minus.tmp2
rm \$hisat2_out_bam_control_input_minus.tmp[12]
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_control_ip_plus.tmp1 \$hisat2_out_bam_control_ip
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_control_ip_plus.tmp2 \$hisat2_out_bam_control_ip
samtools merge -f \$hisat2_out_bam_control_ip_plus \$hisat2_out_bam_control_ip_plus.tmp1 \$hisat2_out_bam_control_ip_plus.tmp2
rm \$hisat2_out_bam_control_ip_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_control_ip_minus.tmp1 \$hisat2_out_bam_control_ip
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_control_ip_minus.tmp2 \$hisat2_out_bam_control_ip
samtools merge -f \$hisat2_out_bam_control_ip_minus \$hisat2_out_bam_control_ip_minus.tmp1 \$hisat2_out_bam_control_ip_minus.tmp2
rm \$hisat2_out_bam_control_ip_minus.tmp[12]
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_treatment_input_plus.tmp1 \$hisat2_out_bam_treatment_input
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_treatment_input_plus.tmp2 \$hisat2_out_bam_treatment_input
samtools merge -f \$hisat2_out_bam_treatment_input_plus \$hisat2_out_bam_treatment_input_plus.tmp1 \$hisat2_out_bam_treatment_input_plus.tmp2
rm \$hisat2_out_bam_treatment_input_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_treatment_input_minus.tmp1 \$hisat2_out_bam_treatment_input
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_treatment_input_minus.tmp2 \$hisat2_out_bam_treatment_input
samtools merge -f \$hisat2_out_bam_treatment_input_minus \$hisat2_out_bam_treatment_input_minus.tmp1 \$hisat2_out_bam_treatment_input_minus.tmp2
rm \$hisat2_out_bam_treatment_input_minus.tmp[12]
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_treatment_ip_plus.tmp1 \$hisat2_out_bam_treatment_ip
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_treatment_ip_plus.tmp2 \$hisat2_out_bam_treatment_ip
samtools merge -f \$hisat2_out_bam_treatment_ip_plus \$hisat2_out_bam_treatment_ip_plus.tmp1 \$hisat2_out_bam_treatment_ip_plus.tmp2
rm \$hisat2_out_bam_treatment_ip_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_treatment_ip_minus.tmp1 \$hisat2_out_bam_treatment_ip
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_treatment_ip_minus.tmp2 \$hisat2_out_bam_treatment_ip
samtools merge -f \$hisat2_out_bam_treatment_ip_minus \$hisat2_out_bam_treatment_ip_minus.tmp1 \$hisat2_out_bam_treatment_ip_minus.tmp2
rm \$hisat2_out_bam_treatment_ip_minus.tmp[12]";
			}
			$peak_calling_3peakSuite_divbam_samtools.="
wait
samtools index \$hisat2_out_bam_control_input_plus
samtools index \$hisat2_out_bam_control_input_minus
samtools index \$hisat2_out_bam_control_ip_plus
samtools index \$hisat2_out_bam_control_ip_minus
samtools index \$hisat2_out_bam_treatment_input_plus
samtools index \$hisat2_out_bam_treatment_input_minus
samtools index \$hisat2_out_bam_treatment_ip_plus
samtools index \$hisat2_out_bam_treatment_ip_minus
wait
";
			$peak_calling_3peakSuite_divbam_rm="
rm \$hisat2_out_bam_control_input_plus
rm \$hisat2_out_bam_control_input_minus
rm \$hisat2_out_bam_control_ip_plus
rm \$hisat2_out_bam_control_ip_minus
rm \$hisat2_out_bam_control_input_plus.bai
rm \$hisat2_out_bam_control_input_minus.bai
rm \$hisat2_out_bam_control_ip_plus.bai
rm \$hisat2_out_bam_control_ip_minus.bai
rm \$hisat2_out_bam_treatment_input_plus
rm \$hisat2_out_bam_treatment_input_minus
rm \$hisat2_out_bam_treatment_ip_plus
rm \$hisat2_out_bam_treatment_ip_minus
rm \$hisat2_out_bam_treatment_input_plus.bai
rm \$hisat2_out_bam_treatment_input_minus.bai
rm \$hisat2_out_bam_treatment_ip_plus.bai
rm \$hisat2_out_bam_treatment_ip_minus.bai
#mv \$peak_out_dir/\$compare/\$tool/plus/\$tool.Rdata \$peak_out_dir/\$compare/\$tool/\$tool.plus.Rdata
#mv \$peak_out_dir/\$compare/\$tool/minus/\$tool.Rdata \$peak_out_dir/\$compare/\$tool/\$tool.minus.Rdata
#rm -rf \$peak_out_dir/\$compare/\$tool/plus
#rm -rf \$peak_out_dir/\$compare/\$tool/minus
";
		} else {
			$peak_calling_3peakSuite_divbam="";
			$peak_calling_3peakSuite_divbam_samtools="";
			$peak_calling_3peakSuite_divbam_rm="";
		}
	} else {
		$bam_setting.="
hisat2_out_bam_control_input=\$hisat2_out_dir/\$$controlInputLabel/\$$controlInputLabel.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.bam
hisat2_out_bam_control_ip=\$hisat2_out_dir/\$$controlIPLabel/\$$controlIPLabel.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.bam
hisat2_out_bam_treatment_input=\$hisat2_out_dir/\$$treatmentInputLabel/\$$treatmentInputLabel.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.bam
hisat2_out_bam_treatment_ip=\$hisat2_out_dir/\$$treatmentIPLabel/\$$treatmentIPLabel.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.bam
";
		if ($strandness eq "R" || $strandness eq "F" || $strandness eq "RF" || $strandness eq "FR") {
			$peak_calling_3peakSuite_divbam="
#divide bam files into plus and minus ones
hisat2_out_bam_control_input_plus=\$hisat2_out_dir/\$$controlInputLabel/\$$controlInputLabel.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.plus.bam
hisat2_out_bam_control_input_minus=\$hisat2_out_dir/\$$controlInputLabel/\$$controlInputLabel.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.minus.bam
hisat2_out_bam_control_ip_plus=\$hisat2_out_dir/\$$controlIPLabel/\$$controlIPLabel.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.plus.bam
hisat2_out_bam_control_ip_minus=\$hisat2_out_dir/\$$controlIPLabel/\$$controlIPLabel.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.minus.bam
hisat2_out_bam_treatment_input_plus=\$hisat2_out_dir/\$$treatmentInputLabel/\$$treatmentInputLabel.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.plus.bam
hisat2_out_bam_treatment_input_minus=\$hisat2_out_dir/\$$treatmentInputLabel/\$$treatmentInputLabel.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.minus.bam
hisat2_out_bam_treatment_ip_plus=\$hisat2_out_dir/\$$treatmentIPLabel/\$$treatmentIPLabel.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.plus.bam
hisat2_out_bam_treatment_ip_minus=\$hisat2_out_dir/\$$treatmentIPLabel/\$$treatmentIPLabel.\$hisat2_index_name.${rmRep}align.${rmDup}sorted.minus.bam";
			if ($strandness eq "R" || $strandness eq "F") {
				$peak_calling_3peakSuite_divbam_samtools="
if [[ ! -d \$peak_out_dir/\$compare ]];then mkdir -p \$peak_out_dir/\$compare;fi
samtools view -hub -f 0 -F 16 \$hisat2_out_bam_control_input | samtools sort -T \$hisat2_out_dir/\$$controlInputLabel -@ \$ncpus -o \$hisat2_out_bam_control_input_plus -
samtools view -hub -f 16 \$hisat2_out_bam_control_input | samtools sort -T \$hisat2_out_dir/\$$controlInputLabel -@ \$ncpus -o \$hisat2_out_bam_control_input_minus -
samtools view -hub -f 0 -F 16 \$hisat2_out_bam_control_ip | samtools sort -T \$hisat2_out_dir/\$$controlIPLabel -@ \$ncpus -o \$hisat2_out_bam_control_ip_plus -
samtools view -hub -f 16 \$hisat2_out_bam_control_ip | samtools sort -T \$hisat2_out_dir/\$$controlIPLabel -@ \$ncpus -o \$hisat2_out_bam_control_ip_minus -
samtools view -hub -f 0 -F 16 \$hisat2_out_bam_treatment_input | samtools sort -T \$hisat2_out_dir/\$$treatmentInputLabel -@ \$ncpus -o \$hisat2_out_bam_treatment_input_plus -
samtools view -hub -f 16 \$hisat2_out_bam_treatment_input | samtools sort -T \$hisat2_out_dir/\$$treatmentInputLabel -@ \$ncpus -o \$hisat2_out_bam_treatment_input_minus -
samtools view -hub -f 0 -F 16 \$hisat2_out_bam_treatment_ip | samtools sort -T \$hisat2_out_dir/\$$treatmentIPLabel -@ \$ncpus -o \$hisat2_out_bam_treatment_ip_plus -
samtools view -hub -f 16 \$hisat2_out_bam_treatment_ip | samtools sort -T \$hisat2_out_dir/\$$treatmentIPLabel -@ \$ncpus -o \$hisat2_out_bam_treatment_ip_minus -";
			} elsif ($strandness eq "RF") {
				$peak_calling_3peakSuite_divbam_samtools="
if [[ ! -d \$peak_out_dir/\$compare ]];then mkdir -p \$peak_out_dir/\$compare;fi
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_control_input_plus.tmp1 \$hisat2_out_bam_control_input
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_control_input_plus.tmp2 \$hisat2_out_bam_control_input
samtools merge -f \$hisat2_out_bam_control_input_plus \$hisat2_out_bam_control_input_plus.tmp1 \$hisat2_out_bam_control_input_plus.tmp2
rm \$hisat2_out_bam_control_input_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_control_input_minus.tmp1 \$hisat2_out_bam_control_input
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_control_input_minus.tmp2 \$hisat2_out_bam_control_input
samtools merge -f \$hisat2_out_bam_control_input_minus \$hisat2_out_bam_control_input_minus.tmp1 \$hisat2_out_bam_control_input_minus.tmp2
rm \$hisat2_out_bam_control_input_minus.tmp[12]
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_control_ip_plus.tmp1 \$hisat2_out_bam_control_ip
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_control_ip_plus.tmp2 \$hisat2_out_bam_control_ip
samtools merge -f \$hisat2_out_bam_control_ip_plus \$hisat2_out_bam_control_ip_plus.tmp1 \$hisat2_out_bam_control_ip_plus.tmp2
rm \$hisat2_out_bam_control_ip_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_control_ip_minus.tmp1 \$hisat2_out_bam_control_ip
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_control_ip_minus.tmp2 \$hisat2_out_bam_control_ip
samtools merge -f \$hisat2_out_bam_control_ip_minus \$hisat2_out_bam_control_ip_minus.tmp1 \$hisat2_out_bam_control_ip_minus.tmp2
rm \$hisat2_out_bam_control_ip_minus.tmp[12]
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_treatment_input_plus.tmp1 \$hisat2_out_bam_treatment_input
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_treatment_input_plus.tmp2 \$hisat2_out_bam_treatment_input
samtools merge -f \$hisat2_out_bam_treatment_input_plus \$hisat2_out_bam_treatment_input_plus.tmp1 \$hisat2_out_bam_treatment_input_plus.tmp2
rm \$hisat2_out_bam_treatment_input_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_treatment_input_minus.tmp1 \$hisat2_out_bam_treatment_input
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_treatment_input_minus.tmp2 \$hisat2_out_bam_treatment_input
samtools merge -f \$hisat2_out_bam_treatment_input_minus \$hisat2_out_bam_treatment_input_minus.tmp1 \$hisat2_out_bam_treatment_input_minus.tmp2
rm \$hisat2_out_bam_treatment_input_minus.tmp[12]
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_treatment_ip_plus.tmp1 \$hisat2_out_bam_treatment_ip
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_treatment_ip_plus.tmp2 \$hisat2_out_bam_treatment_ip
samtools merge -f \$hisat2_out_bam_treatment_ip_plus \$hisat2_out_bam_treatment_ip_plus.tmp1 \$hisat2_out_bam_treatment_ip_plus.tmp2
rm \$hisat2_out_bam_treatment_ip_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_treatment_ip_minus.tmp1 \$hisat2_out_bam_treatment_ip
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_treatment_ip_minus.tmp2 \$hisat2_out_bam_treatment_ip
samtools merge -f \$hisat2_out_bam_treatment_ip_minus \$hisat2_out_bam_treatment_ip_minus.tmp1 \$hisat2_out_bam_treatment_ip_minus.tmp2
rm \$hisat2_out_bam_treatment_ip_minus.tmp[12]";
			} else {
				$peak_calling_3peakSuite_divbam_samtools="
if [[ ! -d \$peak_out_dir/\$compare ]];then mkdir -p \$peak_out_dir/\$compare;fi
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_control_input_plus.tmp1 \$hisat2_out_bam_control_input
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_control_input_plus.tmp2 \$hisat2_out_bam_control_input
samtools merge -f \$hisat2_out_bam_control_input_plus \$hisat2_out_bam_control_input_plus.tmp1 \$hisat2_out_bam_control_input_plus.tmp2
rm \$hisat2_out_bam_control_input_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_control_input_minus.tmp1 \$hisat2_out_bam_control_input
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_control_input_minus.tmp2 \$hisat2_out_bam_control_input
samtools merge -f \$hisat2_out_bam_control_input_minus \$hisat2_out_bam_control_input_minus.tmp1 \$hisat2_out_bam_control_input_minus.tmp2
rm \$hisat2_out_bam_control_input_minus.tmp[12]
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_control_ip_plus.tmp1 \$hisat2_out_bam_control_ip
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_control_ip_plus.tmp2 \$hisat2_out_bam_control_ip
samtools merge -f \$hisat2_out_bam_control_ip_plus \$hisat2_out_bam_control_ip_plus.tmp1 \$hisat2_out_bam_control_ip_plus.tmp2
rm \$hisat2_out_bam_control_ip_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_control_ip_minus.tmp1 \$hisat2_out_bam_control_ip
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_control_ip_minus.tmp2 \$hisat2_out_bam_control_ip
samtools merge -f \$hisat2_out_bam_control_ip_minus \$hisat2_out_bam_control_ip_minus.tmp1 \$hisat2_out_bam_control_ip_minus.tmp2
rm \$hisat2_out_bam_control_ip_minus.tmp[12]
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_treatment_input_plus.tmp1 \$hisat2_out_bam_treatment_input
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_treatment_input_plus.tmp2 \$hisat2_out_bam_treatment_input
samtools merge -f \$hisat2_out_bam_treatment_input_plus \$hisat2_out_bam_treatment_input_plus.tmp1 \$hisat2_out_bam_treatment_input_plus.tmp2
rm \$hisat2_out_bam_treatment_input_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_treatment_input_minus.tmp1 \$hisat2_out_bam_treatment_input
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_treatment_input_minus.tmp2 \$hisat2_out_bam_treatment_input
samtools merge -f \$hisat2_out_bam_treatment_input_minus \$hisat2_out_bam_treatment_input_minus.tmp1 \$hisat2_out_bam_treatment_input_minus.tmp2
rm \$hisat2_out_bam_treatment_input_minus.tmp[12]
samtools view -@ \$ncpus -hub -f 144 -o \$hisat2_out_bam_treatment_ip_plus.tmp1 \$hisat2_out_bam_treatment_ip
samtools view -@ \$ncpus -hub -f 64 -F 16 -o \$hisat2_out_bam_treatment_ip_plus.tmp2 \$hisat2_out_bam_treatment_ip
samtools merge -f \$hisat2_out_bam_treatment_ip_plus \$hisat2_out_bam_treatment_ip_plus.tmp1 \$hisat2_out_bam_treatment_ip_plus.tmp2
rm \$hisat2_out_bam_treatment_ip_plus.tmp[12]
samtools view -@ \$ncpus -hub -f 128 -F 16 -o \$hisat2_out_bam_treatment_ip_minus.tmp1 \$hisat2_out_bam_treatment_ip
samtools view -@ \$ncpus -hub -f 80 -o \$hisat2_out_bam_treatment_ip_minus.tmp2 \$hisat2_out_bam_treatment_ip
samtools merge -f \$hisat2_out_bam_treatment_ip_minus \$hisat2_out_bam_treatment_ip_minus.tmp1 \$hisat2_out_bam_treatment_ip_minus.tmp2
rm \$hisat2_out_bam_treatment_ip_minus.tmp[12]";
			}
			$peak_calling_3peakSuite_divbam_samtools.="
wait
samtools index \$hisat2_out_bam_control_input_plus
samtools index \$hisat2_out_bam_control_input_minus
samtools index \$hisat2_out_bam_control_ip_plus
samtools index \$hisat2_out_bam_control_ip_minus
samtools index \$hisat2_out_bam_treatment_input_plus
samtools index \$hisat2_out_bam_treatment_input_minus
samtools index \$hisat2_out_bam_treatment_ip_plus
samtools index \$hisat2_out_bam_treatment_ip_minu
wait
";
			$peak_calling_3peakSuite_divbam_rm="
rm \$hisat2_out_bam_control_input_plus
rm \$hisat2_out_bam_control_input_minus
rm \$hisat2_out_bam_control_ip_plus
rm \$hisat2_out_bam_control_ip_minus
rm \$hisat2_out_bam_control_input_plus.bai
rm \$hisat2_out_bam_control_input_minus.bai
rm \$hisat2_out_bam_control_ip_plus.bai
rm \$hisat2_out_bam_control_ip_minus.bai
rm \$hisat2_out_bam_treatment_input_plus
rm \$hisat2_out_bam_treatment_input_minus
rm \$hisat2_out_bam_treatment_ip_plus
rm \$hisat2_out_bam_treatment_ip_minus
rm \$hisat2_out_bam_treatment_input_plus.bai
rm \$hisat2_out_bam_treatment_input_minus.bai
rm \$hisat2_out_bam_treatment_ip_plus.bai
rm \$hisat2_out_bam_treatment_ip_minus.bai
#mv \$peak_out_dir/\$compare/\$tool/plus/\$tool.Rdata \$peak_out_dir/\$compare/\$tool/\$tool.plus.Rdata
#mv \$peak_out_dir/\$compare/\$tool/minus/\$tool.Rdata \$peak_out_dir/\$compare/\$tool/\$tool.minus.Rdata
#rm -rf \$peak_out_dir/\$compare/\$tool/plus
#rm -rf \$peak_out_dir/\$compare/\$tool/minus
";
		} else {
			$peak_calling_3peakSuite_divbam="";
			$peak_calling_3peakSuite_divbam_samtools="";
			$peak_calling_3peakSuite_divbam_rm="";
		}
	}
	my $peak_calling_exomePeak="";
	my $peak_calling_MeTDiff="";
	my $peak_calling_macs2="";
	if ($peak_tool eq "all") {
		$peak_calling_macs2=&macs2_diff_peak($strandness);
		$peak_calling_exomePeak=&three_peak_suite_diff_peak("exomePeak",$strandness);
		$peak_calling_MeTDiff=&three_peak_suite_diff_peak("MeTDiff",$strandness);
		$bash_link=$bam_setting.$peak_calling_3peakSuite_divbam.$peak_calling_3peakSuite_divbam_samtools.$peak_calling_macs2.$peak_calling_exomePeak.$peak_calling_MeTDiff.$peak_calling_3peakSuite_divbam_rm;
		return $bash_link;
	} else {
		my @tools=split /,/,$peak_tool;
		my %tool_count=();
		@tools=grep {++$tool_count{$_}<2} @tools;
		my $peak_calling_3peakSuite_link="";
		for (my $i = 0; $i < @tools; $i++) {
			if ($tools[$i] eq "exomePeak") {
				$peak_calling_exomePeak=&three_peak_suite_diff_peak("exomePeak",$strandness);
			} elsif ($tools[$i] eq "MeTDiff") {
				$peak_calling_MeTDiff=&three_peak_suite_diff_peak("MeTDiff",$strandness);
			} elsif ($tools[$i] eq "MACS2") {
				$peak_calling_macs2=&macs2_diff_peak($strandness);
			}
		}
		unless ($peak_calling_exomePeak eq "" && $peak_calling_MeTDiff eq "" && $peak_calling_macs2 eq "" ) {
			$bash_link=$bam_setting.$peak_calling_3peakSuite_divbam.$peak_calling_3peakSuite_divbam_samtools.$peak_calling_macs2.$peak_calling_exomePeak.$peak_calling_MeTDiff.$peak_calling_3peakSuite_divbam_rm;
			return $bash_link;
		} else {
			die "unrecognized tools!\n";
		}
	}
}

sub macs2_diff_peak {
	my $strandness=shift;
	my $macs2_pipeline;
	if ($strandness eq "R" || $strandness eq "F" || $strandness eq "RF" || $strandness eq "FR") {
		$macs2_pipeline="
echo -e \"\\npeak calling -- MACS2\\n\"
echo -e \"\\n\$compare -- control plus strand\\n\"
if [[ ! -d \$peak_out_dir/\$compare/macs2/plus ]];then mkdir -p \$peak_out_dir/\$compare/macs2/plus;fi
macs2 callpeak -t \$hisat2_out_bam_control_ip_plus -c \$hisat2_out_bam_control_input_plus -f BAM -g \$tx_size -n \$peak_out_dir/\$compare/macs2/plus/control \\
	-B --SPMR --nomodel --tsize \$control_read_length --extsize \$control_fragment_length --keep-dup all > \$peak_out_dir/\$compare/macs2/run.plus.log 2>&1
tags=`egrep \"total tags in treatment|total tags in control\" \$peak_out_dir/\$compare/macs2/plus/control_peaks.xls`
echo \$tags >> \$peak_out_dir/\$compare/macs2/run.plus.log
d1=`echo \$tags | awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {split(\$0,a,/(: )|( #)/);t=a[2];c=a[4]} END {if(t>c) {print c} else {print t}}'`
echo -e \"\\n\$compare -- treatment plus strand\\n\"
macs2 callpeak -t \$hisat2_out_bam_treatment_ip_plus -c \$hisat2_out_bam_treatment_input_plus -f BAM -g \$tx_size -n \$peak_out_dir/\$compare/macs2/plus/treatment \\
	-B --SPMR --nomodel --tsize \$treatment_read_length --extsize \$treatment_fragment_length --keep-dup all >> \$peak_out_dir/\$compare/macs2/run.plus.log 2>&1
tags=`egrep \"total tags in treatment|total tags in control\" \$peak_out_dir/\$compare/macs2/plus/treatment_peaks.xls`
echo \$tags >> \$peak_out_dir/\$compare/macs2/run.plus.log
d2=`echo \$tags | awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {split(\$0,a,/(: )|( #)/);t=a[2];c=a[4]} END {if(t>c) {print c} else {print t}}'`
macs2 bdgdiff --t1 \$peak_out_dir/\$compare/macs2/plus/control_treat_pileup.bdg --c1 \$peak_out_dir/\$compare/macs2/plus/control_control_lambda.bdg \\
	--t2 \$peak_out_dir/\$compare/macs2/plus/treatment_treat_pileup.bdg --c2 \$peak_out_dir/\$compare/macs2/plus/treatment_control_lambda.bdg \\
	--d1 \$d1 --d2 \$d2 -g \$read_length -l \$fragment_length --o-prefix \$peak_out_dir/\$compare/macs2/plus/\$compare >> \$peak_out_dir/\$compare/macs2/run.plus.log 2>&1

echo -e \"\\n\$compare -- control minus strand\\n\"
if [[ ! -d \$peak_out_dir/\$compare/macs2/minus ]];then mkdir -p \$peak_out_dir/\$compare/macs2/minus;fi
macs2 callpeak -t \$hisat2_out_bam_control_ip_minus -c \$hisat2_out_bam_control_input_minus -f BAM -g \$tx_size -n \$peak_out_dir/\$compare/macs2/minus/control \\
	-B --SPMR --nomodel --tsize \$control_read_length --extsize \$control_fragment_length --keep-dup all > \$peak_out_dir/\$compare/macs2/run.minus.log 2>&1
tags=`egrep \"total tags in treatment|total tags in control\" \$peak_out_dir/\$compare/macs2/minus/control_peaks.xls`
echo \$tags >> \$peak_out_dir/\$compare/macs2/run.minus.log
d1=`echo \$tags | awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {split(\$0,a,/(: )|( #)/);t=a[2];c=a[4]} END {if(t>c) {print c} else {print t}}'`
echo -e \"\\n\$compare -- treatment minus strand\\n\"
macs2 callpeak -t \$hisat2_out_bam_treatment_ip_minus -c \$hisat2_out_bam_treatment_input_minus -f BAM -g \$tx_size -n \$peak_out_dir/\$compare/macs2/minus/treatment \\
	-B --SPMR --nomodel --tsize \$treatment_read_length --extsize \$treatment_fragment_length --keep-dup all >> \$peak_out_dir/\$compare/macs2/run.minus.log 2>&1
tags=`egrep \"total tags in treatment|total tags in control\" \$peak_out_dir/\$compare/macs2/minus/treatment_peaks.xls`
echo \$tags >> \$peak_out_dir/\$compare/macs2/run.minus.log
d2=`echo \$tags | awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {split(\$0,a,/(: )|( #)/);t=a[2];c=a[4]} END {if(t>c) {print c} else {print t}}'`
macs2 bdgdiff --t1 \$peak_out_dir/\$compare/macs2/minus/control_treat_pileup.bdg --c1 \$peak_out_dir/\$compare/macs2/minus/control_control_lambda.bdg \\
	--t2 \$peak_out_dir/\$compare/macs2/minus/treatment_treat_pileup.bdg --c2 \$peak_out_dir/\$compare/macs2/minus/treatment_control_lambda.bdg \\
	--d1 \$d1 --d2 \$d2 -g \$read_length -l \$fragment_length --o-prefix \$peak_out_dir/\$compare/macs2/minus/\$compare >> \$peak_out_dir/\$compare/macs2/run.minus.log 2>&1";
		if ($strandness eq "F" || $strandness eq "FR" || $strandness eq "RF") {
			$macs2_pipeline.="
awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {if(ARGIND==1) {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_plus_\"a[length(a)],\$5,\"+\"} else {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_minus_\"a[length(a)],\$5,\"-\"}}' \\
	\$peak_out_dir/\$compare/macs2/plus/control_summits.bed \$peak_out_dir/\$compare/macs2/minus/control_summits.bed | sort -k 1,1V -k 2,2n > \$peak_out_dir/\$compare/macs2/control_summits.bed
awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {if(ARGIND==1) {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_plus_\"a[length(a)],\$5,\"+\",\$7,\$8,\$9,\$10} else {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_minus_\"a[length(a)],\$5,\"-\",\$7,\$8,\$9,\$10}}' \\
	\$peak_out_dir/\$compare/macs2/plus/control_peaks.narrowPeak \$peak_out_dir/\$compare/macs2/minus/control_peaks.narrowPeak | sort -k 1,1V -k 2,2n > \$peak_out_dir/\$compare/macs2/control_peaks.narrowPeak
awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {if(ARGIND==1) {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_plus_\"a[length(a)],\$5,\"+\"} else {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_minus_\"a[length(a)],\$5,\"-\"}}' \\
	\$peak_out_dir/\$compare/macs2/plus/treatment_summits.bed \$peak_out_dir/\$compare/macs2/minus/treatment_summits.bed | sort -k 1,1V -k 2,2n > \$peak_out_dir/\$compare/macs2/treatment_summits.bed
awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {if(ARGIND==1) {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_plus_\"a[length(a)],\$5,\"+\",\$7,\$8,\$9,\$10} else {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_minus_\"a[length(a)],\$5,\"-\",\$7,\$8,\$9,\$10}}' \\
	\$peak_out_dir/\$compare/macs2/plus/treatment_peaks.narrowPeak \$peak_out_dir/\$compare/macs2/minus/treatment_peaks.narrowPeak | sort -k 1,1V -k 2,2n > \$peak_out_dir/\$compare/macs2/treatment_peaks.narrowPeak
cat \$peak_out_dir/\$compare/macs2/plus/\${compare}_c3.0_cond[12].bed | awk -F '\\t' 'BEGIN {OFS=\"\\t\"} \$0 !~ /^track/ {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_plus_\"a[length(a)],\$5,\"+\"}' > \\
	\$peak_out_dir/\$compare/macs2/plus/\${compare}_c3.0_cond12.bed
cat \$peak_out_dir/\$compare/macs2/minus/\${compare}_c3.0_cond[12].bed | awk -F '\\t' 'BEGIN {OFS=\"\\t\"} \$0 !~ /^track/ {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_minus_\"a[length(a)],\$5,\"-\"}' > \\
	\$peak_out_dir/\$compare/macs2/minus/\${compare}_c3.0_cond12.bed";
		} else {
			$macs2_pipeline.="
awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {if(ARGIND==1) {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_plus_\"a[length(a)],\$5,\"-\"} else {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_minus_\"a[length(a)],\$5,\"+\"}}' \\
	\$peak_out_dir/\$compare/macs2/plus/control_summits.bed \$peak_out_dir/\$compare/macs2/minus/control_summits.bed | sort -k 1,1V -k 2,2n > \$peak_out_dir/\$compare/macs2/control_summits.bed
awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {if(ARGIND==1) {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_plus_\"a[length(a)],\$5,\"-\",\$7,\$8,\$9,\$10} else {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_minus_\"a[length(a)],\$5,\"+\",\$7,\$8,\$9,\$10}}' \\
	\$peak_out_dir/\$compare/macs2/plus/control_peaks.narrowPeak \$peak_out_dir/\$compare/macs2/minus/control_peaks.narrowPeak | sort -k 1,1V -k 2,2n > \$peak_out_dir/\$compare/macs2/control_peaks.narrowPeak
awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {if(ARGIND==1) {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_plus_\"a[length(a)],\$5,\"-\"} else {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_minus_\"a[length(a)],\$5,\"+\"}}' \\
	\$peak_out_dir/\$compare/macs2/plus/treatment_summits.bed \$peak_out_dir/\$compare/macs2/minus/treatment_summits.bed | sort -k 1,1V -k 2,2n > \$peak_out_dir/\$compare/macs2/treatment_summits.bed
awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {if(ARGIND==1) {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_plus_\"a[length(a)],\$5,\"-\",\$7,\$8,\$9,\$10} else {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_minus_\"a[length(a)],\$5,\"+\",\$7,\$8,\$9,\$10}}' \\
	\$peak_out_dir/\$compare/macs2/plus/treatment_peaks.narrowPeak \$peak_out_dir/\$compare/macs2/minus/treatment_peaks.narrowPeak | sort -k 1,1V -k 2,2n > \$peak_out_dir/\$compare/macs2/treatment_peaks.narrowPeak
cat \$peak_out_dir/\$compare/macs2/plus/\${compare}_c3.0_cond[12].bed | awk -F '\\t' 'BEGIN {OFS=\"\\t\"} \$0 !~ /^track/ {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_plus_\"a[length(a)],\$5,\"-\"}' > \\
	\$peak_out_dir/\$compare/macs2/plus/\${compare}_c3.0_cond12.bed
cat \$peak_out_dir/\$compare/macs2/minus/\${compare}_c3.0_cond[12].bed | awk -F '\\t' 'BEGIN {OFS=\"\\t\"} \$0 !~ /^track/ {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_minus_\"a[length(a)],\$5,\"+\"}' > \\
	\$peak_out_dir/\$compare/macs2/minus/\${compare}_c3.0_cond12.bed";
		}
		$macs2_pipeline.="
cat \$peak_out_dir/\$compare/macs2/run.plus.log \$peak_out_dir/\$compare/macs2/run.minus.log > \$peak_out_dir/\$compare/macs2/run.log
cat \$peak_out_dir/\$compare/macs2/plus/\${compare}_c3.0_common.bed \$peak_out_dir/\$compare/macs2/minus/\${compare}_c3.0_common.bed | sort -k 1,1V -k 2,2n > \$peak_out_dir/\$compare/macs2/\${compare}_c3.0_common.bed
cat \$peak_out_dir/\$compare/macs2/plus/\${compare}_c3.0_cond12.bed \$peak_out_dir/\$compare/macs2/minus/\${compare}_c3.0_cond12.bed | sort -k 1,1V -k 2,2n > \$peak_out_dir/\$compare/macs2/\${compare}_c3.0_cond12.bed
awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {s=int(\$2+(\$3-\$2)/2);print \$1,s,s+1,\$4,\$5,\$6}' \$peak_out_dir/\$compare/macs2/\${compare}_c3.0_cond12.bed > \$peak_out_dir/\$compare/macs2/\${compare}_c3.0_cond12_midpoint.bed";
	} elsif ($strandness eq "U") {
		$macs2_pipeline="
echo -e \"\\npeak calling -- MACS2\\n\"
echo -e \"\\n\$compare -- control\\n\"
if [[ ! -d \$peak_out_dir/\$compare/macs2 ]];then mkdir -p \$peak_out_dir/\$compare/macs2;fi
macs2 callpeak -t \$hisat2_out_bam_control_ip -c \$hisat2_out_bam_control_input -f BAM -g \$tx_size -n \$peak_out_dir/\$compare/macs2/control \\
	-B --SPMR --nomodel --tsize \$control_read_length --extsize \$control_fragment_length --keep-dup all > \$peak_out_dir/\$compare/macs2/run.log 2>&1
awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_\"a[length(a)],\$5,\".\"}' \$peak_out_dir/\$compare/macs2/control_summits.bed > \$peak_out_dir/\$compare/macs2/control_summits.bed.tmp
mv \$peak_out_dir/\$compare/macs2/control_summits.bed.tmp \$peak_out_dir/\$compare/macs2/control_summits.bed
awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_\"a[length(a)],\$5,\$6,\$7,\$8,\$9,\$10}' \$peak_out_dir/\$compare/macs2/control_peaks.narrowPeak > \$peak_out_dir/\$compare/macs2/control_peaks.narrowPeak.tmp
mv \$peak_out_dir/\$compare/macs2/control_peaks.narrowPeak.tmp \$peak_out_dir/\$compare/macs2/control_peaks.narrowPeak
tags=`egrep \"total tags in treatment|total tags in control\" \$peak_out_dir/\$compare/macs2/control_peaks.xls`
echo \$tags >> \$peak_out_dir/\$compare/macs2/run.log
d1=`echo \$tags | awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {split(\$0,a,/(: )|( #)/);t=a[2];c=a[4]} END {if(t>c) {print c} else {print t}}'`

echo -e \"\\n\$compare -- treatment\\n\"
macs2 callpeak -t \$hisat2_out_bam_treatment_ip -c \$hisat2_out_bam_treatment_input -f BAM -g \$tx_size -n \$peak_out_dir/\$compare/macs2/treatment \\
	-B --SPMR --nomodel --tsize \$treatment_read_length --extsize \$treatment_fragment_length --keep-dup all >> \$peak_out_dir/\$compare/macs2/run.log 2>&1
awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_\"a[length(a)],\$5,\".\"}' \$peak_out_dir/\$compare/macs2/treatment_summits.bed > \$peak_out_dir/\$compare/macs2/treatment_summits.bed.tmp
mv \$peak_out_dir/\$compare/macs2/treatment_summits.bed.tmp \$peak_out_dir/\$compare/macs2/treatment_summits.bed
awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_\"a[length(a)],\$5,\$6,\$7,\$8,\$9,\$10}' \$peak_out_dir/\$compare/macs2/treatment_peaks.narrowPeak > \$peak_out_dir/\$compare/macs2/treatment_peaks.narrowPeak.tmp
mv \$peak_out_dir/\$compare/macs2/treatment_peaks.narrowPeak.tmp \$peak_out_dir/\$compare/macs2/treatment_peaks.narrowPeak
tags=`egrep \"total tags in treatment|total tags in control\" \$peak_out_dir/\$compare/macs2/treatment_peaks.xls`
echo \$tags >> \$peak_out_dir/\$compare/macs2/run.log
d2=`echo \$tags | awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {split(\$0,a,/(: )|( #)/);t=a[2];c=a[4]} END {if(t>c) {print c} else {print t}}'`
macs2 bdgdiff --t1 \$peak_out_dir/\$compare/macs2/control_treat_pileup.bdg --c1 \$peak_out_dir/\$compare/macs2/control_control_lambda.bdg \\
	--t2 \$peak_out_dir/\$compare/macs2/treatment_treat_pileup.bdg --c2 \$peak_out_dir/\$compare/macs2/treatment_control_lambda.bdg \\
	--d1 \$d1 --d2 \$d2 -g \$read_length -l \$fragment_length --o-prefix \$peak_out_dir/\$compare/macs2/\$compare >> \$peak_out_dir/\$compare/macs2/run.log 2>&1
cat \$peak_out_dir/\$compare/macs2/\${compare}_c3.0_cond[12].bed | awk -F '\\t' 'BEGIN {OFS=\"\\t\"} \$0 !~ /^track/ {split(\$4,a,\"_\");print \$1,\$2,\$3,a[length(a)-1]\"_\"a[length(a)],\$5,\".\"}' > \\
	\$peak_out_dir/\$compare/macs2/\${compare}_c3.0_cond12_unannot.bed
intersectBed -wo -f 0.5 -a \$peak_out_dir/\$compare/macs2/\${compare}_c3.0_cond12_unannot.bed -b \$bed12_file | awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {print \$1,\$2,\$3,\$4,\$5,\$12}' | uniq | sort -k 1,1V -k 2,2n > \\
	\$peak_out_dir/\$compare/macs2/\${compare}_c3.0_cond12.bed
awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {s=int(\$2+(\$3-\$2)/2);print \$1,s,s+1,$4,$5,$6}' \$peak_out_dir/\$compare/macs2/\${compare}_c3.0_cond12.bed > \$peak_out_dir/\$compare/macs2/\${compare}_c3.0_cond12_midpoint.bed
";
	} else {
		die "exp $runName: strandness must be F, R, FR, RF or U!\n";
	}
	$macs2_pipeline.="
common_peak=`cat \$peak_out_dir/\$compare/macs2/\${compare}_c3.0_common.bed | wc -l`
DE_peak1=`grep cond1 \$peak_out_dir/\$compare/macs2/\${compare}_c3.0_cond12.bed | wc -l`
DE_peak2=`grep cond2 \$peak_out_dir/\$compare/macs2/\${compare}_c3.0_cond12.bed | wc -l`
echo -e \"\\n`expr \$common_peak - 1` common peaks were found\" >> \$peak_out_dir/\$compare/macs2/run.log
echo -e \"\\n\$DE_peak1 cond1 peaks were found.\\n\" >> \$peak_out_dir/\$compare/macs2/run.log
echo -e \"\\n\$DE_peak2 cond2 peaks were found.\\n\" >> \$peak_out_dir/\$compare/macs2/run.log
slopBed -b 25 -g \$chrom_size -i \$peak_out_dir/\$compare/macs2/\${compare}_c3.0_cond12_midpoint.bed > \$peak_out_dir/\$compare/macs2/\${compare}_c3.0_cond12_mid50.bed
fastaFromBed -name+ -split -s -fi \$genome_fa -bed \$peak_out_dir/\$compare/macs2/\${compare}_c3.0_cond12_mid50.bed > \$peak_out_dir/\$compare/macs2/\${compare}_c3.0_cond12_mid50.fa
shuffleBed -incl \$bed12_file -seed 12345 -noOverlapping -i \$peak_out_dir/\$compare/macs2/\${compare}_c3.0_cond12_mid50.bed -g \$chrom_size > \$peak_out_dir/\$compare/macs2/\${compare}_c3.0_cond12_mid50_random.bed
fastaFromBed -name+ -split -s -fi \$genome_fa -bed \$peak_out_dir/\$compare/macs2/\${compare}_c3.0_cond12_mid50_random.bed > \$peak_out_dir/\$compare/macs2/\${compare}_c3.0_cond12_mid50_random.fa
findMotifs.pl \$peak_out_dir/\$compare/macs2/\${compare}_c3.0_cond12_mid50.fa fasta \$peak_out_dir/\$compare/macs2/homer -fasta \$peak_out_dir/\$compare/macs2/\${compare}_c3.0_cond12_mid50_random.fa \\
	-p \$ncpus -len 5,6,7,8 -S 10 -rna -dumpFasta > \$peak_out_dir/\$compare/macs2/homer.run.log 2>&1
echo -e \"\\nmacs2 done!\\n\"
";
}

sub three_peak_suite_diff_peak {
	my $tool=shift;
	my $strandness=shift;
	my $peak_calling_3peakSuite;
	if ($strandness eq "R" || $strandness eq "F" || $strandness eq "RF" || $strandness eq "FR") {
		$peak_calling_3peakSuite="
echo -e \"\\npeak calling -- $tool\\n\"
echo -e \"\\n\$compare\\n\"
if [[ ! -d \$peak_out_dir/\$compare/$tool/plus ]];then mkdir -p \$peak_out_dir/\$compare/$tool/plus;fi
if [[ ! -d \$peak_out_dir/\$compare/$tool/minus ]];then mkdir -p \$peak_out_dir/\$compare/$tool/minus;fi
3peakSuite.R $tool \$gtf_file \$peak_out_dir/\$compare/$tool plus \$fragment_length \$read_length \$hisat2_out_bam_control_input_plus \$hisat2_out_bam_control_ip_plus \\
	\$hisat2_out_bam_treatment_input_plus \$hisat2_out_bam_treatment_ip_plus > \$peak_out_dir/\$compare/$tool/run.plus.log 2>&1 &
3peakSuite.R $tool \$gtf_file \$peak_out_dir/\$compare/$tool minus \$fragment_length \$read_length \$hisat2_out_bam_control_input_minus \$hisat2_out_bam_control_ip_minus \\
	\$hisat2_out_bam_treatment_input_minus \$hisat2_out_bam_treatment_ip_minus > \$peak_out_dir/\$compare/$tool/run.minus.log 2>&1
wait
";
		if ($tool eq "exomePeak") {
			$peak_calling_3peakSuite.="resfiles=(con_sig_diff_peak.bed con_sig_diff_peak.xls diff_peak.bed diff_peak.xls sig_diff_peak.bed sig_diff_peak.xls)"
		} else {
			$peak_calling_3peakSuite.="resfiles=(peak.bed peak.xls diff_peak.bed diff_peak.xls)"
		}
		if ($strandness eq "F" || $strandness eq "FR" || $strandness eq "RF") {
			$peak_calling_3peakSuite.="
for i in \${resfiles[@]}; do
	awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {if(ARGIND==1) {if(\$6 == \"+\") print \$0} else {if(\$6 == \"-\") print \$0}}' \$peak_out_dir/\$compare/$tool/plus/\$i \\
		\$peak_out_dir/\$compare/$tool/minus/\$i > \$peak_out_dir/\$compare/$tool/\$i
done";
		} else {
			$peak_calling_3peakSuite.="
for i in \${resfiles[@]}; do
	awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {if(ARGIND==1) {if(\$6 == \"-\") print \$0} else {if(\$6 == \"+\") print \$0}}' \$peak_out_dir/\$compare/$tool/plus/\$i \\
		\$peak_out_dir/\$compare/$tool/minus/\$i > \$peak_out_dir/\$compare/$tool/\$i
done";
		}
	} elsif ($strandness eq "U") {
		$peak_calling_3peakSuite="
echo -e \"\\npeak calling -- $tool\\n\"
echo -e \"\\n\$compare\\n\"
if [[ ! -d \$peak_out_dir/\$compare/$tool ]];then mkdir -p \$peak_out_dir/\$compare/$tool;fi
3peakSuite.R $tool \$gtf_file \$peak_out_dir/\$compare $tool \$fragment_length \$read_length \$hisat2_out_bam_control_input \$hisat2_out_bam_control_ip \\
	\$hisat2_out_bam_treatment_input \$hisat2_out_bam_treatment_ip > \$peak_out_dir/\$compare/$tool/run.log 2>&1
";
	} else {
		die "exp $runName: strandness must be F, R, FR, RF or U!\n";
	}
	$peak_calling_3peakSuite.="
fastaFromBed -name+ -split -s -fi \$genome_fa -bed \$peak_out_dir/\$compare/$tool/diff_peak.bed > \$peak_out_dir/\$compare/$tool/diff_peak.fa
shuffleBed -incl \$bed12_file -seed 12345 -noOverlapping -i \$peak_out_dir/\$compare/$tool/diff_peak.bed -g \$chrom_size > \$peak_out_dir/\$compare/$tool/random_diff_peak.bed
fastaFromBed -name+ -split -s -fi \$genome_fa -bed \$peak_out_dir/\$compare/$tool/random_diff_peak.bed > \$peak_out_dir/\$compare/$tool/random_diff_peak.fa
findMotifs.pl \$peak_out_dir/\$compare/$tool/diff_peak.fa fasta \$peak_out_dir/\$compare/$tool/homer -fasta \$peak_out_dir/\$compare/$tool/random_diff_peak.fa \\
	-p \$ncpus -len 5,6,7,8 -S 10 -rna -dumpFasta > \$peak_out_dir/\$compare/$tool/homer_run.log 2>&1
echo -e \"\\n$tool done!\\n\"
";
}


__END__


=head1 NAME

easym6A - Creat a straightword workflow in a bash script. 

=head1 SYNOPSIS

single job mode (call peaks):

$ easym6A.pl -s <sampleList.txt> -c <configure.txt> -n <ID,ID> -e <run_name> [options]

single job mode (call differential peaks):

$ easym6A.pl -s <sampleList.txt> -c <configure.txt> -n <ID,ID,ID,ID> -e <run_name> [options]

batch job mode:

$ easym6A.pl -s <sampleList.txt> -b <batchList.txt> -x <num> -y <num> [options]

 Options:
   -h/--help                   brief help message.
   --man                       full documentation.
   -s/--samplelist <file>      a table file that records sample infomation. (Required)
   -c/--configure <file>       a configuration file that records paths of required files and output. (Requried)
   -n/--sampleno <string>      a comma-seperated ID list of samples. (Required)
   -e/--runname <string>       a user-defined job name. (Required)
   -t/--threads <num>          number of CPU threads to run the pipeline. (Default: 1)
   -l/--parallel               reduce the analysis time in parallel mode. (Default: off)
   -m/--method <string>        name(s) of peak calling tool(s) that are used. (Default: all)
   -b/--batch <file>           a table file that records batch job information.
   -x/--bstart <int>           batch job start ID. Required when -b/--batch is set. Used together with -y/--bend.
   -y/--bend <int>             batch job end ID. Required when -b/--batch is set. Used together with -x/--bstart.
   -a/--onlybam                run gene expression analysis only. (Default: off)
   -u/--daq                    calculate gene expression in both of samples. (Default: off)
   -k/--onlypeak               run peak calling analysis only. (Default: off)
   -p/--rmrep                  remove reads from repetitive elements. (Default: off)
   -d/--rmdup                  remove reads from PCR duplication. (Default: off)
   -f/--keeptmp                keep intermediated files. (Default: off)
   -r/--run                    run bash script(s) generated by the pipeline. (Default: off)

=head1 OPTIONS

=over 4

=item B<-h/--help>

Print a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

=item B<-s/--samplelist> 

Required. The path of a table file that includes m6A/MeRIP-seq sample information. Created in Step 2.

=item B<-c/--configure> 

Required when the single job mode is used. The path of a configuration file that includes paths of required files and output. Created in Step 3.

=item B<-n/--sampleno> 

Required when the single job mode is used. A comma-seperated ID list of samples. Sample IDs are indicated in the table file of Step 2. Note that the 
current version of easym6A can not analyze samples with replicates. Therefore, the option is always `ID,ID` (input,IP) for performing peak detection or `ID,ID,ID,ID` (contorl input,control IP,treatment input,treatment IP) for finding differential peaks. For instance, `1,2` represents that easym6A analyzes samples whose ID are 1 (input sample) and 2 (IP sample); `1,2,3,4` represents that easym6A analyzes samples whose ID are 1 (control input sample), 2 (control IP sample), 3 (treatment input sample) and 4 (treatment IP sample).

=item B<-e/--runname> 

Required when the single job mode is used. Used for naming the bash script and log file..

=item B<-t/--threads> 

Optional. The number of CPU threads is used to run the pipeline. Default: `1`.

=item B<-l/--parallel> 

Optional. The current version of easym6A does not support the option.

=item B<-m/--method> 

Optional. Accepted values are `exomePeak`, `MeTPeak`, 'MeTDiff' `MACS2` or `all`. A comma-seperated list is available to tell easym6A which tools are 
used together.  Note that `exomePeak`, `MeTPeak` and `MACS2`are used in peak detection. `exomePeak` and `MetDiff` are used in finding differential peaks. Default: `all`.

=item B<-b/--batch> 

Required when the batch job mode is used. The path of a table file that includes batch job information. See below.

=item B<-x/--bstart> 

Required when the batch job mode is used. The starting ID of a batch job. Used together with the option `-y/--bend`.

=item B<-y/--bend> 

Required when the batch job mode is used. The ending ID of a batch job. Used together with the option `-x/--bstart`.

=item B<-a/--onlybam> 

Optional. If it is specified, easym6A only generates alignment and gene quantification results, exactly as what RNA-seq data was processed. When 
easym6A was used to process RNA-seq data, `-a/--onlybam` should be specified together with `-u/--daq`. Default: off.

=item B<-u/--daq> 

Optional. Originally, easym6A is designed to calculate gene quantification in the input sample but not in the IP sample of m6A/MeRIP-seq data. When `
-u/--daq` is specified, gene expression is calculated in both of samples. It is useful when easym6A is regarded as a RNA-seq data processing pipeline. Default: off.

=item B<-k/--onlypeak> 

Optional. If it is specified, easym6A only performs the step of peak calling by using exisiting bam files whose paths are in the sample information 
file. Default: off.

=item B<-p/--rmrep> 

Optional. Remove reads mapped to repetitive elements first before they are mapped to genomes. When `-p/--rmrep` is specified, a HISAT2 index of 
repetitive elements is required. Default: off.

=item B<-d/--rmdup> 

Optional. Remove reads from PCR duplication detected by Picard. Default: off.

=item B<-f/--keeptmp> 

Optional. Keep intermediated files generated during data processing. Default: off.

=item B<-r/--run> 

Run bash scripts generated by easym6A. Default: off.

=back

=head1 DESCRIPTION

B<easym6A>: Process m6A/MeRIP-seq data in single or batch job mode.

=cut
