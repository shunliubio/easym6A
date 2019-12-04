#!/usr/bin/env /software/R-3.3-el7-x86_64/bin/Rscript

start <- Sys.time()
cat("Starting Time is ",format(start,format="%Y-%m-%d %H:%M:%S"),"\n\n",sep="")

argv <- commandArgs(TRUE)

.libPaths("/software/R-3.3-el7-x86_64/lib64/R/library")
options(scipen=5)

tool <- argv[1]
gtf <- argv[2]
out_dir <- argv[3]
exp_name <- argv[4]
frag_length <- as.numeric(argv[5])
read_length <- as.numeric(argv[6])
input_bam <- argv[7]
ip_bam <- argv[8]
peak_cutoff_pvalue <- 1e-3
peak_cutoff_fdr <- 0.05
window_width <- 50
sliding_step <- 10
minimal_mapq <- 20
fold_enrichment <- 2
diff_peak_cutoff_fdr <- 0.1
diff_peak_abs_fold_change <- 1.5
diff_peak_consistent_abs_fold_change <- 1.5
dup <- F

work_dir <- paste(out_dir,exp_name,sep="/")
setwd(work_dir)
cat("work directory changed to ",work_dir,"\n\n",sep="")

cat("tool is ",tool,"\n",sep="")
cat("gtf is ",gtf,"\n",sep="")
cat("out_dir is ",out_dir,"\n",sep="")
cat("exp_name is ",exp_name,"\n",sep="")
cat("frag_length is ",frag_length,"\n",sep="")
cat("read_length is ",read_length,"\n",sep="")
cat("input_bam is ",input_bam,"\n",sep="")
cat("ip_bam is ",ip_bam,"\n",sep="")
cat("peak_cutoff_pvalue is ",peak_cutoff_pvalue,"\n",sep="")
cat("peak_cutoff_fdr is ",peak_cutoff_fdr,"\n",sep="")
cat("window_width is ",window_width,"\n",sep="")
cat("sliding_step is ",sliding_step,"\n",sep="")
cat("minimal_mapq is ",minimal_mapq,"\n",sep="")
cat("fold_enrichment is ",fold_enrichment,"\n",sep="")
cat("diff_peak_cutoff_fdr is ",diff_peak_cutoff_fdr,"\n",sep="")
cat("diff_peak_abs_fold_change is ",diff_peak_abs_fold_change,"\n",sep="")
cat("diff_peak_consistent_abs_fold_change is ",diff_peak_consistent_abs_fold_change,"\n",sep="")
cat("deduplication is ",dup,"\n\n",sep="")

input_rep <- regexpr(",",input_bam)[1]
ip_rep <- regexpr(",",ip_bam)[1]
if (input_rep != -1) {
	input_bam <- unlist(strsplit(input_bam,","))
}
if (ip_rep != -1) {
	ip_bam <- unlist(strsplit(ip_bam,","))
}

if (length(argv) == 8) {
	if (tool == "exomePeak") {
		library(exomePeak)
		exomepeak_result=exomepeak(GENE_ANNO_GTF=gtf,OUTPUT_DIR=out_dir,EXPERIMENT_NAME=exp_name,FRAGMENT_LENGTH=frag_length,READ_LENGTH=read_length,IP_BAM=ip_bam,INPUT_BAM=input_bam,
			PEAK_CUTOFF_FDR=peak_cutoff_fdr,WINDOW_WIDTH=window_width,SLIDING_STEP=sliding_step,MINIMAL_MAPQ=minimal_mapq,FOLD_ENRICHMENT=fold_enrichment,REMOVE_LOCAL_TAG_ANOMALITIES=dup)
	} else if (tool == "MeTPeak") {
		library(MeTPeak)
		metpeak_result=metpeak(GENE_ANNO_GTF=gtf,OUTPUT_DIR=out_dir,EXPERIMENT_NAME=exp_name,FRAGMENT_LENGTH=frag_length,READ_LENGTH=read_length,IP_BAM=ip_bam,INPUT_BAM=input_bam,
			PEAK_CUTOFF_FDR=peak_cutoff_fdr,WINDOW_WIDTH=window_width,SLIDING_STEP=sliding_step,MINIMAL_MAPQ=minimal_mapq,FOLD_ENRICHMENT=fold_enrichment,REMOVE_LOCAL_TAG_ANOMALITIES=dup)
	} else {
		cat("Either exomePeak or MeTPeak is used in peak calling mode!!!\n")
	}
} else if (length(argv) == 10) {
	treated_input_bam <- argv[9]
	treated_ip_bam <- argv[10]
	cat("treated_input_bam is ",treated_input_bam,"\n",sep="")
	cat("treated_ip_bam is ",treated_ip_bam,"\n\n",sep="")
	treated_input_rep <- regexpr(",",treated_input_bam)[1]
	treated_ip_rep <- regexpr(",",treated_ip_bam)[1]
	if (treated_input_rep != -1) {
		treated_input_bam <- unlist(strsplit(treated_input_bam,","))
	}
	if (treated_ip_rep != -1) {
		treated_ip_bam <- unlist(strsplit(treated_ip_bam,","))
	}
	if (tool == "exomePeak") {
		library(exomePeak)
		exomepeak_result=exomepeak(GENE_ANNO_GTF=gtf,OUTPUT_DIR=out_dir,EXPERIMENT_NAME=exp_name,FRAGMENT_LENGTH=frag_length,READ_LENGTH=read_length,IP_BAM=ip_bam,INPUT_BAM=input_bam,
			TREATED_IP_BAM=treated_ip_bam,TREATED_INPUT_BAM=treated_input_bam,PEAK_CUTOFF_PVALUE=peak_cutoff_pvalue,WINDOW_WIDTH=window_width,SLIDING_STEP=sliding_step,
			MINIMAL_MAPQ=minimal_mapq,FOLD_ENRICHMENT=fold_enrichment,DIFF_PEAK_CUTOFF_FDR=diff_peak_cutoff_fdr,DIFF_PEAK_ABS_FOLD_CHANGE=diff_peak_abs_fold_change,
			DIFF_PEAK_CONSISTENT_ABS_FOLD_CHANGE=diff_peak_consistent_abs_fold_change,REMOVE_LOCAL_TAG_ANOMALITIES=dup)
	} else if (tool == "MeTDiff") {
		library(MeTDiff)
		metdiff_result=metdiff(GENE_ANNO_GTF=gtf,OUTPUT_DIR=out_dir,EXPERIMENT_NAME=exp_name,FRAGMENT_LENGTH=frag_length,READ_LENGTH=read_length,IP_BAM=ip_bam,INPUT_BAM=input_bam,
			TREATED_IP_BAM=treated_ip_bam,TREATED_INPUT_BAM=treated_input_bam,PEAK_CUTOFF_PVALUE=peak_cutoff_pvalue,WINDOW_WIDTH=window_width,SLIDING_STEP=sliding_step,
			MINIMAL_MAPQ=minimal_mapq,FOLD_ENRICHMENT=fold_enrichment,DIFF_PEAK_CUTOFF_FDR=diff_peak_cutoff_fdr,DIFF_PEAK_ABS_FOLD_CHANGE=diff_peak_abs_fold_change,
			REMOVE_LOCAL_TAG_ANOMALITIES=dup)
	} else {
		cat("Either exomePeak or MeTDiff is used in differential peak calling mode!!!\n")
	}
} else {
	cat("The number of arguments is either 8 or 10!!!\n")
	q(status = 1)
}
warnings()
save.image()

end <- Sys.time()
cat("Ending Time is",format(end,format="%Y-%m-%d %H:%M:%S"),"\n")
runtime <- difftime(end,start,unit="mins")
cat("Used Time is",runtime,"mins\n")
