#! /usr/bin/env nextflow

// IARC - C. VOEGELE - Last update 22 Apr 2016 //

// Usage : ./GATK_alignment.nf --fastq_folder FASTQ/ --project MY_PROJECT --bed MY_BED --cluster_options MY_OPTIONS 

// Set parameters
params.extension = "_sequence.fq"
params.suffix1 = "_1"
params.suffix2 = "_2"
params.cluster_options = ""

if (params.help) {
    log.info ''
    log.info '--------------------------------------------------'
    log.info 'NEXTFLOW: ALIGNMENT FOLLOWING GATK BEST PRACTICES'
    log.info '--------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow GATK_alignment.nf --fastq_folder FASTQ/ [--extension fq] [--suffix _]'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --fastq_folder   FOLDER                  Folder containing FASTQ paired files'
    log.info '    --project        STRING                  Project name for bam @RG ID'
    log.info '    --bed            STRING                  Bed file (path)'
    log.info 'Options:'
    log.info '    --extension      STRING                  Extension of fastq files (default : fq)'
    log.info '    --suffix1         STRING                 Suffix of fastq files 1 (default : _1)'
    log.info '    --suffix2         STRING                 Suffix of fastq files 2 (default : _2)'
    log.info '    --cluster_options	STRING		   Specific options for cluster scheduler'
    log.info ''
    exit 1
}

// Check folder with fastq files
assert file(params.fastq_folder).listFiles().findAll { it.name ==~ /.*${params.extension}/ }.size() > 0 : "FASTQ folder contains no files"

keys1 = file(params.fastq_folder).listFiles().findAll { it.name ==~ /.*${params.suffix1}${params.extension}/ }.collect { it.getName() }
                                                                                                               .collect { it.replace("${params.suffix1}${params.extension}",'') }
keys2 = file(params.fastq_folder).listFiles().findAll { it.name ==~ /.*${params.suffix2}${params.extension}/ }.collect { it.getName() }
                                                                                                               .collect { it.replace("${params.suffix2}${params.extension}",'') }
if ( !(keys1.containsAll(keys2)) || !(keys2.containsAll(keys1)) ) {println "\n ERROR : There is at least one fastq without its mate, please check your fastq files."; System.exit(0)}

// Gather files ending with _1 suffix
reads1 = Channel
    .fromPath( params.fastq_folder+'/*'+params.suffix1+params.extension )
    .map {  path -> [ path.name.replace("${params.suffix1}${params.extension}",""), path ] }

// Gather files ending with _2 suffix
reads2 = Channel
    .fromPath( params.fastq_folder+'/*'+params.suffix2+params.extension )
    .map {  path -> [ path.name.replace("${params.suffix2}${params.extension}",""), path ] }

// Match the pairs on two channels having the same 'key' (name) and emit a new pair containing the expected files
readPairs = reads1
    .phase(reads2)
    .map { pair1, pair2 -> [ pair1[1], pair2[1] ] }


////////// STEP 01 ################### Alignment => Generates a SAM file containing aligned reads

process creation_sam {
	echo "creation_sam"
	tag { pair_tag }
	cpus 8
	memory '12GB'
	clusterOptions = params.cluster_options
    	//clusterOptions = '-m cn12 -R "rusage[mem=12000]" -M 12000'// PUT NODES AS PARAM	
input:
	file pair from readPairs
output:
	set val(pair_tag), file("${pair_tag}_aligned.sam") into sam_files
shell:
    	pair_tag = pair[0].name.replace("${params.suffix1}${params.extension}","")
    	'''
    	!{params.bwa} mem -M -t 8 -R '@RG\\tID:!{params.project}\\tSM:!{pair[0]}\\tPL:illumina\\tLB:!{pair[0]}\\tPU:unit1' !{params.genome_ref} !{pair[0]} !{pair[1]} > !{pair_tag}_aligned.sam
    	'''
}

////////// STEP 02 ################### Conversion to bam and sort

process creation_bam {
   	echo "creation_bam"
	tag { pair_tag }
	cpus 8
	memory '12GB'
	clusterOptions = params.cluster_options
input:
	set val(pair_tag), file("${pair_tag}_aligned.sam") from sam_files
output:
	set val(pair_tag), file("${pair_tag}_aligned.bam") into bam_files
shell:
   	'''
   	!{params.java17} -Xmx2g -Djava.io.tmpdir=/data/tmp -jar !{params.picard} SortSam INPUT=!{pair_tag}_aligned.sam OUTPUT=!{pair_tag}_aligned.bam TMP_DIR=/data/tmp SORT_ORDER=coordinate
   	'''
}

////////// STEP 03 ################### Mark Duplicates

process creation_md {
	echo "creation_md"
	tag { pair_tag }
	cpus 6
	publishDir 'MD_output',mode:'move',pattern:'*.txt'
	memory '24GB'
	clusterOptions = params.cluster_options
input:
	set val(pair_tag), file("${pair_tag}_aligned.bam") from bam_files
output:
	set val(pair_tag), file("${pair_tag}_aligned_MD.bam") into md_files
shell:
	'''
	!{params.java17} -jar !{params.picard} MarkDuplicates INPUT=!{pair_tag}_aligned.bam OUTPUT=!{pair_tag}_aligned_MD.bam METRICS_FILE=!{pair_tag}_md_metrics.txt
	'''
}

////////// STEP 04 ################### Index MD bam file

process creation_md_bai {
	echo "creation_md_bai"
	tag { pair_tag }
	memory '8GB'
	clusterOptions = params.cluster_options
input:
	set val(pair_tag), file("${pair_tag}_aligned_MD.bam") from md_files
output:
	set val(pair_tag), file("${pair_tag}_aligned_MD.bam"), file("${pair_tag}_aligned_MD.bai") into bam_bai_files, bam_bai_files2
shell:
	'''
	!{params.java17} -jar !{params.picard} BuildBamIndex INPUT=!{pair_tag}_aligned_MD.bam			
	'''
}

////////// STEP 05 ################### Local realignment around indels - create target list to be realigned

process creation_indel_realign_target {
	echo "creation_indel_realign_target"
	tag { pair_tag }
	cpus 6
	memory '16GB'
	clusterOptions = params.cluster_options
input:
	set val(pair_tag), file("${pair_tag}_aligned_MD.bam"), file("${pair_tag}_aligned_MD.bai") from bam_bai_files
output:
	set val(pair_tag), file("${pair_tag}_target_intervals.list") into indel_realign_target_files
shell:
	'''
	!{params.java17} -jar !{params.gatk}GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 6 -R !{params.genome_ref} -I !{pair_tag}_aligned_MD.bam -known !{params.Mills_indels} -known !{params.ThousandG_indels} -o !{pair_tag}_target_intervals.list
	'''
}

////////// STEP 06 ################### Local realignment around indels - perform realignment

process creation_indel_realign_bam {
	echo "creation_indel_realign_bam"
	tag { pair_tag }
	cpus 6
	memory '12GB'
	clusterOptions = params.cluster_options
input:
	set val(pair_tag), file("${pair_tag}_aligned_MD.bam"), file("${pair_tag}_aligned_MD.bai") from bam_bai_files2
	set val(pair_tag), file("${pair_tag}_target_intervals.list") from indel_realign_target_files
output:
	set val(pair_tag), file("${pair_tag}_realigned.bam"), file("${pair_tag}_realigned.bai") into realign_bam_bai_files, realign_bam_bai_files2, realign_bam_bai_files3
shell:
	'''
	!{params.java17} -jar !{params.gatk}GenomeAnalysisTK.jar -T IndelRealigner -R !{params.genome_ref} -I !{pair_tag}_aligned_MD.bam -targetIntervals !{pair_tag}_target_intervals.list -known !{params.Mills_indels} -known !{params.ThousandG_indels} -o !{pair_tag}_realigned.bam			
	'''
}

////////// STEP 07 ################### Base quality score recalibration - Analyze patterns of covariation in the sequence dataset

process creation_recal_table {
	echo "creation_recal_table"
	tag { pair_tag }
	cpus 6
	memory '16GB'
	clusterOptions = params.cluster_options
input:
	set val(pair_tag), file("${pair_tag}_realigned.bam"), file("${pair_tag}_realigned.bai") from realign_bam_bai_files
output:
	set val(pair_tag), file("${pair_tag}_recal.table") into recal_table_files, recal_table_files2, recal_table_files3
shell:
	'''
	!{params.java17} -jar !{params.gatk}GenomeAnalysisTK.jar -T BaseRecalibrator -nct 6 -R !{params.genome_ref} -I !{pair_tag}_realigned.bam -knownSites !{params.dbsnp} -knownSites !{params.Mills_indels} -knownSites !{params.ThousandG_indels} -L !{params.bed} -o !{pair_tag}_recal.table			
	'''
}

////////// STEP 08 ################### Base quality score recalibration - Do a second pass to analyze covariation remaining after recalibration

process creation_recal_table_post {
	echo "creation_recal_table_post"
	tag { pair_tag }
	cpus 6
	memory '16GB'
	clusterOptions = params.cluster_options
input:
    	set val(pair_tag), file("${pair_tag}_realigned.bam"), file("${pair_tag}_realigned.bai") from realign_bam_bai_files2
    	set val(pair_tag), file("${pair_tag}_recal.table") from recal_table_files
output:
    	set val(pair_tag), file("${pair_tag}_post_recal.table") into recal_table_post_files
shell:
	'''
	!{params.java17} -jar !{params.gatk}GenomeAnalysisTK.jar -T BaseRecalibrator -nct 6 -R !{params.genome_ref} -I !{pair_tag}_realigned.bam -knownSites !{params.dbsnp} -knownSites !{params.Mills_indels} -knownSites !{params.ThousandG_indels} -BQSR !{pair_tag}_recal.table -L !{params.bed} -o !{pair_tag}_post_recal.table		
	'''
}

////////// STEP 09 ################### Base quality score recalibration - Generate before/after plots

process creation_recal_plots {
	echo "creation_recal_plots"
	tag { pair_tag }
	cpus 6
	clusterOptions = params.cluster_options
input:
	set val(pair_tag), file("${pair_tag}_recal.table") from recal_table_files2
 	set val(pair_tag), file("${pair_tag}_post_recal.table") from recal_table_post_files
output:
	set val(pair_tag), file("${pair_tag}_recalibration_plots.pdf") into recal_plots_files
shell:
	'''
	!{params.java17} -jar !{params.gatk}GenomeAnalysisTK.jar -T AnalyzeCovariates -R !{params.genome_ref} -before !{pair_tag}_recal.table -after !{pair_tag}_post_recal.table -plots !{pair_tag}_recalibration_plots.pdf			
	'''
}

////////// STEP 10 ################### Base quality score recalibration - Apply the recalibration to your sequence data

process creation_recal_bam {
	echo "creation_recal_bam"
	tag { pair_tag }
	cpus 6
	memory '16GB'
	clusterOptions = params.cluster_options
input:
	set val(pair_tag), file("${pair_tag}_realigned.bam"), file("${pair_tag}_realigned.bai") from realign_bam_bai_files3
	set val(pair_tag), file("${pair_tag}_recal.table") from recal_table_files3
output:
	set val(pair_tag), file("${pair_tag}_realigned_recal.bam") into recal_bam_files
shell:
	'''
	!{params.java17} -jar !{params.gatk}GenomeAnalysisTK.jar -T PrintReads -nct 6 -R !{params.genome_ref} -I !{pair_tag}_realigned.bam -BQSR !{pair_tag}_recal.table -L !{params.bed} -o !{pair_tag}_realigned_recal.bam			
	'''
}