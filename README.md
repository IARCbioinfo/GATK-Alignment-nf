# Alignment-GATK-BP-nf
Performs bwa alignment and pre-processing (realignment and recalibration) following GATK best practices

Before using the pipeline, specify in your config file the paths to the following files and softwares:
       genome_ref = '/path/ucsc.hg19.fasta'
		   dbsnp = '/path/GATK_Bundle/dbsnp_138.hg19.vcf'
		   Mills_indels='/path/GATK_Bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf'
		   1000G_indels='/path/GATK_Bundle/1000G_phase1.indels.hg19.sites.vcf'
		   bwa = "/appli57/bwa/bwa/-0-7-12/bwa"
		   java17 = "java"
		   picard = "/appli57/picard/picard-tools-1.131/picard.jar"
		   gatk = "/appli57/GenomeAnalysisTK/GATK-3.4-0/"
		   samtools = "/appli57/samtools/samtools-1.2/samtools"


		   

