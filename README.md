# GATK-Alignment-nf
Performs bwa alignment and pre-processing (realignment and recalibration) following GATK best practices

The whole pipeline is made up of 10 steps: 

		1. Alignment generating a sam file

		2. Conversion to bam and sorting

		3. Marking of duplicates

		4. Indexing of bam

		5. & 6. Local realignment around indels

		7. to 10. Base quality score recalibration


Before using the pipeline, specify in your config file the paths to the following files and softwares:

		   genome_ref = '/appli/reference/GATK_Bundle/ucsc.hg19.fasta'
		   
		   dbsnp = '/appli/reference/GATK_Bundle/dbsnp_138.hg19.vcf'
		   
		   Mills_indels = '/appli/reference/GATK_Bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf'
		   
		   ThousandG_indels = '/appli/reference/GATK_Bundle/1000G_phase1.indels.hg19.sites.vcf'
		   
		   cosmic = '/appli/mutect/references/Cosmic_v73.hg19_noMT.vcf'
		   
		   bwa = '/appli/bwa/bwa-0-7-12/bwa'
		   
		   java17 = 'java'
		   
		   picard = '/appli/picard/picard-tools-1.131/picard.jar'
		   
		   gatk = '/appli/GenomeAnalysisTK/GATK-3.4-0/'
		   
		   samtools = '/appli/samtools/samtools-1.2/samtools'



		   

