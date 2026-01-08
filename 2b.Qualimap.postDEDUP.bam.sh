#!/bin/bash -l

################################################################################################################################
###		Bos taurus 30 samples PE illumina 2024 
###		Andreas E Lundberg; andreas.e.lundberg@slu.se; anlu4433@student.uu.se 
### 	
### 	Qualimap 
###		post MarkDuplicatesSpark, post alignment
###		ARS-UCD1.2_bosTau9
###		
###		
###		
################################################################################################################################

# exit 1 
# REMOVE BEFORE RUN SCRIPT


####### Expected paths under project root 
#			.
#			├── bin
#			├── data
#			│   ├── fastq
#			│   └── refseq
#			├── results
#			│   ├── bam
#			│   ├── stats
#			│   └── vcf
#			└── tmp
#
# symlink tmp directory to fast storage SSD/NVMe disk


####### Variables & paths
# absolute path to base directory of project
projectPath="/home/andreas/Bos_taurus.30bulls.illuminaPE"
binPath="/home/andreas/bin"


# make nicer and use 1.5x threads of node 
# threadsNode="60" 	# 1.5x 32 threads
# fastpThreads="6"
# sortThreads="14"
# bwaThreads="40" 		# "((${threadsNode}-${fastpThreads}-${sortThreads}))"

# fasta refeq file 
# refSeq="$projectPath/data/refseq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fasta"
			# refSeq="$projectPath/data/bwa-meme.refseq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fasta"


# indexPrefix="Bos_taurus.ARS-UCD1.2"

# fastqPath="$projectPath/data/fastq" 
bamPath="$projectPath/results/bam" 
# vcfPath="$projectPath/results/vcf"
statsPath="$projectPath/results/stats"
tmpPath="$projectPath/tmp" # use faster NVMe disk 

mkdir -p $tmpPath ${statsPath}/qualimap.postDEDUP.bam

################################################################################################################################
## MarkDuplicates - Spark parallelization 


##	For loop commands per Sample
for sampleName in `ls ${bamPath}/*.QC.SORTED.BWA-MEME.MarkDuplicatesSpark.DEDUP.bam | awk -F"/" '{print $NF}' | awk -F"." '{print $1}'`; do
	# sampleName variable example is: SH-2270-SE-402_S8_L003

	mkdir -p ${statsPath}/qualimap.postDEDUP.bam/${sampleName}

	# echo -e "${sampleName}
	# debug log output of variables
	# echo -e "${sampleName}\t$bamPath/${sampleName}" 2>&1 | tee -a $statsPath/${sampleName}.debug.log
	# echo -e "\t$fastqPath/${sampleName}_R1_001.fastq\t$bamPath/${sampleName}.sam" 2>&1 | tee -a $statsPath/${sampleName}.debug.log

	${binPath}/qualimap bamqc \
		-bam ${bamPath}/${sampleName}.QC.SORTED.BWA-MEME.MarkDuplicatesSpark.DEDUP.bam \
		-outformat PDF \
		-outdir ${statsPath}/qualimap.postDEDUP.bam/${sampleName} \
		-outfile ${sampleName}.QC.SORTED.BWA-MEME.MarkDuplicatesSpark.DEDUP.bam.postDEDUP.Qualimap.pdf \
		--paint-chromosome-limits \
		--java-mem-size=320G

#qualimap bamqc 
#-bam SH-2272-SE-432_S139_L003.QC.SORTED.BWA-MEME.MarkDuplicatesSpark.DEDUP.bam 
#-outformat PDF 
#-outdir ../SH-2272-SE-432_S139_L003.QC.SORTED.BWA-MEME.MarkDuplicatesSpark.DEDUP.bam.qualimap 
#-outfile SH-2272-SE-432_S139_L003.QC.SORTED.BWA-MEME.MarkDuplicatesSpark.DEDUP.bam.qualimap.pdf
# --java-mem-size=320G

done


